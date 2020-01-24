/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Author
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "codeRules.H"
#include "WENOBase.H"
#include "WENOPolynomial.H"
#include "geometryWENO.H"
#include "SVD.H"
#include "processorFvPatch.H"
#include "labelListIOList.H"
#include "OFstream.H"
#include "IFstream.H"

#include <iostream>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::WENOBase::splitStencil
(
    const fvMesh& mesh,
    const label cellI,
    label& nStencilsI
)
{
    const pointField& pts = mesh.points();
    const cell& faces = mesh.cells()[cellI];

    List<List<scalarSquareMatrix> > JacobiInvQ(nStencilsI - 1);

    label exludeFace = 0;

    for
    (
        label stencilI = 1;
        stencilI < cellToPatchMap_[cellI].size();
        stencilI ++
    )
    {
        cellToPatchMap_[cellI][stencilI].setSize(1);
        cellToPatchMap_[cellI][stencilI][0] = -1;
    }

    // Fill lists with inverse jacobians J_Q for each subsector

    forAll(faces, faceI)
    {
        if (faces[faceI] < mesh.nInternalFaces())
        {
            List<tetIndices> faceTets =
                polyMeshTetDecomposition::faceTetIndices
                (
                    mesh,
                    faces[faceI],
                    cellI
                );

            triFaceList triFaces(faceTets.size());

            forAll(faceTets, cTI)
            {
                triFaces[cTI] = faceTets[cTI].faceTriIs(mesh);
            }

            JacobiInvQ[faceI - exludeFace].setSize(triFaces.size());

            forAll(triFaces, i)
            {
                const triFace& tri(triFaces[i]);

                scalarSquareMatrix J = geometryWENO::jacobi(
                    mesh.C()[cellI][0], mesh.C()[cellI][1],
                    mesh.C()[cellI][2],
                    pts[tri[0]][0], pts[tri[0]][1], pts[tri[0]][2],
                    pts[tri[1]][0], pts[tri[1]][1], pts[tri[1]][2],
                    pts[tri[2]][0], pts[tri[2]][1], pts[tri[2]][2]
                );

                JacobiInvQ[faceI - exludeFace][i] = Foam::geometryWENO::JacobiInverse(J);
            }
        }
        else
        {
            exludeFace++;
        }
    }

    // Distribute cells from big central stencil to the appropriate sector

    for (label cellJ = 1; cellJ <stencilsID_[cellI][0].size(); cellJ++)
    {
        bool attached = false;

        label actualFace = 0;

        while (attached == false)
        {
            if (stencilsID_[cellI][actualFace + 1].size() < 2.1*nDvt_)
            {
                forAll(JacobiInvQ[actualFace], triangleI)
                {
                    point transCenterJ = pTraits<point>::zero;

                    if (cellToPatchMap_[cellI][0][cellJ] == -1)
                    {
                        transCenterJ =
                            Foam::geometryWENO::transformPoint
                            (
                                JacobiInvQ[actualFace][triangleI],
                                mesh.C()[stencilsID_[cellI][0][cellJ]],
                                mesh.C()[cellI]
                            );
                    }
                    else if (cellToPatchMap_[cellI][0][cellJ] > -1)
                    {
                        transCenterJ =
                            Foam::geometryWENO::transformPoint
                            (
                                JacobiInvQ[actualFace][triangleI],
                                haloCenters_[cellToPatchMap_[cellI][0][cellJ]]
                                    [stencilsID_[cellI][0][cellJ]],
                                mesh.C()[cellI]
                            );
                    }

                    if
                    (
                        sign(transCenterJ.x()) > 0
                     && sign(transCenterJ.y()) > 0
                     && sign(transCenterJ.z()) > 0
                    )
                    {
                        stencilsID_[cellI][actualFace + 1].append
                        (
                            stencilsID_[cellI][0][cellJ]
                        );

                        cellToPatchMap_[cellI][actualFace + 1].append
                        (
                            cellToPatchMap_[cellI][0][cellJ]
                        );

                        attached = true;
                        break;
                    }
                }
            }

            if (actualFace < JacobiInvQ.size() - 1)
            {
                actualFace++;
            }
            else
            {
                break;
            }
        }
    }

    // Reject sectors without enough cells
    // and cut the stencils to the necessary size

    const scalar necSize = 2.0*nDvt_ + 1;

    forAll(stencilsID_[cellI], stencilI)
    {
        if (stencilsID_[cellI][stencilI].size() >= necSize)
        {
            stencilsID_[cellI][stencilI].resize(necSize);
            cellToPatchMap_[cellI][stencilI].resize(necSize);
        }
        else
        {
            stencilsID_[cellI][stencilI].resize(1);
            cellToPatchMap_[cellI][stencilI].resize(1);
            stencilsID_[cellI][stencilI][0] = -1;
            cellToPatchMap_[cellI][stencilI][0] = -4;

            nStencilsI--;
        }
    }
}


void Foam::WENOBase::extendStencils
(
    const fvMesh& mesh,
    const label cellI,
    labelList& lastNeighboursI,
    label& minStencilSize
)
{
    const label startEntry = stencilsID_[cellI][0].size() - lastNeighboursI[0];
    const label lastEntry = stencilsID_[cellI][0].size();

    lastNeighboursI[0] = 0;

    for (label nI = startEntry; nI < lastEntry; nI++)
    {
        const labelList& ngbhC = mesh.cellCells()[stencilsID_[cellI][0][nI]];

        forAll(ngbhC, I)
        {
            bool add = true;

            for (label J = stencilsID_[cellI][0].size() - 1; J > -1; J--)
            {
                if (stencilsID_[cellI][0][J] == ngbhC[I])
                {
                    add = false;
                    break;
                }
            }

            if (add == true)
            {
                 stencilsID_[cellI][0].append(ngbhC[I]);
                 lastNeighboursI[0] += 1;
            }
        }
    }

    lastNeighboursI.setSize(lastNeighboursI[0] + 1);

    for (label i = lastEntry; i < lastEntry+lastNeighboursI[0]; i++)
    {
        lastNeighboursI[i + 1 - lastEntry] = stencilsID_[cellI][0][i];
    }

    minStencilSize = stencilsID_[cellI][0].size();
}



void Foam::WENOBase::sortStencil
(
    const fvMesh& mesh,
    const label cellI,
    const label maxSize
)
{
    point transCcellI =
        Foam::geometryWENO::transformPoint
        (
            JInv_[cellI],
            mesh.C()[stencilsID_[cellI][0][0]],
            refPoint_[cellI]
        );

    scalarField distField(stencilsID_[cellI][0].size(), 0.0);
    scalarField numberField(stencilsID_[cellI][0].size(), cellI);
    scalarField mapField(stencilsID_[cellI][0].size(), -1);

    for (label i = 1; i < stencilsID_[cellI][0].size(); i++)
    {
        if (cellToPatchMap_[cellI][0][i] == -1)
        {
            point transCJ =
                Foam::geometryWENO::transformPoint
                (
                    JInv_[cellI],
                    mesh.C()[stencilsID_[cellI][0][i]],
                    refPoint_[cellI]
                );

            distField[i] = mag(transCJ - transCcellI);
        }
        else
        {
            point transCJ =
                Foam::geometryWENO::transformPoint
                (
                    JInv_[cellI],
                    haloCenters_[cellToPatchMap_[cellI][0][i]]
                        [stencilsID_[cellI][0][i]],
                    refPoint_[cellI]
                );

            distField[i] = mag(    transCJ - transCcellI);
        }

        numberField[i] = stencilsID_[cellI][0][i];
        mapField[i] = cellToPatchMap_[cellI][0][i];
    }

    // Bubblesort algorithm
    for (label p = 0; p < distField.size(); p++)
    {
        bool swapped = false;
        for (label j = 0; j < distField.size() - (p + 1); j++)
        {
            if (distField[j] > distField[j + 1])
            {
                scalar tempDist = distField[j];
                scalar tempID = numberField[j];
                scalar tempMap = mapField[j];
                distField[j] = distField[j + 1];
                numberField[j] = numberField[j + 1];
                mapField[j] = mapField[j + 1];
                distField[j + 1] = tempDist;
                numberField[j + 1] = tempID;
                mapField[j + 1] = tempMap;

                swapped = true;
            }
        }
        if (!swapped) break;
    }

    // Cut stencil to necessary size

    stencilsID_[cellI][0].resize(min(maxSize, numberField.size()));
    cellToPatchMap_[cellI][0].resize(stencilsID_[cellI][0].size());

    for (label i = 0; i < stencilsID_[cellI].size(); i++)
    {
        stencilsID_[cellI][0][i] = numberField[i];
        cellToPatchMap_[cellI][0][i] = mapField[i];
    }
}


void Foam::WENOBase::distributeStencils
(
    const fvMesh& mesh,
    labelListList& haloCells,
    List<List<List<point> > >& haloTriFaceCoord
)
{
#ifdef FOAM_PSTREAM_COMMSTYPE_IS_ENUMCLASS 
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
#else
    PstreamBuffers pBufs(Pstream::nonBlocking);
#endif

    // Collect centres of halo cells
    forAll(haloCells, patchI)
    {
        forAll(haloCells[patchI], cellI)
        {
            haloCenters_[patchI].append
            (
                mesh.C()[haloCells[patchI][cellI]]
            );
        }
    }

    // Distribute halo cell ID's
    // Assigning new ID, starting at 0
    forAll(patchToProcMap_, patchI)
    {
        if (patchToProcMap_[patchI] != -1)
        {
            UOPstream toBuffer(patchToProcMap_[patchI], pBufs);
            toBuffer << haloCells[patchI];
        }
    }

    pBufs.finishedSends();

    forAll(patchToProcMap_, patchI)
    {
        haloCells[patchI].clear();

        if (patchToProcMap_[patchI] != -1)
        {
            UIPstream fromBuffer(patchToProcMap_[patchI], pBufs);
            fromBuffer >> haloCells[patchI];
        }
    }

    forAll(haloCells, patchI)
    {
        if (patchToProcMap_[patchI] != -1)
        {
            label newID = 0;
            forAll(haloCells[patchI], ID)
            {
                haloCells[patchI][ID] = newID;
                newID++;
            }
        }
    }


    // Distribute halo cell center coordinates
    forAll(patchToProcMap_, patchI)
    {
        if (patchToProcMap_[patchI] != -1)
        {
            UOPstream toBuffer(patchToProcMap_[patchI], pBufs);
            toBuffer << haloCenters_[patchI];
        }
    }

    pBufs.finishedSends();

    forAll(patchToProcMap_, patchI)
    {
        haloCenters_[patchI].clear();

        if (patchToProcMap_[patchI] != -1)
        {
            UIPstream fromBuffer(patchToProcMap_[patchI], pBufs);
            fromBuffer >> haloCenters_[patchI];
        }
    }


    // Distribute coordinates of halo faces for calculating volume integrals
    forAll(patchToProcMap_, patchI)
    {
        if (patchToProcMap_[patchI] != -1)
        {
            UOPstream toBuffer(patchToProcMap_[patchI], pBufs);
            toBuffer << haloTriFaceCoord[patchI];
        }
    }

    pBufs.finishedSends();

    forAll(patchToProcMap_, patchI)
    {
        haloTriFaceCoord[patchI].clear();

        if (patchToProcMap_[patchI] != -1)
        {
            UIPstream fromBuffer(patchToProcMap_[patchI], pBufs);
            fromBuffer >> haloTriFaceCoord[patchI];
        }
    }
}


Foam::scalarRectangularMatrix Foam::WENOBase::calcMatrix
(
    const fvMesh& mesh,
    const label cellI,
    const label stencilI,
    const List<List<List< point> > >& haloTriFaceCoord
)
{
    const label stencilSize = stencilsID_[cellI][stencilI].size();

    scalarRectangularMatrix A
    (
        stencilSize - 1,
        nDvt_,
        scalar(0.0)
    );

    point transCenterI = Foam::geometryWENO::transformPoint
    (
        JInv_[cellI],
        mesh.C()[cellI],
        refPoint_[cellI]
    );

    volIntegralType volIntegralsIJ = volIntegralsList_[cellI];

    // Add one line per cell
    for (label cellJ = 1; cellJ < stencilSize; cellJ++)
    {
        if (cellToPatchMap_[cellI][stencilI][cellJ] == -1)
        {
            point transCenterJ =
                Foam::geometryWENO::transformPoint
                (
                    JInv_[cellI],
                    mesh.C()[stencilsID_[cellI][stencilI][cellJ]],
                    refPoint_[cellI]
                );

            volIntegralType transVolMom =
                Foam::geometryWENO::transformIntegral
                (
                    mesh,
                    stencilsID_[cellI][stencilI][cellJ],
                    transCenterJ,
                    polOrder_,
                    JInv_[cellI],
                    refPoint_[cellI],
                    refDet_[cellI]
                );

            for (label n = 0; n <= dimList_[cellI][0]; n++)
            {
                for (label m = 0; m <= dimList_[cellI][1]; m++)
                {
                    for (label l = 0; l <= dimList_[cellI][2]; l++)
                    {
                        if ((n + m + l) <= polOrder_ && (n + m + l) > 0)
                        {
                            volIntegralsIJ[n][m][l] =
                                calcGeom
                                (
                                    transCenterJ - transCenterI,
                                    n,
                                    m,
                                    l,
                                    transVolMom,
                                    volIntegralsList_[cellI]
                                );
                        }
                    }
                }
            }

        }
        else if (cellToPatchMap_[cellI][stencilI][cellJ] > -1)
        {
            point transCenterJ =
                Foam::geometryWENO::transformPoint
                (
                    JInv_[cellI],
                    haloCenters_[cellToPatchMap_[cellI][stencilI][cellJ]]
                        [stencilsID_[cellI][stencilI][cellJ]],
                    refPoint_[cellI]
                );

            volIntegralType transVolMom =
                Foam::geometryWENO::getHaloMoments
                (
                    mesh,
                    transCenterJ,
                    haloTriFaceCoord[cellToPatchMap_[cellI][stencilI][cellJ]]
                        [stencilsID_[cellI][stencilI][cellJ]],
                    polOrder_,
                    JInv_[cellI],
                    refPoint_[cellI]
                );

            for (label n = 0; n <= dimList_[cellI][0]; n++)
            {
                for (label m = 0; m <= dimList_[cellI][1]; m++)
                {
                    for (label l = 0; l <= dimList_[cellI][2]; l++)
                    {
                        if ((n + m + l) <= polOrder_ && (n + m + l) > 0)
                        {
                            volIntegralsIJ[n][m][l] =
                                calcGeom
                                (
                                    transCenterJ - transCenterI,
                                    n,
                                    m,
                                    l,
                                    transVolMom,
                                    volIntegralsList_[cellI]
                                );
                        }
                    }
                }
            }
        }

        WENOPolynomial::addCoeffs
        (
            A[cellJ - 1],
            polOrder_,
            dimList_[cellI],
            volIntegralsIJ
        );
    }

    // Returning pseudoinverse using SVD

    SVD svd(A, 1e-5);

    return svd.VSinvUt();
}


Foam::scalar Foam::WENOBase::calcGeom
(
    const vector x_ij,
    const label n,
    const label m,
    const label o,
    const volIntegralType& volMomJ,
    const volIntegralType& volMomI
)
{
    scalar geom = 0.0;

    for (label k = 0; k <= n; k++)
    {
        for (label l = 0; l <= m; l++)
        {
            for (label j = 0; j <= o; j++)
            {
                geom +=
                    factorial(n)/(factorial(k)*factorial((n - k)))
                   *factorial(m)/(factorial(l)*factorial((m - l)))
                   *factorial(o)/(factorial(j)*factorial((o - j)))
                   *pow(x_ij.x(), k)*pow(x_ij.y(), l)*pow(x_ij.z(), j)
                   *volMomJ[n - k][m - l][o - j];
            }
        }
    }

    return (geom - volMomI[n][m][o]);
}


Foam::WENOBase::WENOBase
(
    const fvMesh& mesh,
    const label polOrder
)
{
    polOrder_ = polOrder;

    Dir_ = mesh.time().path()/"constant"/"WENOBase" + Foam::name(polOrder_);

    // 3D version
    if (mesh.nSolutionD() == 3)
    {
        nDvt_ = (polOrder_ + 1.0)*(polOrder_ + 2.0)*(polOrder_ + 3.0)/6.0 - 1.0;
    }
    else // 2D version (always only one cell in z-direction)
    {
        nDvt_ = (polOrder_ + 1.0)*(polOrder_ + 2.0)/2.0 - 1.0;
    }

    // Set the dimList
    labelList dummyList(3,0);
    dimList_.setSize(mesh.nCells(),dummyList);
    
    // Vector with valid dimensions
    vector dimMesh = mesh.solutionD();
    
    for (label i = 0; i < mesh.nCells(); i++)
    {
        dimList_[i][0] = (dimMesh[0] == 1 ? polOrder_ : 0);
        dimList_[i][1] = (dimMesh[1] == 1 ? polOrder_ : 0);
        dimList_[i][2] = (dimMesh[2] == 1 ? polOrder_ : 0);
    }


    // Check for existing lists
    bool listExist = readList(mesh);

    // Create new lists if necessary
    if (listExist == false)
    {
        // Read expert factor

        IOdictionary WENODict
        (
            IOobject
            (
                "WENODict",
                mesh.time().caseSystem(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );

        const scalar extendRatio =
            WENODict.lookupOrDefault<scalar>("extendRatio", 2.5);

        // Get big central stencils

        stencilsID_.setSize(mesh.nCells());
        cellToPatchMap_.setSize(mesh.nCells());

        labelList nStencils(mesh.nCells(),0);

        volIntegralType volIntegrals;

        volIntegrals.resize((polOrder_ + 1));

        for (label i = 0; i < (polOrder_+1); i++)
        {
            volIntegrals[i].resize((polOrder_+ 1)-i);

            for (label j = 0; j < ((polOrder_+1)-i); j++)
            {
                volIntegrals[i][j].resize((polOrder_ + 1)-i, 0.0);
            }
        }


        volIntegralsList_.setSize(mesh.nCells(), volIntegrals);

        labelListList lastNeighbours(mesh.nCells());

        const fvPatchList& patches = mesh.boundary();

        patchToProcMap_.setSize(patches.size(), -1);

        labelListList haloCells(patches.size());
        List<List<List<point> > > haloTriFaceCoord(patches.size());

        haloCenters_.setSize(patches.size());

        ownHalos_ = haloCells;

        JInv_.setSize(mesh.nCells());
        refPoint_.setSize(mesh.nCells());
        refDet_.setSize(mesh.nCells());

        for (label cellI = 0; cellI < mesh.nCells(); cellI++)
        {
            const cell& faces = mesh.cells()[cellI];

            nStencils[cellI] = 1;

            forAll(faces, faceI)
            {
                if (faces[faceI] < mesh.nInternalFaces())
                {
                    nStencils[cellI]++;
                }
            }

            stencilsID_[cellI].setSize(nStencils[cellI]);
            cellToPatchMap_[cellI].setSize(nStencils[cellI]);

            forAll(stencilsID_[cellI],stencilI)
            {
                stencilsID_[cellI][stencilI].append(cellI);
            }
            stencilsID_[cellI][0].append(mesh.cellCells()[cellI]);

            lastNeighbours[cellI].setSize(stencilsID_[cellI][0].size());
            lastNeighbours[cellI][0] = lastNeighbours[cellI].size() - 1;

            for (label i = 1; i < stencilsID_[cellI][0].size(); i++)
            {
                lastNeighbours[cellI][i] = stencilsID_[cellI][0][i];
            }

            Foam::geometryWENO::initIntegrals
            (
                mesh,
                cellI,
                polOrder_,
                volIntegralsList_[cellI],
                JInv_[cellI],
                refPoint_[cellI],
                refDet_[cellI]
            );


            // Extend central stencil to neccessary size

            label minStencilSize = 0;

            while (minStencilSize < 1.2*extendRatio*nDvt_*nStencils[cellI])
            {
                extendStencils
                (
                    mesh,
                    cellI,
                    lastNeighbours[cellI],
                    minStencilSize
                );
            }

            // Sort and cut stencil

            labelList dummyLabels(stencilsID_[cellI][0].size(),-1);
            cellToPatchMap_[cellI][0] = dummyLabels;

            sortStencil(mesh,cellI, extendRatio*nDvt_*nStencils[cellI]);
        }

        // Extension to halo cells, if neccessary

        if(Pstream::parRun())
        {
            labelListList stencilNeedsHalo(mesh.nCells());

            forAll(stencilNeedsHalo, i)
            {
                stencilNeedsHalo[i].setSize(patches.size(), -1);
            }

            forAll(patches, patchI)
            {
                if(isA<processorFvPatch>(patches[patchI]))
                {
                    labelList faceCells = patches[patchI].faceCells();

                    patchToProcMap_[patchI] =
                        refCast<const processorFvPatch>
                        (patches[patchI]).neighbProcNo();

                    forAll(faceCells, cellI)
                    {
                        // Add halo cells and mark stencils needing halo cells
                        forAll(stencilsID_[faceCells[cellI]][0], cellJ)
                        {
                            label haloCell =
                                stencilsID_[faceCells[cellI]][0][cellJ];

                            stencilNeedsHalo[haloCell][patchI] = patchI;

                            bool add = true;

                            forAll(haloCells[patchI], cellK)
                            {
                                if (haloCells[patchI][cellK] == haloCell)
                                {
                                    add = false;
                                    break;
                                }
                            }
                            if (add == true)
                            {
                                haloCells[patchI].append(haloCell);
                            }
                        }
                    }
                }
            }

            forAll(haloTriFaceCoord, patchI)
            {
                haloTriFaceCoord[patchI].setSize(haloCells[patchI].size());

                forAll(haloCells[patchI], cellI)
                {
                    haloTriFaceCoord[patchI][cellI] =
                        Foam::geometryWENO::getTriFaces
                        (
                            mesh,
                            haloCells[patchI][cellI]
                        );
                }
            }


            // Distribute halo cells, triFaces and centres
            // New cell ID's begin behind the local cells
            ownHalos_ = haloCells;

            distributeStencils
            (
                mesh,
                haloCells,
                haloTriFaceCoord
            );

            // Add halo cells to stencils
            forAll(stencilsID_, stencilI)
            {
                forAll(stencilNeedsHalo[stencilI], patchI)
                {
                    if (stencilNeedsHalo[stencilI][patchI] != -1)
                    {
                        scalar radius =
                            mag
                            (
                                mesh.C()[stencilsID_[stencilI][0][0]]
                              - mesh.C()[stencilsID_[stencilI][0]
                                    [stencilsID_[stencilI][0].size()-1]]
                            );

                        labelList haloLayer =
                            haloCells[stencilNeedsHalo[stencilI][patchI]];

                        List<point> centers =
                            haloCenters_[stencilNeedsHalo[stencilI][patchI]];

                        forAll(haloLayer,i)
                        {
                            scalar radiusI =
                                mag
                                (
                                    mesh.C()[stencilsID_[stencilI][0][0]]
                                  - centers[i]
                                );

                            if (radiusI <= radius)
                            {
                                stencilsID_[stencilI][0].append(haloLayer[i]);
                                cellToPatchMap_[stencilI][0].append(patchI);
                            }
                        }
                    }
                }
            }

            // Get final big central stencils
            for (label cellI = 0; cellI < mesh.nCells(); cellI++)
            {
                sortStencil(mesh, cellI, (extendRatio*nDvt_)*nStencils[cellI]);
            }
        }

        // Split the stencil in several sectorial stencils
        for (label cellI = 0; cellI < mesh.nCells(); cellI++)
        {
            splitStencil(mesh, cellI, nStencils[cellI]);
        }


        // Get the least squares matrices and their pseudoinverses

        LSmatrix_.setSize(mesh.nCells());

        for (label cellI = 0; cellI < mesh.nCells(); cellI++)
        {
            label excludeFace = 0;

            LSmatrix_[cellI].setSize(nStencils[cellI]);

            forAll(stencilsID_[cellI], stencilI)
            {
                if (stencilsID_[cellI][stencilI][0] != -1)
                {
                    LSmatrix_[cellI][stencilI - excludeFace] =
                        calcMatrix
                        (
                            mesh,
                            cellI,
                            stencilI,
                            haloTriFaceCoord
                        );
                }
                else
                {
                    excludeFace++;
                }
            }

        }

        // Get the smoothness indicator matrices

        B_.setSize(mesh.nCells());

        for(label cellI = 0; cellI < mesh.nCells(); cellI++)
        {
            B_[cellI] =
                Foam::geometryWENO::getB
                (
                    mesh,
                    cellI,
                    polOrder_,
                    nDvt_,
                    JInv_[cellI],
                    refPoint_[cellI],
                    dimList_[cellI]
                );

        }

        // Get surface integrals over basis functions in transformed coordinates

        intBasTrans_.setSize(mesh.nFaces());

        scalarList dummy(2,0.0);
        refFacAr_.setSize(mesh.nFaces(),dummy);

        for (label faceI = 0; faceI < mesh.nFaces(); faceI++)
        {
            intBasTrans_[faceI].setSize(2);
            intBasTrans_[faceI][0] = volIntegrals;
            intBasTrans_[faceI][1] = volIntegrals;
        }

        Foam::geometryWENO::surfIntTrans
        (
            mesh,
            polOrder_,
            volIntegralsList_,
            JInv_,
            refPoint_,
            intBasTrans_,
            refFacAr_
        );

        // Write Lists to constant folder
        writeList
        (
            mesh
        );
    }
}


bool Foam::WENOBase::readList
(
    const fvMesh& mesh
)
{
    if(isDir(Dir_))
    {
        Info<< "\nRead existing lists from constant folder \n" << endl;

        IFstream isDL(Dir_/"DimLists");
        dimList_.setSize(mesh.nCells());

        forAll(dimList_, cellI)
        {
            isDL >> dimList_[cellI];
        }

        IFstream isSID(Dir_/"StencilIDs");
        stencilsID_.setSize(mesh.nCells());
        scalar nEntries;

        for (label cellI = 0; cellI < mesh.nCells(); cellI++)
        {
            isSID >> nEntries;

            stencilsID_[cellI].setSize(nEntries);

            for (label stencilI = 0; stencilI < nEntries; stencilI++)
            {
                isSID >> stencilsID_[cellI][stencilI];
            }
        }

        IFstream isCToP(Dir_/"CellToPatchMaps");
        cellToPatchMap_.setSize(mesh.nCells());

        for (label cellI = 0; cellI < mesh.nCells(); cellI++)
        {
            isCToP >> nEntries;

            cellToPatchMap_[cellI].setSize(nEntries);

            for (label stencilI = 0; stencilI < nEntries; stencilI++)
            {
                isCToP >> cellToPatchMap_[cellI][stencilI];
            }
        }

        IFstream isLS(Dir_/"Pseudoinverses");
        LSmatrix_.setSize(mesh.nCells());

        for (label cellI = 0; cellI < mesh.nCells(); cellI++)
        {
            isLS >> nEntries;

            LSmatrix_[cellI].setSize(nEntries);

            for (label stencilI = 0; stencilI < nEntries; stencilI++)
            {
                isLS >> LSmatrix_[cellI][stencilI];
            }
        }

        IFstream isB(Dir_/"B");
        B_.setSize(mesh.nCells());

        forAll(B_, cellI)
        {
            isB >> B_[cellI];
        }

        const fvPatchList& patches = mesh.boundary();

        patchToProcMap_.setSize(patches.size());
        IFstream isPToP(Dir_/"PatchToProcMaps");

        forAll(patchToProcMap_, patchI)
        {
            isPToP >> patchToProcMap_[patchI];
        }

        ownHalos_.setSize(patches.size());
        IFstream isOH(Dir_/"OwnHalos");

        forAll(ownHalos_, patchI)
        {
            isOH >> nEntries;

            ownHalos_[patchI].setSize(nEntries);

            forAll(ownHalos_[patchI], cellI)
            {
                isOH >> ownHalos_[patchI][cellI];
            }
        }

        haloCenters_.setSize(patches.size());
        IFstream isHalo(Dir_/"HaloCenters");

        forAll(haloCenters_, patchI)
        {
            isHalo >> nEntries;

            haloCenters_[patchI].setSize(nEntries);

            forAll(haloCenters_[patchI], cellI)
            {
                isHalo >> haloCenters_[patchI][cellI];
            }
        }

        // Calculating volume integrals in transformed coordinates,
        // faster than writting and reading


        volIntegralType volIntegrals;

        volIntegrals.resize((polOrder_+1));

        for (label i = 0; i < (polOrder_+1); i++)
        {
            volIntegrals[i].resize((polOrder_+ 1)-i);

            for (label j = 0; j < ((polOrder_+1)-i); j++)
            {
                volIntegrals[i][j].resize((polOrder_ + 1)-i, 0.0);
            }
        }

        volIntegralsList_.setSize(mesh.nCells(),volIntegrals);
        JInv_.setSize(mesh.nCells());
        refPoint_.setSize(mesh.nCells());
        refDet_.setSize(mesh.nCells());

        for (label cellI = 0; cellI < mesh.nCells(); cellI++)
        {
            Foam::geometryWENO::initIntegrals
            (
                mesh,
                cellI,
                polOrder_,
                volIntegralsList_[cellI],
                JInv_[cellI],
                refPoint_[cellI],
                refDet_[cellI]
            );
        }

        // Get surface integrals in transformed coordinates

        intBasTrans_.setSize(mesh.nFaces());

        scalarList dummy(2,0.0);
        refFacAr_.setSize(mesh.nFaces(),dummy);

        for (label faceI = 0; faceI < mesh.nFaces(); faceI++)
        {
            intBasTrans_[faceI].setSize(2);
            intBasTrans_[faceI][0] = volIntegrals;
            intBasTrans_[faceI][1] = volIntegrals;
        }

        Foam::geometryWENO::surfIntTrans
        (
            mesh,
            polOrder_,
            volIntegralsList_,
            JInv_,
            refPoint_,
            intBasTrans_,
            refFacAr_
        );

        return true;
    }
    else
    {
        Info<< "Create new lists \n" << endl;

        mkDir(Dir_);

        return false;
    }
}


void Foam::WENOBase::writeList
(
    const fvMesh& mesh
)
{
    Info<< "Write created lists to constant folder \n" << endl;

    OFstream osPToP(Dir_/"PatchToProcMaps");

    forAll(patchToProcMap_, i)
    {
        osPToP<< patchToProcMap_[i] << endl;
    }

    OFstream osDL(Dir_/"DimLists");

    forAll(dimList_, cellI)
    {
        osDL<< dimList_[cellI] << endl;
    }

    OFstream osSID(Dir_/"StencilIDs");

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        osSID<< stencilsID_[cellI].size() << endl;

        forAll(stencilsID_[cellI], stenciI)
        {
            osSID<< stencilsID_[cellI][stenciI] << endl;
        }
    }

    OFstream osCToP(Dir_/"CellToPatchMaps");

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        osCToP<< cellToPatchMap_[cellI].size() << endl;

        forAll(cellToPatchMap_[cellI], cellJ)
        {
            osCToP<< cellToPatchMap_[cellI][cellJ] << endl;
        }
    }

    OFstream osLS(Dir_/"Pseudoinverses");
    osLS.precision(10);

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        osLS<< LSmatrix_[cellI].size() << endl;

        forAll(LSmatrix_[cellI], stenciI)
        {
            osLS<< LSmatrix_[cellI][stenciI] << endl;
        }
    }

    OFstream osHalo(Dir_/"HaloCenters");
    osHalo.precision(10);

    forAll(haloCenters_, patchI)
    {
        osHalo<< haloCenters_[patchI].size() << endl;

        forAll(haloCenters_[patchI], cellI)
        {
            osHalo<< haloCenters_[patchI][cellI] << endl;
        }
    }

    OFstream osOH(Dir_/"OwnHalos");
    osOH.precision(10);

    forAll(ownHalos_, patchI)
    {
        osOH<< ownHalos_[patchI].size() << endl;

        forAll(ownHalos_[patchI], cellI)
        {
            osOH<< ownHalos_[patchI][cellI] << endl;
        }
    }

    OFstream osB(Dir_/"B");
    osB.precision(10);

    forAll(B_, cellI)
    {
        osB<< B_[cellI] << endl;
    }


    /*
    label fx = 1;
    OFstream::streamFormat ofmt = OFstream::BINARY;
    OFstream osBbin(Dir_/"Bbinary", ofmt);
    osBbin.write(fx);


    label test = 22;

    std::ofstream outStream("yourFile", std::ios::binary);
    binary_write(outStream, test);
    outStream.close();

    std::ifstream inStream("yourFile", std::ios::binary);
    binary_read(inStream, test);
    inStream.close();
    Info<<test<<endl;
    */
}


// ************************************************************************* //
