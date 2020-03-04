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
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "codeRules.H"
#include "WENOBase.H"
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
    const fvMesh& globalMesh,
    const fvMesh& localMesh,
    const label localCellI,
    const label globalCellI,
    label& nStencilsI
)
{
    const pointField& pts = globalMesh.points();
    const cell& faces = globalMesh.cells()[globalCellI];

    List<List<scalarSquareMatrix> > JacobiInvQ(nStencilsI - 1);

    label exludeFace = 0;

    for
    (
        label stencilI = 1;
        stencilI < cellToProcMap_[localCellI].size();
        stencilI ++
    )
    {
        cellToProcMap_[localCellI][stencilI].setSize(1);
        cellToProcMap_[localCellI][stencilI][0] = int(Cell::local);
    }

    // Fill lists with inverse jacobians J_Q for each subsector
    forAll(faces, faceI)
    {
        if (faces[faceI] < globalMesh.nInternalFaces())
        {
            List<tetIndices> faceTets =
                polyMeshTetDecomposition::faceTetIndices
                (
                    globalMesh,
                    faces[faceI],
                    globalCellI
                );

            triFaceList triFaces(faceTets.size());

            forAll(faceTets, cTI)
            {
                triFaces[cTI] = faceTets[cTI].faceTriIs(globalMesh);
            }

            JacobiInvQ[faceI - exludeFace].setSize(triFaces.size());

            forAll(triFaces, i)
            {
                const triFace& tri(triFaces[i]);

                scalarSquareMatrix J = geometryWENO::jacobi(
                    globalMesh.C()[globalCellI][0],
                    globalMesh.C()[globalCellI][1],
                    globalMesh.C()[globalCellI][2],
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
    for (label cellJ = 1; cellJ <stencilsID_[localCellI][0].size(); cellJ++)
    {
        bool attached = false;

        label actualFace = 0;

        while (attached == false)
        {
            if (stencilsID_[localCellI][actualFace + 1].size() < 2.1*nDvt_)
            {
                forAll(JacobiInvQ[actualFace], triangleI)
                {
                    point transCenterJ = pTraits<point>::zero;

                    if (cellToProcMap_[localCellI][0][cellJ] == int(Cell::local))
                    {
                        transCenterJ =
                            Foam::geometryWENO::transformPoint
                            (
                                JacobiInvQ[actualFace][triangleI],
                                localMesh.C()[stencilsID_[localCellI][0][cellJ]],
                                globalMesh.C()[globalCellI]
                            );
                    }
                    else if (cellToProcMap_[localCellI][0][cellJ] != int(Cell::local) 
                          && cellToProcMap_[localCellI][0][cellJ] != int(Cell::deleted))
                    {
                        transCenterJ =
                            Foam::geometryWENO::transformPoint
                            (
                                JacobiInvQ[actualFace][triangleI],
                                haloCenters_[cellToProcMap_[localCellI][0][cellJ]]
                                    [stencilsID_[localCellI][0][cellJ]],
                                globalMesh.C()[globalCellI]
                            );
                    }

                    if
                    (
                        sign(transCenterJ.x()) > 0
                     && sign(transCenterJ.y()) > 0
                     && sign(transCenterJ.z()) > 0
                    )
                    {
                        stencilsID_[localCellI][actualFace + 1].append
                        (
                            stencilsID_[localCellI][0][cellJ]
                        );

                        cellToProcMap_[localCellI][actualFace + 1].append
                        (
                            cellToProcMap_[localCellI][0][cellJ]
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
    //Info << "stencil["<<localCellI<<"] and necSize "<<necSize<<": "<<stencilsID_[localCellI] <<nl<< endl;
    forAll(stencilsID_[localCellI], stencilI)
    {
        if (stencilsID_[localCellI][stencilI].size() >= necSize)
        {
            stencilsID_[localCellI][stencilI].resize(necSize);
            cellToProcMap_[localCellI][stencilI].resize(necSize);
        }
        else
        {
            stencilsID_[localCellI][stencilI].resize(1);
            cellToProcMap_[localCellI][stencilI].resize(1);
            stencilsID_[localCellI][stencilI][0] = int(Cell::deleted);
            cellToProcMap_[localCellI][stencilI][0] = int(Cell::deleted);

            nStencilsI--;
        }
    }
}


void Foam::WENOBase::extendStencils
(
    const fvMesh& mesh,
    const label cellI,
    std::map<int,bool>& stencilCells,
    label& minStencilSize
)
{
    std::map<int,bool> temp = stencilCells;
    
    // Iterate over all current neighbours and add there cells to the stencil
    for (auto it : stencilCells)
    {
        // If it is a neighbour cell
        if (it.second == true)
        {
            
            // Reference for cellIndex
            const int& cellInd = it.first;
            
            const labelList& ngbhC = mesh.cellCells()[cellInd];
            
            forAll(ngbhC,i)
            {
                if (temp.find(ngbhC[i]) == temp.end())
                {
                    stencilsID_[cellI][0].append(ngbhC[i]);
                    temp.insert(std::pair<int,bool>(ngbhC[i],true));
                }
            }
            
            // Set iterator to false for the next iteration
            auto itTemp = temp.find(it.first);
            itTemp->second = false;
        }
    }
    
    stencilCells.clear();
    stencilCells = temp;

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
    scalarField mapField(stencilsID_[cellI][0].size(), int(Cell::local));

    for (label i = 1; i < stencilsID_[cellI][0].size(); i++)
    {
        if (cellToProcMap_[cellI][0][i] == int(Cell::local))
        {
            //Pout << "cellID "<< stencilsID_[cellI][0][i] 
            //     << "  mesh.C(): " << mesh.C()[stencilsID_[cellI][0][i]]<<endl;
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
                    haloCenters_[cellToProcMap_[cellI][0][i]]
                        [stencilsID_[cellI][0][i]],
                    refPoint_[cellI]
                );

            distField[i] = mag(transCJ - transCcellI);
        }

        numberField[i] = stencilsID_[cellI][0][i];
        mapField[i] = cellToProcMap_[cellI][0][i];
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
    cellToProcMap_[cellI][0].resize(stencilsID_[cellI][0].size());

    for (label i = 0; i < stencilsID_[cellI].size(); i++)
    {
        stencilsID_[cellI][0][i] = numberField[i];
        cellToProcMap_[cellI][0][i] = mapField[i];
    }
}


void Foam::WENOBase::distributeStencils
(
    labelListList& haloCells
)
{
    // Get the procList of all other processors
    List<labelList> allValues(Pstream::nProcs());

    allValues[Pstream::myProcNo()] = receiveProcList_;

    Pstream::gatherList(allValues);
    
    Pstream::scatterList(allValues);

    // Clear old list and fill with -1
    sendProcList_.setSize(Pstream::nProcs());
    forAll (sendProcList_,i)
    {
        sendProcList_[i] = -1;
    }

    // Loop over ProcList
    forAll(allValues,procI)
    {
        if (procI == Pstream::myProcNo())
            continue;
        
        forAll(allValues[procI],i)
        {
            if (allValues[procI][i] == Pstream::myProcNo())
                sendProcList_[procI] = procI;
        }
    }
    
    
    
    #ifdef FOAM_PSTREAM_COMMSTYPE_IS_ENUMCLASS 
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    #else
        PstreamBuffers pBufs(Pstream::nonBlocking);
    #endif

    // Distribute halo cell ID's
    // Assigning new ID, starting at 0
    forAll(sendProcList_, procI)
    {
        if (sendProcList_[procI] != -1)
        {
            UOPstream toBuffer(sendProcList_[procI], pBufs);
            toBuffer << haloCells[procI];
        }
    }

    pBufs.finishedSends();

    forAll(receiveProcList_, procI)
    {
        haloCells[procI].clear();
        if (receiveProcList_[procI] != -1)
        {
            UIPstream fromBuffer(receiveProcList_[procI], pBufs);
            fromBuffer >> haloCells[procI];
        }
    }
    
}


Foam::scalarRectangularMatrix Foam::WENOBase::calcMatrix
(
    const fvMesh& globalMesh,
    const fvMesh& localMesh,
    const label   localCellI,
    const label stencilI,
    const List<List<List< point> > >& haloTriFaceCoord
)
{
    const label stencilSize = stencilsID_[localCellI][stencilI].size();

    scalarRectangularMatrix A
    (
        stencilSize - 1,
        nDvt_,
        scalar(0.0)
    );

    point transCenterI = Foam::geometryWENO::transformPoint
    (
        JInv_[localCellI],
        localMesh.C()[localCellI],
        refPoint_[localCellI]
    );

    volIntegralType volIntegralsIJ = volIntegralsList_[localCellI];

    // Add one line per cell
    for (label cellJ = 1; cellJ < stencilSize; cellJ++)
    {
        if (cellToProcMap_[localCellI][stencilI][cellJ] == int(Cell::local))
        {
            point transCenterJ =
                Foam::geometryWENO::transformPoint
                (
                    JInv_[localCellI],
                    localMesh.C()[stencilsID_[localCellI][stencilI][cellJ]],
                    refPoint_[localCellI]
                );

            volIntegralType transVolMom =
                Foam::geometryWENO::transformIntegral
                (
                    localMesh,
                    stencilsID_[localCellI][stencilI][cellJ],
                    transCenterJ,
                    polOrder_,
                    JInv_[localCellI],
                    refPoint_[localCellI],
                    refDet_[localCellI]
                );

            for (label n = 0; n <= dimList_[localCellI][0]; n++)
            {
                for (label m = 0; m <= dimList_[localCellI][1]; m++)
                {
                    for (label l = 0; l <= dimList_[localCellI][2]; l++)
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
                                    volIntegralsList_[localCellI]
                                );
                        }
                    }
                }
            }

        }
        else if (cellToProcMap_[localCellI][stencilI][cellJ] != int(Cell::local)
              && cellToProcMap_[localCellI][stencilI][cellJ] != int(Cell::deleted))
        {
            point transCenterJ =
                Foam::geometryWENO::transformPoint
                (
                    JInv_[localCellI],
                    haloCenters_[cellToProcMap_[localCellI][stencilI][cellJ]]
                        [stencilsID_[localCellI][stencilI][cellJ]],
                    refPoint_[localCellI]
                );
                

            volIntegralType transVolMom =
                Foam::geometryWENO::getHaloMoments
                (
                    transCenterJ,
                    haloTriFaceCoord[cellToProcMap_[localCellI][stencilI][cellJ]]
                        [stencilsID_[localCellI][stencilI][cellJ]],
                    polOrder_,
                    JInv_[localCellI],
                    refPoint_[localCellI]
                );

            for (label n = 0; n <= dimList_[localCellI][0]; n++)
            {
                for (label m = 0; m <= dimList_[localCellI][1]; m++)
                {
                    for (label l = 0; l <= dimList_[localCellI][2]; l++)
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
                                    volIntegralsList_[localCellI]
                                );
                        }
                    }
                }
            }
        }
        
        // Populate the matrix A
        addCoeffs(A,cellJ,polOrder_,dimList_[localCellI],volIntegralsIJ);
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


void Foam::WENOBase::addCoeffs
(
    scalarRectangularMatrix& A,
    const label cellj,
    const label polOrder,
    const labelList& dim,
    const volIntegralType& volIntegralsIJ 
)
{
    label currIdx = 0;
    for (label n=0 ; n <= dim[0] ; n++)
    {
        for (label m=0 ; m <= dim[1] ; m++)
        {
            for (label l=0 ; l <= dim[2] ; l++)
            {
                if( (n+m+l) <= polOrder_ && (n+m+l) > 0 )
                {
                    A[cellj-1][currIdx++] = volIntegralsIJ[n][m][l];
                }
            }
        }
    }
}

// ---------------------------- Constructor ------------------------------------

Foam::WENOBase::WENOBase
(
    const fvMesh& mesh,
    const label polOrder
)
{
    /**************************** General Note ********************************\
    Collecting Stencils:
      At first all neighbouring cells are collected within one list which is 
      always saved in stencilID_[cellI][0]. After all cells have been collected
      and potentially been corrected in case of a prallel run the cells are
      distributed to each sector, where a sector is one face of the cell.
      Therefore each cell has n number of stencils where n is the number of
      faces of the cell plus 1 for the collection of all the cells before
      splitting.
      
      This is also why until the splitStencil() function only the 
      stencilID_[cellI][0] list is manipulated and no other list. 
    \**************************************************************************/
    
    const WENO::globalfvMesh globalfvMesh(mesh);

    // Note the local mesh is the mesh of the processor, the global mesh is the
    // reconstructed mesh from all processors 
    const fvMesh& localMesh = globalfvMesh.localMesh();
    const fvMesh& globalMesh = globalfvMesh();

    polOrder_ = polOrder;

    Dir_ = localMesh.time().path()/"constant"/"WENOBase" + Foam::name(polOrder_);

    // Calculate the degrees of freedom and sets the dimensions 
    setDegreeOfFreedom(localMesh);

    // Create new lists if necessary
    if (!readList(localMesh))
    {
        // Read expert factor
        IOdictionary WENODict
        (
            IOobject
            (
                "WENODict",
                localMesh.time().caseSystem(),
                localMesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );

        const scalar extendRatio =
            WENODict.lookupOrDefault<scalar>("extendRatio", 2.5);

        // ------------- Initialize Lists --------------------------------------

        stencilsID_.setSize(localMesh.nCells());
        cellToProcMap_.setSize(localMesh.nCells());

        labelList nStencils(localMesh.nCells(),0);

        haloCenters_.setSize(Pstream::nProcs());
        
        ownHalos_.setSize(Pstream::nProcs());

        List<List<List<point> > > haloTriFaceCoord(Pstream::nProcs());

        // ------------------ Start Processing ---------------------------------
        
        // Initialize the volume integrals 
        volIntegralType volIntegrals;   // Dummy variable for volumeIntegral of one cell
        initVolIntegrals(globalfvMesh,volIntegrals);

        Info << "\t1) Create local stencils..." << endl;
        createStencilID(globalMesh,globalfvMesh.localToGlobalCellID(),nStencils,extendRatio);

        
        // Extension to halo cells, if neccessary
        if(Pstream::parRun())
        {
            Info << "\t2) Create haloCells ... " << endl;
            correctParallelRun(globalfvMesh,haloTriFaceCoord,nStencils,extendRatio);
        }

        Info << "\t3) Split stencil ... " << endl;
        // Split the stencil in several sectorial stencils
        const labelList& localToGlobalCellID = globalfvMesh.localToGlobalCellID();
        for
        (
            label localCellI = 0, globalCellI = localToGlobalCellID[localCellI];
            localCellI < localToGlobalCellID.size();
            localCellI++, globalCellI=localToGlobalCellID[localCellI < localToGlobalCellID.size() ? localCellI : 0]
        )
        {
            splitStencil(globalMesh, localMesh, localCellI, globalCellI, nStencils[localCellI]);
        }

        Info << "\t4) Calculate LS matrix ..." << endl;
        // Get the least squares matrices and their pseudoinverses
        LSmatrix_.setSize(localMesh.nCells());

        for (label cellI = 0; cellI < localMesh.nCells(); cellI++)
        {
            label excludeFace = 0;

            LSmatrix_[cellI].setSize(nStencils[cellI]);

            forAll(stencilsID_[cellI], stencilI)
            {
                if (stencilsID_[cellI][stencilI][0] != int(Cell::deleted))
                {
                    LSmatrix_[cellI][stencilI - excludeFace] =
                        calcMatrix
                        (
                            globalMesh,
                            localMesh,
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

        Info << "\t5) Calcualte smoothness indicator B..."<<endl;
        // Get the smoothness indicator matrices
        B_.setSize(localMesh.nCells());

        for(label cellI = 0; cellI < localMesh.nCells(); cellI++)
        {
            B_[cellI] =
                Foam::geometryWENO::getB
                (
                    localMesh,
                    cellI,
                    polOrder_,
                    nDvt_,
                    JInv_[cellI],
                    refPoint_[cellI],
                    dimList_[cellI]
                );

        }

        // Get surface integrals over basis functions in transformed coordinates

        intBasTrans_.setSize(localMesh.nFaces());
        
        refFacAr_.setSize(localMesh.nFaces());

        for (label faceI = 0; faceI < localMesh.nFaces(); faceI++)
        {
            intBasTrans_[faceI][0] = volIntegrals;
            intBasTrans_[faceI][1] = volIntegrals;
        }

        Foam::geometryWENO::surfIntTrans
        (
            localMesh,
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
            localMesh
        );
    }
    
    
    #ifdef FULLDEBUG
        volScalarField excludedStencils
        (
          IOobject
          (
           "excludeStencil",
           localMesh.time().timeName(),
           localMesh,
           IOobject::NO_READ,
           IOobject::NO_WRITE
          ),
          localMesh,
          dimensioned<scalar>("alphaSu", dimless, 0)
        );

        forAll(stencilsID_,cellI)
        {
            forAll(stencilsID_[cellI],stencilI)
            {
                if (stencilsID_[cellI][stencilI][0] == int(Cell::deleted))
                    excludedStencils[cellI] = excludedStencils[cellI]+1;
            }

        }

        excludedStencils.write();
    #endif


    // Clear all unwanted fields:
    volIntegralsList_.clear();

    JInv_.clear();

    refDet_.clear();

    refPoint_.clear();
}


void Foam::WENOBase::createStencilID
(
    const fvMesh& globalMesh,         // here the global mesh
    const labelList& cellID,
    labelList& nStencils,
    const scalar extendRatio
    
)
{
    for
    (
        label cellI = 0, globalCellI = cellID[cellI];
        cellI < cellID.size();
        cellI++,globalCellI=cellID[cellI < cellID.size() ? cellI : 0]
    )
    {
        // Note: local variables as nStencils or stencilID_ are accessed with 
        //       cellI. Global mesh values are accessed with globalCellI
        //       At first the stencil is populated with the globalCellI and is 
        //       later corrected
        const cell& faces = globalMesh.cells()[globalCellI];

        nStencils[cellI] = 1;

        forAll(faces, faceI)
        {
            if (faces[faceI] < globalMesh.nInternalFaces())
            {
                nStencils[cellI]++;
            }
        }

        stencilsID_[cellI].setSize(nStencils[cellI]);
        cellToProcMap_[cellI].setSize(nStencils[cellI]);

        forAll(stencilsID_[cellI],stencilI)
        {
            stencilsID_[cellI][stencilI].append(globalCellI);
        }
        stencilsID_[cellI][0].append(globalMesh.cellCells()[globalCellI]);

        
        // Store all cells within a map 
        // Neighbour cells are marked as true
        std::map<int,bool> stencilCells;
        
        stencilCells.insert(std::pair<int,bool>(stencilsID_[cellI][0][0],false));
        forAll(stencilsID_[cellI][0],i)
        {
            stencilCells.insert(std::pair<int,bool>(stencilsID_[cellI][0][i],true));
        }
        


        // Extend central stencil to neccessary size
        label minStencilSize = 0;
        // Maximum number of iterations for extendRatio
        const label maxIter = 100;
        label iter = 0;
        while (minStencilSize < 1.2*extendRatio*nDvt_*nStencils[cellI])
        {
            extendStencils
            (
                globalMesh,
                cellI,
                stencilCells,
                minStencilSize
            );
            iter++;
            if (iter > maxIter)
            {
                Pout << "ExtendStencil failed to reach criteria " 
                     << minStencilSize << " < " << 1.2*extendRatio*nDvt_*nStencils[cellI]
                     << "  for cell: " << cellI << nl
                     << "Maximum iteration reached. Continue with this stencil size..."<<endl;
                break;
            }
        }

        // Sort and cut stencil
        labelList dummyLabels(stencilsID_[cellI][0].size(),static_cast<int>(Cell::local));
        cellToProcMap_[cellI][0] = dummyLabels;

        sortStencil(globalMesh,cellI, extendRatio*nDvt_*nStencils[cellI]);
    }
}


void Foam::WENOBase::correctParallelRun
(
    const WENO::globalfvMesh& globalfvMesh,
    List<List<List<point> > >& haloTriFaceCoord,
    const labelList& nStencils,
    const scalar extendRatio
)
{
    const fvMesh& localMesh = globalfvMesh.localMesh();
    
    
    // labelList of halo cells as their cellID in the local mesh
    // of their processor domain
    labelListList haloProcessorCellID(Pstream::nProcs());
    
    // labelList of halo cells as their cellID in the global mesh
    labelListList haloGlobalCellID(Pstream::nProcs());

    // List of new halo cell ID starting at zero
    List<std::map<int,int> > stencilNewHaloID(Pstream::nProcs());
    
    labelList haloCellsPerProcessor(Pstream::nProcs(),0);

    // Default is -1 and if a processor is needed it is set to the processorID
    receiveProcList_.setSize(Pstream::nProcs());
    forAll(receiveProcList_,procI)
    {
        receiveProcList_[procI] = -1;
    }

    // Checking for FULLDEBUG mode if the cell centers match before and 
    // after correcting the stencilID list
    #ifdef FULLDEBUG
        List<List<vector > > centers;
        centers.setSize(stencilsID_.size());
        forAll(centers,cellI)
        {
            centers[cellI].setSize(stencilsID_[cellI][0].size());
            forAll(centers[cellI],j)
            {
                centers[cellI][j] = globalfvMesh().C()[stencilsID_[cellI][0][j]];
            }
        }
    #endif



    
    // Loop over all stencil and check if the cells are local or halo
    forAll(stencilsID_,cellI)
    {
        forAll(stencilsID_[cellI][0],i)
        {
            if (!globalfvMesh.isLocalCell(stencilsID_[cellI][0][i]))
            {
                int procID = globalfvMesh.getProcID(stencilsID_[cellI][0][i]);
                
                cellToProcMap_[cellI][0][i] = procID;
                
                receiveProcList_[procID] = procID;

                // If the cell has not been added jet add it to the halo Cell 
                auto it = stencilNewHaloID[procID].find(stencilsID_[cellI][0][i]);
                
                if (it == stencilNewHaloID[procID].end())
                {
                    haloProcessorCellID[procID].append
                    (
                        globalfvMesh.processorCellID(stencilsID_[cellI][0][i])
                    );
                    
                    
                    haloGlobalCellID[procID].append
                    (
                        stencilsID_[cellI][0][i]
                    );
                    
                    haloCenters_[procID].append
                    (
                        globalfvMesh().C()[stencilsID_[cellI][0][i]]
                    );
                    
                    // Create entry in map
                    stencilNewHaloID[procID].insert
                    (
                        std::pair<int,int>
                        (
                            stencilsID_[cellI][0][i],
                            haloCellsPerProcessor[procID]++
                        )
                    );
                    stencilsID_[cellI][0][i] 
                        = stencilNewHaloID[procID][stencilsID_[cellI][0][i]];
                }
                else
                {
                    stencilsID_[cellI][0][i] 
                        = stencilNewHaloID[procID][stencilsID_[cellI][0][i]];
                }
                     
                
            }
            else
            {
                // If cell is a local cell the stencilID has to be changed to a 
                // local cellID 
                stencilsID_[cellI][0][i] = 
                    globalfvMesh.processorCellID(stencilsID_[cellI][0][i]);
                    
                cellToProcMap_[cellI][0][i] = int(Cell::local);
            }
            
        }
    }

    // Calculate and store TriFace coordinates for volume integral of global mesh
    forAll(haloTriFaceCoord, procI)
    {
        haloTriFaceCoord[procI].setSize(haloGlobalCellID[procI].size());

        forAll(haloGlobalCellID[procI], cellI)
        {
            haloTriFaceCoord[procI][cellI] =
                Foam::geometryWENO::getTriFaces
                (
                    globalfvMesh(),
                    haloGlobalCellID[procI][cellI]
                );
        }
    }


    // Distribute halo cells
    distributeStencils
    (
        haloProcessorCellID
    );

    // Store the local cellID of your own halos
    ownHalos_ = haloProcessorCellID;
    
    #ifdef FULLDEBUG
        // Check centers:
        forAll(centers,cellI)
        {
            forAll(centers[cellI],j)
            {
                if (cellToProcMap_[cellI][0][j] == int(Cell::local))
                {
                    if (mag(centers[cellI][j] - localMesh.C()[stencilsID_[cellI][0][j]])>1E-15)
                        FatalError << "Local cell Center does not match "<<nl
                                   << "local: "<<localMesh.C()[stencilsID_[cellI][0][j]] 
                                   << "  global: "<<centers[cellI][j]<< exit(FatalError);
                }
                else
                {
                    
                    if (centers[cellI][j] != haloCenters_[cellToProcMap_[cellI][0][j]][stencilsID_[cellI][0][j]])
                        FatalError << "Halo cell Center does not match " << nl
                                   << "local: "<<haloCenters_[cellToProcMap_[cellI][0][j]][stencilsID_[cellI][0][j]]
                                   << "  global: "<<centers[cellI][j]<< exit(FatalError);
                }
            }
        }
    #endif
    
    // Get final big central stencils
    for (label cellI = 0; cellI < localMesh.nCells(); cellI++)
    {
        sortStencil(localMesh, cellI, (extendRatio*nDvt_)*nStencils[cellI]);
    }
}



void Foam::WENOBase::setDegreeOfFreedom(const fvMesh& localMesh)
{
    // 3D version
    if (localMesh.nSolutionD() == 3)
    {
        nDvt_ = (polOrder_ + 1.0)*(polOrder_ + 2.0)*(polOrder_ + 3.0)/6.0 - 1.0;
    }
    else // 2D version (always only one cell in z-direction)
    {
        nDvt_ = (polOrder_ + 1.0)*(polOrder_ + 2.0)/2.0 - 1.0;
    }

    // Set the dimList
    labelList dummyList(3,0);
    dimList_.setSize(localMesh.nCells(),dummyList);
    
    // Vector with valid dimensions
    vector dimMesh = localMesh.solutionD();
    
    for (label i = 0; i < localMesh.nCells(); i++)
    {
        dimList_[i][0] = (dimMesh[0] == 1 ? polOrder_ : 0);
        dimList_[i][1] = (dimMesh[1] == 1 ? polOrder_ : 0);
        dimList_[i][2] = (dimMesh[2] == 1 ? polOrder_ : 0);
    }
}


void Foam::WENOBase::initVolIntegrals
(
    const WENO::globalfvMesh& globalfvMesh,
    volIntegralType& volIntegrals
)
{
    // local mesh
    const fvMesh& localMesh = globalfvMesh.localMesh();
    const fvMesh& globalMesh = globalfvMesh();

    volIntegrals.resize((polOrder_ + 1));

    for (label i = 0; i < (polOrder_+1); i++)
    {
        volIntegrals[i].resize((polOrder_+ 1)-i);

        for (label j = 0; j < ((polOrder_+1)-i); j++)
        {
            volIntegrals[i][j].resize((polOrder_ + 1)-i, 0.0);
        }
    }

    volIntegralsList_.setSize(localMesh.nCells(), volIntegrals);

    JInv_.setSize(localMesh.nCells());
    refPoint_.setSize(localMesh.nCells());
    refDet_.setSize(localMesh.nCells());
    
    for
    (
        label cellI = 0, globalCellI = globalfvMesh.localToGlobalCellID()[cellI];
        cellI < globalfvMesh.localToGlobalCellID().size();
        cellI++, globalCellI=globalfvMesh.localToGlobalCellID()[cellI < globalfvMesh.localToGlobalCellID().size() ? cellI : 0]
    )
    {
        // Create the volume integral of each cell
        Foam::geometryWENO::initIntegrals
        (
            globalMesh,
            globalCellI,
            polOrder_,
            volIntegralsList_[cellI],
            JInv_[cellI],
            refPoint_[cellI],
            refDet_[cellI]
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

        IFstream isDL(Dir_/"DimLists",IFstream::streamFormat::BINARY);
        dimList_.setSize(mesh.nCells());

        forAll(dimList_, cellI)
        {
            isDL >> dimList_[cellI];
        }

        IFstream isSID(Dir_/"StencilIDs",IFstream::streamFormat::BINARY);
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

        IFstream isCToP(Dir_/"CellToProcMap",IFstream::streamFormat::BINARY);
        cellToProcMap_.setSize(mesh.nCells());

        for (label cellI = 0; cellI < mesh.nCells(); cellI++)
        {
            isCToP >> nEntries;

            cellToProcMap_[cellI].setSize(nEntries);

            for (label stencilI = 0; stencilI < nEntries; stencilI++)
            {
                isCToP >> cellToProcMap_[cellI][stencilI];
            }
        }

        IFstream isLS(Dir_/"Pseudoinverses",IFstream::streamFormat::BINARY);
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

        IFstream isB(Dir_/"B",IFstream::streamFormat::BINARY);
        B_.setSize(mesh.nCells());

        forAll(B_, cellI)
        {
            isB >> B_[cellI];
        }


        sendProcList_.setSize(Pstream::nProcs());
        IFstream isPToPSend(Dir_/"sendProcList",IFstream::streamFormat::BINARY);

        forAll(sendProcList_, procI)
        {
            isPToPSend >> sendProcList_[procI];
        }

        receiveProcList_.setSize(Pstream::nProcs());
        IFstream isPToPReceive(Dir_/"receiveProcList",IFstream::streamFormat::BINARY);

        forAll(receiveProcList_, procI)
        {
            isPToPReceive >> receiveProcList_[procI];
        }

        ownHalos_.setSize(Pstream::nProcs());
        IFstream isOH(Dir_/"OwnHalos",IFstream::streamFormat::BINARY);

        forAll(ownHalos_, procI)
        {
            isOH >> nEntries;

            ownHalos_[procI].setSize(nEntries);

            forAll(ownHalos_[procI], cellI)
            {
                isOH >> ownHalos_[procI][cellI];
            }
        }

        haloCenters_.setSize(Pstream::nProcs());
        IFstream isHalo(Dir_/"HaloCenters",IFstream::streamFormat::BINARY);

        forAll(haloCenters_, procI)
        {
            isHalo >> nEntries;

            haloCenters_[procI].setSize(nEntries);

            forAll(haloCenters_[procI], cellI)
            {
                isHalo >> haloCenters_[procI][cellI];
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

        refFacAr_.setSize(mesh.nFaces());

        for (label faceI = 0; faceI < mesh.nFaces(); faceI++)
        {
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
        return false;
    }
}


void Foam::WENOBase::writeList
(
    const fvMesh& mesh
)
{
    Info<< "Write created lists to constant folder \n" << endl;

    mkDir(Dir_);

    OFstream osPToPSend(Dir_/"sendProcList",OFstream::streamFormat::BINARY);

    forAll(sendProcList_, i)
    {
        osPToPSend << sendProcList_[i] << endl;
    }

    OFstream osPToPReceive(Dir_/"receiveProcList",OFstream::streamFormat::BINARY);

    forAll(receiveProcList_, i)
    {
        osPToPReceive << receiveProcList_[i] << endl;
    }

    OFstream osDL(Dir_/"DimLists",OFstream::streamFormat::BINARY);

    forAll(dimList_, cellI)
    {
        osDL<< dimList_[cellI] << endl;
    }

    OFstream osSID(Dir_/"StencilIDs",OFstream::streamFormat::BINARY);

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        osSID<< stencilsID_[cellI].size() << endl;

        forAll(stencilsID_[cellI], stenciI)
        {
            osSID<< stencilsID_[cellI][stenciI] << endl;
        }
    }

    OFstream osCToP(Dir_/"CellToProcMap",OFstream::streamFormat::BINARY);

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        osCToP<< cellToProcMap_[cellI].size() << endl;

        forAll(cellToProcMap_[cellI], cellJ)
        {
            osCToP<< cellToProcMap_[cellI][cellJ] << endl;
        }
    }

    OFstream osLS(Dir_/"Pseudoinverses",OFstream::streamFormat::BINARY);
    osLS.precision(10);

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        osLS<< LSmatrix_[cellI].size() << endl;

        forAll(LSmatrix_[cellI], stenciI)
        {
            osLS<< LSmatrix_[cellI][stenciI] << endl;
        }
    }

    OFstream osHalo(Dir_/"HaloCenters",OFstream::streamFormat::BINARY);
    osHalo.precision(10);

    forAll(haloCenters_, patchI)
    {
        osHalo<< haloCenters_[patchI].size() << endl;

        forAll(haloCenters_[patchI], cellI)
        {
            osHalo<< haloCenters_[patchI][cellI] << endl;
        }
    }

    OFstream osOH(Dir_/"OwnHalos",OFstream::streamFormat::BINARY);
    osOH.precision(10);

    
    forAll(ownHalos_, procI)
    {
        osOH<< ownHalos_[procI].size() << endl;

        forAll(ownHalos_[procI], cellI)
        {
            osOH<< ownHalos_[procI][cellI] << endl;
        }
    }

    OFstream osB(Dir_/"B",OFstream::streamFormat::BINARY);
    osB.precision(10);

    forAll(B_, cellI)
    {
        osB<< B_[cellI] << endl;
    }
}


//void Foam::WENOBase::printStatistics()
//{



    //// Get the information of all stencilIDs of all processors 
    //// Check if parallel 
    //if (Pstream::parRun())
    //{
        //Info << "Processor Statistics:"<<endl;
        
        //// Copy the field to avoid overwriting it
        //List<labelListList> stencilIDList(stencilsID_);
        
        //List<List<labelListList> > allValues(Pstream::nProcs());

        //allValues[Pstream::myProcNo()] = stencilIDList;

        //Pstream::gatherList(allValues);

        //if (Pstream::master())
        //{
            //// Loop over all processor entries and check how many stencils
            //// have been excluded
            //forAll(allValues,i)
            //{
                //List<labelListList>& list = allValues[i];
                
                //labelList excludedStencils(list.size(),0);

                //forAll(list,cellI)
                //{
                    //forAll(list[cellI],stencilI)
                    //{
                        //if (list[cellI][stencilI][0] == int(Cell::deleted))
                        //{
                            //excludedStencils[cellI]++;
                        //}
                    //}
                //}
                //Info << "\t" << "Processor"<<i<<":"<<nl
                     //<< "\t" << "\tExcludedStencils: "<< sum(excludedStencils) << endl;

            //}
        //}
        //Info << endl;
    //}
    
    
//}

// ************************************************************************* //
