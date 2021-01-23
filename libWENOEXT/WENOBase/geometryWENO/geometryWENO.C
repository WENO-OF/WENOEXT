/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "geometryWENO.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


void Foam::geometryWENO::initIntegrals
(
    const fvMesh& mesh,
    const label cellI,
    const label polOrder,
    volIntegralType& volIntegrals,
    scalarSquareMatrix& JInvI,
    point& refPointI,
    scalar& refDetI
)
{
    const pointField& pts = mesh.points();
    const faceList& fcs = mesh.faces();
    const cell & cc = mesh.cells()[cellI];

    const labelList pLabels(cc.labels(fcs));
    const labelList pEdge = mesh.pointPoints()[pLabels[0]];
    
    const scalar cellVolume = mesh.V()[cellI];
    

    // Create reference frame of new space

    labelList referenceFrame(1,pLabels[0]);

    forAll(pEdge, i)
    {
        forAll(pLabels, j)
        {
            if (pEdge[i] == pLabels[j])
            {
                referenceFrame.append(pEdge[i]);
            }
        }
    }

    refPointI = pts[referenceFrame[0]];

    label k = 1;

    // Check the quality of the chosen frame and change if necessary
    while
    (
        referenceFrame.size() < 4 
     || mag(det(jacobi(pts,referenceFrame)))/cellVolume < 1E-10 // cell normalized determinante 
    )
    {
        if (k >= pLabels.size())
            WarningIn("geometryWENO::initIntegrals calculate reference frame") 
                << "Determinante of Jacobian matrix smaller than 1E-10 ("
                <<det(jacobi(pts,referenceFrame))<<")"<<endl;
            
        const labelList pEdgeMod = mesh.pointPoints()[pLabels[k]];

        labelList modrefFrame(1, pLabels[k]);

        forAll(pEdgeMod, i)
        {
            forAll(pLabels, j)
            {
                if (pEdgeMod[i] == pLabels[j])
                {
                    modrefFrame.append(pEdgeMod[i]);
                }
            }
        }

        refPointI = pts[modrefFrame[0]];

        k++;


        // For some meshs it is possible that a point connects only two edges 
        // and thus does not span a tetrahedar
        if (modrefFrame.size() > 3)
            referenceFrame = modrefFrame;
    }

    scalarSquareMatrix J = jacobi(pts,referenceFrame);

    // We have to live with a copy assignment as Foam::Matrix() does not have
    // a move assignment operator... 
    JInvI = JacobiInverse(J);

    refDetI = det(JInvI);

    const point refPointTrans =
        Foam::geometryWENO::transformPoint
        (
            JInvI,
            mesh.cellCentres()[cellI],
            refPointI
        );

    // Triangulate the faces of the cell
    List<tetIndices> cellTets =
        polyMeshTetDecomposition::cellTetIndices(mesh, cellI);

    triFaceList triFaces(cellTets.size());

    forAll(cellTets, cTI)
    {
        triFaces[cTI] = cellTets[cTI].faceTriIs(mesh);
    }

    // Evaluate volume integral using surface integrals over triangulated faces
    // Initialize size of integral list
    volIntegrals.resize((polOrder + 1),(polOrder + 1),(polOrder + 1));
    // Set value to zero
    volIntegrals.setZero();

    forAll(triFaces, i)
    {
        const triFace& tri(triFaces[i]);

        vector v0 = transformPoint(JInvI, pts[tri[0]], refPointI);
        vector v1 = transformPoint(JInvI, pts[tri[1]], refPointI);
        vector v2 = transformPoint(JInvI, pts[tri[2]], refPointI);

        vector vn = (v1 - v0) ^ (v2 - v0);

        scalar area = 0.5*mag(vn);

        if (sign(vn & (v0 - refPointTrans)) < 0.0)
        {
             vn *= -1.0/mag(vn);
        }
        else
        {
            vn /= mag(vn);
        }

        // Evaluate integral using Gaussian quadratures
        for (label n = 0; n <= polOrder; n++)
        {
            for (label m = 0; m <= polOrder; m++)
            {
                for (label l = 0; l <= polOrder; l++)
                {
                    if ((n + m + l) <= polOrder && n > 0)
                    {
                        volIntegrals(n,m,l) +=
                            1.0/(n + 1)*area*vn.x()
                            *gaussQuad(n + 1, m, l, refPointTrans, v0, v1, v2);
                    }
                    else if ((n + m + l) <= polOrder && m > 0)
                    {
                        volIntegrals(n,m,l) +=
                            1.0/(m + 1)*area*vn.y()
                            *gaussQuad(n, m + 1, l, refPointTrans, v0, v1, v2);
                    }
                    else if ((n + m + l) <= polOrder)
                    {
                        volIntegrals(n,m,l) += 
                            1.0/(l + 1)*area*vn.z()
                            *gaussQuad(n, m, l + 1, refPointTrans, v0, v1, v2);
                    }
                }
            }
        }
    }

    for (label n = 0; n <= polOrder; n++)
    {
        for (label m = 0; m <= polOrder; m++)
        {
            for (label l = 0; l <= polOrder; l++)
            {
                if ((n + m + l) <= polOrder)
                {
                    volIntegrals(n,m,l) *=
                        1.0/(mag(refDetI)*mesh.cellVolumes()[cellI]);
                }
            }
        }
    }
}


void Foam::geometryWENO::transformIntegral
(
    const fvMesh& mesh,
    const label cellJ,
    const point transCenterJ,
    const label polOrder,
    const scalarSquareMatrix& JInvI,
    const point refPointI,
    const scalar refDetI,
    volIntegralType& transVolMom
)
{
    const pointField& pts = mesh.points();


    transVolMom.resize((polOrder + 1),(polOrder + 1),(polOrder + 1));

    // Initialize with zero
    transVolMom.setZero();

    // Triangulate the faces of the cell
    List<tetIndices> cellTets =
        polyMeshTetDecomposition::cellTetIndices(mesh, cellJ);

    triFaceList triFaces(cellTets.size());

    forAll(cellTets, cTI)
    {
        triFaces[cTI] = cellTets[cTI].faceTriIs(mesh);
    }

    // Evaluate volume integral using surface integrals over triangulated faces

    forAll(triFaces, i)
    {
        const triFace& tri(triFaces[i]);

        vector v0 = transformPoint(JInvI, pts[tri[0]], refPointI);
        vector v1 = transformPoint(JInvI, pts[tri[1]], refPointI);
        vector v2 = transformPoint(JInvI, pts[tri[2]], refPointI);

        vector vn = (v1 - v0) ^ (v2 - v0);

        scalar area = 0.5*mag(vn);

        if (sign(vn & (v0 - transCenterJ)) < 0.0)
        {
             vn *= -1.0/mag(vn);
        }
        else
        {
            vn /= mag(vn);
        }

        // Evaluate integral using Gaussian quadratures
        for (label n = 0; n <= polOrder; n++)
        {
            for (label m = 0; m <= polOrder; m++)
            {
                for (label l = 0; l <= polOrder; l++)
                {
                    if ((n + m + l) <= polOrder && n > 0)
                    {
                        transVolMom(n,m,l) +=
                            1.0/(n + 1)*area*vn.x()
                           *gaussQuad(n + 1, m, l, transCenterJ, v0, v1, v2);
                    }
                    else if ((n + m + l) <= polOrder && m > 0)
                    {
                        transVolMom(n,m,l) +=
                            1.0/(m + 1)*area*vn.y()
                           *gaussQuad(n, m + 1, l, transCenterJ, v0, v1, v2);
                    }
                    else if ((n + m + l) <= polOrder)
                    {
                        transVolMom(n,m,l) +=
                            1.0/(l + 1)*area*vn.z()
                            * gaussQuad(n, m, l + 1, transCenterJ, v0, v1, v2);
                    }
                }
            }
        }
    }

    for (label n = 0; n <= polOrder; n++)
    {
        for (label m = 0; m <= polOrder; m++)
        {
            for (label l = 0; l <= polOrder; l++)
            {
                if ((n + m + l) <= polOrder)
                {
                    transVolMom(n,m,l) *=
                        1.0/(mag(refDetI)*mesh.cellVolumes()[cellJ]);
                }
            }
        }
    }
}


Foam::scalarSquareMatrix Foam::geometryWENO::JacobiInverse
(
    const scalarSquareMatrix& J
)
{
    scalarSquareMatrix JacobiInv(3,0.0);

    scalar det = Foam::det(J);

    JacobiInv[0][0] = (J(1,1)*J(2,2)-J(1,2)*J(2,1))/det;

    JacobiInv[0][1] = (J(0,2)*J(2,1)-J(0,1)*J(2,2))/det;

    JacobiInv[0][2] = (J(0,1)*J(1,2)-J(0,2)*J(1,1))/det;

    JacobiInv[1][0] = (J(1,2)*J(2,0)-J(1,0)*J(2,2))/det;

    JacobiInv[1][1] = (J(0,0)*J(2,2)-J(0,2)*J(2,0))/det;

    JacobiInv[1][2] = (J(0,2)*J(1,0)-J(0,0)*J(1,2))/det;

    JacobiInv[2][0] = (J(1,0)*J(2,1)-J(1,1)*J(2,0))/det;

    JacobiInv[2][1] = (J(0,1)*J(2,0)-J(0,0)*J(2,1))/det;

    JacobiInv[2][2] = (J(0,0)*J(1,1)-J(0,1)*J(1,0))/det;

    return JacobiInv;
}


Foam::point Foam::geometryWENO::transformPoint
(
    const scalarSquareMatrix& Jinv,
    const point xP,
    const point x0
)
{
    point xiP(0,0,0);

    for (label q = 0; q < 3; q++)
    {
        for (label l = 0; l < 3; l++)
        {
            xiP[q] += Jinv[q][l]*(xP[l] - x0[l]);
        }
    }

    return xiP;
}


Foam::scalar Foam::geometryWENO::gaussQuad
(
        const scalar n,
        const scalar m,
        const scalar l,
        const point xi0,
        const vector v0,
        const vector v1,
        const vector v2
)
{
    // For integer power it is much faster to do an integer multiplication
    // This depends on the compiler used! For portability it is explicitly defined
    // here.
    auto intPow = [](const scalar base,const unsigned int exponent) -> scalar
    {
        scalar temp = 1.0;
        for (unsigned int i =0; i<exponent;i++)
        {
            temp *= base;
        }
        return temp;
    };

    // Sum up over Gaussian points with transformation on projected triangle

    scalar sum = 0.0;

    for (label j = 0; j < 13; j++)
    {
        scalar xi =
            v0.x()* (1 - gaussCoeff[j][0] - gaussCoeff[j][1])
          + v1.x()* gaussCoeff[j][0] + v2.x()*gaussCoeff[j][1] - xi0.x();
        scalar eta =
            v0.y()* (1- gaussCoeff[j][0]- gaussCoeff[j][1])
          + v1.y()* gaussCoeff[j][0] +v2.y()* gaussCoeff[j][1] - xi0.y();
        scalar zeta =
            v0.z()* (1- gaussCoeff[j][0]- gaussCoeff[j][1])
          + v1.z()* gaussCoeff[j][0] +v2.z()* gaussCoeff[j][1] - xi0.z();
    
        sum += gaussCoeff[j][2]* intPow(xi,n)*intPow(eta,m)*intPow(zeta,l);
    }

    return sum;
}



void Foam::geometryWENO::smoothIndIntegrals
(
    const fvMesh& mesh,
    const label cellI,
    const label polOrder,
    const scalarSquareMatrix& JInvI,
    const point refPointI,
    volIntegralType& smoothVolIntegral
)
{
    const pointField& pts = mesh.points();
    const label maxOrder = 2*polOrder - 2;

    smoothVolIntegral.resize((maxOrder + 1),(maxOrder + 1),(maxOrder + 1));
    smoothVolIntegral.setZero();


    List<tetIndices> cellTets =
        polyMeshTetDecomposition::cellTetIndices(mesh, cellI);

    point transCenterI =
        Foam::geometryWENO::transformPoint
        (
            JInvI,
            mesh.cellCentres()[cellI],
            refPointI
        );

    triFaceList triFaces(cellTets.size());

    forAll(cellTets, cTI)
    {
        triFaces[cTI] = cellTets[cTI].faceTriIs(mesh);
    }

    forAll(triFaces, i)
    {
        const triFace& tri(triFaces[i]);

        vector v0 = transformPoint(JInvI,pts[tri[0]],refPointI);
        vector v1 = transformPoint(JInvI,pts[tri[1]],refPointI);
        vector v2 = transformPoint(JInvI,pts[tri[2]],refPointI);

        vector vn = (v1 - v0) ^ (v2 - v0);

        scalar area = 0.5*mag(vn);

        if (sign(vn & (v0 - transCenterI)) < 0.0)
        {
             vn *= -1.0/mag(vn);
        }
        else
        {
            vn /= mag(vn);
        }

        for (label potXi = 0; potXi <= maxOrder; potXi++)
        {
            for (label potEta = 0; potEta <= maxOrder; potEta++)
            {
                for (label potZeta = 0; potZeta <= maxOrder; potZeta++)
                {
                    if ((potXi + potEta + potZeta) <= maxOrder && potXi > 0)
                    {
                        smoothVolIntegral(potXi,potEta,potZeta) +=

                            1.0/(potXi + 1)*area*vn.x()
                           *gaussQuad
                            (
                                potXi + 1,
                                potEta,
                                potZeta,
                                transCenterI,
                                v0,
                                v1,
                                v2
                            );
                    }
                    else if
                    (
                        (potXi + potEta + potZeta) <= maxOrder
                     && potEta > 0
                    )
                    {
                        smoothVolIntegral(potXi,potEta,potZeta) +=
                            1.0/(potEta + 1)*area*vn.y()
                           *gaussQuad
                            (
                                potXi,
                                potEta + 1,
                                potZeta,
                                transCenterI,
                                v0,
                                v1,
                                v2
                            );
                    }
                    else if ((potXi + potEta + potZeta) <= maxOrder)
                    {
                        smoothVolIntegral(potXi,potEta,potZeta) +=
                            1.0/(potZeta + 1)*area*vn.z()
                           *gaussQuad
                            (
                                potXi,
                                potEta,
                                potZeta + 1,
                                transCenterI,
                                v0,
                                v1,
                                v2
                            );
                    }
                }
            }
        }
    }
}


Foam::geometryWENO::DynamicMatrix Foam::geometryWENO::getB
(
    const fvMesh& mesh,
    const label cellI,
    const label polOrder,
    const label nDvt,
    const scalarSquareMatrix& JInvI,
    const point refPointI,
    const labelList& dim
)
{
    DynamicMatrix B
    (
        nDvt,
        nDvt,
        scalar(0.0)
    );

    // Get all necessary volume integrals
    volIntegralType intB;
    
    smoothIndIntegrals
    (
        mesh,
        cellI,
        polOrder,
        JInvI,
        refPointI,
        intB
    );

    label p = 0;
    label q = 0;

    for (label n1 = 0; n1 <= dim[0]; n1++)
    {
        for (label m1 = 0; m1 <= dim[1]; m1++)
        {
            for (label l1 = 0; l1 <= dim[2]; l1++)
            {
                if ((n1 + m1 + l1) <= polOrder && (n1 + m1 + l1) > 0)
                {
                    for (label n2 = 0; n2 <= dim[0]; n2++)
                    {
                        for (label m2 = 0; m2 <= dim[1]; m2++)
                        {
                            for (label l2 = 0; l2 <= dim[2]; l2++)
                            {
                                if ((n2+m2+l2) <= polOrder && (n2+m2+l2) > 0)
                                {
                                    for
                                    (
                                        label lambda = 1;
                                        lambda <= polOrder;
                                        lambda++
                                    )
                                    {
                                        for
                                        (
                                            label alpha = 0;
                                            alpha <= lambda;
                                            alpha++
                                        )
                                        {
                                            for
                                            (
                                                label beta = 0;
                                                beta <= (lambda - alpha);
                                                beta++
                                            )
                                            {
                                                label gamma =
                                                    lambda - alpha - beta;

                                                scalar K =
                                                    (
                                                        pos0(n1 - alpha)
                                                       *pos0(n2 - alpha)
                                                       *pos0(m1 - beta)
                                                       *pos0(m2 - beta)
                                                       *pos0(l1 - gamma)
                                                       *pos0(l2 - gamma)
                                                       *Fac(n1)*Fac(m1)*Fac(l1)
                                                       *Fac(n2)*Fac(m2)*Fac(l2)
                                                    )
                                                   /(
                                                        Fac(n1 - alpha)
                                                       *Fac(n2 - alpha)
                                                       *Fac(m1 - beta)
                                                       *Fac(m2 - beta)
                                                       *Fac(l1 - gamma)
                                                       *Fac(l2 - gamma)
                                                    );

                                                if (K != 0)
                                                {
                                                    B(p,q) +=
                                                       K*intB
                                                       (
                                                        n1 + n2 - 2*alpha,
                                                        m1 + m2 - 2*beta,
                                                        l1 + l2 - 2*gamma
                                                       );
                                                }
                                            }
                                        }
                                    }

                                    q++;
                                }
                            }
                        }
                    }

                    p++;
                    q = 0;
                }
            }
        }
    }

    return B;
}


void Foam::geometryWENO::surfIntTrans
(
    const fvMesh& mesh,
    const label polOrder,
    const List<volIntegralType>& volIntegralsList,
    const List<scalarSquareMatrix>& JInv,
    const List<point>& refPoint,
    List<Pair<volIntegralType>>& intBasTrans,
    List<scalar>& refFacAr
)
{
    const pointField& pts = mesh.points();
    const labelUList& N = mesh.neighbour();

    // Initialize fields
    forAll(intBasTrans,cellI)
    {
        intBasTrans[cellI][0].resize(polOrder+1,polOrder+1,polOrder+1);
        intBasTrans[cellI][0].setZero();
        intBasTrans[cellI][1].resize(polOrder+1,polOrder+1,polOrder+1);
        intBasTrans[cellI][1].setZero();
    }

    refFacAr.resize(mesh.nFaces(),0);

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        point refPointTrans =
            Foam::geometryWENO::transformPoint
            (
                JInv[cellI],
                mesh.cellCentres()[cellI],
                refPoint[cellI]
            );

        const cell& faces = mesh.cells()[cellI];

        for (label faceI = 0; faceI < faces.size(); faceI++)
        {
            // If face is neither in owner or neighbour it is at the boundary
            // and thus an owner 
            label OwnNeighIndex = 0;
            
            if (faces[faceI] < N.size() && cellI == N[faces[faceI]])
            {
                OwnNeighIndex = 1;
            }
            
            // Triangulate the faces
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

            scalar area = 0;

            // Evaluate surface integral using Gaussian quadratures
            forAll(triFaces, i)
            {
                const triFace& tri(triFaces[i]);

                vector v0 =
                    Foam::geometryWENO::transformPoint
                    (
                        JInv[cellI],
                        pts[tri[0]],
                        refPoint[cellI]
                    );
                vector v1 =
                    Foam::geometryWENO::transformPoint
                    (
                        JInv[cellI],
                        pts[tri[1]],
                        refPoint[cellI]
                    );
                vector v2 =
                    Foam::geometryWENO::transformPoint
                    (
                        JInv[cellI],
                        pts[tri[2]],
                        refPoint[cellI]
                    );

                vector vn = (v1 - v0) ^ (v2 - v0);

                area = 0.5*mag(vn);

                /**************************************************************\
                Note: The face is the same for the neighbour and the owner
                      Therefore integration is the same and looping over all
                      owner faces will include all faces.
                \**************************************************************/
                if (OwnNeighIndex == 0)                
                    refFacAr[faces[faceI]] += area;

                if (sign(vn & (v0 - refPointTrans)) < 0.0)
                {
                     vn *= -1.0/mag(vn);
                }
                else
                {
                    vn /= mag(vn);
                }

                for (label n = 0; n <= polOrder; n++)
                {
                    for (label m = 0; m <= polOrder; m++)
                    {
                        for (label l = 0; l <= polOrder; l++)
                        {
                            if ((n + m + l) <= polOrder)
                            {
                                intBasTrans[faces[faceI]][OwnNeighIndex](n,m,l) +=
                                    area
                                   *geometryWENO::gaussQuad
                                    (
                                        n,
                                        m,
                                        l,
                                        refPointTrans,
                                        v0,
                                        v1,
                                        v2
                                    );
                            }
                        }
                    }
                }
            }

            // Subtract volume integrals
            for (label n = 0; n <= polOrder; n++)
            {
                for (label m = 0; m <= polOrder; m++)
                {
                    for (label l = 0; l <= polOrder; l++)
                    {
                        if ((n + m + l) <= polOrder)
                        {
                            intBasTrans[faces[faceI]][OwnNeighIndex](n,m,l) -=
                            (
                                area*volIntegralsList[cellI](n,m,l)
                            );
                        }
                    }
                }
            }
        }
    }
}


Foam::scalar Foam::geometryWENO::Fac(label x)
{
    if (x <= 0)
    { 
        return 1;
    } 
    else
    {
        return Foam::factorial(x);
    }
}


Foam::vector Foam::geometryWENO::compCheck
(
    const label n,
    const label m,
    const label l,
    const volIntegralType& intBasisfI
)
{
    vector result(0.0,0.0,0.0);

    if (n > 0) result[0] = n*intBasisfI(n - 1,m,l);
    if (m > 0) result[1] = m*intBasisfI(n,m - 1,l);
    if (l > 0) result[2] = l*intBasisfI(n,m,l - 1);

    return result;
}


Foam::geometryWENO::scalarSquareMatrix Foam::geometryWENO::jacobi
(
    const pointField& pts,
    const labelList& referenceFrame
)
{
    scalarSquareMatrix J(3,0);
    
    for (int i = 0;i<3;i++)
    {
        for (int j = 0; j<3;j++)
        {
            J(i,j) = pts[referenceFrame[j+1]][i] - pts[referenceFrame[0]][i];
        }
    }
    
    return J;
}


Foam::geometryWENO::scalarSquareMatrix Foam::geometryWENO::jacobi
(
    const scalar x0, const scalar y0, const scalar z0,
    const scalar x1, const scalar y1, const scalar z1,
    const scalar x2, const scalar y2, const scalar z2,
    const scalar x3, const scalar y3, const scalar z3
)
{
    scalarSquareMatrix J(3,0);
    
    J(0,0) = x1-x0;
    J(0,1) = x2-x0;
    J(0,2) = x3-x0;
    
    J(1,0) = y1-y0;
    J(1,1) = y2-y0;
    J(1,2) = y3-y0;
    
    J(2,0) = z1-z0;
    J(2,1) = z2-z0;
    J(2,2) = z3-z0;
    
    return J;
}

// ************************************************************************* //
