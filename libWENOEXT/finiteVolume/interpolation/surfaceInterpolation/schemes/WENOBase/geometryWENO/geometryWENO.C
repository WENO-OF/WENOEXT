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
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "geometryWENO.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


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
    // Points and weights for the Gaussian quadrature of a standard triangle
    // - 1.row: x-values
    // - 2.row: y-values
    // - 3.row: weights
    scalar xw[7][3];

    xw[0][0] = 0.33333333333333;
    xw[0][1] = 0.33333333333333;
    xw[0][2] = 0.22500000000000;
    xw[1][0] = 0.47014206410511;
    xw[1][1] = 0.47014206410511;
    xw[1][2] = 0.13239415278851;
    xw[2][0] = 0.47014206410511;
    xw[2][1] = 0.05971587178977;
    xw[2][2] = 0.13239415278851;
    xw[3][0] = 0.05971587178977;
    xw[3][1] = 0.47014206410511;
    xw[3][2] = 0.13239415278851;
    xw[4][0] = 0.10128650732346;
    xw[4][1] = 0.10128650732346;
    xw[4][2] = 0.12593918054483;
    xw[5][0] = 0.10128650732346;
    xw[5][1] = 0.79742698535309;
    xw[5][2] = 0.12593918054483;
    xw[6][0] = 0.79742698535309;
    xw[6][1] = 0.10128650732346;
    xw[6][2] = 0.12593918054483;

    // Sum up over Gaussian points with transformation on projected triangle

    scalar sum = 0.0;

    for (label j = 0; j < 7; j++)
    {
        scalar xi =
            v0.x()*(1 - xw[j][0] - xw[j][1])
          + v1.x()*xw[j][0] + v2.x()*xw[j][1];
        scalar eta =
            v0.y()*(1 - xw[j][0] - xw[j][1])
          + v1.y()*xw[j][0] + v2.y()*xw[j][1];
        scalar zeta =
            v0.z()*(1 - xw[j][0] - xw[j][1])
          + v1.z()*xw[j][0] + v2.z()*xw[j][1];

        sum +=
            xw[j][2]*pow((eta - xi0.y()), m)
           *pow((zeta - xi0.z()), l)*pow((xi - xi0.x()), n);
    }

    return sum;
}


void Foam::geometryWENO::initIntegrals
(
    const fvMesh& mesh,
    const label cellI,
    const label polOrder,
    scalarMatrix& Integral,
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

    const scalar cellDimension = cc.mag(pts,fcs);



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
        checkRefFrame
        (
            jacobi(pts,referenceFrame),
            cellDimension
        )
        != true
        && (k < pLabels.size())
    )
    {
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

        // Initialize size of integral list
        Integral.resize((polOrder + 1));

        for (label i = 0; i <= polOrder; i++)
        {
            Integral[i].resize((polOrder+ 1));

            for (label j = 0; j <= polOrder; j++)
            {
                Integral[i][j].resize((polOrder + 1), 0.0);
            }
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
                        Integral[n][m][l] +=
                            1.0/(n + 1)*area*vn.x()
                            *gaussQuad(n + 1, m, l, refPointTrans, v0, v1, v2);
                    }
                    else if ((n + m + l) <= polOrder && m > 0)
                    {
                        Integral[n][m][l] +=
                            1.0/(m + 1)*area*vn.y()
                            *gaussQuad(n, m + 1, l, refPointTrans, v0, v1, v2);
                    }
                    else if ((n + m + l) <= polOrder)
                    {
                        Integral[n][m][l] +=
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
                    Integral[n][m][l] *=
                        1.0/(mag(refDetI)*mesh.cellVolumes()[cellI]);
                }
            }
        }
    }
}


Foam::List<Foam::point> Foam::geometryWENO::getTriFaces
(
    const fvMesh& mesh,
    const label cellI
)
{
    const pointField& pts = mesh.points();

    List<tetIndices> cellTets =
        polyMeshTetDecomposition::cellTetIndices(mesh, cellI);

    triFaceList triFaces(cellTets.size());

    forAll(cellTets, cTI)
    {
        triFaces[cTI] = cellTets[cTI].faceTriIs(mesh);
    }

    List<point> triFaceCoord(triFaces.size()*3, pTraits<vector>::zero);

    label k = 0;

    forAll(triFaces, i)
    {
        const triFace& tri(triFaces[i]);

        triFaceCoord[k] = pts[tri[0]];
        triFaceCoord[k+1] = pts[tri[1]];
        triFaceCoord[k+2] = pts[tri[2]];

        k = k + 3;
    }

    return triFaceCoord;
}


Foam::geometryWENO::scalarMatrix Foam::geometryWENO::getHaloMoments
(
    const fvMesh& mesh,
    const point transCenterJ,
    const List<point>& triFaceCoord,
    const label polOrder,
    const scalarSquareMatrix& JInvI,
    const point refPointI
)
{
    scalarMatrix Integral;

    Integral.resize((polOrder + 1));
    for (label i = 0; i < (polOrder + 1); i++)
    {
        Integral[i].resize((polOrder + 1));
        for(label j = 0; j < (polOrder + 1); j++)
        {
            Integral[i][j].resize((polOrder + 1), 0.0);
        }
    }

    label nTriFaces = triFaceCoord.size()/3.0;
    label k = 0;

    // Evaluate volume integral using surface integrals over triangulated faces

    for (label i = 0; i < nTriFaces; i++)
    {
        vector v0 = transformPoint(JInvI, triFaceCoord[k], refPointI);
        vector v1 = transformPoint(JInvI, triFaceCoord[k+1], refPointI);
        vector v2 = transformPoint(JInvI, triFaceCoord[k+2], refPointI);

        k = k + 3;

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
                        Integral[n][m][l] +=
                            1.0/(n + 1)*area*vn.x()
                           *gaussQuad(n + 1, m, l, transCenterJ, v0, v1, v2);
                    }
                    else if ((n + m + l) <= polOrder && m > 0)
                    {
                        Integral[n][m][l] +=
                            1.0/(m + 1)*area*vn.y()
                           *gaussQuad(n, m + 1, l, transCenterJ, v0, v1, v2);
                    }
                    else if ((n + m + l) <= polOrder)
                    {
                        Integral[n][m][l] +=
                            1.0/(l + 1)*area*vn.z()
                           *gaussQuad(n, m, l + 1, transCenterJ, v0, v1, v2);
                    }
                }
            }
        }
    }

    const scalar cellVolume = Integral[0][0][0];

    for (label n = 0; n <= polOrder; n++)
    {
        for (label m = 0; m <= polOrder; m++)
        {
            for (label l = 0; l <= polOrder; l++)
            {
                if ((n + m + l) <= polOrder)
                {
                    Integral[n][m][l] *= 1.0/(mag(cellVolume));
                }
            }
        }
    }

    return Integral;
}


Foam::geometryWENO::scalarMatrix Foam::geometryWENO::transformIntegral
(
    const fvMesh& mesh,
    const label cellJ,
    const point transCenterJ,
    const label polOrder,
    const scalarSquareMatrix& JInvI,
    const point refPointI,
    const scalar refDetI
)
{
    const pointField& pts = mesh.points();

    scalarMatrix Integral;

    Integral.resize((polOrder + 1));
    for (label i = 0; i <= polOrder; i++)
    {
        Integral[i].resize((polOrder + 1));
        for (label j = 0; j <= polOrder; j++)
        {
            Integral[i][j].resize((polOrder + 1), 0.0);
        }
    }

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
                        Integral[n][m][l] +=
                            1.0/(n + 1)*area*vn.x()
                           *gaussQuad(n + 1, m, l, transCenterJ, v0, v1, v2);
                    }
                    else if ((n + m + l) <= polOrder && m > 0)
                    {
                        Integral[n][m][l] +=
                            1.0/(m + 1)*area*vn.y()
                           *gaussQuad(n, m + 1, l, transCenterJ, v0, v1, v2);
                    }
                    else if ((n + m + l) <= polOrder)
                    {
                        Integral[n][m][l] +=
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
                    Integral[n][m][l] *=
                        1.0/(mag(refDetI)*mesh.cellVolumes()[cellJ]);
                }
            }
        }
    }

    return Integral;
}


bool Foam::geometryWENO::checkRefFrame
(
    const scalarSquareMatrix&& J,
    const scalar cellDimension
)
{
    // Calculate determinante of Jacobian matrix 
    // Determinante has to be greater than zero to calculate the inverse
    if
    (
        det(J) < cellDimension
    )
    {
        return false;
    }
    else
    {
        return true;
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


Foam::scalar Foam::geometryWENO::gaussQuadB
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
    // Points and weights for the Gaussian quadrature of a standard triangle
    // - 1.row: x-values
    // - 2.row: y-values
    // - 3.row: weights
    scalar xw[13][3];


    xw[0][0] = 0.33333333333333;
    xw[0][1] = 0.33333333333333;
    xw[0][2] = -0.14957004446768;
    xw[1][0] = 0.26034596607904;
    xw[1][1] = 0.26034596607904;
    xw[1][2] = 0.17561525743321;
    xw[2][0] = 0.26034596607904;
    xw[2][1] = 0.47930806784192;
    xw[2][2] = 0.17561525743321;
    xw[3][0] = 0.47930806784192;
    xw[3][1] = 0.26034596607904;
    xw[3][2] = 0.17561525743321;
    xw[4][0] = 0.06513010290222;
    xw[4][1] = 0.06513010290222;
    xw[4][2] = 0.05334723560884;
    xw[5][0] = 0.06513010290222;
    xw[5][1] = 0.86973979419557;
    xw[5][2] = 0.05334723560884;
    xw[6][0] = 0.86973979419557;
    xw[6][1] = 0.06513010290222;
    xw[6][2] = 0.05334723560884;
    xw[7][0] = 0.31286549600487;
    xw[7][1] = 0.63844418856981;
    xw[7][2] = 0.07711376089026;
    xw[8][0] = 0.63844418856981;
    xw[8][1] = 0.04869031542532;
    xw[8][2] = 0.07711376089026;
    xw[9][0] = 0.04869031542532;
    xw[9][1] = 0.31286549600487;
    xw[9][2] = 0.07711376089026;
    xw[10][0] = 0.63844418856981;
    xw[10][1] = 0.31286549600487;
    xw[10][2] = 0.07711376089026;
    xw[11][0] = 0.31286549600487;
    xw[11][1] = 0.04869031542532;
    xw[11][2] = 0.07711376089026;
    xw[12][0] = 0.04869031542532;
    xw[12][1] = 0.63844418856981;
    xw[12][2] = 0.07711376089026;


    // Sum up over Gaussian points with transformation on projected triangle

    scalar sum = 0.0;

    for (label j = 0; j < 13; j++)
    {
        scalar xi =
            v0.x()*(1 - xw[j][0] - xw[j][1])
          + v1.x()*xw[j][0] + v2.x()*xw[j][1];
        scalar eta =
            v0.y()* (1- xw[j][0]- xw[j][1])
          + v1.y()* xw[j][0] +v2.y()* xw[j][1] ;
        scalar zeta =
            v0.z()* (1- xw[j][0]- xw[j][1])
          + v1.z()* xw[j][0] +v2.z()* xw[j][1] ;

        sum +=
            xw[j][2]*pow(xi - xi0.x(), n)
           *pow(eta - xi0.y(), m)*pow(zeta - xi0.z(), l);
    }

    return sum;
}



Foam::geometryWENO::scalarMatrix Foam::geometryWENO::smoothIndIntegrals
(
    const fvMesh& mesh,
    const label cellI,
    const label polOrder,
    const scalarSquareMatrix& JInvI,
    const point refPointI
)
{
    const pointField& pts = mesh.points();
    const label maxOrder = 2*polOrder - 2;

    scalarMatrix Integral;
    Integral.setSize((maxOrder + 1));

    for (label i = 0; i < (maxOrder + 1); i++)
    {
        Integral[i].setSize((maxOrder + 1));

        for (label j = 0; j < (maxOrder + 1); j++)
        {
            Integral[i][j].setSize((maxOrder + 1), 0.0);
        }
    }

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
                        Integral[potXi][potEta][potZeta] +=
                            1.0/(potXi + 1)*area*vn.x()
                           *gaussQuadB
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
                        Integral[potXi][potEta][potZeta] +=
                            1.0/(potEta + 1)*area*vn.y()
                           *gaussQuadB
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
                        Integral[potXi][potEta][potZeta] +=
                            1.0/(potZeta + 1)*area*vn.z()
                           *gaussQuadB
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

    return Integral;
}


Foam::scalarRectangularMatrix Foam::geometryWENO::getB
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
    scalarRectangularMatrix B
    (
        nDvt,
        nDvt,
        scalar(0.0)
    );

    // Get all necessary volume integrals
    scalarMatrix intB =
        smoothIndIntegrals
        (
            mesh,
            cellI,
            polOrder,
            JInvI,
            refPointI
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
                                                        Pos(n1 - alpha)
                                                       *Pos(n2 - alpha)
                                                       *Pos(m1 - beta)
                                                       *Pos(m2 - beta)
                                                       *Pos(l1 - gamma)
                                                       *Pos(l2 - gamma)
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
                                                    B[p][q] +=
                                                       K*intB[n1 + n2 - 2*alpha]
                                                       [m1 + m2 - 2*beta]
                                                       [l1 + l2 - 2*gamma];
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
    const List<scalarMatrix>& volMom,
    const List<scalarSquareMatrix>& JInv,
    const List<point>& refPoint,
    List<List<scalarMatrix> >& intBasTrans,
    List<scalarList>& refFacAr
)
{
    const pointField& pts = mesh.points();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

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
            // initialize 
            label OwnNeighIndex = -1;

            if (faces[faceI] < P.size() && cellI == P[faces[faceI]])
            {
                OwnNeighIndex = 0;
            }
            else if (faces[faceI] < N.size() && cellI == N[faces[faceI]])
            {
                OwnNeighIndex = 1;
            }
            else
            {
                // If face is neither in owner or neighbour it is at the boundary
                // and thus an owner 
                OwnNeighIndex = 0;
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

                scalar area = 0.5*mag(vn);

                refFacAr[faces[faceI]][OwnNeighIndex] += area;

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
                                intBasTrans[faces[faceI]][OwnNeighIndex][n][m][l] +=
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
                            intBasTrans[faces[faceI]][OwnNeighIndex][n][m][l] -=
                            (
                                refFacAr[faces[faceI]][OwnNeighIndex]
                               *volMom[cellI][n][m][l]
                            );
                        }
                    }
                }
            }
        }
    }
}


Foam::scalar Foam::geometryWENO::Pos(scalar x)
{
    if (x >= 0) return 1.0;
    else
    {
        return 0.0;
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
    const scalarMatrix& intBasisfI
)
{
    vector result(0.0,0.0,0.0);

    if (n > 0) result[0] = n*intBasisfI[n - 1][m][l];
    if (m > 0) result[1] = m*intBasisfI[n][m - 1][l];
    if (l > 0) result[2] = l*intBasisfI[n][m][l - 1];

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
