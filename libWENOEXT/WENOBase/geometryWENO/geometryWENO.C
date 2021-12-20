/*---------------------------------------------------------------------------*\
       ██╗    ██╗███████╗███╗   ██╗ ██████╗     ███████╗██╗  ██╗████████╗
       ██║    ██║██╔════╝████╗  ██║██╔═══██╗    ██╔════╝╚██╗██╔╝╚══██╔══╝
       ██║ █╗ ██║█████╗  ██╔██╗ ██║██║   ██║    █████╗   ╚███╔╝    ██║   
       ██║███╗██║██╔══╝  ██║╚██╗██║██║   ██║    ██╔══╝   ██╔██╗    ██║   
       ╚███╔███╔╝███████╗██║ ╚████║╚██████╔╝    ███████╗██╔╝ ██╗   ██║   
        ╚══╝╚══╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝     ╚══════╝╚═╝  ╚═╝   ╚═╝   
-------------------------------------------------------------------------------                                                                                                                                                      
License
    This file is part of WENO Ext.

    WENO Ext is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    WENO Ext is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with  WENO Ext.  If not, see <http://www.gnu.org/licenses/>.

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "geometryWENO.H"
#include "mathFunctionsWENO.H"
#include "codeRules.H"

#if defined(__AVX__)
    #include <immintrin.h>
#endif

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
    

    labelList referenceFrame;
    
    // Loop over all points of the cell to find best conditioned reference system
    scalar bestCond = 1E+10;
    forAll(pLabels,k)
    {
        labelList modRefFrame(1,pLabels[k]);
        const labelList& pEdge = mesh.pointPoints(pLabels[k]);
        
        forAll(pEdge, i)
        {
            forAll(pLabels, j)
            {
                if (pEdge[i] == pLabels[j])
                {
                    modRefFrame.append(pEdge[i]);
                }
            }
        }
   
           
        if (modRefFrame.size() > 3 && cond(jacobi(pts,modRefFrame)) < bestCond)
        {
            bestCond = cond(jacobi(pts,modRefFrame));
            referenceFrame = modRefFrame;
        }
        else if (modRefFrame.size() > 3 && k == pLabels.size()-1 && bestCond == 1E+10)
        { 
            referenceFrame = modRefFrame;
            WarningInFunction
                << "Cannot calculate the reference frame with good condition ("
                << cond(jacobi(pts,modRefFrame)) << ") for cell: "<<cellI
                << " with coordinates "<<mesh.C()[cellI]<<endl;
            break;
        }
        else if (k == pLabels.size()-1 && bestCond == 1E+10)
        {
            FatalError << "Could not calculate reference frame in "
                       << "geometryWENO::initIntegrals() for point "<<mesh.C()[cellI]
                       << exit(FatalError);
        }
    }

    refPointI = pts[referenceFrame[0]];

    scalarSquareMatrix J = jacobi(pts,referenceFrame);
    
    calculateInverseJacobi(J,JInvI);

    // Use math functions WENO as blaze::det could be subject
    // to fatal floating point cancellation: 
    // see also: https://bitbucket.org/blaze-lib/blaze/issues/425/
    refDetI = mathFunctionsWENO::det(JInvI);
    if (mag(refDetI) < SMALL)
    {
        WarningInFunction
            << "Determinante of referenceFrame is too small ("
            << refDetI << ") for cell: "<<cellI
            << " with coordinates "<<mesh.C()[cellI]<<endl;
        refDetI = refDetI < 0 ? -1.0*SMALL : SMALL;
    }
        
    
    
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
        
        // For some mesh the triangle has a close to zero area. These triangles 
        // are excluded
        if (area < ROOTVSMALL)
            continue;

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
    const point& transCenterJ,
    const label polOrder,
    const scalarSquareMatrix& JInvI,
    const point& refPointI,
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
        
        // For some mesh the triangle has a close to zero area. These triangles 
        // are excluded
        if (area < ROOTVSMALL)
            continue;

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


Foam::point Foam::geometryWENO::transformPoint
(
    const scalarSquareMatrix& Jinv,
    const point& xP,
    const point& x0
)
{
    blaze::StaticVector<scalar,3UL,blaze::columnVector> v;
    v[0] = xP[0]-x0[0];
    v[1] = xP[1]-x0[1];
    v[2] = xP[2]-x0[2];

    auto vT = Jinv * v;

    return point(vT[0],vT[1],vT[2]);
}


Foam::scalar Foam::geometryWENO::gaussQuad
(
        const int n,
        const int m,
        const int l,
        const point& xi0,
        const vector& v0,
        const vector& v1,
        const vector& v2
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
#if defined(__AVX__)
    auto intPowAVX = [](__m256d& base,const unsigned int exponent) -> void
    {
        if (exponent == 0)
        {
            base = _mm256_set1_pd(1.0);
            return;
        }
        
        const __m256d temp = base;
        
        for (unsigned int i=1; i<exponent; i++)
        {
            base = _mm256_mul_pd(temp,base);
        }
    };
    
    {
    scalar xi[4];
    scalar eta[4];
    scalar zeta[4];
    // go through 1 to 12 as they align in memory
    for (label j = 0; j < 12; j=j+4)
    {
        for (size_t k=0; k < 4; k++)
        {
            xi[k] =
                v0.x()* (1 - gaussCoeff[j+k][0] - gaussCoeff[j+k][1])
              + v1.x()* gaussCoeff[j+k][0] + v2.x()*gaussCoeff[j+k][1] - xi0.x();
            eta[k] =
                v0.y()* (1- gaussCoeff[j+k][0]- gaussCoeff[j+k][1])
              + v1.y()* gaussCoeff[j+k][0] +v2.y()* gaussCoeff[j+k][1] - xi0.y();
            zeta[k] =
                v0.z()* (1- gaussCoeff[j+k][0]- gaussCoeff[j+k][1])
              + v1.z()* gaussCoeff[j+k][0] +v2.z()* gaussCoeff[j+k][1] - xi0.z();
        }
        __m256d mxi   = _mm256_set_pd(  xi[0],   xi[1],   xi[2],   xi[3]);
        __m256d meta  = _mm256_set_pd( eta[0],  eta[1],  eta[2],  eta[3]);
        __m256d mzeta = _mm256_set_pd(zeta[0], zeta[1], zeta[2], zeta[3]);
        
        intPowAVX(mxi,n);
        intPowAVX(meta,m);
        intPowAVX(mzeta,l);

        __m256d mgaussCoeff = _mm256_set_pd(gaussCoeff[j][2],gaussCoeff[j+1][2],
                                            gaussCoeff[j+2][2],gaussCoeff[j+3][2]);

        __m256d temp;
        temp = _mm256_mul_pd(mxi,meta);
        temp = _mm256_mul_pd(temp,mzeta);
        temp = _mm256_mul_pd(mgaussCoeff,temp);

        sum += temp[0]+temp[1]+temp[2]+temp[3];
    }
    }
    // Add the last iteration
    size_t j = 12;
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
    
    return sum;
#else
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
#endif
}



void Foam::geometryWENO::smoothIndIntegrals
(
    const fvMesh& mesh,
    const label cellI,
    const label polOrder,
    const scalarSquareMatrix& JInvI,
    const point& refPointI,
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

        // For some mesh the triangle has a close to zero area. These triangles 
        // are excluded
        if (area < ROOTVSMALL)
            continue;

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
    const point& refPointI,
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
    const blazeList& JInv,
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
    
    // Clear and initialize with zero
    refFacAr.clear();
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
                const triFace& tri = triFaces[i];

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
    scalarSquareMatrix J;
    
    for (int i = 0;i<3;i++)
    {
        for (int j = 0; j<3;j++)
        {
            J(i,j) = pts[referenceFrame[j+1]][i] - pts[referenceFrame[0]][i];
        }
    }
    return J;
}


Foam::scalar Foam::geometryWENO::cond(const scalarSquareMatrix& J)
{
    #ifdef USE_LAPACK
    // Calcualte the eigenvalues of the Jacobi matrix 
        auto sigma = eigen(J);
    #else
        auto sigma = mathFunctionsWENO::eigen(J); 
    #endif
    
    if (min(abs(sigma)) < ROOTVSMALL)
        return GREAT;
    
    return max(abs(sigma))/min(abs(sigma));
}


Foam::geometryWENO::scalarSquareMatrix Foam::geometryWENO::jacobi
(
    const scalar x0, const scalar y0, const scalar z0,
    const scalar x1, const scalar y1, const scalar z1,
    const scalar x2, const scalar y2, const scalar z2,
    const scalar x3, const scalar y3, const scalar z3
)
{
    scalarSquareMatrix J;
    
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


void Foam::geometryWENO::calculateInverseJacobi
(
    const scalarSquareMatrix& J, 
    scalarSquareMatrix& JInvI
)
{
    // Use Blaze functions if LAPACK is available
    #ifdef USE_LAPACK 
        blaze::StaticMatrix<scalar,3UL,3UL,blaze::columnMajor> temp = blaze::inv(J);
       
        // Use LU decomposition to improve the Jacobi Matrix
        blaze::DynamicMatrix<double,blaze::columnMajor> L, U, P;
        lu(temp,L,U,P);
        JInvI = L*U;
    #else
        JInvI = mathFunctionsWENO::inv(J); 
        mathFunctionsWENO::pivot(JInvI);
    #endif
}

// ************************************************************************* //
