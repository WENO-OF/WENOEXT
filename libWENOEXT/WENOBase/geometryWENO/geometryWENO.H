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

Class
    Foam::geometryWENO

Description
    Class for geometric calculations of WENO reconstructions

SourceFiles
    geometryWENO.C

Author
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020

\*---------------------------------------------------------------------------*/

#ifndef geometryWENO_H
#define geometryWENO_H

#include "triSurface.H"
#include "polyMeshTetDecomposition.H"
#include "fvMesh.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class geometryWENO Declaration
\*---------------------------------------------------------------------------*/

namespace geometryWENO
{

    //- Typedef for 3D scalar matrix
    using volIntegralType = List< List< List<scalar> > > ;
    using scalarSquareMatrix = SquareMatrix<scalar>;

    // Member Functions

        //- Evaluate the surface integral using Gaussian quadrature
        scalar gaussQuad
        (
            const scalar n,
            const scalar m,
            const scalar l,
            const point x0,
            const vector v0,
            const vector v1,
            const vector v2
        );

        //- Calculate integral of Eq. (3.20) in Master Thesis 
        //  for any combination of l,m,n with the restrain '(l+m+n) < r'
        void initIntegrals
        (
            const fvMesh& mesh,
            const label cellI,
            const label polOrder,
            volIntegralType& volIntegrals,
            scalarSquareMatrix& JInvI,
            point& refPointI,
            scalar& refDetI
        );

        //- Calculate the volume integrals of stencil cells in reference space
        //- of owner cell
        volIntegralType transformIntegral
        (
            const fvMesh& mesh,
            const label cellJ,
            const point transCenterJ,
            const label polOrder,
            const scalarSquareMatrix& JInvI,
            const point refPointI,
            const scalar refDetI
        );

        //- Calculate inverse of Jacobian matrix
        scalarSquareMatrix JacobiInverse
        (
            const scalarSquareMatrix& J
        );

        //- Transform an arbitrary point into reference space of owner cell
        point transformPoint
        (
            const scalarSquareMatrix& Jinv,
            const point xP,
            const point x0
        );


        //- Calculate the volume integrals for the smoothness indicator matrix
        volIntegralType smoothIndIntegrals
        (
            const fvMesh& mesh,
            const label cellI,
            const label polOrder,
            const scalarSquareMatrix& JacobianI,
            const point refPointI
        );

        //- Calculate the smoothness indicator matrices
        scalarRectangularMatrix getB
        (
            const fvMesh& mesh,
            const label cellI,
            const label polOrder,
            const label nDvt,
            const scalarSquareMatrix& JInvI,
            const point refPointI,
            const labelList& dim
        );

        //- Calculate factorials of variable
        scalar Fac(label x);

        //- Calculation of surface integrals for convective terms
        void surfIntTrans
        (
            const fvMesh& mesh,
            const label polOrder,
            const List<volIntegralType>& volIntegralsList,
            const List<scalarSquareMatrix>& JInv,
            const List<point>& refPoint,
            List<Pair<volIntegralType> >& intBasTrans,
            List<scalar>& refFacAr
        );

        vector compCheck
        (
            const label n,
            const label m,
            const label l,
            const volIntegralType& intBasisfI
        );
        
        //- Create Jacobi matrix from pointField and label
        scalarSquareMatrix jacobi
        (
            const pointField& pts,
            const labelList& referenceFrame
        );
        
        //- Create Jacobi matrix from entry list
        scalarSquareMatrix jacobi
        (
            const scalar x0, const scalar y0, const scalar z0,
            const scalar x1, const scalar y1, const scalar z1,
            const scalar x2, const scalar y2, const scalar z2,
            const scalar x3, const scalar y3, const scalar z3
        );

} // End namespace geometryWENO

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
