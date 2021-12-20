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

Class
    Foam::geometryWENO

Description
    Class for geometric calculations of WENO reconstructions

SourceFiles
    geometryWENO.C
    
    [1] T. Pringuey, “Large Eddy Simulation of Primary Breakup,” 2012.

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
#include "blaze/Math.h"
#include "BlazeIO.H"
#include "List3D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class geometryWENO Declaration
\*---------------------------------------------------------------------------*/

namespace geometryWENO
{

    //- Typedef for 3D scalar matrix
    using volIntegralType = List3D<scalar>;
    using scalarSquareMatrix = blaze::StaticMatrix<scalar,3UL,3UL,blaze::rowMajor>;
    using DynamicMatrix = blaze::DynamicMatrix<scalar>;
    
    //- For static matrices a custom allocator is required when dynamic
    //  allocation is used. OpenFOAM does not provide an option for an 
    //  allocator with a list. 
    using blazeList = std::vector<geometryWENO::scalarSquareMatrix,blaze::AlignedAllocator<scalar>>;
    
    // Member Functions

        //- Evaluate the surface integral using Gaussian quadrature
        scalar gaussQuad
        (
            const int n,
            const int m,
            const int l,
            const point& x0,
            const vector& v0,
            const vector& v1,
            const vector& v2
        );

        //- Calculate integral of Eq. (6.23) in Pringuey PhD thesis [1]
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
        void transformIntegral
        (
            const fvMesh& mesh,
            const label cellJ,
            const point& transCenterJ,
            const label polOrder,
            const scalarSquareMatrix& JInvI,
            const point& refPointI,
            const scalar refDetI,
            volIntegralType& transVolMom
        );
        

        //- Transform an arbitrary point into reference space of owner cell
        point transformPoint
        (
            const scalarSquareMatrix& Jinv,
            const point& xP,
            const point& x0
        );


        //- Calculate the volume integrals for the smoothness indicator matrix
        void smoothIndIntegrals
        (
            const fvMesh& mesh,
            const label cellI,
            const label polOrder,
            const scalarSquareMatrix& JacobianI,
            const point& refPointI,
            volIntegralType& smoothVolIntegral
        );

        //- Calculate the smoothness indicator matrices
        DynamicMatrix getB
        (
            const fvMesh& mesh,
            const label cellI,
            const label polOrder,
            const label nDvt,
            const scalarSquareMatrix& JInvI,
            const point& refPointI,
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
            const blazeList& JInv,
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

        //- Calculate the matrix condition of matrix J
        scalar cond(const scalarSquareMatrix& J);

        // Const array for Gauss quadratur 7th order
        // Coefficients can be found in 
        //   Quadrature Formulas in Two DimensionsMath 5172 
        //   - Finite Element MethodSection 001, Spring 2010
        // Points and weights for the Gaussian quadrature of a standard triangle
        // - 1.row: x-values
        // - 2.row: y-values
        // - 3.row: weights
        const static scalar gaussCoeff[13][3] = 
        {
            {0.33333333333333, 0.33333333333333, -0.14957004446768}, // 0
            {0.26034596607904, 0.26034596607904,  0.17561525743321}, // 1
            {0.26034596607904, 0.47930806784192,  0.17561525743321}, // 2
            {0.47930806784192, 0.26034596607904,  0.17561525743321}, // 3
            {0.06513010290222, 0.06513010290222,  0.05334723560884}, // 4
            {0.06513010290222, 0.86973979419557,  0.05334723560884}, // 5
            {0.86973979419557, 0.06513010290222,  0.05334723560884}, // 6
            {0.31286549600487, 0.63844418856981,  0.07711376089026}, // 7
            {0.63844418856981, 0.04869031542532,  0.07711376089026}, // 8
            {0.04869031542532, 0.31286549600487,  0.07711376089026}, // 9
            {0.63844418856981, 0.31286549600487,  0.07711376089026}, // 10
            {0.31286549600487, 0.04869031542532,  0.07711376089026}, // 11
            {0.04869031542532, 0.63844418856981,  0.07711376089026} // 12
        };
        // Coeffcients for 5th order
        //const static scalar gaussCoeff[7][3] = 
        //{
            //{0.33333333333333, 0.33333333333333, 0.22500000000000},
            //{0.47014206410511, 0.47014206410511, 0.13239415278851},
            //{0.47014206410511, 0.05971587178977, 0.13239415278851},
            //{0.05971587178977, 0.47014206410511, 0.13239415278851},
            //{0.10128650732346, 0.10128650732346, 0.12593918054483},
            //{0.10128650732346, 0.79742698535309, 0.12593918054483},
            //{0.79742698535309, 0.10128650732346, 0.12593918054483}
        //};
        
        //- Calculate the inverse of the Jacobi Matrix J and perform
        //  additional checks and modifications to improve the condition
        void calculateInverseJacobi
        (
            const scalarSquareMatrix& J, 
            scalarSquareMatrix& JInvI
        );

} // End namespace geometryWENO

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
