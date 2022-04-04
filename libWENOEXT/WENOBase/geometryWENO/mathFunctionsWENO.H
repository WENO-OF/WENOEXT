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

Description
    An additional math library for the calculation of matrix properties
    if no LAPACK library is present 

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020

\*---------------------------------------------------------------------------*/

#ifndef mathFunctions_H 
#define mathFunctions_H

#include "geometryWENO.H"
#include "cubicEqn.H"
#include "mathematicalConstants.H"
#include <math.h>

#ifdef USE_FMA
#if !defined(FP_FAST_FMA) && defined(__AVX2__)
    #define FP_FAST_FMA
#endif
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Namespace for WENO math functions 
\*---------------------------------------------------------------------------*/

namespace mathFunctionsWENO
{
    //- Calculate eigen values of 3x3 matrix
    //- see: https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices
    blaze::DynamicVector<double,blaze::columnVector> eigen
    (
        const geometryWENO::scalarSquareMatrix& A
    );

    // Calculate the 3x3 determinant
    double det(const blaze::DynamicMatrix<double>& A);
    
    //- Calculate the 2x2 determinant
    //  solves the system A = a b
    //                        c d 
    //  with det2 = a*d - c*b
    double det2 (
        const double a, 
        const double d,
        const double c,
        const double b
    );

    //- Stablize a 3x3 matrix by calculating the pivoted matrix A' = P*A
    void pivot(geometryWENO::scalarSquareMatrix& A);

    //- Invert a 3x3 matrix
    geometryWENO::scalarSquareMatrix inv(const geometryWENO::scalarSquareMatrix& A);

} // End namespace mathFunctionsWENO 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
blaze::DynamicVector<double,blaze::columnVector> 
Foam::mathFunctionsWENO::eigen(const geometryWENO::scalarSquareMatrix& A)
{
    // Eigen values are computed in the same manner as OpenFOAM 
    // tensor.C does for eigenValues()
    
    // Create blaze vector with three entries
    blaze::DynamicVector<double,blaze::columnVector> lambda(3); 
	
	// Coefficients of the characteristic cubic polynomial (a = 1)
    const scalar b =
       - A(0,0) - A(1,1) - A(2,2);
    const scalar c =
         A(0,0)*A(1,1) + A(0,0)*A(2,2) + A(1,1)*A(2,2)
       - A(0,1)*A(1,0) - A(1,2)*A(2,1) - A(2,0)*A(0,2);
    const scalar d =
       - A(0,0)*A(1,1)*A(2,2)
       - A(0,1)*A(1,2)*A(2,0) - A(0,2)*A(2,1)*A(1,0)
       + A(0,0)*A(1,2)*A(2,1) + A(1,1)*A(2,0)*A(0,2) + A(2,2)*A(0,1)*A(1,0);
 
     // Solve
     Roots<3> roots = cubicEqn(1, b, c, d).roots();
 
     // Check the root types
     forAll(roots, i)
     {
         switch (roots.type(i))
         {
             case roots::real:
                 lambda[i] = roots[i];
                 break;
             case roots::complex:
                 lambda[i] = 0;
                 break;
             case roots::posInf:
                 lambda[i] = 1E+15;
                 break;
             case roots::negInf:
                 lambda[i] = -1E+15;
                 break;
             case roots::nan:
                 FatalErrorInFunction
                     << "Eigenvalue calculation failed for tensor: " 
                     << exit(FatalError);
         }
     }
     
     return lambda;
}


double
Foam::mathFunctionsWENO::det(const blaze::DynamicMatrix<double>& A)
{
    return   det2
             (
                A(0,0),det2(A(1,1),A(2,2),A(1,2),A(2,1)),
                A(0,1),det2(A(1,0),A(2,2),A(1,2),A(2,0))
             )
           + A(0,2)*det2(A(1,0),A(2,1),A(1,1),A(2,0));
}


double
Foam::mathFunctionsWENO::det2
(
    const double a, 
    const double d,
    const double c,
    const double b
)
{
    // Solves the determinant for a 2x2 matrix of the form
    // A = a  b
    //     c  d
    // with det2 = a*d - c*b
    // see also: https://pharr.org/matt/blog/2019/11/03/difference-of-floats
    double cb = c * b;
    #ifdef FP_FAST_FMA
        double err = std::fma(-c, b, cb);
        double dop = std::fma(a, d, -cb);
        return dop + err;
    #else
        return (a*d) - cb;
    #endif
}


void Foam::mathFunctionsWENO::pivot(geometryWENO::scalarSquareMatrix& A)
{
    // create the identity matrix
    constexpr geometryWENO::scalarSquareMatrix I
    {
        {1,0,0},
        {0,1,0},
        {0,0,1}
    };

    // Find the largest element in the first column
    int maxRowIndex = 0;
    double Amax = -1E+15;
    for (unsigned int i=0; i < A.rows(); i++)
    {
        if (std::abs(A(i,0)) > Amax)
        {
            maxRowIndex = i;
            Amax = std::abs(A(i,0));
        }
    }

    // create pivot matrix
    geometryWENO::scalarSquareMatrix P1 = I;
    // create a copy of the first row
    blaze::DynamicVector<double,blaze::rowVector> row0;
    row0 = blaze::row<0>(P1);
    auto rowI = blaze::row(P1,maxRowIndex);
    row(P1,0) = rowI;
    row(P1,maxRowIndex) = row0;

    geometryWENO::scalarSquareMatrix PA = P1*A;

    // execute firs step of the gauss elemination process
    geometryWENO::scalarSquareMatrix M1 = I;
    M1(1,0) = -1.*PA(1,0)/PA(0,0);
    M1(2,0) = -1.*PA(2,0)/PA(0,0);

    PA = M1*PA;

    // Find largest element in the second column 
    Amax = -1E+15;
    for (unsigned int i=1; i < A.rows(); i++)
    {
        if (std::abs(A(i,1)) > Amax)
        {
            maxRowIndex = i;
            Amax = std::abs(A(i,1));
        }
    }

    // Permutation matrix is only required if maxRowIndex != 1
    if (maxRowIndex != 1)
    {
        geometryWENO::scalarSquareMatrix P2 = I;
        blaze::DynamicVector<double,blaze::rowVector> row1;
        row1 = blaze::row<1>(P2);
        auto row2 = blaze::row<2>(P2);
        row(P2,1) = row2;
        row(P2,2) = row1;

        P1 = P2*P1;
    }

    // permutate the matrix A
    A = P1*A;
}

    
Foam::geometryWENO::scalarSquareMatrix 
Foam::mathFunctionsWENO::inv(const geometryWENO::scalarSquareMatrix& A)
{
    // Based on blaze::invertGeneral3x3
    geometryWENO::scalarSquareMatrix AInv;
    
    AInv(0, 0) = det2(A(1, 1),A(2, 2),A(2, 1),A(1, 2));
    AInv(1, 0) = det2(A(1, 2),A(2, 0),A(1, 0),A(2, 2));
    AInv(2, 0) = det2(A(1, 0),A(2, 1),A(2, 0),A(1, 1));
    
    const double Det = A(0,0)*AInv(0,0) + A(0,1)*AInv(1,0) + A(0,2)*AInv(2,0);
    
    #ifdef FULLDEBUG
        if (Det == 0)
            FatalErrorInFunction 
                << "Determinant is zero! Division by zero occurs"
                << exit(FatalError);
    #endif

    AInv(0, 1) = det2(A(0, 2),A(2, 1),A(0, 1),A(2, 2));
    AInv(0, 2) = det2(A(0, 1),A(1, 2),A(0, 2),A(1, 1));
    AInv(1, 1) = det2(A(0, 0),A(2, 2),A(0, 2),A(2, 0));
    AInv(1, 2) = det2(A(1, 0),A(0, 2),A(0, 0),A(1, 2));
    AInv(2, 1) = det2(A(2, 0),A(0, 1),A(0, 0),A(2, 1));
    AInv(2, 2) = det2(A(0, 0),A(1, 1),A(1, 0),A(0, 1));

    return AInv/Det;
}


#endif

// ************************************************************************* //
