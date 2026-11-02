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
#include "realEigenValues.H"
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
    //  Based on Numerical Recipes 3rd Edition Chapter 11.6
    Eigen::VectorXcd eigen
    (
        const geometryWENO::scalarSquareMatrix& A
    );

    // Calculate the 3x3 determinant
    double det(const geometryWENO::scalarSquareMatrix& A);
    
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
Eigen::VectorXcd
Foam::mathFunctionsWENO::eigen(const geometryWENO::scalarSquareMatrix& A)
{
    realEigenValues eigenVal(A);

    return eigenVal.eig();
}


double
Foam::mathFunctionsWENO::det(const geometryWENO::scalarSquareMatrix& A)
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
    const geometryWENO::scalarSquareMatrix I =
        geometryWENO::scalarSquareMatrix::Identity();

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
    // swap the first row with the row of the largest element
    P1.row(0).swap(P1.row(maxRowIndex));

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
        P2.row(1).swap(P2.row(2));

        P1 = P2*P1;
    }

    // permutate the matrix A
    A = P1*A;
}

    
Foam::geometryWENO::scalarSquareMatrix 
Foam::mathFunctionsWENO::inv(const geometryWENO::scalarSquareMatrix& A)
{
    // Based on the cofactor (adjugate) formula of a 3x3 matrix
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
