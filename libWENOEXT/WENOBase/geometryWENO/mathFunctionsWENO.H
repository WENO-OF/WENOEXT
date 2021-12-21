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
    //  see: https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices
    blaze::DynamicVector<double,blaze::columnVector> eigen
    (
        const geometryWENO::scalarSquareMatrix& A
    );

    // Calculate the 3x3 determinant
    double det(const blaze::DynamicMatrix<double>& A);
    
    // Calculate the 2x2 determinant
    // solves the system A = a b
    //                       c d 
    // with det2 = a*d - c*b
    double det2 (
        const double a, 
        const double d,
        const double c,
        const double b
    );

    // Stablize a 3x3 matrix by calculating the pivoted matrix A' = P*A
    void pivot(geometryWENO::scalarSquareMatrix& A);

    geometryWENO::scalarSquareMatrix inv(const geometryWENO::scalarSquareMatrix& A);

} // End namespace mathFunctionsWENO 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
blaze::DynamicVector<double,blaze::columnVector> 
Foam::mathFunctionsWENO::eigen(const geometryWENO::scalarSquareMatrix& A)
{
    // Create blaze vector with three entries
    blaze::DynamicVector<double,blaze::columnVector> eig(3);

    // Create identity matrix
    geometryWENO::scalarSquareMatrix I;
    I = blaze::declid(I);
    
    const double p1 = A(0,1)*A(0,1) + A(0,2)*A(0,2) + A(1,2)*A(1,2);
    if (p1 == 0)
    {
        eig[0] = A(0,0);
        eig[1] = A(1,1);
        eig[2] = A(2,2);
    }
    else
    {
        const double q  = blaze::trace(A)/3.0;
        const double p2 =   std::pow((A(0,0) - q),2)
                          + std::pow((A(1,1) - q),2)
                          + std::pow((A(2,2) - q),2)
                          + 2.*p1;
        const double p  = sqrt(p2/6.0);
        const auto   B  = (A - (q*I))/p;
        const double r  = 0.5*det(B);
        // In exact arithmetic for a symmetric matrix  -1 <= r <= 1
        // but computation error can leave it slightly outside this range
        double phi;
        if (r <= -1)
            phi = M_PI/3.0;
        else if (r >= 1)
            phi = 0;
        else
            phi = std::acos(r)/3.0;

        eig[0] = q + 2.0*p*std::cos(phi);
        eig[1] = q + 2.0*p*std::cos(phi + (2.0*M_PI/3.0));
        eig[2] = 3*q - eig[0] - eig[1];
    }
    return eig;
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
    const double invDet = 1.0/det(A);

    geometryWENO::scalarSquareMatrix AInv;

    AInv(0, 0) = det2(A(1, 1),A(2, 2),A(2, 1),A(1, 2)) * invDet;
    AInv(0, 1) = det2(A(0, 2),A(2, 1),A(0, 1),A(2, 2)) * invDet;
    AInv(0, 2) = det2(A(0, 1),A(1, 2),A(0, 2),A(1, 1)) * invDet;
    AInv(1, 0) = det2(A(1, 2),A(2, 0),A(1, 0),A(2, 2)) * invDet;
    AInv(1, 1) = det2(A(0, 0),A(2, 2),A(0, 2),A(2, 0)) * invDet;
    AInv(1, 2) = det2(A(1, 0),A(0, 2),A(0, 0),A(1, 2)) * invDet;
    AInv(2, 0) = det2(A(1, 0),A(2, 1),A(2, 0),A(1, 1)) * invDet;
    AInv(2, 1) = det2(A(2, 0),A(0, 1),A(0, 0),A(2, 1)) * invDet;
    AInv(2, 2) = det2(A(0, 0),A(1, 1),A(1, 0),A(0, 1)) * invDet;

    return AInv;
}


#endif

// ************************************************************************* //
