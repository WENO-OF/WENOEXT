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


Application
    mathFunctionsWENO-Test

Description
    Test the math functions
    
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Catch::Approx is needed when floats are compared

#include "mathFunctionsWENO.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Compare complex eigen values with expected real values
// Uses an absolute tolerance similar to the former relaxed Blaze comparison
void checkEigenValues
(
    const Eigen::VectorXcd& eigVal,
    const std::vector<double>& expected
)
{
    REQUIRE(eigVal.size() == Eigen::Index(expected.size()));
    for (Eigen::Index i = 0; i < eigVal.size(); i++)
    {
        REQUIRE(eigVal[i].real() == Catch::Approx(expected[i]).margin(1E-8));
        REQUIRE(eigVal[i].imag() == Catch::Approx(0).margin(1E-8));
    }
}


TEST_CASE("mathFunctionsWENO","[baseTest],[mathFunctions]")
{
    SECTION("determinant of 3x3 matrix")
    {
        // Check the determinante of 3x3 matrix
        const geometryWENO::scalarSquareMatrix A
        {
            {2,-5,3},
            {0,7,-2},
            {-1,4,1}
        };
        REQUIRE(mathFunctionsWENO::det(A) == 41);
    }

    // Check the eigen values 
    SECTION("eigenvalues of 3x3 matrix")
    {
        {
            INFO("With valid eigen values");
            // Check the determinante of 3x3 matrix
            const geometryWENO::scalarSquareMatrix A
            {
                { 5, 2, 0},
                { 2, 5, 0},
                {-3, 4, 6}
            };

            checkEigenValues(mathFunctionsWENO::eigen(A), {7,6,3});
        }
        
        {
            INFO("With one zero entry");
            // Check for a badly conditioned matrix system
            const geometryWENO::scalarSquareMatrix B
            {
                { 1, 1, 0.1},
                { 0, 2, 0},
                { 1, 0, 0.1}
            };

            checkEigenValues(mathFunctionsWENO::eigen(B), {2,1.1,0});
        }
        
        {
            INFO("With two zero eigen values");
            // Check the determinante of 3x3 matrix
            const geometryWENO::scalarSquareMatrix A
            {
                {-5.6555245e-05,  5.6558097e-05, 2.873191e-05},
                { 4.92445437e-05, -4.9247027e-05, 3.27478754e-05},
                { 0     ,  0     ,-3.46944695e-18}
            };
            
            checkEigenValues
            (
                mathFunctionsWENO::eigen(A),
                {
                    1.70008e-14,
                    -3.46945e-18,
                    -0.000105802
                }
            );
        }
    }

    SECTION("inverse of 3x3 matrix")
    {
        const geometryWENO::scalarSquareMatrix A
        {
            { 0,-3,-2},
            { 1,-4,-2},
            {-3, 4, 1}
        };

        const geometryWENO::scalarSquareMatrix AInv
        {
            { 4,-5,-2},
            { 5,-6,-2},
            {-8, 9, 3}
        };

        REQUIRE(mathFunctionsWENO::inv(A) == AInv); 
    }
}

