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
        // Check the determinante of 3x3 matrix
        const geometryWENO::scalarSquareMatrix A
        {
            { 5, 2, 0},
            { 2, 5, 0},
            {-3, 4, 6}
        };

        blaze::DynamicVector<double,blaze::columnVector> eigVal{7,6,3};

        CHECK(mathFunctionsWENO::eigen(A) == eigVal);
        
        
        // Check for a badly conditioned matrix system
        const geometryWENO::scalarSquareMatrix B
        {
            { 1, 1, 0.1},
            { 0, 2, 0},
            { 1, 0, 0.1}
        };

        blaze::DynamicVector<double,blaze::columnVector> eigVal2{2,1.1,0};

        REQUIRE(mathFunctionsWENO::eigen(B) == eigVal2);
        
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

