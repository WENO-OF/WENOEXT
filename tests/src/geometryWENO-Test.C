/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    geometryWENO-Test

Description
    Test geometryWENO class with Catch2
    
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file


#include <catch2/catch.hpp>

#include "fvCFD.H"
#include "geometryWENO.H"
#include <cmath>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TEST_CASE("geometryWENO Jacobi Class")
{
    using scalarSquareMatrix = SquareMatrix<scalar>;
    
    pointField pts(4,vector(0,0,0));
    
    // Populate points
    pts[0] =vector(0,0,0);
    pts[1] =vector(1,0,0);
    pts[2] =vector(0,2,0);
    pts[3] =vector(0,0,3);
    
    labelList referenceFrame(4);
    std::iota(referenceFrame.begin(), referenceFrame.end(),0);
    
    scalarSquareMatrix J = geometryWENO::jacobi(pts,referenceFrame);
    
    for (int i = 0; i<J.n();i++)
    {
        for (int j = 0; j < J.n(); j++)
        {
            if (i==j)
                CHECK(Approx(J(i,j)) == i+1);
            else
                CHECK(Approx(J(i,j)) == 0);
        }
    }
    
    // Check the determinante
    CHECK(Approx(det(J)) == 6);
    
    // Check the inverse
    scalarRectangularMatrix JInv = geometryWENO::JacobiInverse(J);
    
        for (int i = 0; i<J.n();i++)
    {
        for (int j = 0; j < J.n(); j++)
        {
            if (i==j)
                CHECK(Approx(JInv(i,j)) == 1.0/(i+1));
            else
                CHECK(Approx(JInv(i,j)) == 0);
        }
    }
    
}
