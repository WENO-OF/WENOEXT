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


#include "catch.hpp"

#include "fvCFD.H"
#include "geometryWENO.H"
#include <cmath>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TEST_CASE("geometryWENO: Jacobi Matrix")
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
    
    // Construct empty geometryWENO object
    geometryWENO geo;
    
    scalarSquareMatrix J = geo.jacobi(pts,referenceFrame);
    
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
    
    // Create Jacobi from components
    scalarSquareMatrix J2 = geo.jacobi
    (
        pts[referenceFrame[0]][0], pts[referenceFrame[0]][1],
        pts[referenceFrame[0]][2], pts[referenceFrame[1]][0],
        pts[referenceFrame[1]][1], pts[referenceFrame[1]][2],
        pts[referenceFrame[2]][0], pts[referenceFrame[2]][1],
        pts[referenceFrame[2]][2], pts[referenceFrame[3]][0],
        pts[referenceFrame[3]][1], pts[referenceFrame[3]][2]
    );
    
    for (int i = 0; i<J2.n();i++)
    {
        for (int j = 0; j < J2.n(); j++)
        {
            if (i==j)
                CHECK(Approx(J2(i,j)) == i+1);
            else
                CHECK(Approx(J2(i,j)) == 0);
        }
    }
    
    
    // Check the determinante
    CHECK(Approx(det(J)) == 6);
    
    // Check the inverse
    scalarSquareMatrix JInv = geo.JacobiInverse(J);
    
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
    
    
    // Check transform point
    const point x0(1,0,0);
    const point xp(2,2,3);
    const point res = geo.transformPoint(JInv,xp,x0);
    
    CHECK(Approx(res[0])==1);
    CHECK(Approx(res[1])==1);
    CHECK(Approx(res[2])==1);

}



TEST_CASE("geometryWENO: Integration")
{
    // Replace setRootCase.H for Catch2   
    int argc = 1;
    char **argv = static_cast<char**>(malloc(sizeof(char*)));
    char executable[] = {'m','a','i','n'};
    argv[0] = executable;
    Foam::argList args(argc, argv,false,false,false);
        
        
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    
    geometryWENO geo;
    
    using scalarMatrix = List< List< List<scalar> > > ;
    using scalarSquareMatrix = SquareMatrix<scalar>;
    
    //const fvMesh& mesh,
    const label cellI = 33;
    const label polOrder = 2;
    scalarMatrix Integral;
    scalarSquareMatrix JInvI;
    point refPointI;
    scalar refDetI;
    
    geo.initIntegrals(mesh,cellI,polOrder,Integral,JInvI,refPointI,refDetI);
    
    Info << Integral << endl;
    
}


