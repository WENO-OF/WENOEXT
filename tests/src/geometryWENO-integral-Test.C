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

#include "catch.hpp"

#include "fvCFD.H"
#include "geometryWENO.H"
#include <cmath>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TEST_CASE("geometryWENO: Integration","[3D]")
{
    // Replace setRootCase.H for Catch2   
    int argc = 1;
    char **argv = static_cast<char**>(malloc(sizeof(char*)));
    char executable[] = {'m','a','i','n'};
    argv[0] = executable;
    Foam::argList args(argc, argv,false,false,false);
        
    // create the mesh from case file
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    
    using volIntegralType = List< List< List<scalar> > > ;
    using scalarSquareMatrix = SquareMatrix<scalar>;
    
    //const fvMesh& mesh,
    const label cellI = 33;
    const label polOrder = 2;
    volIntegralType Integral;
    scalarSquareMatrix JInvI;
    point refPointI;
    scalar refDetI;
    
    geometryWENO::initIntegrals(mesh,cellI,polOrder,Integral,JInvI,refPointI,refDetI);
    
    
    
}


