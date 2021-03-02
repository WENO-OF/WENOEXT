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
    WENOUpwindFit boundary condition test 
    
Description
    Test following boundary conditions
    
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#include "catch.hpp"
#include "fvCFD.H"
#include "WENOUpwindFit.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TEST_CASE("WENOUpwindFit Boundary Test","[bCyclic]")
{
    // ------------------------------------------------------------------------
    //                          OpenFOAM Start-Up 
    // ------------------------------------------------------------------------
    // Replace setRootCase.H for Catch2   
    int argc = 1;
    char **argv = static_cast<char**>(malloc(sizeof(char*)));
    char executable[] = {'m','a','i','n'};
    argv[0] = executable;
    Foam::argList args(argc, argv,false,false,false);
    #include "createTime.H"        // create the time object
    #include "createMesh.H"        // create the mesh object

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    // Create velocity field 
    volVectorField U
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    
    volScalarField psi
    (
        IOobject
        (
            "psi",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    
    
    // Get the surface fields
    surfaceScalarField phi = fvc::flux(U);

    // Create WENOUpwindFit object
    WENOUpwindFit<scalar> WENO(mesh,phi,mesh.divScheme("WENO"));
    
    auto corrField = WENO.correction(psi);
    
    // Check the cyclic patches left/right
    const auto bField = corrField().boundaryField();
    
    const fvPatchList& patches = mesh.boundary();
    
    forAll(bField,patchI)
    {
        if (isA<cyclicFvPatch>(patches[patchI]))
        {
            forAll(bField[patchI],faceI)
            {
                CHECK(bField[patchI][faceI] == 0); 
            }
        }
    }
}


TEST_CASE("WENOUpwindFit Boundary Test cyclicAMI","[bCyclicAMI]")
{
    // ------------------------------------------------------------------------
    //                          OpenFOAM Start-Up 
    // ------------------------------------------------------------------------
    // Replace setRootCase.H for Catch2   
    int argc = 1;
    char **argv = static_cast<char**>(malloc(sizeof(char*)));
    char executable[] = {'m','a','i','n'};
    argv[0] = executable;
    Foam::argList args(argc, argv,false,false,false);
    #include "createTime.H"        // create the time object
    #include "createMesh.H"        // create the mesh object

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    // Create velocity field 
    volVectorField U
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    
    volScalarField psi
    (
        IOobject
        (
            "psi",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    
    // Get the surface fields
    surfaceScalarField phi = fvc::flux(U);

    // Create WENOUpwindFit object
    WENOUpwindFit<scalar> WENO(mesh,phi,mesh.divScheme("WENO"));
    
    auto corrField = WENO.correction(psi);
    
    // Check the cyclic patches left/right
    const auto bField = corrField().boundaryField();
    
    const fvPatchList& patches = mesh.boundary();
    
    forAll(bField,patchI)
    {
        if (isA<cyclicAMIFvPatch>(patches[patchI]))
        {
            forAll(bField[patchI],faceI)
            {
                INFO("Checking patch field "<<patches[patchI].name())
                REQUIRE(bField[patchI][faceI] == 0); 
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

