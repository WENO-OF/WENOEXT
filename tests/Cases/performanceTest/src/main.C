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
    WENOUpwindFit interpolation test

Description
    Test the performance of the WENOUpwindFit scheme. 
    1. Section: Build up of the required matrix list 
    2. Section: Run time during execution
    
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include <chrono>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main()
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

    const vectorField& centre = mesh.C();

    // ------------------------------------------------------------------------
    //                          Create Input Data 
    // ------------------------------------------------------------------------
    
    
    // Mesh has two patches, outlet with fixed value and empty direction
    wordList patchTypes(2);
    patchTypes[mesh.boundary()["outlet"].index()] = "fixedValue"; 
    patchTypes[mesh.boundary()["topBottom"].index()] = "fixedValue";
    

    auto psiFunc = [](const double x, const double y) -> double
    {
        return std::sin(x);
    };

    // Create a volScalarField with the sinus curve
    volScalarField psi
    (
        IOobject
        (
            "psi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0),
        patchTypes
    );

    // populate psi
    forAll(mesh.C(),celli)
    {
        psi[celli] = psiFunc(centre[celli].x(),centre[celli].y());
    }
    
    {
        auto& PatchField = psi.boundaryFieldRef();
        forAll(PatchField,patchi)
        {
            const vectorField& faceCenters = PatchField[patchi].patch().Cf();
            forAll(PatchField[patchi],facei)
            {
                PatchField[patchi][facei] = psiFunc(faceCenters[facei].x(),faceCenters[facei].y());
            }
        }
    }
    psi.write();
    
    // Create velocity field 
    volVectorField U
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimVelocity, vector(0.0,0.0,0.0)),
        patchTypes
    );
    
    forAll(mesh.C(),celli)
    {
        U[celli] = vector(1,1,0);
    }
    // Correct the boundary
    // Find patch outlet
    {
        const label outletIndex = mesh.boundary().findPatchID("outlet");
        auto& PatchField = U.boundaryFieldRef()[outletIndex];
    
        forAll(PatchField,facei)
        {
            PatchField[facei] = vector(1,1,0);
        }
    }
    
    U.write();
    
    // Get the surface fields
    surfaceScalarField phi = fvc::flux(U);


    auto t1 = std::chrono::high_resolution_clock::now();
    volScalarField divWENO("divWENO",fvc::div(phi,psi,"div(WENO)"));
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    Info << "Duration Build-Up: "<<duration/1E+6<<" seconds"<<endl;
    
    
    t1 = std::chrono::high_resolution_clock::now();
    for (int i=0; i < 1000; ++i)
    {
        volScalarField divWENO("divWENO",fvc::div(phi,psi,"div(WENO)"));
    }
    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    Info << "Duration Run-Time: "<<duration/1E+6<<" seconds"<< endl;
    
    
    t1 = std::chrono::high_resolution_clock::now();
    for (int i=0; i < 1000; ++i)
    {
        volScalarField divWENO("divWENO",fvc::div(phi,psi,"div(Linear)"));
    }
    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    Info << "Duration Run-Time Linear: "<<duration/1E+6<<" seconds"<<endl;

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

