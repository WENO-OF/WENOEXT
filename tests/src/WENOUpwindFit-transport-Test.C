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
    Test the transport with the WENO scheme within a steady state convection
    diffusion problem with the equation:
    \verbatim
        div(phi,Y) == laplacian(D,Y)
    \endverbatim
    
    This has the analytical solution:
    \verbatim
        Y(X) = Y(0) + Y(L) - Y(0)*(exp(phi/D*X)-1)/(exp(phi/D*L)-1)
    \endverbatim
    
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#include "catch.hpp"

#include "fvCFD.H"
#include "gaussConvectionScheme.H"
#include "WENOUpwindFit.H"
#include "upwind.H"
#include "IStringStream.H"
#include <cmath>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TEST_CASE("WENOUpwindFit Transport","[1D][transport]")
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
    

    // Diffusion Coefficient
    scalar Gamma = 0.1;
    // Maximum length of the domain is 1m
    scalar L = 1.0;
    
    
    // Read the field 
    volScalarField Y
    (
        IOobject
        (
            "Y",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    // Copy the field to solve it with linear scheme 
    volScalarField YLinear("YLinear",Y);
    
    // Copy the field to solve it with upwind scheme 
    volScalarField YUpwind("YUpwind",Y);
    
    volScalarField D
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("D",pow(dimLength,2)/dimTime,Gamma)
    );
    
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    
    scalar flux = U[0][0];
    
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux(U)
    );
    

    
    while (runTime.run())
    {
        // Create Stream to read pol Order from
        IStringStream stream("1 0");
        
        runTime++;
        solve
        (
            fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                WENOUpwindFit<scalar>(mesh, phi,stream)
            ).fvmDiv(phi, Y)
          ==
            fvm::laplacian(D,Y)
        );
        
        solve
        (
            fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                linear<scalar>(mesh)
            ).fvmDiv(phi, YLinear)
          ==
            fvm::laplacian(D,YLinear)
        );
        
        
        solve
        (
            fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, YUpwind)
          ==
            fvm::laplacian(D,YUpwind)
        );
        
        runTime.write();
    }
    YLinear.write();
    YUpwind.write();
    
    // Calculate the analytical solution on the field 
    volScalarField analyticalSolu("analyticalSolu",Y);
    
    const label patchIDInlet = mesh.boundaryMesh().findPatchID("inlet");
    const label patchIDOutlet = mesh.boundaryMesh().findPatchID("outlet");
    
    forAll(analyticalSolu,celli)
    {
        // Get the cell center position:
        scalar X = mesh.C()[celli] & vector(1,0,0);
        analyticalSolu[celli] =   Y.boundaryField()[patchIDOutlet][0]
                                + Y.boundaryField()[patchIDInlet][0]
                                *(
                                    1.0
                                  - (std::exp(flux/Gamma*X)-1)
                                    /(std::exp(flux/Gamma*L)-1)
                                 );
    }
    analyticalSolu.write();
    
    // The upwind solution is the lower bound 
    // Check that the Y solution is always equal or above
    forAll(Y,celli)
    {
        REQUIRE(Y[celli] >= Approx(YUpwind[celli]).margin(1e-3));
        //REQUIRE(Y[celli] <= Approx(analyticalSolu[celli]).margin(0.001));
    }
    
    // Calculate the maximum and the mean error to the analytical solution
    double maxError = 0;
    double meanErrorWENO = 0;
    double meanErrorUpwind = 0;
    double meanErrorLinear = 0;
    forAll(Y,celli)
    {
        double errorWENO = std::abs(Y[celli] - analyticalSolu[celli])/analyticalSolu[celli];
        double errorUpwind = std::abs(YUpwind[celli] - analyticalSolu[celli])/analyticalSolu[celli];
        double errorLinear = std::abs(YLinear[celli] - analyticalSolu[celli])/analyticalSolu[celli];
        if (errorWENO > maxError)
            maxError = errorWENO;
        meanErrorWENO += errorWENO;
        meanErrorUpwind += errorUpwind;
        meanErrorLinear += errorLinear;
    }
    meanErrorWENO /= Y.size();
    meanErrorUpwind /= Y.size();
    meanErrorLinear /= Y.size();
    
    Info << nl << "***********************************************"<<endl;
    Info << "Maximum error WENO: \t" << maxError*100.0 <<"[%]"<< endl;
    Info << "Mean error WENO:  \t"  << meanErrorWENO*100.0   <<"[%]"<< endl;
    Info << "Mean error Upwind: \t" << meanErrorUpwind*100.0 <<"[%]"<< endl;
    Info << "Mean error Linear: \t" << meanErrorLinear*100.0 <<"[%]"<< endl;
    Info << "***********************************************"<<endl;
}


