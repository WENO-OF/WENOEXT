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
    Test the WENOUpwindFit scheme by using the interpolation function 
    
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#include "catch.hpp"

#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TEST_CASE("WENOUpwindFit 2D Test","[2D]")
{
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

    // -------------------------------------------------------------------------
    //                          Create Input Data 
    // -------------------------------------------------------------------------
        
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
        "zeroGradient"
    );

    auto calcSinus = [](const double x, const double y) -> double
    {
        return std::sin(std::sqrt(x*x+y*y));
    };

    // populate psi
    forAll(mesh.C(),celli)
    {
        psi[celli] = calcSinus(centre[celli].x(),centre[celli].y());
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
        "zeroGradient"
    );
    
    forAll(mesh.C(),celli)
    {
        const double x = centre[celli].x();
        const double y = centre[celli].y();
        double mag = (x*x+y*y);
        U[celli] = vector(x/mag,y/mag,0);
    }
    
    U.write();
    
    surfaceScalarField phi = fvc::flux(U);
    
    // Surface field of psi
    surfaceScalarField psiSf = fvc::interpolate(psi,phi,"interpolate(psiLinear)");
    // Correct
    auto surfCenters = psiSf.mesh().Cf();
    forAll(psiSf,facei)
    {
        psiSf[facei] = calcSinus(surfCenters[facei].x(),surfCenters[facei].y());
    }
        
            
    SECTION("Interpolate Function")
    {
        // ---------------------------------------------------------------------
        //                     Run Test Interpolate
        // ---------------------------------------------------------------------
        // Test the surface interpolation scheme
        surfaceScalarField interpWENO = fvc::interpolate(psi,phi,"interpolate(psiWENO)");
        surfaceScalarField interpLinear = fvc::interpolate(psi,phi,"interpolate(psiLinear)");

        dimensionedScalar dimSmall("dimSmall",dimless,SMALL);
        surfaceScalarField errorWENO = (interpWENO-psiSf)/(psiSf+SMALL);
        surfaceScalarField errorLinear = (interpLinear-psiSf)/(psiSf+SMALL);
        
        double meanErrorWENO = 0;
        double meanErrorLinear = 0;
        double maxErrorWENO = 0;
        double maxErrorLinear = 0;
        
        // Calculate the mean of the result and exclude first and second entry
        for (int facei = 1; facei < (errorWENO.size()-1); facei++)
        {
            meanErrorWENO += std::abs(errorWENO[facei]);
            meanErrorLinear += std::abs(errorLinear[facei]);
            if (std::abs(errorWENO[facei]) > maxErrorWENO)
                maxErrorWENO = std::abs(errorWENO[facei]);
            if (std::abs(errorLinear[facei]) > maxErrorLinear)
                maxErrorLinear = std::abs(errorLinear[facei]);
        }
        meanErrorWENO /= (errorWENO.size()-2);
        meanErrorLinear /= (errorLinear.size()-2);
        
        Info << "---------------------------\n"
             << "       Interpolate         \n"
             << "---------------------------\n"
             << "Mean Error Linear: "<<meanErrorLinear<<nl
             << "Mean Error WENO:   "<<meanErrorWENO<<nl<<nl
             << "Max Error Linear:  "<<maxErrorLinear<<nl
             << "Max Error WENO:    "<<maxErrorWENO << nl
             << "---------------------------" << endl; 
        
        REQUIRE(meanErrorLinear > meanErrorWENO);
        CHECK_NOFAIL(maxErrorLinear > maxErrorWENO);
    }
    
    SECTION("Divergence Function")
    {
        // ---------------------------------------------------------------------
        //                     Run Test Divergence
        // ---------------------------------------------------------------------
        // Test the surface interpolation scheme
        volScalarField divWENO = fvc::div(phi,psi,"div(WENO)");
        volScalarField divLinear = fvc::div(phi,psi,"div(Linear)");
        volScalarField divLimitedLinear = fvc::div(phi,psi,"div(LimitedLinear)");

        // Analytical sultion
        auto analyticalSoultion = [](const double x, const double y)
        {
            auto L = std::sqrt(x*x+y*y);
            return std::cos(L)/(L);
        };
        
        volScalarField analSolu
        (
            IOobject
            (
                "analyticalSolution",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("0", dimensionSet(0,0,-1,0,0), 0.0),
            "zeroGradient"
        );

        forAll(mesh.C(),celli)
        {
            const double x = centre[celli].x();
            const double y = centre[celli].y();
            analSolu[celli] = analyticalSoultion(x,y);
        }
        analSolu.write();

        dimensionedScalar dimSmall("dimSmall",dimensionSet(0,0,-1,0,0),SMALL);
        volScalarField errorWENO = (divWENO-analSolu)/(analSolu+dimSmall);
        volScalarField errorLinear = (divLinear-analSolu)/(analSolu+dimSmall);
        volScalarField errorLimitedLinear = (divLimitedLinear-analSolu)/(analSolu+dimSmall);
        
        double meanErrorWENO = 0;
        double meanErrorLinear = 0;
        double meanErrorLimitedLinear = 0;
        double maxErrorWENO = 0;
        double maxErrorLinear = 0;
        double maxErrorLimitedLinear = 0;
        
        // Calculate the mean of the result and exclude first and second entry
        forAll(errorWENO,celli)
        {
            meanErrorWENO += std::abs(errorWENO[celli]);
            meanErrorLinear += std::abs(errorLinear[celli]);
            meanErrorLimitedLinear += std::abs(errorLimitedLinear[celli]);
            if (std::abs(errorWENO[celli]) > maxErrorWENO)
                maxErrorWENO = std::abs(errorWENO[celli]);
            if (std::abs(errorLinear[celli]) > maxErrorLinear)
                maxErrorLinear = std::abs(errorLinear[celli]);
            if (std::abs(errorLimitedLinear[celli]) > maxErrorLimitedLinear)
                maxErrorLimitedLinear = std::abs(errorLimitedLinear[celli]);
        }
        meanErrorWENO /= (errorWENO.size());
        meanErrorLinear /= (errorLinear.size());
        meanErrorLimitedLinear /= (errorLimitedLinear.size());
        
        Info << "---------------------------\n"
             << "       Divergence          \n"
             << "---------------------------\n"
             << "Mean Error Linear:        "<<meanErrorLinear<<nl
             << "Mean Error LimitedLinear: "<<meanErrorLimitedLinear<<nl
             << "Mean Error WENO:          "<<meanErrorWENO<<nl<<nl
             << "Max Error Linear:         "<<maxErrorLinear<<nl
             << "Max Error LimitedLinear:  "<<maxErrorLimitedLinear<<nl
             << "Max Error WENO:           "<<maxErrorWENO << nl
             << "---------------------------" << endl; 
        
        REQUIRE(meanErrorLinear > meanErrorWENO);
        CHECK_NOFAIL(maxErrorLinear > maxErrorWENO);
    }
}

//SECTION("Step Function")
    //{
     //// ---------------------------------------------------------------------
        ////                   Create Input Data 
        //// ---------------------------------------------------------------------
            
            //// Create a volScalarField with the sinus curve
            //volScalarField psi
            //(
                //IOobject
                //(
                    //"psi",
                    //mesh.time().timeName(),
                    //mesh,
                    //IOobject::NO_READ,
                    //IOobject::NO_WRITE
                //),
                //mesh,
                //dimensionedScalar("0", dimless, 0.0),
                //"zeroGradient"
            //);

            //auto calcStepFunction = [](const double x, const double y) -> double
            //{
                //return std::sin(std::sqrt(x*x+y*y)) > 0 ? 1.0 : -0.5;
            //};

            //// populate psi
            //forAll(mesh.C(),celli)
            //{
                //psi[celli] = calcStepFunction(centre[celli].x(),centre[celli].y());
            //}

            //psi.write();
            

            
            
            //// Create velocity field 
            //volVectorField U
            //(
                //IOobject
                //(
                    //"U",
                    //mesh.time().timeName(),
                    //mesh,
                    //IOobject::NO_READ,
                    //IOobject::NO_WRITE
                //),
                //mesh,
                //dimensionedVector("0", dimless, vector(0.0,0.0,0.0)),
                //"zeroGradient"
            //);
            
            //forAll(mesh.C(),celli)
            //{
                //const double x = centre[celli].x();
                //const double y = centre[celli].y();
                //double mag = (x*x+y*y);
                //U[celli] = vector(x/mag,y/mag,0);
            //}
            
            //U.write();
            
            //surfaceScalarField phi = fvc::flux(U);
            
            //// Surface field of psi
            //surfaceScalarField psiSf = fvc::interpolate(psi,phi,"interpolate(psiLinear)");
            //// Correct
            //auto surfCenters = psiSf.mesh().Cf();
            //forAll(psiSf,facei)
            //{
                //psiSf[facei] = calcStepFunction(surfCenters[facei].x(),surfCenters[facei].y());
            //}
        //// ---------------------------------------------------------------------
        ////                     Run Test 
        //// ---------------------------------------------------------------------
        //// Test the surface interpolation scheme
        //surfaceScalarField interpWENO = fvc::interpolate(psi,phi,"interpolate(psiWENO)");
        //surfaceScalarField interpLinear = fvc::interpolate(psi,phi,"interpolate(psiLinear)");

        //dimensionedScalar dimSmall("dimSmall",dimless,SMALL);
        //surfaceScalarField errorWENO = (interpWENO-psiSf)/(psiSf+SMALL);
        //surfaceScalarField errorLinear = (interpLinear-psiSf)/(psiSf+SMALL);
        
        //double meanErrorWENO = 0;
        //double meanErrorLinear = 0;
        //double maxErrorWENO = 0;
        //double maxErrorLinear = 0;
        
        //// Calculate the mean of the result and exclude first and second entry
        //for (int facei = 1; facei < (errorWENO.size()-1); facei++)
        //{
            //meanErrorWENO += std::abs(errorWENO[facei]);
            //meanErrorLinear += std::abs(errorLinear[facei]);
            //if (std::abs(errorWENO[facei]) > maxErrorWENO)
                //maxErrorWENO = std::abs(errorWENO[facei]);
            //if (std::abs(errorLinear[facei]) > maxErrorLinear)
                //maxErrorLinear = std::abs(errorLinear[facei]);
        //}
        //meanErrorWENO /= (errorWENO.size()-2);
        //meanErrorLinear /= (errorLinear.size()-2);
        
        //Info << "---------------------------\n"
             //<< "Mean Error Linear: "<<meanErrorLinear<<nl
             //<< "Mean Error WENO:   "<<meanErrorWENO<<nl<<nl
             //<< "Max Error Linear:  "<<maxErrorLinear<<nl
             //<< "Max Error WENO:    "<<maxErrorWENO<<endl;
        
        //REQUIRE(meanErrorLinear > meanErrorWENO);
        //CHECK(maxErrorLinear > maxErrorWENO);
    //}
