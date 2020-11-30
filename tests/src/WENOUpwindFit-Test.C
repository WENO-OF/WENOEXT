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
#include "WENOBase.H"
#include "fvCFD.H"
#include "writeToFile.H"

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
    
    
    // Mesh has two patches, outlet with fixed value and empty direction
    wordList patchTypes(2);
    patchTypes[mesh.boundary()["outlet"].index()] = "fixedValue"; 
    patchTypes[mesh.boundary()["empty"].index()] = "empty";
    

    auto psiFunc = [](const double x, const double y) -> double
    {
        //return std::sin(std::sqrt(2.0)/2.0*(x+y));
        return std::sin(x);
    };

    // Analytical sultion
    auto divPsiFunc = [](const double x, const double y)
    {
        //return std::sqrt(2.0)*std::cos(std::sqrt(2.0)/2.0*(x+y));
        return std::cos(x);
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
        "fixedValue"
    );

    forAll(mesh.C(),celli)
    {
        const double x = centre[celli].x();
        const double y = centre[celli].y();
        analSolu[celli] = divPsiFunc(x,y);
    }
    
    {
        auto& PatchField = analSolu.boundaryFieldRef();
        forAll(PatchField,patchi)
        {
            const vectorField& faceCenters = PatchField[patchi].patch().Cf();
            forAll(PatchField[patchi],facei)
            {
                PatchField[patchi][facei] = divPsiFunc(faceCenters[facei].x(),faceCenters[facei].y());
            }
        }
    }
    
        
    analSolu.write();
    
    
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
    surfaceScalarField psiSf = fvc::interpolate(psi,phi,"interpolate(psiLinear)");
    
    // Correct surface fields with analytical solution
    auto surfCenters = psiSf.mesh().Cf();
    auto surfCentersPhi = phi.mesh().Cf();
    auto Sf = phi.mesh().Sf();
    forAll(psiSf,facei)
    {
        const scalar x = surfCenters[facei].x();
        const scalar y = surfCenters[facei].y();
        psiSf[facei] = psiFunc(x,y);
        phi[facei] = Sf[facei] & vector(1,1,0);
    }
    
    // Correct boundary
    forAll(phi.boundaryField(),patchi)
    {
        auto& psiSfPatchField = psiSf.boundaryFieldRef()[patchi];
        auto& phiPatchField = phi.boundaryFieldRef()[patchi];
        const vectorField& faceCenters = phiPatchField.patch().Cf();
        auto Sf = phiPatchField.patch().Sf();
        forAll(phiPatchField,facei)
        {
            const scalar x = faceCenters[facei].x();
            const scalar y = faceCenters[facei].y();
            psiSfPatchField[facei] = psiFunc(x,y);
            phiPatchField[facei] = Sf[facei] & vector(1,1,0);
        }
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
        surfaceScalarField errorWENO = (interpWENO-psiSf);
        surfaceScalarField errorLinear = (interpLinear-psiSf);
        
        double meanErrorWENO = 0;
        double meanErrorLinear = 0;
        double maxErrorWENO = 0;
        double maxErrorLinear = 0;
        
        for (int facei = 0; facei < errorWENO.size(); facei++)
        {
            meanErrorWENO += std::fabs(errorWENO[facei]);
            meanErrorLinear += std::fabs(errorLinear[facei]);
            if (std::fabs(errorWENO[facei]) > maxErrorWENO)
                maxErrorWENO = std::fabs(errorWENO[facei]);
            if (std::fabs(errorLinear[facei]) > maxErrorLinear)
                maxErrorLinear = std::fabs(errorLinear[facei]);
        }
        meanErrorWENO /= (errorWENO.size());
        meanErrorLinear /= (errorLinear.size());
        
        Info << "---------------------------\n"
             << "       Interpolate         \n"
             << "---------------------------\n"
             << "Mean Error Linear: "<<meanErrorLinear<<nl
             << "Mean Error WENO:   "<<meanErrorWENO<<nl<<nl
             << "Max Error Linear:  "<<maxErrorLinear<<nl
             << "Max Error WENO:    "<<maxErrorWENO << nl
             << "---------------------------" << endl; 
        
        CHECK(meanErrorLinear > meanErrorWENO);
        CHECK(maxErrorLinear > maxErrorWENO);
    }
    
    SECTION("Divergence Function")
    {
        // ---------------------------------------------------------------------
        //                     Run Test Divergence
        // ---------------------------------------------------------------------
        // Test the surface interpolation scheme
        volScalarField divWENO("divWENO",fvc::div(phi,psi,"div(WENO)"));
        volScalarField divLinear("divLinear",fvc::div(phi,psi,"div(Linear)"));
        volScalarField divLimitedLinear("divLimitedLinear",fvc::div(phi,psi,"div(LimitedLinear)"));

        divWENO.write();
        divLinear.write();
        divLimitedLinear.write();
        
        volScalarField errorWENO = (divWENO-analSolu);
        volScalarField errorLinear = (divLinear-analSolu);
        volScalarField errorLimitedLinear = (divLimitedLinear-analSolu);
        
        double meanErrorWENO = 0;
        double meanErrorLinear = 0;
        double meanErrorLimitedLinear = 0;
        double maxErrorWENO = 0;
        double maxErrorLinear = 0;
        double maxErrorLimitedLinear = 0;
        

        forAll(errorWENO,celli)
        {
            meanErrorWENO += std::fabs(errorWENO[celli]);
            meanErrorLinear += std::fabs(errorLinear[celli]);
            meanErrorLimitedLinear += std::fabs(errorLimitedLinear[celli]);
            if (std::fabs(errorWENO[celli]) > maxErrorWENO)
                maxErrorWENO = std::fabs(errorWENO[celli]);
            if (std::fabs(errorLinear[celli]) > maxErrorLinear)
                maxErrorLinear = std::fabs(errorLinear[celli]);
            if (std::fabs(errorLimitedLinear[celli]) > maxErrorLimitedLinear)
                maxErrorLimitedLinear = std::fabs(errorLimitedLinear[celli]);
        }
        meanErrorWENO /= errorWENO.size();
        meanErrorLinear /= errorLinear.size();
        meanErrorLimitedLinear /= errorLimitedLinear.size();
        
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
        
        CHECK(meanErrorLinear > meanErrorWENO);
        CHECK(maxErrorLinear > maxErrorWENO);
        
        // Write Results to file
        std::vector<double> data = 
        {
            meanErrorLinear,meanErrorLimitedLinear,meanErrorWENO,
            maxErrorLinear,maxErrorLimitedLinear,maxErrorWENO
        };
        writeToFile(data);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

