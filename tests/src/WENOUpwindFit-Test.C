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
#include <math.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TEST_CASE("WENOUpwindFit Test","[upwindFitTest]")
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
    // Create patch field list:
    const wordList patches = mesh.boundaryMesh().types();
    wordList patchTypes(patches.size());
    
    forAll(patches,patchI)
    {
        if (patches[patchI] == "patch")
            patchTypes[patchI] = "fixedValue";
        else if (patches[patchI] == "cyclic")
            patchTypes[patchI] = "cyclic";
        else if (patches[patchI] == "cyclicAMI")
            patchTypes[patchI] = "cyclicAMI";
        else if (patches[patchI] == "empty")
            patchTypes[patchI] = "empty";
        else
            patchTypes[patchI] = "processor";
    }

    // Boolean value to check if it is 2D 
    bool sim2D = mesh.nSolutionD() == 2 ? true : false;

    auto psiFunc = [&sim2D](const double x, const double y, const double z) -> double
    {
        //return std::sin(std::sqrt(2.0)/2.0*(x+y));
        if (sim2D)
            return std::sin(M_PI*x)+std::sin(M_PI*y);
        return std::sin(M_PI*x)+std::sin(M_PI*y)+std::sin(M_PI*z);
        
    };

    // Analytical sultion
    auto divPsiFunc = [&sim2D](const double x, const double y, const double z)
    {
        //return std::sqrt(2.0)*std::cos(std::sqrt(2.0)/2.0*(x+y));
        if (sim2D)
            return M_PI*std::cos(M_PI*x)+M_PI*std::cos(M_PI*y);
        return M_PI*std::cos(M_PI*x)+M_PI*std::cos(M_PI*y)+M_PI*std::cos(M_PI*z);
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
        psi[celli] = psiFunc(centre[celli].x(),centre[celli].y(),centre[celli].z());
    }
    
    {
        auto& PatchField = psi.boundaryFieldRef();
        forAll(PatchField,patchi)
        {
            const vectorField& faceCenters = PatchField[patchi].patch().Cf();
            forAll(PatchField[patchi],facei)
            {
                PatchField[patchi][facei] = psiFunc(faceCenters[facei].x(),faceCenters[facei].y(),faceCenters[facei].z());
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
        const double z = centre[celli].z();
        analSolu[celli] = divPsiFunc(x,y,z);
    }
    
    {
        auto& PatchField = analSolu.boundaryFieldRef();
        forAll(PatchField,patchi)
        {
            const vectorField& faceCenters = PatchField[patchi].patch().Cf();
            forAll(PatchField[patchi],facei)
            {
                PatchField[patchi][facei] = divPsiFunc
                (faceCenters[facei].x(),faceCenters[facei].y(),faceCenters[facei].z());
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
        U[celli] = vector(1,1,1);
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
        const scalar z = surfCenters[facei].z();
        psiSf[facei] = psiFunc(x,y,z);
        phi[facei] = Sf[facei] & vector(1,1,1);
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
            const scalar z = faceCenters[facei].z();
            psiSfPatchField[facei] = psiFunc(x,y,z);
            phiPatchField[facei] = Sf[facei] & vector(1,1,1);
        }
    }
    

    
            
    SECTION("Interpolate Function")
    {
        // ---------------------------------------------------------------------
        //                     Run Test Interpolate
        // ---------------------------------------------------------------------
        // Test the surface interpolation scheme
        surfaceScalarField interpWENO("interpWENO",fvc::interpolate(psi,phi,"interpolate(psiWENO)"));
        surfaceScalarField interpLinear = fvc::interpolate(psi,phi,"interpolate(psiLinear)");


        dimensionedScalar dimSmall("dimSmall",dimless,SMALL);
        surfaceScalarField errorWENO("errorWENO",interpWENO-psiSf);
        surfaceScalarField errorLinear("errorLinear",interpLinear-psiSf);
        
        errorWENO.write();
        errorLinear.write();
        
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
        
        // Check boundary fields
        forAll(psiSf.boundaryField(),patchi)
        {
            auto& psiSfPatchField = psiSf.boundaryField()[patchi];
            auto& errorWENOPatchField = errorWENO.boundaryField()[patchi];
            auto& errorLinearPatchField = errorLinear.boundaryField()[patchi];
            //auto& WENOPatchField = interpWENO.boundaryField()[patchi];
            //auto& linearPatchField = interpLinear.boundaryField()[patchi];
            forAll(psiSfPatchField,facei)
            {
                meanErrorWENO += std::fabs(errorWENOPatchField[facei]);
                meanErrorLinear += std::fabs(errorLinearPatchField[facei]);
                if (std::fabs(errorWENOPatchField[facei]) > maxErrorWENO)
                    maxErrorWENO = std::fabs(errorWENOPatchField[facei]);
                if (std::fabs(errorLinearPatchField[facei]) > maxErrorLinear)
                    maxErrorLinear = std::fabs(errorLinearPatchField[facei]);
            }
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

