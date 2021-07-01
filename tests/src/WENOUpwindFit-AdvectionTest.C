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
    WENOUpwindFit advection test

Description
    Test the WENOUpwindFit scheme by using the test case,
    the rotation of a slotted disk, designed by Zalesak
    
    [1] Zalesak, S.T. Fully Multidimensional Flux-Corrected Algorithms
        for Fluids. J. Comput. Phys. 1979, 31, 335–362.
    
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#include "catch.hpp"
#include "WENOBase.H"
#include "fvCFD.H"
#include "EulerDdtScheme.H"
#include "backwardDdtScheme.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TEST_CASE("WENOUpwindFit 2D Advection Test","[Advection]")
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

    // -------------------------------------------------------------------------
    //                          Create Input Data 
    // -------------------------------------------------------------------------
    
    
    // Mesh has two patches, outlet with fixed value and empty direction
    wordList patchTypes(2);
    patchTypes[0] = "fixedValue"; 
    patchTypes[1] = "empty";
    
    // Create a volScalarField to be transported
    volScalarField psiWENO
    (
        IOobject
        (
            "psiWENO",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0),
        patchTypes
    );

    // create the disk 
    const vectorField& centre = mesh.C();
    forAll(mesh.C(),celli)
    {
        const scalar R = 0.3;   // Radius of the disk
        const vector diskOrigin(0,0.5,0);
        const scalar x = centre[celli].x();
        const scalar y = centre[celli].y();
        const scalar dist = mag(vector(x,y,0)-diskOrigin);
        if 
        (
            dist < R
         && (fabs(x) > 0.05 || y > 0.7)
        )
        psiWENO[celli] = 1.0;
    }
    
    volScalarField psiLimitedLinear
    (
        IOobject
        (
            "psiLimitedLinear",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        psiWENO
    );

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
        psiWENO
    );

    psi.write();
    psiWENO.write();
    psiLimitedLinear.write();
    
    

    
    
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
    const double PI = Foam::constant::mathematical::pi;
    forAll(mesh.C(),celli)
    {
        const double x = centre[celli].x();
        const double y = centre[celli].y();
        U[celli] = vector(-2.0*PI*y,2.0*PI*x,0);
    }
    
    U.write();
    
    // Get the surface fields
    surfaceScalarField phi = fvc::flux(U);
    
    auto surfCenters = phi.mesh().Cf();
    auto Sf = phi.mesh().Sf();
    forAll(phi,facei)
    {
        const scalar x = surfCenters[facei].x();
        const scalar y = surfCenters[facei].y();
        phi[facei] = Sf[facei] & vector(-2.0*PI*y,2.0*PI*x,0);
    }
    
    // Correct boundary
    forAll(phi.boundaryField(),patchi)
    {
        auto& phiPatchField = phi.boundaryFieldRef()[patchi];
        const vectorField& faceCenters = phiPatchField.patch().Cf();
        auto Sf = phiPatchField.patch().Sf();
        forAll(phiPatchField,facei)
        {
            const scalar x = faceCenters[facei].x();
            const scalar y = faceCenters[facei].y();
            phiPatchField[facei] = Sf[facei] & vector(-2.0*PI*y,2.0*PI*x,0);
        }
    }

    // ---------------------------------------------------------------------
    //                      Run Advection Test 
    // ---------------------------------------------------------------------
    
    SECTION("Euler Time Discretization")
    {
    
        // Set time to zero
        runTime.setTime(0,0);
        // Set end time 
        runTime.setEndTime(1.0);
        runTime.setDeltaT(2.0e-04);       // Co < 0.3 for 300 cells
        
        while (runTime.run())
        {
            runTime++;
            Info<< "Time = " << runTime.timeName()<<endl;
            solve(fv::EulerDdtScheme<scalar>(mesh).fvmDdt(psiWENO) + fvm::div(phi,psiWENO,"div(WENO)"));
            solve(fv::EulerDdtScheme<scalar>(mesh).fvmDdt(psiLimitedLinear) + fvm::div(phi,psiLimitedLinear,"div(LimitedLinear)"));
        }
        psiWENO.write();
        psiLimitedLinear.write();
        psi.write();
        
        // accepted tolerance
        const double tol = 1e-3;
        
        // Check that WENO scheme is in bounds
        INFO("Check with Euler scheme and tolerance "<<tol<<" failed");
        forAll(mesh.C(),celli)
        {
            REQUIRE(psiWENO[celli] < (1.0+tol));
            REQUIRE(psiWENO[celli] > (0.0-tol));
        }
    }
    SECTION("Backward Time Discretization")
    {
    
        // Set time to zero
        runTime.setTime(0,0);
        // Set end time 
        runTime.setEndTime(1.0);
        runTime.setDeltaT(2.0e-04);       // Co < 0.3 for 300 cells
        
        while (runTime.run())
        {
            runTime++;
            Info<< "Time = " << runTime.timeName()<<endl;
            solve(fv::backwardDdtScheme<scalar>(mesh).fvmDdt(psiWENO) + fvm::div(phi,psiWENO,"div(WENO)"));
            solve(fv::backwardDdtScheme<scalar>(mesh).fvmDdt(psiLimitedLinear) + fvm::div(phi,psiLimitedLinear,"div(LimitedLinear)"));
        }
        psiWENO.write();
        psiLimitedLinear.write();
        psi.write();
        
        // accepted tolerance
        const double tol = 1e-2;
        
        INFO("Check with backward scheme and tolerance "<<tol<<" failed");
        // Check that WENO scheme is in bounds
        forAll(mesh.C(),celli)
        {
            REQUIRE(psiWENO[celli] < (1.0+tol));
            REQUIRE(psiWENO[celli] > (0.0-tol));
        }
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

