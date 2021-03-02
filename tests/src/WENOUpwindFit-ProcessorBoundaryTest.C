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

TEST_CASE("WENOUpwindFit Processor Boundary Test","[procBoundary]")
{
    // Replace setRootCase.H for Catch2   
    int argc = 2;
    char executable[] = "main";
    char parallel[] = "-parallel";
    
    char **argv = static_cast<char**>(malloc(argc*sizeof(char*)));
    argv[0] = executable;
    argv[1] = parallel;
    Foam::argList args(argc, argv);
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }
    #include "createTime.H"        // create the time object
    #include "createMesh.H"        // create the mesh object

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // -------------------------------------------------------------------------
    //                          Create Input Data 
    // -------------------------------------------------------------------------
    
    // Create patch field list:
    const wordList patches = mesh.boundaryMesh().types();
    wordList patchTypes(patches.size());
    forAll(patches,patchI)
    {
        if (patches[patchI] == "patch")
            patchTypes[patchI] = "fixedValue";
        else
            patchTypes[patchI] = "processor";
    }
    
    // -------------------------------------------------------------------------
    //                          Create Fields
    // -------------------------------------------------------------------------
    
    const vectorField& centre = mesh.C();
    
    // Create a volScalarField to be transported
    volScalarField psi
    (
        IOobject
        (
            "psi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 1.0),
        patchTypes
    );

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


    // -------------------------------------------------------------------------
    //                              Test Sections
    // -------------------------------------------------------------------------


    SECTION("Simple Test")
    {
        forAll(mesh.C(),celli)
        {
            const double x = centre[celli].x();
            const double y = centre[celli].y();
            const double z = centre[celli].z();
            psi[celli] = std::cos(x)+std::sin(y);
        }

        {
            auto& PatchField = psi.boundaryFieldRef();
            forAll(PatchField,patchi)
            {
                const vectorField& faceCenters = PatchField[patchi].patch().Cf();
                forAll(PatchField[patchi],facei)
                {
                    PatchField[patchi][facei] = std::cos(faceCenters[facei].x())+std::sin(faceCenters[facei].y());
                }
            }
        }

        psi.write();
        
        
        forAll(mesh.C(),celli)
        {
            U[celli] = vector(1,1,0);
        }
        
        U.write();
        
        // Get the surface fields
        surfaceScalarField phi("phi",fvc::flux(U));
        
        auto surfCenters = phi.mesh().Cf();
        auto Sf = phi.mesh().Sf();
        forAll(phi,facei)
        {
            phi[facei] = Sf[facei] & vector(1,1,0);
        }
        
        // Correct boundary
        forAll(phi.boundaryField(),patchi)
        {
            auto& phiPatchField = phi.boundaryFieldRef()[patchi];
            const vectorField& faceCenters = phiPatchField.patch().Cf();
            auto Sf = phiPatchField.patch().Sf();
            forAll(phiPatchField,facei)
            {
                phiPatchField[facei] = Sf[facei] & vector(1,1,0);
            }
        }

        // ---------------------------------------------------------------------
        //                      Run Processor Boundary Test
        // ---------------------------------------------------------------------
        
        volScalarField divPsiAnalytical
        (
            IOobject
            (
                "divPsiAnalytical",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("0", dimless, 1.0)
        );

        
        volScalarField invalidCellsField
        (
            IOobject
            (
                "invalidCells",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("0", dimless, 0.0)
        );
        
        forAll(divPsiAnalytical,celli)
        {
            const double x = centre[celli].x();
            const double y = centre[celli].y();
            const double z = centre[celli].z();
            divPsiAnalytical[celli] = -std::sin(x)+std::cos(y);
        }

        divPsiAnalytical.write();
        
        psi.correctBoundaryConditions();
        
        // Calculate divergence field
        volScalarField divPsi("divPsi",fvc::div(phi,psi));
        int invalidCells = 0;
        forAll(divPsi,i)
        {
            if (abs(divPsi[i]-divPsiAnalytical[i]) > 1E-6)
            {
                invalidCellsField[i] = 1;
                invalidCells++;
            }
            
        }
        invalidCellsField.write();
        divPsi.write();
        
        INFO(
            "Divergence is not calculated correctly in the parallel case\n"
         << "Field of invalid cells is written out and can be opened with paraview");
        REQUIRE(invalidCells == 0);
    }


    //forAll(mesh.C(),celli)
    //{
        //const double x = centre[celli].x();
        //const double y = centre[celli].y();
        //const double z = centre[celli].z();
        //psi[celli] = x*y*z;
    //}

    //{
        //auto& PatchField = psi.boundaryFieldRef();
        //forAll(PatchField,patchi)
        //{
            //const vectorField& faceCenters = PatchField[patchi].patch().Cf();
            //forAll(PatchField[patchi],facei)
            //{
                //PatchField[patchi][facei] = faceCenters[facei].x()*faceCenters[facei].y()*faceCenters[facei].z();
            //}
        //}
    //}

    //psi.write();
    


    

    
    //// get velocity field by taking the roation of the u = rot(w,x) 
    //const vector omega(1.0,1.0,1.0);   // rotational velocity
    //forAll(mesh.C(),celli)
    //{
        //const double x = centre[celli].x();
        //const double y = centre[celli].y();
        //const double z = centre[celli].z();
        //U[celli] = vector(omega.y()*z-omega.z()*y, omega.z()*x - omega.x()*z,omega.x()*y-omega.y()*x);
    //}
    
    //U.write();
    
    //// Get the surface fields
    //surfaceScalarField phi("phi",fvc::flux(U));
    
    //auto surfCenters = phi.mesh().Cf();
    //auto Sf = phi.mesh().Sf();
    //forAll(phi,facei)
    //{
        //const scalar x = surfCenters[facei].x();
        //const scalar y = surfCenters[facei].y();
        //const scalar z = surfCenters[facei].z();
        //phi[facei] = Sf[facei] & vector(omega.y()*z-omega.z()*y, omega.z()*x - omega.x()*z,omega.x()*y-omega.y()*x);
    //}
    
    //// Correct boundary
    //forAll(phi.boundaryField(),patchi)
    //{
        //auto& phiPatchField = phi.boundaryFieldRef()[patchi];
        //const vectorField& faceCenters = phiPatchField.patch().Cf();
        //auto Sf = phiPatchField.patch().Sf();
        //forAll(phiPatchField,facei)
        //{
            //const scalar x = faceCenters[facei].x();
            //const scalar y = faceCenters[facei].y();
            //const scalar z = faceCenters[facei].z();
            //phiPatchField[facei] = Sf[facei] & vector(omega.y()*z-omega.z()*y, omega.z()*x - omega.x()*z,omega.x()*y-omega.y()*x);
        //}
    //}

    //// ---------------------------------------------------------------------
    ////                      Run Processor Boundary Test
    //// ---------------------------------------------------------------------
    
    //// As the velocity field is divergence free the div(phi,psi) can be calculated
    //// with: div(phi,psi) = psi*div(phi) + phi*grad(psi)
    //// where the first term is zero. 
    //// Analytical solution to grad(psi)
    //volScalarField divPsi
    //(
        //IOobject
        //(
            //"divPsiAnalytical",
            //mesh.time().timeName(),
            //mesh,
            //IOobject::NO_READ,
            //IOobject::AUTO_WRITE
        //),
        //mesh,
        //dimensionedScalar("0", dimless, 1.0)
    //);

    
    //volScalarField invalidCellsField
    //(
        //IOobject
        //(
            //"invalidCells",
            //mesh.time().timeName(),
            //mesh,
            //IOobject::NO_READ,
            //IOobject::AUTO_WRITE
        //),
        //mesh,
        //dimensionedScalar("0", dimless, 0.0)
    //);
    
    //forAll(divPsi,celli)
    //{
        //const double x = centre[celli].x();
        //const double y = centre[celli].y();
        //const double z = centre[celli].z();
        //divPsi[celli] = U[celli] & vector(y*z,x*z,x*y);
    //}

    //divPsi.write();
    
    //psi.correctBoundaryConditions();
    
    //// Calculate divergence field
    //volScalarField field("divPsi",fvc::div(phi,psi));
    //int invalidCells = 0;
    //forAll(field,i)
    //{
        //if (abs(field[i]-divPsi[i]) > 1E-9)
        //{
            //invalidCellsField[i] = 1;
            //invalidCells++;
        //}
        
    //}
    //invalidCellsField.write();
    //field.write();
    
    //INFO(
        //"Divergence is not calculated correctly in the parallel case\n"
     //<< "Field of invalid cells is written out and can be opened with paraview");
    //REQUIRE(invalidCells == 0);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

