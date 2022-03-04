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
    globalFvMesh Test 

Description
    Test the if the global mesh and functions are implemented correctly for a 
    3D cube with 5.12 million cells and decomposed into 512 domains with 8x8x8 processor
    
    This means the first 10000 cells are of processor0 and the second of processor1
    and so on... 
    
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Catch::Approx is needed when floats are compared

#include "fvCFD.H"
#include "globalfvMesh.H"

#include "globalFoamArgs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



TEST_CASE("globalFvMesh Test","[parallel]")
{
    Foam::argList& args = getFoamArgs();
    #include "createTime.H"
    #include "createMesh.H"
    

        
        
    // Create global mesh
    WENO::globalfvMesh globalfvMesh(mesh);
    
    // Check the cellID 
    int myProc = Pstream::myProcNo();

    // Check some basic functions
    for(int localCellI=0;localCellI<mesh.nCells();localCellI++)
    {
        if (globalfvMesh.isLocalCell(localCellI) != true)
            FatalError << "isLocalCell() error "
                   << "localCellI: "<<localCellI <<" not found as local Cell "<< endl;
                   
        if (globalfvMesh.getProcID(localCellI) != myProc)
            FatalError << "getProcID() error "
                   << "localCellI: "<<localCellI<< " not found correct processor ID"<<endl;
        
        if (mag(globalfvMesh().C()[globalfvMesh.localToGlobalCellID()[localCellI]] - globalfvMesh.localMesh().C()[localCellI])>1E-15)
        FatalError << "Mesh location error "
                   << "localCellI: "<<localCellI <<"  globalCellI: "<<globalfvMesh.localToGlobalCellID()[localCellI]<< nl
                   << globalfvMesh().C()[globalfvMesh.localToGlobalCellID()[localCellI]]<<" != " <<globalfvMesh.localMesh().C()[localCellI]
                   <<exit(FatalError);
    }
    
    // -------------------------------------------------------------------------
    //      Check that the reconstructed mesh gives correct cell centers
    // -------------------------------------------------------------------------
    
    // Get all processor cell centers
    
      
    if (Pstream::parRun())
    {
        
        // distribute the list
        List<List<vector>> allCellCenters(Pstream::nProcs());

        vectorField localCellCenter = mesh.C();

        allCellCenters[Pstream::myProcNo()] = localCellCenter;

        Pstream::gatherList(allCellCenters);
        
        Pstream::scatterList(allCellCenters);

        // Loop over reconstructed global mesh and check cell centers
        
        for(int globalCellI=0;globalCellI < globalfvMesh().nCells();globalCellI++)
        {
            // Get processor ID
            int procID = globalfvMesh.getProcID(globalCellI);
            
            // get global cell center
            const point globalPoint = globalfvMesh().C()[globalCellI];
            
            // Get processor cell center of allCellCenter List
            const List<vector>& cellCenterList = allCellCenters[procID];
            
            // Check that it contains the searched point
            bool found = false;
            forAll(cellCenterList,celli)
            {
                if (mag(cellCenterList[celli]-globalPoint) < 1E-9)
                {
                    found = true;
                    break;
                }
            }
            
            INFO("Did not find correct cell center location:"<<nl
                 << "For point p: ["<<globalPoint.x()<<","<<globalPoint.y()
                 << ","<<globalPoint.z()<<"] in processor "<<procID);
            REQUIRE(found == true);            
        }
    }
}


