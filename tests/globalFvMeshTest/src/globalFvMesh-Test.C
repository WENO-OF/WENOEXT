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
    3D cube with 50x50x50 cells and decomposed into 8 domains with 2x2x2 processor
    
    This means the first 125 cells are of processor0 and the second of processor1
    and so on... 
    
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "globalfvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    

        
        
    // Create global mesh
    WENO::globalfvMesh globalfvMesh(mesh);
    
    // Check the cellID 
    int myProc = Pstream::myProcNo();

    // Check local to global cellID
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
        
    Info << "END RUN 2"<<endl;
    return 0;
}


