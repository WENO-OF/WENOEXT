/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020.

\*---------------------------------------------------------------------------*/


#include "globalfvMesh.H"

Foam::WENO::globalfvMesh::globalfvMesh(const fvMesh& mesh)
:
    globalMeshPtr_
    (
        [](const fvMesh& mesh) -> fvMesh*
        {
            if (Pstream::parRun())
            {
                return new fvMesh
                (
                    IOobject
                    (
                        fvMesh::defaultRegion,
                        mesh.time().rootPath()/mesh.time().globalCaseName(),
                        Time
                        (
                            Foam::Time::controlDictName,
                            mesh.time().rootPath(),
                            mesh.time().globalCaseName()
                        ),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                );
            }
            return nullptr;
        }(mesh) 
    ),
    globalMesh_
    (   
        globalMeshPtr_.valid() ? globalMeshPtr_() : mesh
    ),
    localMesh_(mesh),
    localToGlobalCellID_
    (
        [this]()
        {
            // Fill cellID list
            labelList localToGlobalCellID(localMesh_.nCells(),-1);
            
            // user defined tolerance 
            const double tol = 1E-15;
            
            int startCell = 0;
            
            const vectorField& localCellCenters  = localMesh_.C();
            
            const vectorField& globalCellCenters = globalMesh_.C();
            
            const int nCells = globalMesh_.nCells();
            
            // Loop over all local cells and look for cell center
            forAll(localCellCenters,cellI)
            {
                // Search through global mesh
                for (int j=0; j < nCells; j++)
                {
                    // Loop forward
                    int globalCellIndex =  
                        (startCell + j) >= nCells ?  (startCell + j - nCells) : (startCell + j);
                        
                    if (mag(localCellCenters[cellI] - globalCellCenters[globalCellIndex]) < tol)
                    {
                        localToGlobalCellID[cellI] = globalCellIndex;
                        startCell = globalCellIndex;
                        break;
                    }
                    
                    // loop backwards
                    globalCellIndex =  
                        (startCell - j) < 0 ? (j - startCell) : (startCell - j);
                        
                    if (mag(localCellCenters[cellI] - globalCellCenters[globalCellIndex]) < tol)
                    {
                        localToGlobalCellID[cellI] = globalCellIndex;
                        startCell = globalCellIndex;
                        break;
                    }

                }
            }
            
            return localToGlobalCellID;
        }()
    ),
    globalToLocalCellID_
    (
        [this]() -> labelList
        {
            labelList globalToLocalCellID(globalMesh_.nCells(),-1);
            
            // Check if parallel 
            if (Pstream::parRun())
            {
                List<labelList> allValues(Pstream::nProcs());

                allValues[Pstream::myProcNo()] = localToGlobalCellID_;

                Pstream::gatherList(allValues);
                
                Pstream::scatterList(allValues);

                // Loop over all lists
                forAll(allValues,procI)
                {
                    forAll(allValues[procI],cellI)
                    {
                        globalToLocalCellID[allValues[procI][cellI]] = cellI;
                    }
                }
                
            }
            // If it is running in serial
            else
            {
                forAll(globalToLocalCellID,cellI)
                {
                    globalToLocalCellID[cellI] = cellI;
                }
            }

            return globalToLocalCellID;
        }()
    ),
    procList_
    (
        [this]()
        {
            labelList procList(globalMesh_.nCells(),-1);
            
            // Check if parallel 
            if (Pstream::parRun())
            {
                List<labelList> allValues(Pstream::nProcs());

                allValues[Pstream::myProcNo()] = localToGlobalCellID_;

                Pstream::gatherList(allValues);
                
                Pstream::scatterList(allValues);

                // Loop over all lists
                forAll(allValues,procI)
                {
                    forAll(allValues[procI],cellI)
                    {
                        procList[allValues[procI][cellI]] = procI;
                    }
                }
            }
            // If it is running in serial
            else
            {
                forAll(procList,cellI)
                {
                    procList[cellI] = 0;
                }
            }
            
            return procList;
        }()
    )
{}



bool Foam::WENO::globalfvMesh::isLocalCell(const int globalCellID) const
{
    if (int(procList_[globalCellID]) == Pstream::myProcNo())
        return true;
    else
        return false;
}


int Foam::WENO::globalfvMesh::getProcID(const int globalCellID) const
{
    return procList_[globalCellID];
}


int Foam::WENO::globalfvMesh::processorCellID(const int globalCellID) const
{
    return globalToLocalCellID_[globalCellID];
}
