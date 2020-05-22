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
#include "reconstructRegionalMesh.H"

Foam::WENO::globalfvMesh::globalfvMesh(const fvMesh& mesh)
:
    neighborProcessor_
    (
        [](const fvMesh& mesh)
        {
            /******************************************************************\
            Neighbour Processor List:
                Contains the neighbouring processors and its neighbours. 
                With this all processors around the central processor are found.
                
                            *------*------*------*
                            |  CD  |  D   |  DA  |
                            |      |      |      |
                            *------*------*------*
                            |  C   | Main |  A   |
                            |      |      |      |
                            *------*------*------*
                            |  CB  |  B   |  AB  |
                            |      |      |      |
                            *------*------*------*
                
                The processors A,B,C and D are direct neighbours. To get the 
                processors CD, DA, AB and CB the neighbour of the neighbours 
                are collected. Here called second neighbours.
                
                **IMPORTANT** Reconstructing the mesh reads the processor 
                direcotries and creates an MPI blocking when the lists are not
                equally sized for all processors!
            \******************************************************************/
            
            
            if (Pstream::parRun())
            {
                const fvPatchList& patches = mesh.boundary();
                
                labelList myNeighbourProc;
                
                // Keep track of processor IDs already added to avoid duplicate
                // entries
                std::set<label> addedProcessor;
            
                // Add your own processor ID
                myNeighbourProc.append(Pstream::myProcNo());
                addedProcessor.insert(Pstream::myProcNo());
                
                forAll(patches, patchI)
                {
                    if(isA<processorFvPatch>(patches[patchI]))
                    {
                        label procID = 
                            (
                                refCast<const processorFvPatch>
                                (patches[patchI]).neighbProcNo()
                            );
                        if (addedProcessor.find(procID) == addedProcessor.end())
                            myNeighbourProc.append(procID);
                            
                        addedProcessor.insert(procID);
                    }
                }
                
                // distribute the list
                List<labelList> allNeighbours(Pstream::nProcs());

                allNeighbours[Pstream::myProcNo()] = myNeighbourProc;

                Pstream::gatherList(allNeighbours);
                
                Pstream::scatterList(allNeighbours);

                // Get the neighbour and second neighbour list for your processor
                labelList secondNeighborList;
                forAll(myNeighbourProc,procI)
                {
                    const label neighbourProcI = myNeighbourProc[procI];
                    secondNeighborList.append(allNeighbours[neighbourProcI]);
                }
                
                forAll(secondNeighborList,i)
                {
                    if ( addedProcessor.find(secondNeighborList[i]) == addedProcessor.end())
                    {
                        myNeighbourProc.append(secondNeighborList[i]);
                        addedProcessor.insert(secondNeighborList[i]);
                    }
                }
                
                // Get the maximum of myNeighbourProc cells and make everyone equal
                List<labelList> allValues(Pstream::nProcs());
                allValues[Pstream::myProcNo()] = myNeighbourProc;

                Pstream::gatherList(allValues);
                
                Pstream::scatterList(allValues);
                
                label maxProcNumber = 0;
                forAll(allValues,i)
                {
                    if (allValues[i].size()>maxProcNumber)
                        maxProcNumber = allValues[i].size();
                }
                
                while (myNeighbourProc.size() < maxProcNumber)
                {
                    // Get the neighbour and second neighbour list for your processor
                    secondNeighborList.clear();
                    forAll(myNeighbourProc,procI)
                    {
                        const label neighbourProcI = myNeighbourProc[procI];
                        secondNeighborList.append(allNeighbours[neighbourProcI]);
                    }
                    
                    forAll(secondNeighborList,i)
                    {
                        if ( addedProcessor.find(secondNeighborList[i]) == addedProcessor.end())
                        {
                            myNeighbourProc.append(secondNeighborList[i]);
                            addedProcessor.insert(secondNeighborList[i]);
                        }
                    }
                }
                
                // resize list
                myNeighbourProc.resize(maxProcNumber);

                #ifdef FULLDEBUG
                    Info << "Number of neighbour processor: "<<maxProcNumber<<endl;
                #endif
                return myNeighbourProc;
            }
            else
            {
                return labelList(0);
            }
        }(mesh)
    ),
    sendToProcessor_
    (
        [this](const fvMesh& mesh) -> labelList
        {
            if (Pstream::parRun())
            {
                labelList sendToProcessor;
                
                // distribute the list
                List<labelList> allValues(Pstream::nProcs());

                allValues[Pstream::myProcNo()] = neighborProcessor_;

                Pstream::gatherList(allValues);
                
                Pstream::scatterList(allValues);

                forAll(allValues,procI)
                {
                    if (procI != Pstream::myProcNo())
                    {
                        forAll(allValues[procI],i)
                        {
                            if (allValues[procI][i] == Pstream::myProcNo())
                            {
                                sendToProcessor.append(procI);
                                break;
                            }
                        }
                    }
                }
                
                stableSort(sendToProcessor);
               
                return sendToProcessor;
            }
            else
            {
                return labelList(0);
            }
        }(mesh)
    ),
    globalMeshPtr_
    (
        [this](const fvMesh& mesh) -> autoPtr<fvMesh>
        {
            if (Pstream::parRun())
            {
                return reconstructRegionalMesh::reconstruct(neighborProcessor_,mesh);
            }
            return autoPtr<fvMesh>(nullptr);
        }(mesh) 
    ),
    globalMesh_
    (   
        globalMeshPtr_.valid() ? globalMeshPtr_() : mesh
    ),
    localMesh_(mesh),
    procList_(),
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
            
            procList_.setSize(globalMesh_.nCells(),-1);

            // Check if parallel 
            if (Pstream::parRun())
            {
                       
                PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
              
                // Distribute the cell centers of the local mesh
                forAll(sendToProcessor_, procI)
                {
                    // Have to send tmpList because vectorField cannot use 
                    // overload of operator<<() for reading the list from buffer
                    UOPstream toBuffer(sendToProcessor_[procI], pBufs);
                    List<vector> tmpList(localMesh_.C());
                    toBuffer << tmpList;
                }
                
                pBufs.finishedSends();
                
                // Neighbour processor list contains the own processor itself.
                List<vectorField> cellCentersList(neighborProcessor_.size());
                
                forAll(neighborProcessor_, procI)
                {
                    if (neighborProcessor_[procI] != Pstream::myProcNo())
                    {
                        List<vector> tmpList;
                        UIPstream fromBuffer(neighborProcessor_[procI], pBufs);
                        fromBuffer >> tmpList;
                        cellCentersList[procI] = vectorField(tmpList);
                    }
                    else
                    {
                        cellCentersList[procI] = localMesh_.C();
                    }
                }
    
                const vectorField& globalCellCenters = globalMesh_.C();
                    
                const int nCells = globalMesh_.nCells();
                
                // Loop over all global cells and find corresponding cell center
                forAll(cellCentersList,procI)
                {
                    // user defined tolerance 
                    const double tol = 1E-15;
                    
                    int startCell = 0;
                    
                    const vectorField& localCellCentersI = cellCentersList[procI];
                    
                    forAll(localCellCentersI,cellI)
                    {
                        // Search through global mesh
                        for (int j=0; j < nCells; j++)
                        {
                            // Loop forward
                            int globalCellIndex =  
                                (startCell + j) >= nCells ?  (startCell + j - nCells) : (startCell + j);
                                
                            if (mag(localCellCentersI[cellI] - globalCellCenters[globalCellIndex]) < tol)
                            {
                                globalToLocalCellID[globalCellIndex] = cellI;
                                procList_[globalCellIndex] = neighborProcessor_[procI];
                                startCell = globalCellIndex;
                                break;
                            }
                            
                            // loop backwards
                            globalCellIndex =  
                                (startCell - j) < 0 ? (j - startCell) : (startCell - j);
                                
                            if (mag(localCellCentersI[cellI] - globalCellCenters[globalCellIndex]) < tol)
                            {
                                globalToLocalCellID[globalCellIndex] = cellI;
                                procList_[globalCellIndex] = neighborProcessor_[procI];
                                startCell = globalCellIndex;
                                break;
                            }

                        }
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
    )
{
    #ifdef FULLDEBUG
    Pout << "Reconstructed Mesh for processor "<<Pstream::myProcNo() << endl;
    #endif
}



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
