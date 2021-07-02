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

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020.

\*---------------------------------------------------------------------------*/


#include "globalfvMesh.H"
#include "reconstructRegionalMesh.H"
#include "meshSearch.H"

Foam::WENO::globalfvMesh::globalfvMesh(const fvMesh& mesh)
:
    neighborProcessor_(neighborProcessorList(mesh)),
    sendToProcessor_(sendToProcessorList(mesh)),
    globalMeshPtr_(createGlobalMesh(mesh)),
    localMeshPtr_(createLocalMesh(mesh)),
    globalMesh_
    (   
        globalMeshPtr_.valid() ? globalMeshPtr_() : mesh
    ),
    localMesh_
    (   
        localMeshPtr_.valid() ? localMeshPtr_() : mesh
    ),
    procList_(),
    localToGlobalCellID_(localToGlobalCellIDList()),
    globalToLocalCellID_(globalToLocalCellIDList())
{
    #ifdef FULLDEBUG
    Pout << "Reconstructed Mesh for processor "<<Pstream::myProcNo() << endl;
    #endif
}


Foam::labelList Foam::WENO::globalfvMesh::neighborProcessorList(const fvMesh& mesh)
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
        labelListList secondNeighbourList(myNeighbourProc.size()-1);
        
        // Jump first entry as it is its own processor
        for (label procI = 1; procI < myNeighbourProc.size(); ++procI)
        {
            const label neighbourProcI = myNeighbourProc[procI];
            secondNeighbourList[procI-1] = allNeighbours[neighbourProcI];
        }
        
        forAll(secondNeighbourList,procI)
        {
            forAll(secondNeighbourList[procI],i)
            {
                if ( addedProcessor.find(secondNeighbourList[procI][i]) == addedProcessor.end())
                {
                    // Get processor list of this neighbour 
                    const labelList& procsSecNeighbour = allNeighbours[secondNeighbourList[procI][i]];

                    // Get list of direct neighbours of own processor
                    const labelList& directNeighbour = allNeighbours[Pstream::myProcNo()];

                    // count the number of common processor interfaces
                    label commonProcessorInterfaces = 0;
                    forAll(procsSecNeighbour,k)
                    {
                        forAll(directNeighbour,m)
                        {
                            if (procsSecNeighbour[k] == directNeighbour[m])
                                commonProcessorInterfaces++;
                        }
                    }
                    
                    // Only add to list if it connects with at least 1 other direct neighbour 
                    if (commonProcessorInterfaces > 1)
                    {
                        myNeighbourProc.append(secondNeighbourList[procI][i]);
                        addedProcessor.insert(secondNeighbourList[procI][i]);
                    }
                }
            }
        }
        
        #ifdef FULLDEBUG
            Pout << "Number of neighbour processor: "<<myNeighbourProc.size()<<endl;
        #endif
        return myNeighbourProc;
    }

    // If not parallel return label list 0
    return labelList(0);
}
        
        
Foam::labelList Foam::WENO::globalfvMesh::sendToProcessorList(const fvMesh& mesh)
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

    // If not parallel return label list 0
    return labelList(0);
}
        
Foam::autoPtr<Foam::fvMesh> Foam::WENO::globalfvMesh::createGlobalMesh(const fvMesh& mesh)
{
    if (Pstream::parRun())
    {
        return reconstructRegionalMesh::reconstruct(neighborProcessor_,mesh);
    }
    return autoPtr<fvMesh>(nullptr);
}

Foam::autoPtr<Foam::fvMesh> Foam::WENO::globalfvMesh::createLocalMesh(const fvMesh& mesh)
{
    if (mesh.topoChanging())
        FatalErrorInFunction << "Mesh cannot change topology if used with WENO scheme"
                             << exit(FatalError);
    if (mesh.moving())
    {
        // If the mesh is moving it needs to be constructed from the same time field the 
        // global mesh is constructed from. 
        if (Pstream::parRun())
        {
            labelList list(1);
            list[0] = Pstream::myProcNo();
            return reconstructRegionalMesh::reconstruct(list,mesh);
        }
    }

    return autoPtr<fvMesh>(nullptr);
}


        
Foam::labelList Foam::WENO::globalfvMesh::localToGlobalCellIDList()
{
    // Fill cellID list
    labelList localToGlobalCellID(localMesh_.nCells(),-1);
       
    if (Pstream::parRun())
    {
        const vectorField& localCellCenters  = localMesh_.C();

        meshSearch mSearch(globalMesh_);

        label oldCellIndex = 0;

        // Loop over all local cells and look for cell center
        forAll(localCellCenters,cellI)
        {
            // Find cell
            label globalCellIndex = mSearch.findCell(localCellCenters[cellI],oldCellIndex);
            if (globalCellIndex == -1)
                globalCellIndex = mSearch.findCell(localCellCenters[cellI]);
                
            oldCellIndex = globalCellIndex;
            
            if (globalCellIndex >= 0)
                localToGlobalCellID[cellI] = globalCellIndex;
            else
            {
                globalMesh_.write();
                FatalErrorInFunction << "Could not find local cell in global mesh" <<nl
                                     << "Local cellID: "<<cellI<<"  position: "<< localCellCenters[cellI] << nl
                                     << "Global cell Index: "<<globalCellIndex<<nl
                                     << "Writing out global mesh" 
                                     << exit(FatalError);
            }
        }
    }
    // If it is running in serial
    else
    {
        forAll(localToGlobalCellID,cellI)
        {
            localToGlobalCellID[cellI] = cellI;
        }
    }
    
    return localToGlobalCellID;
}
        
Foam::labelList Foam::WENO::globalfvMesh::globalToLocalCellIDList()
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
        
        meshSearch mSearch(globalMesh_);
        
        // Loop over all global cells and find corresponding cell center
        forAll(cellCentersList,procI)
         { 
            const vectorField& localCellCentersI = cellCentersList[procI];
            
            label oldCellIndex = 0;
            
            forAll(localCellCentersI,cellI)
            { 
                // Find cell
                label globalCellIndex = mSearch.findCell(localCellCentersI[cellI],oldCellIndex);
                if (globalCellIndex == -1)
                    globalCellIndex = mSearch.findCell(localCellCentersI[cellI]);
                    
                oldCellIndex = globalCellIndex;
                    
                if (globalCellIndex >= 0)
                {
                    globalToLocalCellID[globalCellIndex] = cellI;
                    procList_[globalCellIndex] = neighborProcessor_[procI];
                }
                else
                    FatalErrorInFunction << "Processor local point not found in global mesh" <<nl
                               << "For local cell "<<localCellCentersI[cellI]
                               << " of processor "<<procI
                               << " for global mesh of processor "
                               <<Pstream::myProcNo()<<nl
                               << exit(FatalError);
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
}


// ***************************************************************************
//                   Public Member Functions



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
