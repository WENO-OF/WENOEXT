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
    cellDist_
    (
        [](const fvMesh& mesh)
        {
            if (Pstream::parRun())
            {
                IOobject cellDistLocalIO
                (
                    "cellDist",
                    (mesh.time().findClosestTime(mesh.time().startTime().value())).name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                );
                
                if (!cellDistLocalIO.typeHeaderOk<volScalarField>())
                    FatalError << "field 'cellDist' not found. "
                               << "Did you decompose with: decomposePar -cellDist? " 
                               << endl;
                
                volScalarField cellDistLocal
                (
                    cellDistLocalIO,
                    mesh
                );
                
                
                scalarField combinedField(cellDistLocal);
                // Check if parallel 
                if (Pstream::parRun())
                {
                    List<scalarField> allValues(Pstream::nProcs());

                    allValues[Pstream::myProcNo()] = combinedField;

                    Pstream::gatherList(allValues);
                    
                    Pstream::scatterList(allValues);

                    combinedField =
                        ListListOps::combine<scalarField>
                        (
                            allValues,
                            accessOp<scalarField>()
                        );
                }
                return combinedField;
                
            }
            
            return scalarField(1,1);
        }(mesh)
    )
{
    // List of processor with a count of cells already found
    labelList procCount(Pstream::nProcs(),0);
    
    forAll(cellDist_,cellI)
    {
        mapGlobalCellIDToLocal_.insert(std::pair<int,int>(cellI,procCount[cellDist_[cellI]]++));
    } 
    
    // Fill cellID list
    cellID_.resize(localMesh_.nCells(),-1);
    if (Pstream::parRun())
    {
        // Loop over cellDist and check if it matches with your processorID
        label j = 0;
        forAll(cellDist_,celli)
        {
            if (int(cellDist_[celli]) == Pstream::myProcNo())
            {
                cellID_[j] = celli;
                j++;
            }
        }
    }
    else
    {
        forAll(cellID_,cellI)
        {
            cellID_[cellI] = cellI;
        }
    }
}


const Foam::labelList& Foam::WENO::globalfvMesh::cellID() const
{
    return cellID_;
}

bool Foam::WENO::globalfvMesh::isLocalCell(const int globalCellID) const
{
    if (int(cellDist_[globalCellID]) == Pstream::myProcNo())
        return true;
    else
        return false;
}


int Foam::WENO::globalfvMesh::getProcID(const int globalCellID)
{
    return cellDist_[globalCellID];
}


int Foam::WENO::globalfvMesh::processorCellID(const int globalCellID)
{
    auto it = mapGlobalCellIDToLocal_.find(globalCellID);
    if (it != mapGlobalCellIDToLocal_.end())
        return it->second;
    else
        FatalErrorInFunction
        << "Global cell ID cannot be mapped to a local processor cell ID"
        << exit(FatalError);
    return -1;
}
