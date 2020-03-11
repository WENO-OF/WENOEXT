/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020

\*---------------------------------------------------------------------------*/

#include "reconstructRegionalMesh.H"

Foam::autoPtr<Foam::fvMesh> Foam::reconstructRegionalMesh::reconstruct
(
    const labelList processorList,
    const fvMesh& localMesh
)
{
    word regionName = polyMesh::defaultRegion;
    word regionDir = word::null;

    scalar mergeTol = 1E-7;
    
    label nProcs = processorList.size();
    
    // Read all time databases
    PtrList<Time> databases(nProcs);
    
    forAll(databases, proci)
    {
       // Pout << "proci: "<<proci<<endl;
        databases.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                localMesh.time().rootPath(),
                localMesh.time().globalCaseName()/fileName(word("processor") + name(processorList[proci]))
            )
        );
    }
    
    // Read point on individual processors to determine merge tolerance
    // (otherwise single cell domains might give problems)
    const boundBox bb = procBounds(databases, regionDir);
    const scalar mergeDist = mergeTol*bb.mag();

    //Info<< "Overall mesh bounding box  : " << bb << nl
        //<< "Relative tolerance         : " << mergeTol << nl
        //<< "Absolute matching distance : " << mergeDist << nl
        //<< endl;


    // Addressing from processor to reconstructed case
    labelListList cellProcAddressing(nProcs);
    labelListList faceProcAddressing(nProcs);
    labelListList pointProcAddressing(nProcs);
    labelListList boundaryProcAddressing(nProcs);

    List<fvMesh*> masterMesh(nProcs);

    for (label proci=0; proci<nProcs; proci++)
    {
        masterMesh[proci] = 
        (
            new fvMesh
            (
                IOobject
                (
                    regionName,
                    localMesh.time().timeName(),
                    localMesh.time(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                xferCopy(pointField()),
                xferCopy(faceList()),
                xferCopy(cellList())
            )
        );

        
        // Mesh cannot be constructed with boundaries as this calls the 
        // updateMesh() function of coupled processors and causes an MPI error
        // --> Solution construct without boundary
        
        pointField points = static_cast<pointField>
        (
            pointIOField
            (
                 IOobject
                 (
                    "points",
                    "constant/" + polyMesh::meshSubDir,
                    databases[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                 )
            )
         );

        faceList faces = static_cast<faceList>
        (
            faceCompactIOList
            (
                IOobject
                (
                    "faces",
                    "constant/" + polyMesh::meshSubDir,
                    databases[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
         );

        labelList owner = static_cast<labelList>
        (
            labelIOList
            (
                IOobject
                (
                    "owner",
                    "constant/" + polyMesh::meshSubDir,
                    databases[proci],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );
        
        labelList neighbour = static_cast<labelList>
        (
            labelIOList
            (
                IOobject
                (
                    "neighbour",
                    "constant/" + polyMesh::meshSubDir,
                    databases[proci],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
         );
            
        fvMesh meshToAdd
        (
            IOobject
            (
                regionName,
                databases[proci].timeName(),
                databases[proci],
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            xferMove(points),
            xferMove(faces),
            xferMove(owner),
            xferMove(neighbour)
        );

        // Now add the boundaries by creating a polyPatch with a type patch 
        // for generic value
        
        // Read polyPatchList
        const polyBoundaryMesh& polyMeshRef = meshToAdd.boundaryMesh();
         autoPtr<ISstream> isPtr = fileHandler().readStream
         (
            const_cast<polyBoundaryMesh&>(polyMeshRef),
            databases[proci].path()/fileName("constant/polyMesh/boundary"),
            "polyBoundaryMesh"
        );
 
        ISstream& is = isPtr();
 
         PtrList<entry> patchEntries(is);
         
         List<polyPatch*> patches(patchEntries.size());
 
         forAll(patches, patchi)
         {
             patches[patchi] = 
                (polyPatch::New
                (
                   "patch",
                    patchEntries[patchi].keyword(),
                    patchEntries[patchi].dict(),
                    patchi,
                    polyMeshRef
                )).ptr();
         }

        meshToAdd.addPatches(patches,false);

        // Initialize its addressing
        cellProcAddressing[proci] = identity(meshToAdd.nCells());
        faceProcAddressing[proci] = identity(meshToAdd.nFaces());
        pointProcAddressing[proci] = identity(meshToAdd.nPoints());
        boundaryProcAddressing[proci] =
            identity(meshToAdd.boundaryMesh().size());

        // Find geometrically shared points/faces.
        autoPtr<faceCoupleInfo> couples
        (
            new faceCoupleInfo
            (
                *masterMesh[proci],
                meshToAdd,
                mergeDist,      // Absolute merging distance
                true            // Matching faces identical
            )
        );
        
        //autoPtr<faceCoupleInfo> couples = determineCoupledFaces
        //(
            //fullMatch,
            //proci,
            //proci,
            //*masterMesh[proci],
            //proci,
            //proci,
            //meshToAdd,
            //mergeDist
        //);

        // Add elements to mesh
        autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
        (
            *masterMesh[proci],
            meshToAdd,
            couples
        );

        // Added processor
        renumber(map().addedCellMap(), cellProcAddressing[proci]);
        renumber(map().addedFaceMap(), faceProcAddressing[proci]);
        renumber(map().addedPointMap(), pointProcAddressing[proci]);
        renumber(map().addedPatchMap(), boundaryProcAddressing[proci]);
    }

    for (label step=2; step<nProcs*2; step*=2)
    {
        for (label proci=0; proci<nProcs; proci+=step)
        {
            label next = proci + step/2;
            if(next >= nProcs)
            {
                continue;
            }

            Info<< "Merging mesh " << proci << " with " << next << endl;

            // Find geometrically shared points/faces.
            autoPtr<faceCoupleInfo> couples
            (
                new faceCoupleInfo
                (
                    *masterMesh[proci],
                    *masterMesh[next],
                    mergeDist,      // Absolute merging distance
                    true            // Matching faces identical
                )
            );

            // Add elements to mesh
            autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
            (
                *masterMesh[proci],
                *masterMesh[next],
                couples
            );

            // Processors that were already in masterMesh
            for (label mergedI=proci; mergedI<next; mergedI++)
            {
                renumber
                (
                    map().oldCellMap(),
                    cellProcAddressing[mergedI]
                );

                renumber
                (
                    map().oldFaceMap(),
                    faceProcAddressing[mergedI]
                );

                renumber
                (
                    map().oldPointMap(),
                    pointProcAddressing[mergedI]
                );

                // Note: boundary is special since can contain -1.
                renumber
                (
                    map().oldPatchMap(),
                    boundaryProcAddressing[mergedI]
                );
            }

            // Added processor
            for
            (
                label addedI=next;
                addedI<min(proci+step, nProcs);
                addedI++
            )
            {
                renumber
                (
                    map().addedCellMap(),
                    cellProcAddressing[addedI]
                );

                renumber
                (
                    map().addedFaceMap(),
                    faceProcAddressing[addedI]
                );

                renumber
                (
                    map().addedPointMap(),
                    pointProcAddressing[addedI]
                );

                renumber
                (
                    map().addedPatchMap(),
                    boundaryProcAddressing[addedI]
                );
            }
        }
    }

    // See if any points on the mastermesh have become connected
    // because of connections through processor meshes.
    mergeSharedPoints(mergeDist, *masterMesh[0], pointProcAddressing);

    // delete all pointers except master
    for (int i=1;i<masterMesh.size();i++)
    {
        delete masterMesh[i];
    }

    fvMesh* meshPtr = masterMesh[0];

    return autoPtr<fvMesh>(meshPtr);
}



void Foam::reconstructRegionalMesh::renumber
(
    const labelList& map,
    labelList& elems
)
{
    forAll(elems, i)
    {
        if (elems[i] >= 0)
        {
            elems[i] = map[elems[i]];
        }
    }
}



Foam::autoPtr<Foam::mapPolyMesh> Foam::reconstructRegionalMesh::mergeSharedPoints
(
    const scalar mergeDist,
    polyMesh& mesh,
    labelListList& pointProcAddressing
)
{
    // Find out which sets of points get merged and create a map from
    // mesh point to unique point.
    Map<label> pointToMaster
    (
        fvMeshAdder::findSharedPoints
        (
            mesh,
            mergeDist
        )
    );

    Info<< "mergeSharedPoints : detected " << pointToMaster.size()
        << " points that are to be merged." << endl;

    if (returnReduce(pointToMaster.size(), sumOp<label>()) == 0)
    {
        return autoPtr<mapPolyMesh>(nullptr);
    }

    polyTopoChange meshMod(mesh);

    fvMeshAdder::mergePoints(mesh, pointToMaster, meshMod);

    // Change the mesh (no inflation). Note: parallel comms allowed.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

    // Update fields. No inflation, parallel sync.
    mesh.updateMesh(map);

    // pointProcAddressing give indices into the master mesh so adapt them
    // for changed point numbering.

    // Adapt constructMaps for merged points.
    forAll(pointProcAddressing, proci)
    {
        labelList& constructMap = pointProcAddressing[proci];

        forAll(constructMap, i)
        {
            label oldPointi = constructMap[i];

            // New label of point after changeMesh.
            label newPointi = map().reversePointMap()[oldPointi];

            if (newPointi < -1)
            {
                constructMap[i] = -newPointi-2;
            }
            else if (newPointi >= 0)
            {
                constructMap[i] = newPointi;
            }
            else
            {
                FatalErrorInFunction
                    << "Problem. oldPointi:" << oldPointi
                    << " newPointi:" << newPointi << abort(FatalError);
            }
        }
    }

    return map;
}


Foam::boundBox Foam::reconstructRegionalMesh::procBounds
(
    const PtrList<Time>& databases,
    const word& regionDir
)
{
    boundBox bb = boundBox::invertedBox;

    forAll(databases, proci)
    {
        pointIOField points
        (
             IOobject
             (
                "points",
                "constant/" + polyMesh::meshSubDir,
                databases[proci],
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
             )
        );

        boundBox domainBb(points, false);

        bb.min() = min(bb.min(), domainBb.min());
        bb.max() = max(bb.max(), domainBb.max());
    }

    return bb;
}
