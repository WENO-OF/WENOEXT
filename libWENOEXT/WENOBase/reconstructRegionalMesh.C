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
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020

\*---------------------------------------------------------------------------*/

#include "reconstructRegionalMesh.H"
#include "masterUncollatedFileOperation.H"



Foam::autoPtr<Foam::fvMesh> Foam::reconstructRegionalMesh::reconstruct
(
    const labelList& processorList,
    const labelList& sendToProcessorList,
    const fvMesh& localMesh
)
{
    // Get points, faces, neighbour, owner lists from processors 
    List<faceList>   allFaces(Pstream::nProcs());
    List<pointField> allPoints(Pstream::nProcs());
    List<labelList>  allOwners(Pstream::nProcs());
    List<labelList>  allNeighbours(Pstream::nProcs());


    getMeshInformation(allFaces,allPoints,allOwners,allNeighbours,
        processorList,sendToProcessorList,localMesh);


    List<wordList> allPatchTypes;
    List<labelList> allPatchSize;
    List<labelList> allStart;
    List<labelList> allIndex;
    
    getPatchInformation(allPatchTypes,allPatchSize,allStart,allIndex,
        processorList,sendToProcessorList,localMesh);

    // Get the dimensions of the current processor domain and calculate the 
    // merging tolerance
    word regionName = polyMesh::defaultRegion;
    word regionDir = word::null;

    scalar mergeTol = 1E-7;
    
    // Read point on individual processors to determine merge tolerance
    // (otherwise single cell domains might give problems)
    const boundBox bb = procBounds(processorList,allPoints);
    const scalar mergeDist = mergeTol*bb.mag();

    //Info<< "Overall mesh bounding box  : " << bb << nl
        //<< "Relative tolerance         : " << mergeTol << nl
        //<< "Absolute matching distance : " << mergeDist << nl
        //<< endl;
        
    fvMesh* masterMesh = 
        (
            new fvMesh
            (
                IOobject
                (
                    regionName,
                    localMesh.time().constant(),
                    localMesh.time(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                #ifdef FOAM_MOVE_CONSTRUCTOR
                    pointField(),
                    faceList(),
                    labelList(),
                    labelList()
                #else
                    xferCopy(pointField()),
                    xferCopy(faceList()),
                    xferCopy(labelList()),
                    xferCopy(labelList())
                #endif
            )
        );

    for (label procI : processorList)
    {        
        // Mesh cannot be constructed with boundaries as this calls the 
        // updateMesh() function of coupled processors and causes an MPI error
        // --> Solution construct without boundary

        fvMesh meshToAdd
        (
            IOobject
            (
                regionName,
                localMesh.time().constant(),
                localMesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            #ifdef FOAM_MOVE_CONSTRUCTOR
                std::move(allPoints[procI]),
                std::move(allFaces[procI]),
                std::move(allOwners[procI]),
                std::move(allNeighbours[procI]),
            #else
                xferMove(allPoints[procI]),
                xferMove(allFaces[procI]),
                xferMove(allOwners[procI]),
                xferMove(allNeighbours[procI]),
            #endif
            false // Do not synchronize --> in polyMesh() bounds_ calls an mpirecv 
        );

        // Now add the boundaries by creating a polyPatch with a type patch 
        // for generic value

        // Read polyPatchList
        const polyBoundaryMesh& polyMeshRef = meshToAdd.boundaryMesh();

        IFstream is(localPath(localMesh,procI,"boundary"));
        readHeader(is);
        PtrList<entry> patchEntries(is);
         
        List<polyPatch*> patches(patchEntries.size());
 
        forAll(patches, patchi)
        {
            patches[patchi] = 
                (
                    polyPatch::New
                    (
                        "patch",
                        patchEntries[patchi].keyword(),
                        patchEntries[patchi].dict(),
                        patchi,
                        polyMeshRef
                    )
                ).ptr();
        }

        meshToAdd.addPatches(patches,false);








        // // Access to boundary mesh
        // const polyBoundaryMesh& polyMeshRef = meshToAdd.boundaryMesh();

        // List<polyPatch*> patches(allPatchTypes[procI].size());
 
        // forAll(patches, patchi)
        // {
        //     Pout << "patchi: "<<patchi<<"  index "<<allIndex[procI][patchi]<<endl;
        //     patches[patchi] = 
        //         (
        //             polyPatch::New
        //             (
        //                 "patch",        // use generic patch instead of patchType
        //                 "patchName",
        //                 allPatchSize[procI][patchi],
        //                 allStart[procI][patchi],
        //                 allIndex[procI][patchi],
        //                 polyMeshRef
        //             )
        //         ).ptr();
        // }

        // meshToAdd.addPatches(patches,false);


        // Find geometrically shared points/faces.
        autoPtr<faceCoupleInfo> couples
        (
            new faceCoupleInfo
            (
                *masterMesh,
                meshToAdd,
                mergeDist,      // Absolute merging distance
                true            // Matching faces identical
            )
        );

        // Add elements to mesh
        autoPtr<mapAddedPolyMesh> map = add
        (
            *masterMesh,
            meshToAdd,
            couples,
            false // do not synchronize
        );
        
    }

    return autoPtr<fvMesh>(masterMesh);
}


Foam::autoPtr<Foam::mapAddedPolyMesh> Foam::reconstructRegionalMesh::add
(
    fvMesh& mesh0,
    const fvMesh& mesh1,
    const faceCoupleInfo& coupleInfo,
    const bool validBoundary
)
{
    // -----------------------------------
    // Function taken from fvMeshAddr::add
    // -----------------------------------
    
    mesh0.clearOut();
    
    // Resulting merged mesh (polyMesh only!)
    autoPtr<mapAddedPolyMesh> mapPtr
    (
        polyMeshAdder::add
        (
            mesh0,
            mesh1,
            coupleInfo,
            validBoundary
        )
    );
    
    return mapPtr;
}


Foam::boundBox Foam::reconstructRegionalMesh::procBounds
(
    const labelList processorList,
    const List<pointField>& allPoints
)
{
    boundBox bb = boundBox::invertedBox;

    forAll(processorList, procI)
    {
        boundBox domainBb(allPoints[processorList[procI]], false);

        bb.min() = min(bb.min(), domainBb.min());
        bb.max() = max(bb.max(), domainBb.max());
    }

    return bb;
}


void Foam::reconstructRegionalMesh::getMeshInformation
(
    List<faceList>& allFaces,
    List<pointField>& allPoints,
    List<labelList>& allOwner,
    List<labelList>& allNeighbours,
    const labelList& processorList,
    const labelList& sendToProcessorList,
    const fvMesh& mesh
)
{
    // Add points and faces for own processor
    allFaces[Pstream::myProcNo()] = mesh.faces();
    allPoints[Pstream::myProcNo()] = mesh.points();
    allOwner[Pstream::myProcNo()] = mesh.faceOwner();
    allNeighbours[Pstream::myProcNo()] = mesh.faceNeighbour();

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    for (label procI : sendToProcessorList)
    {
        // Check for own processor ID as sendToProcessor contains its own procI
        if (procI == Pstream::myProcNo())
            continue;

        UOPstream send(procI,pBufs);
        send << allFaces[Pstream::myProcNo()];
        send << allPoints[Pstream::myProcNo()];
        send << allOwner[Pstream::myProcNo()];
        send << allNeighbours[Pstream::myProcNo()];
    }

    pBufs.finishedSends();

    for (label procI : processorList)
    {
        if (procI == Pstream::myProcNo())
            continue;
        
        UIPstream recv(procI,pBufs);
        recv >> allFaces[procI];
        recv >> allPoints[procI];
        recv >> allOwner[procI];
        recv >> allNeighbours[procI];
    }
}



void Foam::reconstructRegionalMesh::getPatchInformation
(
    List<wordList>& allPatchTypes,
    List<labelList>& allPatchSize,
    List<labelList>& allStart,
    List<labelList>& allIndex,
    const labelList& processorList,
    const labelList& sendToProcessorList,
    const fvMesh& mesh  
)
{
    // Add points and faces for own processor
    allPatchTypes.resize(Pstream::nProcs());
    allPatchSize.resize(Pstream::nProcs());
    allStart.resize(Pstream::nProcs());
    allIndex.resize(Pstream::nProcs());

    // Get local boundary mesh
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    allPatchTypes[Pstream::myProcNo()].resize(bMesh.size());
    allPatchSize[Pstream::myProcNo()].resize(bMesh.size());
    allStart[Pstream::myProcNo()].resize(bMesh.size());
    allIndex[Pstream::myProcNo()].resize(bMesh.size());

    // Loop over all patches
    forAll(bMesh,i)
    {
        const polyPatch& patch = bMesh[i];
        allPatchTypes[Pstream::myProcNo()][i] = patch.type();
        allPatchSize[Pstream::myProcNo()][i] = patch.size();
        allStart[Pstream::myProcNo()][i] = patch.start();
        allIndex[Pstream::myProcNo()][i] = patch.index();
    }

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    for (label procI : sendToProcessorList)
    {
        // Check for own processor ID as sendToProcessor contains its own procI
        if (procI == Pstream::myProcNo())
            continue;

        UOPstream send(procI,pBufs);

        send << allPatchTypes[Pstream::myProcNo()];
        send << allPatchSize[Pstream::myProcNo()];
        send << allStart[Pstream::myProcNo()];
        send << allIndex[Pstream::myProcNo()];
    }

    pBufs.finishedSends();

    for (label procI : processorList)
    {
        if (procI == Pstream::myProcNo())
            continue;
        
        UIPstream recv(procI,pBufs);

        // Create temporary lists
        recv >> allPatchTypes[procI];
        recv >> allPatchSize[procI];
        recv >> allStart[procI];
        recv >> allIndex[procI];
    }
}


// * * * * * * * * * * * Depricated Functions * * * * * * * * * * * * * * * * *

Foam::fileName Foam::reconstructRegionalMesh::localPath
(
    const fvMesh& localMesh,
    const label proci,
    const fileName file
)
{        
    // Create start time value
    const scalar startTimeValue = (localMesh.time().startTime()).value();

    // Fist check if a mesh is present in the time directory
    fileName pathToTimeDir = localMesh.time().path().path()
          / fileName("processor" + name(proci) + "/" + localMesh.time().timeName(startTimeValue) + "/" + polyMesh::meshSubDir)
          / file;

    if (exists(pathToTimeDir))
        return pathToTimeDir;

    return  localMesh.time().path().path()
          / fileName("processor" + name(proci) + "/constant/" + polyMesh::meshSubDir)
          / file;
}


void Foam::reconstructRegionalMesh::readHeader(Istream& is)
{
     // Check Istream not already bad
     if (!is.good())
     {
        FatalIOErrorInFunction(is)
             << " Stream is not good"  
             << exit(FatalIOError);
     }
 
     token firstToken(is);
 
     if
     (
         is.good()
      && firstToken.isWord()
      && firstToken.wordToken() == "FoamFile"
     )
     {
         dictionary headerDict(is);
 
         #if (OF_FORK_VERSION >= 1912 && OF_FORK_VERSION < 2006)
            is.version(headerDict.lookup("version"));
            is.format(headerDict.get<word>("format"));
            const word headerClassName = headerDict.get<word>("class");
            const word headerObject(headerDict.get<word>("object"));
         #elif (OF_FORK_VERSION >= 2006 )
            is.version(headerDict.get<token>("version"));
            is.format(headerDict.get<word>("format"));
            const word headerClassName = headerDict.get<word>("class");
            const word headerObject(headerDict.get<word>("object"));
         #else
            is.version(headerDict.lookup("version"));
            is.format(headerDict.lookup("format"));
            const word headerClassName = word(headerDict.lookup("class"));
            const word headerObject(headerDict.lookup("object"));
         #endif
 
        //if (headerClassName == "faceCompactList")
            //FatalIOErrorInFunction(is)
             //<< " stream failure while "
             //<< " reading file " << is.name() << nl
             //<< "\tFile is in binary format. Mesh has to be saved in ascii format"
             //<< exit(FatalIOError);
 
         
         #ifdef FULLDEBUG
             if (headerObject != headerClassName)
             {
                 IOWarningInFunction(is)
                     << " object renamed from "
                     << headerClassName << " to " << headerObject
                     << " for file " << is.name() << endl;
             }
         #endif
     }
     else
     {
        FatalIOErrorInFunction(is)
             << " stream failure while reading header"
             << " on line " << is.lineNumber()
             << " of file " << is.name()
             << exit(FatalIOError);
     }
 
     // Check stream is still OK
     if (!is.good())
     {
         FatalIOErrorInFunction(is)
             << " stream failure while reading header"
             << " on line " << is.lineNumber()
             << " of file " << is.name()    
             << exit(FatalIOError);
     }
}


Foam::List<Foam::face> Foam::reconstructRegionalMesh::readFaceList
(
    const fileName path
)
{
    // Create an IFStream object
    IFstream is(path);
    readHeader(is);
    
    if (is.format() == IOstream::streamFormat::BINARY)
    {
        List<face> L;
        
        // Read compact
        const labelList start(is);
        const List<label> elems(is);
    
        // Convert
        L.setSize(start.size()-1);
    
        forAll(L, i)
        {
            face& subList = L[i];
    
            label index = start[i];
            subList.setSize(start[i+1] - index);
    
            forAll(subList, j)
            {
                subList[j] = elems[index++];
            }
        }
        
        return L;
    }
    else
        return List<face>(is);
}

void Foam::reconstructRegionalMesh::fileHandlerControl::setUncollated()
{
    // Store old file handler type
    if (isA<fileOperations::masterUncollatedFileOperation>(fileHandler()))
        oldFileHandlerType = "collated";
    else
        oldFileHandlerType = "uncollated";

    // For the generation of the partial mesh the uncollated file opteration 
    // has to be used, as otherwise the MPI gets stuck calling an allGather()
    fileOperation::fileHandlerPtr_ = fileOperation::New("uncollated",false);
    #if (OF_FORK_VERSION >= 2006 )
        fileHandler(std::move(fileOperation::fileHandlerPtr_));
    #else
        fileHandler(fileOperation::fileHandlerPtr_);
    #endif
}


void Foam::reconstructRegionalMesh::fileHandlerControl::reset()
{
    // Reset file handler to default
    fileOperation::fileHandlerPtr_ = fileOperation::New(oldFileHandlerType,false);
    #if (OF_FORK_VERSION >= 2006 )
        fileHandler(std::move(fileOperation::fileHandlerPtr_));
    #else
        fileHandler(fileOperation::fileHandlerPtr_);
    #endif
}

Foam::reconstructRegionalMesh::fileHandlerControl::~fileHandlerControl()
{
    // call reset before destruction
    reset();
}
