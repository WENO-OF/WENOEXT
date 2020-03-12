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

Description
    Reconstruct regional mesh.
    Bases on reconstructParMesh.C of OpenFOAM 

SourceFiles
    reconstructRegionalMesh.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020

\*---------------------------------------------------------------------------*/

#include "IOobjectList.H"
#include "labelIOList.H"
#include "processorPolyPatch.H"
#include "mapAddedPolyMesh.H"
#include "polyMeshAdder.H"
#include "faceCoupleInfo.H"
#include "fvMeshAdder.H"
#include "polyTopoChange.H"
#include "extrapolatedCalculatedFvPatchFields.H"


namespace Foam
{

namespace reconstructRegionalMesh
{

    //- Reconsturct mesh depending on processor list
    autoPtr<fvMesh> reconstruct
    (
        const labelList processorList,
        const fvMesh& localMesh
    );

    void renumber
    (
        const labelList& map,
        labelList& elems
    );
    
    autoPtr<faceCoupleInfo> determineCoupledFaces
    (
        const bool fullMatch,
        const label masterMeshProcStart,
        const label masterMeshProcEnd,
        const polyMesh& masterMesh,
        const label meshToAddProcStart,
        const label meshToAddProcEnd,
        const polyMesh& meshToAdd,
        const scalar mergeDist
    );
    
    void mergeSharedPoints
    (
        const scalar mergeDist,
        polyMesh& mesh
    );
    
    boundBox procBounds
    (
        const PtrList<Time>& databases,
        const word& regionDir
    );

    bool findProci(const labelList& processorList,label proci);

}

}