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

Class
    Foam::globalMesh

Description
    Return the global mesh of the case by constructing a new mesh field if 
    it is a prallel run 

SourceFiles
    globalMesh.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020

\*---------------------------------------------------------------------------*/

#ifndef globalMesh_H
#define globalMesh_H

#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace WENO
{
/*---------------------------------------------------------------------------*\
                            Class WENOBase Declaration
\*---------------------------------------------------------------------------*/

class globalMesh
{
    protected:
        
        
        
        //- Pointer to the global mesh
        autoPtr<fvMesh> globalMeshPtr_;
        
        const fvMesh& globalMesh_;
        
        const fvMesh& localMesh_;
        
        //- Global cellDist map 
        const scalarField cellDist_;
        
        //- List of cellID corresponding to cellID in global mesh
        labelList cellID_;
        
    public:
    
        globalMesh(const fvMesh& mesh);
        
        //- Return cellID for processor cells in global mesh
        const labelList& cellID();
        
        const fvMesh& operator()() {return globalMesh_;}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace WENO
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
