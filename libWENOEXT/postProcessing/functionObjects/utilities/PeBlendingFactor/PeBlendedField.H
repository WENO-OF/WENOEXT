/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::PeBlendedField

SourceFiles
    PeBlendedField.C
    IOPeBlendedField.H

Description
    Blending factor for hybrid RANS-LES computations on low-quality meshes
    Takes into account not only the RANS-LES subdivision in the computational
    domain, but also the mesh quality (optional) name for mesh quality field
    can be specified

\*---------------------------------------------------------------------------*/

#ifndef PeBlendedField_H
#define PeBlendedField_H

#include "functionObjectFile.H"
#include "pointFieldFwd.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;
class polyMesh;
class mapPolyMesh;
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class PeBlendedField Declaration
\*---------------------------------------------------------------------------*/

class PeBlendedField
{
    // Private data

        //- Name of this set of PeBlendedField objects
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;
        bool writeBlendingField_;

        scalar maxFirstSchemePercent_  ;
        scalar minFirstSchemePercent_  ;
                
        // get access to the mesh
        const fvMesh& mesh_;
        //- 
        
        surfaceScalarField StreletsBlendingFactor_;
        
    // Private Member Functions

        //- Return true if the pressure field corresponds to kinematic pressure
//        bool isKinematicPressure();

        //- Disallow default bitwise copy construct
        PeBlendedField(const PeBlendedField&);

        //- Disallow default bitwise assignment
        void operator=(const PeBlendedField&);


public:

    //- Runtime type information
    TypeName("PeBlendingFactor");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        PeBlendedField
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    // Destructor

        virtual ~PeBlendedField();


    // Member Functions

        //- Return name of the set of PeBlendedField
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the PeBlendedField data
        virtual void read(const dictionary&);

        //- Execute
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Write the PeBlendedField at specified write time
        virtual void write();
        
        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Update for changes of mesh, not needed here
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh, not needed here
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
