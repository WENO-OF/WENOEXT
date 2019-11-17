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

\*---------------------------------------------------------------------------*/

#include "PeBlendedField.H"
//#include "volFields.H"
#include "dictionary.H"
#include "fvCFD.H"
//#include  "turbulenceModel.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PeBlendedField, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

/*bool Foam::PeBlendedField::isKinematicPressure()
{
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

    return p.dimensions() == sqr(dimLength)/sqr(dimTime);
}*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PeBlendedField::PeBlendedField
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    writeBlendingField_(dict.lookupOrDefault<bool>("write",false)),
    // Max 100%
    maxFirstSchemePercent_(dict.lookupOrDefault<scalar>("maxValue",1)),
    // Min 0%
    minFirstSchemePercent_(dict.lookupOrDefault<scalar>("minValue",0)),
    mesh_(refCast<const fvMesh>(obr_)),    
    StreletsBlendingFactor_
    (
        IOobject
        (
            "PeBlendingFactor",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("null",dimless,0.)
    )
{

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PeBlendedField::~PeBlendedField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PeBlendedField::read(const dictionary& dict)
{

}


void Foam::PeBlendedField::execute()
{

    Info << "Updating PeBlendedField" << endl;
    //~ const volScalarField& rans = obr_.lookupObject<volScalarField>("rans");
    //~ const surfaceScalarField 
    
    const Foam::incompressible::turbulenceModel& turbMod =
         obr_.lookupObject<Foam::incompressible::turbulenceModel>("turbulenceModel");
    const surfaceScalarField& phi = turbMod.phi();


    surfaceScalarField Pe
    (
        IOobject
        (
            "Pe",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ
        ),
        mag(phi)
       /(
            mesh_.magSf()
          * mesh_.surfaceInterpolation::deltaCoeffs()
          * fvc::interpolate(turbMod.nuEff()))
    );

    // localBlendedBy:  bf * scheme1 + (1-bf)*scheme2
    // let say, you write linearUpwindV grad(U) linear
    // using rans a blending factor you get linearUpwind in RANS region
    // and linear in LES region
    
    StreletsBlendingFactor_ = neg(scalar(2) - Pe)*maxFirstSchemePercent_ +
    (1 - neg(scalar(2) - Pe))*minFirstSchemePercent_;
    
    Info << "Average Pe Blending Factor" << gAverage(StreletsBlendingFactor_) << endl;
    //Info << min(StreletsBlendingFactor_) << endl;
    
    //~ if(obr_.foundObject<surfaceScalarField>(qualityBlendingFieldName_))
    //~ {
        //~ StreletsBlendingFactor_ = max
        //~ (
            //~ StreletsBlendingFactor_,
            //~ obr_.lookupObject<surfaceScalarField>(qualityBlendingFieldName_)
        //~ );
    //~ }

    //~ StreletsBlendingFactor_.min(maxFirstSchemePercent_);
    //~ StreletsBlendingFactor_.max(minFirstSchemePercent_);
}


void Foam::PeBlendedField::end()
{
    // Do nothing - only valid on write

}

void Foam::PeBlendedField::timeSet()
{
    // Do nothing
}



void Foam::PeBlendedField::write()
{
    // write PeBlendedField
    if (active_ && writeBlendingField_)
    {
        volScalarField StreletsBlendingFactorVol("SBF",fvc::average(StreletsBlendingFactor_));
        StreletsBlendingFactorVol.write();
        
    }
}


// ************************************************************************* //
