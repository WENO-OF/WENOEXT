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
    Foam::WENOSensor

Description
    WENOSensor to detect discontinuities
    Idea is based on the work of Zhao [1]. The sensor proposed there does not
    work as we have different definition of the stencil weights. 
    
    The shock sensor is based on the maximum of the smoothness indicator
    function.
    
    References:
    \verbatim
      [1] Guo-Yan Zhao, Ming-Bo Sun, Sergio Pirozzoli,
          "On shock sensors for hybrid compact/WENO schemes",
          Computers & Fluids,Volume 199,2020,104439,ISSN 0045-7930,
          https://doi.org/10.1016/j.compfluid.2020.104439.
     
    \endverbatim

SourceFiles
    WENOSensor.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020

\*---------------------------------------------------------------------------*/

#ifndef WENOSensor_H
#define WENOSensor_H

#include "DynamicField.H"
#include "WENOCoeff.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class WENOSensor Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class WENOSensor
:
    public WENOCoeff<Type>
{
    private:
    // Private Member Functions

        //- Exponent for shock sensor
        scalar theta_;

        using coeffType = blaze::DynamicVector<Type>;

        //- Disallow default bitwise copy construct
        WENOSensor(const WENOSensor&);

        //- Disallow default bitwise assignment
        void operator=(const WENOSensor&);
        
        //- Get weighted combination for any other type
        virtual void calcWeight
        (
            Field<Type>& coeffsWeightedI,
            const label cellI,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const DynamicList<coeffType>& coeffsI
        ) const override;

public:

    // Constructor

        WENOSensor
        (
            const fvMesh& mesh,
            const label polOrder
        );
        
        GeometricField<Type, fvPatchField, volMesh>& getShockSensor() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WENOSensor.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
