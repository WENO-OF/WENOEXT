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
    Foam::WENOHybrid

Description
    Hybrid WENO scheme using linear as default and using WENO at 
    discontinuities. Detection of the discontinuity is based on WENO weights.
    
    The work is based on Guo-Yan 2020 [1]
    
    References:
    \verbatim
       [1] Guo-Yan Zhao, Ming-Bo Sun, Sergio Pirozzoli,
           "On shock sensors for hybrid compact/WENO schemes",
           Computers & Fluids,Volume 199,2020,104439,ISSN 0045-7930,
           https://doi.org/10.1016/j.compfluid.2020.104439.
    
    \endverbatim
     

SourceFiles
    WENOHybrid.C
    
See also
    WENOSensor.H

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020

\*---------------------------------------------------------------------------*/

#ifndef WENOHybrid_H
#define WENOHybrid_H

#include "surfaceInterpolationScheme.H"
#include "WENOSensor.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class WENOHybrid Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class WENOHybrid
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data

        //- Typedef for 3D scalar matrix
        using volIntegralType = List3D<scalar> ;
        
        //- C++11 typedef for squareMatrix
        //  This is used for Jacobian matrix
        using scalarSquareMatrix = SquareMatrix<scalar>;

        //- SurfaceScalarField of U() & Sf()
        const surfaceScalarField& faceFlux_;

        //- Polynomial order
        //  User defined parameter
        const scalar polOrder_;

        //- Limiting factor
        //  User defined:
        //  -  0: not limited
        //  - +1: limited
        const scalar limFac_;

        //- Limiting factor default is one 
        const scalar limiter_;

        //- Reference to WENOSensor class
        WENOSensor<Type> WENOSensor_;

        //- Reference to WENOBase
        const WENOBase& WENOBase_;

        // Store the correction field 
        mutable autoPtr<GeometricField<Type, fvsPatchField, surfaceMesh> > correctionPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        WENOHybrid(const WENOHybrid&);

        //- Disallow default bitwise assignment
        void operator=(const WENOHybrid&);

        //- Distribute data through coupled patches and solve Riemann problem
        void swapData
        (
            const fvMesh& mesh,
            typename GeometricField<Type, fvsPatchField, surfaceMesh>::
                Boundary& btsf
        )   const;

        void coupledRiemannSolver
        (
            const fvMesh& mesh,
            GeometricField<Type, fvsPatchField, surfaceMesh>& tsfP,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const Field<Field<Type> >& coeffsWeighted
        )   const;

        //- Calculating the face flux values
        Type sumFlux
        (
            const labelList& dim,
            const Field<Type>& coeffcI,
            const volIntegralType& intBasiscIfI
        ) const;

        //- Calculating the polynomial limiters
        void calcLimiter
        (
            const fvMesh& mesh,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            GeometricField<Type, fvsPatchField, surfaceMesh>& tsfP,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& tsfN
        ) const;


        //- Return a zero surfaceScalarField
        //  Required for construct from mesh constructor
        tmp<surfaceScalarField> zeroFlux() const
        {
            tmp<surfaceScalarField> zeroFlux
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "zeroFlux",
                        this->mesh().time().timeName(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar("zeroFlux",dimless,0.0)
                )
            );
            return zeroFlux;
        }
        
        // calculate correction flux
        void calcCorrectionFlux
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
        
        
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > calcShockSensor() const;

public:

    //- Runtime type information
    TypeName("WENOHybrid");


    // Constructors

        //- Construct from mesh (not used)
        WENOHybrid
        (
            const fvMesh& mesh,
            const label polOrder
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(zeroFlux()),
            polOrder_(polOrder),
            limFac_(0),
            limiter_(1),
            WENOSensor_(mesh,polOrder_),
            WENOBase_(WENOSensor_.WENOBaseRef())
        {}

        //- Construct from mesh and Istream
        WENOHybrid
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            ) ,
            polOrder_(readScalar(is)),
            limFac_(readScalar(is)),
            limiter_(readScalar(is)),
            WENOSensor_(mesh,polOrder_),
            WENOBase_(WENOSensor_.WENOBaseRef())
        {}

        //- Construct from mesh, faceFlux and Istream
        WENOHybrid
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux),
            polOrder_(readScalar(is)),
            limFac_(readScalar(is)),
            limiter_(readScalar(is)),
            WENOSensor_(mesh,polOrder_),
            WENOBase_(WENOSensor_.WENOBaseRef())
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const override;
        
        
        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WENOHybrid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
