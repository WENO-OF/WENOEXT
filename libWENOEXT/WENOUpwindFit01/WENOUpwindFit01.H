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

Class
    Foam::WENOUpwindFit01

Description
    Specialized WENO interpolation scheme for bounded transport between 0 
    and 1.
    
SourceFiles
    WENOUpwindFit01.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020

\*---------------------------------------------------------------------------*/

#ifndef WENOUpwindFit01_H
#define WENOUpwindFit01_H

#include "codeRules.H"
#include "surfaceInterpolationScheme.H"
#include "WENOCoeff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class WENOUpwindFit01 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class WENOUpwindFit01
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data

        //- Typedef for 3D scalar matrix
        using volIntegralType = List3D<scalar>;
        
        //- C++11 typedef for squareMatrix
        //  This is used for Jacobian matrix
        using scalarSquareMatrix = SquareMatrix<scalar>;

        //- SurfaceScalarField of U() & Sf()
        const surfaceScalarField& faceFlux_;

        //- Polynomial order
        //  User defined parameter
        const scalar polOrder_;

        //- Reference to WENOCoeff class
        WENOCoeff<Type> WENOCoeff_;

        //- Reference to WENOBase
        const WENOBase& WENOBase_;


    // Private Member Functions

        //- Distribute data through coupled patches and solve Riemann problem
        void swapData
        (
            const fvMesh& mesh,
            typename GeometricField<Type, fvsPatchField, surfaceMesh>::
            #ifdef FOAM_NEW_GEOMFIELD_RULES
                Boundary& btsf
            #else 
                GeometricBoundaryField& btsf
            #endif
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
            const Field<Field<Type> >& coeffsWeighted,
            GeometricField<Type, fvsPatchField, surfaceMesh>& tsfP
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


public:

    //- Runtime type information
    TypeName("WENOUpwindFit01");


    // Constructors

        //- Construct from mesh (not used)
        WENOUpwindFit01
        (
            const fvMesh& mesh,
            const label polOrder
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(zeroFlux()),
            polOrder_(polOrder),
            WENOCoeff_(mesh,polOrder_),
            WENOBase_(WENOCoeff_.WENOBaseRef())
        {}

        //- Construct from mesh and Istream
        WENOUpwindFit01
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
            WENOCoeff_(mesh,polOrder_),
            WENOBase_(WENOCoeff_.WENOBaseRef())
        {}

        //- Construct from mesh, faceFlux and Istream
        WENOUpwindFit01
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux),
            polOrder_(readScalar(is)),
            WENOCoeff_(mesh,polOrder_),
            WENOBase_(WENOCoeff_.WENOBaseRef())
        {}


        //- Disallow default bitwise copy construct
        WENOUpwindFit01(const WENOUpwindFit01&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const WENOUpwindFit01&) = delete;

    // Member Functions

        //- Return the interpolation weighting factors for implicit part
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            return pos(faceFlux_);
        }

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
#   include "WENOUpwindFit01.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
