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
    Foam::WENOLinearFit

Description
    WENO interpolation scheme class using an exact Riemann solver. The 
    explicit correction vector is used to determine the implicit neighbor
    weights.
    Essentially this scheme is similar to limitedLinear, utilizing the 
    WENO polynoms as limiter function. 

    Usage:

    div(phi,U)  Gauss WENOLinear Fit 3 0.8 1;
    

    First index -- here 3:
        Sets the polynomial order
    Second index -- here 0.8:
        If it is set to 0, it is essentially a linear scheme, if it is set
        to 1, it uses the WENO weights. E.g., if it is set to 0.9 it means 
        that the WENO weight must differ by more than 10% to the linear 
        weight to be used. Otherwise the central differencing 
        (linear scheme) weights are used.
    Third index -- here 1:
        Boolean value to select if the explicit correction is used
        If it is used, the explicit correcion uses a cell limited 
        approach

SourceFiles
    WENOLinearFit.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2026

\*---------------------------------------------------------------------------*/

#ifndef WENOLinearFit_H
#define WENOLinearFit_H

#include "codeRules.H"
#include "surfaceInterpolationScheme.H"
#include "WENOCoeff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class WENOLinearFit Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class WENOLinearFit
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

        //- Store the correction value to avoid recomputation
        mutable autoPtr<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfCorrP_;

        //- Polynomial order
        //  User defined parameter
        const scalar polOrder_;

        //- Limiting factor
        //  User defined:
        //  -  0: not limited -- fully linear scheme
        //  - +1: limited -- use only WENO weights
        const scalar limFac_;

        //- Use explicit correction for upwind terms
        const bool explicitCorrection_;

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

        //- Calculate implicit weight based on explicit correction
        scalar calcWeight
        (
            const Type& corr,
            const Type& psi,
            const Type& psiN
        ) const;


        //- Limit the explicit correction to cell values
        void cellLimitedCorrection
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
    TypeName("WENOLinearFit");


    // Constructors

        //- Construct from mesh (not used)
        WENOLinearFit
        (
            const fvMesh& mesh,
            const label polOrder
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(zeroFlux()),
            polOrder_(polOrder),
            limFac_(1.0),
            explicitCorrection_(false),
            WENOCoeff_(mesh,polOrder_),
            WENOBase_(WENOCoeff_.WENOBaseRef())
        {}

        //- Construct from mesh and Istream
        WENOLinearFit
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
            explicitCorrection_(readBool(is)),
            WENOCoeff_(mesh,polOrder_),
            WENOBase_(WENOCoeff_.WENOBaseRef())
        {}

        //- Construct from mesh, faceFlux and Istream
        WENOLinearFit
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
            explicitCorrection_(readBool(is)),
            WENOCoeff_(mesh,polOrder_),
            WENOBase_(WENOCoeff_.WENOBaseRef())
        {}


        //- Disallow default bitwise copy construct
        WENOLinearFit(const WENOLinearFit&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const WENOLinearFit&) = delete;

    // Member Functions

        //- Return the interpolation weighting factors for implicit part
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const override
        {
            return explicitCorrection_;
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
#   include "WENOLinearFit.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
