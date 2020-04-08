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
    Foam::WENOCentredFit

Description
     WENO interpolation scheme class using a centered scheme for the calculation
     of the face flux
 
SourceFiles
    WENOCentredFit.C

Author
    Jan Wilhelm Gärtner, <jan.gaertner@outlook.de> 

\*---------------------------------------------------------------------------*/

#ifndef WENOCentredFit_H
#define WENOCentredFit_H

#include "codeRules.H"
#include "surfaceInterpolationScheme.H"
#include "WENOCoeff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class WENOCentredFit Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class WENOCentredFit
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data
    
        //- Typedef for 3D scalar matrix
        using volIntegralType = List< List< List<scalar> > > ;
        
        //- C++11 typedef for squareMatrix
        //  This is used for Jacobian matrix
        using scalarSquareMatrix = SquareMatrix<scalar>;

        //- Polynomial order
        //  User defined parameter
        const scalar polOrder_;

        //- Reference to WENOCoeff class
        WENOCoeff<Type> WENOCoeff_;

        //- Reference to WENOBase
        const WENOBase& WENOBase_;


    // Private Member Functions

        void coupledRiemannSolver
        (
            const fvMesh& mesh,
            GeometricField<Type, fvsPatchField, surfaceMesh>& tsf,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const surfaceScalarField& weights,
            const Field<Field<Type> >& coeffsWeighted
        ) const; 
        
        
        //- Calculating the face flux values
        Type sumFlux    
        (            
            const labelList& dim,
            const Field<Type>& coeffcI,
            const volIntegralType& intBasiscIfI
        ) const;     
        
        
        void swapData
        (
            const fvMesh& mesh,
            typename GeometricField<Type, fvsPatchField, surfaceMesh>::
            #ifdef FOAM_NEW_GEOMFIELD_RULES
                Boundary& btsf
            #else 
                GeometricBoundaryField& btsf
            #endif
        ) const;

public:

    //- Runtime type information
    TypeName("WENOCentredFit");


    // Constructors

        //- Construct from mesh
        WENOCentredFit
        (
            const fvMesh& mesh, 
            const label polOrder
        )                    
        :                                                            
            surfaceInterpolationScheme<Type>(mesh),
            polOrder_(polOrder),
            WENOCoeff_(mesh,polOrder_),
            WENOBase_(WENOCoeff_.WENOBaseRef())        
        {}

        //- Construct from mesh and Istream
        WENOCentredFit
        (
            const fvMesh& mesh, 
            Istream& is            
        )                    
        :                                                            
            surfaceInterpolationScheme<Type>(mesh),
            polOrder_(readScalar(is)),
            WENOCoeff_(mesh,polOrder_),
            WENOBase_(WENOCoeff_.WENOBaseRef())
        {}

        //- Construct from mesh, faceFlux and Istream
        WENOCentredFit
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            polOrder_(readScalar(is)),
            WENOCoeff_(mesh,polOrder_),
            WENOBase_(WENOCoeff_.WENOBaseRef())
        {}


        //- Disallow default bitwise copy construct
        WENOCentredFit(const WENOCentredFit&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const WENOCentredFit&) = delete;

    // Member Functions           
    
        //- Return the interpolation weighting factors for implicit part
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const 
        {
            return this->mesh().surfaceInterpolation::weights();
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
        ) const ;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#ifdef NoRepository
#   include "WENOCentredFit.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
