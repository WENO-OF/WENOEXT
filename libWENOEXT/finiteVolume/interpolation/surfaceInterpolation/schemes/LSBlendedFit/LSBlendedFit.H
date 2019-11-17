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
    Foam::LSBlendedFit

Description
 	LS interpolation scheme class using a blended combination of linear and
    upwind. Suitable for linearised convection terms.
 
SourceFiles
    LSBlendedFit.C

Author
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef LSBlendedFit_H
#define LSBlendedFit_H

#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class LSBlendedFit Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LSBlendedFit
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data
           
        //- Typedef for 3D scalar matrix    
		typedef List< List< List<scalar> > > scalarMatrix;
				
		//- Lists of surface integrals of basis functions 
		//  Calculated in the reference space		
		List<List<scalarMatrix> >** intBasTrans_;	
	
		//- List of face areas in the reference space
		List<scalarList>** refFacAr_;		
	
		//- Dimensionality of the geometry
		//  Individual for each stencil
        labelListList** dimList_;	
        
        //- SurfaceScalarField of U() & Sf()	
		const surfaceScalarField& faceFlux_;			
		
        //- Polynomial order
        //  User defined parameter        
        const scalar polOrder_;    
        
        //- Name of blending field
        //  Field should be available at runtime
        word blendingFactorName_;         
        	
    // Private Member Functions

        //- Disallow default bitwise copy construct
        LSBlendedFit(const LSBlendedFit&);

        //- Disallow default bitwise assignment
        void operator=(const LSBlendedFit&);

		//- Distribute polynomials through coupled patches
        void swapData
        (
			const fvMesh& mesh,
		    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
				GeometricBoundaryField& btsfCd
        )   const;   
          
        void coupledRiemannSolver
        (
			const fvMesh& mesh,
		    GeometricField<Type, fvsPatchField, surfaceMesh>& tsf,
		    const GeometricField<Type, fvPatchField, volMesh>& vf,
		    const surfaceScalarField& weights,
		    const surfaceScalarField& blendingFactor,
			Field<Field<Type> > coeffsWeighted
        )   const; 
        

public:

    //- Runtime type information
    TypeName("LSBlendedFit");


    // Constructors

        //- Construct from mesh
        LSBlendedFit
        (
			const fvMesh& mesh, 
			const label polOrder
		)					
        :															
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(zeroFlux()),
            polOrder_(polOrder),
            blendingFactorName_("StreletsBlendingFactor") // Dummy name
		{}


        //- Construct from mesh and Istream
        LSBlendedFit
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
            ),
            polOrder_(readScalar(is)),
            blendingFactorName_(is)
		{}


        //- Construct from mesh, faceFlux and Istream
        LSBlendedFit
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux),
            polOrder_(readScalar(is)),
            blendingFactorName_(is)
        {}


    // Member Functions	       
    

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }
        
        tmp<surfaceScalarField> weights
        (
        const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const ;
        
        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        )	const ;  
		
		//- Calculating the face flux values
		Type sumFlux	
		(        	
			const labelList& dim,
			const Field<Type>& coeffcI,
			const scalarMatrix intBasiscIfI
		) 	const;   	
         
        surfaceScalarField zeroFlux() const 
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
		} ; 
        
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#ifdef NoRepository
#   include "LSBlendedFit.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
