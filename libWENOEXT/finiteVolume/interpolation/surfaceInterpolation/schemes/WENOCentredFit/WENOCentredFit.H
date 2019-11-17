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
 	WENO interpolation scheme class using a linear combination. Suitable
 	for linearised convection terms.
 
SourceFiles
    WENOCentredFit.C

Author
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef WENOCentredFit_H
#define WENOCentredFit_H

#include "surfaceInterpolationScheme.H"

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
		typedef List< List< List<scalar> > > scalarMatrix;
				
		//- Lists of surface integrals of basis functions 
		//  Calculated in the reference space		
		List<List<scalarMatrix> >** intBasTrans_;	
	
		//- List of face areas in the reference space
		List<scalarList>** refFacAr_;		
	
		//- Dimensionality of the geometry
		//  Individual for each stencil
        labelListList** dimList_;				
		
        //- Polynomial order
        //  User defined parameter        
        const scalar polOrder_;      
        

    // Private Member Functions

        //- Disallow default bitwise copy construct
        WENOCentredFit(const WENOCentredFit&);

        //- Disallow default bitwise assignment
        void operator=(const WENOCentredFit&);

		//- Distribute polynomials through coupled patches
        void swapData
        (
			const fvMesh& mesh,
		    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
				GeometricBoundaryField& btsf
        )   const;   
          
        void coupledRiemannSolver
        (
			const fvMesh& mesh,
		    GeometricField<Type, fvsPatchField, surfaceMesh>& tsf,
		    const GeometricField<Type, fvPatchField, volMesh>& vf,
		    const surfaceScalarField& weights,
			Field<Field<Type> > coeffsWeighted
        )   const; 
        
 
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
            polOrder_(polOrder)        
		{}

        //- Construct from mesh and Istream
        WENOCentredFit
        (
			const fvMesh& mesh, 
			Istream& is			
		)					
        :															
            surfaceInterpolationScheme<Type>(mesh),
            polOrder_(readScalar(is))
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
            polOrder_(readScalar(is))
        {}


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
        )	const ;  

		//- Calculating the face flux values
		Type sumFlux	
		(        	
			const labelList& dim,
			const Field<Type>& coeffcI,
			const scalarMatrix intBasiscIfI
		) 	const;          
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#ifdef NoRepository
#   include "WENOCentredFit.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
