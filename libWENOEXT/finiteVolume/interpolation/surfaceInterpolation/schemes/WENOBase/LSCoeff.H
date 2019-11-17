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
    Foam::LSCoeff

Description
    LS base class for runtime operations of LS schemes

SourceFiles
    LSCoeff.C
     
Author
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.    

\*---------------------------------------------------------------------------*/

#ifndef LSCoeff_H
#define LSCoeff_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LSCoeff Declaration                         
\*---------------------------------------------------------------------------*/

template<class Type>
class LSCoeff
{
private:

    // Private Data
           
        //- Typedef for 3D scalar matrix        
		typedef List< List< List<scalar> > > scalarMatrix;		

		//- WENO weighting factors
		scalar p_;
		scalar dm_;

		//- Dimensionality of the geometry
		//  Individual for each stencil
        labelListList* dimList_;		
        
		//- Order of polynomials 
        const scalar polOrder_;								
		
		//- Number of derivatives        
        label nDvt_;
		
		//- Lists of central and sectorial stencil ID's for each cell 
		List<labelListList>* stencilsID_;
		
		//- Information about cell property
		//  Could be extended to boundary conditions 
		//  - -1 : local cell 
		//  - >-1: halo cell
		List<labelListList>* cellToPatchMap_;
		
		//- Information about processor neighbours of patches 			
		labelList* patchToProcMap_;	
		
		//- Lists of surface integrals of basis functions 
		//  Calculated in the reference space		
		List<List<scalarMatrix> >* intBasTrans_;			

		//- List of face areas in the reference space
		List<scalarList>* refFacAr_;	
		
		//- Lists of field values of halo cells						
		List<List<Type> > haloData_;	
		
		//- Lists of center coordinates of halo cells		
		List<List<point> >* haloCenters_;
		
		//- Lists of ID's of halo cells
		labelListList* ownHalos_;		
		
		//- Lists of inverse Jacobians for each stencil
		List<scalarRectangularMatrix>* JInv_;			
		
		//- Lists of pseudoinverses for each stencil of each cell	
		List<List<scalarRectangularMatrix> >* LSmatrix_;		


    // Private Member Functions

        //- Disallow default bitwise copy construct
        LSCoeff(const LSCoeff&);

        //- Disallow default bitwise assignment
        void operator=(const LSCoeff&);

		//- Distribute data if multiple processors are involved
        void collectData(const fvMesh& mesh);   


public:

    // Constructor

        LSCoeff
        (
			const fvMesh& mesh,
			const label polOrder
		) 
        : 
			polOrder_(polOrder)           
		{											
			// 3D version
			if (mesh.nGeometricD() == 3)
			{
				nDvt_ = 
					(polOrder_ + 1.0)*(polOrder_ + 2.0)*(polOrder_ + 3.0) 
				   /6.0 - 1.0;
				
				Info<< "Reconstruction using WENO"
				    << polOrder_ << " (3D version)" << endl;				
			}
			else // 2D version (always only one cell in z-direction)								
			{	
				nDvt_ = (polOrder_ + 1.0)*(polOrder_ + 2.0)/2.0 - 1.0;
				
				Info<< "Reconstruction using WENO"
				    << polOrder_ << " (2D version)" << endl;			
			}
		}

    // Member Functions	       
  
        //- Calling function from different schemes 
        Field<Field<Type> > getWENOPol
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        )	;  
  
		//- Calculating the coefficients for each stencil of each cell 
		void calcCoeff
		(		
			const label cellI,
			const GeometricField<Type, fvPatchField, volMesh>& dataField,
			Field<Type>& dvtI,
			const label stencilI,
			const label excludeStencils
        )	;
		
		//- Get necessary lists for runtime operations
        inline List<List<scalarMatrix> >** getPointerIntBasTrans()
        {
            return &intBasTrans_;
        }; 		  
        inline List<scalarList>** getPointerRefFacAr()
        {
            return &refFacAr_;
        }; 	
        inline labelListList** getPointerDimList()
        {
            return &dimList_;
        }; 	  
        inline List<scalarRectangularMatrix>** getPointerJInv()
        {
            return &JInv_;
        }; 	                             
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LSCoeff.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
