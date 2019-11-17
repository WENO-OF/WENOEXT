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
    Foam::fv::LSGrad

Description
    LS gradient scheme returns a volField for cell centers

SourceFiles
    LSGrad.C

Author
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef LSGrad_H
#define LSGrad_H

#include "gradScheme.H"
#include "processorFvPatch.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class LSGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LSGrad
:
    public fv::gradScheme<Type>
{
    // Private Data

        //- Typedef for resulting gradField            
		typedef typename outerProduct<vector, Type>::type GradType;           	         
		
		//- Dimensionality of the geometry
		//  Individual for each stencil
        labelListList** dimList_;			
		
		//- Lists of inverse Jacobians for each stencil
		List<scalarRectangularMatrix>** JInv_;			
		
        //- polynomial order
        const scalar polOrder_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        LSGrad(const LSGrad&);

        //- Disallow default bitwise assignment
        void operator=(const LSGrad&);


public:

    //- Runtime type information
    TypeName("LSGrad");


    // Constructors

        //- Construct from Istream
        LSGrad(const fvMesh& mesh, Istream& is)
        :
            gradScheme<Type>(mesh),
            polOrder_(readScalar(is))
        {}


    // Member Functions

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;     


        vector matTVecProd					
        (
			const scalarRectangularMatrix& JInvI,
			const Field<scalar>& gradXi
        )	const
        {			
			vector result = pTraits<vector>::zero;
		
			for (label i=0; i < JInvI.n(); i++)
			{			
				for (label j=0; j < JInvI.m(); j++)
				{		
					result[i] += JInvI[j][i] * gradXi[j];				
				}				
			}																		
			return result;
		}		 
		
        tensor matTVecProd					
        (
			const scalarRectangularMatrix& JInvI,
			const Field<vector>& gradXi
        )	const
        {						
			tensor result = pTraits<tensor>::zero;
		
			result.xx() = gradXi[0][0] * JInvI[0][0] 
				+ gradXi[1][0] * JInvI[1][0] 
				+ gradXi[2][0] * JInvI[2][0]; 
			result.xy() = gradXi[0][1] * JInvI[0][0] 
				+ gradXi[1][1] * JInvI[1][0] 
				+ gradXi[2][1] * JInvI[2][0]; 
			result.xz() = gradXi[0][2] * JInvI[0][0] 
				+ gradXi[1][2] * JInvI[1][0] 
				+ gradXi[2][2] * JInvI[2][0]; 
			
			result.yx() = gradXi[0][0] * JInvI[0][1] 
				+ gradXi[1][0] * JInvI[1][1] 
				+ gradXi[2][0] * JInvI[2][1]; 
			result.yy() = gradXi[0][1] * JInvI[0][1] 
				+ gradXi[1][1] * JInvI[1][1] 
				+ gradXi[2][1] * JInvI[2][1]; 
			result.yz() = gradXi[0][2] * JInvI[0][1] 
				+ gradXi[1][2] * JInvI[1][1] 
				+ gradXi[2][2] * JInvI[2][1]; 
			
			result.zx() = gradXi[0][0] * JInvI[0][2] 
				+ gradXi[1][0] * JInvI[1][2] 
				+ gradXi[2][0] * JInvI[2][2]; 
			result.zy() = gradXi[0][1] * JInvI[0][2] 
				+ gradXi[1][1] * JInvI[1][2] 
				+ gradXi[2][1] * JInvI[2][2];
			result.zz() = gradXi[0][2] * JInvI[0][2] 
				+ gradXi[1][2] * JInvI[1][2] 
				+ gradXi[2][2] * JInvI[2][2];			
											
			return result;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LSGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
