/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Author
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "WENOGrad.H"
#include "gaussGrad.H"

#include "WENOCoeff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::WENOGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vtf,
    const word& name
) const
{
	const fvMesh& mesh = this->mesh();	

	// Get degrees of freedom from WENOCoeff class

	Foam::WENOCoeff<Type> getWeights(mesh, polOrder_);
	Field<Field<Type> > coeffsWeighted = getWeights.getWENOPol(vtf);

	WENOGrad *ptr = const_cast<WENOGrad*>(this);
	
	ptr->dimList_ = getWeights.getPointerDimList();	
	ptr->JInv_	  = getWeights.getPointerJInv();		
	
    tmp<GeometricField<GradType, fvPatchField, volMesh> > twenoGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                name,
                vtf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "zero",
                vtf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& wenoGrad = twenoGrad();    


	//- calculate gradients in cell centers in physical coordinates	
	//- depends on dimension of domain; gradXi = gradXi, gradEta, gradZeta
	
	Field<Type> gradXi(3,pTraits<Type>::zero);	
	
	if (   (**dimList_)[0][0] == (**dimList_)[0][1] 
		&& (**dimList_)[0][0] == (**dimList_)[0][2])	// 3D
	{
		for (label cellI = 0; cellI<mesh.nCells(); cellI++)
		{
			gradXi[2] = coeffsWeighted[cellI][0];
			gradXi[1] = coeffsWeighted[cellI][polOrder_];

			if (polOrder_ == 1) gradXi[0] = coeffsWeighted[cellI][2];			
			else if (polOrder_ == 2) gradXi[0] = coeffsWeighted[cellI][5];		
			else if (polOrder_ == 3) gradXi[0] = coeffsWeighted[cellI][9];		
			else if (polOrder_ == 4) gradXi[0] = coeffsWeighted[cellI][14];		
			else
			{
				Info<<"not implemented yet"<<endl;
			}
		
			wenoGrad[cellI] = matTVecProd((**JInv_)[cellI], gradXi);
		}
	}
	else // 2D
	{
		for (label cellI = 0; cellI<mesh.nCells(); cellI++)
		{
			if ( (**dimList_)[cellI][0] == 0)	// xi = 0
			{
				gradXi[2] = coeffsWeighted[cellI][0]; 		
				gradXi[1] = coeffsWeighted[cellI][polOrder_];	
				gradXi[0] = pTraits<Type>::zero;	 	
					
				wenoGrad[cellI] = matTVecProd((**JInv_)[cellI], gradXi);
			}
			else if ( (**dimList_)[cellI][1] == 0)	// eta = 0
			{
				gradXi[2] = coeffsWeighted[cellI][0];			
				gradXi[1] = pTraits<Type>::zero; 	
				gradXi[0] = coeffsWeighted[cellI][polOrder_];	
					
				wenoGrad[cellI] = matTVecProd((**JInv_)[cellI], gradXi);
			}
			else if ( (**dimList_)[cellI][2] == 0)	// zeta = 0
			{
				gradXi[2] = pTraits<Type>::zero;			
				gradXi[1] = coeffsWeighted[cellI][0];	
				gradXi[0] = coeffsWeighted[cellI][polOrder_];

				wenoGrad[cellI] = matTVecProd((**JInv_)[cellI], gradXi);
			}						
		}
	}

	// Coupled boundary faces

    // Correct the boundary conditions
    // Currently linear reconstruction
    wenoGrad.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vtf, wenoGrad);			

    return twenoGrad;
}


// ************************************************************************* //
