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

#include "LSCoeff.H"
#include "LSBase.H"

#include "processorFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void Foam::LSCoeff<Type>::calcCoeff			
(
	const label cellI,
	const GeometricField<Type, fvPatchField, volMesh>& vf,
	Field<Type>& coeffI,
	const label stencilI,
	const label excludeStencils	
)															 
{ 	
	const List<label>& stencilsIDI = 
	    (*stencilsID_)[cellI][stencilI + excludeStencils];
	const scalarRectangularMatrix& A = 
		(*LSmatrix_)[cellI][stencilI];
	const List<label>& cellToPatchMapI = 
		(*cellToPatchMap_)[cellI][stencilI + excludeStencils];
		
    // Calculate degrees of freedom of stencil as a matrix vector product
    // First line is always constraint line 
	
	Type bJ = pTraits<Type>::zero;
	
	for (label j = 1; j < stencilsIDI.size(); j++)									
    {		
		// Distinguish between local and halo cells
		if (cellToPatchMapI[j] == -1)
		{	
			bJ = vf[stencilsIDI[j]] - vf[cellI];	
			
			for (label i = 0; i < A.n(); i++)
			{	
				coeffI[i] += A[i][j-1]*bJ;
			}	
		}
		else if(cellToPatchMapI[j] > -4)
		{
			bJ = 
			    haloData_[cellToPatchMapI[j]][stencilsIDI[j]] 
			  - vf[cellI];
			
			for (label i = 0; i < A.n(); i++)
			{	
				coeffI[i] += A[i][j-1]*bJ;
			}			
		}		
	}	
}


template<class Type>
void Foam::LSCoeff<Type>::collectData(const fvMesh& mesh)	
{
	PstreamBuffers pBufs(Pstream::nonBlocking);							
			
	// Distribute data
	forAll((*patchToProcMap_), patchI)
	{
		if ((*patchToProcMap_)[patchI] > -1)								
		{
			UOPstream toBuffer((*patchToProcMap_)[patchI], pBufs);			
			toBuffer << haloData_[patchI];															
		}
	}
        
	pBufs.finishedSends();										
			
	// Collect data				
	forAll((*patchToProcMap_), patchI)
	{		
		if ((*patchToProcMap_)[patchI] > -1)								
		{
			haloData_[patchI].clear();
			
			UIPstream fromBuffer((*patchToProcMap_)[patchI], pBufs);		
			fromBuffer >> haloData_[patchI];							   									
		}
		else if ((*patchToProcMap_)[patchI] == -2)	
		{
			const cyclicPolyPatch& cycPatch = 
				refCast<const cyclicPolyPatch>(mesh.boundaryMesh()[patchI]);
		
			if (cycPatch.owner() == true)
			{			
				const label neighbPatchID = cycPatch.neighbPatchID();

				const List<Type> ownHaloData = haloData_[patchI];
				
				haloData_[patchI] = haloData_[neighbPatchID];
				haloData_[neighbPatchID] = ownHaloData;			
			}
		}	
	}
}


template<class Type> 
Foam::Field<Foam::Field<Type> > 
Foam::LSCoeff<Type>::getWENOPol
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) 
{
    const fvMesh& mesh = vf.mesh();	
	const fvPatchList& patches = mesh.boundary();
	
	// Get preprocessing lists from LSBase class
	
	Foam::LSBase& init = 
	    LSBase::instance
		(
			mesh,
			polOrder_
		);	

    stencilsID_ = init.getPointerStencilID();
    cellToPatchMap_ = init.getPointerCellToPatchMap();
    patchToProcMap_ = init.getPointerPatchToProcMap();    
    haloCenters_ = init.getPointerHaloCenters();
    ownHalos_ = init.getPointerOwnHalos();    
    LSmatrix_ = init.getPointerLSmatrix();
	intBasTrans_ = init.getPointerIntBasTrans();	
	refFacAr_ = init.getPointerRefFacAr();	
	dimList_ = init.getPointerDimList();
	JInv_ = init.getPointerJInv();		
	
	// Distribute data to neighbour processors
	
	haloData_.setSize((*ownHalos_).size());
			
	forAll(patches, patchI)
	{
		if (isA<processorFvPatch>(patches[patchI]) ||
			isA<cyclicFvPatch>(patches[patchI])
		   )			
		{
			forAll(haloData_, patchI)
			{
				haloData_[patchI].setSize((*ownHalos_)[patchI].size());
			
				forAll(haloData_[patchI], cellI)
				{
					haloData_[patchI][cellI] = 
						vf.internalField()[(*ownHalos_)[patchI][cellI]];
				}
			}
		}
	}	

	collectData(mesh);

	// Runtime operations
		
	Field<Field<Type> > coeffsLS(mesh.nCells());					
	
	for (label cellI = 0; cellI < mesh.nCells(); cellI++)
	{
		// Calculate degrees of freedom for central stencil of the cell						
		coeffsLS[cellI].setSize(nDvt_,pTraits<Type>::zero);	
		
		calcCoeff												
		( 
			cellI,
			vf,
			coeffsLS[cellI],
			0,
			0
		);	
	}
	
	return coeffsLS;		
} 

// ************************************************************************* //
