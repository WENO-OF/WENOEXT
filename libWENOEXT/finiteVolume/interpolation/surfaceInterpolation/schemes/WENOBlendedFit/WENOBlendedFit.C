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

#include "WENOCoeff.H"
#include "WENOBlendedFit.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::WENOBlendedFit<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
	const fvMesh& mesh = this->mesh();	

	// Get degrees of freedom from WENOCoeff class

	Foam::WENOCoeff<Type> getWeights(mesh, polOrder_);
	Field<Field<Type> > coeffsWeighted = getWeights.getWENOPol(vf);	

	WENOBlendedFit *ptr = const_cast<WENOBlendedFit*>(this);
	
	ptr->intBasTrans_ = getWeights.getPointerIntBasTrans();
	ptr->refFacAr_ = getWeights.getPointerRefFacAr();	
	ptr->dimList_ = getWeights.getPointerDimList();	    
    
    
    // Get the TVD limiters
    
	surfaceScalarField blendingFactor
    (
		tLimitedScheme_().limiter(vf)
    );


	// Calculate the interpolated face values
	
	const labelUList& P = mesh.owner();
	const labelUList& N = mesh.neighbour();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "tvf",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), pTraits<Type>::zero)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& tsf = tsfCorr();
	
	const surfaceScalarField& weights = mesh.surfaceInterpolation::weights();		
		
	// Blended combination of correction polynomials
		
	forAll(P, faceI)	
	{			
		Type valueOwn =
			sumFlux
			(
				(**dimList_)[P[faceI]],
				coeffsWeighted[P[faceI]],
				(**intBasTrans_)[faceI][0]
			)  /(**refFacAr_)[faceI][0];			
		
		Type valueNeigh = 
			sumFlux
			(
				(**dimList_)[N[faceI]],
				coeffsWeighted[N[faceI]],
				(**intBasTrans_)[faceI][1]
			)  /(**refFacAr_)[faceI][1];
			
		Type upwindVal = 
			pos(faceFlux_[faceI])*valueOwn 
		  + (1.0 - pos(faceFlux_[faceI]))*valueNeigh;
			
		Type centredVal =  
			weights[faceI]*valueOwn 
		  + (1.0 - weights[faceI])*valueNeigh;	
			
		tsf[faceI] = blendingFactor[faceI]*centredVal								
		  + (1.0 - blendingFactor[faceI])*upwindVal;
	}
																	
	coupledRiemannSolver
	(
		mesh, 
		tsf, 
		vf, 
		weights, 
		blendingFactor, 
		coeffsWeighted
	);	

	return tsfCorr;		
} 


template<class Type> 
Foam::tmp<Foam::surfaceScalarField> Foam::WENOBlendedFit<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const 
{
	surfaceScalarField blendingFactor
    (
		tLimitedScheme_().limiter(vf)
    );
        
    return 
        (blendingFactor*this->mesh().surfaceInterpolation::weights() 
       +(scalar(1.0) - blendingFactor)*pos(faceFlux_));
}


template<class Type>
Type Foam::WENOBlendedFit<Type>::sumFlux
(        
    const labelList& dim,
    const Field<Type>& coeffcI,
    const scalarMatrix intBasiscIfI
)	const
{
	Type flux = pTraits<Type>::zero;

	label nCoeff = 0;

	for (label n = 0; n <= dim[0]; n++)
	{
		for (label m = 0; m <= dim[1]; m++)
		{
			for (label l = 0; l <= dim[2]; l++)
			{
				if ((n+m+l) <= polOrder_ && (n+m+l) > 0)
				{				
					flux += 
						coeffcI[nCoeff]*intBasiscIfI[n][m][l];
						
					nCoeff++;													
				}
			}
		}
	}	
		
	return flux; 
}	


template<class Type>
void Foam::WENOBlendedFit<Type>::swapData
(
	const fvMesh& mesh,
	typename GeometricField<Type, fvsPatchField, surfaceMesh>::
			GeometricBoundaryField& btsf
)   const
{
	const fvPatchList& patches = mesh.boundary();
	
	PstreamBuffers pBufs(Pstream::nonBlocking);							
		
	// Distribute data
	forAll(btsf, patchI)
	{		
		if (isA<processorFvPatch>(patches[patchI]))
		{			
			UOPstream toBuffer
				(
					refCast<const processorFvPatch>
						(patches[patchI]).neighbProcNo(), 
					pBufs
				);
							
			forAll(btsf[patchI],faceI)
			{
				toBuffer << btsf[patchI][faceI];		
			}													
		}
	}
        
	pBufs.finishedSends();										
			
	// Collect data					
	forAll(btsf, patchI)
	{		
		if (isA<processorFvPatch>(patches[patchI]))
		{
			UIPstream fromBuffer
				(
					refCast<const processorFvPatch>
						(patches[patchI]).neighbProcNo(), 
					pBufs
				);
				
			forAll(btsf[patchI],faceI)
			{						
				fromBuffer >> btsf[patchI][faceI];			
			}				   									
		}
	}
}

template<class Type>
void Foam::WENOBlendedFit<Type>::coupledRiemannSolver
(
	const fvMesh& mesh,
	GeometricField<Type, fvsPatchField, surfaceMesh>& tsf,
	const GeometricField<Type, fvPatchField, volMesh>& vf,
	const surfaceScalarField& weights,
	const surfaceScalarField& blendingFactor,
	Field<Field<Type> > coeffsWeighted
)   const
{
	const fvPatchList& patches = mesh.boundary(); 
		
	typename GeometricField<Type, fvsPatchField, surfaceMesh>::
		GeometricBoundaryField& btsf = tsf.boundaryField();				
					
	tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorrBD
	(
		new GeometricField<Type, fvsPatchField, surfaceMesh>
		(
			IOobject
			(
				"tsfCorrBD",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE,
				false
			),
			mesh,
			dimensioned<Type>
				(vf.name(), vf.dimensions(), pTraits<Type>::zero)
		)
	);
	GeometricField<Type, fvsPatchField, surfaceMesh>& tsfBD = 
		tsfCorrBD();				
				
	typename GeometricField<Type, fvsPatchField, surfaceMesh>::
		GeometricBoundaryField& btsfBD = tsfBD.boundaryField();							
		
	forAll(btsf, patchI)
	{
		fvsPatchField<Type>& pSfCorr = btsf[patchI];					
	
		if (isA<processorFvPatch>(patches[patchI]))
		{					
			const labelUList& pOwner = mesh.boundary()[patchI].faceCells();	
				
			label startFace = patches[patchI].start();
		
			forAll(pOwner, faceI)
			{				
				label own = pOwner[faceI];
					
				btsfBD[patchI][faceI] = 
					sumFlux
					(
						(**dimList_)[own],
						coeffsWeighted[own],
						(**intBasTrans_)[faceI + startFace][0]
					)  /(**refFacAr_)[faceI + startFace][0] ;
						
				pSfCorr[faceI] = btsfBD[patchI][faceI];										
			}
		}
	}

	swapData(mesh, btsfBD);
		
	forAll(btsf, patchI)
	{
		fvsPatchField<Type>& pSfCorr = btsf[patchI];					
	
		if (isA<processorFvPatch>(patches[patchI]))
		{					
			const labelUList& pOwner = mesh.boundary()[patchI].faceCells();

			const scalarField& pFaceFlux = 
				faceFlux_.boundaryField()[patchI];

			forAll(pOwner, faceI)
			{						
				Type upwindVal = 
					pos(pFaceFlux[faceI])*pSfCorr[faceI] 
				   +(1.0 - pos(pFaceFlux[faceI]))*btsfBD[patchI][faceI];
					
				Type centredVal =  
					weights[faceI]*pSfCorr[faceI] 
				  + (1.0 - weights.boundaryField()[patchI][faceI])
				   *btsfBD[patchI][faceI];		
					
				pSfCorr[faceI] = blendingFactor[faceI]*centredVal
				  + (1.0 - blendingFactor[faceI])*upwindVal;   						
			}
		}
		else if (isA<cyclicFvPatch>(patches[patchI]))	
		{
			const cyclicPolyPatch& cycPatch = 
				refCast<const cyclicPolyPatch>(mesh.boundaryMesh()[patchI]);
		
			if (cycPatch.owner() == true)
			{							
				const label neighbPatchID = cycPatch.neighbPatchID();
				
				const scalarField& pFaceFlux = 
					faceFlux_.boundaryField()[patchI];				
				
				const labelUList& pOwner = mesh.boundary()[patchI].faceCells();	
				const labelUList& pNeigh = 
					mesh.boundary()[neighbPatchID].faceCells();
				
				const label startFaceOwn = patches[patchI].start();
				const label startFaceNeigh = patches[neighbPatchID].start();
			
				forAll(pOwner, faceI)
				{				
					pSfCorr[faceI] = 
						sumFlux
						(
							(**dimList_)[pOwner[faceI]],
							coeffsWeighted[pOwner[faceI]],
							(**intBasTrans_)[faceI + startFaceOwn][0]
						)  /(**refFacAr_)[faceI + startFaceOwn][0];				
				
					Type neighValue = 
						sumFlux
						(
							(**dimList_)[pNeigh[faceI]],
							coeffsWeighted[pNeigh[faceI]],
							(**intBasTrans_)[faceI + startFaceNeigh][0]
						)  /(**refFacAr_)[faceI + startFaceNeigh][0];				

					Type upwindVal = 
						pos(pFaceFlux[faceI])*pSfCorr[faceI] 
					   +(1.0 - pos(pFaceFlux[faceI]))*neighValue;
						
					Type centredVal =  
						weights[faceI]*pSfCorr[faceI] 
					  + (1.0 - weights.boundaryField()[patchI][faceI])
					   *neighValue;		
						
					pSfCorr[faceI] = blendingFactor[faceI]*centredVal 
					  + (1.0 - blendingFactor[faceI])*upwindVal; 	
				
					btsf[neighbPatchID][faceI] = pSfCorr[faceI]; 				
				}				
			}
		}		
		
		
	}		
}

// ************************************************************************* //
