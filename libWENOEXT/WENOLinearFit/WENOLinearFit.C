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

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "codeRules.H"
#include "WENOLinearFit.H"
#include "processorFvPatch.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::WENOLinearFit<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Fatal Error if the correction vector has not yet been set
    if (!tsfCorrP_.valid())
        FatalError << "Explicit correction in WENO was not calculated"
                   << exit(FatalError);


    auto& tsfP = tsfCorrP_.ref();
    // Loop over the correction vector and set the values
    forAll(requireExplicitCorrection_(),faceI)
    {
        // If it does not require a correction, it is set to zero
        if (requireExplicitCorrection_()[faceI] == 0.0)
        {
            tsfP[faceI] = pTraits<Type>::zero;
        }
    }
    // Reset to nullptr to trigger if correction is used before weights is called
    requireExplicitCorrection_.reset(nullptr);
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfCorrTmp(tsfCorrP_.release());
    return tsfCorrTmp;
}


//- Return the interpolation weighting factors for implicit part
template<class Type>
Foam::tmp<surfaceScalarField> Foam::WENOLinearFit<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // Get the linear weights
    surfaceScalarField linearWeights(this->mesh().surfaceInterpolation::weights());

    // Create field for the WENO weights and set to upwind
    tmp<surfaceScalarField> WENOWeightsTmp
    (
        new surfaceScalarField
        (
            IOobject
            (
                "WENOWeights",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pos(faceFlux_)
        )
    );

    surfaceScalarField& WENOWeights = WENOWeightsTmp.ref();




    // Calculate WENO explicit correction vector
    // ========================================= 

    // Get degrees of freedom from WENOCoeff class
    tmp<Field<Field<Type> > > coeffsWeightedTmp = WENOCoeff_.getWENOPol(vf);
    const Field<Field<Type> >& coeffsWeighted = coeffsWeightedTmp();

    // Calculate the interpolated face values
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    requireExplicitCorrection_.reset
    (
        new GeometricField<scalar, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "requireExplicitCorrection_"+vf.name(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<scalar>(vf.name(), dimless, 0.0)
        )
    );

    auto& requireExplicitCorrectionRef = requireExplicitCorrection_.ref();

    tsfCorrP_.reset
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "tvfP_"+vf.name(),
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
    GeometricField<Type, fvsPatchField, surfaceMesh>& tsfP = tsfCorrP_.ref();

    // Evaluate the WENO polynom from the upside cell to the face

    // Get the cell center value of this polynome
    // called psi, the neighbour one is called psiN
    Type psi;   // this is always the upstream value
    Type psiN;  // this is always the downstream neighbor
    forAll(P, faceI)
    {    
        if (faceFlux_[faceI] > 0)
        {
            tsfP[faceI] =
                sumFlux
                (
                    WENOBase_.dimList()[P[faceI]],
                    coeffsWeighted[P[faceI]],
                    WENOBase_.intBasTrans()[faceI][0]
                ) / WENOBase_.refFacAr()[faceI];

            // Calculate the weight
            psi  = vf[P[faceI]];
            psiN = vf[N[faceI]];

            WENOWeights[faceI] = calcWeight(tsfP[faceI],psi,psiN);
            if (WENOWeights[faceI] == 1.0)
                requireExplicitCorrectionRef[faceI] = 1.0;

            // If the WENO weight only differs limFac_ fraction 
            // -- set to linear weight
            if (mag(WENOWeights[faceI]-linearWeights[faceI]) < (1.0-limFac_))
                WENOWeights[faceI] = linearWeights[faceI];
        }
        else if (faceFlux_[faceI] < 0)
        {
            tsfP[faceI] =
                sumFlux
                (
                    WENOBase_.dimList()[N[faceI]],
                    coeffsWeighted[N[faceI]],
                    WENOBase_.intBasTrans()[faceI][1]
                )  /WENOBase_.refFacAr()[faceI];

            // Calculate the weight
            psi  = vf[N[faceI]];
            psiN = vf[P[faceI]];

            // Note: Here a weight of 0 denotes the weighting to the upstream
            WENOWeights[faceI] = 1.0-calcWeight(tsfP[faceI],psi,psiN);
            if (WENOWeights[faceI] == 0.0)
                requireExplicitCorrectionRef[faceI] = 1.0;

            // If the WENO weight only differs limFac_ fraction 
            // -- set to linear weight
            if (mag(WENOWeights[faceI]-linearWeights[faceI]) < (1.0-limFac_))
                WENOWeights[faceI] = linearWeights[faceI];
        }
        else
        {
            tsfP[faceI] = pTraits<Type>::zero;
        }
    }
    
    coupledRiemannSolver(mesh, tsfP, vf, coeffsWeighted);

    return WENOWeightsTmp;
}


template<class Type>
Foam::scalar Foam::WENOLinearFit<Type>::calcWeight
(
    const Type& corr,
    const Type& psi,    // Always the upstream cell center point
    const Type& psiN    // Always the downstream cell center point
) const
{
    // Set the weight to upwind -- psi_f = (w*psi + (1-w)*psiN)
    Type weight = pTraits<Type>::one;
    scalar meanWeight = 0.0;

    // Loop over the components
    const label nComp = pTraits<Type>::nComponents;
    for (label cI=0; cI < nComp; cI++)
    {
        const scalar delta = component(psiN - psi,cI);
        if (delta < SMALL)
            continue;

        setComponent(weight,cI) = 1.0 - component(corr,cI)/delta;
        meanWeight += component(weight,cI);
    }

    for (label cI=0; cI < nComp; cI++)
    {
        // If the mean weight differs more by 10% set return 1.0 for fully upwind
        // and use explicit correction
        if 
        (
            meanWeight/component(weight,cI) > 1.1
         || meanWeight/component(weight,cI) < 0.9
        ) return 1.0;
    }
    return meanWeight;
}


// Scalar specialization
template<>
Foam::scalar Foam::WENOLinearFit<Foam::scalar>::calcWeight
(
    const scalar& corr,
    const scalar& psi,    // Always the upstream cell center point
    const scalar& psiN    // Always the downstream cell center point
) const
{
    // Set the weight to upwind -- psi_f = (w*psi + (1-w)*psiN)
    scalar weight = 1.0;
    const scalar delta = psiN - psi;
    if (delta < SMALL)
        return 1.0;

    return 1.0 - corr/delta;
}


template<class Type>
Type Foam::WENOLinearFit<Type>::sumFlux
(
    const labelList& dim,
    const Field<Type>& coeffcI,
    const volIntegralType& intBasiscIfI
)    const
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
                        coeffcI[nCoeff]*intBasiscIfI(n,m,l);

                    nCoeff++;
                }
            }
        }
    }

    return flux;
}


template<class Type>
void Foam::WENOLinearFit<Type>::swapData
(
    const fvMesh& mesh,
    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
    #ifdef FOAM_NEW_GEOMFIELD_RULES
        Boundary& btsf
    #else 
        GeometricBoundaryField& btsf
    #endif
) const
{
    const fvPatchList& patches = mesh.boundary();

    #ifdef FOAM_PSTREAM_COMMSTYPE_IS_ENUMCLASS 
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    #else 
        PstreamBuffers pBufs(Pstream::nonBlocking);
    #endif

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
void Foam::WENOLinearFit<Type>::coupledRiemannSolver
(
    const fvMesh& mesh,
    GeometricField<Type, fvsPatchField, surfaceMesh>& tsfP,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const Field<Field<Type> >& coeffsWeighted
)   const
{
    const fvPatchList& patches = mesh.boundary();

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
    #ifdef FOAM_NEW_GEOMFIELD_RULES
        Boundary& btsfP = tsfP.boundaryFieldRef();
    #else 
        GeometricBoundaryField& btsfP = tsfP.boundaryField();
    #endif

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfUDCoupled
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "tsfUDCoupled",
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
    GeometricField<Type, fvsPatchField, surfaceMesh>& tsfUD =
    #ifdef FOAM_NEW_TMP_RULES
        tsfUDCoupled.ref();
    #else 
        tsfUDCoupled();
    #endif

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
    #ifdef FOAM_NEW_GEOMFIELD_RULES
        Boundary& btsfUD = tsfUD.boundaryFieldRef();
    #else 
        GeometricBoundaryField& btsfUD = tsfUD.boundaryField();
    #endif

    forAll(btsfP, patchI)
    {
        fvsPatchField<Type>& pSfCorr = btsfP[patchI];

        // for all coupled patches the first step is the same
        if ((patches[patchI]).coupled())
        {
            const scalarField& pFaceFlux =
                faceFlux_.boundaryField()[patchI];

            const labelUList& pOwner = mesh.boundary()[patchI].faceCells();

            label startFace = patches[patchI].start();

            forAll(pOwner, faceI)
            {
                if (pFaceFlux[faceI] > 0)
                {
                    label own = pOwner[faceI];

                    btsfUD[patchI][faceI] =
                        sumFlux
                        (
                            WENOBase_.dimList()[own],
                            coeffsWeighted[own],
                            WENOBase_.intBasTrans()[faceI + startFace][0]
                        )  /WENOBase_.refFacAr()[faceI + startFace];

                    pSfCorr[faceI] = btsfUD[patchI][faceI];
                }
            }
        }
    }

    swapData(mesh, btsfUD);

    forAll(btsfP, patchI)
    {
        fvsPatchField<Type>& pSfCorr = btsfP[patchI];

        if (isA<processorFvPatch>(patches[patchI]))
        {
            const scalarField& pFaceFlux =
                faceFlux_.boundaryField()[patchI];

            const labelUList& pOwner = mesh.boundary()[patchI].faceCells();

            forAll(pOwner, faceI)
            {
                if (pFaceFlux[faceI] < 0)
                {
                    pSfCorr[faceI] = btsfUD[patchI][faceI];
                }
            }
        }
        else if (isA<cyclicFvPatch>(patches[patchI]))
        {
            // If coupled the value at the face of the neighbour patch can be 
            // used.
            const scalarField& pFaceFlux =
                faceFlux_.boundaryField()[patchI];

            const labelUList& pOwner = mesh.boundary()[patchI].faceCells();

            forAll(pOwner, faceI)
            {
                #ifdef FOAM_NEW_COUPLED_PATCHES
                const label neighbPatchID = refCast<const cyclicFvPatch>
                        (patches[patchI]).nbrPatchID();
                #else 
                const label neighbPatchID = refCast<const cyclicFvPatch>
                        (patches[patchI]).neighbPatchID();
                #endif
                
                
                if (pFaceFlux[faceI] < 0)
                {
                    pSfCorr[faceI] = btsfP[neighbPatchID][faceI];
                }
            }
        }
        else if (isA<cyclicAMIFvPatch>(patches[patchI]))
        {
            /*************************** NOTE *******************************
            * Currently not used as it is not quite clear how 
            * the interpolation will affect the results 
            ****************************************************************/
            //// If coupled the value at the face of the neighbour patch can be 
            //// used.
            //const scalarField& pFaceFlux =
                //faceFlux_.boundaryField()[patchI];

            //const labelUList& pOwner = mesh.boundary()[patchI].faceCells();

            //forAll(pOwner, faceI)
            //{
                //const label neighbPatchID = refCast<const cyclicAMIFvPatch>
                        //(patches[patchI]).neighbPatchID();
                //// inerpolate results to patch neighbour field
                //tmp<Field<Type>> interpField = refCast<const cyclicAMIFvPatch>
                        //(patches[patchI]).interpolate(btsfUD[neighbPatchID]);
                
                //if (pFaceFlux[faceI] < 0)
                //{
                    //pSfCorr[faceI] = btsfP[neighbPatchID][faceI];
                //}
            //}
        }
    }
}


// ************************************************************************* //
