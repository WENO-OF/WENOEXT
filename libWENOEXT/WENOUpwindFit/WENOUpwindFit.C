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
#include "WENOUpwindFit.H"
#include "processorFvPatch.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::WENOUpwindFit<Type>::calcLimiter
(
    const fvMesh& mesh,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const Field<Field<Type> >& coeffsWeighted,
    GeometricField<Type, fvsPatchField, surfaceMesh>& tsfP
)    const
{
    // Limit the explicit correction if the polynome at the face would exceed
    // the cell center values. This is similar to cellLimited schemes of 
    // OpenFOAM
    // If the face values are calculated with the polynome p,
    // p = psi + \sum{alpha_k * Omega_k}
    // where \sum{alpha_k * Omega_k} is here the field 'tsfP'.
    // Then psi is the implicit part of the upwind scheme and only tsfP is 
    // limited. 

    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    const label nComp = pTraits<Type>::nComponents;


    // --------------- Calculate theta ----------------------------------------
    // ------------------------------------------------------------------------
    
    GeometricField<Type, fvsPatchField, surfaceMesh> theta
    (
        IOobject
        (
            "theta",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensioned<Type>(vf.name(), vf.dimensions(), pTraits<Type>::zero)
    );

    // Calculate theta for internal field
    forAll(P,faceI)
    {
        // Get the cell center value of this polynome
        // called psi, the neighbour one is called psiN
        Type psi;
        Type psiN;
        if (faceFlux_[faceI] > 0)
        {
            psi  = vf[P[faceI]];
            psiN = vf[N[faceI]];
        }
        else
        {
            psi = vf[N[faceI]];
            psiN = vf[P[faceI]];
        }
        
        // Loop over all components
        for (label cI = 0; cI < nComp; cI++)
        {
            // Check that psi is not larger than the difference of psi to psiN
            if (mag(component(tsfP[faceI],cI)) > mag(component(psi-psiN,cI)))
                setComponent(theta[faceI],cI) = 
                    min(mag(component(psi-psiN,cI))/mag(component(tsfP[faceI],cI)),1.0);
            else
                setComponent(theta[faceI],cI) = 1.0;
        }
    }

    // Calculate theta for boundary:
    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
    #ifdef FOAM_NEW_GEOMFIELD_RULES
        Boundary& btheta = theta.boundaryFieldRef();
    #else 
        GeometricBoundaryField& btheta = theta.boundaryField();
    #endif

    forAll(btheta, patchI)
    {
        if ((btheta[patchI]).coupled())
        {
            const scalarField& pFaceFlux =
                faceFlux_.boundaryField()[patchI];

            const labelUList& pOwner = mesh.boundary()[patchI].faceCells();

            const fvsPatchField<Type>& pbtsfP = tsfP.boundaryField()[patchI];

            // Get patch neighbour field
            const Field<Type>& vfN = (vf.boundaryField()[patchI].patchNeighbourField())();

            fvsPatchField<Type>& pbtheta = btheta[patchI];

            forAll(pOwner, faceI)
            {
                label own = pOwner[faceI];

                // Get the cell center value of this polynome
                // called psi, the neighbour one is called psiN
                Type psi;
                Type psiN;
                if (pFaceFlux[faceI] > 0)
                {
                    psi  = vf[own];
                    psiN = vfN[faceI];
                }
                else
                {
                    psi = vfN[faceI];
                    psiN = vf[own];
                }
                
                // Loop over all components
                for (label cI = 0; cI < nComp; cI++)
                {
                    // Check that psi is not larger than the difference of psi to psiN
                    if (mag(component(pbtsfP[faceI],cI)) > mag(component(psi-psiN,cI)))
                        setComponent(pbtheta[faceI],cI) = 
                            min(mag(component(psi-psiN,cI))/mag(component(pbtsfP[faceI],cI)),1.0);
                    else
                        setComponent(pbtheta[faceI],cI) = 1.0;
                }
            }
        }
    }

    // ---------------------- Adjust tsfP field -------------------------------
    // ------------------------------------------------------------------------
    
    // Evaluate the limited internal fluxes

    forAll(P, faceI)
    {
        for (label cI = 0; cI < nComp; cI++)
        {
            setComponent(tsfP[faceI],cI) =
                component(theta[faceI],cI) * component(tsfP[faceI],cI);
        }
    }

    // Adjust boundary
    forAll(tsfP.boundaryField(), patchI)
    {
        const fvPatchList& patches = mesh.boundary();

        fvsPatchField<Type>& pbtsfP =
        #ifdef FOAM_NEW_GEOMFIELD_RULES
            tsfP.boundaryFieldRef()[patchI];
        #else 
            tsfP.boundaryField()[patchI];
        #endif
        
        const fvsPatchField<Type>& pbtheta = theta.boundaryField()[patchI];
        
        if (isA<processorFvPatch>(patches[patchI]))
        {
            const labelUList& pOwner = mesh.boundary()[patchI].faceCells();
            
            forAll(pOwner, faceI)
            {
                for (label cI = 0; cI < nComp; cI++)
                {
                    setComponent(pbtsfP[faceI],cI) =
                        component(pbtheta[faceI],cI) * component(pbtsfP[faceI],cI);
                }
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::WENOUpwindFit<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // Get degrees of freedom from WENOCoeff class
    tmp<Field<Field<Type> > > coeffsWeightedTmp = WENOCoeff_.getWENOPol(vf);
    const Field<Field<Type> >& coeffsWeighted = coeffsWeightedTmp();


    // Calculate the interpolated face values
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorrP
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "tvfP",
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
    GeometricField<Type, fvsPatchField, surfaceMesh>& tsfP =
    #ifdef FOAM_NEW_TMP_RULES
        tsfCorrP.ref();
    #else 
        tsfCorrP();
    #endif


    // Exact Riemann solver at each internal and coupled face
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
        }
        else
        {
            tsfP[faceI] = pTraits<Type>::zero;
        }
    }
    
    coupledRiemannSolver(mesh, tsfP, vf, coeffsWeighted);
    
    if (limFac_)
        calcLimiter(mesh,vf,coeffsWeighted,tsfP);

    return tsfCorrP;
}


template<class Type>
Type Foam::WENOUpwindFit<Type>::sumFlux
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
void Foam::WENOUpwindFit<Type>::swapData
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
void Foam::WENOUpwindFit<Type>::coupledRiemannSolver
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
