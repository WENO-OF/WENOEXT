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

\*---------------------------------------------------------------------------*/

#include "codeRules.H"
#include "WENOUpwindFit01.H"
#include "processorFvPatch.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::WENOUpwindFit01<Type>::calcLimiter
(
    const fvMesh& mesh,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const Field<Field<Type> >& coeffsWeighted,
    GeometricField<Type, fvsPatchField, surfaceMesh>& tsfP
)    const
{
    const Field<Type>& vfI = vf.internalField();

    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    const label nComp = pTraits<Type>::nComponents;

    // Evaluate the limiters
    Field<Type> theta(mesh.nCells(),pTraits<Type>::zero);

    const Type maxVfI = pTraits<Type>::one;
    const Type minVfI = pTraits<Type>::zero;

    Type maxP = pTraits<Type>::zero;
    Type minP = pTraits<Type>::zero;

    scalar argMax = 0.0;
    scalar argMin = 0.0;

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        const cell& faces = mesh.cells()[cellI];

        maxP = vfI[cellI];
        minP = vfI[cellI];

        for (label cI = 0; cI < nComp; cI++)
        {
            // Buffer components:
            auto maxPci = component(maxP,cI);
            auto minPci = component(minP,cI);
            
            forAll(faces, fI)
            {
                if (faces[fI] < mesh.nInternalFaces())
                {
                    // Coefficient to compare:
                    // Set to dummy value so it has the correct type
                    auto tsfPci = maxPci;
                    if (cellI == P[faces[fI]])
                    {
                        if (faceFlux_[faces[fI]] > 0)
                        {
                            // See Eq. (3.47) of master thesis                     
                            tsfPci = 
                                    component(tsfP[faces[fI]],cI)
                                  + component(vfI[P[faces[fI]]],cI);
                        }
                        else
                        {
                            tsfPci = component
                            (   vfI[P[faces[fI]]] + 
                                sumFlux
                                (
                                    WENOBase_.dimList()[P[faces[fI]]],
                                    coeffsWeighted[P[faces[fI]]],
                                    WENOBase_.intBasTrans()[faces[fI]][0]
                                )  /WENOBase_.refFacAr()[faces[fI]]
                            ,cI
                            );
                        }
                    }
                    else if (cellI == N[faces[fI]])
                    {
                        if (faceFlux_[faces[fI]] < 0)
                        {
                            // See Eq. (3.47) of master thesis                     
                            tsfPci = 
                                    component(tsfP[faces[fI]],cI)
                                  + component(vfI[cellI],cI);
                        }
                        else
                        {
                            tsfPci = component
                            (   vfI[cellI] + 
                                sumFlux
                                (
                                    WENOBase_.dimList()[cellI],
                                    coeffsWeighted[cellI],
                                    WENOBase_.intBasTrans()[cellI][1]
                                )  /WENOBase_.refFacAr()[faces[fI]]
                            ,cI
                            );
                        }
                    }
                    
                    
                    if (tsfPci > maxPci)
                        maxPci = tsfPci;
                    else if (tsfPci < minPci)
                        minPci = tsfPci;
                }
            }

            if (mag((maxPci - component(vfI[cellI],cI))) < 1E-9)
            {
                argMax = 1.0;
            }
            else
            {
                argMax =
                    mag((component(maxVfI,cI) - component(vfI[cellI],cI))
                   /(maxPci - component(vfI[cellI],cI)));
            }

            if (mag((minPci - component(vfI[cellI],cI))) < 1E-9)
            {
                argMin = 1.0;
            }
            else
            {
                argMin =
                    mag((component(minVfI,cI)- component(vfI[cellI],cI))
                   /(minPci - component(vfI[cellI],cI)));
            }
            setComponent(theta[cellI],cI) = min(min(argMax, argMin), 1.0);
        }
    }

    // Evaluate the limited internal fluxes

    forAll(P, faceI)
    {
        for (label cI = 0; cI < nComp; cI++)
        {
            if (faceFlux_[faceI] > 0)
                setComponent(tsfP[faceI],cI) =
                        component(theta[P[faceI]],cI)
                      * (component(tsfP[faceI],cI));
             else
                setComponent(tsfP[faceI],cI) =
                        component(theta[N[faceI]],cI)
                      * (component(tsfP[faceI],cI));
                  
        }
    }

    forAll(tsfP.boundaryField(), patchI)
    {
        const fvPatchList& patches = mesh.boundary();

        fvsPatchField<Type>& pbtsfP =
        #ifdef FOAM_NEW_GEOMFIELD_RULES
            tsfP.boundaryFieldRef()[patchI];
        #else 
            tsfP.boundaryField()[patchI];
        #endif
        
        if (isA<processorFvPatch>(patches[patchI]))
        {
            const labelUList& pOwner =
                mesh.boundary()[patchI].faceCells();

            forAll(pOwner, faceI)
            {
                label own = pOwner[faceI];

                for (label cI = 0; cI < nComp; cI++)
                {
                    setComponent(pbtsfP[faceI],cI) =
                            component(theta[own],cI)
                          * component(pbtsfP[faceI],cI);
                }
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::WENOUpwindFit01<Type>::correction
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
    
    calcLimiter(mesh,vf,coeffsWeighted,tsfP);

    return tsfCorrP;
}


template<class Type>
Type Foam::WENOUpwindFit01<Type>::sumFlux
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
void Foam::WENOUpwindFit01<Type>::swapData
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
void Foam::WENOUpwindFit01<Type>::coupledRiemannSolver
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
