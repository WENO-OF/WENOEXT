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
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "codeRules.H"
#include "WENOUpwindFit.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::WENOUpwindFit<Type>::calcLimiter
(
    const fvMesh& mesh,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, fvsPatchField, surfaceMesh>& tsfP,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& tsfN
)    const
{
    const Field<Type>& vfI = vf.internalField();

    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    const label nComp = vfI[0].size();

    // Evaluate the limiters

    Field<Type> theta(mesh.nCells(),pTraits<Type>::zero);

    const Type maxPhi = max(vfI);
    const Type minPhi = min(vfI);

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
            forAll(faces, fI)
            {
                if (faces[fI] < mesh.nInternalFaces())
                {
                    if (cellI == P[faces[fI]])
                    {
                        if (tsfP[faces[fI]][cI] > maxP[cI])
                        {
                            maxP[cI] = tsfP[faces[fI]][cI];
                        }
                        else if (tsfP[faces[fI]][cI] < minP[cI])
                        {
                            minP[cI] = tsfP[faces[fI]][cI];
                        }
                    }
                    else
                    {
                        if (tsfN[faces[fI]][cI] > maxP[cI])
                        {
                            maxP[cI] = tsfN[faces[fI]][cI];
                        }
                        else if (tsfN[faces[fI]][cI] < minP[cI])
                        {
                            minP[cI] = tsfN[faces[fI]][cI];
                        }
                    }
                }
            }

            if (mag(maxP[cI] - vfI[cellI][cI]) < 1e-10)
            {
                argMax = 1.0;
            }
            else
            {
                argMax =
                    mag((maxPhi[cI] - vfI[cellI][cI])
                   /(maxP[cI] - vfI[cellI][cI]));
            }

            if (mag(minP[cI] - vfI[cellI][cI]) < 1e-10)
            {
                argMin = 1.0;
            }
            else
            {
                argMin =
                    mag((minPhi[cI] - vfI[cellI][cI])
                   /(minP[cI] - vfI[cellI][cI]));
            }

            theta[cellI][cI] = min(min(argMax, argMin), 1.0);
        }
    }

    // Evaluate the limited internal fluxes

    forAll(P, faceI)
    {
        if (faceFlux_[faceI] > 0)
        {
            for (label cI = 0; cI < nComp; cI++)
            {
                tsfP[faceI][cI] =
                    limFac_*(theta[P[faceI]][cI]
                   *(tsfP[faceI][cI] - vfI[P[faceI]][cI])
                  + vfI[P[faceI]][cI])
                  + (1.0 - limFac_)*tsfP[faceI][cI];

                tsfP[faceI][cI] -= vfI[P[faceI]][cI];
            }
        }
        else if (faceFlux_[faceI] < 0)
        {
            for (label cI = 0; cI < nComp; cI++)
            {
                tsfP[faceI][cI] =
                    limFac_*(theta[N[faceI]][cI]
                   *(tsfN[faceI][cI] - vfI[N[faceI]][cI])
                  + vfI[N[faceI]][cI])
                  + (1.0 - limFac_)*tsfN[faceI][cI];

                tsfP[faceI][cI] -= vfI[N[faceI]][cI];
            }
        }
        else
        {
            tsfP[faceI] = pTraits<Type>::zero;
        }
    }

    forAll(tsfP.boundaryField(), patchI)
    {
        fvsPatchField<Type>& pbtsfP =
#ifdef FOAM_NEW_GEOMFIELD_RULES
            tsfP.boundaryFieldRef()[patchI];
#else 
            tsfP.boundaryField()[patchI];
#endif
        const fvsPatchField<Type>& pbtsfN =
            tsfN.boundaryField()[patchI];

        if (pbtsfP.coupled())
        {
            const labelUList& pOwner =
                mesh.boundary()[patchI].faceCells();

            const scalarField& pFaceFlux =
                faceFlux_.boundaryField()[patchI];

            const List<Type>& vfN =
                vf.boundaryField()[patchI].patchNeighbourField();

            forAll(pOwner, faceI)
            {
                label own = pOwner[faceI];

                if (pFaceFlux[faceI] > 0)
                {
                    for (label cI = 0; cI < nComp; cI++)
                    {
                        pbtsfP[faceI][cI] =
                            limFac_*(theta[own][cI]
                           *(pbtsfP[faceI][cI] - vfI[own][cI])
                          + vfI[own][cI]) + (1.0 - limFac_)
                           *pbtsfP[faceI][cI];

                        pbtsfP[faceI][cI] -= vfI[own][cI];
                    }
                }
                else if (pFaceFlux[faceI] < 0)
                {
                    for (label cI = 0; cI < nComp; cI++)
                    {
                        pbtsfP[faceI][cI] =                             // unlimited
                            limFac_*(1.0*(pbtsfN[faceI][cI]
                          - vfN[faceI][cI])
                          + vfN[faceI][cI]) + (1.0 - limFac_)
                           *pbtsfN[faceI][cI];

                        pbtsfP[faceI][cI] -= vfN[faceI][cI];
                    }
                }
                else
                {
                    pbtsfP[faceI] = pTraits<Type>::zero;
                }
            }
        }
    }
}



//- Calculating the limiters for scalar fields
template<>
void Foam::WENOUpwindFit<Foam::scalar>::calcLimiter
(
    const fvMesh& mesh,
    const volScalarField& vf,
    surfaceScalarField& tsfP,
    const surfaceScalarField& tsfN
)    const
{
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    const scalarField& vfI = vf.internalField();

    scalarField theta(mesh.nCells(),0.0);

    scalar maxP = 0.0;
    scalar minP = 0.0;
    const scalar maxPhi = max(vfI);
    const scalar minPhi = min(vfI);
    scalar argMax = 0.0;
    scalar argMin = 0.0;

    // Evaluate the internal limiters

    for (label cellI = 0; cellI< mesh.nCells(); cellI++)
    {
        const cell& faces = mesh.cells()[cellI];

        maxP = vfI[cellI];
        minP = vfI[cellI];

        forAll(faces, fI)
        {
            if (faces[fI] < mesh.nInternalFaces())
            {
                if (cellI == P[faces[fI]])
                {
                    if (tsfP[faces[fI]] > maxP)
                    {
                        maxP = tsfP[faces[fI]];
                    }
                    else if (tsfP[faces[fI]] < minP)
                    {
                        minP = tsfP[faces[fI]];
                    }
                }
                else
                {
                    if (tsfN[faces[fI]] > maxP)
                    {
                        maxP = tsfN[faces[fI]];
                    }
                    else if (tsfN[faces[fI]] < minP)
                    {
                        minP = tsfN[faces[fI]];
                    }
                }
            }
        }

        if (mag(maxP - vfI[cellI]) < 1e-10)
        {
            argMax = 1.0;
        }
        else
        {
            argMax = mag((maxPhi - vfI[cellI])/(maxP - vfI[cellI]));
        }

        if (mag(minP - vfI[cellI]) < 1e-10)
        {
            argMin = 1.0;
        }
        else
        {
            argMin = mag((minPhi - vfI[cellI])/(minP - vfI[cellI]));
        }

        theta[cellI] = min(min(argMax, argMin), 1.0);
    }

    // Evaluate the limited fluxes

    forAll(P, faceI)
    {
        if (faceFlux_[faceI] > 0)
        {
            tsfP[faceI] =
                limFac_*(theta[P[faceI]]*(tsfP[faceI] - vfI[P[faceI]])
              + vfI[P[faceI]]) + (1.0 - limFac_)*tsfP[faceI];

            tsfP[faceI] -= vfI[P[faceI]];
        }
        else if (faceFlux_[faceI] < 0)
        {
            tsfP[faceI] =
                limFac_*(theta[N[faceI]]*(tsfN[faceI] -    vfI[N[faceI]])
              + vfI[N[faceI]]) + (1.0 - limFac_)*tsfN[faceI];

            tsfP[faceI] -= vfI[N[faceI]];
        }
        else
        {
            tsfP[faceI] =  0.0;
        }
    }

    forAll(tsfP.boundaryField(), patchI)
    {
        fvsPatchField<scalar>& pbtsfP =
#ifdef FOAM_NEW_GEOMFIELD_RULES
            tsfP.boundaryFieldRef()[patchI];
#else 
            tsfP.boundaryField()[patchI];
#endif
        const fvsPatchField<scalar>& pbtsfN = tsfN.boundaryField()[patchI];

        if (tsfP.boundaryField()[patchI].coupled())
        {
            const labelUList& pOwner = mesh.boundary()[patchI].faceCells();

            const scalarField& pFaceFlux =
                faceFlux_.boundaryField()[patchI];

            const scalarList& vfN =
                vf.boundaryField()[patchI].patchNeighbourField();

            forAll(pOwner, faceI)
            {
                label own = pOwner[faceI];

                if (pFaceFlux[faceI] > 0)
                {
                    pbtsfP[faceI] =
                        limFac_*(theta[own]*(pbtsfP[faceI] - vfI[own])
                      + vfI[own]) + (1.0 - limFac_)*pbtsfP[faceI];

                    pbtsfP[faceI] -= vfI[own];
                }
                else if (pFaceFlux[faceI] < 0)
                {
                    pbtsfP[faceI] =
                        limFac_*(1.0*(pbtsfN[faceI] - vfN[faceI])                 // unlimited
                      + vfN[faceI]) + (1.0 - limFac_)*pbtsfN[faceI];

                    pbtsfP[faceI] -= vfN[faceI];
                }
                else
                {
                    pbtsfP[faceI] = 0.0;
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
    
    Field<Field<Type> > coeffsWeighted = WENOCoeff_.getWENOPol(vf);


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


    // Unlimited polynomial
    if (limFac_ == 0)
    {
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
                    ) / WENOBase_.refFacAr()[faceI][0];
            }
            else if (faceFlux_[faceI] < 0)
            {
                tsfP[faceI] =
                    sumFlux
                    (
                        WENOBase_.dimList()[N[faceI]],
                        coeffsWeighted[N[faceI]],
                        WENOBase_.intBasTrans()[faceI][1]
                    )  /WENOBase_.refFacAr()[faceI][1];
            }
            else
            {
                tsfP[faceI] = pTraits<Type>::zero;
            }
        }

        coupledRiemannSolver(mesh, tsfP, vf, coeffsWeighted);
    }
    // Limited polynomials
    else
    {
        const fvPatchList& patches = mesh.boundary();

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorrN
            (
                new GeometricField<Type, fvsPatchField, surfaceMesh>
                (
                    IOobject
                    (
                        "tsfCorrN",
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
        GeometricField<Type, fvsPatchField, surfaceMesh>& tsfN =
#ifdef FOAM_NEW_TMP_RULES
        tsfCorrN.ref();
#else 
        tsfCorrN();
#endif

        typename GeometricField<Type, fvsPatchField, surfaceMesh>::
#ifdef FOAM_NEW_GEOMFIELD_RULES
            Boundary& btsfN =
#else 
            GeometricBoundaryField& btsfN =
#endif
#ifdef FOAM_NEW_GEOMFIELD_RULES
        tsfN.boundaryFieldRef();
#else 
        tsfN.boundaryField();
#endif

        typename GeometricField<Type, fvsPatchField, surfaceMesh>::
#ifdef FOAM_NEW_GEOMFIELD_RULES
        Boundary& btsfP = tsfP.boundaryFieldRef();
#else 
        GeometricBoundaryField& btsfP = tsfP.boundaryField();
#endif

        // Calculating face fluxes from both sides

        forAll(P, faceI)
        {
            tsfP[faceI] =
                vf[P[faceI]] + sumFlux
                (
                    WENOBase_.dimList()[P[faceI]],
                    coeffsWeighted[P[faceI]],
                    WENOBase_.intBasTrans()[faceI][0]
                )  /WENOBase_.refFacAr()[faceI][0];

            tsfN[faceI] =
                vf[N[faceI]] + sumFlux
                (
                    WENOBase_.dimList()[N[faceI]],
                    coeffsWeighted[N[faceI]],
                    WENOBase_.intBasTrans()[faceI][1]
                )  /WENOBase_.refFacAr()[faceI][1];
        }

        forAll(btsfN, patchI)
        {
            fvsPatchField<Type>& pbtsfP = btsfP[patchI];
            fvsPatchField<Type>& pbtsfN = btsfN[patchI];

            if (isA<processorFvPatch>(patches[patchI]))
            {
                const labelUList& pOwner = mesh.boundary()[patchI].faceCells();

                label startFace = patches[patchI].start();

                forAll(pOwner, faceI)
                {
                    label own = pOwner[faceI];

                    pbtsfN[faceI] =
                        vf[own] + sumFlux
                        (
                            WENOBase_.dimList()[own],
                            coeffsWeighted[own],
                            WENOBase_.intBasTrans()[faceI + startFace][0]
                        )  /WENOBase_.refFacAr()[faceI + startFace][0] ;

                    pbtsfP[faceI] = pbtsfN[faceI];
                }
            }
        }

        swapData(mesh, btsfN);

        // Limiting the polynomials and evaluating the upwind fluxes

        calcLimiter(mesh, vf, tsfP, tsfN);
    }

    return tsfCorrP;
}


template<class Type>
Type Foam::WENOUpwindFit<Type>::sumFlux
(
    const labelList& dim,
    const Field<Type>& coeffcI,
    const volIntegralType intBasiscIfI
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
                        coeffcI[nCoeff]*intBasiscIfI[n][m][l];

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
            )   const
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
    Field<Field<Type> > coeffsWeighted
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

        if (isA<processorFvPatch>(patches[patchI]))
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
                        )  /WENOBase_.refFacAr()[faceI + startFace][0] ;

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
    }

}


// ************************************************************************* //
