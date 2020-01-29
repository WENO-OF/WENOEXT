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
