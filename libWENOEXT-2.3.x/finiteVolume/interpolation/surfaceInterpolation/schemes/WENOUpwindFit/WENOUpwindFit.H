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
    Foam::WENOUpwindFit

Description
     WENO interpolation scheme class using an exact Riemann solver. Suitable
     for linearised convection terms.

SourceFiles
    WENOUpwindFit.C

Author
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef WENOUpwindFit_H
#define WENOUpwindFit_H

#include "codeRules.H"
#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class WENOUpwindFit Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class WENOUpwindFit
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data

        //- Typedef for 3D scalar matrix
        typedef List< List< List<scalar> > > scalarMatrix;

        //- Lists of surface integrals of basis functions
        //  Calculated in the reference space
        List<List<scalarMatrix> >** intBasTrans_;

        //- List of face areas in the reference space
        List<scalarList>** refFacAr_;

        //- Dimensionality of the geometry
        //  Individual for each stencil
        labelListList** dimList_;

        //- SurfaceScalarField of U() & Sf()
        const surfaceScalarField& faceFlux_;

        //- Polynomial order
        //  User defined parameter
        const scalar polOrder_;

        //- Limiting factor
        //  User defined:
        //  -  0: not limited
        //  - +1: limited
        const scalar limFac_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        WENOUpwindFit(const WENOUpwindFit&);

        //- Disallow default bitwise assignment
        void operator=(const WENOUpwindFit&);

        //- Distribute data through coupled patches and solve Riemann problem
        void swapData
        (
            const fvMesh& mesh,
            typename GeometricField<Type, fvsPatchField, surfaceMesh>::
#ifdef FOAM_NEW_GEOMFIELD_RULES
                Boundary& btsf
#else 
                GeometricBoundaryField& btsf
#endif
        )   const;

        void coupledRiemannSolver
        (
            const fvMesh& mesh,
            GeometricField<Type, fvsPatchField, surfaceMesh>& tsfP,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            Field<Field<Type> > coeffsWeighted
        )   const;


public:

    //- Runtime type information
    TypeName("WENOUpwindFit");


    // Constructors

        //- Construct from mesh    (not used)
        WENOUpwindFit
        (
            const fvMesh& mesh,
            const label polOrder
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(zeroFlux()),
            polOrder_(polOrder),
            limFac_(0)
        {}

        //- Construct from mesh and Istream
        WENOUpwindFit
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            ) ,
            polOrder_(readScalar(is)),
            limFac_(readScalar(is))
        {}

        //- Construct from mesh, faceFlux and Istream
        WENOUpwindFit
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux),
            polOrder_(readScalar(is)),
            limFac_(readScalar(is))
        {}


    // Member Functions

        //- Return a zero surfaceScalarField
        tmp<surfaceScalarField> zeroFlux() const
        {
            tmp<surfaceScalarField> zeroFlux
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "zeroFlux",
                        this->mesh().time().timeName(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar("zeroFlux",dimless,0.0)
                )
            );
            return zeroFlux;
        }

        //- Return the interpolation weighting factors for implicit part
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            return pos(faceFlux_);
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        )    const ;

        //- Calculating the face flux values
        Type sumFlux
        (
            const labelList& dim,
            const Field<Type>& coeffcI,
            const scalarMatrix intBasiscIfI
        )     const;

        //- Calculating the polynomial limiters
        void calcLimiter
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

                if (tsfP.boundaryField()[patchI].coupled())
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
};


//- Calculating the limiters for scalar fields
template<>
void WENOUpwindFit<scalar>::calcLimiter
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WENOUpwindFit.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
