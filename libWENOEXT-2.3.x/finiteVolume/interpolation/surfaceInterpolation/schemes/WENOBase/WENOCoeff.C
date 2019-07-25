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
#include "WENOBase.H"

#include "processorFvPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void Foam::WENOCoeff<Type>::calcCoeff
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
void Foam::WENOCoeff<Type>::collectData(const fvMesh& mesh)
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
    }
}


template<class Type>
Foam::Field<Foam::Field<Type> >
Foam::WENOCoeff<Type>::getWENOPol
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();
    const fvPatchList& patches = mesh.boundary();

    // Get preprocessing lists from WENOBase class

    Foam::WENOBase& init =
        WENOBase::instance
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
    B_ = init.getPointerB();
    intBasTrans_ = init.getPointerIntBasTrans();
    refFacAr_ = init.getPointerRefFacAr();
    dimList_ = init.getPointerDimList();

    // Distribute data to neighbour processors

    haloData_.setSize((*ownHalos_).size());

    forAll(patches, patchI)
    {
        if (isA<processorFvPatch>(patches[patchI]))
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

    Field<Field<Type> > coeffsWeighted(mesh.nCells());

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        const label nStencilsI = (*LSmatrix_)[cellI].size();
        Field<Field<Type> > coeffsI(nStencilsI);

        coeffsWeighted[cellI].setSize(nDvt_,pTraits<Type>::zero);

        label excludeStencils = 0;
        label stencilI = 0;

        // Calculate degrees of freedom for each stencil of the cell
        while (stencilI < nStencilsI)
        {
            // Offset for deleted stencils
            if ((*stencilsID_)[cellI][stencilI+excludeStencils][0] == -1)
            {
                excludeStencils++;
            }
            else
            {
                coeffsI[stencilI].setSize(nDvt_,pTraits<Type>::zero);

                calcCoeff
                (
                    cellI,
                    vf,
                    coeffsI[stencilI],
                    stencilI,
                    excludeStencils
                );

                stencilI++;
            }
        }

        // Get weighted combination
        calcWeight
        (
            coeffsWeighted[cellI],
            cellI,
            vf,
            coeffsI
        );
    }

    return coeffsWeighted;
}


template<class Type>
void Foam::WENOCoeff<Type>::calcWeightComp
(
    Field<Type>& coeffsWeightedI,
    const label cellI,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const Field<Field<Type> >& coeffsI
)
{
    // Get weighted combination for each component separately

    scalar gamma = 0.0;

    for (label compI = 0; compI < vf[0].size(); compI++)
    {
        scalar gammaSum = 0.0;

        forAll(coeffsI, stencilI)
        {
            const Field<Type> coeffsIsI = coeffsI[stencilI];

            // Get smoothness indicator

            scalar smoothInd = 0.0;

            forAll(coeffsIsI, coeffP)
            {
                scalar sumB = 0.0;

                forAll(coeffsIsI, coeffQ)
                {
                    sumB +=
                    (
                        (*B_)[cellI][coeffP][coeffQ]
                       *coeffsIsI[coeffQ][compI]
                    );
                }

                smoothInd += coeffsIsI[coeffP][compI]*sumB;
            }

            // Calculate gamma for central and sectorial stencils

            if (stencilI == 0)
            {
                gamma = dm_/(pow(10e-6 + smoothInd,p_));
            }
            else
            {
                gamma = 1.0/(pow(10e-6 + smoothInd,p_));
            }

            gammaSum += gamma;

            forAll(coeffsIsI, coeffI)
            {
                coeffsWeightedI[coeffI][compI] +=
                    coeffsIsI[coeffI][compI]*gamma;
            }
        }

        forAll(coeffsWeightedI, coeffI)
        {
            coeffsWeightedI[coeffI][compI] /= gammaSum;
        }
    }
}


// ************************************************************************* //
