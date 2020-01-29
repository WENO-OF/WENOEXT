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
#include "WENOCoeff.H"
#include "DynamicField.H"

#include "processorFvPatch.H"

// * * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * //

template<class Type>
Foam::WENOCoeff<Type>::WENOCoeff
(
    const fvMesh& mesh,
    const label polOrder
)
:
    polOrder_(polOrder),
    WENOBase_
    (
        WENOBase::instance
        (
            mesh,
            polOrder_
        )
    )
{
    // 3D version
    if (mesh.nGeometricD() == 3)
    {
        nDvt_ =
            (polOrder_ + 1.0)*(polOrder_ + 2.0)*(polOrder_ + 3.0)
           /6.0 - 1.0;

        Info<< "Reconstruction using WENO"
            << polOrder_ << " (3D version)" << endl;
    }
    else // 2D version (always only one cell in z-direction)
    {
        nDvt_ = (polOrder_ + 1.0)*(polOrder_ + 2.0)/2.0 - 1.0;

        Info<< "Reconstruction using WENO"
            << polOrder_ << " (2D version)" << endl;
    }

    // Read expert factors
    IOdictionary WENODict
    (
        IOobject
        (
            "WENODict",
            mesh.time().path()/"system",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    p_ = WENODict.lookupOrDefault<scalar>("p", 4.0);
    dm_ = WENODict.lookupOrDefault<scalar>("dm", 1000.0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void Foam::WENOCoeff<Type>::calcCoeff
(
    const label cellI,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    List<Type>& coeff,
    const label stencilI,
    const label excludeStencils
) const
{
    const List<label>& stencilsIDI =
        WENOBase_.stencilsID()[cellI][stencilI + excludeStencils];
    const scalarRectangularMatrix& A =
        WENOBase_.LSmatrix()[cellI][stencilI];
    const List<label>& cellToPatchMapI =
        WENOBase_.cellToPatchMap()[cellI][stencilI + excludeStencils];

    // Calculate degrees of freedom of stencil as a matrix vector product
    // First line is always constraint line
    
    coeff.setSize(nDvt_,pTraits<Type>::zero);

    Type bJ = pTraits<Type>::zero;


    for (label j = 1; j < stencilsIDI.size(); j++)
    {

        // Distinguish between local and halo cells
        if (cellToPatchMapI[j] == -1)
        {
            bJ = vf[stencilsIDI[j]] - vf[cellI];

            for (label i = 0; i < nDvt_; i++)
            {
                coeff[i] += A[i][j-1]*bJ;
            }
        }
        else if(cellToPatchMapI[j] > -4)
        {
            bJ =
                haloData_[cellToPatchMapI[j]][stencilsIDI[j]]
              - vf[cellI];

            for (label i = 0; i < nDvt_; i++)
            {
                coeff[i] += A[i][j-1]*bJ;
            }
        }
    }
}


template<class Type>
void Foam::WENOCoeff<Type>::collectData(const fvMesh& mesh) const
{
#ifdef FOAM_PSTREAM_COMMSTYPE_IS_ENUMCLASS 
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
#else
    PstreamBuffers pBufs(Pstream::nonBlocking);
#endif

    // Distribute data
    forAll(WENOBase_.patchToProcMap(), patchI)
    {
        if (WENOBase_.patchToProcMap()[patchI] > -1)
        {
            UOPstream toBuffer(WENOBase_.patchToProcMap()[patchI], pBufs);
            toBuffer << haloData_[patchI];
        }
    }

    pBufs.finishedSends();

    // Collect data

    forAll(WENOBase_.patchToProcMap(), patchI)
    {
        if (WENOBase_.patchToProcMap()[patchI] > -1)
        {
            haloData_[patchI].clear();

            UIPstream fromBuffer(WENOBase_.patchToProcMap()[patchI], pBufs);
            fromBuffer >> haloData_[patchI];
        }
    }
}


template<class Type>
Foam::Field<Foam::Field<Type> >
Foam::WENOCoeff<Type>::getWENOPol
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = vf.mesh();
    const fvPatchList& patches = mesh.boundary();

    // Distribute data to neighbour processors

    haloData_.setSize(WENOBase_.ownHalos().size());

    forAll(patches, patchI)
    {
        if (isA<processorFvPatch>(patches[patchI]))
        {
            forAll(haloData_, patchI)
            {
                haloData_[patchI].setSize(WENOBase_.ownHalos()[patchI].size());

                forAll(haloData_[patchI], cellI)
                {
                    haloData_[patchI][cellI] =
                        vf.internalField()[WENOBase_.ownHalos()[patchI][cellI]];
                }
            }
        }
    }

    collectData(mesh);



    // Runtime operations
    
    Field<Field<Type> > coeffsWeighted(mesh.nCells());
    
    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        coeffsWeighted[cellI].setSize(nDvt_,pTraits<Type>::zero);

        const label nStencilsI = WENOBase_.LSmatrix()[cellI].size();
        
        List<List<Type> > coeffsI(nStencilsI);
        
	    label excludeStencils = 0;
        label stencilI = 0;


        // Calculate degrees of freedom for each stencil of the cell
        while (stencilI < nStencilsI)
        {
            // Offset for deleted stencils
            if (WENOBase_.stencilsID()[cellI][stencilI+excludeStencils][0] == -1)
            {
                excludeStencils++;
            }
            else
            {
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
void Foam::WENOCoeff<Type>::calcWeight
(
    Field<scalar>& coeffsWeightedI,
    const label cellI,
    const GeometricField<scalar, fvPatchField, volMesh>& vf,
    const List<List<scalar> >& coeffsI
) const
{
    scalar gamma = 0.0;
    scalar gammaSum = 0.0;

    forAll(coeffsI, stencilI)
    {
        const List<scalar> coeffsIsI = coeffsI[stencilI];

        // Get smoothness indicator

        scalar smoothInd = 0.0;

        forAll(coeffsIsI, coeffP)
        {
            scalar sumB = 0.0;

            forAll(coeffsIsI, coeffQ)
            {
                sumB +=
                    WENOBase_.B()[cellI][coeffP][coeffQ]
                   *coeffsIsI[coeffQ];
            }

            smoothInd += coeffsIsI[coeffP]*sumB;
        }

        // Calculate gamma for central and sectorial stencils

        if (stencilI == 0)
        {
            gamma = dm_/(pow(epsilon_ + smoothInd,p_));
        }
        else
        {
            gamma = 1.0/(pow(epsilon_ + smoothInd,p_));
        }

        gammaSum += gamma;

        forAll(coeffsIsI, coeffI)
        {
            coeffsWeightedI[coeffI] += coeffsIsI[coeffI]*gamma;
        }
    }

    forAll(coeffsWeightedI, coeffI)
    {
        coeffsWeightedI[coeffI] /= gammaSum;
    }
}



template<class Type>
void Foam::WENOCoeff<Type>::calcWeightComp
(
    Field<Type>& coeffsWeightedI,
    const label cellI,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const List<List<Type> >& coeffsI
) const 
{
    // Get weighted combination for each component separately

    scalar gamma = 0.0;

    for (label compI = 0; compI < vf[0].size(); compI++)
    {
        scalar gammaSum = 0.0;

        forAll(coeffsI, stencilI)
        {
            const List<Type>& coeffsIsI = coeffsI[stencilI];

            // Get smoothness indicator

            scalar smoothInd = 0.0;

            forAll(coeffsIsI, coeffP)
            {
                scalar sumB = 0.0;

                forAll(coeffsIsI, coeffQ)
                {
                    sumB +=
                    (
                        WENOBase_.B()[cellI][coeffP][coeffQ]
                       *coeffsIsI[coeffQ][compI]
                    );
                }

                smoothInd += coeffsIsI[coeffP][compI]*sumB;
            }

            // Calculate gamma for central and sectorial stencils

            if (stencilI == 0)
            {
                gamma = dm_/(pow(epsilon_ + smoothInd,p_));
            }
            else
            {
                gamma = 1.0/(pow(epsilon_ + smoothInd,p_));
            }

            gammaSum += gamma;



            forAll(coeffsIsI, coeffI)
            {
                coeffsWeightedI[coeffI][compI] += coeffsIsI[coeffI][compI]*gamma;
            }
        }

        forAll(coeffsWeightedI, coeffI)
        {
            coeffsWeightedI[coeffI][compI] /= gammaSum;
        }
    }
}


// ************************************************************************* //
