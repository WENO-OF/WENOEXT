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
#include "version.h"
#include "WENOCoeff.H"
#include "DynamicField.H"
#include "processorFvPatch.H"


// Generate macro for forAll expansion of unsigned integers
#define forAllU(list,i) \
    for (unsigned int i=0; i<(list).size();i++)

// * * * * * * * * * * * * * *  Static Variables * * * * * * * * * * * * * * //
template<class Type>
bool Foam::WENOCoeff<Type>::printWENODict_=false;

template<class Type>
scalar Foam::WENOCoeff<Type>::dm_=1000;

template<class Type>
scalar Foam::WENOCoeff<Type>::epsilon_=1E-40;

template<class Type>
scalar Foam::WENOCoeff<Type>::p_=4;

// * * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * //

template<class Type>
Foam::WENOCoeff<Type>::WENOCoeff
(
    const fvMesh& mesh,
    const label polOrder
)
:
    mesh_(mesh),
    polOrder_(polOrder),
    WENOBase_
    (
        WENOBase::instance
        (
            mesh,
            polOrder_
        )
    ),
    nDvt_(WENOBase_.degreesOfFreedom())
{
    bJ_.reserve(100);
    
    if (!printWENODict_)
    {
        // Read expert factors
        IOdictionary WENODict
        (
            IOobject
            (
                "WENODict",
                mesh.time().caseSystem(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );


        p_ = WENODict.lookupOrAddDefault<scalar>("p", 4.0);
        dm_ = WENODict.lookupOrAddDefault<scalar>("dm", 1000.0);
        epsilon_ = WENODict.lookupOrAddDefault<scalar>("epsilon",1E-40);
        
        // Add pol order for printing 
        WENODict.add<label>("polOrder",polOrder_);

        Info << "WENO Version: "<< word(GIT_BRANCH)<<" "<<word(GIT_COMMIT_HASH) << nl
             << "WENODict:"
             << WENODict << endl;
        printWENODict_=true;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void Foam::WENOCoeff<Type>::calcCoeff
(
    const label cellI,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelListList& stencilList,
    DynamicList<coeffType>& coeffsList
) const
{
    // Set coefficient size to number of stencils
    coeffsList.setSize(stencilList.size(),coeffType(1,pTraits<Type>::zero));
    
    label coeffIndex = 0;
    

    forAll(stencilList,stencilI)
    {
        if (stencilList[stencilI][0] != int(WENOBase::Cell::deleted))
        {
            const List<label>& stencilsIDI =
                WENOBase_.stencilsID()[cellI][stencilI];
            const blaze::DynamicMatrix<double>& A =
                WENOBase_.LSmatrix()[cellI][stencilI]();
            const List<label>& cellToProcMapI =
                WENOBase_.cellToProcMap()[cellI][stencilI];

            bJ_.resize(A.columns());

            // Calculate degrees of freedom of stencil as a matrix vector product
            // First line is always constraint line

            for (label j = 1; j < stencilsIDI.size(); j++)
            {

                // Distinguish between local and halo cells
                if (cellToProcMapI[j] == int(WENOBase::Cell::local))
                {
                    bJ_[j-1] = vf[stencilsIDI[j]] - vf[cellI];
                }
                else if(cellToProcMapI[j] != int(WENOBase::Cell::deleted))
                {
                    bJ_[j-1] =
                        haloData_[cellToProcMapI[j]][stencilsIDI[j]]
                      - vf[cellI];
                }
            }
            
            // calculate coefficients
            coeffsList[coeffIndex] = A*bJ_;
            coeffIndex++;
        }
    }
    // Correct for deleted stencils
    coeffsList.setSize(coeffIndex);
}


template<class Type>
void Foam::WENOCoeff<Type>::collectData
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Distribute data to neighbour processors

    haloData_.setSize(WENOBase_.ownHalos().size());

    forAll(haloData_, procI)
    {
        haloData_[procI].setSize(WENOBase_.ownHalos()[procI].size());

        forAll(haloData_[procI], cellI)
        {
            haloData_[procI][cellI] =
                vf.internalField()[WENOBase_.ownHalos()[procI][cellI]];
        }
    }
    
    
    #ifdef FOAM_PSTREAM_COMMSTYPE_IS_ENUMCLASS 
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    #else
        PstreamBuffers pBufs(Pstream::nonBlocking);
    #endif

    // Distribute data
    forAll(WENOBase_.sendProcList(), procI)
    {
        if (WENOBase_.sendProcList()[procI] != -1)
        {
            UOPstream toBuffer(WENOBase_.sendProcList()[procI], pBufs);
            toBuffer << haloData_[procI];
        }
    }

    pBufs.finishedSends();

    // Collect data

    forAll(WENOBase_.receiveProcList(), procI)
    {
        haloData_[procI].clear();
        if (WENOBase_.receiveProcList()[procI] != -1)
        {
            UIPstream fromBuffer(WENOBase_.receiveProcList()[procI], pBufs);
            fromBuffer >> haloData_[procI];
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Foam::Field<Type> > >
Foam::WENOCoeff<Type>::getWENOPol
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (Pstream::parRun())
        collectData(vf);


    // Runtime operations
    
    tmp<Field<Field<Type> > > coeffsWeightedTmp
    (
        new Field<Field<Type> >
        (
            mesh_.nCells(),
            Field<Type>(nDvt_,pTraits<Type>::zero)
        )
    );
    
    Field<Field<Type> >& coeffsWeighted = coeffsWeightedTmp.ref();
    
    DynamicList<coeffType> coeffsI(10,coeffType(1,pTraits<Type>::zero));


    for (label cellI = 0; cellI < mesh_.nCells(); cellI++)
    {
        const labelListList& stencilList = WENOBase_.stencilsID()[cellI];
        
        
        // If no valid stencils are given return zero list of weighted 
        // coefficients
        if (stencilList[0][0] == int(WENOBase::Cell::empty))
            continue;
        
        calcCoeff
        (
            cellI,
            vf,
            stencilList,
            coeffsI
        );
        
        // Get weighted combination
        calcWeight
        (
            coeffsWeighted[cellI],
            cellI,
            vf,
            coeffsI
        );
    }

    return coeffsWeightedTmp;
}


// Specialisation for scalar
template<>
inline void Foam::WENOCoeff<Foam::scalar>::calcWeight
(
    Field<scalar>& coeffsWeightedI,
    const label cellI,
    const GeometricField<scalar, fvPatchField, volMesh>& vf,
    const DynamicList<coeffType>& coeffsList
) const
{
    scalar gamma = 0.0;
    scalar gammaSum = 0.0;

    forAll(coeffsList, stencilI)
    {
        const auto& coeffs = coeffsList[stencilI];
        
        const scalar smoothInd = trans(coeffs) * (WENOBase_.B()[cellI]*coeffs);

        // Calculate gamma for central and sectorial stencils

        if (stencilI == 0)
        {
            gamma = dm_/(intPow(epsilon_ + smoothInd,p_));
        }
        else
        {
            gamma = 1.0/(intPow(epsilon_ + smoothInd,p_));
        }

        gammaSum += gamma;

        forAllU(coeffs, coeffI)
        {
            coeffsWeightedI[coeffI] += coeffs[coeffI]*gamma;
        }
    }

    forAll(coeffsWeightedI, coeffI)
    {
        coeffsWeightedI[coeffI] /= gammaSum;
    }       
}



template<class Type>
void Foam::WENOCoeff<Type>::calcWeight
(
    Field<Type>& coeffsWeightedI,
    const label cellI,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const DynamicList<coeffType>& coeffsList
) const 
{
    scalar gamma = 0.0;

    for (label compI = 0; compI < vf[0].size(); compI++)
    {
        scalar gammaSum = 0.0;

        forAll(coeffsList, stencilI)
        {
            const auto& coeffs = coeffsList[stencilI];

            // Get smoothness indicator

            scalar smoothInd = 0.0;
            auto sumB =  WENOBase_.B()[cellI] * coeffs;
            forAllU(coeffs, coeffP)
            {
                smoothInd += coeffs[coeffP][compI]*sumB[coeffP][compI];
            }

            // Calculate gamma for central and sectorial stencils

            if (stencilI == 0)
            {
                gamma = dm_/(intPow(epsilon_ + smoothInd,p_));
            }
            else
            {
                gamma = 1.0/(intPow(epsilon_ + smoothInd,p_));
            }

            gammaSum += gamma;

            forAllU(coeffs, coeffI)
            {
                coeffsWeightedI[coeffI][compI] += coeffs[coeffI][compI]*gamma;
            }
        }

        forAll(coeffsWeightedI, coeffI)
        {
            coeffsWeightedI[coeffI][compI] /= gammaSum;
        }
        
    }
}


template<class Type>
Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& Foam::WENOCoeff<Type>::storeOrRetrieve
(
    const word fieldName
) const
{
    if (!mesh_.objectRegistry::foundObject<GeometricField<Type, fvPatchField, volMesh> >(fieldName))
    {
        GeometricField<Type, fvPatchField, volMesh>* fPtr
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                IOobject
                (
                    fieldName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensioned<Type>("sensor",dimless,Type())
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return mesh_.objectRegistry::lookupObjectRef<GeometricField<Type, fvPatchField, volMesh> >(fieldName);
}


// ************************************************************************* //
