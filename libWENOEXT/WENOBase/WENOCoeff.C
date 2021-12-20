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
#include "version.h"
#include "WENOCoeff.H"
#include "DynamicField.H"
#include "processorFvPatch.H"

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
    coeffsList.setSize(stencilList.size());
    
    label coeffIndex = 0;
    

    forAll(stencilList,stencilI)
    {
        if (stencilList[stencilI][0] != int(WENOBase::Cell::deleted))
        {
            const List<label>& stencilsIDI =
                WENOBase_.stencilsID()[cellI][stencilI];
            const auto& A =
                WENOBase_.LSmatrix()[cellI][stencilI]();
            const List<label>& cellToProcMapI =
                WENOBase_.cellToProcMap()[cellI][stencilI];

            // Storage for bJ matrix needed in calcCoeff
            // Note: resize also pre reserves the space by default, see blaze wiki
            const label nComp = pTraits<Type>::nComponents;
            bJ_.resize(A.columns(),nComp);

            // Calculate degrees of freedom of stencil as a matrix vector product
            // First line is always constraint line

            for (label j = 1; j < stencilsIDI.size(); j++)
            {
                // Distinguish between local and halo cells
                if (cellToProcMapI[j] == int(WENOBase::Cell::local))
                {
                    // Loop over the components
                    for (label compI = 0; compI < nComp; compI++)
                        bJ_(j-1,compI) = component(vf[stencilsIDI[j]],compI) 
                                        - component(vf[cellI],compI);
                }
                else if(cellToProcMapI[j] != int(WENOBase::Cell::deleted))
                {
                    // Loop over the components
                    for (label compI = 0; compI < nComp; compI++)
                        bJ_(j-1,compI) =
                            component(receiveHaloData_[cellToProcMapI[j]][stencilsIDI[j]],compI)
                          - component(vf[cellI],compI);
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
void Foam::WENOCoeff<Type>::calcWeight
(
    Field<Type>& coeffsWeightedI,
    const label cellI,
    const DynamicList<coeffType>& coeffsList
) const 
{
    const label nComp = pTraits<Type>::nComponents;

    // Get smoothness indicator matrix B
    const auto& B = WENOBase_.B()[cellI];

    scalar gammaSum[nComp] = {0.0};
    scalar gamma[nComp] = {0.0};

    forAll(coeffsList, stencilI)
    {
        const auto& coeffs = coeffsList[stencilI];
        
        // First multiply B with coefficients
        // See Eq. (6.37)  in [1]
        auto sumB =  B*coeffs;
        for (label compI=0; compI < nComp; compI++)
        {
            // Get reference to column
            const auto colCoeff = column(coeffs,compI);
            const auto colSumB = column(sumB,compI);
            const scalar smoothInd = trans(colCoeff)*colSumB;
            
            if (stencilI == 0)
            {
                gamma[compI] = dm_/(intPow(epsilon_ + smoothInd,p_));
            }
            else
            {
                gamma[compI] = 1.0/(intPow(epsilon_ + smoothInd,p_));
            }
            
            gammaSum[compI] += gamma[compI];
         
            forAll(coeffsWeightedI, coeffI)
            {
                setComponent(coeffsWeightedI[coeffI],compI) += gamma[compI] * colCoeff[coeffI];
            }
        }
    }

    forAll(coeffsWeightedI, coeffI)
    {
        for (label compI = 0; compI < nComp; compI++)
            setComponent(coeffsWeightedI[coeffI],compI) /= gammaSum[compI];
    }
}


template<class Type>
void Foam::WENOCoeff<Type>::collectData
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Distribute data to neighbour processors
    receiveHaloData_.setSize(WENOBase_.receiveHaloSize().size());
    sendHaloData_.setSize(WENOBase_.sendHaloCellIDList().size());
    
    
    outstandingRecvRequest_.setSize(WENOBase_.receiveHaloSize().size(),-1);
    
    // store current request index
    const label nReq = UPstream::nRequests();
    
    // This represents initEvaluate of processorFvPatchField.C
    forAll(receiveHaloData_, procI)
    {
        receiveHaloData_[procI].setSize(WENOBase_.receiveHaloSize()[procI]);
        // Make const references for easier access
        const labelList& sendHaloCellIDs = WENOBase_.sendHaloCellIDList()[procI];
        const label sendProcID = WENOBase_.sendProcList()[procI];
        const label receiveProcID = WENOBase_.receiveProcList()[procI];
        sendHaloData_[procI].setSize(sendHaloCellIDs.size());
        
        if (sendProcID == -1 && receiveProcID == -1)
            continue;

        // Fill halo data to send to other processors
        forAll(sendHaloData_[procI], cellI)
        {
            sendHaloData_[procI][cellI] =
                vf.internalField()[sendHaloCellIDs[cellI]];
        }
        
        if (receiveProcID != -1)
        {
            outstandingRecvRequest_[procI] = UPstream::nRequests();
            // UIPstream from processorFvPatchField.C
            UIPstream::read
            (
                Pstream::commsTypes::nonBlocking,
                receiveProcID,
                reinterpret_cast<char*>(receiveHaloData_[procI].data()), // The data to read into
                receiveHaloData_[procI].byteSize(),
                UPstream::msgType(),   // this is UPstream::msgType() from processorFvPatch.H
                mesh_.comm()   // this is the communicator stored e.g. in the mesh object
            );
        }
        
        if (sendProcID != -1)
        {
            // UIPstream from processorFvPatchField.C
            UOPstream::write
            (
                Pstream::commsTypes::nonBlocking,
                sendProcID,
                reinterpret_cast<char*>(sendHaloData_[procI].data()), // The data to read into
                sendHaloData_[procI].byteSize(),
                UPstream::msgType(),   // this is UPstream::msgType() from processorFvPatch.H
                mesh_.comm()   // this is the communicator stored e.g. in the mesh object
            );
        }
            
    }
    
    // Note: For large scale simulations it appears that one 
    // MPI_Waitall is better than using MPI_Wait for each request.
    // Also using MPI_Test in calcCoeff() for a real non-blocking communication
    // increased communication time for large number of processors (>2000)
    UPstream::waitRequests(nReq);
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
    
    // Construct list with default 10 elements
    DynamicList<coeffType> coeffsI(10);


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
            coeffsI
        );
    }

    return coeffsWeightedTmp;
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
