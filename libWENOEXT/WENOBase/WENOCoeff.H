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
    Foam::WENOCoeff

Description
    WENO base class for runtime operations of WENO schemes

SourceFiles
    WENOCoeff.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020
    Tobias Martin, <tobimartin2@googlemail.com>

\*---------------------------------------------------------------------------*/

#ifndef WENOCoeff_H
#define WENOCoeff_H

#include "DynamicField.H"
#include "WENOBase.H"
#include "blaze/Math.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class WENOCoeff Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class WENOCoeff
{
protected:

    // Type definitions
    
        //- Typedef for 3D scalar matrix
        using volIntegralType = List3D<scalar> ;
        
        //- C++11 typedef for squareMatrix
        //  This is used for Jacobian matrix
        using scalarSquareMatrix = SquareMatrix<scalar>;
        
        //- Typedef for coefficient type with blaze object
        using coeffType = blaze::DynamicVector<Type>;


    // Private Data

        //- Store reference to mesh 
        const fvMesh& mesh_;
        
        //- Static variable to print WENODict just once during runtime
        static bool printWENODict_;
        
        //- Order of polynomials
        const unsigned int polOrder_;
        
        //- Reference to WENOBase class
        const WENOBase&  WENOBase_;

        //- Number of derivates
        const label nDvt_;

        //- WENO weighting factors
        static scalar p_;
        static scalar dm_;
        static scalar epsilon_;

    // Allocate storage for dynamic variables

        //- Lists of field values of halo cells
        //  Has to be mutable so getWENOPol is const 
        mutable List<List<Type> > haloData_;

        //- Storage for bJ vector needed in calcCoeff
        mutable blaze::DynamicVector<Type> bJ_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        WENOCoeff(const WENOCoeff&);

        //- Disallow default bitwise assignment
        void operator=(const WENOCoeff&);

        //- Distribute data if multiple processors are involved
        void collectData
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
        
        //- Calculating the coefficients for each stencil of each cell
        void calcCoeff
        (
            const label cellI,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const labelListList& stencilList,
            DynamicList<coeffType>& coeffsList
        ) const;


        //- Get weighted combination for any other type
        virtual void calcWeight
        (
            Field<Type>& coeffsWeightedI,
            const label cellI,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const DynamicList<coeffType>& coeffsI
        ) const;


        // For integer power it is much faster to do an integer multiplication
        // This depends on the compiler used! For portability it is explicitly defined
        // here
        inline scalar intPow(const scalar base, const unsigned int exponent) const
        {
            scalar temp = 1.0;
            for (unsigned int i =0; i<exponent;++i)
            {
                temp *= base;
            }
            return temp;
        };

public:

    // Constructor

        WENOCoeff
        (
            const fvMesh& mesh,
            const label polOrder
        );
        
    // Destructor
        virtual ~WENOCoeff() = default;

    // Accessor to WENO Base 
    
        //- Return reference to WENOBase
        const WENOBase& WENOBaseRef() const
        {
            return WENOBase_;
        }

    // Member Functions

        //- Calling function from different schemes
        tmp<Field<Field<Type> > > getWENOPol
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        
        //- Function to store or retrieve fields from the database 
        GeometricField<Type, fvPatchField, volMesh>& storeOrRetrieve
        (
            const word fieldName
        ) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WENOCoeff.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
