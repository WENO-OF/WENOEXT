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
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef WENOCoeff_H
#define WENOCoeff_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class WENOCoeff Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class WENOCoeff
{
private:

    // Private Data

        //- Typedef for 3D scalar matrix
		typedef List< List< List<scalar> > > scalarMatrix;

		//- WENO weighting factors
		scalar p_;
		scalar dm_;

		//- Dimensionality of the geometry
		//  Individual for each stencil
        labelListList* dimList_;

		//- Order of polynomials
        const scalar polOrder_;

		//- Number of derivatives
        label nDvt_;

		//- Lists of central and sectorial stencil ID's for each cell
		List<labelListList>* stencilsID_;

		//- Information about cell property
		//  Could be extended to boundary conditions
		//  - -1 : local cell
		//  - >-1: halo cell
		List<labelListList>* cellToPatchMap_;

		//- Information about processor neighbours of patches
		labelList* patchToProcMap_;

		//- Lists of surface integrals of basis functions
		//  Calculated in the reference space
		List<List<scalarMatrix> >* intBasTrans_;

		//- List of face areas in the reference space
		List<scalarList>* refFacAr_;

		//- Lists of field values of halo cells
		List<List<Type> > haloData_;

		//- Lists of center coordinates of halo cells
		List<List<point> >* haloCenters_;

		//- Lists of ID's of halo cells
		labelListList* ownHalos_;

		//- Lists of pseudoinverses for each stencil of each cell
		List<List<scalarRectangularMatrix> >* LSmatrix_;

		//- Lists of oscillation matrices for each stencil of each cell
		List<scalarRectangularMatrix>* B_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        WENOCoeff(const WENOCoeff&);

        //- Disallow default bitwise assignment
        void operator=(const WENOCoeff&);

		//- Distribute data if multiple processors are involved
        void collectData(const fvMesh& mesh);


public:

    // Constructor

        WENOCoeff
        (
			const fvMesh& mesh,
			const label polOrder
		)
        :
			polOrder_(polOrder)
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


    // Member Functions

        //- Calling function from different schemes
        Field<Field<Type> > getWENOPol
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        )	;

		//- Calculating the coefficients for each stencil of each cell
		void calcCoeff
		(
			const label cellI,
			const GeometricField<Type, fvPatchField, volMesh>& dataField,
			Field<Type>& dvtI,
			const label stencilI,
			const label excludeStencils
        )	;

		//- Get weighted combination for scalar fields
        void calcWeight
        (
            Field<scalar>& coeffsWeightedI,
            const label cellI,
            const GeometricField<scalar, fvPatchField, volMesh>& vf,
            const Field<Field<scalar> >& coeffsI
        )
        {
			scalar gamma = 0.0;
			scalar gammaSum = 0.0;

			forAll(coeffsI, stencilI)
			{
				const scalarField coeffsIsI = coeffsI[stencilI];

				// Get smoothness indicator

				scalar smoothInd = 0.0;

				forAll(coeffsIsI, coeffP)
				{
					scalar sumB = 0.0;

					forAll(coeffsIsI, coeffQ)
					{
						sumB +=
						    (*B_)[cellI][coeffP][coeffQ]
						   *coeffsIsI[coeffQ];
					}

					smoothInd += coeffsIsI[coeffP]*sumB;
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
					coeffsWeightedI[coeffI] += coeffsIsI[coeffI]*gamma;
				}
			}

			forAll(coeffsWeightedI, coeffI)
			{
				coeffsWeightedI[coeffI] /= gammaSum;
			}
		}

		//- Get weighted combination for non-scalar fields
        void calcWeightComp
        (
            Field<Type>& coeffsWeightedI,
            const label cellI,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const Field<Field<Type> >& coeffsI
        )	;

        void calcWeight
        (
            Field<vector>& coeffsWeightedI,
            const label cellI,
            const GeometricField<vector, fvPatchField, volMesh>& vf,
            const Field<Field<vector> >& coeffsI
        )
        {
			calcWeightComp
			(
				coeffsWeightedI,
				cellI,
				vf,
				coeffsI
			);
		}

        void calcWeight
        (
            Field<tensor>& coeffsWeightedI,
            const label cellI,
            const GeometricField<tensor, fvPatchField, volMesh>& vf,
            const Field<Field<tensor> >& coeffsI
        )
        {
			calcWeightComp
			(
				coeffsWeightedI,
				cellI,
				vf,
				coeffsI
			);
		}

        void calcWeight
        (
            Field<symmTensor>& coeffsWeightedI,
            const label cellI,
            const GeometricField<symmTensor, fvPatchField, volMesh>& vf,
            const Field<Field<symmTensor> >& coeffsI
        )
        {
			calcWeightComp
			(
				coeffsWeightedI,
				cellI,
				vf,
				coeffsI
			);
		}

        void calcWeight
        (
            Field<sphericalTensor>& coeffsWeightedI,
            const label cellI,
            const GeometricField<sphericalTensor, fvPatchField, volMesh>& vf,
            const Field<Field<sphericalTensor> >& coeffsI
        )
        {
			calcWeightComp
			(
				coeffsWeightedI,
				cellI,
				vf,
				coeffsI
			);
		}

		//- Get necessary lists for runtime operations
        inline List<List<scalarMatrix> >** getPointerIntBasTrans()
        {
            return &intBasTrans_;
        };
        inline List<scalarList>** getPointerRefFacAr()
        {
            return &refFacAr_;
        };
        inline labelListList** getPointerDimList()
        {
            return &dimList_;
        };
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
