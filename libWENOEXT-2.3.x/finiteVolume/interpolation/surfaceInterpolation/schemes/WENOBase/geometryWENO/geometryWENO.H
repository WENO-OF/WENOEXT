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

Class
    Foam::geometryWENO

Description
    Class for geometric calculations of WENO reconstructions

SourceFiles
    geometryWENO.C

Author
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef geometryWENO_H
#define geometryWENO_H

#include "triSurface.H"
#include "polyMeshTetDecomposition.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class geometryWENO Declaration
\*---------------------------------------------------------------------------*/

class geometryWENO
{
private:

    //- Typedef for 3D scalar matrix
	typedef List< List< List<scalar> > > scalarMatrix;


public:

    // Member Functions

		//- Evaluate the surface integral using Gaussian quadrature
		static scalar gaussQuad
		(
			const scalar n,
			const scalar m,
			const scalar l,
			const point x0,
			const vector v0,
			const vector v1,
			const vector v2
		);

		//- Build the transformation and calculate the volume integral for
		//- each cell
		static void initIntegrals
		(
			const fvMesh& mesh,
			const label cellI,
			const label polOrder,
			scalarMatrix& Integral,
			scalarRectangularMatrix& JInvI,
			point& refPointI,
			scalar& refDetI
		);

		//- Calculate triangulation of halo faces
		static List<point> getTriFaces
		(
			const fvMesh& mesh,
			const label cellI
		);

		//- Calculate volume integrals of halo cells in reference space of
		//- owner cell
		static scalarMatrix getHaloMoments
		(
			const fvMesh& mesh,
			const point transCenterJ,
			const List<point>& triFaces,
			const label polOrder,
			const scalarRectangularMatrix& JInvI,
			const point refPointI
		);

		//- Calculate the volume integrals of stencil cells in reference space
		//- of owner cell
		static scalarMatrix transformIntegral
		(
			const fvMesh& mesh,
			const label cellJ,
			const point transCenterJ,
			const label polOrder,
			const scalarRectangularMatrix& JInvI,
			const point refPointI,
			const scalar refDetI
		);

		//- Check validity of chosen reference frame
		static bool checkRefFrame
		(
			const scalar x0, const scalar y0, const scalar z0,
			const scalar x1, const scalar y1, const scalar z1,
			const scalar x2, const scalar y2, const scalar z2,
			const scalar x3, const scalar y3, const scalar z3
		);

		//- Calculate inverse of Jacobian matrix
		static scalarRectangularMatrix JacobiInverse
		(
			const scalar x0, const scalar y0, const scalar z0,
			const scalar x1, const scalar y1, const scalar z1,
			const scalar x2, const scalar y2, const scalar z2,
			const scalar x3, const scalar y3, const scalar z3
		);

		//- Calculate determinant of inverse Jacobian matrix
		static scalar determinantJacInv
		(
			const scalar x0, const scalar y0, const scalar z0,
			const scalar x1, const scalar y1, const scalar z1,
			const scalar x2, const scalar y2, const scalar z2,
			const scalar x3, const scalar y3, const scalar z3
		);

		//- Transform an arbitrary point into reference space of owner cell
		static point transformPoint
		(
			const scalarRectangularMatrix& Jinv,
			const point xP,
			const point x0
		);

		//- Evaluate the surface integrals for smoothness indicator
		static scalar gaussQuadB
		(
			const scalar n,
			const scalar m,
			const scalar l,
			const point x0,
			const vector v0,
			const vector v1,
			const vector v2
		);

		//- Calculate the volume integrals for the smoothness indicator matrix
		static scalarMatrix smoothIndIntegrals
		(
			const fvMesh& mesh,
			const label cellI,
			const label polOrder,
			const scalarRectangularMatrix& JacobianI,
			const point refPointI
		);

		//- Calculate the smoothness indicator matrices
		static scalarRectangularMatrix getB
		(
			const fvMesh& mesh,
			const label cellI,
			const label polOrder,
			const label nDvt,
			const scalarRectangularMatrix& JInvI,
			const point refPointI,
			const labelList& dim
		);

		//- Check sign of variable
		static scalar Pos(scalar x);

		//- Calculate factorials of variable
		static scalar Fac(scalar x);

		//- Calculation of surface integrals for convective terms
		static void surfIntTrans
		(
			const fvMesh& mesh,
			const label polOrder,
			const List<scalarMatrix>& volMom,
			const List<scalarRectangularMatrix>& JInv,
			const List<point>& refPoint,
			List<List<scalarMatrix> >& intBasTrans,
			List<scalarList>& refFacAr
		);

		static vector compCheck
		(
			const label n,
			const label m,
			const label l,
			const scalarMatrix& intBasisfI
		);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
