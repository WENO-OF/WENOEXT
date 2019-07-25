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
    Foam::WENOPolynomial

Description
    Cell centred Taylor polynomial for interpolation fitting.
    Can also be used for boundary condition extension.

Author
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef WENOPolynomial_H
#define WENOPolynomial_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class WENOPolynomial Declaration
\*---------------------------------------------------------------------------*/

class WENOPolynomial
{
private:

    //- Typedef for 3D scalar matrix
	typedef List< List< List<scalar> > > scalarMatrix;


public:

    // Member functions

		//- Fill matrix with coefficients
        static void addCoeffs
        (
            scalar* coeffs,
			const label polOrder,
			const labelList& dim,
            const scalarMatrix& x_circ_
        )
        {
			register label curIdx = 0;

			for (label n=0 ; n <= dim[0] ; n++)
			{
				for (label m=0 ; m <= dim[1] ; m++)
				{
					for (label l=0 ; l <= dim[2] ; l++)
					{
						if( (n+m+l) <= polOrder && (n+m+l) > 0 )
						{
							coeffs[curIdx++] = x_circ_[n][m][l];
						}
					}
				}
			}
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
