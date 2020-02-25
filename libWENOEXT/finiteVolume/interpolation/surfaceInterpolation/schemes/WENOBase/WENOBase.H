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
    Foam::WENOBase

Description
    WENO base class for preprocessing operations of WENO schemes

SourceFiles
    WENOBase.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020de>
    Tobias Martin, <tobimartin2@googlemail.com>.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef WENOBase_H
#define WENOBase_H

#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class WENOBase Declaration
\*---------------------------------------------------------------------------*/

class WENOBase
{
    private:

        //- Enumerator for halo cells
        //  halo cells have positive IDs that correspond with their cellID
        enum class Cell
        {
            local    = -1,   // is a local cell, in e.g.: cellToPatchMap
            deleted  = -4    // was deleted in splitStencil
        };


    //- Constructors

        WENOBase
        (
            const fvMesh& mesh,
            const label polOrder
        );

       //- Disallow default bitwise copy construct
       WENOBase(const WENOBase&);

       //- Disallow default bitwise assignment
       WENOBase& operator=(const WENOBase&);


        //- Destructor
        ~WENOBase(){};


    //- Private Data

        //- Typedef for 3D scalar matrix
        using volIntegralType = List< List< List<scalar> > > ;
        
        //- C++11 typedef for squareMatrix
        //  This is used for Jacobian matrix
        using scalarSquareMatrix = SquareMatrix<scalar>;

        //- Path to lists in constant folder
        fileName Dir_;

        //- Dimensionality of the geometry
        //  Individual for each stencil
        labelListList dimList_;

        //- Order of polynomials
        label polOrder_;

        //- Number of derivatives
        label nDvt_;

        //- Lists of central and sectorial stencil ID's for each cell
        List<labelListList> stencilsID_;

        //- Lists of center coordinates of halo cells
        List<List<point> > haloCenters_;

        //- Lists of ID's of halo cells
        labelListList ownHalos_;

        //- Information about cell property
        //  Could be extended to boundary conditions
        //  - -1 : local cell
        //  - >-1: halo cell
        //  - -4:  if the stencil is deleted see splitStencil()
        List<labelListList> cellToPatchMap_;

        //- Information about processor neighbours of patches
        labelList patchToProcMap_;

        //- List of volume integrals over basis functions
        //  Calculated in the reference space
        List<volIntegralType> volIntegralsList_;

        //- Lists of surface integrals of basis functions
        //  Calculated in the reference space stored for the owner and neighbour
        //  side. Owner side is the first entry and neighbour side the second
        List<Pair<volIntegralType> > intBasTrans_;

        //- List of face areas in the reference space
        List<scalar> refFacAr_;

        //- Lists of inverse Jacobians for each cell
        List<scalarSquareMatrix> JInv_;

        //- Lists of determinants of inverse Jacobian for each cell
        List<scalar> refDet_;

        //- Lists of reference points of the reference space for each cell
        List<point> refPoint_;

        //- Lists of pseudoinverses for each stencil of each cell
        List<List<scalarRectangularMatrix> > LSmatrix_;

        //- Lists of oscillation matrices for each stencil of each cell
        List<scalarRectangularMatrix> B_;


    //- Private member functions

        //- Split big central stencil into sectorial stencils
        void splitStencil
        (
            const fvMesh& mesh,
            const label cellI,
            label& nStencilsI
        );

        //- Extend stencil by one layer
        void extendStencils
        (
            const fvMesh& mesh,
            const label cellI,
            labelList& lastNeighboursI,
            label& minStencilSize
        );

        //- Sort list of stencil cells from nearest to farest
        //- and cut to necessary size
        void sortStencil
        (
            const fvMesh& mesh,
            const label cellI,
            const label maxSize
        );

        //- Distribute data between processors
        void distributeStencils
        (
            const fvMesh& mesh,
            labelListList& haloCells,
            List<List<List<point> > >& haloTriFaces
        );

        //- Distribute data between local boundaries
        void distributeLocalStencils
        (
            const fvMesh& mesh,
            labelListList& haloCells,
            List<List<List<point> > >& haloTriFaces
        );

        //- Fill the least squares matrices and calculate the
        //- pseudoinverses for each cell
        scalarRectangularMatrix calcMatrix
        (
            const fvMesh& mesh,
            const label cellI,
            const label stencilI,
            const List<List<List<point> > >& haloTriFaces
        );

        //- Calculate the entries of the least squares matrices
        scalar calcGeom
        (
            const vector x_ij,
            const label m,
            const label n,
            const label o,
            const volIntegralType& integralsj,
            const volIntegralType& integralsi
        );

        //- Check for existing lists in constant folder and read them
        bool readList(const fvMesh& mesh);

        //- Write lists to constant folder
        void writeList(const fvMesh& mesh);

        //- Draw final stencils for postprocessing
        void drawStencils
        (
            const List<labelListList>& stencilsID,
            const List<List<point> >& haloCenters,
            const List<labelListList>& cellToPatchMap,
            const fvMesh& mesh,
            const fileName DirBase
        );
        
        //- Add the coefficients to the matrix A for each row
        void addCoeffs
        (
            scalarRectangularMatrix& A,
            const label cellj,
            const label polOrder,
            const labelList& dim,
            const volIntegralType& volIntegralsIJ
        );


public:

    // Member Functions

        static WENOBase& instance
        (
            const fvMesh& mesh,
            const label polOrder
        )
        {
            static WENOBase instance_
            (
                mesh,
                polOrder
            );

            return instance_;
        }

    // Accessor functions for member variables as const reference

        //- Get necessary lists for runtime operations
        inline const List<labelListList>& stencilsID() const
        {
            return stencilsID_;
        }
        
        inline const List<labelListList>& cellToPatchMap() const 
        {
            return cellToPatchMap_;
        }
        
        inline const labelList& patchToProcMap() const
        {
            return patchToProcMap_;
        }
        
        inline const List<List<point> > haloCenters() const 
        {
            return haloCenters_;
        }
        
        inline const labelListList& ownHalos() const 
        {
            return ownHalos_;
        }
        
        inline const List<List<scalarRectangularMatrix> >& LSmatrix() const
        {
            return LSmatrix_;
        }
        
        inline const List<scalarRectangularMatrix>& B() const 
        {
            return B_;
        }
        
        inline const List<Pair<volIntegralType>>& intBasTrans() const 
        {
            return intBasTrans_;
        }
        
        inline const List<scalar>& refFacAr() const 
        {
            return refFacAr_;
        }
        
        inline const labelListList& dimList() const 
        {
            return dimList_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
