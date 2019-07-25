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
        typedef List< List< List<scalar> > > scalarMatrix;

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
        List<labelListList> cellToPatchMap_;

        //- Information about processor neighbours of patches
        labelList patchToProcMap_;

        //- List of volume integrals over basis functions
        //  Calculated in the reference space
        List<scalarMatrix> volMom_;

        //- Lists of surface integrals of basis functions
        //  Calculated in the reference space
        List<List<scalarMatrix> > intBasTrans_;

        //- List of face areas in the reference space
        List<scalarList> refFacAr_;

        //- Lists of inverse Jacobians for each cell
        List<scalarRectangularMatrix> JInv_;

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
            const scalar m,
            const scalar n,
            const scalar o,
            const scalarMatrix& integralsj,
            const scalarMatrix& integralsi
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

        //- Get necessary lists for runtime operations
        inline List<labelListList>* getPointerStencilID()
        {
            return &stencilsID_;
        };
        inline List<labelListList>* getPointerCellToPatchMap()
        {
            return &cellToPatchMap_;
        };
        inline labelList* getPointerPatchToProcMap()
        {
            return &patchToProcMap_;
        };
        inline List<List<point> >* getPointerHaloCenters()
        {
            return &haloCenters_;
        } ;
        inline labelListList* getPointerOwnHalos()
        {
            return &ownHalos_;
        } ;
        inline List<List<scalarRectangularMatrix> >* getPointerLSmatrix()
        {
            return &LSmatrix_;
        } ;
        inline List<scalarRectangularMatrix>* getPointerB()
        {
            return &B_;
        }
        inline List<List<scalarMatrix> >* getPointerIntBasTrans()
        {
            return &intBasTrans_;
        };
        inline List<scalarList>* getPointerRefFacAr()
        {
            return &refFacAr_;
        };
        inline labelListList* getPointerDimList()
        {
            return &dimList_;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
