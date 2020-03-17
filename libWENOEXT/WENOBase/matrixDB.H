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
    Foam::matrixDB

Description
    Databank for the inverse matrices
    Used to decrease storage demands. Similar matrices are linked such that 
    the information is only stored once 

SourceFiles
    matrixDB.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020

\*---------------------------------------------------------------------------*/



#ifndef matrixDB_H
#define matrixDB_H

#include "linear.H"
#include <map>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class matrixDB Declaration
\*---------------------------------------------------------------------------*/

class matrixDB
{
    class scalarRectangularMatrixPtr
    {
        private:
            //- Const iterator to the scalar matrix
            //  An iterator is storead as the standard states:
            //  All Associative Containers: The insert and emplace members shall
            //  not affect the validity of iterators and references
            //  to the container [26.2.6/9]
            std::map<scalar, scalarRectangularMatrix>::const_iterator itr_;
            
            //- Reference to parent object
            matrixDB* matrixDB_ = nullptr;
        
        public: 
        
        // Constructor
        
            //- Default construtor
            scalarRectangularMatrixPtr() = default;
        
            //- Constructor
            scalarRectangularMatrixPtr(matrixDB*);
            
        // Public Member
        
            //- add a new element
            void add(const scalarRectangularMatrix&& A);
            
            //- Dereference the pointer
            //  Throw an execption if called for a nullptr
            const scalarRectangularMatrix& operator()();
    };
    
    
    private:
    
        //- Tolerance to accept similar matrices 
        //  Default value 1E-9
        const scalar tol_;
        
        //- counter to store the number of saved matrices through pointer
        int counter_ = 0;
    
        //- Lists of pseudoinverses for each stencil of each cell
        //  Stored as a pointer to the underlying data structure
        List<List<scalarRectangularMatrixPtr> > LSmatrix_;

        //- Map of stored matrices as the underlying data bank
        //  The key is the sum of the elements
        std::map<scalar, scalarRectangularMatrix> DB_;


    //- Private member functions
        
        //- Check if a matrix already exist in databank and otherwise add matrix
        //  returns iterator to this matrix
        std::map<scalar,scalarRectangularMatrix>::iterator similar
        (
            const scalarRectangularMatrix&& A
        );

public:

    // Constructors 
    
        //- Default Constructor
        matrixDB(const scalar tol = 1E-9) : tol_(tol) {}
        
        // delete copy constructors
        matrixDB(const matrixDB&) = delete;
    
    // Public Member Functions
    
        //- Resize LSmatrix list
        void resize(const label size)
        {
            LSmatrix_.resize(size);
        }
    
        //- Set size of stencil sub list 
        void resizeSubList(const label cellI, const label size);
        
        //- Access an element
        const List<scalarRectangularMatrixPtr>& operator[](const label celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
