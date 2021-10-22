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
#include "Ostream.H"
#include "blaze/Math.h"
#include <map>
#include "geometryWENO.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class matrixDB;

Istream& operator >> (Istream&, matrixDB&);

Ostream& operator << (Ostream&, const matrixDB&);

/*---------------------------------------------------------------------------*\
                            Class matrixDB Declaration
\*---------------------------------------------------------------------------*/

class matrixDB
{
    using DynamicMatrix = blaze::DynamicMatrix<double>;
    using keyType = int32_t;
    
    using iterType = std::multimap<keyType, DynamicMatrix>::iterator;
    using constIterType = std::multimap<keyType, DynamicMatrix>::const_iterator;
    
    class MatrixPtr
    {
        private:
            //- Reference to parent object
            matrixDB* matrixDB_ = nullptr;
            
            //- Iterator to the scalar matrix
            //  An iterator is storead as the standard states:
            //  All Associative Containers: The insert and emplace members shall
            //  not affect the validity of iterators and references
            //  to the container [26.2.6/9]
            iterType itr_;
            
        
        public: 
        
        // Constructor
        
            //- Default construtor
            MatrixPtr() = default;
        
            //- Constructor
            MatrixPtr(matrixDB*);
            
        // Public Member
        
            //- set the iterator
            void set(iterType it) {itr_ = it;}
        
            //- add a new element
            void add(const scalarRectangularMatrix&& A);
            
            //- Dereference the pointer
            //  Throw an execption if called for a nullptr
            const DynamicMatrix& operator()() const;
            
            //- Return const iterator
            const constIterType iterator() const {return itr_;}
            
            //- Return non-const iterator
            iterType iterator() {return itr_;}
            
            //- Check if the container is valid
            bool valid() const;
    };
    
    
    private:
    
        //- Tolerance to accept similar matrices 
        //  Default value 1E-9
        const scalar epsilon_;

        //- Number of key neighbours to check
        const label checkRange_ = 5;
        
        //- counter to store the number of saved matrices through pointer
        int counter_ = 0;
    
        //- Lists of pseudoinverses for each stencil of each cell
        //  Stored as a pointer to the underlying data structure
        List<List<MatrixPtr> > LSmatrix_;

        //- Map of stored matrices as the underlying data bank
        std::multimap<keyType, DynamicMatrix> DB_;


    //- Private member functions
        
        //- Check if a matrix already exist in databank and otherwise add matrix
        //  returns iterator to this matrix
        iterType similar
        (
            const scalarRectangularMatrix&& A
        );
        
        //- Create a hash value of a matrix. 
        keyType hashMatrix(const scalarRectangularMatrix& A);
        
        //- Convert OpenFOAM Matrix to Blaze dynamic matrix
        DynamicMatrix convertToBlaze(const scalarRectangularMatrix& A);

public:

    // Constructors 
    
        //- Default Constructor
        matrixDB(const scalar epsilon = 1E-9) : epsilon_(epsilon) {}
        
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
        inline const List<MatrixPtr>& operator[](const label celli) const 
        {return LSmatrix_[celli];}
        
        inline List<MatrixPtr>& operator[](const label celli) 
        {return LSmatrix_[celli];}
        
    // Access
        
        // get size of LSmatrix list
        label size() const {return LSmatrix_.size();}
        
        //- Print information to screen 
        void info();
        
    // matrixDB IO
    
        void write(Ostream& os) const;
        
        void read(Istream& is);
        
        friend Istream& operator>>(Istream& is, matrixDB&);
        
        friend Ostream& operator<<(Ostream& is, const matrixDB&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
