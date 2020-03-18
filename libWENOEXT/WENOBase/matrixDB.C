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

\*---------------------------------------------------------------------------*/

#include "matrixDB.H"
#include <stdint.h>
#include <inttypes.h>
// * * * * * * * * * * *  ScalarRectangularMatrixPtr * * * * * * * * * * * * //

Foam::matrixDB::scalarRectangularMatrixPtr::scalarRectangularMatrixPtr(matrixDB* db)
: 
    matrixDB_(db),
    itr_(db->DB_.end())
{}


void Foam::matrixDB::scalarRectangularMatrixPtr::add
(
    const scalarRectangularMatrix&& A
)
{
    // search the databank for a similar matrix and return iterator
    itr_ = matrixDB_->similar(std::move(A));
}


const Foam::scalarRectangularMatrix& 
Foam::matrixDB::scalarRectangularMatrixPtr::operator()() const
{
    return itr_->second;
}


bool Foam::matrixDB::scalarRectangularMatrixPtr::valid()
{
    if ((matrixDB_ != nullptr) && (itr_ != matrixDB_->DB_.end()))
        return true;

    return false;
}


// * * * * * * * * * * * * * * * matrixDB  * * * * * * * * * * * * * * * * * //

std::multimap<Foam::scalar,Foam::scalarRectangularMatrix>::const_iterator
Foam::matrixDB::similar
(
    const scalarRectangularMatrix&& A
)
{
    scalar key = 0;
    for (int i = 0; i < A.m(); i++)
    {
        for (int j = 0; j < A.n(); j++)
        {
            key += A[i][j]*i*j;
        }
    }

    if (DB_.size() == 0)
    {
        auto it = DB_.emplace(key,A);
        return it;
    }
    
    // Get lower bound and upper bound of entry 
    auto bound = DB_.equal_range(key);
    
    // Loop over the bound
    for (auto it = bound.first; it != bound.second; ++it)
    {    
        // Check if it is within the tolerance
        const scalarRectangularMatrix& cmpA = it->second;
        
        scalar diff = 0;
        if (cmpA.size() == A.size())
        {
            for (int i = 0; i < A.m(); i++)
            {
                for (int j = 0; j < A.n(); j++)
                {
                    diff += mag(cmpA[i][j] - A[i][j]);
                }
            }
            
            if (diff < tol_)
            {
                counter_++;
                return it;
            }
        }        
    }

    // Use insert with hint 
    auto it = DB_.insert
    (
        bound.first,
        std::pair<scalar,scalarRectangularMatrix>(key,A)
    );
    
    return it;
}


int Foam::matrixDB::hashMatrix
(
    const scalarRectangularMatrix& A
)
{
    /* Returns a representation of the specified floating-point value
     * according to the IEEE 754 floating-point "double
     * format" bit layout, preserving Not-a-Number (NaN) values.
     * 
     * Bit 63 (the bit that is selected by the mask 0x8000000000000000L) 
     * represents the sign of the floating-point number. Bits 62-52 
     * (the bits that are selected by the mask 0x7ff0000000000000L represent 
     * the exponent. Bits 51-0 (the bits that are selected by the mask 
     * 0x000fffffffffffffL) represent the significand (sometimes called the 
     * mantissa) of the floating-point number. 
     *
     * If the argument is positive infinity, the result is
     * 0x7ff0000000000000L.
     * 
     * If the argument is negative infinity, the result is
     * 0xfff0000000000000L.
     * 
     * If the argument is NaN, the result is the long
     * integer representing the actual NaN value.  
     */
    auto doubleToRawBits = [](double x) -> uint64_t
    {
        uint64_t bits;
        memcpy(&bits, &x, sizeof bits);
        return bits;
    };
    
    // Hash function of Java Arrays.hashCode 
    // (http://developer.classpath.org/doc/java/util/Arrays-source.html)
    int key = 1;
    for (int i = 0; i < A.m(); i++)
    {
        for (int j = 0; j < A.n(); j++)
        {
            long long int l = doubleToRawBits(A[i][j]);
            int elt = static_cast<int>(l ^ (l >> 32));
            key = 31*key + elt;
        }
    }
    return key;
}


const Foam::List<Foam::matrixDB::scalarRectangularMatrixPtr>& 
Foam::matrixDB::operator[](const label celli) const
{
    return LSmatrix_[celli];
}


Foam::List<Foam::matrixDB::scalarRectangularMatrixPtr>& 
Foam::matrixDB::operator[](const label celli)
{
    return LSmatrix_[celli];
}


void Foam::matrixDB::resizeSubList(const label cellI, const label size)
{
    LSmatrix_[cellI].resize(size,scalarRectangularMatrixPtr(this));
}


void Foam::matrixDB::info()
{
    scalar numElements = 0;
    forAll(LSmatrix_,celli)
    {
        forAll(LSmatrix_[celli],stencilI)
        {
            if (LSmatrix_[celli][stencilI].valid())
                numElements++;
        }
    }
    
    Pout << "\tMatrix Database Statistics: "<<nl
         << "\t\tTotal Number of matrices: "<< numElements << nl
         << "\t\tNumber matrices stored: "<<DB_.size() <<nl
         << "\t\tCounter: "<<counter_<< endl;
}


void Foam::matrixDB::write(Ostream& os) const
{    
    os << LSmatrix_.size()<<endl;
    forAll(LSmatrix_,cellI)
    {
        os << LSmatrix_[cellI].size()<<endl;
        forAll(LSmatrix_[cellI],stencilI)
        {
            os << LSmatrix_[cellI][stencilI].iterator()->first<<endl;
            os << LSmatrix_[cellI][stencilI].iterator()->second;
        }
    }
}


void Foam::matrixDB::read(Istream& is)
{
    scalarRectangularMatrix matrix;
    scalar key;
    // Read in the LSMatrix list
    label size;
    is >> size;
    LSmatrix_.resize(size);
    
    forAll(LSmatrix_,cellI)
    {
        is >> size;
        LSmatrix_[cellI].resize(size,scalarRectangularMatrixPtr(this));
        
        forAll(LSmatrix_[cellI],stencilI)
        {
            is >> key;
            is >> matrix;
            DB_.emplace(key,matrix);
            LSmatrix_[cellI][stencilI].add(std::move(matrix));
        }
    }
}
