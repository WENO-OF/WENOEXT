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

std::map<Foam::scalar,Foam::scalarRectangularMatrix>::iterator
Foam::matrixDB::similar
(
    const scalarRectangularMatrix&& A
)
{
    // Calculate the unique identifier of A 
    scalar sum = 0;
    for (int i = 0; i < A.m(); i++)
    {
        for (int j = 0; j < A.n(); j++)
        {
            sum += A[i][j]*i*j;
        }
    }
    
    // Get lower bound and upper bound of entry 
    auto it = DB_.find(sum);
    
    if (it == DB_.end())
    {
        auto upper = DB_.upper_bound(sum);
        
        if (upper == DB_.end() && DB_.size() == 0)
        {
            auto pair = DB_.emplace(sum,A);
            return pair.first;
        }
        else if (upper == DB_.end() &&  DB_.size() == 1)
        {
            it = upper--;
        }
        else
        {
            if (upper == DB_.end())
                upper--;
            
            scalar sumUpper = upper->first;
        
            // As the lower bound points to either the exact key or above
            // it has to be decreased by one for the value below it
            auto lower = upper--;
            scalar sumLower = lower->first;
            
            // get the closer value of the key 
            it = mag(sumLower-sum) < mag(sumUpper-sum) ? lower : upper;
        }
    }
    
    
    // Check if it is within the tolerance
    const scalarRectangularMatrix& cmpA = it->second;
    
    if (cmpA.size() == A.size())
    {
        scalar diff = 0;
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
    
    auto pair = DB_.emplace(sum,A);
    
    #ifdef FULLDEBUG
        // Check that it was added 
        if (pair.second == false)
        {
            FatalErrorInFunction 
                << "Matrix was not added to database!" << nl
                << "Matrix to Add: "<<A << nl
                << "------------------------------" << nl
                << "Found matrix: " << DB_[sum] << nl
                << exit(FatalError);
        }
    #endif
    return pair.first;
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
