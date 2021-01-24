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
    Foam::List3D

Description
    A three dimensional list.
    Substitute for List<List<List<Type>>>

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020>

\*---------------------------------------------------------------------------*/

#ifndef List3D_H
#define List3D_H


#include <type_traits>
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                            Class List3D
\*---------------------------------------------------------------------------*/

template<class Type>
class List3D
{
    static_assert(std::is_trivial<Type>(), "A List3D can only hold classes that are trivial");
    
    private: 
        int sizeX_{0};
        int sizeY_{0};
        int sizeZ_{0};
        
        Type* data_{nullptr};

    public:
    
        // Constructors
            List3D() = default;
            
            //- Constructor with size
            List3D(const int sizeX, const int sizeY, const int sizeZ)
            : sizeX_(sizeX),sizeY_(sizeY),sizeZ_(sizeZ)
            {
                data_ = new Type[sizeX_*sizeY_*sizeZ_];
            }
            
            //- Destructor
            ~List3D() {delete[] data_;}
            
            //- Copy Constructor
            List3D(const List3D& other)
            {
                sizeX_ = other.sizeX_;
                sizeY_ = other.sizeY_;
                sizeZ_ = other.sizeZ_;
                if (data_ !=  nullptr)
                    delete[] data_;
                    
                int fullSize = sizeX_*sizeY_*sizeZ_;
                data_ = new Type[fullSize];
                
                for(int i=0; i < fullSize; ++i)
                {
                    data_[i] = other.data_[i];
                }
            }
            
            //- Copy assignment
            List3D& operator=(const List3D& other)
            {
                sizeX_ = other.sizeX_;
                sizeY_ = other.sizeY_;
                sizeZ_ = other.sizeZ_;
                if (data_ !=  nullptr)
                    delete[] data_;
                    
                int fullSize = sizeX_*sizeY_*sizeZ_;
                data_ = new Type[fullSize];
                
                for(int i=0; i < fullSize; ++i)
                {
                    data_[i] = other.data_[i];
                }
                return *this;
            }
            

            //- Move Constructor
            List3D(List3D&& other)
            {
                sizeX_ = other.sizeX_;
                sizeY_ = other.sizeY_;
                sizeZ_ = other.sizeZ_;
                if (data_ !=  nullptr)
                    delete[] data_;
                    
                data_ = other.data_;
                other.data_ = nullptr;
            }
            
            //- Move assignment 
            List3D& operator=(List3D&& other)
            {
                sizeX_ = other.sizeX_;
                sizeY_ = other.sizeY_;
                sizeZ_ = other.sizeZ_;
                if (data_ !=  nullptr)
                    delete[] data_;
                    
                this->data_ = other.data_;
                other.data_ = nullptr;
                return *this;
            }

        // Access
        inline const Type& operator()(const int i, const int j, const int k) const;
        
        inline Type& operator()(const int i, const int j, const int k);

        
        void resize(const int sizeX, const int sizeY, const int sizeZ)
        {
            const auto oldSize = sizeX_*sizeY_*sizeZ_;
            sizeX_ = sizeX;
            sizeY_ = sizeY;
            sizeZ_ = sizeZ;
            
            const auto newSize = sizeX_*sizeY_*sizeZ_;
            if (newSize > oldSize)
            {
                Type* tmp = new Type[newSize];
                for (int i=0; i < oldSize; ++i)
                {
                    tmp[i] = data_[i];
                }
                delete[] data_;
                data_ = tmp;
            }
        }
        
        //- Initialize content to zero
        void setZero()
        {
            const auto fullSize = sizeX_*sizeY_*sizeZ_;
            for (int i=0; i < fullSize;++i)
            {
                data_[i] = 0;
            }
        }
        
        int sizeX() const {return sizeX_;}
        int sizeY() const {return sizeY_;}
        int sizeZ() const {return sizeZ_;}
};


template<class Type>
inline Type& Foam::List3D<Type>::operator()(const int i, const int j, const int k)
{
    #ifdef FULLDEBUG
        if ((i*(sizeY_*sizeZ_) + j*(sizeZ_) + k) > (sizeX_*sizeY_*sizeZ_))
        {
            std::cerr << "List3D: Access element out of bounds: "
                      << (i*(sizeY_*sizeZ_) + j*(sizeZ_) + k) <<  " of "
                      << (sizeX_*sizeY_*sizeZ_)<<std::endl;
            std::abort();
        }
    #endif
    
    return data_[i*(sizeY_*sizeZ_) + j*(sizeZ_) + k];
}


template<class Type>
inline const Type& Foam::List3D<Type>::operator()(const int i, const int j, const int k) const
{
    #ifdef FULLDEBUG
        if ((i*(sizeY_*sizeZ_) + j*(sizeZ_) + k) > (sizeX_*sizeY_*sizeZ_))
        {
            std::cerr << "List3D: Access element out of bounds: "
                      << (i*(sizeY_*sizeZ_) + j*(sizeZ_) + k) <<  " of "
                      << (sizeX_*sizeY_*sizeZ_)<<std::endl;
            std::abort();
        }
    #endif
    
    return data_[i*(sizeY_*sizeZ_) + j*(sizeZ_) + k];
}




} // End Namespace Foam
#endif
