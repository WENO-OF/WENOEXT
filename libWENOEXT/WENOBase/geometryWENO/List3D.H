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
#include "Ostream.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                            Class List3D
\*---------------------------------------------------------------------------*/

// Forward declaration for I/O functions
template<class Type>
class List3D;

template<class Type>
Istream& operator >> (Istream&, List3D<Type>&);

template<class Type>
Ostream& operator << (Ostream&, const List3D<Type>&);



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
            List3D(const int sizeX, const int sizeY, const int sizeZ);
            
            //- Construct with ISstream
            List3D(Istream& is);
            
            //- Destructor
            ~List3D() {delete[] data_;}
            
            //- Copy Constructor
            List3D(const List3D& other);
            
            //- Copy assignment
            List3D& operator=(const List3D& other);
            
            //- Move Constructor
            List3D(List3D&& other);
            
            //- Move assignment 
            List3D& operator=(List3D&& other);

        // Access
            inline const Type& operator()(const int i, const int j, const int k) const;
            
            inline Type& operator()(const int i, const int j, const int k);

            //- Resize data array
            void resize(const int sizeX, const int sizeY, const int sizeZ);
            
            //- Initialize content to zero
            void setZero();
            
            int size() const  {return sizeX_*sizeY_*sizeZ_;}
            int sizeX() const {return sizeX_;}
            int sizeY() const {return sizeY_;}
            int sizeZ() const {return sizeZ_;}
            
        // I/O Functions
        void write(Ostream& os) const;
        
        void read(Istream& is);
        
        template<class Type2>
        friend Istream& operator>>(Istream& is, List3D<Type2>&);
        
        template<class Type2>
        friend Ostream& operator<<(Ostream& is, const List3D<Type2>&);
};
} // End Namespace Foam



// * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * * *


template<class Type>
Foam::List3D<Type>::List3D(const int sizeX, const int sizeY, const int sizeZ)
: sizeX_(sizeX),sizeY_(sizeY),sizeZ_(sizeZ)
{
    data_ = new Type[sizeX_*sizeY_*sizeZ_];
}


template<class Type>
Foam::List3D<Type>::List3D(Istream& is)
{
    this->read(is);
}



template<class Type>
Foam::List3D<Type>::List3D(const List3D& other)
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


template<class Type>
Foam::List3D<Type>& Foam::List3D<Type>::operator=(const List3D& other)
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


template<class Type>
Foam::List3D<Type>::List3D(List3D&& other)
{
    sizeX_ = other.sizeX_;
    sizeY_ = other.sizeY_;
    sizeZ_ = other.sizeZ_;
    if (data_ !=  nullptr)
        delete[] data_;
        
    data_ = other.data_;
    other.data_ = nullptr;
}


template<class Type>
Foam::List3D<Type>& Foam::List3D<Type>::operator=(List3D&& other)
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



// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * 

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


template<class Type>
void Foam::List3D<Type>::resize(const int sizeX, const int sizeY, const int sizeZ)
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


template<class Type>
void Foam::List3D<Type>::setZero()
{
    const auto fullSize = sizeX_*sizeY_*sizeZ_;
    for (int i=0; i < fullSize;++i)
    {
        data_[i] = 0;
    }
}


// * * * * * * * * * * * * * * * * I/O Functions * * * * * * * * * * * * * * *

template<class Type>
void Foam::List3D<Type>::write(Ostream& os) const
{
    if (os.format() == IOstream::ASCII)
    {
        os << sizeX_<<endl;
        os << sizeY_<<endl;
        os << sizeZ_<<endl;
        auto fullSize = sizeX_*sizeY_*sizeZ_;
        for (int i=0; i<fullSize; i++)
        {
            os << data_[i] << endl;
        }
    }
    else
    {
        os << sizeX_<<endl;
        os << sizeY_<<endl;
        os << sizeZ_<<endl;
        auto fullSize = sizeX_*sizeY_*sizeZ_;
        os.write(reinterpret_cast<const char*>(data_),fullSize*sizeof(Type));
        os.flush();
    }
}


template<class Type>
void Foam::List3D<Type>::read(Istream& is)
{
    if (is.format() == IOstream::ASCII)
    {
        std::cout << "ASCII is executed" <<std::endl;
        is >> sizeX_;
        is >> sizeY_;
        is >> sizeZ_;
        auto fullSize = sizeX_*sizeY_*sizeZ_;
        if (data_ !=  nullptr)
            delete[] data_;
        
        data_ = new Type[fullSize];
        
        for (int i=0; i<fullSize; i++)
        {
            is >> data_[i];
        }
    }
    else
    {
        is >> sizeX_;
        is >> sizeY_;
        is >> sizeZ_;
        auto fullSize = sizeX_*sizeY_*sizeZ_;
        if (data_ != nullptr)
            delete[] data_;
        data_ = new Type[fullSize];
        is.read(reinterpret_cast<char*>(data_),fullSize*sizeof(Type));
    }
}


template<class Type>
Foam::Istream& Foam::operator >>(Istream& is, List3D<Type>& List)
{
    List.read(is);
    return is;
}

template<class Type>
Foam::Ostream& Foam::operator <<(Ostream& os,const List3D<Type>& List)
{
    List.write(os);
    return os;
}

#endif
