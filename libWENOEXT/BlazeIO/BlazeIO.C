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
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020

\*---------------------------------------------------------------------------*/

#include "BlazeIO.H"

Foam::Istream& Foam::operator >>
(
    Istream& is, 
    blaze::DynamicMatrix<double,blaze::rowMajor>& M
)
{
    // Frist write out the size 
    unsigned int rows;
    unsigned int columns;
    
    is >> rows;
    is >> columns;
    
    M.resize(rows,columns);
    
    if (is.format() == IOstream::ASCII)
    {
        for (unsigned int i=0; i<M.rows(); i++)
        {
            for (unsigned int j=0; j<M.columns(); j++)
            {
                is >> M(i,j);
            }
        }
    }
    else
    {
        unsigned int spacing;
        is >> spacing;
        is.read(reinterpret_cast<char*>(M.data()),rows*spacing*sizeof(double));
    }
    
    return is;
}


Foam::Istream& Foam::operator >>
(
    Istream& is, 
    blaze::DynamicMatrix<double,blaze::columnMajor>& M
)
{
    // Frist write out the size 
    unsigned int rows;
    unsigned int columns;
    
    is >> rows;
    is >> columns;
    
    M.resize(rows,columns);
    
    if (is.format() == IOstream::ASCII)
    {
        for (unsigned int i=0; i<M.rows(); i++)
        {
            for (unsigned int j=0; j<M.columns(); j++)
            {
                is >> M(i,j);
            }
        }
    }
    else
    {
        unsigned int spacing;
        is >> spacing;
        is.read(reinterpret_cast<char*>(M.data()),columns*spacing*sizeof(double));
    }
    
    return is;
}


Foam::Ostream& Foam::operator <<
(
    Ostream& os,
    const blaze::DynamicMatrix<double,blaze::rowMajor>& M
)
{
    // Frist write out the size 
    os << M.rows() << endl;
    os << M.columns() << endl;
    
    if (os.format() == IOstream::ASCII)
    {
        for (unsigned int i=0; i<M.rows(); i++)
        {
            for (unsigned int j=0; j<M.columns(); j++)
            {
                os << M(i,j)<<" ";
            }
            os << endl;
        }
    }
    else
    {
        // Matrix can be padded for alignment. Rows and columns does not give 
        // the spacing
        os << M.spacing()<<endl;
        os.write(reinterpret_cast<const char*>(M.data()),(M.spacing()*M.rows()* sizeof(double)));
        os.flush();
    }
    
    return os;
}


Foam::Ostream& Foam::operator <<
(
    Ostream& os,
    const blaze::DynamicMatrix<double,blaze::columnMajor>& M
)
{
    // Frist write out the size 
    os << M.rows() << endl;
    os << M.columns() << endl;
    
    if (os.format() == IOstream::ASCII)
    {
        for (unsigned int i=0; i<M.rows(); i++)
        {
            for (unsigned int j=0; j<M.columns(); j++)
            {
                os << M(i,j)<<" ";
            }
            os << endl;
        }
    }
    else
    {
        // Matrix can be padded for alignment. Rows and columns does not give 
        // the spacing
        os << M.spacing()<<endl;
        os.write(reinterpret_cast<const char*>(M.data()),(M.spacing()*M.columns()* sizeof(double)));
        os.flush();
    }
    
    return os;
}


Foam::Ostream& Foam::operator <<
(
    Ostream& os,
    const blaze::StaticMatrix<double,3UL,3UL,blaze::rowMajor>& M
)
{
    for (size_t i = 0; i < M.rows(); i++)
    {
        os << M(i,0);
        for (size_t j = 1; j < M.columns(); j++)
            os <<" "<< M(i,j);
        os << "\n";
    }
    return os;
}


Foam::Ostream& Foam::operator <<
(
    Ostream& os,
    const blaze::StaticMatrix<double,3UL,3UL,blaze::columnMajor>& M
)
{
    for (size_t i = 0; i < M.rows(); i++)
    {
        os << M(i,0);
        for (size_t j = 1; j < M.columns(); j++)
            os <<" "<< M(i,j);
        os << "\n";
    }
    return os;
}


Foam::Ostream& Foam::operator << (Ostream& os, const blaze::DynamicVector<double,blaze::columnVector>& V)
{
    for (size_t i = 0; i < V.size() - 1; i++)
    {
        os << V[i] << ",";
    }
    os << V[V.size()-1]<<"\n";
    return os;

}
Foam::Ostream& Foam::operator << (Ostream& os, const blaze::DynamicVector<double,blaze::rowVector>& V)
{
    for (size_t i = 0; i < V.size() - 1; i++)
    {
        os << V[i] << ",";
    }
    os << V[V.size()-1]<<"\n";
    return os;
}
