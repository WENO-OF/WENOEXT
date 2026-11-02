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

#include "EigenIO.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

template<class MatrixType>
void readMatrix(Foam::Istream& is, MatrixType& M)
{
    Foam::label rows;
    Foam::label columns;

    is >> rows;
    is >> columns;

    M.resize(rows,columns);

    if (is.format() == Foam::IOstream::ASCII)
    {
        for (Eigen::Index i=0; i<M.rows(); i++)
        {
            for (Eigen::Index j=0; j<M.cols(); j++)
            {
                is >> M(i,j);
            }
        }
    }
    else
    {
        // Spacing is the size of the inner dimension in memory
        Foam::label spacing;
        is >> spacing;

        const Eigen::Index inner = M.innerSize();
        const Eigen::Index outer = M.outerSize();

        if (spacing == inner)
        {
            is.read
            (
                reinterpret_cast<char*>(M.data()),
                outer*inner*sizeof(double)
            );
        }
        else
        {
            // Data written by the Blaze based version can be padded,
            // read the full buffer and skip the padding
            std::vector<double> buffer(outer*spacing);
            is.read
            (
                reinterpret_cast<char*>(buffer.data()),
                outer*spacing*sizeof(double)
            );

            for (Eigen::Index o=0; o<outer; o++)
            {
                for (Eigen::Index i=0; i<inner; i++)
                {
                    M.data()[o*inner + i] = buffer[o*spacing + i];
                }
            }
        }
    }
}


template<class MatrixType>
void writeMatrix(Foam::Ostream& os, const MatrixType& M)
{
    // Frist write out the size
    os << Foam::label(M.rows()) << Foam::endl;
    os << Foam::label(M.cols()) << Foam::endl;

    if (os.format() == Foam::IOstream::ASCII)
    {
        for (Eigen::Index i=0; i<M.rows(); i++)
        {
            for (Eigen::Index j=0; j<M.cols(); j++)
            {
                os << M(i,j)<<" ";
            }
            os << Foam::endl;
        }
    }
    else
    {
        // Eigen matrices are not padded, the spacing is the size
        // of the inner dimension
        os << Foam::label(M.innerSize()) << Foam::endl;
        os.write
        (
            reinterpret_cast<const char*>(M.data()),
            M.size()*sizeof(double)
        );
        os.flush();
    }
}


template<class MatrixType>
void writeFixedMatrix(Foam::Ostream& os, const MatrixType& M)
{
    for (Eigen::Index i = 0; i < M.rows(); i++)
    {
        os << M(i,0);
        for (Eigen::Index j = 1; j < M.cols(); j++)
            os <<" "<< M(i,j);
        os << "\n";
    }
}


template<class VectorType>
void writeVector(Foam::Ostream& os, const VectorType& V)
{
    for (Eigen::Index i = 0; i < V.size() - 1; i++)
    {
        os << V[i] << ",";
    }
    os << V[V.size()-1]<<"\n";
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator >> (Istream& is, EigenIO::rowMajorMatrix& M)
{
    readMatrix(is, M);
    return is;
}


Foam::Istream& Foam::operator >> (Istream& is, EigenIO::colMajorMatrix& M)
{
    readMatrix(is, M);
    return is;
}


Foam::Ostream& Foam::operator << (Ostream& os, const EigenIO::rowMajorMatrix& M)
{
    writeMatrix(os, M);
    return os;
}


Foam::Ostream& Foam::operator << (Ostream& os, const EigenIO::colMajorMatrix& M)
{
    writeMatrix(os, M);
    return os;
}


Foam::Ostream& Foam::operator << (Ostream& os, const EigenIO::rowMajorMatrix3& M)
{
    writeFixedMatrix(os, M);
    return os;
}


Foam::Ostream& Foam::operator << (Ostream& os, const EigenIO::colMajorMatrix3& M)
{
    writeFixedMatrix(os, M);
    return os;
}


Foam::Ostream& Foam::operator << (Ostream& os, const Eigen::VectorXd& V)
{
    writeVector(os, V);
    return os;
}


Foam::Ostream& Foam::operator << (Ostream& os, const Eigen::RowVectorXd& V)
{
    writeVector(os, V);
    return os;
}


// ************************************************************************* //
