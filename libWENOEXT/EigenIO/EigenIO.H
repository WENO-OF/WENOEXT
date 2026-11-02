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
    Foam::EigenIO

Description
    Functions to use OpenFOAM stream objects with Eigen matrices

    Binary format of dynamic matrices: rows, columns, spacing, raw data.
    The spacing is the size of the inner (storage order) dimension. Files
    written with the former Blaze based version can contain padded data
    (spacing larger than the inner dimension) and are still readable.

SourceFiles
    EigenIO.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020

\*---------------------------------------------------------------------------*/



#ifndef EigenIO_H
#define EigenIO_H

// Include Eigen before OpenFOAM headers to avoid macro clashes
#include <Eigen/Dense>
#include "linear.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace EigenIO
{
    //- Dynamic row major matrix
    using rowMajorMatrix =
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;

    //- Dynamic column major matrix
    using colMajorMatrix =
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>;

    //- 3x3 row major matrix
    using rowMajorMatrix3 = Eigen::Matrix<double,3,3,Eigen::RowMajor>;

    //- 3x3 column major matrix
    using colMajorMatrix3 = Eigen::Matrix<double,3,3,Eigen::ColMajor>;
}

// Read operators for Eigen matrix
Istream& operator >> (Istream&, EigenIO::rowMajorMatrix&);
Istream& operator >> (Istream&, EigenIO::colMajorMatrix&);

// Write operators for Eigen matrix
Ostream& operator << (Ostream&, const EigenIO::rowMajorMatrix&);
Ostream& operator << (Ostream&, const EigenIO::colMajorMatrix&);

Ostream& operator << (Ostream&, const EigenIO::rowMajorMatrix3&);
Ostream& operator << (Ostream&, const EigenIO::colMajorMatrix3&);

// Write operators for Eigen vector
Ostream& operator << (Ostream&, const Eigen::VectorXd&);
Ostream& operator << (Ostream&, const Eigen::RowVectorXd&);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
