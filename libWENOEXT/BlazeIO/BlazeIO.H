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
    Foam::BlazeIO

Description
    Functions to use OpenFOAM stream objects with Blaze matrices 

SourceFiles
    BlazeIO.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020

\*---------------------------------------------------------------------------*/



#ifndef BlazeIO_H
#define BlazeIO_H

#include "linear.H"
#include "Ostream.H"
#include "blaze/Math.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Read operators for blaze matrix
Istream& operator >> (Istream&, blaze::DynamicMatrix<double,blaze::rowMajor>&);
Istream& operator >> (Istream&, blaze::DynamicMatrix<double,blaze::columnMajor>&);

// Write operators for blaze matr
Ostream& operator << (Ostream&, const blaze::DynamicMatrix<double,blaze::rowMajor>&);
Ostream& operator << (Ostream&, const blaze::DynamicMatrix<double,blaze::columnMajor>&);

Ostream& operator << (Ostream&, const blaze::StaticMatrix<double,3UL,3UL,blaze::rowMajor>&);
Ostream& operator << (Ostream&, const blaze::StaticMatrix<double,3UL,3UL,blaze::columnMajor>&);

// Write operators for blaze vector
Ostream& operator << (Ostream&, const blaze::DynamicVector<double,blaze::columnVector>&);
Ostream& operator << (Ostream&, const blaze::DynamicVector<double,blaze::rowVector>&);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
