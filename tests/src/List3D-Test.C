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

Application
    List3D-Test

Description
    Test List3D class 
    
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#include "catch.hpp"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include "List3D.H"
#include "OFstream.H"
#include "IFstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


TEST_CASE("List3D Test","[baseTest]")
{
    // Create random number generator
    std::srand(std::time(nullptr));
    
    
    // Create 3D Matrix with 10x9x8 elements 
    Foam::List3D<int> matrix(10,9,8);
    
    // Create a 3D vector matrix to compare
    std::vector<std::vector<std::vector<int>>> vecMatrix;
    vecMatrix.resize(10);
    for (int i = 0; i< vecMatrix.size(); ++i)
    {
        vecMatrix[i].resize(9);
        for (int j=0; j < vecMatrix[i].size(); j++)
        {
            vecMatrix[i][j].resize(8,0);
        }
    }
    
    // Fill matrices
    for (int i = 0; i < 10; ++i)
    {
        for (int j=0; j < 9; ++j)
        {
            for (int k = 0; k < 8; k++)
            {
                int randomNumber = std::rand();
                matrix(i,j,k) = randomNumber;
                vecMatrix[i][j][k] = randomNumber;
            }
        }
    }
    
    SECTION("I/O Constructor")
    {
        INFO("I/O Constructor");
        Foam::OFstream osIntBasTrans("List3D.dat",Foam::OFstream::streamFormat::BINARY);
        matrix.write(osIntBasTrans);
        
        Foam::IFstream isIntBasTrans("List3D.dat",Foam::IFstream::streamFormat::BINARY);
        Foam::List3D<int> copyMatrix(isIntBasTrans);
        
        REQUIRE(copyMatrix.size() == matrix.size());
        REQUIRE(copyMatrix.sizeX() == matrix.sizeX());
        REQUIRE(copyMatrix.sizeY() == matrix.sizeY());
        REQUIRE(copyMatrix.sizeZ() == matrix.sizeZ());
        
        for (int i = 0; i < 10; ++i)
        {
            for (int j=0; j < 9; ++j)
            {
                for (int k = 0; k < 8; k++)
                {
                    REQUIRE(matrix(i,j,k) == copyMatrix(i,j,k));
                }
            }
        }
    }
    
    
    SECTION("Copy Constructor")
    {
        Foam::List3D<int> copyMatrix(matrix);
        for (int i = 0; i < 10; ++i)
        {
            for (int j=0; j < 9; ++j)
            {
                for (int k = 0; k < 8; k++)
                {
                    REQUIRE(vecMatrix[i][j][k] == copyMatrix(i,j,k));
                }
            }
        }
    }
    
    SECTION("Copy Assignment Operator")
    {
        Foam::List3D<int> copyMatrix(10,9,8);
        copyMatrix = matrix;   
        for (int i = 0; i < 10; ++i)
        {
            for (int j=0; j < 9; ++j)
            {
                for (int k = 0; k < 8; k++)
                {
                    REQUIRE(vecMatrix[i][j][k] == copyMatrix(i,j,k));
                }
            }
        }
    }
    
    
    SECTION("Move Constructor")
    {
        Foam::List3D<int> temp = matrix;
        Foam::List3D<int> copyMatrix(std::move(temp));
        for (int i = 0; i < 10; ++i)
        {
            for (int j=0; j < 9; ++j)
            {
                for (int k = 0; k < 8; k++)
                {
                    REQUIRE(vecMatrix[i][j][k] == copyMatrix(i,j,k));
                }
            }
        }
    }
    
    
    SECTION("Move Assignment")
    {
        Foam::List3D<int> temp = matrix;
        Foam::List3D<int> copyMatrix = std::move(temp);
        for (int i = 0; i < 10; ++i)
        {
            for (int j=0; j < 9; ++j)
            {
                for (int k = 0; k < 8; k++)
                {
                    REQUIRE(vecMatrix[i][j][k] == copyMatrix(i,j,k));
                }
            }
        }
    }
    
    
    SECTION("Access Operator () Constant")
    {
        const Foam::List3D<int> copyMatrix(matrix);
        for (int i = 0; i < 10; ++i)
        {
            for (int j=0; j < 9; ++j)
            {
                for (int k = 0; k < 8; k++)
                {
                    REQUIRE(vecMatrix[i][j][k] == copyMatrix(i,j,k));
                }
            }
        }
    }
    
    
    SECTION("Access Operator ()")
    {
        Foam::List3D<int> copyMatrix(matrix);
        for (int i = 0; i < 10; ++i)
        {
            for (int j=0; j < 9; ++j)
            {
                for (int k = 0; k < 8; k++)
                {
                    REQUIRE(vecMatrix[i][j][k] == copyMatrix(i,j,k));
                }
            }
        }
    }
    
    SECTION("Check += operator")
    {
        Foam::List3D<int> copyMatrix(matrix);
        int temp = copyMatrix(1,1,1);
        
        copyMatrix(1,1,1) += 1;
        REQUIRE(copyMatrix(1,1,1) == (temp+1));
        
    }


    SECTION("I/O Operator ASCII")
    {
        Foam::OFstream osIntBasTrans("List3D-ascii.dat",Foam::OFstream::streamFormat::ASCII);
        matrix.write(osIntBasTrans);
        
        Foam::IFstream isIntBasTrans("List3D-ascii.dat",Foam::IFstream::streamFormat::ASCII);
        Foam::List3D<int> copyMatrix;
        copyMatrix.read(isIntBasTrans);
        
        REQUIRE(copyMatrix.size() == matrix.size());
        REQUIRE(copyMatrix.sizeX() == matrix.sizeX());
        REQUIRE(copyMatrix.sizeY() == matrix.sizeY());
        REQUIRE(copyMatrix.sizeZ() == matrix.sizeZ());
        
        for (int i = 0; i < 10; ++i)
        {
            for (int j=0; j < 9; ++j)
            {
                for (int k = 0; k < 8; k++)
                {
                    REQUIRE(matrix(i,j,k) == copyMatrix(i,j,k));
                }
            }
        }
    }

    SECTION("I/O Operator Binary")
    {
        INFO("Binary Read/Write");
        Foam::OFstream osIntBasTrans("List3D.dat",Foam::OFstream::streamFormat::BINARY);
        matrix.write(osIntBasTrans);
        
        Foam::IFstream isIntBasTrans("List3D.dat",Foam::IFstream::streamFormat::BINARY);
        Foam::List3D<int> copyMatrix;
        copyMatrix.read(isIntBasTrans);
        
        REQUIRE(copyMatrix.size() == matrix.size());
        REQUIRE(copyMatrix.sizeX() == matrix.sizeX());
        REQUIRE(copyMatrix.sizeY() == matrix.sizeY());
        REQUIRE(copyMatrix.sizeZ() == matrix.sizeZ());
        
        for (int i = 0; i < 10; ++i)
        {
            for (int j=0; j < 9; ++j)
            {
                for (int k = 0; k < 8; k++)
                {
                    REQUIRE(matrix(i,j,k) == copyMatrix(i,j,k));
                }
            }
        }
    }
    

}
    

