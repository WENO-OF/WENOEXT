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
    matrixDB-Test
    
Description
    Test case for matrixDB
    
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#include "catch.hpp"

#include "fvCFD.H"
#include "matrixDB.H"
#include "OFstream.H"
#include "IFstream.H"
#include "blaze/Math.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TEST_CASE("matrixDB Test Case","[2DMesh][singleCore]")
{
    // Replace setRootCase.H for Catch2   
    int argc = 1;
    char **argv = static_cast<char**>(malloc(sizeof(char*)));
    char executable[] = {'m','a','i','n'};
    argv[0] = executable;
    Foam::argList args(argc, argv,false,false,false);
        
    // create the mesh from case file
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    
    
    
    // -------------------------- Helper Functions -----------------------------
    
    // Function to create matrix
    auto createMatrix = [](const label cellI, const label stencilI) -> scalarRectangularMatrix
    {
        scalarRectangularMatrix A
            (
                5,
                10,
                scalar((cellI*stencilI) % 1000)
            );
        return A;
    };
    
    // function to check the sum 
    auto compareMatrix = [](const scalarRectangularMatrix& A, const blaze::DynamicMatrix<double>& B) -> void
    {
        for (int i = 0; i < A.m(); i++)
        {
            for (int j = 0; j < A.n(); j++)
            {
                REQUIRE(A(i,j) == B(i,j));
            }
        }
    };
    
    
    // ------------------------- Start of Testing ------------------------------
    
    // Create matrixDB object
    matrixDB matrixDataBank;
    
    // Create a list of matrices
    List<List<scalarRectangularMatrix> > LSmatrix;
    LSmatrix.resize(1000);
    
    // Set size of matrix data bank and check size
    matrixDataBank.resize(LSmatrix.size());
    REQUIRE(matrixDataBank.size() == LSmatrix.size());
    
    
    forAll(LSmatrix, cellI)
    {
        LSmatrix[cellI].resize(10);
        matrixDataBank.resizeSubList(cellI,LSmatrix[cellI].size());
        
        // Check the sub list size
        REQUIRE(matrixDataBank[cellI].size() == LSmatrix[cellI].size());
        
        forAll(LSmatrix[cellI],stencilI)
        {
            LSmatrix[cellI][stencilI] = createMatrix(cellI,stencilI);
            matrixDataBank[cellI][stencilI].add(createMatrix(cellI,stencilI));
        }
    }

    
    // Display Info 
    matrixDataBank.info();
    
    // Check that all matrix values are correct
    
    forAll(LSmatrix,cellI)
    {
        forAll(LSmatrix[cellI],stencilI)
        {
            compareMatrix(LSmatrix[cellI][stencilI],matrixDataBank[cellI][stencilI]());
        }
    }


    // ------------------------- Check IO Functions ----------------------------

    fileName path = mesh.time().path()/"constant/matrixDataBankTest";
    
    OFstream osLS(path,OFstream::streamFormat::BINARY);
    matrixDataBank.write(osLS);
    
    matrixDB newMatrixDB;
    IFstream isLS(path,IFstream::streamFormat::BINARY);
    newMatrixDB.read(isLS);
    
    // Check that all values are read in correctly:
    forAll(LSmatrix,cellI)
    {
        forAll(LSmatrix[cellI],stencilI)
        {
            compareMatrix(LSmatrix[cellI][stencilI],newMatrixDB[cellI][stencilI]());
        }
    }
}
