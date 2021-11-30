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
    WENOBase I/O Test 
    
Description
    Test the I/O function of WENOBase
    
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#include "catch.hpp"
#include "WENOBase.H"
#include "fvCFD.H"
#include <chrono>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ListType>
void checkList(const ListType& L1, const ListType& L2)
{
    forAll(L1,index)
    {
        checkList(L1[index],L2[index]);
    }
}


template<>
void checkList(const List<label>& L1, const List<label>& L2)
{
    forAll(L1,i)
    {
        REQUIRE(L1[i] == L2[i]);
    }
}

template<>
void checkList(const List<scalar>& L1, const List<scalar>& L2)
{
    forAll(L1,i)
    {
        REQUIRE(L1[i] == L2[i]);
    }
}

template<>
void checkList(const geometryWENO::DynamicMatrix& M1, const geometryWENO::DynamicMatrix& M2)
{
    for (unsigned int i=0; i<M1.rows(); i++)
    {
        for (unsigned int j=0; j<M1.columns(); j++)
        {
            REQUIRE(M1(i,j) == M2(i,j));
        }
    }
}


void checkVolIntegral(const geometryWENO::volIntegralType& L1, const geometryWENO::volIntegralType& L2)
{
    for (int l = 0; l < L1.sizeX(); ++l)
    {
        for (int m = 0; m < L1.sizeY(); ++m)
        {
            for (int n = 0; n < L1.sizeZ(); ++n)
            {
                REQUIRE(Approx(L1(l,m,n)) == L2(l,m,n));
            }
        }
    }
}

template<>
void checkList(const Pair<geometryWENO::volIntegralType>& L1, const Pair<geometryWENO::volIntegralType>& L2)
{
    checkVolIntegral(L1.first(),L2.first());
    checkVolIntegral(L1.second(),L2.second());
}





TEST_CASE("WENOBase IO Test","[2DMesh][singleCore][IOTest]")
{
    // ------------------------------------------------------------------------
    //                          OpenFOAM Start-Up 
    // ------------------------------------------------------------------------
    // Replace setRootCase.H for Catch2   
    int argc = 1;
    char **argv = static_cast<char**>(malloc(sizeof(char*)));
    char executable[] = {'m','a','i','n'};
    argv[0] = executable;
    Foam::argList args(argc, argv,false,false,false);
    #include "createTime.H"        // create the time object
    #include "createMesh.H"        // create the mesh object

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    WENOBase& WENO = WENOBase::instance(mesh,3);
    
    WENO.writeList(mesh);
    
    // Store current Lists
    List<labelListList> stencilID = WENO.stencilsID();
    List<labelListList> cellToProcMap = WENO.cellToProcMap();
    labelList receiveProcList = WENO.receiveProcList();
    labelList sendProcList = WENO.sendProcList();
    labelListList sendHaloCellIDList = WENO.sendHaloCellIDList();
    List<geometryWENO::DynamicMatrix> B = WENO.B();
    List<Pair<geometryWENO::volIntegralType>> intBasTrans = WENO.intBasTrans();
    List<scalar> refFacAr =WENO.refFacAr();
    labelListList dimList = WENO.dimList();
    label degreesOfFreedom = WENO.degreesOfFreedom();
    List<List<geometryWENO::DynamicMatrix>> LSMatrix;
    LSMatrix.resize(WENO.LSmatrix().size());
    forAll(LSMatrix,cellI)
    {
        LSMatrix[cellI].resize(WENO.LSmatrix()[cellI].size());
        forAll(LSMatrix[cellI],stencilI)
        {
            if (WENO.LSmatrix()[cellI][stencilI].valid())
                LSMatrix[cellI][stencilI] = WENO.LSmatrix()[cellI][stencilI]();
        }
    }
    
    
    // Read the data
    WENO.readList(mesh);
    
    // Check that the entries are the same
    
    INFO("Check stencilID ...");
    checkList(stencilID,WENO.stencilsID());
    INFO("Check cellToProcMap ...");
    checkList(cellToProcMap,WENO.cellToProcMap());
    INFO("Check receiveProcList ...");
    checkList(receiveProcList,WENO.receiveProcList());
    INFO("Check sendProcList ...");
    checkList(sendProcList,WENO.sendProcList());
    INFO("Check sendHaloCellIDList ...");
    checkList(sendHaloCellIDList,WENO.sendHaloCellIDList());
    INFO("Check B List ...")
    checkList(B,WENO.B());
    INFO("Check intBasTrans ...")
    checkList(intBasTrans,WENO.intBasTrans());
    INFO("Check refFacAr ...");
    checkList(refFacAr,WENO.refFacAr());
    INFO("Check dimList ...")
    checkList(dimList,WENO.dimList());
    INFO("Check LSMatrix ...")
    forAll(LSMatrix,cellI)
    {
        forAll(LSMatrix[cellI],stencilI)
        {
            if (WENO.LSmatrix()[cellI][stencilI].valid())
                checkList(LSMatrix[cellI][stencilI],WENO.LSmatrix()[cellI][stencilI]());
        }
    }
    REQUIRE(degreesOfFreedom == WENO.degreesOfFreedom());
    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

