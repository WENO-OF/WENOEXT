/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    WENOBase Test for parallel runs
    
    Checking if parallel stencil lists are calculated correctly by comparing
    to the single core solution
    
Description
    Test the WENOUpwindFit scheme by using the test case,
    the rotation of a slotted disk, designed by Zalesak
    
    [1] Zalesak, S.T. Fully Multidimensional Flux-Corrected Algorithms
        for Fluids. J. Comput. Phys. 1979, 31, 335–362.
    
Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#include "catch.hpp"
#include "WENOBase.H"
#include "codeRules.H"
#include "fvCFD.H"
#include "KDTree/KDTree.hpp"
#include <vector>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


std::vector<std::vector<double>> convertOpenFOAMPoints(const Foam::List<Foam::point>& coordinates)
{
    std::vector<std::vector<double>> converted;
    converted.resize(coordinates.size());
    
    forAll(coordinates,cellI)
    {
        std::vector<double> point(3,0);
        for (int i = 0; i<3; i++)
        {
            point[i] = coordinates[cellI][i];
        }
        converted[cellI] = point;
    }
    
    return std::move(converted);
}


std::vector<double> convertPoint(const Foam::point& OFPoint)
{
    std::vector<double> point(3,0);
    for (int i = 0; i<3; i++)
    {
        point[i] = OFPoint[i];
    }
    
    return std::move(point);
}

namespace Foam
{
void collectData
(
    const fvMesh& mesh,
    const WENOBase& WENO,
    List<List<point>>& processorData
)
{
    // Distribute data to neighbour processors

    processorData.setSize(WENO.ownHalos().size());

    forAll(processorData, procI)
    {
        processorData[procI].setSize(WENO.ownHalos()[procI].size());

        forAll(processorData[procI], cellI)
        {
            processorData[procI][cellI] =
                mesh.C()[WENO.ownHalos()[procI][cellI]];
        }
    }
    
    
    #ifdef FOAM_PSTREAM_COMMSTYPE_IS_ENUMCLASS 
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    #else
        PstreamBuffers pBufs(Pstream::nonBlocking);
    #endif

    // Distribute data
    forAll(WENO.sendProcList(), procI)
    {
        if (WENO.sendProcList()[procI] != -1)
        {
            UOPstream toBuffer(WENO.sendProcList()[procI], pBufs);
            toBuffer << processorData[procI];
        }
    }

    pBufs.finishedSends();

    // Collect data

    forAll(WENO.receiveProcList(), procI)
    {
        processorData[procI].clear();
        if (WENO.receiveProcList()[procI] != -1)
        {
            UIPstream fromBuffer(WENO.receiveProcList()[procI], pBufs);
            fromBuffer >> processorData[procI];
        }
    }
}
};



TEST_CASE("WENOBase Parallel Test","[3Dparallel]")
{
    // Single core stencil list
    List<List<List<vector>>> stencilCoordinates;
    List<vector> singleCoreMeshCells;
    // Map to store coordinate to cell index
    KDTree* cellIndexTree;
    
    // First create single core results:
    {
        // Replace setRootCase.H for Catch2   
        int argc = 1;
        char executable[] = "main";
        
        char **argv = static_cast<char**>(malloc(argc*sizeof(char*)));
        argv[0] = executable;
        Foam::argList args(argc, argv);
        if (!args.checkRootCase())
        {
            Foam::FatalError.exit();
        }
        #include "createTime.H"        // create the time object
        #include "createMesh.H"        // create the mesh object

        autoPtr<WENOBase> WENOPtr(WENOBase::nonStaticInstance(mesh,3));

        List<labelListList> stencilID = WENOPtr().stencilsID();
        

        stencilCoordinates.resize(stencilID.size());
        singleCoreMeshCells.resize(mesh.nCells());
        
        // convert OpenFOAM point list to std::vector list
        
        cellIndexTree =  new KDTree(convertOpenFOAMPoints(mesh.C()));
        
        
        forAll(stencilID,cellI)
        {
            auto coords = mesh.C()[cellI];
            
            singleCoreMeshCells[cellI] = mesh.C()[cellI];
    
            stencilCoordinates[cellI].resize(stencilID[cellI].size());
            forAll(stencilID[cellI],stencilI)
            {
                stencilCoordinates[cellI][stencilI].resize(stencilID[cellI][stencilI].size());
                forAll(stencilID[cellI][stencilI],j)
                {
                    coords = mesh.C()[stencilID[cellI][stencilI][j]];
                    stencilCoordinates[cellI][stencilI][j] = coords;
                }
            }
        }
    }
    

    // Replace setRootCase.H for Catch2   
    int argc = 2;
    char executable[] = "main";
    char parallel[] = "-parallel";
    
    char **argv = static_cast<char**>(malloc(argc*sizeof(char*)));
    argv[0] = executable;
    argv[1] = parallel;
    Foam::argList args(argc, argv);
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }
    #include "createTime.H"        // create the time object
    #include "createMesh.H"        // create the mesh object

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    // Create the WENO instance
    WENOBase& WENO = WENOBase::instance(mesh,3);

    List<List<point>> processorData;
    collectData(mesh,WENO,processorData);


    // Get parallel stencilList
    auto& parallelStencilList = WENO.stencilsID();
    
    forAll(parallelStencilList,cellI)
    {
        // Find corresponding single core cell index
        const label singleCoreCellI = cellIndexTree->nearest_index(convertPoint(mesh.C()[cellI]));
        
        // Check that it is the same cell
        REQUIRE(singleCoreMeshCells[singleCoreCellI].x() == Approx(mesh.C()[cellI].x()));
        REQUIRE(singleCoreMeshCells[singleCoreCellI].y() == Approx(mesh.C()[cellI].y()));
        REQUIRE(singleCoreMeshCells[singleCoreCellI].z() == Approx(mesh.C()[cellI].z()));
        
        
        //forAll(parallelStencilList[cellI],stencilI)
        label stencilI = 0;
        {
            if (parallelStencilList[cellI][stencilI][0] != int(WENOBase::Cell::deleted))
            {
                const List<label>& cellToProcMapI =
                    WENO.cellToProcMap()[cellI][stencilI];

                List<point> singleCoreCoords = stencilCoordinates[singleCoreCellI][stencilI];
                List<point> parallelCoords(parallelStencilList[cellI][stencilI].size());
                
                // Scope INFO statement
                {
                    INFO
                    (
                        "Processor: "<<Pstream::myProcNo()
                        << " Failed check that single core and parallel stencil list have equal size\n"
                        << "Stencil: "<<stencilI<<" Size singleCore "<<singleCoreCoords.size()
                        << " parallelCore: "<<parallelCoords.size()
                    );
                    
                    REQUIRE(singleCoreCoords.size() == parallelCoords.size());
                }
                
                forAll(parallelCoords,j)
                {
                    if (cellToProcMapI[j] == int(WENOBase::Cell::local))
                        parallelCoords[j] = mesh.C()[parallelStencilList[singleCoreCellI][stencilI][j]];
                    else if(cellToProcMapI[j] != int(WENOBase::Cell::deleted))
                        parallelCoords[j] = processorData[cellToProcMapI[j]][parallelStencilList[cellI][stencilI][j]];
                }
                
                
                auto func = [](const point& p1, const point& p2)
                {
                    if (p1.x() == Approx(p2.x()))
                    {
                        if (p1.y() == Approx(p2.y()))
                        {
                            return p1.z() < p2.z();
                        }
                        else
                        {
                            return p1.y() < p2.y();
                        }
                    }
                    return p1.x() < p2.x();
                };
                
                
                
                // Sort the list
                std::stable_sort
                (
                    singleCoreCoords.begin(),
                    singleCoreCoords.end(),
                    func
                );
                
                std::stable_sort
                (
                    parallelCoords.begin(),
                    parallelCoords.end(),
                    func
                );
                
                forAll(parallelCoords,j)
                {
                    for (int compI=0; compI<3; compI++)
                    {
                        INFO
                        (
                            "Processor: "<<Pstream::myProcNo()
                         << " Checking coordinate of cell: "<<cellI
                         << " stencilI: "<<stencilI << " index: "<<j
                         << " and component "<<compI 
                         << " compared to single core cell " << singleCoreCellI
                         << " Coords: ("<<parallelCoords[j].x()<<", "<<parallelCoords[j].y()<<", "
                         << parallelCoords[j].z()<<")"
                         << " vs ("<<singleCoreCoords[j].x()<<", "
                         << singleCoreCoords[j].y()<<", "
                         << singleCoreCoords[j].z()<<")\n"
                         << "Cell Type "<< cellToProcMapI[j]
                        );
                            
                            
                        CHECKED_ELSE(singleCoreCoords[j][compI] == Approx(parallelCoords[j][compI]))
                        {
                            // Write out stencil to visualize in paraview
                            volScalarField stencilFIELDsingleCore
                            (
                                IOobject
                                (
                                    "stencilListSingleCore",
                                    mesh.time().timeName(),
                                    mesh,
                                    IOobject::NO_READ,
                                    IOobject::AUTO_WRITE
                                ),
                                mesh,
                                dimensionedScalar("0", dimless, 0.0)
                            );
                            volScalarField stencilFIELDParallel
                            (
                                IOobject
                                (
                                    "stencilListParallel",
                                    mesh.time().timeName(),
                                    mesh,
                                    IOobject::NO_READ,
                                    IOobject::AUTO_WRITE
                                ),
                                mesh,
                                dimensionedScalar("0", dimless, 0.0)
                            );
                            
                            // Built tree
                            KDTree meshSearch(convertOpenFOAMPoints(mesh.C()));
                            
                            forAll(parallelStencilList[cellI][stencilI],j)
                            {
                                if (cellToProcMapI[j] == int(WENOBase::Cell::local))
                                    stencilFIELDParallel[parallelStencilList[cellI][stencilI][j]]
                                    = 1.0;

                                stencilFIELDsingleCore
                                [
                                    meshSearch.nearest_index
                                    (
                                        convertPoint(singleCoreCoords[j])
                                    )
                                ] = 2;
                                if (cellToProcMapI[j] == int(WENOBase::Cell::local))
                                    Info << "parallel Cell: ["<<j<<"]  "
                                         <<"(" << parallelCoords[j].x() << ", "
                                         << parallelCoords[j].y() << ", "
                                         << parallelCoords[j].z()<<")  "
                                         <<"(" << singleCoreCoords[j].x() << ", "
                                         << singleCoreCoords[j].y() << ", "
                                         << singleCoreCoords[j].z()<<")  "<<cellToProcMapI[j]<<endl;
                                else
                                    Info << "parallel Cell: ["<<j<<"]  "
                                         <<"(" << parallelCoords[j].x() << ", "
                                         << parallelCoords[j].y() << ", "
                                         << parallelCoords[j].z()<<")  "
                                         <<"(" << singleCoreCoords[j].x() << ", "
                                         << singleCoreCoords[j].y() << ", "
                                         << singleCoreCoords[j].z()<<")  Processor"<<endl;
                                
                            }
                            
                            stencilFIELDParallel[cellI] = 10;
                            stencilFIELDsingleCore[cellI] = 10;
                            stencilFIELDParallel.write();
                            stencilFIELDsingleCore.write(); 
                            
                            REQUIRE(singleCoreCoords[j][compI] == Approx(parallelCoords[j][compI]));
                        }
                    }
                    
                }
                
                
                

                //forAll(parallelStencilList[cellI][stencilI],j)
                //{
                    //auto& coords 
                        //= mesh.C()[parallelStencilList[cellI][stencilI][j]];
                            
                    //auto singleCoreCoords = stencilCoordinates[singleCoreCellI][stencilI][j];
                    
                    //Pout << "j: "<<j<<" ("<< coords.x()<<", "<<coords.y()<<", "<<coords.z() <<")"
                         //<< " vs ("<< singleCoreCoords.x()<<", "<<singleCoreCoords.y()<<", "<<singleCoreCoords.z() <<")"<<endl;
                    
                //}
                //Pout << "---------------------------------------"<<endl;

            }
        }
    }
    
    delete cellIndexTree;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

