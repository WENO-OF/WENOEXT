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

Author: Jan Wilhelm Gärtner
\*---------------------------------------------------------------------------*/

#include "vtkWriter.H"

void checkFileName(std::string& fileName)
{
    std::stringstream ss(fileName);
    std::string token;
    while (std::getline(ss,token,'.'))
    
    if (token != "vtk")
        fileName  = fileName + ".vtk";
}


void writeVTKPoints(std::string fileName,const List<point>& points)
{
    checkFileName(fileName);
    
    
    std::ofstream file(fileName);
    file << "# vtk DataFile Version 2.0\n"
         << fileName << "\n"
         << "ASCII\n"
         << "DATASET POLYDATA\n"
         << "POINTS "<<points.size() << " double\n"; 
    
    forAll(points,i)
    {
        file << points[i].x() << " " << points[i].y() << " " << points[i].z() <<"\n";
    }
    file.flush();
    file.close();
}




