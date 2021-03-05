/*---------------------------------------------------------------------------*\
                            Write VTK Files

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




