#include <fstream>
#include <vector>
#include <string>

#ifndef writeToFile_H
#define writeToFile_H



template<class T>
void writeToFile(const std::vector<T> data)
{
    std::ofstream file("PLOT/results.dat",std::ofstream::app);
    
    const std::string SEP = "\t";
    
    for (unsigned int i = 0; i < (data.size() -1);i++)
    {
        file << data.at(i) << SEP;
    }
    file << (*--data.end())<< "\n";
    file.close();
};

#endif
