


#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"



#include <iostream>
#include <fstream>

void read_phase_space_factor(const char* buffer, double &value)
{

    const char* p{buffer};
    while(*p != '0')
    {
        ++ p;
    }
    std::string value_string;
    while((*p != '\0') && (*p != ' '))
    {
        value_string.push_back(*p);
        ++ p;
    }

    value = std::stod(value_string);

}

int main()
{

    // ifstream objects
    std::ifstream ifs_nEqNull("nEqNull.dat", std::ios::ate);
    std::ifstream ifs_nEqTwo("nEqTwo.dat", std::ios::ate);

    std::ifstream ifs_psiN0("psiN0.txt", std::ios::ate);
    std::ifstream ifs_psiN2("psiN2.txt", std::ios::ate);

    // get size
    std::streampos ifs_nEqNull_size{ifs_nEqNull.tellg()};
    std::streampos ifs_nEqTwo_size{ifs_nEqTwo.tellg()};
    
    std::streampos ifs_psiN0_size{ifs_psiN0.tellg()};
    std::streampos ifs_psiN2_size{ifs_psiN2.tellg()};

    // allocate buffers
    char* buf_nEqNull{new char[ifs_nEqNull_size + 1]};
    char* buf_nEqTwo{new char[ifs_nEqTwo_size + 1]};
    
    char* buf_psiN0{new char[ifs_psiN0_size + 1]};
    char* buf_psiN2{new char[ifs_psiN2_size + 1]};

    // seek
    ifs_nEqNull.seekg(0);
    ifs_nEqTwo.seekg(0);

    ifs_psiN0.seekg(0);
    ifs_psiN2.seekg(0);

    // read
    ifs_nEqNull.read(buf_nEqNull, ifs_nEqNull_size);
    ifs_nEqTwo.read(buf_nEqTwo, ifs_nEqTwo_size);

    ifs_psiN0.read(buf_psiN0, ifs_psiN0_size);
    ifs_psiN2.read(buf_psiN2, ifs_psiN2_size);

    // close streams
    ifs_nEqNull.close();
    ifs_nEqTwo.close();

    ifs_psiN0.close();
    ifs_psiN2.close();

    // allocate data storage
    std::vector<std::vector<double>> data_nEqNull;
    std::vector<std::vector<double>> data_nEqTwo;

    double psiN0;
    double psiN2;

    // extract phase space factors
    //std::cout << buf_psiN0 << std::endl;
    //std::cout << buf_psiN2 << std::endl;
    read_phase_space_factor(buf_psiN0, psiN0);
    read_phase_space_factor(buf_psiN2, psiN2);

    std::cout << psiN0 << std::endl;
    std::cout << psiN2 << std::endl;

    // read data
    read_data(buf_nEqNull, data_nEqNull);
    read_data(buf_nEqTwo, data_nEqTwo);
    
    return 0;

}
