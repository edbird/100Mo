


#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
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

template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<std::vector<T>> v)
{
    typename std::vector<std::vector<T>>::const_iterator it{v.cbegin()};
    for(; it != v.cend(); ++ it)
    {
        os << *it << "\n";
    }
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> v)
{
    typename std::vector<T>::const_iterator it{v.cbegin()};
    for(; it != v.cend(); ++ it)
    {
        os << *it << " ";
    }
    return os;
}

void read_data(const char* buffer, std::vector<std::vector<double>>& data)
{
    
    data.clear();
    data.emplace_back(std::vector<double>());

    //std::cout << buffer << std::endl;

    const char* p{buffer};
    char* p_end;
    while(*p != '\0')
    {
        while(*p == ' ') ++ p;
        double value{strtod(p, &p_end)};
        if(p == p_end)
        {
            data.pop_back();
            break;
        }
        data.back().emplace_back(value);
        p = p_end;
        if(*p == '\n')
        {
            //std::cout << data.at(data.size() - 1) << std::endl;
            data.emplace_back(std::vector<double>());
        }
    }

    //std::cout << "data read" << std::endl;
    //std::cout << "dimension x: " << data.size() << std::endl;
    //for(std::size_t i{0}; i < data.size(); ++ i)
    //{
    //    std::cout << "x=" << i << " dimension y: " << data[i].size() << std::endl;
    //    std::cin.get();
    //}
    //std::cout << "dimension y: " << data.at(0).size() << std::endl;

}

int main()
{

    // ifstream objects
    std::ifstream ifs_nEqNull("nEqNull.dat", std::ios::ate);
    //std::ifstream ifs_nEqNull("test.dat", std::ios::ate);
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
    std::cout << "Finished reading " << "nEqNull.dat" << std::endl;
    read_data(buf_nEqTwo, data_nEqTwo);
    std::cout << "Finished reading " << "nEqTwo.dat" << std::endl;

    //std::cout << buf_nEqNull << std::endl;
    //std::cout << buf_nEqTwo << std::endl;

    std::vector<std::vector<double>> ratio;
    ratio.resize(data_nEqNull.size());
    for(std::size_t i{0}; i < ratio.size(); ++ i)
    {
        ratio[i].resize(data_nEqNull[i].size());
    }

    const double epsilon_31{0.8};
    for(std::size_t i{0}; i < ratio.size(); ++ i)
    {
        for(std::size_t j{0}; j < 2 /*ratio[i].size()*/; ++ j)
        {
            ratio[i][j] = data_nEqNull[i][j];
        }
        ratio[i][2] = (data_nEqNull[i][2] + epsilon_31 * data_nEqTwo[i][2]) / data_nEqNull[i][2];
    }
    std::cout << "Finished constructing ratio" << std::endl;

    //std::cout << ratio << std::endl;

    Int_t dimension_xy{1001};
    TH2D *h_nEqNull = new TH2D("h_nEqNull", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    TH2D *h_nEqTwo = new TH2D("h_nEqTwo", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    TH2D *h_ratio = new TH2D("h_ratio", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    h_ratio->SetStats(0);
    h_nEqNull->SetStats(0);
    h_nEqTwo->SetStats(0);

    for(std::size_t i{0}; i < dimension_xy; ++ i)
    {
        for(std::size_t j{0}; j < dimension_xy; ++ j)
        {
            h_nEqNull->SetBinContent(i, j, data_nEqNull.at(i * dimension_xy + j)[2]);
            h_nEqTwo->SetBinContent(i, j, data_nEqTwo.at(i * dimension_xy + j)[2]);
            //if(i < dimension_xy - j)
            if(i + j < dimension_xy - 1)
            {
                h_ratio->SetBinContent(i, j, ratio.at(i * dimension_xy + j)[2]);
            }
        }
    }
    std::cout << "Finished constructing histogram" << std::endl;

    TCanvas *c_ratio = new TCanvas("c_ratio", "", 4000, 3000);
    h_ratio->Draw("colz");
    c_ratio->SaveAs("c_ratio.png");
    c_ratio->SaveAs("c_ratio.pdf");
    c_ratio->SaveAs("c_ratio.C");
    delete c_ratio;
    
    TCanvas *c_nEqNull = new TCanvas("c_nEqNull", "", 4000, 3000);
    h_nEqNull->Draw("colz");
    c_nEqNull->SaveAs("c_nEqNull.png");
    c_nEqNull->SaveAs("c_nEqNull.pdf");
    c_nEqNull->SaveAs("c_nEqNull.C");
    delete c_nEqNull;

    TCanvas *c_nEqTwo = new TCanvas("c_nEqTwo", "", 4000, 3000);
    h_nEqTwo->Draw("colz");
    c_nEqTwo->SaveAs("c_nEqTwo.png");
    c_nEqTwo->SaveAs("c_nEqTwo.pdf");
    c_nEqTwo->SaveAs("c_nEqTwo.C");
    delete c_nEqTwo;

    return 0;

}
