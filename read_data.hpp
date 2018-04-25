#ifndef READDATA_HPP
#define READDATA_HPP

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

void read_data_helper(const std::string& filename_nEqNull, const std::string& filename_nEqTwo,
                      const std::string& filename_psiN0, const std::string& filename_psiN2,
                      std::vector<std::vector<double>> &data_nEqNull_ret, std::vector<std::vector<double>> &data_nEqTwo_ret,
                      double& psiN0_ret, double& psiN2_ret)
{

    ////////////////////////////////////////////////////////////////////////////
    // READ DATA IN FROM FILES
    ////////////////////////////////////////////////////////////////////////////

    // ifstream objects
    //std::ifstream ifs_nEqNull("nEqNull.dat", std::ios::ate);
    std::ifstream ifs_nEqNull(filename_nEqNull.c_str(), std::ios::ate);
    //std::ifstream ifs_nEqNull("test.dat", std::ios::ate);
    //std::ifstream ifs_nEqTwo("nEqTwo.dat", std::ios::ate);
    std::ifstream ifs_nEqTwo(filename_nEqTwo.c_str(), std::ios::ate);

    //std::ifstream ifs_psiN0("psiN0.txt", std::ios::ate);
    std::ifstream ifs_psiN0(filename_psiN0.c_str(), std::ios::ate);
    //std::ifstream ifs_psiN2("psiN2.txt", std::ios::ate);
    std::ifstream ifs_psiN2(filename_psiN2.c_str(), std::ios::ate);

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
    
    // return data
    data_nEqNull_ret = data_nEqNull;
    data_nEqTwo_ret = data_nEqTwo;
    psiN0_ret = psiN0;
    psiN2_ret = psiN2;
    

}

void convert_data_to_histogram_format(const std::vector<std::vector<double>> &data_nEqNull,
                                      const std::vector<std::vector<double>> &data_nEqTwo,
                                      TH2D * h_nEqNull_return,
                                      TH2D * h_nEqTwo_return)
{

    ////////////////////////////////////////////////////////////////////////////
    // CONVERT INPUT DATA TO HISTOGRAM FORMAT
    // Note: Added 2018-04-23 (After INTERMEDIATE DATA below)
    // Note: These histograms do NOT have the phase space variable included
    ////////////////////////////////////////////////////////////////////////////

    const Int_t dimension_xy{1001};
    // don't plot raw data
    TH2D *h_nEqNull = new TH2D("h_nEqNull", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    TH2D *h_nEqTwo = new TH2D("h_nEqTwo", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_ratio = new TH2D("h_ratio", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    h_nEqNull->SetStats(0);
    h_nEqTwo->SetStats(0);
    //h_ratio->SetStats(0);

    for(std::size_t i{0}; i < dimension_xy; ++ i)
    {
        for(std::size_t j{0}; j < dimension_xy; ++ j)
        {
            h_nEqNull->SetBinContent(i, j, data_nEqNull.at(i * dimension_xy + j)[2]);
            h_nEqTwo->SetBinContent(i, j, data_nEqTwo.at(i * dimension_xy + j)[2]);
            //if(i < dimension_xy - j)
            //if(i + j < dimension_xy - 1)
            //{
                // TODO: move above lines to inside this if
                //h_ratio->SetBinContent(i, j, ratio.at(i * dimension_xy + j)[2]);
            //}
        }
    }
    std::cout << "Finished constructing input data histograms" << std::endl;
    
    h_nEqNull_return = h_nEqNull;
    h_nEqTwo_return = h_nEqTwo;

}

void create_data_with_phase_space_factor(std::vector<std::vector<double>> &data_ret, const double epsilon,
                                         std::vector<std::vector<double>> &data_nEqNull,
                                         std::vector<std::vector<double>> &data_nEqTwo,
                                         const double psiN0, const double psiN2)
{
    
    std::vector<std::vector<double>> data;
    data.resize(data_nEqNull.size());
    for(std::size_t i{0}; i < data.size(); ++ i)
    {
        data[i].resize(data_nEqNull[i].size());
    }

    // epsilon_31 = 
    const double epsilon_31{epsilon};
    for(std::size_t i{0}; i < data.size(); ++ i)
    {
        for(std::size_t j{0}; j < 2; ++ j)
        {
            data[i][j] = data_nEqNull[i][j]; // NOTE: this is not used
        }
        double phase_space_factor{1.0 / (psiN0 + epsilon_31 * psiN2)};
        data[i][2] = phase_space_factor * (data_nEqNull[i][2] + epsilon_31 * data_nEqTwo[i][2]);
    }
    
    data_ret = data;
    
}

#endif
