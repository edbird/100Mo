


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

// electron energies: T1, T2 (in units of Q value)
// theory parameter: epsilon
// data: G0, G2 (units of Q value)
//Double_t ReWeight(const Double_t T1, const Double_t T2, const Double_t epsilon,
//                  const std::vector<std::vector<double>>& data_nEqNull,
//                  const std::vector<std::vector<double>>& data_nEqTwo)
Double_t ReWeight(const Double_t T1, const Double_t T2, const Double_t epsilon,
                  const TH2D* const h_nEqNull,
                  const TH2D* const h_nEqTwo,
                  const Double_t psiN0, const Double_t psiN2)
{

    // TODO: NO INTERPOLATION DONE YET

    /*
    std::size_t l{data_nEqNull.size()};
    Int_t sub_length{1001};

    // "primary" index
    Int_t index1{(Int_t)(Double_t)(sub_length) * T1}; index1 *= sub_length;
    Int_t index3{(Int_t)(Double_t)(sub_length) * T2};
    //Int_t index1{(Int_t)((Double_t)(sub_length * sub_length) * T1 + (Double_t)(sub_length) * T2)};
    Int_t index2{index1 + 1};
    //Int_t index3{index1 + sub_length};
    //Int_t index4{index1 + sub_length + 1};
    Int_t index4{index3 + 1};

    Double_t energy1{data_nEqNull.at(index1).at(0)};
    Double_t energy2{data_nEqNull.at(index2).at(1)};
    Double_t energy3{data_nEqNull.at(index3).at(0)};
    Double_t energy4{data_nEqNull.at(index4).at(1)};

    Double_t probability1{data_nEqNull.at(index1).at(2)};
    Double_t probability2{data_nEqNull.at(index2).at(2)};
    Double_t probability3{data_nEqNull.at(index3).at(2)};
    Double_t probability4{data_nEqNull.at(index4).at(2)};
    */

    // interpolate in 2D

    // convert from input energy T1 to "higher and lower energies" as used in
    // histogram data_nEqNull / data_nEqTwo
    
    std::cout << "T1=" << T1 << std::endl;
    std::cout << "T2=" << T2 << std::endl;
    Int_t bin_x{h_nEqNull->GetXaxis()->FindBin(T1)};
    Int_t bin_y{h_nEqNull->GetYaxis()->FindBin(T2)};
    std::cout << "bin_x=" << bin_x << std::endl;
    std::cout << "bin_y=" << bin_y << std::endl;
    Double_t bin_x_low{h_nEqNull->GetXaxis()->GetBinLowEdge(bin_x)};
    Double_t bin_x_high{h_nEqNull->GetXaxis()->GetBinLowEdge(bin_x) + h_nEqNull->GetXaxis()->GetBinWidth(bin_x)};
    Double_t bin_y_low{h_nEqNull->GetXaxis()->GetBinLowEdge(bin_y)};
    Double_t bin_y_high{h_nEqNull->GetXaxis()->GetBinLowEdge(bin_y) + h_nEqNull->GetXaxis()->GetBinWidth(bin_y)};
    std::cout << bin_x_low << " < " << T1 << " < " << bin_x_high << std::endl;
    std::cout << bin_y_low << " < " << T2 << " < " << bin_y_high << std::endl;

    // get the weight for this T1, T2
    // the input data is for epsilon = 0.0
    Double_t phase_1{1.0 / psiN0};
    Double_t weight_1{phase_1 * h_nEqNull->GetBinContent(bin_x, bin_y)};
    std::cout << "weight_1=" << weight_1 << std::endl;

    // get the weight for this T1, T2
    // the input data is for epsilon = some arbitary value
    Double_t phase_2{1.0 / (psiN0 + epsilon * psiN2)};
    Double_t weight_2{phase_2 * (h_nEqNull->GetBinContent(bin_x, bin_y) + epsilon * h_nEqTwo->GetBinContent(bin_x, bin_y))};
    std::cout << "weight_2=" << weight_2 << std::endl;

    std::cout << "reweight factor: " << weight_2 / weight_1 << std::endl;
    return weight_2 / weight_1;

    //std::cout << "T1=" << T1 << " T2=" << T2 << " energy1=" << energy1 << " energy2=" << energy2 << " energy3=" << energy3 << " energy4=" << energy4 << std::endl;
    //std::cout << "T1=" << T1 << " probability1=" << probability1 << " probability2=" << probability2 << " probability3=" << probability3 << " probability4=" << probability4 << std::endl;

}

int main()
{

    ////////////////////////////////////////////////////////////////////////////
    // READ DATA IN FROM FILES
    ////////////////////////////////////////////////////////////////////////////

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

    ////////////////////////////////////////////////////////////////////////////
    // TODO: apply phase space factor to data
    ////////////////////////////////////////////////////////////////////////////

    // ...

    ////////////////////////////////////////////////////////////////////////////
    // CONVERT INPUT DATA TO HISTOGRAM FORMAT
    // Note: Added 2018-04-23 (After INTERMEDIATE DATA below)
    // Note: These histograms do NOT have the phase space variable included
    ////////////////////////////////////////////////////////////////////////////

    Int_t dimension_xy{1001};
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

    ////////////////////////////////////////////////////////////////////////////
    // CREATE INTERMEDIATE DATA
    // Format: Nx3 array, each array a different value of epsilon
    ////////////////////////////////////////////////////////////////////////////

    // create data array for "complete data"
    // (with phase space factors)
    
    std::vector<std::vector<double>> data_0;
    data_0.resize(data_nEqNull.size());
    for(std::size_t i{0}; i < data_0.size(); ++ i)
    {
        data_0[i].resize(data_nEqNull[i].size());
    }

    std::vector<std::vector<double>> data_1;
    data_1.resize(data_nEqNull.size());
    for(std::size_t i{0}; i < data_1.size(); ++ i)
    {
        data_1[i].resize(data_nEqNull[i].size());
    }

    std::vector<std::vector<double>> data_2;
    data_2.resize(data_nEqNull.size());
    for(std::size_t i{0}; i < data_2.size(); ++ i)
    {
        data_2[i].resize(data_nEqNull[i].size());
    }

    // epsilon_31 = 0.0
    const double epsilon_31_0{0.0};
    for(std::size_t i{0}; i < data_0.size(); ++ i)
    {
        for(std::size_t j{0}; j < 2; ++ j)
        {
            data_0[i][j] = data_nEqNull[i][j]; // NOTE: this is not used
        }
        double phase_space_factor{1.0 / (psiN0 + epsilon_31_0 * psiN2)};
        data_0[i][2] = phase_space_factor * (data_nEqNull[i][2] + epsilon_31_0 * data_nEqTwo[i][2]);
    }
    
    // epsilon_31 = 0.4
    const double epsilon_31_1{0.4};
    for(std::size_t i{0}; i < data_1.size(); ++ i)
    {
        for(std::size_t j{0}; j < 2; ++ j)
        {
            data_1[i][j] = data_nEqNull[i][j]; // NOTE: this is not used
        }
        double phase_space_factor{1.0 / (psiN0 + epsilon_31_1 * psiN2)};
        data_1[i][2] = phase_space_factor * (data_nEqNull[i][2] + epsilon_31_1 * data_nEqTwo[i][2]);
    }

    // epsilon_31 = 0.8
    const double epsilon_31_2{0.8};
    for(std::size_t i{0}; i < data_2.size(); ++ i)
    {
        for(std::size_t j{0}; j < 2; ++ j)
        {
            data_2[i][j] = data_nEqNull[i][j]; // NOTE: this is not used
        }
        double phase_space_factor{1.0 / (psiN0 + epsilon_31_2 * psiN2)};
        data_2[i][2] = phase_space_factor * (data_nEqNull[i][2] + epsilon_31_2 * data_nEqTwo[i][2]);
    }


    // create data array for ratio
    /*
    std::vector<std::vector<double>> ratio;
    ratio.resize(data_nEqNull.size());
    for(std::size_t i{0}; i < ratio.size(); ++ i)
    {
        ratio[i].resize(data_nEqNull[i].size());
    }

    // epsilon_31 = 0.8
    const double epsilon_31_2{0.8};
    for(std::size_t i{0}; i < ratio.size(); ++ i)
    {
        for(std::size_t j{0}; j < 2 /*ratio[i].size()*//*; ++ j)
        {
            ratio[i][j] = data_nEqNull[i][j];
        }
        ratio[i][2] = (psiN0 / (psiN0 + epsilon_31 * psiN2)) * ((data_nEqNull[i][2] + epsilon_31 * data_nEqTwo[i][2]) / data_nEqNull[i][2]);
    }
    std::cout << "Finished constructing ratio" << std::endl;
    */

    std::cout << "Finished constructing intermediate data" << std::endl;

    //std::cout << ratio << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    // CREATE OUTPUT HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////

    //Int_t dimension_xy{1001};
    // don't plot raw data
    //TH2D *h_nEqNull = new TH2D("h_nEqNull", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_nEqTwo = new TH2D("h_nEqTwo", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    TH2D *h_data_0 = new TH2D("h_data_0", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    TH2D *h_data_1 = new TH2D("h_data_1", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    TH2D *h_data_2 = new TH2D("h_data_2", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_ratio = new TH2D("h_ratio", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //h_nEqNull->SetStats(0);
    //h_nEqTwo->SetStats(0);
    h_data_0->SetStats(0);
    h_data_1->SetStats(0);
    h_data_2->SetStats(0);
    //h_ratio->SetStats(0);

    for(std::size_t i{0}; i < dimension_xy; ++ i)
    {
        for(std::size_t j{0}; j < dimension_xy; ++ j)
        {
            //h_nEqNull->SetBinContent(i, j, data_nEqNull.at(i * dimension_xy + j)[2]);
            //h_nEqTwo->SetBinContent(i, j, data_nEqTwo.at(i * dimension_xy + j)[2]);
            h_data_0->SetBinContent(i, j, data_0.at(i * dimension_xy + j)[2]);
            h_data_1->SetBinContent(i, j, data_1.at(i * dimension_xy + j)[2]);
            h_data_2->SetBinContent(i, j, data_2.at(i * dimension_xy + j)[2]);
            //if(i < dimension_xy - j)
            if(i + j < dimension_xy - 1)
            {
                // TODO: move above lines to inside this if
                //h_ratio->SetBinContent(i, j, ratio.at(i * dimension_xy + j)[2]);
            }
        }
    }
    std::cout << "Finished constructing histograms" << std::endl;

    /*
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
    */
    
    #define PRINT_1D_CANVAS 0
    #if PRINT_1D_CANVAS

        TCanvas *c_data_0 = new TCanvas("c_data_0", "", 4000, 3000);
        h_data_0->Draw("colz");
        c_data_0->SaveAs("c_data_0.png");
        c_data_0->SaveAs("c_data_0.pdf");
        c_data_0->SaveAs("c_data_0.C");
        delete c_data_0;

        TCanvas *c_data_1 = new TCanvas("c_data_1", "", 4000, 3000);
        h_data_1->Draw("colz");
        c_data_1->SaveAs("c_data_1.png");
        c_data_1->SaveAs("c_data_1.pdf");
        c_data_1->SaveAs("c_data_1.C");
        delete c_data_1;

        TCanvas *c_data_2 = new TCanvas("c_data_2", "", 4000, 3000);
        h_data_2->Draw("colz");
        c_data_2->SaveAs("c_data_2.png");
        c_data_2->SaveAs("c_data_2.pdf");
        c_data_2->SaveAs("c_data_2.C");
        delete c_data_2;

    #endif

    /*
    TCanvas *c_ratio = new TCanvas("c_ratio", "", 4000, 3000);
    h_ratio->Draw("colz");
    c_ratio->SaveAs("c_ratio.png");
    c_ratio->SaveAs("c_ratio.pdf");
    c_ratio->SaveAs("c_ratio.C");
    delete c_ratio;
    */

    // epsilon_31 = 0
    TH1D *h_single_electron_0 = h_data_0->ProjectionX("h_single_electron_0");
    h_single_electron_0->Scale(1.0 / (Double_t)dimension_xy);
    /*for(Int_t i{1}; i <= h_single_electron_0->GetNbinsX(); ++ i)
    {
        Double_t F{1.0 / (Double_t)dimension_xy};
        h_single_electron_0->SetBinContent(i, F * h_single_electron_0->GetBinContent(i));
    }*/
    h_single_electron_0->SetStats(0);
    h_single_electron_0->SetLineColor(2);
    h_single_electron_0->SetMarkerColor(2);
    
    // epsilon_31 = 0.4
    TH1D *h_single_electron_1 = h_data_1->ProjectionX("h_single_electron_1");
    h_single_electron_1->Scale(1.0 / (Double_t)dimension_xy);
    /*for(Int_t i{1}; i <= h_single_electron_1->GetNbinsX(); ++ i)
    {
        Double_t F{1.0 / (Double_t)dimension_xy};
        h_single_electron_1->SetBinContent(i, F * h_single_electron_1->GetBinContent(i));
    }*/
    h_single_electron_1->SetStats(0);
    h_single_electron_1->SetLineColor(3);
    h_single_electron_0->SetMarkerColor(2);
    
    // epsilon_31 = 0.8
    TH1D *h_single_electron_2 = h_data_2->ProjectionX("h_single_electron_2");
    h_single_electron_2->Scale(1.0 / (Double_t)dimension_xy);
    /*for(Int_t i{1}; i <= h_single_electron_2->GetNbinsX(); ++ i)
    {
        Double_t F{1.0 / (Double_t)dimension_xy};
        h_single_electron_2->SetBinContent(i, F * h_single_electron_2->GetBinContent(i));
    }*/
    h_single_electron_2->SetStats(0);
    h_single_electron_2->SetLineColor(4);
    h_single_electron_0->SetMarkerColor(4);

    /*
    TH1D *h_single_electron_0 = new TH1D("h_single_electron_0", "", dimension_xy, 0.0, 1.0);
    h_single_electron_0->SetStats(0);

    for(Int_t i{1}; i <= h_ratio->GetNbinsX(); ++ i)
    {
        Double_t content{0.0};
        for(Int_t j{1}; j <= h_ratio->GetNbinsY(); ++ j)
        {
            content += h_ratio->GetBinContent(i, j);
        }
        h_single_electron_0->SetBinContent(i, content);
    }
    */

    TCanvas *c_single_electron_0 = new TCanvas("c_single_electron_0", "", 4000, 3000);
    h_single_electron_0->Draw("hist");
    c_single_electron_0->SaveAs("c_single_electron_0.png");
    c_single_electron_0->SaveAs("c_single_electron_0.pdf");
    c_single_electron_0->SaveAs("c_single_electron_0.C");
    delete c_single_electron_0;
    
    TCanvas *c_single_electron_1 = new TCanvas("c_single_electron_1", "", 4000, 3000);
    h_single_electron_1->Draw("hist");
    c_single_electron_1->SaveAs("c_single_electron_1.png");
    c_single_electron_1->SaveAs("c_single_electron_1.pdf");
    c_single_electron_1->SaveAs("c_single_electron_1.C");
    delete c_single_electron_1;
    
    TCanvas *c_single_electron_2 = new TCanvas("c_single_electron_2", "", 4000, 3000);
    h_single_electron_2->Draw("hist");
    c_single_electron_2->SaveAs("c_single_electron_2.png");
    c_single_electron_2->SaveAs("c_single_electron_2.pdf");
    c_single_electron_2->SaveAs("c_single_electron_2.C");
    delete c_single_electron_2;

    TCanvas *c_single_electron = new TCanvas("c_single_electron", "", 4000, 3000);
    h_single_electron_0->Draw("hist");
    h_single_electron_1->Draw("histsame");
    h_single_electron_2->Draw("histsame");
    c_single_electron->SaveAs("c_single_electron.png");
    c_single_electron->SaveAs("c_single_electron.pdf");
    c_single_electron->SaveAs("c_single_electron.C");
    delete c_single_electron;

    ////////////////////////////////////////////////////////////////////////////
    // ELECTRON ENERGY CONVERSION
    ////////////////////////////////////////////////////////////////////////////

    TH2D *h_convert_from{nullptr};
    TH2D *h_convert_to{nullptr};

    // convert from epsilon_31 = 0.0 to epsilon_31 = 0.8 
    h_convert_from = h_data_0;
    h_convert_to = h_data_2;

    // electron energy in units of Q value
    // the sum must be less than or equal to 1.0
    // (1.0 - E1 - E2 = neutrino energy)
    Double_t electron_energy_1{0.33};
    Double_t electron_energy_2{0.33};

    // get weights
    Double_t weight_in{h_convert_from->GetBinContent(h_convert_from->FindBin(electron_energy_1, electron_energy_2))};
    Double_t weight_out{h_convert_to->GetBinContent(h_convert_to->FindBin(electron_energy_1, electron_energy_2))};
    Double_t weight_factor{weight_out / weight_in};

    std::cout << "Weight factor is " << weight_factor << std::endl;

    TH1D *h_single_electron_test = new TH1D("h_single_electron_test", "", dimension_xy, 0.0, 1.0);
    h_single_electron_test->SetStats(0);
    for(Int_t i{1}; i <= h_data_0->GetNbinsX(); ++ i)
    {
        Double_t energy_1{h_data_0->GetXaxis()->GetBinCenter(i)};
        Double_t weight_integral_1{0.0};
        Double_t weight_integral_2{0.0};
        for(Int_t j{1}; j <= h_data_0->GetNbinsY(); ++ j)
        {
            Double_t energy_2{h_data_0->GetXaxis()->GetBinCenter(j)};
            if(energy_1 + energy_2 <= 1.0)
            {
                weight_integral_1 += h_data_0->GetBinContent(i, j);
                weight_integral_2 += h_data_2->GetBinContent(i, j);
            }
        }
        if(weight_integral_1 > 0.0)
        {
            Double_t weight_factor{weight_integral_2 / weight_integral_1};
            Double_t input_weight{h_single_electron_0->GetBinContent(i)};
            Double_t output_weight{input_weight * weight_factor};
            h_single_electron_test->SetBinContent(i, output_weight);
        }
    }

 
    // test output, with histogram created from reweighting
    TCanvas *c_single_electron_test = new TCanvas("c_single_electron_test", "", 4000, 3000);
    h_single_electron_test->Draw("hist");
    c_single_electron_test->SaveAs("c_single_electron_test.png");
    c_single_electron_test->SaveAs("c_single_electron_test.pdf");
    c_single_electron_test->SaveAs("c_single_electron_test.C");
    delete c_single_electron_test;

    // compare above histogram (created from reweighting)
    // with histograms created from profile
    TCanvas *c_single_electron_with_test = new TCanvas("c_single_electron_with_test", "", 4000, 3000);
    h_single_electron_0->Draw("hist");
    h_single_electron_1->Draw("histsame");
    h_single_electron_2->Draw("histsame");
    h_single_electron_test->Draw("histsame");
    c_single_electron_with_test->SaveAs("c_single_electron_with_test.png");
    c_single_electron_with_test->SaveAs("c_single_electron_with_test.pdf");
    c_single_electron_with_test->SaveAs("c_single_electron_with_test.C");
    delete c_single_electron_with_test;

    //ReWeight(0.5, 0.5, 0.0, data_nEqNull, data_nEqTwo); 
    //ReWeight(0.45, 0.3, 0.8, h_nEqNull, h_nEqTwo, psiN0, psiN2); 
    
    TFile *f = new TFile("

    ReWeight(0.45, 0.3, 0.8, h_nEqNull, h_nEqTwo, psiN0, psiN2); 

    return 0;

}
