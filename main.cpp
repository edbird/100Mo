


#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"



#include <iostream>


#include "ReWeight.hpp"
#include "read_data.hpp"

// 0.0 MeV to 4.0 MeV = 4.0 MeV range
// num_bins keV bin width: 4.0 MeV / 0.1 MeV = 40 bins
Int_t num_bins{40};

int main(int argc, char* argv[])
{
    // Q value of decay
    // MeV
    Double_t bb_Q{3.034};

    ////////////////////////////////////////////////////////////////////////////
    // READ DATA IN FROM FILES
    ////////////////////////////////////////////////////////////////////////////

    // allocate data storage
    std::vector<std::vector<double>> data_nEqNull;
    std::vector<std::vector<double>> data_nEqTwo;

    double psiN0;
    double psiN2;

    // read
    read_data_helper("nEqNull.dat", "nEqTwo.dat", "psiN0.txt", "psiN2.txt", data_nEqNull, data_nEqTwo, psiN0, psiN2);

    std::cout << psiN0 << std::endl;
    std::cout << psiN2 << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    // READ DATA FOR ELECTRON ENERGY SINGLE AND SUM HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////

    // allocate data storage
    std::vector<std::vector<double>> data_electron_energy_single_0_0;
    std::vector<std::vector<double>> data_electron_energy_single_0_4;
    std::vector<std::vector<double>> data_electron_energy_single_0_8;
    std::vector<std::vector<double>> data_electron_energy_sum_0_0;
    std::vector<std::vector<double>> data_electron_energy_sum_0_4;
    std::vector<std::vector<double>> data_electron_energy_sum_0_8;

    // read
    read_data_helper_2("./data/mo100/single0.0.dat", data_electron_energy_single_0_0);
    read_data_helper_2("./data/mo100/single0.4.dat", data_electron_energy_single_0_4);
    read_data_helper_2("./data/mo100/single0.8.dat", data_electron_energy_single_0_8);
    read_data_helper_2("./data/mo100/sum0.0.dat", data_electron_energy_sum_0_0);
    read_data_helper_2("./data/mo100/sum0.4.dat", data_electron_energy_sum_0_4);
    read_data_helper_2("./data/mo100/sum0.8.dat", data_electron_energy_sum_0_8);

    ////////////////////////////////////////////////////////////////////////////
    // CONVERT INPUT DATA TO GRAPH FORMAT
    ////////////////////////////////////////////////////////////////////////////

    Int_t data_size{data_electron_energy_single_0_0.size()};

    if(data_electron_energy_single_0_4.size() != data_size)
        throw "error: size";
    
    if(data_electron_energy_single_0_8.size() != data_size)
        throw "error: size";
    
    if(data_electron_energy_sum_0_0.size() != data_size)
        throw "error: size";
    
    if(data_electron_energy_sum_0_4.size() != data_size)
        throw "error: size";

    if(data_electron_energy_sum_0_8.size() != data_size)
        throw "error: size";

    Double_t *data_x = new Double_t[data_size];

    Double_t *data_single_0 = new Double_t[data_size];
    Double_t *data_single_1 = new Double_t[data_size];
    Double_t *data_single_2 = new Double_t[data_size];

    Double_t *data_sum_0 = new Double_t[data_size];
    Double_t *data_sum_1 = new Double_t[data_size];
    Double_t *data_sum_2 = new Double_t[data_size];

    for(std::size_t i{0}; i < data_size; ++ i)
        data_x[i] = bb_Q * data_electron_energy_single_0_0[i][0];

    for(std::size_t i{0}; i < data_size; ++ i)
        data_single_0[i] = (1.0 / bb_Q) * data_electron_energy_single_0_0[i][1];

    for(std::size_t i{0}; i < data_size; ++ i)
        data_single_1[i] = (1.0 / bb_Q) * data_electron_energy_single_0_4[i][1];
    
    for(std::size_t i{0}; i < data_size; ++ i)
        data_single_2[i] = (1.0 / bb_Q) * data_electron_energy_single_0_8[i][1];

    for(std::size_t i{0}; i < data_size; ++ i)
        data_sum_0[i] = (1.0 / bb_Q) * data_electron_energy_sum_0_0[i][1];

    for(std::size_t i{0}; i < data_size; ++ i)
        data_sum_1[i] = (1.0 / bb_Q) * data_electron_energy_sum_0_4[i][1];
    
    for(std::size_t i{0}; i < data_size; ++ i)
        data_sum_2[i] = (1.0 / bb_Q) * data_electron_energy_sum_0_8[i][1];

    TGraph *g_el_energy_single_0 = new TGraph(data_size, data_x, data_single_0);
    TGraph *g_el_energy_single_1 = new TGraph(data_size, data_x, data_single_1);
    TGraph *g_el_energy_single_2 = new TGraph(data_size, data_x, data_single_2);

    TGraph *g_el_energy_sum_0 = new TGraph(data_size, data_x, data_sum_0);
    TGraph *g_el_energy_sum_1 = new TGraph(data_size, data_x, data_sum_1);
    TGraph *g_el_energy_sum_2 = new TGraph(data_size, data_x, data_sum_2);

    // compute integral to test
    Double_t data_single_0_integral{0.0};
    for(Int_t i{0 + 1}; i < data_size - 1; ++ i)
    {
        data_single_0_integral += data_single_0[i];
    }
    data_single_0_integral *= 2.0;
    data_single_0_integral += data_single_0[0];
    data_single_0_integral += data_single_0[data_size - 1];
    data_single_0_integral *= (data_x[1] - data_x[0]) / 2.0;
    std::cout << "data_single_0_integral=" << data_single_0_integral << std::endl;


    ////////////////////////////////////////////////////////////////////////////
    // TODO: apply phase space factor to data
    ////////////////////////////////////////////////////////////////////////////

    // ...

    ////////////////////////////////////////////////////////////////////////////
    // CONVERT INPUT DATA TO HISTOGRAM FORMAT
    // Note: Added 2018-04-23 (After INTERMEDIATE DATA below)
    // Note: These histograms do NOT have the phase space variable included
    ////////////////////////////////////////////////////////////////////////////

    const Int_t dimension_xy{1001};
    // don't plot raw data
    TH2D *h_nEqNull = nullptr; //= new TH2D("h_nEqNull", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    TH2D *h_nEqTwo = nullptr; //= new TH2D("h_nEqTwo", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_ratio = new TH2D("h_ratio", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //h_nEqNull->SetStats(0);
    //h_nEqTwo->SetStats(0);
    //h_ratio->SetStats(0);
    convert_data_to_histogram_format(data_nEqNull, data_nEqTwo, h_nEqNull, h_nEqTwo);

    ////////////////////////////////////////////////////////////////////////////
    // CREATE INTERMEDIATE DATA
    // Format: Nx3 array, each array a different value of epsilon
    ////////////////////////////////////////////////////////////////////////////

    // create data array for "complete data"
    // (with phase space factors)
    
    // epsilon_31 = 0.0
    std::vector<std::vector<double>> data_0;
    create_data_with_phase_space_factor(data_0, 0.0, data_nEqNull, data_nEqTwo, psiN0, psiN2);

    // epsilon_31 = 0.4
    std::vector<std::vector<double>> data_1;
    create_data_with_phase_space_factor(data_1, 0.4, data_nEqNull, data_nEqTwo, psiN0, psiN2);
    
    // epsilon_31 = 0.8
    std::vector<std::vector<double>> data_2;
    create_data_with_phase_space_factor(data_2, 0.8, data_nEqNull, data_nEqTwo, psiN0, psiN2);

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
            h_data_0->SetBinContent(i + 1, j + 1, data_0.at(i * dimension_xy + j)[2]);
            h_data_1->SetBinContent(i + 1, j + 1, data_1.at(i * dimension_xy + j)[2]);
            h_data_2->SetBinContent(i + 1, j + 1, data_2.at(i * dimension_xy + j)[2]);
            //if(i < dimension_xy - j)
            // TODO: -1 or -0 ?
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
    
    #define PRINT_DATA_CANVAS 0
    #if PRINT_DATA_CANVAS

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


/*
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
*/

    // these histograms created from projection of the h_data_? histograms
    // along the X direction
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

    // same as the single electron histogram,
    // however it is made using an explicit integration method
    // results appear to be the same
    /*
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
    */


    // compare above histogram (created from reweighting)
    // with histograms created from profile
    // NOTE: removed as this code is exactly the same as c_single_electron
    /*
    TCanvas *c_single_electron_with_test = new TCanvas("c_single_electron_with_test", "", 4000, 3000);
    h_single_electron_0->Draw("hist");
    h_single_electron_1->Draw("histsame");
    h_single_electron_2->Draw("histsame");
    h_single_electron_test->Draw("histsame");
    c_single_electron_with_test->SaveAs("c_single_electron_with_test.png");
    c_single_electron_with_test->SaveAs("c_single_electron_with_test.pdf");
    c_single_electron_with_test->SaveAs("c_single_electron_with_test.C");
    delete c_single_electron_with_test;
    */

    //ReWeight(0.5, 0.5, 0.0, data_nEqNull, data_nEqTwo); 
    //ReWeight(0.45, 0.3, 0.8, h_nEqNull, h_nEqTwo, psiN0, psiN2); 
    

    // load from NEMO-3 data (MC), the reconstructed electron energy (2 in same histogram)
    TH1D *h_el_energy_original = new TH1D("h_el_energy_original", "", num_bins, 0.0, 4.0);
    h_el_energy_original->SetStats(0);
    h_el_energy_original->SetLineColor(2);
    h_el_energy_original->SetMarkerColor(2);
    // same as above but re-weighted
    TH1D *h_el_energy_reweight = new TH1D("h_el_energy_reweight", "", num_bins, 0.0, 4.0);
    h_el_energy_reweight->SetStats(0);
    h_el_energy_reweight->SetLineColor(3);
    h_el_energy_reweight->SetMarkerColor(3);

    // load from NEMO-3 data (mc), the sum of the two reconstructed electron energies
    TH1D *h_el_energy_sum_original = new TH1D("h_el_energy_sum_original", "", num_bins, 0.0, 4.0);
    h_el_energy_sum_original->SetStats(0);
    h_el_energy_sum_original->SetLineColor(2);
    h_el_energy_sum_original->SetMarkerColor(2);
    // same as above but re-weighted
    TH1D *h_el_energy_sum_reweight = new TH1D("h_el_energy_sum_reweight", "", num_bins, 0.0, 4.0);
    h_el_energy_sum_reweight->SetStats(0);
    h_el_energy_sum_reweight->SetLineColor(3);
    h_el_energy_sum_reweight->SetMarkerColor(3);

    // test histograms
    // NEMO-3 data/MC: truth electron energy x2 (T1, T2) (2 in same histo)
    TH1D *h_test_single_original = new TH1D("h_test_single_original", "", num_bins, 0.0, 4.0);
    h_test_single_original->SetStats(0);
    h_test_single_original->SetLineColor(2);
    h_test_single_original->SetMarkerColor(2);
    // re-weighted
    TH1D *h_test_single_reweight = new TH1D("h_test_single_reweight", "", num_bins, 0.0, 4.0);
    h_test_single_reweight->SetStats(0);
    h_test_single_reweight->SetLineColor(3);
    h_test_single_reweight->SetMarkerColor(3);

    // test histograms
    // NEMO-3 data/MC: reconstructed electron energy x2 (2 in same histo)
    TH1D *h_test_sum_original = new TH1D("h_test_sum_original", "", num_bins, 0.0, 4.0);
    h_test_sum_original->SetStats(0);
    h_test_sum_original->SetLineColor(2);
    h_test_sum_original->SetMarkerColor(2);
    // re-weighted
    TH1D *h_test_sum_reweight = new TH1D("h_test_sum_reweight", "", num_bins, 0.0, 4.0);
    h_test_sum_reweight->SetStats(0);
    h_test_sum_reweight->SetLineColor(3);
    h_test_sum_reweight->SetMarkerColor(3);


    // write data histograms to file
    // intermediate file
    TFile *f_histogram = new TFile("f_histogram.root", "recreate");
    h_data_0->Write();
    h_data_1->Write();
    h_data_2->Write();
    f_histogram->Close();


    // NEMO-3 data/MC read from tree
    // input file

    TH2D *h_gen_weight = new TH2D("h_gen_weight", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_gen_weight->SetStats(0);
    
    // read file using program arguments
    std::string filename_default("NewElectronNtuplizerExe_Int_ManDB_output.root");
    std::string filename(filename_default);
    bool gen_weight_enable{false};
    if(argc == 2)
    {
        filename = std::string(argv[1]);
        std::cout << "Specified filename: " << filename << std::endl;
        if(filename != filename_default)
        {
            gen_weight_enable = true;
        }
    }
    //TFile *f = new TFile("NewElectronNtuplizerExe_Int_ManDB_output.root");
    TFile *f = new TFile(filename.c_str());
    TTree *t = (TTree*)f->Get("NewElectronNtuplizer/NewElectronNtuplizer");

    Int_t nElectrons;
    Double_t trueT1;
    Double_t trueT2;
    Double_t el_energy_[2];
    Double_t gen_weight{1.0};

    t->SetBranchAddress("nElectrons", &nElectrons);
    t->SetBranchAddress("trueT1", &trueT1);
    t->SetBranchAddress("trueT2", &trueT2);
    t->SetBranchAddress("el_energy_", el_energy_);
    // new method requires weight to be saved to tree
    if(gen_weight_enable == true)
    {
        t->SetBranchAddress("gen_weight", &gen_weight);
    }


    std::cout << "Processing data" << std::endl;
    Long64_t prog_c{-1};
    for(Long64_t ix{0}; ix < t->GetEntries(); ++ ix)
    {

        t->GetEntry(ix);

        if(nElectrons != 2) continue;


        //std::cout << "trueT1=" << trueT1 << " trueT2=" << trueT2;
        //std::cout << " -> " << ReWeight(trueT1 / bb_Q, trueT2 / bb_Q, 0.8, h_nEqNull, h_nEqTwo, psiN0, psiN2) << std::endl;

        const Double_t epsilon_31{0.8};

        Double_t T1{trueT1 / bb_Q};
        Double_t T2{trueT2 / bb_Q};

        Double_t weight{ReWeight(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};

        h_el_energy_original->Fill(el_energy_[0], 1.0 * gen_weight);
        h_el_energy_original->Fill(el_energy_[1], 1.0 * gen_weight);
        
        h_el_energy_sum_original->Fill(el_energy_[0] + el_energy_[1], 1.0 * gen_weight);

        // test histograms
        h_test_single_original->Fill(trueT1, 1.0 * gen_weight);
        h_test_single_original->Fill(trueT2, 1.0 * gen_weight);

        h_test_sum_original->Fill(trueT1 + trueT2, 1.0 * gen_weight);

        /*
        if(!(weight < 3.0 && weight > 0.001))
        {
            std::cout << "ix=" << ix << std::endl;
            std::cout << "trueT1=" << trueT1 << " trueT2=" << trueT2 << " sum=" << trueT1 + trueT2 << std::endl;

            ReWeight(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true");
            
            std::cout << "weight=" << weight << std::endl;
            std::cin.get();
        }
        */

        h_el_energy_reweight->Fill(el_energy_[0], weight * gen_weight);
        h_el_energy_reweight->Fill(el_energy_[1], weight * gen_weight);

        h_el_energy_sum_reweight->Fill(el_energy_[0] + el_energy_[1], weight * gen_weight);

        // test histograms
        h_test_single_reweight->Fill(trueT1, weight * gen_weight);
        h_test_single_reweight->Fill(trueT2, weight * gen_weight);

        h_test_sum_reweight->Fill(trueT1 + trueT2, weight * gen_weight);

        h_gen_weight->Fill(trueT1, trueT2, gen_weight);

        //if(ix % 10 == 9) std::cin.get();

        if((ix * 10) / t->GetEntries() > prog_c)
        {
            prog_c = (ix * 10) / t->GetEntries();
            std::cout << 10 * (ix * 10) / t->GetEntries() << " %" << std::endl;
        }

    }
    
    // rescale histograms pre other re-scaling
    h_el_energy_original->Scale((1.0 / 0.1) / h_el_energy_original->Integral());
    h_el_energy_reweight->Scale((1.0 / 0.1) / h_el_energy_reweight->Integral());
    h_el_energy_sum_original->Scale((1.0 / 0.1) / h_el_energy_sum_original->Integral());
    h_el_energy_sum_reweight->Scale((1.0 / 0.1) / h_el_energy_sum_reweight->Integral());

    h_test_single_original->Scale((1.0 / 0.1) / h_test_single_original->Integral());
    h_test_single_reweight->Scale((1.0 / 0.1) / h_test_single_reweight->Integral());
    h_test_sum_original->Scale((1.0 / 0.1) / h_test_sum_original->Integral());
    h_test_sum_reweight->Scale((1.0 / 0.1) / h_test_sum_reweight->Integral());

    std::cout << "h_el_energy_original_integral=" << h_el_energy_original->Integral() << std::endl;

    // calculate integral manually
    Double_t h_el_energy_original_integral{0.0};
    for(Int_t i{1}; i <= h_el_energy_original->GetNbinsX(); ++ i)
    {
        h_el_energy_original_integral += h_el_energy_original->GetBinContent(i);
    }
    std::cout << "h_el_energy_original_integral=" << h_el_energy_original_integral << std::endl;

    TCanvas *c_gen_weight = new TCanvas("c_gen_weight", "", 800, 600);
    h_gen_weight->Draw("colz");
    c_gen_weight->SaveAs("c_gen_weight.C");
    c_gen_weight->SaveAs("c_gen_weight.png");
    c_gen_weight->SaveAs("c_gen_weight.pdf");
    delete c_gen_weight;

    // scale the green histogram to match the red one (for sum energy histo)
    Double_t integral_1{h_el_energy_sum_original->Integral()};
    Double_t integral_2{h_el_energy_sum_reweight->Integral()};
    h_el_energy_sum_reweight->Scale(integral_1 / integral_2);
    h_el_energy_reweight->Scale(integral_1 / integral_2);
    Double_t chi_square{chi_square_test(h_el_energy_reweight, h_el_energy_original)};
    std::cout << "chi_square=" << chi_square << std::endl;

    /*
    TCanvas *c_el_energy_original = new TCanvas("c_el_energy_original", "c_el_energy_original", 800, 600);
    h_el_energy_original->Draw("E");
    c_el_energy_original->SaveAs("c_el_energy_original.C");
    c_el_energy_original->SaveAs("c_el_energy_original.png");
    c_el_energy_original->SaveAs("c_el_energy_original.pdf");
    delete c_el_energy_original;

    TCanvas *c_el_energy_reweight = new TCanvas("c_el_energy_reweight", "c_el_energy_reweight", 800, 600);
    h_el_energy_reweight->Draw("E");
    c_el_energy_reweight->SaveAs("c_el_energy_reweight.C");
    c_el_energy_reweight->SaveAs("c_el_energy_reweight.png");
    c_el_energy_reweight->SaveAs("c_el_energy_reweight.pdf");
    delete c_el_energy_reweight; 
    */

    // print single electron distribution
    TCanvas *c_el_energy_both = new TCanvas("e_el_energy_both", "e_el_energy_both", 800, 600);
    h_el_energy_original->Draw("E");
    h_el_energy_reweight->Draw("Esame");
    c_el_energy_both->SaveAs("c_el_energy_both.C");
    c_el_energy_both->SaveAs("c_el_energy_both.png");
    c_el_energy_both->SaveAs("c_el_energy_both.pdf");
    delete c_el_energy_both;

    // print summed distribution
    TCanvas *c_el_energy_sum_both = new TCanvas("e_el_energy_sum_both", "e_el_energy_sum_both", 800, 600);
    h_el_energy_sum_original->Draw("E");
    h_el_energy_sum_reweight->Draw("Esame");
    c_el_energy_sum_both->SaveAs("c_el_energy_sum_both.C");
    c_el_energy_sum_both->SaveAs("c_el_energy_sum_both.png");
    c_el_energy_sum_both->SaveAs("c_el_energy_sum_both.pdf");
    delete c_el_energy_sum_both; 

    // print single electron distribution test histograms
    TCanvas *c_test_single = new TCanvas("c_test_single", "c_test_single", 800, 600);
    h_test_single_original->Draw("E");
    h_test_single_reweight->Draw("Esame");
    g_el_energy_single_0->Draw("same");
    g_el_energy_single_1->Draw("same");
    g_el_energy_single_2->Draw("same");
    //h_test_single_reweight->Draw("E");
    c_test_single->SaveAs("c_test_single.C");
    c_test_single->SaveAs("c_test_single.png");
    c_test_single->SaveAs("c_test_single.pdf");
    delete c_test_single;
    
    // print summed distribution test histograms
    TCanvas *c_test_sum = new TCanvas("c_test_sum", "c_test_sum", 800, 600);
    h_test_sum_original->Draw("E");
    h_test_sum_reweight->Draw("Esame");
    g_el_energy_sum_0->Draw("same");
    g_el_energy_sum_1->Draw("same");
    g_el_energy_sum_2->Draw("same");
    //h_test_sum_reweight->Draw("E");
    c_test_sum->SaveAs("c_test_sum.C");
    c_test_sum->SaveAs("c_test_sum.png");
    c_test_sum->SaveAs("c_test_sum.pdf");
    delete c_test_sum;




        

    // print entries
    std::cout << "Number of entries in each histogram: h_el_energy_original: " << h_el_energy_original->GetEntries() << " h_el_energy_reweight: " << h_el_energy_reweight->GetEntries() << std::endl;



    return 0;

}
