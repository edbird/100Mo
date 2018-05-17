

#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"



#include <iostream>


#include "ReWeight.hpp"
#include "read_data.hpp"

#include "aux.hpp"

#include "program_arguments.hpp"

#include "CanvasFactory.hpp"

//#include "HistogramWrapper.hpp"
// TODO: make some nice histograms
// X/Y AXIS: LABEL: TEXT, FONT SIZE, FONT
// X/Y AXIS: TICK NUMBERS: FONT SIZE, FONT
// X/Y AXIS LIMITS, RANGE, PARTICULARLY Y MAX
// GRAPH COLORS



// 0.0 MeV to 4.0 MeV = 4.0 MeV range
// num_bins keV bin width: 4.0 MeV / 0.1 MeV = 40 bins
Int_t num_bins{40};

int main(int argc, char* argv[])
{
    
    ////////////////////////////////////////////////////////////////////////////
    // PROCESS PROGRAM ARGUMENTS
    ////////////////////////////////////////////////////////////////////////////

    // new program arguments processing method
    ProgramArguments pa;
    pa.Add("help", "--help", "false");
    pa.Add("filename", "--filename", "NewElectronNtuplizerExe_Int_ManDB_output.root");
    pa.Add("epsilon", "--epsilon", "0.0");
    pa.Add("batch_mode", "--batch-mode", "false");
    //pa.Add("300 keV energy cut", "--energy-cut", "false");
    pa.Add("energy_cut", "--energy-cut", "false");
    pa.Add("fit_subrange", "--fit-subrange", "false");
    pa.Add("log_mode", "--log-mode", "true");
    pa.Add("output_filename", "--output-file", "of_data.txt");
    pa.Print();
    pa.Process(argc, argv);

    // this argument is a string
    std::string help{pa.Get("help")};
    // this argument is a string
    std::string filename{pa.Get("filename")};
    // this argument is to be converted to a double
    std::string arg_epsilon_31{pa.Get("epsilon")};
    // this argument is to be converted to a bool
    std::string arg_batch_mode{pa.Get("batch_mode")};
    std::string arg_energy_cut{pa.Get("energy_cut")};
    std::string arg_fit_subrange{pa.Get("fit_subrange")};
    std::string arg_log_mode{pa.Get("log_mode")};
    std::string arg_output_filename{pa.Get("output_filename")};

    double epsilon_31{std::stod(arg_epsilon_31)};
    bool batch_mode{false};
    if(arg_batch_mode == std::string("true"))
    {
        batch_mode = true;
        std::cout << "[ INFO ] : Batch mode: " << "true" << std::endl;
    }
    else
    {
        std::cout << "[ INFO ] : Batch mode: " << "false" << std::endl;
    }

    if(arg_energy_cut == std::string("true"))
    {
        std::cout << "[ INFO ] : Energy cut: " << "true" << std::endl;
    }
    else
    {
        std::cout << "[ INFO ] : Energy cut: " << "false" << std::endl;
    }

    if(arg_fit_subrange == std::string("true"))
    {
        std::cout << "[ INFO ] : Fit subrange: " << "true" << std::endl;
    }
    else
    {
        std::cout << "[ INFO ] : Fit subrange: " << "false" << std::endl;
    }

    if(arg_log_mode == std::string("true"))
    {
        std::cout << "[ INFO ] : Log mode: " << "true" << std::endl;
    }
    else
    {
        std::cout << "[ INFO ] : Log mode: " << "false" << std::endl;
    }

    std::cout << "Writing data to file " << arg_output_filename << " append mode" << std::endl;

    // process gathered argument data
    //filename = arg_filename;
    //std::string filename_default(arg_default_filename); // required for later if test
    bool gen_weight_enable{false};
    //if(arg_filename != filename_default)
    if(filename != pa.GetDefault("filename"))
    {
        gen_weight_enable = true;
    }
    // TODO: check filename exists!

    //epsilon_31 = std::stod(arg_epsilon_31);
    // todo: check valid
    if(0.0 <= epsilon_31 && epsilon_31 <= 10.0)
    {
        std::cout << "[ INFO ] : Set epsilon_31 = " << epsilon_31 << std::endl;
    }
    else
    {
        std::cout << "invalid epsilon_31 value" << std::endl;
        throw "invalid epsilon_31 value";
    }

    /*
    if(arg_batch_mode == std::string("true"))
    {
        batch_mode = true;
        std::string batch_mode_enable_string("false");
        if(batch_mode)
        {
            batch_mode_enable_string = std::string("true");
        }
        std::cout << "[ INFO ] : Batch mode: " << batch_mode_enable_string << std::endl;
    }
    */




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

    // SCALE THE DATA SO THAT X RUNS FROM 0.0 TO Q_BB RATHER THAN 0.0 TO 1.0
    // SCALE THE DATA SO THAT THE INTEGRAL IS UNCHANGED, IE: INTEGRAL = 1.0
    // THIS IS DONE BY SCALING THE Y VALUES IN THE INVERSE WAY TO THE SCALING
    // OF THE X VALUES

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
    g_el_energy_single_0->SetLineColor(2);

    TGraph *g_el_energy_single_1 = new TGraph(data_size, data_x, data_single_1);
    g_el_energy_single_1->SetLineColor(3);
    
    TGraph *g_el_energy_single_2 = new TGraph(data_size, data_x, data_single_2);
    g_el_energy_single_2->SetLineColor(4);

    TGraph *g_el_energy_sum_0 = new TGraph(data_size, data_x, data_sum_0);
    g_el_energy_sum_0->SetLineColor(2);
    
    TGraph *g_el_energy_sum_1 = new TGraph(data_size, data_x, data_sum_1);
    g_el_energy_sum_1->SetLineColor(3);
    
    TGraph *g_el_energy_sum_2 = new TGraph(data_size, data_x, data_sum_2);
    g_el_energy_sum_2->SetLineColor(4);

    // compute integral to test
    // trapezium rule
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
    // Note: Do not do this

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

        if(batch_mode == false)
        {
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
        }

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
    h_single_electron_1->SetMarkerColor(3);
    
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
    h_single_electron_2->SetMarkerColor(4);

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
    ->SetMaximum();
    h_single_electron_0->Draw("hist");
    c_single_electron_0->SaveAs("c_single_electron_0.png");
    c_single_electron_0->SaveAs("c_single_electron_0.pdf");
    c_single_electron_0->SaveAs("c_single_electron_0.C");
    delete c_single_electron_0;
    
    TCanvas *c_single_electron_1 = new TCanvas("c_single_electron_1", "", 4000, 3000);
    ->SetMaximum();
    h_single_electron_1->Draw("hist");
    c_single_electron_1->SaveAs("c_single_electron_1.png");
    c_single_electron_1->SaveAs("c_single_electron_1.pdf");
    c_single_electron_1->SaveAs("c_single_electron_1.C");
    delete c_single_electron_1;
    
    TCanvas *c_single_electron_2 = new TCanvas("c_single_electron_2", "", 4000, 3000);
    ->SetMaximum();
    h_single_electron_2->Draw("hist");
    c_single_electron_2->SaveAs("c_single_electron_2.png");
    c_single_electron_2->SaveAs("c_single_electron_2.pdf");
    c_single_electron_2->SaveAs("c_single_electron_2.C");
    delete c_single_electron_2;
*/

    // these histograms created from projection of the h_data_? histograms
    // along the X direction
    #if PRINT_DATA_CANVAS
        if(batch_mode == false)
        {
            TCanvas *c_single_electron = new TCanvas("c_single_electron", "", 4000, 3000);
            h_single_electron_0->SetMaximum(3.5);
            h_single_electron_0->GetXaxis()->SetTitle("Energy [MeV]");
            h_single_electron_0->GetYaxis()->SetTitle("Events");
            h_single_electron_0->Draw("hist");
            h_single_electron_1->Draw("histsame");
            h_single_electron_2->Draw("histsame");
            c_single_electron->SaveAs("c_single_electron.png");
            c_single_electron->SaveAs("c_single_electron.pdf");
            c_single_electron->SaveAs("c_single_electron.C");
            delete c_single_electron;
        }
    #endif

    ////////////////////////////////////////////////////////////////////////////
    // ELECTRON ENERGY CONVERSION
    ////////////////////////////////////////////////////////////////////////////

    // This is a test code

    /*
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
    */

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
        h_single_electron_0->GetYaxis()->SetMaximum();
        h_single_electron_0->GetXaxis()->SetTitle("Energy [MeV]");
        h_single_electron_0->GetYaxis()->SetTitle("Events");
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
        h_single_electron_0->GetYaxis()->SetMaximum();
        h_single_electron_0->GetXaxis()->SetTitle("Energy [MeV]");
        h_single_electron_0->GetYaxis()->SetTitle("Events");
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


    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS FOR SINGLE AND SUM ENERGY HISTOGRAMS FOR BASELINE (EPS=BASE)
    // AND RESCALED HISTOGRAMS (EPS=SOME OTHER VALUE)
    ////////////////////////////////////////////////////////////////////////////
    
    // load from NEMO-3 data (MC), the reconstructed electron energy (2 in same histogram)
    TH1D *h_el_energy_original = new TH1D("h_el_energy_original", "", num_bins, 0.0, 4.0);
    //h_el_energy_original->SetStats(0);
    //h_el_energy_original->SetLineColor(2);
    //h_el_energy_original->SetMarkerColor(2);
    // same as above but re-weighted
    TH1D *h_el_energy_reweight = new TH1D("h_el_energy_reweight", "", num_bins, 0.0, 4.0);
    //h_el_energy_reweight->SetStats(0);
    //h_el_energy_reweight->SetLineColor(4);
    //h_el_energy_reweight->SetMarkerColor(4);

    // load from NEMO-3 data (mc), the sum of the two reconstructed electron energies
    TH1D *h_el_energy_sum_original = new TH1D("h_el_energy_sum_original", "", num_bins, 0.0, 4.0);
    h_el_energy_sum_original->SetStats(0);
    h_el_energy_sum_original->SetLineColor(2);
    h_el_energy_sum_original->SetMarkerColor(2);
    // same as above but re-weighted
    TH1D *h_el_energy_sum_reweight = new TH1D("h_el_energy_sum_reweight", "", num_bins, 0.0, 4.0);
    h_el_energy_sum_reweight->SetStats(0);
    h_el_energy_sum_reweight->SetLineColor(4);
    h_el_energy_sum_reweight->SetMarkerColor(4);

    // test histograms
    // NEMO-3 data/MC: truth electron energy x2 (T1, T2) (2 in same histo)
    TH1D *h_test_single_original = new TH1D("h_test_single_original", "", num_bins, 0.0, 4.0);
    h_test_single_original->SetStats(0);
    h_test_single_original->SetLineColor(2);
    h_test_single_original->SetMarkerColor(2);
    // re-weighted
    TH1D *h_test_single_reweight = new TH1D("h_test_single_reweight", "", num_bins, 0.0, 4.0);
    h_test_single_reweight->SetStats(0);
    h_test_single_reweight->SetLineColor(4);
    h_test_single_reweight->SetMarkerColor(4);

    // test histograms
    // NEMO-3 data/MC: reconstructed electron energy x2 (2 in same histo)
    TH1D *h_test_sum_original = new TH1D("h_test_sum_original", "", num_bins, 0.0, 4.0);
    h_test_sum_original->SetStats(0);
    h_test_sum_original->SetLineColor(2);
    h_test_sum_original->SetMarkerColor(2);
    // re-weighted
    TH1D *h_test_sum_reweight = new TH1D("h_test_sum_reweight", "", num_bins, 0.0, 4.0);
    h_test_sum_reweight->SetStats(0);
    h_test_sum_reweight->SetLineColor(4);
    h_test_sum_reweight->SetMarkerColor(4);


    // write data histograms to file
    // intermediate file
    // used for creating test histogram
    // with createtree code
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
    //std::string filename_default("NewElectronNtuplizerExe_Int_ManDB_output.root");
    //std::string filename(filename_default);
    //bool gen_weight_enable{false};
    //if(argc == 2)
    //{
    //    filename = std::string(argv[1]);
    //    std::cout << "Specified filename: " << filename << std::endl;
    //    if(filename != filename_default)
    //    {
    //        gen_weight_enable = true;
    //    }
    //}
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
    //const Double epsilon_31{0.0};
    for(Long64_t ix{0}; ix < t->GetEntries(); ++ ix)
    {

        // TODO remove
        //if(ix == 1000) break;

        t->GetEntry(ix);

        if(nElectrons != 2) continue;

        // cut both electrons > 300 keV
        if(arg_energy_cut == std::string("true"))
        {
            if(el_energy_[0] < 0.3) continue;
            if(el_energy_[1] < 0.3) continue;
        }

        //std::cout << "trueT1=" << trueT1 << " trueT2=" << trueT2;
        //std::cout << " -> " << ReWeight(trueT1 / bb_Q, trueT2 / bb_Q, 0.8, h_nEqNull, h_nEqTwo, psiN0, psiN2) << std::endl;

        //const Double_t epsilon_31{0.8};

        Double_t T1{trueT1 / bb_Q};
        Double_t T2{trueT2 / bb_Q};

        // ReWeight = baseline 0.0, ReWeight2 = baseline = 0.382
        Double_t weight{ReWeight2(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};

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
            //std::cout << 10 * (ix * 10) / t->GetEntries() << " %" << std::endl;
        }

    }
    
    // rescale histograms pre other re-scaling
    // NOTE: to get the same >>> amplitude <<<
    // histograms must be scaled by the bin width
    // the bin width for one set is 100 keV
    // the bin width for the other set is 1 MeV / 1000 = 1 keV (?)
    // factor 1000 discrepency?
    // TODO: resolve this issue
    //h_el_energy_original->Scale(1.0 / h_el_energy_original->Integral());
    //h_el_energy_reweight->Scale(1.0 / h_el_energy_reweight->Integral());
    //h_el_energy_sum_original->Scale(1.0 / h_el_energy_sum_original->Integral());
    //h_el_energy_sum_reweight->Scale(1.0 / h_el_energy_sum_reweight->Integral());

    // NOTE: must scale all by same amount here
    // TODO: removed scaling
    /*
    Double_t scale_factor{1.0 / h_el_energy_sum_original->Integral()};
    h_el_energy_original->Scale(scale_factor);
    h_el_energy_reweight->Scale(scale_factor);
    h_el_energy_sum_original->Scale(scale_factor);
    h_el_energy_sum_reweight->Scale(scale_factor);
    */
    
    //TH1D *h_el_energy_sum_reweight_clone = (TH1D*)h_el_energy_sum_reweight->Clone();
    //h_el_energy_sum_reweight_clone->SetLineColor(6);
    //h_el_energy_sum_reweight_clone->SetMarkerColor(6);
    
    //TH1D *h_el_energy_sum_original_clone = (TH1D*)h_el_energy_sum_original->Clone();
    //h_el_energy_sum_original_clone->SetLineColor(7);
    //h_el_energy_sum_original_clone->SetMarkerColor(7);

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

    if(batch_mode == false)
    {
        TCanvas *c_gen_weight = new TCanvas("c_gen_weight", "", 800, 600);
        h_gen_weight->Draw("colz");
        c_gen_weight->SaveAs("c_gen_weight.C");
        c_gen_weight->SaveAs("c_gen_weight.png");
        c_gen_weight->SaveAs("c_gen_weight.pdf");
        delete c_gen_weight;
    }

    // scale the green histogram to match the red one (for sum energy histo) 
    #define FIT_METHOD_1 0
    #if FIT_METHOD_1
        Double_t integral_1{h_el_energy_sum_original->Integral()};
        Double_t integral_2{h_el_energy_sum_reweight->Integral()};
        h_el_energy_sum_reweight->Scale(integral_1 / integral_2);
        h_el_energy_reweight->Scale(integral_1 / integral_2);
        Double_t chi_square{chi_square_test(h_el_energy_reweight, h_el_energy_original)};
        std::cout << "chi_square=" << chi_square << std::endl;
    #endif

    Double_t sensitivity_chisquare{0.0};
    bool fit_subrange{false};
    if(arg_fit_subrange == std::string("true")) fit_subrange = true;
    // number of degrees of freedom
    Int_t non_empty_bins{0};
    #define FIT_METHOD_2 1
    #if FIT_METHOD_2
        TF1 *f_el_energy_sum_original = new TF1("f_el_energy_sum_original", fit_function, 0.0, 4.0, 1 + 2 * (h_el_energy_sum_reweight->GetNbinsX() + 1));
        f_el_energy_sum_original->SetParameter(0, 1.0); // set initial amplitude parameter to 1.0
        f_el_energy_sum_original->SetNpx(100000);
        // set the "parameters" (constants)
        // these are the values from the histogram
        {
            Int_t ix{0};
            Int_t jx{0};
            Double_t par_ix{0.0};
            Double_t par_jx{0.0};

            Int_t nbx{h_el_energy_sum_reweight->GetNbinsX()};

            for(Int_t i{0}; i <= nbx; ++ i)
            {
                ix = 2 * i + 1;
                jx = 2 * i + 2;

                par_ix = h_el_energy_sum_reweight->GetXaxis()->GetBinLowEdge(i + 1);
                par_jx = h_el_energy_sum_reweight->GetBinContent(i + 1);

                //std::cout << "x=" << par_ix << " y=" << par_jx << std::endl;

                f_el_energy_sum_original->FixParameter(ix, par_ix);
                f_el_energy_sum_original->FixParameter(jx, par_jx);
            }

            ix = 2 * nbx + 3;
            ix = 2 * nbx + 4;
            
            par_ix = h_el_energy_sum_reweight->GetXaxis()->GetBinLowEdge(nbx + 1) + h_el_energy_sum_reweight->GetBinWidth(1);
            par_jx = 0.0; // content doesn't make sense here, don't use "overflow" by accident

            // final 2, marks the end point of final "bin"
            f_el_energy_sum_original->FixParameter(ix, par_ix);
            f_el_energy_sum_original->FixParameter(jx, par_jx);
            
            // TODO: am i putting the right information into the fitting function and am
            // i fitting the correct function?
            if(fit_subrange == false)
            {
                h_el_energy_sum_original->Fit("f_el_energy_sum_original", "0");
            }
            else if(fit_subrange == true)
            {
                h_el_energy_sum_original->Fit("f_el_energy_sum_original", "0", "", 2.0, 4.0);
            }

        }
        std::cout << "Amplitude parameter: " << f_el_energy_sum_original->GetParameter(0) << std::endl;
        std::cout << "                err: " << f_el_energy_sum_original->GetParError(0) << std::endl;
        //std::cout << "         chi square: " << f_el_energy_sum_original->GetChisquare() / h_el_energy_sum_original->GetNbinsX() << std::endl;
        std::cout << "         chi square: " << f_el_energy_sum_original->GetChisquare() << std::endl;
        for(Int_t i{1}; i <= h_el_energy_sum_original->GetNbinsX(); ++ i)
        {
            if(h_el_energy_sum_original->GetBinContent(i) != 0.0) ++ non_empty_bins;
        }
        std::cout << " degrees of freedom: " << non_empty_bins << std::endl;
        std::cout << " chi square reduced: " << f_el_energy_sum_original->GetChisquare() / (Double_t)non_empty_bins << std::endl;

        // scale the red histogram using the amplitude parameter
        h_el_energy_sum_reweight->Scale(f_el_energy_sum_original->GetParameter(0));
        h_el_energy_reweight->Scale(f_el_energy_sum_original->GetParameter(0));

        // get chi-square for single electron histograms
        if(fit_subrange == false)
        {
            sensitivity_chisquare = chi_square_test(h_el_energy_reweight, h_el_energy_original);
            std::cout << "chi square of single electron: " << sensitivity_chisquare << std::endl;
            std::cout << " degrees of freedom: " << non_empty_bins << std::endl;
            std::cout << " chi square reduced: " << sensitivity_chisquare / (Double_t)non_empty_bins << std::endl;

        }
        else if(fit_subrange == true)
        {
            sensitivity_chisquare = chi_square_test(h_el_energy_reweight, h_el_energy_original, 2.0, 4.0);
            std::cout << "chi square of single electron, 2.0 MeV - 4.0 MeV: " << sensitivity_chisquare << std::endl;
            std::cout << " degrees of freedom: " << non_empty_bins << std::endl;
            std::cout << " chi square reduced: " << sensitivity_chisquare / (Double_t)non_empty_bins << std::endl;
        }
    #endif
    
    // Note: never got this working
    #define FIT_METHOD_3 0
    #if FIT_METHOD_3
        Double_t amplitude{1.0};
        driver(h_el_energy_sum_original, h_el_energy_sum_reweight, amplitude);
    #endif


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
    
    // enable/disable printing with log canvas
    Bool_t log_mode{false};
    if(arg_log_mode == std::string("true")) log_mode = true;

    if(batch_mode == false)
    {
        //HistogramWrapper2 hwrap_el_energy("el_energy", "x_axis_label_text", "y_axis_label_text", "E");
        //hwrap_el_energy.SetLogMode(log_mode);
        //hwrap_el_energy.Add("h_el_energy_original", h_el_energy_original);
        //hwrap_el_energy.Canvas();
        // TODO: make some nice histograms
        // X/Y AXIS: LABEL: TEXT, FONT SIZE, FONT
        // X/Y AXIS: TICK NUMBERS: FONT SIZE, FONT
        // X/Y AXIS LIMITS, RANGE, PARTICULARLY Y MAX
        // GRAPH COLORS


        // print single electron distribution
        //TCanvas *c_el_energy_both = new TCanvas("e_el_energy_both", "e_el_energy_both", 800, 600);
        //c_el_energy_both->SetLogy(log_mode);
        //h_el_energy_original->SetMaximum(220.0e3); // MC data
        //h_el_energy_original->SetMaximum(1000.0e3); // MC data log mode
        //h_el_energy_original->GetXaxis()->SetTitle("Energy [MeV]");
        //h_el_energy_original->GetYaxis()->SetTitle("Events");
        //h_el_energy_reweight->SetMaximum(1000.0e3);
        //h_el_energy_reweight->GetXaxis()->SetTitle("Energy [MeV]");
        //h_el_energy_reweight->GetYaxis()->SetTitle("Events");
        //h_el_energy_original->Draw("E");
        //h_el_energy_reweight->Draw("Esame");
        //c_el_energy_both->SaveAs("c_el_energy_both.C");
        //c_el_energy_both->SaveAs("c_el_energy_both.png");
        //c_el_energy_both->SaveAs("c_el_energy_both.pdf");
        //delete c_el_energy_both;
        const Double_t max_log_mode{1000.0e3};
        const Double_t max_nolog_mode{220.0e3};
        Double_t max{0.0};
        Double_t min{0.0};
        const std::string dir_log_mode("el_energy_log_dir");
        const std::string dir_nolog_mode("el_energy_nolog_dir");
        std::string dir;
        if(log_mode) { min = 0.1; max = max_log_mode; dir = dir_log_mode; }
        else { min = 0.0; max = max_nolog_mode; dir = dir_nolog_mode; }
        CanvasFactorySettings settings("Energy [MeV]", "Events", min, max, log_mode);
        settings.SetDrawOption("E");
        CanvasFactory factory(settings);
        factory.Canvas("el_energy", dir, h_el_energy_original, "Baseline", h_el_energy_reweight, "Reweighted");


        // print summed distribution
        //TCanvas *c_el_energy_sum_both = new TCanvas("e_el_energy_sum_both", "e_el_energy_sum_both", 800, 600);
        //c_el_energy_sum_both->SetLogy(log_mode);
        //h_el_energy_sum_original->SetMaximum(80.0e3); // MC data
        //h_el_energy_sum_original->SetMaximum(100.0e3); // MC data log mode
        //h_el_energy_sum_original->GetXaxis()->SetTitle("Energy [MeV]");
        //h_el_energy_sum_original->GetYaxis()->SetTitle("Events");
        //h_el_energy_sum_original->Draw("E");
        //h_el_energy_sum_reweight->Draw("Esame");
        //h_el_energy_sum_original_clone->Draw("Esame");
        //h_el_energy_sum_reweight_clone->Draw("Esame");
        //c_el_energy_sum_both->SaveAs("c_el_energy_sum_both.C");
        //c_el_energy_sum_both->SaveAs("c_el_energy_sum_both.png");
        //c_el_energy_sum_both->SaveAs("c_el_energy_sum_both.pdf");
        //delete c_el_energy_sum_both; 
        const Double_t max_log_mode_2{100.0e3};
        const Double_t max_nolog_mode_2{80.0e3};
        const std::string dir_log_mode_2("el_energy_sum_log_dir");
        const std::string dir_nolog_mode_2("el_energy_sum_nolog_dir");
        if(log_mode) { min = 0.1; max = max_log_mode_2; dir = dir_log_mode_2; }
        else { min = 0.0; max = max_nolog_mode_2; dir = dir_nolog_mode_2; }
        CanvasFactorySettings settings_2("Energy [MeV]", "Events", min, max, log_mode);
        settings_2.SetDrawOption("E");
        CanvasFactory factory_2(settings_2);
        factory_2.Canvas("el_energy_sum", dir, h_el_energy_sum_original, "Baseline", h_el_energy_sum_reweight, "Reweighted");
        // TODO: settings should be passed to canvas in case settings should be changed?
        // or setsettings function should be provided


        // print single electron distribution test histograms
        TCanvas *c_test_single = new TCanvas("c_test_single", "c_test_single", 800, 600);
        c_test_single->SetLogy(log_mode);
        //h_test_single_original->SetMaximum(1.4); // MC data
        //h_test_single_original->SetMaximum(1.2); // createtree data
        h_test_single_original->SetMaximum(10.0); // createtree data log mode
        h_test_single_original->GetXaxis()->SetTitle("Energy [MeV]");
        h_test_single_original->GetYaxis()->SetTitle("Events");
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
        // compute chi-square
        //std::cout << "chi1: " << chi_square_test(g_el_energy_single_0, h_test_single_original) / h_test_single_original->GetNbinsX() << std::endl;
        //std::cout << "chi2: " << chi_square_test(g_el_energy_single_2, h_test_single_reweight) / h_test_single_reweight->GetNbinsX() << std::endl;
        
        // print summed distribution test histograms
        TCanvas *c_test_sum = new TCanvas("c_test_sum", "c_test_sum", 800, 600);
        c_test_sum->SetLogy(log_mode);
        //h_test_sum_original->SetMaximum(1.0); // MC data
        //h_test_sum_original->SetMaximum(0.8); // createtree data
        h_test_sum_original->SetMaximum(1.0); // createtree data log mode
        h_test_sum_original->GetXaxis()->SetTitle("Energy [MeV]");
        h_test_sum_original->GetYaxis()->SetTitle("Events");
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
        // compute chi-square
        //std::cout << "chi3: " << chi_square_test(g_el_energy_sum_0, h_test_sum_original) / h_test_sum_original->GetNbinsX() << std::endl;
        //std::cout << "chi4: " << chi_square_test(g_el_energy_sum_2, h_test_sum_reweight) / h_test_sum_reweight->GetNbinsX() << std::endl;

    }

    // distribution of log-likelihoods
    std::vector<Double_t> vec_ll;
    TRandom3 gen;
    
    const Int_t number_of_pseudo_experiments{1000000};
    for(Int_t count{0}; count < number_of_pseudo_experiments; ++ count)
    {

        // log likelihood method
        // create "data" histogram - rounded original histogram
        std::string h_name{std::string("h_el_energy_data_") + std::to_string(count)};
        TH1I *h_el_energy_data = new TH1I(h_name.c_str(), "", num_bins, 0.0, 4.0);
        for(Int_t ix{1}; ix <= h_el_energy_sum_original->GetNbinsX(); ++ ix)
        {
            Double_t content{h_el_energy_sum_original->GetBinContent(ix)};
            //content = std::round(content);
            // Note: moved to pseudorandom data rather than rounding
            Int_t poisson_result{gen.Poisson(content)};
            //h_el_energy_data->SetBinContent(ix, content);
            h_el_energy_data->SetBinContent(ix, poisson_result);
        }

        // compute poisson likelihood for each bin
        Double_t likelihood{1.0};
        for(Int_t ix{1}; ix <= h_el_energy_sum_original->GetNbinsX(); ++ ix)
        {
            Double_t lambda{h_el_energy_sum_reweight->GetBinContent(ix)}; // reweighted
            Int_t data{std::lround(h_el_energy_data->GetBinContent(ix))}; // rounded baseline
            Double_t poi{TMath::Poisson(data, lambda)};
            //std::cout << "bin index = " << ix << ", poisson = " << poi << ", lambda = " << lambda << ", data = " << data << std::endl;
            likelihood *= poi;
        }
        //std::cout << "likelihood = " << likelihood << std::endl;
        Double_t log_likelihood{std::log(likelihood)};

        vec_ll.push_back(-2.0 * log_likelihood);
    
    }
    
    //Double_t min{std::min_element(vec_ll.begin(), vec_ll.end())};
    //Double_t max{std::max_element(vec_ll.begin(), vec_ll.end())};
    std::pair<std::vector<Double_t>::iterator, std::vector<Double_t>::iterator> min_max_pair{std::minmax_element(vec_ll.begin(), vec_ll.end())};
    Double_t min{*min_max_pair.first};
    Double_t max{*min_max_pair.second};
    TH1D *h_ll = new TH1D("h_ll", "", 100, min, max);
    h_ll->GetXaxis()->SetTitle("Log Likelihood Value");
    h_ll->GetYaxis()->SetTitle("Number of Pseudo Experiments");
    // fill the histogram
    for(std::vector<Double_t>::const_iterator it{vec_ll.cbegin()}; it != vec_ll.cend(); ++ it)
    {
        h_ll->Fill(*it);
    }

    /*
    TCanvas *c_ll = new TCanvas("c_ll", "", 800, 600);
    h_ll->Draw("E");
    c_ll->SaveAs("c_ll.png");
    delete c_ll;
    */


    // print entries
    //std::cout << "Number of entries in each histogram: h_el_energy_original: " << h_el_energy_original->GetEntries() << " h_el_energy_reweight: " << h_el_energy_reweight->GetEntries() << std::endl;


    // add data to data output file
    std::ofstream of_data(arg_output_filename.c_str(), std::ios::app);
    if(of_data.tellp() == 0)
    {
        of_data << "epsilon_31, chisquare (fit), degrees of freedom, chisquare (fit reduced), chisquare (sensitivity), chisquare (sensitivity reduced), -2log(l)" << std::endl;
    }
    of_data << epsilon_31 << ','
            << f_el_energy_sum_original->GetChisquare() << ','
            << non_empty_bins << ','
            << f_el_energy_sum_original->GetChisquare() / (Double_t)non_empty_bins << ','
            << sensitivity_chisquare << ','
            << sensitivity_chisquare / (Double_t)non_empty_bins; //<< ','
            //<< -2.0 * log_likelihood
            // TODO what goes here
    for(std::vector<Double_t>::const_iterator it{vec_ll.cbegin()}; it != vec_ll.cend(); ++ it)
    {
        of_data << ',' << *it;
    }
    of_data << std::endl;



    return 0;

}
