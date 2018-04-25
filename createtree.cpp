

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"



#include <iostream>
#include <fstream>


#include "ReWeight.hpp"
#include "read_data.hpp"


int main()
{

    // TODO: remove this section, need to figure out why 
    // if statements with continue are not working
    
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


    // section
    
    //TFile *f = new TFile("NewElectronNtuplizerExe_Int_ManDB_output.root");
    TFile *f = new TFile("f_histogram.root");
    
    TH2D *h_data_0 = (TH2D*)f->Get("h_data_0");

    // T1 reweight single and 2 electron
    // 
    // use c_data_0, to create new output->input tree, feed into this code
    // and check single/sum histograms to see if they are the same as 
    // the ones from the paper
    // TODO: interpolation

    TFile *f_out = new TFile("NewElectronNtuplizerExe_Int_ManDB_output.output.root", "recreate");
    TDirectory *t_dir = f_out->mkdir("NewElectronNtuplizer");
    f_out->cd("NewElectronNtuplizer");
    TTree *t_out = new TTree("NewElectronNtuplizer", "NewElectronNtuplizer");

    Double_t output_T1;
    Double_t output_T2;
    Double_t output_el_energy_[2]{0.0, 0.0};
    Int_t output_nElectrons{2};
    Double_t gen_weight{1.0};

    t_out->Branch("nElectrons", &output_nElectrons);
    t_out->Branch("trueT1", &output_T1);
    t_out->Branch("trueT2", &output_T2);
    t_out->Branch("el_energy_", output_el_energy_);
    t_out->Branch("gen_weight", &gen_weight);

    std::cout << "Processing..." << std::endl;
    const Double_t bb_Q{3.034};
    for(Int_t j{1}; j <= h_data_0->GetNbinsY(); ++ j)
    {
        //std::cout << "processing: j=" << j << std::endl;

        for(Int_t i{1}; i <= h_data_0->GetNbinsX(); ++ i)
        {

            // the input weight
            Double_t weight{h_data_0->GetBinContent(i, j)};

            // the "T1"
            //Double_t input_T1{h_data_0->GetXaxis()->GetBinLowEdge(i)};
            Double_t input_T1{h_data_0->GetXaxis()->GetBinCenter(i)};

            // the "T2"
            //Double_t input_T2{h_data_0->GetYaxis()->GetBinLowEdge(j)};
            Double_t input_T2{h_data_0->GetYaxis()->GetBinCenter(j)};
            
            // included to prevent tree filling with impossible values
            // which cannot be processed by the ReWeight function
            //if(input_T1 + input_T2 > 1.0)
            //const Int_t dimension_xy{1001};
            //if(i - 1 + j - 1 < dimension_xy - 1)
            //{
            //    continue;
            //}
            // TODO: this didn't work, why?

            
            // TODO: this also didn't work
            //if(input_T1 + input_T2 > 1.0)
            if(std::isnan(ReWeight(input_T1, input_T2, 0.8, h_nEqNull, h_nEqTwo, psiN0, psiN2, "false")))
            {
                continue;   
            }
            

            Double_t epsilon_31{0.8}; // TODO: this value follows previous value
            // TODO: need to change when previous value is changed (MANUALLY!)
            // TODO: use same variable

            // the output weight
            //Double_t weight{ReWeight(input_T1, input_T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "false")};
            // set output variables
            output_T1 = input_T1 * bb_Q;
            output_T2 = input_T2 * bb_Q;

            // TODO: how should the weight information be put into the output file?
            gen_weight = weight;
            //Double_t weight_scale{1.0e7};
            //Double_t weight_scale{1.0e3};

            // fill tree with weight
            //Long64_t fill_count{(Long64_t)(weight * weight_scale)};
            //for(Long64_t count{0}; count < fill_count; ++ count)
            //{
                t_out->Fill();
            //}


        }
    }

    t_out->Write();
    f_out->Close();

    std::cout << "Done" << std::endl;

    return 0;

}
