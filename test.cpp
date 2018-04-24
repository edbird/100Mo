

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"



#include <iostream>
#include <fstream>

int main()
{
    
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

    for(Int_t j{1}; j <= h_data_0->GetNbinsY(); ++ j)
    {
        std::cout << "processing: j=" << j << std::endl;

        for(Int_t i{1}; i <= h_data_0->GetNbinsX(); ++ i)
        {
            // the input weight
            Double_t content{h_data_0->GetBinContent(i, j)};

            // the "T1"
            Double_t input_T1{h_data_0->GetXaxis()->GetBinLowEdge(i)};

            // the "T2"
            Double_t input_T2{h_data_0->GetYaxis()->GetBinLowEdge(j)};

            Double_t epsilon_31{0.8}; // TODO: this value follows previous value
            // TODO: need to change when previous value is changed (MANUALLY!)
            // TODO: use same variable

            // the output weight
            //Double_t weight{ReWeight(input_T1, input_T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "false")};
            Double_t weight{content};

            // set output variables
            output_T1 = input_T1;
            output_T2 = input_T2;

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


    return 0;

}
