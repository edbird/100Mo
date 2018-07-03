// self header
#include "Analysis.hpp"


// local headers
#include "ReWeight.hpp"



////////////////////////////////////////////////////////////////////////////////
// CREATE HISTOGRAMS REQUIRED FOR EVENT LOOP
////////////////////////////////////////////////////////////////////////////////


void Analysis::InitEventLoopTree()
{

    ////////////////////////////////////////////////////////////////////////////
    // initialize file and branch addresses
    ////////////////////////////////////////////////////////////////////////////

    // read file using program arguments
    f = new TFile(filename.c_str());
    t = (TTree*)f->Get("NewElectronNtuplizer/NewElectronNtuplizer");

    t->SetBranchAddress("nElectrons", &nElectrons);
    t->SetBranchAddress("trueT1", &trueT1);
    t->SetBranchAddress("trueT2", &trueT2);
    t->SetBranchAddress("el_energy_", el_energy_);
    // new method requires weight to be saved to tree
    gen_weight = 1.0;
    if(_gen_weight_enable_ == true)
    {
        t->SetBranchAddress("gen_weight", &gen_weight);
    }

}


void Analysis::InitEventLoopHistogram()
{

    
    delete h_el_energy_original;
    delete h_el_energy_reweight;
    delete h_el_energy_sum_original;
    delete h_el_energy_sum_reweight;
    delete h_el_energy_2d_original;
    delete h_el_energy_2d_reweight;
    delete h_el_energy_2d_diff;
    delete h_el_energy_2d_pull;
    delete h_el_energy_2d_diff_data_rw;
    delete h_el_energy_2d_diff_data_orig;
    delete h_gen_weight;


    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS FOR SINGLE AND SUM ENERGY HISTOGRAMS FOR BASELINE (EPS=BASE)
    // AND RESCALED HISTOGRAMS (EPS=SOME OTHER VALUE)
    ////////////////////////////////////////////////////////////////////////////
    
    // load from NEMO-3 data (MC), the reconstructed electron energy (2 in same histogram)
    h_el_energy_original = new TH1D("h_el_energy_original", "", num_bins, 0.0, 4.0);
    //h_el_energy_original->SetStats(0);
    //h_el_energy_original->SetLineColor(2);
    //h_el_energy_original->SetMarkerColor(2);
    // same as above but re-weighted
    h_el_energy_reweight = new TH1D("h_el_energy_reweight", "", num_bins, 0.0, 4.0);
    //h_el_energy_reweight->SetStats(0);
    //h_el_energy_reweight->SetLineColor(4);
    //h_el_energy_reweight->SetMarkerColor(4);

    // load from NEMO-3 data (mc), the sum of the two reconstructed electron energies
    h_el_energy_sum_original = new TH1D("h_el_energy_sum_original", "", num_bins, 0.0, 4.0);
    h_el_energy_sum_original->SetStats(0);
    h_el_energy_sum_original->SetLineColor(2);
    h_el_energy_sum_original->SetMarkerColor(2);
    // same as above but re-weighted
    h_el_energy_sum_reweight = new TH1D("h_el_energy_sum_reweight", "", num_bins, 0.0, 4.0);
    h_el_energy_sum_reweight->SetStats(0);
    h_el_energy_sum_reweight->SetLineColor(4);
    h_el_energy_sum_reweight->SetMarkerColor(4); // TODO: don't think specific colors are required anymore

    // 2d fit version of h_el_energy_sum_original
    // instead of using a 1d histogram containing the sum of the 2 electron
    // energies, use a 2d histogram containing both energies
    // TODO: should I sort by energy?
    h_el_energy_2d_original = new TH2D("h_el_energy_2d_original", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_el_energy_2d_original->SetStats(0);
    
    h_el_energy_2d_reweight = new TH2D("h_el_energy_2d_reweight", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_el_energy_2d_reweight->SetStats(0);



    // NEMO-3 data/MC read from tree
    // input file

    h_gen_weight = new TH2D("h_gen_weight", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_gen_weight->SetStats(0);

}


////////////////////////////////////////////////////////////////////////////////
// EVENT LOOP
////////////////////////////////////////////////////////////////////////////////

void Analysis::EventLoop()
{

    std::cout << "Processing data" << std::endl;
    Long64_t prog_c{-1};
    //const Double epsilon_31{0.0};
    for(Long64_t ix{0}; ix < t->GetEntries(); ++ ix)
    {

        t->GetEntry(ix);

        if(nElectrons != 2) continue;

        // cut both electrons > 300 keV
        if(_energy_cut_enable_)
        {
            if(el_energy_[0] < 0.3) continue;
            if(el_energy_[1] < 0.3) continue;
        }

        Double_t T1{trueT1 / bb_Q};
        Double_t T2{trueT2 / bb_Q};

        // ReWeight = baseline 0.0, ReWeight2 = baseline = 0.382
        Double_t weight{ReWeight2(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};
        //Double_t weight{ReWeight(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};

        h_el_energy_original->Fill(el_energy_[0], 1.0 * gen_weight);
        h_el_energy_original->Fill(el_energy_[1], 1.0 * gen_weight);
        
        h_el_energy_sum_original->Fill(el_energy_[0] + el_energy_[1], 1.0 * gen_weight);
        
        if(el_energy_[0] <= el_energy_[1])
        {
            h_el_energy_2d_original->Fill(el_energy_[0], el_energy_[1], 1.0 * gen_weight);
        }
        else
        {
            h_el_energy_2d_original->Fill(el_energy_[1], el_energy_[0], 1.0 * gen_weight);
        }
        
        // test histograms
        if(h_test_single_original != nullptr)
        {
            h_test_single_original->Fill(trueT1, 1.0 * gen_weight);
            h_test_single_original->Fill(trueT2, 1.0 * gen_weight);
        }
        if(h_test_sum_original != nullptr)
        {
            h_test_sum_original->Fill(trueT1 + trueT2, 1.0 * gen_weight);
        }

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

        if(el_energy_[0] <= el_energy_[1])
        {
            h_el_energy_2d_reweight->Fill(el_energy_[0], el_energy_[1], weight * gen_weight);
        }
        else
        {
            h_el_energy_2d_reweight->Fill(el_energy_[1], el_energy_[0], weight * gen_weight);
        }

        // test histograms
        if(h_test_single_reweight != nullptr)
        {
            h_test_single_reweight->Fill(trueT1, weight * gen_weight);
            h_test_single_reweight->Fill(trueT2, weight * gen_weight);
        }
        if(h_test_sum_reweight != nullptr)
        {
            h_test_sum_reweight->Fill(trueT1 + trueT2, weight * gen_weight);
        }

        h_gen_weight->Fill(trueT1, trueT2, gen_weight);

        //if(ix % 10 == 9) std::cin.get();

        if((ix * 10) / t->GetEntries() > prog_c)
        {
            prog_c = (ix * 10) / t->GetEntries();
            //std::cout << 10 * (ix * 10) / t->GetEntries() << " %" << std::endl;
        }

    }

}


////////////////////////////////////////////////////////////////////////////////
// POST PROCESSING
////////////////////////////////////////////////////////////////////////////////

void Analysis::PostProcess()
{

    
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
    
    if(h_test_single_original != nullptr)
    {
        h_test_single_original->Scale((1.0 / 0.1) / h_test_single_original->Integral());
    }
    if(h_test_single_reweight != nullptr)
    {
        h_test_single_reweight->Scale((1.0 / 0.1) / h_test_single_reweight->Integral());
    }
    if(h_test_sum_original != nullptr)
    {
        h_test_sum_original->Scale((1.0 / 0.1) / h_test_sum_original->Integral());
    }
    if(h_test_sum_reweight != nullptr)
    {
        h_test_sum_reweight->Scale((1.0 / 0.1) / h_test_sum_reweight->Integral());
    }

    std::cout << "h_el_energy_original_integral=" << h_el_energy_original->Integral() << std::endl;

    // calculate integral manually
    /*
    Double_t h_el_energy_original_integral{0.0};
    for(Int_t i{1}; i <= h_el_energy_original->GetNbinsX(); ++ i)
    {
        h_el_energy_original_integral += h_el_energy_original->GetBinContent(i);
    }
    std::cout << "h_el_energy_original_integral=" << h_el_energy_original_integral << std::endl;
    */

    if(_batch_mode_ == false)
    {
        c_gen_weight = new TCanvas("c_gen_weight", "", 800, 600);
        c_gen_weight->SetRightMargin(0.15);
        h_gen_weight->Draw("colz");
        c_gen_weight->SaveAs("c_gen_weight.C");
        c_gen_weight->SaveAs("c_gen_weight.png");
        c_gen_weight->SaveAs("c_gen_weight.pdf");
        delete c_gen_weight;
    }


}
