// self header
#include "Analysis.hpp"


// local headers
#include "ReWeight.hpp"


void Analysis::InitEventLoopTree()
{

    ////////////////////////////////////////////////////////////////////////////
    // initialize file and branch addresses
    ////////////////////////////////////////////////////////////////////////////

    // NEMO-3 data/MC read from tree
    // input file

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
    
    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->SetGenWeightEnabled(_gen_weight_enable_);
    }

}


////////////////////////////////////////////////////////////////////////////////
// CREATE HISTOGRAMS REQUIRED FOR EVENT LOOP
////////////////////////////////////////////////////////////////////////////////

void Analysis::InitEventLoopHistogram()
{
    
    delete h_gen_weight;
    
    delete _subanalysis_systematic_default_;
    delete _subanalysis_systematic_low_;
    delete _subanalysis_systematic_high_;
    
    _subanalysis_systematic_default_ = nullptr;
    _subanalysis_systematic_low_ = nullptr;
    _subanalysis_systematic_high_ = nullptr;


    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS FOR SINGLE AND SUM ENERGY HISTOGRAMS FOR BASELINE (EPS=BASE)
    // AND RESCALED HISTOGRAMS (EPS=SOME OTHER VALUE)
    ////////////////////////////////////////////////////////////////////////////

    h_gen_weight = new TH2D("h_gen_weight", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_gen_weight->SetStats(0);
 
 
 
    ////////////////////////////////////////////////////////////////////////////
    // INIT SUBCLASS
    ////////////////////////////////////////////////////////////////////////////
 

     _subanalysis_systematic_default_ = new SubAnalysis("_default",
                                                        output_filename,
                                                        &_sensitivity_record_map_,
                                                        epsilon_31,
                                                        systematic_energy_mult, systematic_energy_mult_enabled,
                                                        systematic_energy_offset,
                                                        systematic_efficiency,
                                                        h_nEqNull, h_nEqTwo, psiN0, psiN2,
                                                        &nElectrons, &trueT1, &trueT2, el_energy_, &gen_weight,
                                                        &gen);
                                                       
    _subanalysis_systematic_low_ = new SubAnalysis("_low",
                                                   output_filename,
                                                   &_sensitivity_record_map_,
                                                   epsilon_31,
                                                   systematic_energy_mult_low, systematic_energy_mult_enabled,
                                                   systematic_energy_offset_low,
                                                   systematic_efficiency_low,
                                                   h_nEqNull, h_nEqTwo, psiN0, psiN2,
                                                   &nElectrons, &trueT1, &trueT2, el_energy_, &gen_weight,
                                                   &gen);
                                                   
    _subanalysis_systematic_high_ = new SubAnalysis("_high",
                                                    output_filename,
                                                    &_sensitivity_record_map_,
                                                    epsilon_31,
                                                    systematic_energy_mult_high, systematic_energy_mult_enabled,
                                                    systematic_energy_offset_high,
                                                    systematic_efficiency_high,
                                                    h_nEqNull, h_nEqTwo, psiN0, psiN2,
                                                    &nElectrons, &trueT1, &trueT2, el_energy_, &gen_weight,
                                                    &gen);
 
    _subanalysis_.clear();
    _subanalysis_.push_back(_subanalysis_systematic_default_);
    _subanalysis_.push_back(_subanalysis_systematic_low_);
    _subanalysis_.push_back(_subanalysis_systematic_high_);
 
    // set mode for energy correction systematic 
    _subanalysis_systematic_default_->SetEnergyCorrectionSystematicMode(0);
    _subanalysis_systematic_low_->SetEnergyCorrectionSystematicMode(-1);
    _subanalysis_systematic_high_->SetEnergyCorrectionSystematicMode(1);

    // set mode
    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        // MC / DATA flag
        (*it)->SetMode(m_mode);

        // set g_energy_correction, regardless if used
        // data for energy degradation
        (*it)->Set_g_energy_correction(g_energy_correction);
        // TODO: move g_energy_correction_[high/low] (systematic graph data) to here
        // do not need copies of them for each sub analysis class

        // enable/disable energy correction systematic
        (*it)->SetEnergyCorrectionSystematicEnabled(b_energy_correction_systematic_enabled);
    }


    // set energy calibration
    // call: RunEnergyCalibration() before this function
    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        Double_t a{energy_calibration_a};
        Double_t b{energy_calibration_b};
        (*it)->SetEnergyCalibrationConstants(a, b);
    }


    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->InitEventLoopHistogram();
    }

    // at this point, subanalysis class InitEventLoopHistogram functions have
    // been called, therefore h_el_energy_original pointer has changed
    TH1D *h_temp{_subanalysis_systematic_default_->GetElEnergyOriginalHistoPointer()};
    _subanalysis_systematic_default_->SetBaselineHistoPointer(h_temp);
    _subanalysis_systematic_low_->SetBaselineHistoPointer(h_temp);
    _subanalysis_systematic_high_->SetBaselineHistoPointer(h_temp);

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
        
        for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
        {
            (*it)->EventLoop();
        }

        // note: no systematic energy shift in this class
        //Double_t el_energy_0{el_energy_[0] * systematic_energy_mult};
        //Double_t el_energy_1{el_energy_[1] * systematic_energy_mult};
        Double_t el_energy_0{el_energy_[0]};
        Double_t el_energy_1{el_energy_[1]};


        // analysis only valid for 2 electron events
        if(nElectrons != 2) continue;
        

        // true energy
        // TODO: presumably this does not exist for data so need to search for
        // all instances of the trueT1, trueT2 variable and remove/replace
        Double_t T1{trueT1 / bb_Q};
        Double_t T2{trueT2 / bb_Q};

        // if MC apply energy degradation (correction)
        if(m_mode == MODE_FLAG::MODE_MC)
        {
            Double_t visible_true_ratio_1{g_energy_correction->Eval(T1)};
            Double_t visible_true_ratio_2{g_energy_correction->Eval(T2)};

            el_energy_0 = el_energy_0 * visible_true_ratio_1;
            el_energy_1 = el_energy_0 * visible_true_ratio_2;
        }

        // TODO: not sure if this should be removed
        // cut both electrons > 300 keV
        if(_energy_cut_enable_)
        {
            if(el_energy_0 < 0.3) continue;
            if(el_energy_1 < 0.3) continue;
        }


        // ReWeight = baseline 0.0, ReWeight2 = baseline = 0.382
        Double_t weight{ReWeight2(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};
        //Double_t weight{ReWeight(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};
        
        // TODO: do the test histograms stay here? they are not dependent on systematic shift?
        // so should stay here?
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

        /*
        if((ix * 10) / t->GetEntries() > prog_c)
        {
            prog_c = (ix * 10) / t->GetEntries();
            //std::cout << 10 * (ix * 10) / t->GetEntries() << " %" << std::endl;
        }
        */

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

    //std::cout << "h_el_energy_original_integral=" << h_el_energy_original->Integral() << std::endl;

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

    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->PostProcess();
    }

}
