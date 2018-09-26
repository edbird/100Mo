// self header
#include "SubAnalysis.hpp"

#include "ReWeight.hpp"


////////////////////////////////////////////////////////////////////////////////
// CREATE HISTOGRAMS REQUIRED FOR EVENT LOOP
////////////////////////////////////////////////////////////////////////////////

void SubAnalysis::InitEventLoopHistogram()
{

    
    delete h_el_energy_sum_original;
    delete h_el_energy_sum_reweight;
    delete h_el_energy_original;
    delete h_el_energy_reweight;
    delete h_el_energy_diff;
    delete h_el_energy_pull;
    delete h_el_energy_ratio;
    delete h_el_energy_2d_original;
    delete h_el_energy_2d_reweight;
    delete h_el_energy_2d_diff;
    delete h_el_energy_2d_pull;
    delete h_el_energy_2d_ratio;
    
    h_el_energy_sum_original = nullptr;
    h_el_energy_sum_reweight = nullptr;
    h_el_energy_original = nullptr;
    h_el_energy_reweight = nullptr;
    h_el_energy_diff = nullptr;
    h_el_energy_pull = nullptr;
    h_el_energy_ratio = nullptr;
    h_el_energy_2d_original = nullptr;
    h_el_energy_2d_reweight = nullptr;
    h_el_energy_2d_diff = nullptr;
    h_el_energy_2d_pull = nullptr;
    h_el_energy_2d_ratio = nullptr;


    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS FOR SINGLE AND SUM ENERGY HISTOGRAMS FOR BASELINE (EPS=BASE)
    // AND RESCALED HISTOGRAMS (EPS=SOME OTHER VALUE)
    ////////////////////////////////////////////////////////////////////////////
    
    // load from NEMO-3 data (MC), the reconstructed electron energy (2 in same histogram)
    h_el_energy_original = new TH1D((std::string("h_el_energy_original") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
    //h_el_energy_original->SetStats(0);
    //h_el_energy_original->SetLineColor(2);
    //h_el_energy_original->SetMarkerColor(2);
    // same as above but re-weighted
    h_el_energy_reweight = new TH1D((std::string("h_el_energy_reweight") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
    //h_el_energy_reweight->SetStats(0);
    //h_el_energy_reweight->SetLineColor(4);
    //h_el_energy_reweight->SetMarkerColor(4);

    // load from NEMO-3 data (mc), the sum of the two reconstructed electron energies
    h_el_energy_sum_original = new TH1D((std::string("h_el_energy_sum_original") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
    h_el_energy_sum_original->SetStats(0);
    h_el_energy_sum_original->SetLineColor(2);
    h_el_energy_sum_original->SetMarkerColor(2);
    // same as above but re-weighted
    h_el_energy_sum_reweight = new TH1D((std::string("h_el_energy_sum_reweight") + h_name_append).c_str(), "", num_bins, 0.0, 4.0);
    h_el_energy_sum_reweight->SetStats(0);
    h_el_energy_sum_reweight->SetLineColor(4);
    h_el_energy_sum_reweight->SetMarkerColor(4); // TODO: don't think specific colors are required anymore

    // 2d fit version of h_el_energy_sum_original
    // instead of using a 1d histogram containing the sum of the 2 electron
    // energies, use a 2d histogram containing both energies
    // TODO: should I sort by energy?
    h_el_energy_2d_original = new TH2D((std::string("h_el_energy_2d_original") + h_name_append).c_str(), "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_el_energy_2d_original->SetStats(0);
    
    h_el_energy_2d_reweight = new TH2D((std::string("h_el_energy_2d_reweight") + h_name_append).c_str(), "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_el_energy_2d_reweight->SetStats(0);

}


////////////////////////////////////////////////////////////////////////////////
// EVENT LOOP
////////////////////////////////////////////////////////////////////////////////

void SubAnalysis::EventLoop()
{

    Double_t el_energy_0{el_energy_[0]};
    Double_t el_energy_1{el_energy_[1]}; // * systematic_energy_mult


    // TODO: not required
    if(*nElectrons != 2) return;

    
    // true energy
    // TODO: presumably this does not exist for data so need to search for
    // all instances of the trueT1, trueT2 variable and remove/replace
    Double_t T1{*trueT1 / bb_Q};
    Double_t T2{*trueT2 / bb_Q};

    // if MC apply energy degradation (correction)
    if(m_mode == MODE_FLAG::MODE_MC)
    {

        Double_t visible_true_ratio_1{1.0};
        Double_t visible_true_ratio_2{1.0};
        
        if(b_energy_correction_systematic_enabled == true)
        {
            //std::cout << "energy correction systematic is enabled" << std::endl;

            // if energy correction systematic is enabled, choose
            // energy correction value depending on which subanalysis
            // class this is
            // 0 is default
            // 1 is high systematic
            // -1 is low systematic
            if(m_energy_correction_systematic_mode == 0)
            {
                visible_true_ratio_1 = g_energy_correction->Eval(1000.0 * T1);
                visible_true_ratio_2 = g_energy_correction->Eval(1000.0 * T2);
            }
            else if(m_energy_correction_systematic_mode == 1)
            {
                visible_true_ratio_1 = g_energy_correction_systematic_high->Eval(1000.0 * T1);
                visible_true_ratio_2 = g_energy_correction_systematic_high->Eval(1000.0 * T2);
            }
            else if(m_energy_correction_systematic_mode == -1)
            {
                visible_true_ratio_1 = g_energy_correction_systematic_low->Eval(1000.0 * T1);
                visible_true_ratio_2 = g_energy_correction_systematic_low->Eval(1000.0 * T2);
            }
        }
        else
        {
            //std::cout << "energy correction systematic is disabled" << std::endl;

            // if systematics for energy correction are disabled...
            visible_true_ratio_1 = g_energy_correction->Eval(1000.0 * T1);
            visible_true_ratio_2 = g_energy_correction->Eval(1000.0 * T2);
        }

        //std::cout << "visible_true_ratio = " << visible_true_ratio_1 << ", " << visible_true_ratio_2 << std::endl;

        // apply energy correction with systematics if enabled
        el_energy_0 = el_energy_0 * visible_true_ratio_1;
        el_energy_1 = el_energy_0 * visible_true_ratio_2;
    }

    
    /*** SYSTEMATICS **********************************************************/

    // this if statement sorts out the logical problem of having different
    // high/low sysematic energy multipliers for the purpose of using them
    // as labels to address the SubAnalysis entries in the map inside Analysis,
    // and simultaniously allowing the systematic energy mult systematic to be
    // turned off while another systematic is on
    if(systematic_energy_mult_enable == true)
    {
        el_energy_0 = el_energy_0 * systematic_energy_mult;
        el_energy_1 = el_energy_1 * systematic_energy_mult;
    }
        
    // linear energy offset systematic
    el_energy_0 = el_energy_0 + systematic_energy_offset;
    el_energy_1 = el_energy_1 + systematic_energy_offset;

    // efficiency systematic
    // TODO: can remove, as weight_efficiency = systematic_efficiency
    Double_t weight_efficiency = 1.0;
    weight_efficiency = weight_efficiency * systematic_efficiency;

    // generator weight (MC weight) multiplied by weight efficiency
    Double_t aux_weight{*gen_weight};
    aux_weight = aux_weight * weight_efficiency;
    
    // TODO; what happens if the energy shift / systematics move the energy
    // out of a valid range
    // answer: nothing, reweight function depends only on T1 T2
    // TODO: should T1 and T2 be shifted by systematic?


    // TODO: energy degratation systematic

    // cut both electrons > 300 keV
    if(_energy_cut_enable_)
    {
        if(el_energy_0 < 0.3) return;
        if(el_energy_1 < 0.3) return;
    }


    // NOTE: more logical to set variables
    // weight_1, weight_2
    // for baseline and reweighted (now "baseline" and "test" / "universe")
    // then fill each histogram with each different weight
    // ?
    // NOTE: why would we reweight at all, why not use the decay rates from the
    // theorists directly?

    // ReWeight = baseline 0.0, ReWeight2 = baseline = 0.382
    Double_t weight{ReWeight2(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};
    //Double_t weight{ReWeight(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};

    // original
    h_el_energy_original->Fill(el_energy_0, 1.0 * aux_weight);
    h_el_energy_original->Fill(el_energy_1, 1.0 * aux_weight);
    
    h_el_energy_sum_original->Fill(el_energy_0 + el_energy_1, 1.0 * aux_weight);
    
    if(el_energy_0 <= el_energy_1)
    {
        h_el_energy_2d_original->Fill(el_energy_0, el_energy_1, 1.0 * aux_weight);
    }
    else
    {
        h_el_energy_2d_original->Fill(el_energy_1, el_energy_0, 1.0 * aux_weight);
    }

    // reweight
    h_el_energy_reweight->Fill(el_energy_0, weight * aux_weight);
    h_el_energy_reweight->Fill(el_energy_1, weight * aux_weight);

    h_el_energy_sum_reweight->Fill(el_energy_0 + el_energy_1, weight * aux_weight);

    if(el_energy_0 <= el_energy_1)
    {
        h_el_energy_2d_reweight->Fill(el_energy_0, el_energy_1, weight * aux_weight);
    }
    else
    {
        h_el_energy_2d_reweight->Fill(el_energy_1, el_energy_0, weight * aux_weight);
    }

}


////////////////////////////////////////////////////////////////////////////////
// POST PROCESSING
////////////////////////////////////////////////////////////////////////////////

void SubAnalysis::PostProcess()
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


}
