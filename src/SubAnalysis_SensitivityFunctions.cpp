// self header
#include "SubAnalysis.hpp"


// local headers
#include "aux.hpp"


// C++ headers
#include <iostream>
#include <fstream>



////////////////////////////////////////////////////////////////////////////////
// SCALE SUMMED HISTOGRAM TO MATCH AMPLITUDE
////////////////////////////////////////////////////////////////////////////////

void SubAnalysis::SummedEnergyFit()
{

    // reset
    non_empty_bins_fit = 0;

    ////////////////////////////////////////////////////////////////////////
    // 1d histogram fit and sensitivity
    ////////////////////////////////////////////////////////////////////////

    delete f_el_energy_sum_original;
    f_el_energy_sum_original = nullptr;

    f_el_energy_sum_original = new TF1((std::string("f_el_energy_sum_original") + h_name_append).c_str(), fit_function, 0.0, 4.0, 1 + 2 * (h_el_energy_sum_reweight->GetNbinsX() + 1));
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
            h_el_energy_sum_original->Fit((std::string("f_el_energy_sum_original") + h_name_append).c_str(), "0Q");
        }
        else if(fit_subrange == true)
        {
            h_el_energy_sum_original->Fit((std::string("f_el_energy_sum_original") + h_name_append).c_str(), "0Q", "", 2.0, 4.0);
        }

    }
    std::cout << "Amplitude parameter: " << f_el_energy_sum_original->GetParameter(0) << std::endl;
    std::cout << "                err: " << f_el_energy_sum_original->GetParError(0) << std::endl;
    //std::cout << "         chi square: " << f_el_energy_sum_original->GetChisquare() / h_el_energy_sum_original->GetNbinsX() << std::endl;
    std::cout << "         chi square: " << f_el_energy_sum_original->GetChisquare() << std::endl;
    for(Int_t i{1}; i <= h_el_energy_sum_original->GetNbinsX(); ++ i)
    {
        if(h_el_energy_sum_original->GetBinContent(i) != 0.0) ++ non_empty_bins_fit;
    }
    std::cout << " degrees of freedom: " << non_empty_bins_fit << std::endl;
    std::cout << " chi square reduced: " << f_el_energy_sum_original->GetChisquare() / (Double_t)non_empty_bins_fit << std::endl;

    // scale the red histogram using the amplitude parameter
    h_el_energy_sum_reweight->Scale(f_el_energy_sum_original->GetParameter(0));
    h_el_energy_reweight->Scale(f_el_energy_sum_original->GetParameter(0));
    h_el_energy_2d_reweight->Scale(f_el_energy_sum_original->GetParameter(0));


    ////////////////////////////////////////////////////////////////////////////
    // FIT CANVAS OUTPUT SUMMED ENERGY
    ////////////////////////////////////////////////////////////////////////////

    #if 0
    // print summed distribution
    Double_t max{0.0};
    Double_t min{0.0};
    const Double_t max_log_mode_2{100.0e3};
    const Double_t max_nolog_mode_2{80.0e3};
    const std::string dir_log_mode_2("el_energy_sum_log_dir");
    const std::string dir_nolog_mode_2("el_energy_sum_nolog_dir");
    std::string dir; //{std::to_string(epsilon_31)};
    if(log_mode) { min = 0.1; max = max_log_mode_2; dir = dir_log_mode_2; }
    else { min = 0.0; max = max_nolog_mode_2; dir = dir_nolog_mode_2; }
    CanvasFactorySettings settings_2("Sum Electron Energy [MeV]", "Events", min, max, log_mode);
    settings_2.SetDrawOption("E");
    CanvasFactory factory_2(settings_2);
    factory_2.Canvas("el_energy_sum", dir, h_el_energy_sum_original, "Baseline", h_el_energy_sum_reweight, "Reweighted");
    // TODO: settings should be passed to canvas in case settings should be changed?
    // or setsettings function should be provided
    #endif

}


////////////////////////////////////////////////////////////////////////////////
// CHISQUARE FIT BETWEEN REWEIGHTED AND BASELINE (NO SYSTEMATIC)
////////////////////////////////////////////////////////////////////////////////

void SubAnalysis::ChiSquare_BaselineNoSystematic()
{

    if(_baseline_histo_p_ == nullptr)
    {
        std::cout << "Error: Unable to call SubAnalysis::ChiSquare_BaselineNoSystematic(): _baseline_histo_p_ = nullptr" << std::endl;
    }
    else
    {
        // call chisquare function

        // reset
        non_empty_bins_baseline_nosystematic = 0;

        ////////////////////////////////////////////////////////////////////////////
        // SINGLE ELECTRON ENERGY CHISQUARE METHOD
        ////////////////////////////////////////////////////////////////////////////

        // get chi-square for single electron histogram
        if(fit_subrange == false)
        {
            
            for(Int_t i{1}; i <= _baseline_histo_p_->GetNbinsX(); ++ i)
            {
                if(_baseline_histo_p_->GetBinContent(i) != 0.0) ++ non_empty_bins_baseline_nosystematic;
            }
            
            sensitivity_chisquare_baseline_nosystematic = chi_square_test(h_el_energy_reweight, _baseline_histo_p_);
            //std::cout << "chi square of single electron: " << sensitivity_chisquare_baseline_nosystematic << std::endl;
            //std::cout << " degrees of freedom: " << non_empty_bins_baseline_nosystematic << std::endl;
            //std::cout << " chi square reduced: " << sensitivity_chisquare_baseline_nosystematic / (Double_t)non_empty_bins_baseline_nosystematic << std::endl;

        }
        else if(fit_subrange == true)
        {
        
            for(Int_t i{1}; i <= _baseline_histo_p_->GetNbinsX(); ++ i)
            {
                if(_baseline_histo_p_->GetBinCenter(i) >= 2.0)
                {
                    if(_baseline_histo_p_->GetBinContent(i) != 0.0) ++ non_empty_bins_baseline_nosystematic;
                }
            }
            
            sensitivity_chisquare_baseline_nosystematic = chi_square_test(h_el_energy_reweight, _baseline_histo_p_, 2.0, 4.0);
            //std::cout << "chi square of single electron, 2.0 MeV - 4.0 MeV: " << sensitivity_chisquare_baseline_nosystematic << std::endl;
            //std::cout << " degrees of freedom: " << non_empty_bins_baseline_nosystematic << std::endl;
            //std::cout << " chi square reduced: " << sensitivity_chisquare_baseline_nosystematic / (Double_t)non_empty_bins_baseline_nosystematic << std::endl;
        }
    }
}


////////////////////////////////////////////////////////////////////////////////
// PRINT OUTPUT TO FILE
////////////////////////////////////////////////////////////////////////////////

void SubAnalysis::PrintOutputToFile()
{

    // the chisquare data for production of chisquare curve graphs is sent to
    // an intermediate output file here
    // the chisquare code needs to be modified to read in data from each of the
    // subanalysis files (default, high, low) and produce an output file with
    // 3 curves - for each of the subanalysis systematic files

    // ******************************
    // variable sensitivity_chisquare
    // ******************************
    //
    // chisquare value for comparison between reweighted and baseline
    // histograms
    //
    // both the baseline and reweighted histograms for this subclass are
    // constructed with the systematic energy shift for this subclass
    //
    // what we want is to measure the chisquare between the baseline
    // histogram with no systematic and the reweighted histogram with a
    // systematic, this reweighted histogram is therefore dependent on both
    // the systematic parameter and the xi parameter
    //
    // this will be implemented in the Analysis_SensitivityFunctions.cpp file
    // in the function MakeSensitivityCanvas()
    //
    

    // add data to sensitivity record storage
    SensitivityRecord r;
    r.epsilon_31 = epsilon_31;
    r.systematic_energy_mult = systematic_energy_mult;
    if(systematic_energy_mult_enable == true) r.systematic_energy_mult_enable = 1; else r.systematic_energy_mult_enable = 0;
    r.systematic_energy_offset = systematic_energy_offset;
    r.fit_chisquare = f_el_energy_sum_original->GetChisquare();
    r.fit_chisquare_bin_count = non_empty_bins_fit;
    r.sensitivity_chisquare_1d = sensitivity_chisquare;
    r.sensitivity_chisquare_1d_bin_count = non_empty_bins;
    r.sensitivity_chisquare_2d = sensitivity_chisquare_2d;
    r.sensitivity_chisquare_2d_bin_count = non_empty_bins_2d;
    
    r.sensitivity_chisquare_1d_baseline_nosystematic = sensitivity_chisquare_baseline_nosystematic;
    r.sensitivity_chisquare_1d_bin_count_baseline_nosystematic = non_empty_bins_baseline_nosystematic;


    std::map<Double_t, std::map<Double_t, SensitivityRecord>> &srm{*_sensitivity_record_map_};
    //if(srm.count(epsilon_31) == 0)
    //{
        //if(srm.at(epsilon_31).count(systematic_energy_mult) == 0)
        if(srm[epsilon_31].count(systematic_energy_mult) == 0)
        {
            // TODO: does this work? what goes here?
           srm[epsilon_31][systematic_energy_mult] = r;
            //_sensitivity_record_map_.insert(epsilon_31, ???[systematic_energy_mult] = r;
        }
        else
        {
            std::cout << "Error: Failed to enter SensitivityRecord: systematic_energy_mult count > 0 (systematic_energy_mult=" << systematic_energy_mult << ")" << std::endl;
        }
    //}
    //else
    //{
    //    std::cout << "Error: Failed to enter SensitivityRecord: epsilon_31 count > 0, (epsilon_31=" << epsilon_31 << ")" << std::endl;
    //}
    

    // add data to data output file
    std::string filename_chisq(output_filename.substr(0, output_filename.find(".txt")) + std::string("_chisq") + h_name_append + std::string(".txt"));
    std::cout << filename_chisq << std::endl;
    std::ofstream of_data_chisq(filename_chisq.c_str(), std::ios::app);
    if(of_data_chisq.tellp() == 0)
    {
        of_data_chisq << "epsilon_31,fit (summed): chisq,reduced chisq,dof,sensitivity (1d single): chisq,reduced chisq,dof,sensitivity (2d single): chisq,reduced chisq,dof" << std::endl;
    }
    of_data_chisq << epsilon_31 << ','
                  << f_el_energy_sum_original->GetChisquare() << ','
                  << f_el_energy_sum_original->GetChisquare() / (Double_t)non_empty_bins_fit << ','
                  << non_empty_bins_fit << ','
                  << sensitivity_chisquare << ','
                  << sensitivity_chisquare / (Double_t)non_empty_bins << ','
                  << non_empty_bins << ','
                  << sensitivity_chisquare_2d << ','
                  << sensitivity_chisquare_2d / (Double_t)non_empty_bins_2d << ','
                  << non_empty_bins_2d << std::endl;


}
