// self header
#include "Analysis.hpp"


// local headers
#include "aux.hpp"


// C++ headers
#include <iostream>
#include <fstream>



////////////////////////////////////////////////////////////////////////////////
// SCALE SUMMED HISTOGRAM TO MATCH AMPLITUDE
////////////////////////////////////////////////////////////////////////////////

void Analysis::SummedEnergyFit()
{

    // reset
    non_empty_bins_fit = 0;

    ////////////////////////////////////////////////////////////////////////
    // 1d histogram fit and sensitivity
    ////////////////////////////////////////////////////////////////////////

    f_el_energy_sum_original = new TF1("f_el_energy_sum_original", fit_function, 0.0, 4.0, 1 + 2 * (h_el_energy_sum_reweight->GetNbinsX() + 1));
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

}


////////////////////////////////////////////////////////////////////////////////
// PRINT OUTPUT TO FILE
////////////////////////////////////////////////////////////////////////////////

void Analysis::PrintOutputToFile()
{


    // print entries
    //std::cout << "Number of entries in each histogram: h_el_energy_original: " << h_el_energy_original->GetEntries() << " h_el_energy_reweight: " << h_el_energy_reweight->GetEntries() << std::endl;


    {
        // add data to data output file
        std::ofstream of_data(output_filename.c_str(), std::ios::app);
        if(of_data.tellp() == 0)
        {
            of_data << "epsilon_31,chisquare (fit),degrees of freedom,chisquare (fit reduced),chisquare (sensitivity),chisquare (sensitivity reduced),-2log(l)" << std::endl;
        }
        of_data << epsilon_31 << ','
                << f_el_energy_sum_original->GetChisquare() << ','
                << non_empty_bins_fit << ','
                << f_el_energy_sum_original->GetChisquare() / (Double_t)non_empty_bins_fit << ','
                << sensitivity_chisquare << ','
                << sensitivity_chisquare / (Double_t)non_empty_bins; //<< ','
                //<< -2.0 * log_likelihood
                // TODO what goes here
        for(std::vector<Double_t>::const_iterator it{vec_ll.cbegin()}; it != vec_ll.cend(); ++ it)
        {
            of_data << ',' << *it;
        }
        of_data << std::endl;


        // new output format, do not change above - will break chi-square code
        std::string filename_chisq(output_filename.substr(0, output_filename.find(".txt")) + std::string("_chisq") + std::string(".txt"));
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

    {
        // 1d ll data
        std::string filename_ll(output_filename.substr(0, output_filename.find(".txt")) + std::string("_ll") + std::string(".txt"));
        std::cout << filename_ll << std::endl;
        std::ofstream of_data_ll(filename_ll.c_str(), std::ios::app);
        for(std::vector<Double_t>::const_iterator it{vec_ll.cbegin()}; it != vec_ll.cend(); ++ it)
        {
            of_data_ll << *it;
            if(it + 1 != vec_ll.cend()) of_data_ll << ',';
        }
        of_data_ll << std::endl;

    }

    {
        // 2d ll data
        std::string filename_ll_2d(output_filename.substr(0, output_filename.find(".txt")) + std::string("_ll_2d") + std::string(".txt"));
        std::cout << filename_ll_2d << std::endl;
        std::ofstream of_data_ll_2d(filename_ll_2d.c_str(), std::ios::app);
        for(std::vector<Double_t>::const_iterator it{vec_ll_2d.cbegin()}; it != vec_ll_2d.cend(); ++ it)
        {
            of_data_ll_2d << *it;
            if(it + 1 != vec_ll_2d.cend()) of_data_ll_2d << ',';
        }
        of_data_ll_2d << std::endl;
    }


}
