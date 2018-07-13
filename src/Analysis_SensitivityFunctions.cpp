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

    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->SummedEnergyFit();
    }

}


////////////////////////////////////////////////////////////////////////////////
// PRINT OUTPUT TO FILE
////////////////////////////////////////////////////////////////////////////////

void Analysis::PrintOutputToFile()
{

    // TODO: this no longer works because these variables are members of the subclass
    // however instead of printing for each systematic energy shift value, I may just want to
    // process the data for the systematic energy shift here and then print from here

    {
        
        // get refs for default analysis subclass
        Double_t                &epsilon_31                 {_subanalysis_systematic_default_->epsilon_31};
        TF1                     *f_el_energy_sum_original   {_subanalysis_systematic_default_->f_el_energy_sum_original};
        Int_t                   &non_empty_bins_fit         {_subanalysis_systematic_default_->non_empty_bins_fit};
        Double_t                &sensitivity_chisquare      {_subanalysis_systematic_default_->sensitivity_chisquare};
        Int_t                   &non_empty_bins             {_subanalysis_systematic_default_->non_empty_bins};
        Double_t                &sensitivity_chisquare_2d   {_subanalysis_systematic_default_->sensitivity_chisquare_2d};
        Int_t                   &non_empty_bins_2d          {_subanalysis_systematic_default_->non_empty_bins_2d};
        std::vector<Double_t>   &vec_ll                     {_subanalysis_systematic_default_->vec_ll};
    
    
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
    
        std::vector<Double_t> &vec_ll{_subanalysis_systematic_default_->vec_ll};
    
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
    
        std::vector<Double_t> &vec_ll_2d{_subanalysis_systematic_default_->vec_ll_2d};
    
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
    
    
    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->PrintOutputToFile();
    }


}


////////////////////////////////////////////////////////////////////////////////
// Produce canvas with ratio histograms
// 1: ratio: Reweight / Baseline
// 2: ratio: systematic high / Baseline
// 3: ratio: systematic low / Baseline
////////////////////////////////////////////////////////////////////////////////

void Analysis::MakeSensitivityCanvas()
{
    
    const std::string eps_string{std::to_string(epsilon_31)};

    TH1D *h_el_energy_ratio{_subanalysis_systematic_default_->h_el_energy_ratio};
    TH1D *h_el_energy_ratio_energy_systematic_low = nullptr; //new TH1D("h_el_energy_ratio_energy_systematic_low", "", num_bins, 0.0, 4.0);
    TH1D *h_el_energy_ratio_energy_systematic_high = nullptr; //new TH1D("h_el_energy_ratio_energy_systematic_high", "", num_bins, 0.0, 4.0);

    TH1D *h_el_energy_original = _subanalysis_systematic_default_->h_el_energy_original;
    TH1D *h_el_energy_original_energy_systematic_low = _subanalysis_systematic_low_->h_el_energy_original;
    TH1D *h_el_energy_original_energy_systematic_high = _subanalysis_systematic_high_->h_el_energy_original;

    // create ratio histogram, systematic low
    std::string h_el_energy_ratio_energy_systematic_low_name{std::string("h_el_energy_ratio_energy_systematic_low") + std::string("_") + eps_string};
    h_el_energy_ratio_energy_systematic_low = new TH1D(h_el_energy_ratio_energy_systematic_low_name.c_str(), "", num_bins, 0.0, 4.0);
    h_el_energy_ratio_energy_systematic_low->SetStats(0);
    for(Int_t i{1}; i <= h_el_energy_ratio_energy_systematic_low->GetNbinsX(); ++ i)
    {
        Double_t content1{h_el_energy_original->GetBinContent(i)};
        Double_t content2{h_el_energy_original_energy_systematic_low->GetBinContent(i)};
        Double_t error1{h_el_energy_original->GetBinError(i)}; // correct way round?
        // note: do not use error or reweighted
        Double_t error2{0.0 * h_el_energy_original_energy_systematic_low->GetBinError(i)};
        if(content1 != 0.0)
        {
            Double_t content{content2 / content1}; // correct way round?
            Double_t error{std::sqrt(std::pow(error1 * (1.0 / content2), 2.0) + std::pow(error2 * (content1 / (content2 * content2)), 2.0))}; // should error2 be used?
            h_el_energy_ratio_energy_systematic_low->SetBinContent(i, content);
            h_el_energy_ratio_energy_systematic_low->SetBinError(i, error);
            // what to do if content1 = 0.0 -> div by zero error
        }
    }
    
    // create ratio histogram, systematic high
    std::string h_el_energy_ratio_energy_systematic_high_name{std::string("h_el_energy_ratio_energy_systematic_high") + std::string("_") + eps_string};
    h_el_energy_ratio_energy_systematic_high = new TH1D(h_el_energy_ratio_energy_systematic_high_name.c_str(), "", num_bins, 0.0, 4.0);
    h_el_energy_ratio_energy_systematic_high->SetStats(0);
    for(Int_t i{1}; i <= h_el_energy_ratio_energy_systematic_high->GetNbinsX(); ++ i)
    {
        Double_t content1{h_el_energy_original->GetBinContent(i)};
        Double_t content2{h_el_energy_original_energy_systematic_high->GetBinContent(i)};
        Double_t error1{h_el_energy_original->GetBinError(i)}; // correct way round?
        // note: do not use error or reweighted
        Double_t error2{0.0 * h_el_energy_original_energy_systematic_high->GetBinError(i)};
        if(content1 != 0.0)
        {
            Double_t content{content2 / content1}; // correct way round?
            Double_t error{std::sqrt(std::pow(error1 * (1.0 / content2), 2.0) + std::pow(error2 * (content1 / (content2 * content2)), 2.0))}; // should error2 be used?
            h_el_energy_ratio_energy_systematic_high->SetBinContent(i, content);
            h_el_energy_ratio_energy_systematic_high->SetBinError(i, error);
            // what to do if content1 = 0.0 -> div by zero error
        }
    }

    h_el_energy_ratio_energy_systematic_low->SetMarkerColor(3);
    h_el_energy_ratio_energy_systematic_low->SetLineColor(3);
    h_el_energy_ratio_energy_systematic_high->SetMarkerColor(2);
    h_el_energy_ratio_energy_systematic_high->SetLineColor(2);

    h_el_energy_ratio->SetMaximum(2.0);
    h_el_energy_ratio->SetMinimum(0.0);

    TCanvas* c_ratio_sensitivity = new TCanvas("c_ratio_sensitivity", "", 800, 600);
    h_el_energy_ratio->Draw("E");
    h_el_energy_ratio_energy_systematic_low->Draw("Esame");
    h_el_energy_ratio_energy_systematic_high->Draw("Esame");
    std::string c_ratio_sensitivity_name(std::string("c_ratio_sensitivity_name") + std::string("_") + eps_string);
    c_ratio_sensitivity->SaveAs((c_ratio_sensitivity_name + std::string(".png")).c_str());
    c_ratio_sensitivity->SaveAs((c_ratio_sensitivity_name + std::string(".pdf")).c_str());
    c_ratio_sensitivity->SaveAs((c_ratio_sensitivity_name + std::string(".eps")).c_str());
    c_ratio_sensitivity->SaveAs((c_ratio_sensitivity_name + std::string(".C")).c_str());

    delete c_ratio_sensitivity;
    c_ratio_sensitivity = nullptr;


    delete h_el_energy_ratio_energy_systematic_low;
    delete h_el_energy_ratio_energy_systematic_high;
    h_el_energy_ratio_energy_systematic_low = nullptr;
    h_el_energy_ratio_energy_systematic_high = nullptr;

}
