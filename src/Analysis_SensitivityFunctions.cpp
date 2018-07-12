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
