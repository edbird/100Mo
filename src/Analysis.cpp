// self header
#include "Analysis.hpp"


// local headers
#include "ReWeight.hpp"
#include "aux.hpp"


Analysis::Analysis(const std::string& filename, const std::string& output_filename)
    : epsilon_31{0.368}
    , systematic_energy_mult{1.0}
    , systematic_energy_mult_high{1.0}
    , systematic_energy_mult_low{1.0}
    , f{nullptr}
    , t{nullptr}
    , num_bins{40}
    , filename{filename}
    , output_filename{output_filename}
    , _batch_mode_{false}
    , log_mode{false}
    , _gen_weight_enable_{false}
    , _energy_cut_enable_{false} // energy cut flag
    , _canvas_enable_raw_data_{false}
    , _canvas_enable_decay_rate_{true}
    , _canvas_enable_single_electron_projection_{false}
    , g_el_energy_single_0{nullptr}
    , g_el_energy_single_1{nullptr}
    , g_el_energy_single_2{nullptr}
    , g_el_energy_sum_0{nullptr}
    , g_el_energy_sum_1{nullptr}
    , g_el_energy_sum_2{nullptr}
    , h_nEqNull{nullptr}
    , h_nEqTwo{nullptr}
    , psiN0{0.0}
    , psiN2{0.0}
    , fit_subrange{false} // fit subrange flag
    , h_data_0{nullptr}
    , h_data_1{nullptr}
    , h_data_2{nullptr}
    , h_single_electron_0{nullptr}
    , h_single_electron_1{nullptr}
    , h_single_electron_2{nullptr}
    , h_test_single_original{nullptr}
    , h_test_single_reweight{nullptr}
    , h_test_sum_original{nullptr}
    , h_test_sum_reweight{nullptr}
    
    , h_gen_weight{nullptr}
    , c_nEqNull{nullptr}
    , c_nEqTwo{nullptr}
    , c_data_0{nullptr}
    , c_data_1{nullptr}
    , c_data_2{nullptr}
    , c_single_electron{nullptr}
    , c_test_single{nullptr}
    , c_test_sum{nullptr}
    , c_gen_weight{nullptr}
    , c_el_energy_diff{nullptr}
    , c_el_energy_pull{nullptr}
    , c_el_energy_2d_original{nullptr}
    , c_el_energy_2d_reweight{nullptr}
    , c_el_energy_2d_diff{nullptr}
    , c_el_energy_2d_pull{nullptr}
    , c_ll{nullptr}
    , c_el_energy_2d_diff_data_rw{nullptr}
    , c_el_energy_2d_diff_data_orig{nullptr}
    , c_el_energy_2d_data{nullptr}
    , c_el_energy_2d_prob{nullptr}
    , c_ll_2d{nullptr}
    
    , _subanalysis_systematic_default_{nullptr}
    , _subanalysis_systematic_low_{nullptr}
    , _subanalysis_systematic_high_{nullptr}
{

}


Analysis::~Analysis()
{

}


////////////////////////////////////////////////////////////////////////////////
// ANALYSIS FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

void Analysis::SetSystematicEnergyMultiplier(const Double_t value)
{
    systematic_energy_mult = value;
}

void Analysis::SetSystematicEnergyMultiplierHighLow(const Double_t high, const Double_t low)
{
    systematic_energy_mult_low = low;
    systematic_energy_mult_high = high;
}

void Analysis::SetEpsilon31(const Double_t epsilon)
{
    epsilon_31 = epsilon;
    
    /*
    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->SetEpsilon31(epsilon);
    }
    */
}

void Analysis::AddEpsilonValue(const Double_t epsilon)
{
    vec_epsilon_31.push_back(epsilon);
}

void Analysis::RunOverEpsilonVector()
{

    ////ReadData();
    ////CanvasDecayRate();
    ////CanvasSingleElectronProjection();
    //CanvasSingleElectronTest(); // must be done after Fill in event loop

    ////InitEventLoopTree();
    //InitEventLoopHistogram();

    std::vector<Double_t>::const_iterator it{vec_epsilon_31.cbegin()};
    for(; it != vec_epsilon_31.cend(); ++ it)
    {

        //epsilon_31 = *it;
        SetEpsilon31(*it); // set for all

        std::cout << "epsilon_31=" << epsilon_31 << std::endl;

        //std::cout << "Press Enter to Init" << std::endl;
        //std::cin.get();
        InitEventLoopHistogram();
        
        //std::cout << "Press Enter to Loop" << std::endl;
        //std::cin.get();
        EventLoop();

        //std::cout << "Press Enter to Postprocess" << std::endl;
        //std::cin.get();
        PostProcess();

        //std::cout << "Press Enter to Fit" << std::endl;
        //std::cin.get();
        SummedEnergyFit();
        
        //std::cout << "Press Enter to Chisquare1" << std::endl;
        //std::cin.get();
        SensitivityMeasurementChisquare1();
        
        //std::cout << "Press Enter to Chisquare2" << std::endl;
        //std::cin.get();
        SensitivityMeasurementChisquare2();
        
        //std::cout << "Press Enter to Loglikelihood1" << std::endl;
        //std::cin.get();
        //SensitivityMeasurementLoglikelihood1();
        
        //std::cout << "Press Enter to Loglikelihood2" << std::endl;
        //std::cin.get();
        //SensitivityMeasurementLoglikelihood2();

        //std::cout << "Press Enter to Print" << std::endl;
        //std::cin.get();
        PrintOutputToFile();

        MakeSensitivityCanvas();
    }


}


////////////////////////////////////////////////////////////////////////////////
// FLAG FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

void Analysis::SetBatchMode(const bool mode)
{
    _batch_mode_ = mode;
    
    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->SetBatchMode(mode);
    }
}

void Analysis::SetLogMode(const bool mode)
{
    log_mode = mode;
    
    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->SetLogMode(mode);
    }
}

void Analysis::SetEnergyCutEnabled(const bool enabled)
{
    _energy_cut_enable_ = enabled;
    
    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->SetEnergyCutEnabled(enabled);
    }
}

void Analysis::SetFitSubrange(const bool flag)
{
    fit_subrange = flag;

    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->SetFitSubrange(flag);
    }
}

void Analysis::SetGenWeightEnabled(const bool enabled)
{
    _gen_weight_enable_ = enabled;
    
    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->SetGenWeightEnabled(enabled);
    }
}

void Analysis::SetCanvasEnableRawData(const bool mode)
{
    _canvas_enable_raw_data_ = mode; // TODO: this might not be needed, or if it is then the one in subclass might not be needed
    
    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->SetCanvasEnableRawData(mode);
    }
}

void Analysis::SetCanvasEnableDecayRate(const bool mode)
{
    _canvas_enable_decay_rate_ = mode; // TODO: this might not be needed, or if it is then the one in subclass might not be needed
    
    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->SetCanvasEnableDecayRate(mode);
    }
}





void Analysis::SetNumberOfPseudoexperiments(const Int_t number_of_pseudo_experiments_1d, const Int_t number_of_pseudo_experiments_2d)
{    
    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->SetNumberOfPseudoexperiments(number_of_pseudo_experiments_1d, number_of_pseudo_experiments_2d);
    }
}


////////////////////////////////////////////////////////////////////////////////
// STATIC CONSTANTS
////////////////////////////////////////////////////////////////////////////////

// Q value of decay, MeV
const Double_t Analysis::bb_Q{3.034};
