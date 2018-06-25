// self header
#include "Analysis.hpp"


// local headers
#include "ReWeight.hpp"
#include "aux.hpp"


Analysis::Analysis(const std::string& filename, const std::string& output_filename)
    : epsilon_31{0.368}
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
    , sensitivity_chisquare{0.0}
    , sensitivity_chisquare_2d{0.0}
    , fit_subrange{true} // fit subrange flag
    , number_of_pseudo_experiments{1}
    , number_of_pseudo_experiments_2d{1}
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
    , h_el_energy_original{nullptr}
    , h_el_energy_reweight{nullptr}
    , h_el_energy_sum_original{nullptr}
    , h_el_energy_sum_reweight{nullptr}
    , h_el_energy_diff{nullptr}
    , h_el_energy_pull{nullptr}
    , h_el_energy_data{nullptr}
    , h_el_energy_prob{nullptr}
    , h_el_energy_diff_data_rw{nullptr}
    , h_el_energy_diff_data_orig{nullptr}
    , h_el_energy_2d_original{nullptr}
    , h_el_energy_2d_reweight{nullptr}
    , h_el_energy_2d_diff{nullptr}
    , h_el_energy_2d_pull{nullptr}
    , h_el_energy_2d_data{nullptr}
    , h_el_energy_2d_prob{nullptr}
    , h_el_energy_2d_diff_data_rw{nullptr}
    , h_el_energy_2d_diff_data_orig{nullptr}
    , h_ll{nullptr}
    , h_ll_2d{nullptr}
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
{


    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////


}


Analysis::~Analysis()
{

}


////////////////////////////////////////////////////////////////////////////////
// ANALYSIS FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

/*
void Analysis::SetEpsilon31(const Double_t epsilon)
{
    epsilon_31 = epsilon;
}
*/

void Analysis::AddEpsilonValue(const Double_t epsilon)
{
    vec_epsilon_31.push_back(epsilon);
}

void Analysis::RunOverEpsilonVector()
{

    ReadData();
    CanvasDecayRate();
    CanvasSingleElectronProjection();
    CanvasSingleElectronTest();

    InitEventLoopTree();

    std::vector<Double_t>::const_iterator it{vec_epsilon_31.cbegin()};
    for(; it != vec_epsilon_31.cend(); ++ it)
    {

        epsilon_31 = *it;

        std::cout << "epsilon_31=" << epsilon_31 << std::endl;

        InitEventLoopHistogram();
        EventLoop();
        PostProcess();
        SummedEnergyFit();
        SensitivityMeasurementChisquare1();
        SensitivityMeasurementChisquare2();
        //SensitivityMeasurementLoglikelihood1();

        PrintOutputToFile();
    }


}


////////////////////////////////////////////////////////////////////////////////
// FLAG FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

void Analysis::SetBatchMode(const bool mode)
{
    _batch_mode_ = mode;
}

void Analysis::SetLogMode(const bool mode)
{
    log_mode = mode;
}

void Analysis::SetEnergyCutEnabled(const bool enabled)
{
    _energy_cut_enable_ = enabled;
}

void Analysis::SetGenWeightEnabled(const bool enabled)
{
    _gen_weight_enable_ = enabled;
}

void Analysis::SetCanvasEnableRawData(const bool mode)
{
    _canvas_enable_raw_data_ = mode;
}

void Analysis::SetCanvasEnableDecayRate(const bool mode)
{
    _canvas_enable_decay_rate_ = mode;
}



void Analysis::SetFitSubrange(const bool flag)
{
    fit_subrange = flag;
}


void Analysis::SetNumberOfPseudoexperiments(const Int_t number_of_pseudo_experiments_1d, const Int_t number_of_pseudo_experiments_2d)
{
    this->number_of_pseudo_experiments = number_of_pseudo_experiments_1d;
    this->number_of_pseudo_experiments_2d = number_of_pseudo_experiments_2d;
}


////////////////////////////////////////////////////////////////////////////////
// STATIC CONSTANTS
////////////////////////////////////////////////////////////////////////////////

// Q value of decay, MeV
const Double_t Analysis::bb_Q{3.034};
