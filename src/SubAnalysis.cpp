// self header
#include "Analysis.hpp"


// local headers
#include "ReWeight.hpp"
#include "aux.hpp"


SubAnalysis::SubAnalysis
    (
        const std::string& h_name_append,
        const std::string& output_filename,
        std::map<Double_t, std::map<Double_t, SensitivityRecord>> *sensitivity_record_map,
        Double_t epsilon_31,
        Double_t systematic_energy_mult, bool systematic_energy_mult_enable,
        Double_t systematic_energy_offset,
        TH2D *h_nEqNull, TH2D *h_nEqTwo, Double_t psiN0, Double_t psiN2,
        Int_t *nElectrons, Double_t *trueT1, Double_t *trueT2, Double_t *el_energy_, Double_t *gen_weight,
        TRandom3 *gen
    )
    : h_name_append{h_name_append}
    , _sensitivity_record_map_{sensitivity_record_map}
    , epsilon_31{epsilon_31}
    , systematic_energy_mult{systematic_energy_mult}
    , systematic_energy_mult_enable{systematic_energy_mult_enable}
    , systematic_energy_offset{systematic_energy_offset}

    , nElectrons{nElectrons}
    , trueT1{trueT1}
    , trueT2{trueT2}
    , el_energy_{el_energy_}
    , gen_weight{gen_weight}
    , gen{gen}
    
    , h_nEqNull{h_nEqNull}
    , h_nEqTwo{h_nEqTwo}
    , psiN0{psiN0}
    , psiN2{psiN2}
    
    , num_bins{40}
    , output_filename{output_filename}
    , _batch_mode_{false}
    , log_mode{false}
    , _gen_weight_enable_{false}
    , _energy_cut_enable_{false} // energy cut flag
    , fit_subrange{false} // fit subrange flag
    , _canvas_enable_raw_data_{false}
    , _canvas_enable_decay_rate_{true}
    , _canvas_enable_single_electron_projection_{false}


    
    , sensitivity_chisquare{0.0}
    , sensitivity_chisquare_2d{0.0}
    , number_of_pseudo_experiments{1}
    , number_of_pseudo_experiments_2d{1}
    , sensitivity_chisquare_baseline_nosystematic{0.0}
    , non_empty_bins_baseline_nosystematic{0}


    , f_el_energy_sum_original{nullptr}
    
    , h_el_energy_original{nullptr}
    , h_el_energy_reweight{nullptr}
    , h_el_energy_sum_original{nullptr}
    , h_el_energy_sum_reweight{nullptr}
    , h_el_energy_diff{nullptr}
    , h_el_energy_pull{nullptr}
    , h_el_energy_ratio{nullptr}
    , h_el_energy_data{nullptr}
    , h_el_energy_prob{nullptr}
    , h_el_energy_diff_data_rw{nullptr}
    , h_el_energy_diff_data_orig{nullptr}
    , h_el_energy_2d_original{nullptr}
    , h_el_energy_2d_reweight{nullptr}
    , h_el_energy_2d_diff{nullptr}
    , h_el_energy_2d_pull{nullptr}
    , h_el_energy_2d_ratio{nullptr}
    , h_el_energy_2d_data{nullptr}
    , h_el_energy_2d_prob{nullptr}
    , h_el_energy_2d_diff_data_rw{nullptr}
    , h_el_energy_2d_diff_data_orig{nullptr}
    , h_ll{nullptr}
    , h_ll_2d{nullptr}
    
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
    
    , _baseline_histo_p_{nullptr}
{

    std::cout << "construct: " << h_name_append << std::endl;

}


SubAnalysis::~SubAnalysis()
{

    std::cout << "delete: " << h_name_append << std::endl;

    delete f_el_energy_sum_original;

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

    delete h_el_energy_data;
    delete h_el_energy_prob;
    delete h_el_energy_diff_data_rw;
    delete h_el_energy_diff_data_orig;

    delete h_el_energy_2d_data;
    delete h_el_energy_2d_prob;
    delete h_el_energy_2d_diff_data_rw;
    delete h_el_energy_2d_diff_data_orig;

    delete h_ll;
    delete h_ll_2d;
    



}


void SubAnalysis::SetBaselineHistoPointer(TH1D* const histo_p)
{
    _baseline_histo_p_ = histo_p;
}

TH1D* const SubAnalysis::GetBaselineHistoPointer()
{
    return _baseline_histo_p_;
}

TH1D* const SubAnalysis::GetElEnergyOriginalHistoPointer()
{
    return h_el_energy_original;
}


////////////////////////////////////////////////////////////////////////////////
// ANALYSIS FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

/*(
void SubAnalysis::SetSystematicEnergyMultiplier(const Double_t value)
{
    systematic_energy_mult = value;
}
*/

/*
void SubAnalysis::SetEpsilon31(const Double_t epsilon)
{
    epsilon_31 = epsilon;
}
*/


////////////////////////////////////////////////////////////////////////////////
// FLAG FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

void SubAnalysis::SetBatchMode(const bool mode)
{
    _batch_mode_ = mode;
}

void SubAnalysis::SetLogMode(const bool mode)
{
    log_mode = mode;
}

void SubAnalysis::SetEnergyCutEnabled(const bool enabled)
{
    _energy_cut_enable_ = enabled;
}

void SubAnalysis::SetGenWeightEnabled(const bool enabled)
{
    _gen_weight_enable_ = enabled;
}

void SubAnalysis::SetCanvasEnableRawData(const bool mode)
{
    _canvas_enable_raw_data_ = mode;
}

void SubAnalysis::SetCanvasEnableDecayRate(const bool mode)
{
    _canvas_enable_decay_rate_ = mode;
}



void SubAnalysis::SetFitSubrange(const bool flag)
{
    fit_subrange = flag;
}


void SubAnalysis::SetNumberOfPseudoexperiments(const Int_t number_of_pseudo_experiments_1d, const Int_t number_of_pseudo_experiments_2d)
{
    this->number_of_pseudo_experiments = number_of_pseudo_experiments_1d;
    this->number_of_pseudo_experiments_2d = number_of_pseudo_experiments_2d;
}


////////////////////////////////////////////////////////////////////////////////
// STATIC CONSTANTS
////////////////////////////////////////////////////////////////////////////////

// Q value of decay, MeV
const Double_t SubAnalysis::bb_Q{3.034};
