#ifndef SUBANALYSIS_HPP
#define SUBANALYSIS_HPP


#include "SensitivityRecord.hpp"
#include "CanvasFactory.hpp"
#include "SubAnalysis.hpp"
#include "Flag.hpp"


#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"


#include <iostream>
#include <map>


class Analysis;

// this class runs for each systematic energy multiplier value
class SubAnalysis
{

    friend Analysis;

    public:
    
    ////////////////////////////////////////////////////////////////////////////
    // CLASS FUNCTIONS
    ////////////////////////////////////////////////////////////////////////////

    SubAnalysis(const std::string&, const std::string&, std::map<Double_t,
                std::map<Double_t, SensitivityRecord>>*,
                Double_t, Double_t, bool, Double_t, Double_t,
                TH2D*, TH2D*, Double_t, Double_t,
                Int_t*, Double_t*, Double_t*, Double_t*, Double_t*, TRandom3*);
    ~SubAnalysis();
    
    
    ////////////////////////////////////////////////////////////////////////////
    // ANALYSIS FUNCTIONS
    ////////////////////////////////////////////////////////////////////////////

    // done by constructor
    //void SetSystematicEnergyMultiplier(const Double_t);
    //void SetEpsilon31(const Double_t);
    
    
    void InitEventLoopHistogram();
    void EventLoop();
    void PostProcess();
    void SummedEnergyFit();
    void SensitivityMeasurementChisquare1(); // 1d
    void SensitivityMeasurementChisquare2(); // 2d
    void SensitivityMeasurementLoglikelihood1();
    void SensitivityMeasurementLoglikelihood2();
    void PrintOutputToFile();

    void ChiSquare_BaselineNoSystematic();
    

    ////////////////////////////////////////////////////////////////////////////
    // CALIBRATION FUNCTIONS
    ////////////////////////////////////////////////////////////////////////////
    
    void SetEnergyCalibrationConstants(const Double_t a, const Double_t b)
    {
        energy_calibration_a = a;
        energy_calibration_b = b;
    }

    void Set_g_energy_correction(TGraphErrors* const);
    void SetMode(const MODE_FLAG);

    void SetEnergyCorrectionSystematicMode(const Int_t);
    void SetEnergyCorrectionSystematicEnabled(const bool);


    ////////////////////////////////////////////////////////////////////////////
    // FLAG FUNCTIONS
    ////////////////////////////////////////////////////////////////////////////

    void SetBatchMode(const bool mode);
    void SetLogMode(const bool mode);
    void SetEnergyCutEnabled(const bool enabled);
    void SetGenWeightEnabled(const bool enabled);
    void SetCanvasEnableRawData(const bool mode);
    void SetCanvasEnableDecayRate(const bool mode);
    void SetCanvasEnableSingleElectronProjection(const bool mode);
    
    
    void SetFitSubrange(const bool);
    void SetNumberOfPseudoexperiments(const Int_t, const Int_t);
    

    void SetBaselineHistoPointer(TH1D* const);
    TH1D* const GetBaselineHistoPointer();
    TH1D* const GetElEnergyOriginalHistoPointer();

    
    ////////////////////////////////////////////////////////////////////////////
    // CONSTANTS
    ////////////////////////////////////////////////////////////////////////////

    static const Double_t bb_Q; // Q value of decay, MeV

    const Int_t dimension_xy{1001}; // dimension of 2d decay rate histograms


    ////////////////////////////////////////////////////////////////////////////
    // FLAGS
    ////////////////////////////////////////////////////////////////////////////

    bool _batch_mode_;
    bool log_mode;
    bool _gen_weight_enable_;
    bool _energy_cut_enable_;

    // canvas draw flags
    bool _canvas_enable_raw_data_;
    bool _canvas_enable_decay_rate_;
    bool _canvas_enable_single_electron_projection_;

    
    ////////////////////////////////////////////////////////////////////////////
    // OUTPUT AND RESULTS
    ////////////////////////////////////////////////////////////////////////////

    // map: first index: epsilon_31 value
    // second index: systematic energy multiplier value
    // SensitivityRecord object contains these values, and they should
    // be equal, if not there is an error
    std::map<Double_t, std::map<Double_t, SensitivityRecord>> *_sensitivity_record_map_;
    
    
    ////////////////////////////////////////////////////////////////////////////
    // DATA
    ////////////////////////////////////////////////////////////////////////////
    
    // string to append to histogram names to avoid name conflict
    std::string h_name_append;

    // don't have f or t (file or tree),
    // instead have pointers to tree variables (below)
    
    
    Double_t epsilon_31;
    
    // systematic energy shift
    Double_t systematic_energy_mult; // (this is the current value for this class)
    bool systematic_energy_mult_enable;
    Double_t systematic_energy_offset; 
    Double_t systematic_efficiency; 

    // energy calibration, Bi 207
    Double_t energy_calibration_a;
    Double_t energy_calibration_b;

    // moved to Analysis
    //Double_t energy_calibration_Bi207_EC_1; // expected Electron Capture peak 1 (keV)
    //Double_t energy_calibration_Bi207_EC_2; // expected Electron Capture peak 2 (keV)
    //Double_t energy_calibration_Bi207_EC_measured_1; // measured EC peak 1 (keV)
    //Double_t energy_calibration_Bi207_EC_measured_2; // measured EC peak 2 (keV)

    // energy correction
    // NOTE: FOR MC ONLY, DO NOT APPLY TO DATA
    TGraphErrors *g_energy_correction;
    TGraph *g_energy_correction_systematic_high;
    TGraph *g_energy_correction_systematic_low;
    MODE_FLAG m_mode;
    Int_t m_energy_correction_systematic_mode;
    //  0 = none,
    //  1 = systematic positive (add uncertainties)
    // -1 = systematic negative (subtract uncertainties)
    bool b_energy_correction_systematic_enabled;
    
    // 0.0 MeV to 4.0 MeV = 4.0 MeV range
    // num_bins keV bin width: 4.0 MeV / 0.1 MeV = 40 bins
    Int_t num_bins;

    std::string output_filename;
    
    // phase space values
    double psiN0;
    double psiN2;
    
    
    // sensitivity measurement
    Double_t sensitivity_chisquare{0.0};
    Double_t sensitivity_chisquare_2d{0.0};
    // number of degrees of freedom
    Int_t non_empty_bins{0}; // 1d single electron energy
    Int_t non_empty_bins_2d{0};
    Int_t non_empty_bins_fit{0}; // summed electron energy

    bool fit_subrange;
    
    // chisquare between reweighted with arbitary xi and systematics
    // and baseline with baseline xi (0.368) and no systematics
    Double_t sensitivity_chisquare_baseline_nosystematic;
    Int_t non_empty_bins_baseline_nosystematic;
    
    
    // random generator
    TRandom3 *gen;
    
    // distribution of log-likelihood
    std::vector<Double_t> vec_ll;
    std::vector<Double_t> vec_ll_2d;
    Int_t number_of_pseudo_experiments; //{1000000};
    Int_t number_of_pseudo_experiments_2d; //{1000000};
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    // EVENT LOOP DATA
    ////////////////////////////////////////////////////////////////////////////

    Int_t *nElectrons;
    Double_t *trueT1;
    Double_t *trueT2;
    Double_t *el_energy_;
    Double_t *gen_weight;
    
    
    ////////////////////////////////////////////////////////////////////////////
    // OUTPUT STREAMS
    ////////////////////////////////////////////////////////////////////////////

    
    ////////////////////////////////////////////////////////////////////////////
    // GRAPHS
    ////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////
    // FUNCTIONS
    ////////////////////////////////////////////////////////////////////////////

    TF1 *f_el_energy_sum_original;
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////

    // data histograms
    TH2D *h_nEqNull;
    TH2D *h_nEqTwo;


    // original electron energy (data)
    TH1D *h_el_energy_original;
    // reweighted electron energy after reweighting by epsilon (model)
    TH1D *h_el_energy_reweight;
    
    // original summed electron energy (data)
    TH1D *h_el_energy_sum_original;
    // reweigted summed electron energy (model)
    TH1D *h_el_energy_sum_reweight;

    // difference: electron energy, original - reweighted
    TH1D *h_el_energy_diff;
    // pull
    TH1D *h_el_energy_pull;
    // ratio
    TH1D *h_el_energy_ratio;

    // 2d single original electron energy
    TH2D *h_el_energy_2d_original;
    // 2d reweighted single electron energy
    TH2D *h_el_energy_2d_reweight;

    // 2d orig - 2d reweight (data - model)
    TH2D *h_el_energy_2d_diff;
    // chisquare pull
    TH2D *h_el_energy_2d_pull;
    // ratio
    TH2D *h_el_energy_2d_ratio;
    
    
    // log likelihood pseudoexperiments method
    TH1I *h_el_energy_data;
    TH1D *h_el_energy_prob;
    // data - reweighted
    TH1D *h_el_energy_diff_data_rw;
    // data - original
    TH1D *h_el_energy_diff_data_orig;

    // log likelihood pseudoexperiments method
    TH2I *h_el_energy_2d_data;
    TH2D *h_el_energy_2d_prob;
    // 2d data (poisson log likelihood method) - 2d reweight
    TH2D *h_el_energy_2d_diff_data_rw;
    // 2d data (poisson log likelihood method) - 2d original
    TH2D *h_el_energy_2d_diff_data_orig;
    
    
    // log likelihood histograms
    TH1D *h_ll;
    TH1D *h_ll_2d;
    
    
    ////////////////////////////////////////////////////////////////////////////
    // CANVAS
    ////////////////////////////////////////////////////////////////////////////
    
    // not sure if I want anything here yet
    
    TCanvas *c_nEqNull;
    TCanvas *c_nEqTwo;
    TCanvas *c_data_0;
    TCanvas *c_data_1;
    TCanvas *c_data_2;
    TCanvas *c_single_electron;
    TCanvas *c_test_single;
    TCanvas *c_test_sum;
    TCanvas *c_gen_weight;
    TCanvas *c_el_energy_diff;
    TCanvas *c_el_energy_pull;
    TCanvas *c_el_energy_diff_data_rw;
    TCanvas *c_el_energy_diff_data_orig;
    TCanvas *c_el_energy_2d_original;
    TCanvas *c_el_energy_2d_reweight;
    TCanvas *c_el_energy_2d_diff;
    TCanvas *c_el_energy_2d_pull;
    TCanvas *c_ll;
    TCanvas *c_el_energy_2d_diff_data_rw;
    TCanvas *c_el_energy_2d_diff_data_orig;
    TCanvas *c_el_energy_2d_data;
    TCanvas *c_el_energy_2d_prob;
    TCanvas *c_ll_2d;


    // pointer to the "baseline histogram"
    // this is the histogram of the class where all systematics are set
    // to default values, and xi is set to the baseline value (0.368)
    TH1D *_baseline_histo_p_;

};


#endif // SUBANALYSIS_HPP
