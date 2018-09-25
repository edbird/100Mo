#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP


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


class Analysis
{


    public:

    ////////////////////////////////////////////////////////////////////////////
    // CLASS FUNCTIONS
    ////////////////////////////////////////////////////////////////////////////

    Analysis(const std::string&, const std::string&);
    ~Analysis();


    ////////////////////////////////////////////////////////////////////////////
    // CHISQUARE OUTPUT FUNCTIONS
    ////////////////////////////////////////////////////////////////////////////

    void MakeChiSquareType1();

    ////////////////////////////////////////////////////////////////////////////
    // ANALYSIS FUNCTIONS
    ////////////////////////////////////////////////////////////////////////////

    void SetSystematicEnergyMultiplier(const Double_t);
    void SetSystematicEnergyMultiplierHighLow(const Double_t, const Double_t);
    void SetSystematicEnergyMultiplierEnabled(const bool);

    void SetSystematicEnergyOffset(const Double_t, const Double_t);
    void SetSystematicEfficiency(const Double_t, const Double_t);

    void SetEnergyCorrectionSystematicEnabled(const bool);

    void SetEpsilon31(const Double_t);
    void AddEpsilonValue(const Double_t epsilon);
    void RunOverEpsilonVector();

    void ReadData();
    //void CanvasRawData();
    void CanvasDecayRate();
    void CanvasSingleElectronProjection(); // TODO: note this creates histograms as well as canvas output, misleading name
    void CanvasSingleElectronTest(); // graphs of data from theory paper 
    void InitSingleElectronTest();
    void InitEventLoopTree();
    void InitEventLoopHistogram();
    void EventLoop();
    void PostProcess();
    void SummedEnergyFit();
    void SensitivityMeasurementChisquare1(); // 1d
    void SensitivityMeasurementChisquare2(); // 2d
    void SensitivityMeasurementLoglikelihood1();
    void SensitivityMeasurementLoglikelihood2();
    void PrintOutputToFile();
    void MakeSensitivityCanvas();


    ////////////////////////////////////////////////////////////////////////////
    // CALIBRATION FUNCTIONS
    ////////////////////////////////////////////////////////////////////////////
    
    void SetEnergyCalibrationPoints(const Double_t expected_Bi207_EC_1 = 481.7, const Double_t expected_Bi207_EC_2 = 975.7)
    {
        // reference: https://www.ezag.com/fileadmin/ezag/user-uploads/isotopes/isotopes/Isotrak/isotrak-pdf/Decay_Schema_Data/Bi-207.pdf

        // set the calibration points
        energy_calibration_Bi207_EC_1 = expected_Bi207_EC_1; // expected Electron Capture peak 1 (keV)
        energy_calibration_Bi207_EC_2 = expected_Bi207_EC_2; // expected Electron Capture peak 2 (keV)
    }

    void RunEnergyCalibration(const Double_t measured_Bi207_EC_1, const Double_t measured_Bi207_EC_2)
    {
        // set the calibration points
        energy_calibration_Bi207_EC_measured_1 = measured_Bi207_EC_1; // measured EC peak 1 (keV)
        energy_calibration_Bi207_EC_measured_2 = measured_Bi207_EC_2; // measured EC peak 2 (keV)

        // assume energy calibration points are set

        Double_t E1{energy_calibration_Bi207_EC_1};
        Double_t E2{energy_calibration_Bi207_EC_2};
        Double_t X1{energy_calibration_Bi207_EC_measured_1};
        Double_t X2{energy_calibration_Bi207_EC_measured_2};

        Double_t a{(E1 - E2) / (X1 - X2)};
        Double_t b{E1 - a * X1};
        if(std::abs((E2 - a * X2) - b) > 1e-4)
        {
            std::cout << "WARNING: CALIBRATION FAILED" << std::endl;
        }

        energy_calibration_a = a;
        energy_calibration_b = b;

        std::cout << "Energy Calibration Values: a = " << a << ", b = " << b << std::endl;
    }


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


    ////////////////////////////////////////////////////////////////////////////
    // CLASSES
    ////////////////////////////////////////////////////////////////////////////
    
    // low and high systematic energy shift versions
    std::vector<SubAnalysis*> _subanalysis_;
    SubAnalysis* _subanalysis_systematic_low_;
    SubAnalysis* _subanalysis_systematic_high_;
    SubAnalysis* _subanalysis_systematic_default_;

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
    
    bool fit_subrange;

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
    std::map<Double_t, std::map<Double_t, SensitivityRecord>> _sensitivity_record_map_;


    ////////////////////////////////////////////////////////////////////////////
    // DATA
    ////////////////////////////////////////////////////////////////////////////

    std::vector<Double_t> vec_epsilon_31;
    Double_t epsilon_31;

    // systematic energy shift (this is the default value, 1.0)
    Double_t systematic_energy_mult;
    // systematic energy shift, low level
    Double_t systematic_energy_mult_low;
    // systematic energy shift, high level
    Double_t systematic_energy_mult_high;
    // enable / disable switch, this is explained in comment in main.cpp
    bool systematic_energy_mult_enabled;

    // systematic energy offset
    Double_t systematic_energy_offset;
    Double_t systematic_energy_offset_low;
    Double_t systematic_energy_offset_high;

    // systematic efficiency
    Double_t systematic_efficiency;
    Double_t systematic_efficiency_low;
    Double_t systematic_efficiency_high;

    // energy calibration, Bi 207
    Double_t energy_calibration_a;
    Double_t energy_calibration_b;
    // TODO: systematics

    Double_t energy_calibration_Bi207_EC_1; // expected Electron Capture peak 1 (keV)
    Double_t energy_calibration_Bi207_EC_2; // expected Electron Capture peak 2 (keV)
    Double_t energy_calibration_Bi207_EC_measured_1; // measured EC peak 1 (keV)
    Double_t energy_calibration_Bi207_EC_measured_2; // measured EC peak 2 (keV)
    // TODO: uncertainties (systematics?)

    // energy correction file
    // NOTE: FOR MC ONLY, DO NOT APPLY TO DATA
    TFile *f_energy_correction;
    TGraphErrors *g_energy_correction;
    MODE_FLAG m_mode; // simulation (MC) or real (DATA)
    bool b_energy_correction_systematic_enabled;

    TFile *f;
    TTree *t;

    // 0.0 MeV to 4.0 MeV = 4.0 MeV range
    // num_bins keV bin width: 4.0 MeV / 0.1 MeV = 40 bins
    Int_t num_bins;

    std::string filename;
    std::string output_filename; // TODO: remove?

    // allocate data storage for raw data
    std::vector<std::vector<double>> data_nEqNull;
    std::vector<std::vector<double>> data_nEqTwo;
    double psiN0;
    double psiN2;

    // for 2d decay rate histograms
    std::vector<std::vector<double>> data_0; // epsilon_31 = 0.0
    std::vector<std::vector<double>> data_1; // epsilon_31 = 0.4
    std::vector<std::vector<double>> data_2; // epsilon_31 = 0.8

    



    // random generator
    TRandom3 gen;


    ////////////////////////////////////////////////////////////////////////////
    // EVENT LOOP DATA
    ////////////////////////////////////////////////////////////////////////////

    Int_t nElectrons;
    Double_t trueT1;
    Double_t trueT2;
    Double_t el_energy_[2];
    Double_t gen_weight;
    
    

    ////////////////////////////////////////////////////////////////////////////
    // OUTPUT STREAMS
    ////////////////////////////////////////////////////////////////////////////

    
    ////////////////////////////////////////////////////////////////////////////
    // GRAPHS
    ////////////////////////////////////////////////////////////////////////////

    // data from theory paper for test_single / test_sum
    TGraph *g_el_energy_single_0;
    TGraph *g_el_energy_single_1;
    TGraph *g_el_energy_single_2;
    TGraph *g_el_energy_sum_0;
    TGraph *g_el_energy_sum_1;
    TGraph *g_el_energy_sum_2;


    ////////////////////////////////////////////////////////////////////////////
    // FUNCTIONS
    ////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////

    // data histograms
    TH2D *h_nEqNull;
    TH2D *h_nEqTwo;
    
    TH2D *h_data_0;
    TH2D *h_data_1;
    TH2D *h_data_2;

    // single electron projection histograms
    TH1D *h_single_electron_0;
    TH1D *h_single_electron_1;
    TH1D *h_single_electron_2;

    // test histograms, curves from paper plotted with original and reweighed
    // these use truth energies rather than reconstructed energy
    // so ignore detector effects
    TH1D *h_test_single_original;
    TH1D *h_test_single_reweight;
    TH1D *h_test_sum_original;
    TH1D *h_test_sum_reweight;


    // general weight histogram
    TH2D *h_gen_weight;


    ////////////////////////////////////////////////////////////////////////////
    // CANVAS
    ////////////////////////////////////////////////////////////////////////////

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



};

#endif // ANALYSIS_HPP
