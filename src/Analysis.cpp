// self header
#include "Analysis.hpp"


// local headers
#include "ReWeight.hpp"
#include "aux.hpp"


#include "TCanvas.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TLegend.h"


Analysis::Analysis(const std::string& filename, const std::string& output_filename)
    : epsilon_31{0.368}
    , systematic_energy_mult{1.0}
    , systematic_energy_mult_high{1.0}
    , systematic_energy_mult_low{1.0}
    , systematic_energy_mult_enabled{true}
    , systematic_energy_offset{0.0}
    , systematic_energy_offset_high{0.0}
    , systematic_energy_offset_low{0.0}
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
// CHISQUARE OUTPUT HISTOGRAM FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

void Analysis::MakeChiSquareType1()
{
    
    // get 1 sigma width for default, and minima for high/low

    // NOTE: chisquare results not stored in subclass
    //Double_t low_min{_subanalysis_systematic_low_->};
    //Double_t high_min{_subanalysis_systematic_high_->};

    // find minimum, low
    // get object
    //std::map<Double_t, SensitivityRecord> &s_rec_low{_sensitivity_record_map_.at()};
    //for(;;)


    // construct data for graphs

    std::size_t data_size_4{vec_epsilon_31.size()};
    std::size_t data_size_5{vec_epsilon_31.size()};
    std::size_t data_size_6{vec_epsilon_31.size()};

    Double_t *data_x_4 = new Double_t[data_size_4];
    Double_t *data_x_5 = new Double_t[data_size_5];
    Double_t *data_x_6 = new Double_t[data_size_6];

    Double_t *data_y_4 = new Double_t[data_size_4];
    Double_t *data_y_5 = new Double_t[data_size_5];
    Double_t *data_y_6 = new Double_t[data_size_6];

    for(std::size_t i{0}; i < data_size_4; ++ i)
    {
        Double_t epsilon_31{vec_epsilon_31.at(i)};
        std::cout << epsilon_31 << std::endl;
        std::cout << systematic_energy_mult << std::endl;

        data_x_4[i] = _sensitivity_record_map_.at(epsilon_31).at(systematic_energy_mult).epsilon_31;
        data_y_4[i] = _sensitivity_record_map_.at(epsilon_31).at(systematic_energy_mult).sensitivity_chisquare_1d_baseline_nosystematic;
    }
    for(std::size_t i{0}; i < data_size_5; ++ i)
    {
        Double_t epsilon_31{vec_epsilon_31.at(i)};
        std::cout << epsilon_31 << std::endl;
        std::cout << systematic_energy_mult << std::endl;

        data_x_5[i] = _sensitivity_record_map_.at(epsilon_31).at(systematic_energy_mult_high).epsilon_31;
        data_y_5[i] = _sensitivity_record_map_.at(epsilon_31).at(systematic_energy_mult_high).sensitivity_chisquare_1d_baseline_nosystematic;

        // TODO: this still doesn't work because we need the chisquare between the default and the systematic reweighted
        // not the systematic default and systematic reweighted
    }
    for(std::size_t i{0}; i < data_size_6; ++ i)
    {
        Double_t epsilon_31{vec_epsilon_31.at(i)};
        std::cout << epsilon_31 << std::endl;
        std::cout << systematic_energy_mult << std::endl;

        data_x_6[i] = _sensitivity_record_map_.at(epsilon_31).at(systematic_energy_mult_low).epsilon_31;
        data_y_6[i] = _sensitivity_record_map_.at(epsilon_31).at(systematic_energy_mult_low).sensitivity_chisquare_1d_baseline_nosystematic;
    }

    TGraph *g_4 = new TGraph(data_size_4, data_x_4, data_y_4); // FULL RANGE NO CUT DEFAULT
    g_4->SetTitle("");
    g_4->GetXaxis()->SetTitle("Parameter #xi_{31}^{2#nu}");
    g_4->GetYaxis()->SetTitle("#chi^{2}");
    //g_4->GetYaxis()->SetTitleOffset(1.75);
    //g_4->GetYaxis()->SetMaxDigits(3);
    TGaxis* g_4_tgaxis = (TGaxis*)g_4->GetYaxis();
    g_4_tgaxis->SetMaxDigits(3);
    g_4->SetLineColor(2);

    TGraph *g_5 = new TGraph(data_size_5, data_x_5, data_y_5); // FULL RANGE NO CUT HIGH
    g_5->SetTitle("");
    g_5->GetXaxis()->SetTitle("Parameter #xi_{31}^{2#nu}");
    g_5->GetYaxis()->SetTitle("#chi^{2}");
    //g_5->GetYaxis()->SetTitleOffset(1.75);
    //g_5->GetYaxis()->SetMaxDigits(3);
    TGaxis* g_5_tgaxis = (TGaxis*)g_5->GetYaxis();
    g_5_tgaxis->SetMaxDigits(3);
    g_5->SetLineColor(4);

    TGraph *g_6 = new TGraph(data_size_6, data_x_6, data_y_6); // FULL RANGE NO CUT LOW
    g_6->SetTitle("");
    g_6->GetXaxis()->SetTitle("Parameter #xi_{31}^{2#nu}");
    g_6->GetYaxis()->SetTitle("#chi^{2}");
    //g_6->GetYaxis()->SetTitleOffset(1.75);
    //g_6->GetYaxis()->SetMaxDigits(3);
    TGaxis* g_6_tgaxis = (TGaxis*)g_6->GetYaxis();
    g_6_tgaxis->SetMaxDigits(3);
    g_6->SetLineColor(6);

    TLegend *l_systematic = new TLegend(0.40, 0.70, 0.60, 0.85);
    l_systematic->AddEntry(g_4, "Default", "l"); // FULL RANGE NO CUT DEFAULT
    l_systematic->AddEntry(g_5, "High", "l"); // HIGH
    l_systematic->AddEntry(g_6, "Low", "l"); // LOW

    TCanvas *c_systematic = new TCanvas("c_systematic", "", 804, 628);
    c_systematic->GetPad(0)->SetTicks(1, 2);
    //g_4->GetYaxis()->SetRangeUser(0.0, 6.0e1); // energy systematic m
    g_4->GetYaxis()->SetRangeUser(0.0, 1.5e3); // energy systematic c
    g_4->Draw("AL");
    g_5->Draw("same");
    g_6->Draw("same");
    l_systematic->Draw();
    c_systematic->SaveAs("c_systematic_out.png");
    c_systematic->SaveAs("c_systematic_out.pdf");
    c_systematic->SaveAs("c_systematic_out.eps");
    c_systematic->SaveAs("c_systematic_out.C");

    TFile *c_systematic_out = new TFile("c_systematic_out.root", "RECREATE");
    g_4->SetName("g_systematic_default");
    g_5->SetName("g_systematic_high");
    g_6->SetName("g_systematic_low");
    g_4->Write();
    g_5->Write();
    g_6->Write();
    c_systematic_out->Close();
    delete c_systematic_out;
    c_systematic_out = nullptr;

    // get minimum, low, high
    Double_t min_low{data_y_6[0]};
    Int_t min_low_ix{0};
    Double_t min_high{data_y_5[0]};
    Int_t min_high_ix{0};
    Double_t min_default{data_y_4[0]};
    Int_t min_default_ix{0};
    Double_t sigma_1_low{0.0};
    Int_t sigma_1_low_ix{0};
    Double_t sigma_1_high{0.0};
    Int_t sigma_1_high_ix{0};

    for(std::size_t i{0}; i < data_size_4; ++ i)
    {
        Double_t val_low{data_y_6[i]};
        Double_t val_high{data_y_5[i]};
        Double_t val_default{data_y_4[i]};

        if(val_low < min_low)
        {
            min_low = val_low;
            min_low_ix = i;
        }
        if(val_high < min_high)
        {
            min_high = val_high;
            min_high_ix = i;
        }
        if(val_default < min_default)
        {
            min_default = val_default;
            min_default_ix = i;
        }

    }

    // find 1 sigma width and interpolate
    Int_t x0_index{min_default_ix};
    Int_t data_size{data_size_4};
    Double_t *data_y{data_y_4};
    Double_t *data_x{data_x_4};
    Double_t y0{min_default};
    Double_t x_high;
    Double_t x_low;
    Double_t x_high_simple;
    Double_t x_low_simple;
    for(Int_t ix{x0_index}; ix < data_size - 1; ++ ix)
    {
        if(data_y[ix] - y0 >= 1.0)
        {
            //std::cout << "found delta >= 1 sigma: delta=" << data_y[ix] - y0 << std::endl;
            x_high_simple = data_x[ix];
            Double_t y_max_1sigma{y0 + 1.0};
            Double_t y_minus{data_y[ix - 1]};
            Double_t y_plus{data_y[ix]};
            Double_t x_minus{data_x[ix - 1]};
            Double_t x_plus{data_x[ix]};
            x_high = x_minus + (x_plus - x_minus) * ((y_max_1sigma - y_minus) / (y_plus - y_minus));
            std::cout << "1 sigma interpolate, high: " << x_minus << " " << x_high << " " << x_plus << " simple=" << x_high_simple << std::endl;
            std::cout << "the index is: " << ix << ", the initial index is: " << x0_index << std::endl;
            sigma_1_high = x_high;
            sigma_1_high_ix = ix; // somewhere in the middle, due to interpolation
            break;
        }
    }
    for(Int_t ix{x0_index}; ix >= 0; -- ix)
    {
        if(data_y[ix] - y0 >= 1.0)
        {
            //std::cout << "found delta >= 1 sigma: delta=" << data_y[ix] - y0 << std::endl;
            x_low_simple = data_x[ix];
            Double_t y_max_1sigma{y0 + 1.0};
            Double_t y_minus{data_y[ix]};
            Double_t y_plus{data_y[ix + 1]};
            Double_t x_minus{data_x[ix]};
            Double_t x_plus{data_x[ix + 1]};
            x_low = x_minus + (x_plus - x_minus) * ((y_max_1sigma - y_minus) / (y_plus - y_minus));
            std::cout << "1 sigma interpolate, low: " << x_minus << " " << x_low << " " << x_plus << " simple=" << x_low_simple << std::endl;
            std::cout << "the index is: " << ix << ", the initial index is: " << x0_index << std::endl;
            sigma_1_low = x_low;
            sigma_1_low_ix = ix;
            break;
        }
    }
    
    std::cout << "width is: " << x_high - x_low << std::endl;
    std::cout << "systematics width is: " << data_x_6[min_low_ix] << " - " << data_x_5[min_high_ix] << " = " << data_x_6[min_low_ix] - data_x_5[min_high_ix] << std::endl;
    // TODO: NEED TO CHECK INTERPOLATION ALGORITHM
    // SEE ALSO: int main() of sensitivity.cpp

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

void Analysis::SetSystematicEnergyMultiplierEnabled(const bool enabled)
{
    systematic_energy_mult_enabled = enabled;
}

void Analysis::SetSystematicEnergyOffset(const Double_t high, const Double_t low)
{
    systematic_energy_offset_low = low;
    systematic_energy_offset_high = high;
}

// TODO: add new functions here, for systematics


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

        for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
        {
            (*it)->ChiSquare_BaselineNoSystematic();
        }
        
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

    MakeChiSquareType1();

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
