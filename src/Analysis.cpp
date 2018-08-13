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

    TLegend *l_systematic = new TLegend(0.45, 0.70, 0.65, 0.85);
    l_systematic->AddEntry(g_4, "Default", "l"); // FULL RANGE NO CUT DEFAULT
    l_systematic->AddEntry(g_5, "High", "l"); // HIGH
    l_systematic->AddEntry(g_6, "Low", "l"); // LOW

    TCanvas *c_systematic = new TCanvas("c_systematic", "", 804, 628);
    c_systematic->GetPad(0)->SetTicks(1, 2);
    g_4->GetYaxis()->SetRangeUser(0.0, 1.6e3);
    g_4->Draw("AL");
    g_5->Draw("same");
    g_6->Draw("same");
    l_systematic->Draw();
    c_systematic->SaveAs("c_systematic_out.png");
    c_systematic->SaveAs("c_systematic_out.pdf");
    c_systematic->SaveAs("c_systematic_out.eps");
    c_systematic->SaveAs("c_systematic_out.C");


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
