


#include "TCanvas.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"



#include <iostream>


#include "read_data.hpp"

#include "aux.hpp"

#include "program_arguments.hpp"

//#include "CanvasFactory.hpp"

//#include "HistogramWrapper.hpp"
// TODO: make some nice histograms
// X/Y AXIS: LABEL: TEXT, FONT SIZE, FONT
// X/Y AXIS: TICK NUMBERS: FONT SIZE, FONT
// X/Y AXIS LIMITS, RANGE, PARTICULARLY Y MAX
// GRAPH COLORS



// 0.0 MeV to 4.0 MeV = 4.0 MeV range
// num_bins keV bin width: 4.0 MeV / 0.1 MeV = 40 bins
Int_t num_bins{40};

int main(int argc, char* argv[])
{
    
    ////////////////////////////////////////////////////////////////////////////
    // PROCESS PROGRAM ARGUMENTS
    ////////////////////////////////////////////////////////////////////////////

    // new program arguments processing method
    ProgramArguments pa;
    pa.Add("help", "--help", "false");
    //pa.Add("filename", "--filename", "");
    pa.Print();
    pa.Process(argc, argv);

    // this argument is a string
    std::string help{pa.Get("help")};
    // this argument is a string
    //std::string filename{pa.Get("filename")};


    // Q value of decay
    // MeV
    Double_t bb_Q{3.034};

    ////////////////////////////////////////////////////////////////////////////
    // READ DATA IN FROM FILES
    ////////////////////////////////////////////////////////////////////////////

    // allocate data storage
    std::vector<std::vector<double>> data_FULLRANGE_CUT;
    std::vector<std::vector<double>> data_FULLRANGE_NOCUT;
    std::vector<std::vector<double>> data_SUBRANGE_CUT;
    std::vector<std::vector<double>> data_SUBRANGE_NOCUT;

    // read
    read_data_helper_2("of_data_FULLRANGE_CUT.txt", data_FULLRANGE_CUT, ',');
    read_data_helper_2("of_data_FULLRANGE_NOCUT.txt", data_FULLRANGE_NOCUT, ',');
    read_data_helper_2("of_data_SUBRANGE_CUT.txt", data_SUBRANGE_CUT, ',');
    read_data_helper_2("of_data_SUBRANGE_NOCUT.txt", data_SUBRANGE_NOCUT, ',');

    std::cout << data_FULLRANGE_CUT.size() << std::endl;
    
    ////////////////////////////////////////////////////////////////////////////
    // CONVERT INPUT DATA TO GRAPH FORMAT
    ////////////////////////////////////////////////////////////////////////////

    Int_t data_size_0{data_FULLRANGE_CUT.size()};
    Int_t data_size_1{data_FULLRANGE_NOCUT.size()};
    Int_t data_size_2{data_SUBRANGE_CUT.size()};
    Int_t data_size_3{data_SUBRANGE_NOCUT.size()};

    Double_t *data_x_0 = new Double_t[data_size_0];
    Double_t *data_x_1 = new Double_t[data_size_1];
    Double_t *data_x_2 = new Double_t[data_size_2];
    Double_t *data_x_3 = new Double_t[data_size_3];

    Double_t *data_y_0 = new Double_t[data_size_0];
    Double_t *data_y_1 = new Double_t[data_size_1];
    Double_t *data_y_2 = new Double_t[data_size_2];
    Double_t *data_y_3 = new Double_t[data_size_3];

    // SCALE THE DATA SO THAT X RUNS FROM 0.0 TO Q_BB RATHER THAN 0.0 TO 1.0
    // SCALE THE DATA SO THAT THE INTEGRAL IS UNCHANGED, IE: INTEGRAL = 1.0
    // THIS IS DONE BY SCALING THE Y VALUES IN THE INVERSE WAY TO THE SCALING
    // OF THE X VALUES

    for(std::size_t i{0}; i < data_size_0; ++ i)
    {
        data_x_0[i] = data_FULLRANGE_CUT[i][1];
        data_y_0[i] = data_FULLRANGE_CUT[i][3];
    }

    for(std::size_t i{0}; i < data_size_1; ++ i)
    {
        data_x_1[i] = data_FULLRANGE_NOCUT[i][1];
        data_y_1[i] = data_FULLRANGE_NOCUT[i][3];
    }

    for(std::size_t i{0}; i < data_size_2; ++ i)
    {
        data_x_2[i] = data_SUBRANGE_CUT[i][1];
        data_y_2[i] = data_SUBRANGE_CUT[i][3];
    }

    for(std::size_t i{0}; i < data_size_3; ++ i)
    {
        data_x_3[i] = data_SUBRANGE_NOCUT[i][1];
        data_y_3[i] = data_SUBRANGE_NOCUT[i][3];
    }

    TGraph *g_0 = new TGraph(data_size_0, data_x_0, data_y_0);
    g_0->SetTitle("");
    g_0->GetXaxis()->SetTitle("Parameter #xi");
    g_0->GetYaxis()->SetTitle("#chi^{2}");
    //g_0->GetYaxis()->SetTitleOffset(1.75);
    //g_0->GetYaxis()->SetMaxDigits(3);
    TGaxis* g_0_tgaxis = (TGaxis*)g_0->GetYaxis();
    g_0_tgaxis->SetMaxDigits(3);
    g_0->SetLineColor(2);

    TGraph *g_1 = new TGraph(data_size_1, data_x_1, data_y_1);
    g_1->SetTitle("");
    g_1->GetXaxis()->SetTitle("Parameter #xi");
    g_1->GetYaxis()->SetTitle("#chi^{2}");
    //g_1->GetYaxis()->SetTitleOffset(1.75);
    //g_1->GetYaxis()->SetMaxDigits(3);
    TGaxis* g_1_tgaxis = (TGaxis*)g_1->GetYaxis();
    g_1_tgaxis->SetMaxDigits(3);
    g_1->SetLineColor(4);

    TGraph *g_2 = new TGraph(data_size_2, data_x_2, data_y_2);
    g_2->SetTitle("");
    g_2->GetXaxis()->SetTitle("Parameter #xi");
    g_2->GetYaxis()->SetTitle("#chi^{2}");
    //g_2->GetYaxis()->SetTitleOffset(1.75);
    //g_2->GetYaxis()->SetMaxDigits(3);
    TGaxis* g_2_tgaxis = (TGaxis*)g_2->GetYaxis();
    g_2_tgaxis->SetMaxDigits(3);
    g_2->SetLineColor(6);

    TGraph *g_3 = new TGraph(data_size_3, data_x_3, data_y_3);
    g_3->SetTitle("");
    g_3->GetXaxis()->SetTitle("Parameter #xi");
    g_3->GetYaxis()->SetTitle("#chi^{2}");
    //g_3->GetYaxis()->SetTitleOffset(1.75);
    //((TGAxis*)(g_3->GetYaxis()))->SetMaxDigits(3);
    TGaxis* g_3_tgaxis = (TGaxis*)g_3->GetYaxis();
    g_3_tgaxis->SetMaxDigits(3);
    g_3->SetLineColor(3);

    std::cout << "Finished constructing intermediate data" << std::endl;

    TLegend *l = new TLegend(0.1, 0.7, 0.3, 0.9);
    l->AddEntry(g_1, "Full range cut", "l");
    l->AddEntry(g_0, "Full range no cut", "l");
    l->AddEntry(g_3, "Sub range cut", "l");
    l->AddEntry(g_2, "Sub range no cut", "l");

    ////////////////////////////////////////////////////////////////////////////
    // CREATE OUTPUT GRAPH
    ////////////////////////////////////////////////////////////////////////////
    
    TCanvas *c = new TCanvas("c", "", 800, 600);
    c->GetPad(0)->SetTicks(1, 2);
    g_1->Draw();
    g_0->Draw("same");
    g_2->Draw("same");
    g_3->Draw("same");
    l->Draw();

    c->SaveAs("c_out.png");
    c->SaveAs("c_out.pdf");
    c->SaveAs("c_out.eps");
    c->SaveAs("c_out.C");
            
    TFile *f_out = new TFile("f_out.root", "RECREATE");
    c->Write();
    g_0->Write();
    g_1->Write();
    g_2->Write();
    g_3->Write();
    f_out->Close();
   
    /*

    // enable/disable printing with log canvas
    Bool_t log_mode{true};

    if(batch_mode == false)
    {
        //HistogramWrapper2 hwrap_el_energy("el_energy", "x_axis_label_text", "y_axis_label_text", "E");
        //hwrap_el_energy.SetLogMode(log_mode);
        //hwrap_el_energy.Add("h_el_energy_original", h_el_energy_original);
        //hwrap_el_energy.Canvas();
        // TODO: make some nice histograms
        // X/Y AXIS: LABEL: TEXT, FONT SIZE, FONT
        // X/Y AXIS: TICK NUMBERS: FONT SIZE, FONT
        // X/Y AXIS LIMITS, RANGE, PARTICULARLY Y MAX
        // GRAPH COLORS


        // print single electron distribution
        //TCanvas *c_el_energy_both = new TCanvas("e_el_energy_both", "e_el_energy_both", 800, 600);
        //c_el_energy_both->SetLogy(log_mode);
        //h_el_energy_original->SetMaximum(220.0e3); // MC data
        //h_el_energy_original->SetMaximum(1000.0e3); // MC data log mode
        //h_el_energy_original->GetXaxis()->SetTitle("Energy [MeV]");
        //h_el_energy_original->GetYaxis()->SetTitle("Events");
        //h_el_energy_reweight->SetMaximum(1000.0e3);
        //h_el_energy_reweight->GetXaxis()->SetTitle("Energy [MeV]");
        //h_el_energy_reweight->GetYaxis()->SetTitle("Events");
        //h_el_energy_original->Draw("E");
        //h_el_energy_reweight->Draw("Esame");
        //c_el_energy_both->SaveAs("c_el_energy_both.C");
        //c_el_energy_both->SaveAs("c_el_energy_both.png");
        //c_el_energy_both->SaveAs("c_el_energy_both.pdf");
        //delete c_el_energy_both;
        const Double_t max_log_mode{1000.0e3};
        const Double_t max_nolog_mode{220.0e3};
        Double_t max{0.0};
        Double_t min{0.0};
        const std::string dir_log_mode("el_energy_log_dir");
        const std::string dir_nolog_mode("el_energy_nolog_dir");
        std::string dir;
        if(log_mode) { min = 0.1; max = max_log_mode; dir = dir_log_mode; }
        else { min = 0.0; max = max_nolog_mode; dir = dir_nolog_mode; }
        CanvasFactorySettings settings("Energy [MeV]", "Events", min, max, log_mode);
        settings.SetDrawOption("E");
        CanvasFactory factory(settings);
        factory.Canvas("el_energy", dir, h_el_energy_original, "Baseline", h_el_energy_reweight, "Reweighted");


        // print summed distribution
        //TCanvas *c_el_energy_sum_both = new TCanvas("e_el_energy_sum_both", "e_el_energy_sum_both", 800, 600);
        //c_el_energy_sum_both->SetLogy(log_mode);
        //h_el_energy_sum_original->SetMaximum(80.0e3); // MC data
        //h_el_energy_sum_original->SetMaximum(100.0e3); // MC data log mode
        //h_el_energy_sum_original->GetXaxis()->SetTitle("Energy [MeV]");
        //h_el_energy_sum_original->GetYaxis()->SetTitle("Events");
        //h_el_energy_sum_original->Draw("E");
        //h_el_energy_sum_reweight->Draw("Esame");
        //h_el_energy_sum_original_clone->Draw("Esame");
        //h_el_energy_sum_reweight_clone->Draw("Esame");
        //c_el_energy_sum_both->SaveAs("c_el_energy_sum_both.C");
        //c_el_energy_sum_both->SaveAs("c_el_energy_sum_both.png");
        //c_el_energy_sum_both->SaveAs("c_el_energy_sum_both.pdf");
        //delete c_el_energy_sum_both; 
        const Double_t max_log_mode_2{100.0e3};
        const Double_t max_nolog_mode_2{80.0e3};
        const std::string dir_log_mode_2("el_energy_sum_log_dir");
        const std::string dir_nolog_mode_2("el_energy_sum_nolog_dir");
        if(log_mode) { min = 0.1; max = max_log_mode_2; dir = dir_log_mode_2; }
        else { min = 0.0; max = max_nolog_mode_2; dir = dir_nolog_mode_2; }
        CanvasFactorySettings settings_2("Energy [MeV]", "Events", min, max, log_mode);
        settings_2.SetDrawOption("E");
        CanvasFactory factory_2(settings_2);
        factory_2.Canvas("el_energy_sum", dir, h_el_energy_sum_original, "Baseline", h_el_energy_sum_reweight, "Reweighted");
        // TODO: settings should be passed to canvas in case settings should be changed?
        // or setsettings function should be provided


        // print single electron distribution test histograms
        TCanvas *c_test_single = new TCanvas("c_test_single", "c_test_single", 800, 600);
        c_test_single->SetLogy(log_mode);
        //h_test_single_original->SetMaximum(1.4); // MC data
        //h_test_single_original->SetMaximum(1.2); // createtree data
        h_test_single_original->SetMaximum(10.0); // createtree data log mode
        h_test_single_original->GetXaxis()->SetTitle("Energy [MeV]");
        h_test_single_original->GetYaxis()->SetTitle("Events");
        h_test_single_original->Draw("E");
        h_test_single_reweight->Draw("Esame");
        g_el_energy_single_0->Draw("same");
        g_el_energy_single_1->Draw("same");
        g_el_energy_single_2->Draw("same");
        //h_test_single_reweight->Draw("E");
        c_test_single->SaveAs("c_test_single.C");
        c_test_single->SaveAs("c_test_single.png");
        c_test_single->SaveAs("c_test_single.pdf");
        delete c_test_single;
        // compute chi-square
        std::cout << "chi1: " << chi_square_test(g_el_energy_single_0, h_test_single_original) / h_test_single_original->GetNbinsX() << std::endl;
        std::cout << "chi2: " << chi_square_test(g_el_energy_single_2, h_test_single_reweight) / h_test_single_reweight->GetNbinsX() << std::endl;
        
        // print summed distribution test histograms
        TCanvas *c_test_sum = new TCanvas("c_test_sum", "c_test_sum", 800, 600);
        c_test_sum->SetLogy(log_mode);
        //h_test_sum_original->SetMaximum(1.0); // MC data
        //h_test_sum_original->SetMaximum(0.8); // createtree data
        h_test_sum_original->SetMaximum(1.0); // createtree data log mode
        h_test_sum_original->GetXaxis()->SetTitle("Energy [MeV]");
        h_test_sum_original->GetYaxis()->SetTitle("Events");
        h_test_sum_original->Draw("E");
        h_test_sum_reweight->Draw("Esame");
        g_el_energy_sum_0->Draw("same");
        g_el_energy_sum_1->Draw("same");
        g_el_energy_sum_2->Draw("same");
        //h_test_sum_reweight->Draw("E");
        c_test_sum->SaveAs("c_test_sum.C");
        c_test_sum->SaveAs("c_test_sum.png");
        c_test_sum->SaveAs("c_test_sum.pdf");
        delete c_test_sum;
        // compute chi-square
        std::cout << "chi3: " << chi_square_test(g_el_energy_sum_0, h_test_sum_original) / h_test_sum_original->GetNbinsX() << std::endl;
        std::cout << "chi4: " << chi_square_test(g_el_energy_sum_2, h_test_sum_reweight) / h_test_sum_reweight->GetNbinsX() << std::endl;

    }

    */

    return 0;

}
