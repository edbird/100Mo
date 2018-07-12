
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"


#include <iostream>


int main()
{


    TFile *f = new TFile("f_test_histos.root", "READ");

    TH1D* h_test_single_original = (TH1D*)f->Get("h_test_single_original");
    TH1D* h_test_single_reweight = (TH1D*)f->Get("h_test_single_reweight");
    TH1D* h_test_sum_original = (TH1D*)f->Get("h_test_sum_original");
    TH1D* h_test_sum_reweight = (TH1D*)f->Get("h_test_sum_reweight");
    TGraph* g_el_energy_single_0 = (TGraph*)f->Get("g_el_energy_single_0");
    TGraph* g_el_energy_single_1 = (TGraph*)f->Get("g_el_energy_single_1");
    TGraph* g_el_energy_single_2 = (TGraph*)f->Get("g_el_energy_single_2");
    TGraph* g_el_energy_sum_0 = (TGraph*)f->Get("g_el_energy_sum_0");
    TGraph* g_el_energy_sum_1 = (TGraph*)f->Get("g_el_energy_sum_1");
    TGraph* g_el_energy_sum_2 = (TGraph*)f->Get("g_el_energy_sum_2");
    
    f->Close();

    std::cout << "data read from file" << std::endl;

    TCanvas *c_test_single, *c_test_sum;
    Bool_t log_mode{kFALSE};


    std::cout << "c_test_single" << std::endl;

    // print single electron distribution test histograms
    c_test_single = new TCanvas("c_test_single", "c_test_single", 800, 600);
    std::cout << "new canvas" << std::endl;
    c_test_single->SetLogy(log_mode);
    std::cout << "set log mode" << std::endl;
    std::cout << h_test_single_original << std::endl;
    //h_test_single_original->SetMaximum(1.4); // MC data
    //h_test_single_original->SetMaximum(1.2); // createtree data
    h_test_single_original->SetMaximum(10.0); // createtree data log mode
    std::cout << "set max" << std::endl;
    h_test_single_original->GetXaxis()->SetTitle("Energy [MeV]");
    h_test_single_original->GetYaxis()->SetTitle("Events");
    std::cout << "predraw" << std::endl;
    h_test_single_original->Draw("E");
    std::cout << "original draw" << std::endl;
    h_test_single_reweight->Draw("Esame");
    std::cout << "reweight draw" << std::endl;
    g_el_energy_single_0->Draw("same");
    g_el_energy_single_1->Draw("same");
    g_el_energy_single_2->Draw("same");
    //h_test_single_reweight->Draw("E");
    c_test_single->SaveAs("./canvas_test_output/c_test_single.C");
    c_test_single->SaveAs("./canvas_test_output/c_test_single.png");
    c_test_single->SaveAs("./canvas_test_output/c_test_single.pdf");
    c_test_single->SaveAs("./canvas_test_output/c_test_single.eps");
    // compute chi-square
    //std::cout << "chi1: " << chi_square_test(g_el_energy_single_0, h_test_single_original) / h_test_single_original->GetNbinsX() << std::endl;
    //std::cout << "chi2: " << chi_square_test(g_el_energy_single_2, h_test_single_reweight) / h_test_single_reweight->GetNbinsX() << std::endl;
    

    std::cout << "c_test_sum" << std::endl;

    // print summed distribution test histograms
    c_test_sum = new TCanvas("c_test_sum", "c_test_sum", 800, 600);
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
    c_test_sum->SaveAs("./canvas_test_output/c_test_sum.C");
    c_test_sum->SaveAs("./canvas_test_output/c_test_sum.png");
    c_test_sum->SaveAs("./canvas_test_output/c_test_sum.pdf");
    c_test_sum->SaveAs("./canvas_test_output/c_test_sum.eps");


}
