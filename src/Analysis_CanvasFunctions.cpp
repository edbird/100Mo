// self header
#include "Analysis.hpp"



////////////////////////////////////////////////////////////////////////////////
// CREATE 2D DECAY RATE HISTOGRAMS
////////////////////////////////////////////////////////////////////////////////

void Analysis::CanvasDecayRate()
{

    ////////////////////////////////////////////////////////////////////////////
    // CREATE OUTPUT HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////

    //Int_t dimension_xy{1001};
    // don't plot raw data
    //TH2D *h_nEqNull = new TH2D("h_nEqNull", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_nEqTwo = new TH2D("h_nEqTwo", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    h_data_0 = new TH2D("h_data_0", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    h_data_1 = new TH2D("h_data_1", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    h_data_2 = new TH2D("h_data_2", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_ratio = new TH2D("h_ratio", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //h_nEqNull->SetStats(0);
    //h_nEqTwo->SetStats(0);
    h_data_0->SetStats(0);
    h_data_1->SetStats(0);
    h_data_2->SetStats(0);
    //h_ratio->SetStats(0);

    for(std::size_t i{0}; i < dimension_xy; ++ i)
    {
        for(std::size_t j{0}; j < dimension_xy; ++ j)
        {
            //h_nEqNull->SetBinContent(i, j, data_nEqNull.at(i * dimension_xy + j)[2]);
            //h_nEqTwo->SetBinContent(i, j, data_nEqTwo.at(i * dimension_xy + j)[2]);
            h_data_0->SetBinContent(i + 1, j + 1, data_0.at(i * dimension_xy + j)[2]);
            h_data_1->SetBinContent(i + 1, j + 1, data_1.at(i * dimension_xy + j)[2]);
            h_data_2->SetBinContent(i + 1, j + 1, data_2.at(i * dimension_xy + j)[2]);
            //if(i < dimension_xy - j)
            // TODO: -1 or -0 ?
            if(i + j < dimension_xy - 1)
            {
                // TODO: move above lines to inside this if
                //h_ratio->SetBinContent(i, j, ratio.at(i * dimension_xy + j)[2]);
            }
        }
    }
    std::cout << "Finished constructing histograms" << std::endl;

    // used for createtree input data
    // write data histograms to file
    // intermediate file
    // used for creating test histogram
    // with createtree code
    TFile *f_histogram = new TFile("f_histogram.root", "recreate");
    h_data_0->Write();
    h_data_1->Write();
    h_data_2->Write();
    f_histogram->Close();

    
    if(_batch_mode_ == false)
    {
        if(_canvas_enable_raw_data_)
        {
            c_nEqNull = new TCanvas("c_nEqNull", "", 4000, 3000);
            h_nEqNull->Draw("colz");
            c_nEqNull->SaveAs("c_nEqNull.png");
            c_nEqNull->SaveAs("c_nEqNull.pdf");
            c_nEqNull->SaveAs("c_nEqNull.C");

            c_nEqTwo = new TCanvas("c_nEqTwo", "", 4000, 3000);
            h_nEqTwo->Draw("colz");
            c_nEqTwo->SaveAs("c_nEqTwo.png");
            c_nEqTwo->SaveAs("c_nEqTwo.pdf");
            c_nEqTwo->SaveAs("c_nEqTwo.C");
        }


        if(_canvas_enable_decay_rate_)
        {
            c_data_0 = new TCanvas("c_data_0", "", 4000, 3000);
            h_data_0->Draw("colz");
            c_data_0->SaveAs("c_data_0.png");
            c_data_0->SaveAs("c_data_0.pdf");
            c_data_0->SaveAs("c_data_0.C");

            c_data_1 = new TCanvas("c_data_1", "", 4000, 3000);
            h_data_1->Draw("colz");
            c_data_1->SaveAs("c_data_1.png");
            c_data_1->SaveAs("c_data_1.pdf");
            c_data_1->SaveAs("c_data_1.C");

            c_data_2 = new TCanvas("c_data_2", "", 4000, 3000);
            h_data_2->Draw("colz");
            c_data_2->SaveAs("c_data_2.png");
            c_data_2->SaveAs("c_data_2.pdf");
            c_data_2->SaveAs("c_data_2.C");
        }
    }

    /*
    TCanvas *c_ratio = new TCanvas("c_ratio", "", 4000, 3000);
    h_ratio->Draw("colz");
    c_ratio->SaveAs("c_ratio.png");
    c_ratio->SaveAs("c_ratio.pdf");
    c_ratio->SaveAs("c_ratio.C");
    delete c_ratio;
    */
}


////////////////////////////////////////////////////////////////////////////////
// CREATE SINGLE ELECTRON CANVAS
////////////////////////////////////////////////////////////////////////////////

void Analysis::CanvasSingleElectronProjection()
{

    // epsilon_31 = 0
    h_single_electron_0 = h_data_0->ProjectionX("h_single_electron_0");
    h_single_electron_0->Scale(1.0 / (Double_t)dimension_xy);
    /*for(Int_t i{1}; i <= h_single_electron_0->GetNbinsX(); ++ i)
    {
        Double_t F{1.0 / (Double_t)dimension_xy};
        h_single_electron_0->SetBinContent(i, F * h_single_electron_0->GetBinContent(i));
    }*/
    h_single_electron_0->SetStats(0);
    h_single_electron_0->SetLineColor(2);
    h_single_electron_0->SetMarkerColor(2);
    
    // epsilon_31 = 0.4
    h_single_electron_1 = h_data_1->ProjectionX("h_single_electron_1");
    h_single_electron_1->Scale(1.0 / (Double_t)dimension_xy);
    /*for(Int_t i{1}; i <= h_single_electron_1->GetNbinsX(); ++ i)
    {
        Double_t F{1.0 / (Double_t)dimension_xy};
        h_single_electron_1->SetBinContent(i, F * h_single_electron_1->GetBinContent(i));
    }*/
    h_single_electron_1->SetStats(0);
    h_single_electron_1->SetLineColor(3);
    h_single_electron_1->SetMarkerColor(3);
    
    // epsilon_31 = 0.8
    h_single_electron_2 = h_data_2->ProjectionX("h_single_electron_2");
    h_single_electron_2->Scale(1.0 / (Double_t)dimension_xy);
    /*for(Int_t i{1}; i <= h_single_electron_2->GetNbinsX(); ++ i)
    {
        Double_t F{1.0 / (Double_t)dimension_xy};
        h_single_electron_2->SetBinContent(i, F * h_single_electron_2->GetBinContent(i));
    }*/
    h_single_electron_2->SetStats(0);
    h_single_electron_2->SetLineColor(4);
    h_single_electron_2->SetMarkerColor(4);

    /*
    TH1D *h_single_electron_0 = new TH1D("h_single_electron_0", "", dimension_xy, 0.0, 1.0);
    h_single_electron_0->SetStats(0);

    for(Int_t i{1}; i <= h_ratio->GetNbinsX(); ++ i)
    {
        Double_t content{0.0};
        for(Int_t j{1}; j <= h_ratio->GetNbinsY(); ++ j)
        {
            content += h_ratio->GetBinContent(i, j);
        }
        h_single_electron_0->SetBinContent(i, content);
    }
    */


    // these histograms created from projection of the h_data_? histograms
    // along the X direction
    if(_batch_mode_ == false)
    {
        if(_canvas_enable_single_electron_projection_)
        {
            c_single_electron = new TCanvas("c_single_electron", "", 4000, 3000);
            h_single_electron_0->SetMaximum(3.5);
            h_single_electron_0->GetXaxis()->SetTitle("Energy [MeV]");
            h_single_electron_0->GetYaxis()->SetTitle("Events");
            h_single_electron_0->Draw("hist");
            h_single_electron_1->Draw("histsame");
            h_single_electron_2->Draw("histsame");
            c_single_electron->SaveAs("c_single_electron.png");
            c_single_electron->SaveAs("c_single_electron.pdf");
            c_single_electron->SaveAs("c_single_electron.C");
        }
    }


    // same as the single electron histogram,
    // however it is made using an explicit integration method
    // results appear to be the same
    /*
    TH1D *h_single_electron_test = new TH1D("h_single_electron_test", "", dimension_xy, 0.0, 1.0);
    h_single_electron_test->SetStats(0);
    for(Int_t i{1}; i <= h_data_0->GetNbinsX(); ++ i)
    {
        Double_t energy_1{h_data_0->GetXaxis()->GetBinCenter(i)};
        Double_t weight_integral_1{0.0};
        Double_t weight_integral_2{0.0};
        for(Int_t j{1}; j <= h_data_0->GetNbinsY(); ++ j)
        {
            Double_t energy_2{h_data_0->GetXaxis()->GetBinCenter(j)};
            if(energy_1 + energy_2 <= 1.0)
            {
                weight_integral_1 += h_data_0->GetBinContent(i, j);
                weight_integral_2 += h_data_2->GetBinContent(i, j);
            }
        }
        if(weight_integral_1 > 0.0)
        {
            Double_t weight_factor{weight_integral_2 / weight_integral_1};
            Double_t input_weight{h_single_electron_0->GetBinContent(i)};
            Double_t output_weight{input_weight * weight_factor};
            h_single_electron_test->SetBinContent(i, output_weight);
        }
    }
    
 
    // test output, with histogram created from reweighting
    TCanvas *c_single_electron_test = new TCanvas("c_single_electron_test", "", 4000, 3000);
        h_single_electron_0->GetYaxis()->SetMaximum();
        h_single_electron_0->GetXaxis()->SetTitle("Energy [MeV]");
        h_single_electron_0->GetYaxis()->SetTitle("Events");
    h_single_electron_test->Draw("hist");
    c_single_electron_test->SaveAs("c_single_electron_test.png");
    c_single_electron_test->SaveAs("c_single_electron_test.pdf");
    c_single_electron_test->SaveAs("c_single_electron_test.C");
    delete c_single_electron_test;
    */

    // compare above histogram (created from reweighting)
    // with histograms created from profile
    // NOTE: removed as this code is exactly the same as c_single_electron
    /*
    TCanvas *c_single_electron_with_test = new TCanvas("c_single_electron_with_test", "", 4000, 3000);
        h_single_electron_0->GetYaxis()->SetMaximum();
        h_single_electron_0->GetXaxis()->SetTitle("Energy [MeV]");
        h_single_electron_0->GetYaxis()->SetTitle("Events");
    h_single_electron_0->Draw("hist");
    h_single_electron_1->Draw("histsame");
    h_single_electron_2->Draw("histsame");
    h_single_electron_test->Draw("histsame");
    c_single_electron_with_test->SaveAs("c_single_electron_with_test.png");
    c_single_electron_with_test->SaveAs("c_single_electron_with_test.pdf");
    c_single_electron_with_test->SaveAs("c_single_electron_with_test.C");
    delete c_single_electron_with_test;
    */
}


////////////////////////////////////////////////////////////////////////////////
// THEORY DATA COMPARISON CANVAS OUTPUT
////////////////////////////////////////////////////////////////////////////////

void Analysis::InitSingleElectronTest()
{

    // test histograms
    // NEMO-3 data/MC: truth electron energy x2 (T1, T2) (2 in same histo)
    h_test_single_original = new TH1D("h_test_single_original", "", num_bins, 0.0, 4.0);
    h_test_single_original->SetStats(0);
    h_test_single_original->SetLineColor(2);
    h_test_single_original->SetMarkerColor(2);
    // re-weighted
    h_test_single_reweight = new TH1D("h_test_single_reweight", "", num_bins, 0.0, 4.0);
    h_test_single_reweight->SetStats(0);
    h_test_single_reweight->SetLineColor(4);
    h_test_single_reweight->SetMarkerColor(4);

    // test histograms
    // NEMO-3 data/MC: reconstructed electron energy x2 (2 in same histo)
    h_test_sum_original = new TH1D("h_test_sum_original", "", num_bins, 0.0, 4.0);
    h_test_sum_original->SetStats(0);
    h_test_sum_original->SetLineColor(2);
    h_test_sum_original->SetMarkerColor(2);
    // re-weighted
    h_test_sum_reweight = new TH1D("h_test_sum_reweight", "", num_bins, 0.0, 4.0);
    h_test_sum_reweight->SetStats(0);
    h_test_sum_reweight->SetLineColor(4);
    h_test_sum_reweight->SetMarkerColor(4);

}


void Analysis::CanvasSingleElectronTest()
{

    // create file with test histogerams and graph data
    TFile *f_test_histos = new TFile("f_test_histos.root", "RECREATE");
    h_test_single_original->Write();
    h_test_single_reweight->Write();
    h_test_sum_original->Write();
    h_test_sum_reweight->Write();
    g_el_energy_single_0->SetName("g_el_energy_single_0");
    g_el_energy_single_1->SetName("g_el_energy_single_1");
    g_el_energy_single_2->SetName("g_el_energy_single_2");
    g_el_energy_single_0->Write();
    g_el_energy_single_1->Write();
    g_el_energy_single_2->Write();
    g_el_energy_sum_0->SetName("g_el_energy_sum_0");
    g_el_energy_sum_1->SetName("g_el_energy_sum_1");
    g_el_energy_sum_2->SetName("g_el_energy_sum_2");
    g_el_energy_sum_0->Write();
    g_el_energy_sum_1->Write();
    g_el_energy_sum_2->Write();
    f_test_histos->Close();


    // NOTE: development of this section has been moved to a different file
    // (program)
    if(_batch_mode_ == false)
    {

        TLegend *l_test_single = new TLegend(0.65, 0.63, 0.85, 0.78);
        l_test_single->SetBorderSize(0);
        l_test_single->AddEntry(g_el_energy_single_0, " #xi = 0.0", "l");
        l_test_single->AddEntry(g_el_energy_single_1, " #xi = 0.4", "l");
        l_test_single->AddEntry(g_el_energy_single_2, " #xi = 0.8", "l");

        // print single electron distribution test histograms
        c_test_single = new TCanvas("c_test_single", "c_test_single", 800, 600);
        c_test_single->SetLogy(log_mode);
        //h_test_single_original->SetMaximum(1.4); // MC data
        h_test_single_original->SetMaximum(1.2); // createtree data
        //h_test_single_original->SetMaximum(10.0); // createtree data log mode
        h_test_single_original->GetXaxis()->SetTitle("Single Electron Energy [MeV]");
        h_test_single_original->GetYaxis()->SetTitle("Events");
        h_test_single_original->Draw("E");
        h_test_single_reweight->Draw("Esame");
        g_el_energy_single_0->Draw("Lsame");
        g_el_energy_single_1->Draw("Lsame");
        g_el_energy_single_2->Draw("Lsame");
        //h_test_single_reweight->Draw("E");
        l_test_single->Draw();
        c_test_single->SaveAs("c_test_single.C");
        c_test_single->SaveAs("c_test_single.png");
        c_test_single->SaveAs("c_test_single.pdf");
        c_test_single->SaveAs("c_test_single.eps");
        // compute chi-square
        //std::cout << "chi1: " << chi_square_test(g_el_energy_single_0, h_test_single_original) / h_test_single_original->GetNbinsX() << std::endl;
        //std::cout << "chi2: " << chi_square_test(g_el_energy_single_2, h_test_single_reweight) / h_test_single_reweight->GetNbinsX() << std::endl;
    }

    if(_batch_mode_ == false)
    {

        TLegend *l_test_sum = new TLegend(0.65, 0.63, 0.85, 0.78);
        l_test_sum->SetBorderSize(0);
        l_test_sum->AddEntry(g_el_energy_sum_0, " #xi = 0.0", "l");
        l_test_sum->AddEntry(g_el_energy_sum_1, " #xi = 0.4", "l");
        l_test_sum->AddEntry(g_el_energy_sum_2, " #xi = 0.8", "l");

        // print summed distribution test histograms
        c_test_sum = new TCanvas("c_test_sum", "c_test_sum", 800, 600);
        //c_test_sum->SetLogy(log_mode);
        //h_test_sum_original->SetMaximum(1.0); // MC data
        h_test_sum_original->SetMaximum(0.8); // createtree data
        //h_test_sum_original->SetMaximum(1.0); // createtree data log mode
        h_test_sum_original->GetXaxis()->SetTitle("Sum Electron Energy [MeV]");
        h_test_sum_original->GetYaxis()->SetTitle("Events");
        h_test_sum_original->Draw("E");
        h_test_sum_reweight->Draw("Esame");
        //g_el_energy_sum_0->Draw("AL");
        g_el_energy_sum_0->Draw("Lsame");
        g_el_energy_sum_1->Draw("Lsame");
        g_el_energy_sum_2->Draw("Lsame");
        //h_test_sum_reweight->Draw("E");
        l_test_sum->Draw();
        c_test_sum->SaveAs("c_test_sum.C");
        c_test_sum->SaveAs("c_test_sum.png");
        c_test_sum->SaveAs("c_test_sum.pdf");
        c_test_sum->SaveAs("c_test_sum.eps");
        // compute chi-square
        //std::cout << "chi3: " << chi_square_test(g_el_energy_sum_0, h_test_sum_original) / h_test_sum_original->GetNbinsX() << std::endl;
        //std::cout << "chi4: " << chi_square_test(g_el_energy_sum_2, h_test_sum_reweight) / h_test_sum_reweight->GetNbinsX() << std::endl;

    }
}
