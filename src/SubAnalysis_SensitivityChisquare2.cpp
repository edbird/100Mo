#include "Analysis.hpp"


#include "aux.hpp"


////////////////////////////////////////////////////////////////////////////////
// 2d histogram chisquare sensitivity
////////////////////////////////////////////////////////////////////////////////

void SubAnalysis::SensitivityMeasurementChisquare2()
{

    const std::string eps_string{std::to_string(epsilon_31)};
    const std::string systematic_energy_mult_string{std::to_string(systematic_energy_mult)};

    // reset
    non_empty_bins_2d = 0;
    
    ////////////////////////////////////////////////////////////////////////////
    // INDEPENDENT (2D) ELECTRON ENERGY CHISQUARE METHOD
    ////////////////////////////////////////////////////////////////////////////

    // get chi-square for 2d electron histogram
    if(fit_subrange == false)
    {

        for(Int_t j{1}; j <= h_el_energy_2d_original->GetNbinsY(); ++ j)
        {
            for(Int_t i{1}; i <= h_el_energy_2d_original->GetNbinsX(); ++ i)
            {
                if(h_el_energy_2d_original->GetBinContent(i, j) != 0.0) ++ non_empty_bins_2d;
            }
        }

        sensitivity_chisquare_2d = chi_square_test(h_el_energy_2d_reweight, h_el_energy_2d_original);
        std::cout << "chi square of 2 electron: " << sensitivity_chisquare_2d << std::endl;
        std::cout << " degrees of freedom (2d): " << non_empty_bins_2d << std::endl;
        std::cout << " chi square reduced: " << sensitivity_chisquare_2d / (Double_t)non_empty_bins_2d << std::endl;

    }
    else if(fit_subrange == true)
    {

        /*
        for(Int_t j{1}; j <= h_el_energy_2d_original->GetNbinsY(); ++ j)
        {
            for(Int_t i{1}; i <= h_el_energy_2d_original->GetNbinsX(); ++ i)
            {
                if(h_el_energy_2d_original->GetYaxis()->GetBinCenter(j) >= 2.0)
                {
                    if(h_el_energy_2d_original->GetXaxis()->GetBinCenter(i) >= 2.0)
                    {
                        if(h_el_energy_2d_original->GetBinContent(i, j) != 0.0) ++ non_empty_bins_2d;
                    }
                }
            }
        }
        */

        sensitivity_chisquare_2d = chi_square_test(h_el_energy_2d_reweight, h_el_energy_2d_original, 2.0, 4.0, non_empty_bins_2d);
        std::cout << "chi square of 2 electron, 2.0 MeV - 2.0 MeV: " << sensitivity_chisquare_2d << std::endl;
        std::cout << " degrees of freedom (2d): " << non_empty_bins_2d << std::endl;
        std::cout << " chi square reduced: " << sensitivity_chisquare_2d / (Double_t)non_empty_bins_2d << std::endl;
        
    }

    // NOTE TO SELF: There is no 2d fit, the 2d histogram is used to evaluate
    // the sensitivity, therefore there is only a chi-square test


    ////////////////////////////////////////////////////////////////////////
    // CREATE DIFFERENCE AND PULL HISTOGRAMS AND CANVAS OUTPUT
    ////////////////////////////////////////////////////////////////////////

    // create difference histogram
    h_el_energy_2d_diff = new TH2D("h_el_energy_2d_diff", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_el_energy_2d_diff->SetStats(0);
    for(Int_t j{1}; j <= h_el_energy_2d_diff->GetNbinsY(); ++ j)
    {
        for(Int_t i{1}; i <= h_el_energy_2d_diff->GetNbinsX(); ++ i)
        {
            Double_t content1{h_el_energy_2d_original->GetBinContent(i, j)};
            Double_t content2{h_el_energy_2d_reweight->GetBinContent(i, j)};
            Double_t error1{h_el_energy_2d_original->GetBinError(i, j)};
            // note: do not use error or reweighted
            Double_t error2{0.0 * h_el_energy_2d_reweight->GetBinError(i, j)};
            Double_t content{content1 - content2};
            Double_t error{std::sqrt(error1 * error1 + error2 * error2)};
            h_el_energy_2d_diff->SetBinContent(i, j, content);
            h_el_energy_2d_diff->SetBinError(i, j, error);
        }
    }
        
    // create pull histogram
    h_el_energy_2d_pull = new TH2D("h_el_energy_2d_pull", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_el_energy_2d_pull->SetStats(0);
    for(Int_t j{1}; j <= h_el_energy_2d_pull->GetNbinsY(); ++ j)
    {
        for(Int_t i{1}; i <= h_el_energy_2d_pull->GetNbinsX(); ++ i)
        {
            Double_t content{h_el_energy_2d_diff->GetBinContent(i, j)};
            Double_t error{h_el_energy_2d_diff->GetBinError(i, j)};
            if(error == 0.0)
            {
                //h_el_energy_2d_pull->SetBinContent(i, j, 0.0);
                h_el_energy_2d_pull->SetBinError(i, j, 0.0);
            }
            else
            {
                h_el_energy_2d_pull->SetBinContent(i, j, content / error);
                h_el_energy_pull->SetBinError(i, j, 0.0);
            }
        }
    }

    // create ratio histogram
    std::string h_el_energy_2d_ratio_name{std::string("h_el_energy_2d_ratio") + std::string("_") + eps_string + std::string("_") + systematic_energy_mult_string};
    h_el_energy_2d_ratio = new TH2D(h_el_energy_2d_ratio_name.c_str(), "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_el_energy_2d_ratio->SetStats(0);
    for(Int_t j{1}; j <= h_el_energy_2d_ratio->GetNbinsY(); ++ j)
    {
        for(Int_t i{1}; i <= h_el_energy_2d_ratio->GetNbinsX(); ++ i)
        {
            Double_t content1{h_el_energy_2d_original->GetBinContent(i, j)};
            Double_t content2{h_el_energy_2d_reweight->GetBinContent(i, j)};
            Double_t error1{h_el_energy_2d_original->GetBinError(i, j)}; // correct way round?
            // note: do not use error or reweighted
            Double_t error2{0.0 * h_el_energy_2d_reweight->GetBinError(i, j)};
            if(content2 != 0.0)
            {
                Double_t content{content1 / content2}; // correct way round?
                Double_t error{std::sqrt(std::pow(error1 * (1.0 / content2), 2.0) + std::pow(error2 * (content1 / (content2 * content2)), 2.0))}; // should error2 be used?
                h_el_energy_2d_ratio->SetBinContent(i, j, content);
                h_el_energy_2d_ratio->SetBinError(i, j, error);
                // what to do if content1 = 0.0 -> div by zero error
            }
        }
    }
    std::string f_el_energy_2d_ratio_name{std::string("h_el_energy_2d_ratio") + std::string("_") + eps_string + std::string(".root")};
    TFile *f_el_energy_2d_ratio = new TFile(f_el_energy_2d_ratio_name.c_str(), "UPDATE");
    h_el_energy_2d_ratio->Write();
    f_el_energy_2d_ratio->Close();
    delete f_el_energy_2d_ratio;
    f_el_energy_2d_ratio = nullptr;


    ////////////////////////////////////////////////////////////////////////////
    // 2D ORIGINAL AND REWEIGHT CANVAS OUTPUT
    ////////////////////////////////////////////////////////////////////////////
        
    #if 0
    if(_batch_mode_ == false)
    {

        // 2d original and reweight histograms
        c_el_energy_2d_original = new TCanvas("c_el_energy_2d_original", "", 800, 600);
        c_el_energy_2d_original->SetRightMargin(0.17);
        //c_el_energy_2d_original->SetLogz();
        h_el_energy_2d_original->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
        h_el_energy_2d_original->GetYaxis()->SetTitle("High Energy Electron [MeV]");
        h_el_energy_2d_original->GetZaxis()->SetTitle("Events");
        h_el_energy_2d_original->GetZaxis()->SetTitleOffset(1.6);
        h_el_energy_2d_original->Draw("colz");
        c_el_energy_2d_original->SaveAs("c_el_energy_2d_original.C");
        c_el_energy_2d_original->SaveAs("c_el_energy_2d_original.png");
        c_el_energy_2d_original->SaveAs("c_el_energy_2d_original.pdf");

        c_el_energy_2d_reweight = new TCanvas("c_el_energy_2d_reweight", "", 800, 600);
        c_el_energy_2d_reweight->SetRightMargin(0.17);
        //c_el_energy_2d_reweight->SetLogz();
        h_el_energy_2d_reweight->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
        h_el_energy_2d_reweight->GetYaxis()->SetTitle("High Energy Electron [MeV]");
        h_el_energy_2d_reweight->GetZaxis()->SetTitle("Events");
        h_el_energy_2d_reweight->GetZaxis()->SetTitleOffset(1.6);
        h_el_energy_2d_reweight->Draw("colz");
        c_el_energy_2d_reweight->SaveAs("c_el_energy_2d_reweight.C");
        c_el_energy_2d_reweight->SaveAs("c_el_energy_2d_reweight.png");
        c_el_energy_2d_reweight->SaveAs("c_el_energy_2d_reweight.pdf");
        
        c_el_energy_2d_diff = new TCanvas("c_el_energy_2d_diff", "", 800, 600);
        c_el_energy_2d_diff->SetRightMargin(0.15);
        //c_el_energy_2d_diff->SetLogz();
        h_el_energy_2d_diff->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
        h_el_energy_2d_diff->GetYaxis()->SetTitle("High Energy Electron [MeV]");
        h_el_energy_2d_diff->GetZaxis()->SetTitle("Events");
        h_el_energy_2d_diff->GetZaxis()->SetTitleOffset(1.2);
        h_el_energy_2d_diff->Draw("colz");
        c_el_energy_2d_diff->SaveAs("c_el_energy_2d_diff.C");
        c_el_energy_2d_diff->SaveAs("c_el_energy_2d_diff.png");
        c_el_energy_2d_diff->SaveAs("c_el_energy_2d_diff.pdf");
        
        c_el_energy_2d_pull = new TCanvas("c_el_energy_2d_pull", "", 800, 600);
        c_el_energy_2d_pull->SetRightMargin(0.15);
        //c_el_energy_2d_pull->SetLogz();
        h_el_energy_2d_pull->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
        h_el_energy_2d_pull->GetYaxis()->SetTitle("High Energy Electron [MeV]");
        h_el_energy_2d_pull->GetZaxis()->SetTitle("Events");
        h_el_energy_2d_pull->GetZaxis()->SetTitleOffset(1.2);
        h_el_energy_2d_pull->Draw("colz");
        c_el_energy_2d_pull->SaveAs("c_el_energy_2d_pull.C");
        c_el_energy_2d_pull->SaveAs("c_el_energy_2d_pull.png");
        c_el_energy_2d_pull->SaveAs("c_el_energy_2d_pull.pdf");

        ////////////////////////////////////////////////////////////////////////
        // SAVE TO ROOT FILE
        ////////////////////////////////////////////////////////////////////////

        TFile *f_el_energy_2d_chisquare = new TFile("f_el_energy_2d_chisquare.root", "RECREATE");
        h_el_energy_2d_original->Write();
        c_el_energy_2d_original->Write();
        h_el_energy_2d_reweight->Write();
        c_el_energy_2d_reweight->Write();
        h_el_energy_2d_diff->Write();
        c_el_energy_2d_diff->Write();
        h_el_energy_2d_pull->Write();
        c_el_energy_2d_pull->Write();
        f_el_energy_2d_chisquare->Close();
        //delete f_el_energy_2d;


        // double check chi-squre result
        // result: OK
        /*
        Double_t chisq{0.0};
        for(Int_t j{1}; j <= h_el_energy_2d_pull->GetNbinsY(); ++ j)
        {
            for(Int_t i{1}; i <= h_el_energy_2d_pull->GetNbinsX(); ++ i)
            {
                Double_t content{h_el_energy_2d_pull->GetBinContent(i, j)};
                chisq += content * content;
            }
        }
        std::cout << chisq << std::endl;
        */

    }
    #endif
 
}
