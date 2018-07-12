#include "Analysis.hpp"


#include "aux.hpp"


////////////////////////////////////////////////////////////////////////////////
// 2d histogram chisquare sensitivity
////////////////////////////////////////////////////////////////////////////////

void Analysis::SensitivityMeasurementChisquare2()
{

    const std::string eps_string{std::to_string(epsilon_31)};
    const std::string systematic_energy_mult_string{std::to_string(systematic_energy_mult)};

    


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
