#include "Analysis.hpp"


#include "aux.hpp"


////////////////////////////////////////////////////////////////////////////////
// 1D SINGLE ELECTRON ENERGY CHISQUARE SENSITIVITY
////////////////////////////////////////////////////////////////////////////////

void Analysis::SensitivityMeasurementChisquare1()
{

    const std::string eps_string{std::to_string(epsilon_31)};
    //const std::string systematic_energy_mult_string{std::to_string(systematic_energy_mult)};
    
    ////////////////////////////////////////////////////////////////////////////
    // SINGLE ELECTRON (1D) CHISQUARE CANVAS OUTPUT
    ////////////////////////////////////////////////////////////////////////////

    /*
    if(_batch_mode_ == false)
    {
    
        // print single electron distribution
        const Double_t max_log_mode{1000.0e3};
        const Double_t max_nolog_mode{220.0e3};
        Double_t max{0.0};
        Double_t min{0.0};
        const std::string dir_log_mode("el_energy_log_dir");
        const std::string dir_nolog_mode("el_energy_nolog_dir");
        std::string dir; //{std::to_string(epsilon_31)};
        if(log_mode) { min = 0.1; max = max_log_mode; dir = dir_log_mode; }
        else { min = 0.0; max = max_nolog_mode; dir = dir_nolog_mode; }
        CanvasFactorySettings settings("Single Electron Energy [MeV]", "Events", min, max, log_mode);
        settings.SetDrawOption("E");
        CanvasFactory factory(settings);
        factory.Canvas("el_energy", dir, h_el_energy_original, "Baseline", h_el_energy_reweight, "Reweighted");

        c_el_energy_diff = new TCanvas("c_el_energy_diff", "", 800, 600);
        //c_el_energy_diff->SetRightMargin(0.12);
        //c_el_energy_diff->SetLogz();
        h_el_energy_diff->GetXaxis()->SetTitle("Single Electron Energy [MeV]");
        h_el_energy_diff->GetYaxis()->SetTitle("Events");
        h_el_energy_diff->Draw("E");
        c_el_energy_diff->SaveAs("c_el_energy_diff.C");
        c_el_energy_diff->SaveAs("c_el_energy_diff.png");
        c_el_energy_diff->SaveAs("c_el_energy_diff.pdf");
        
        c_el_energy_pull = new TCanvas("c_el_energy_pull", "", 800, 600);
        //c_el_energy_pull->SetLogz();
        h_el_energy_pull->GetXaxis()->SetTitle("Single Electron Energy [MeV]");
        h_el_energy_pull->GetYaxis()->SetTitle("Pull #sigma");
        h_el_energy_pull->Draw("hist");
        c_el_energy_pull->SaveAs("c_el_energy_pull.C");
        c_el_energy_pull->SaveAs("c_el_energy_pull.png");
        c_el_energy_pull->SaveAs("c_el_energy_pull.pdf");

    }
    */

    for(std::vector<SubAnalysis*>::iterator it{_subanalysis_.begin()}; it != _subanalysis_.end(); ++ it)
    {
        (*it)->SensitivityMeasurementChisquare1();
    }

}
