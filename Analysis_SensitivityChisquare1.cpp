////////////////////////////////////////////////////////////////////////////////
// 1D SINGLE ELECTRON ENERGY CHISQUARE SENSITIVITY
////////////////////////////////////////////////////////////////////////////////

void Analysis::SensitivityMeasurementChisquare1()
{

    // reset
    non_empty_bins = 0;

    ////////////////////////////////////////////////////////////////////////////
    // SINGLE ELECTRON ENERGY CHISQUARE METHOD
    ////////////////////////////////////////////////////////////////////////////

    // get chi-square for single electron histograms
    if(fit_subrange == false)
    {
        for(Int_t i{1}; i <= h_el_energy_original->GetNbinsX(); ++ i)
        {
            if(h_el_energy_original->GetBinContent(i) != 0.0) ++ non_empty_bins;
        }
        sensitivity_chisquare = chi_square_test(h_el_energy_reweight, h_el_energy_original);
        std::cout << "chi square of single electron: " << sensitivity_chisquare << std::endl;
        std::cout << " degrees of freedom: " << non_empty_bins << std::endl;
        std::cout << " chi square reduced: " << sensitivity_chisquare / (Double_t)non_empty_bins << std::endl;

    }
    else if(fit_subrange == true)
    {
        for(Int_t i{1}; i <= h_el_energy_original->GetNbinsX(); ++ i)
        {
            if(h_el_energy_original->GetBinCenter(i) >= 2.0)
            {
                if(h_el_energy_original->GetBinContent(i) != 0.0) ++ non_empty_bins;
            }
        }
        sensitivity_chisquare = chi_square_test(h_el_energy_reweight, h_el_energy_original, 2.0, 4.0);
        std::cout << "chi square of single electron, 2.0 MeV - 4.0 MeV: " << sensitivity_chisquare << std::endl;
        std::cout << " degrees of freedom: " << non_empty_bins << std::endl;
        std::cout << " chi square reduced: " << sensitivity_chisquare / (Double_t)non_empty_bins << std::endl;
    }


    ////////////////////////////////////////////////////////////////////////
    // CREATE DIFFERENCE AND PULL HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////

    // create difference histogram
    h_el_energy_diff = new TH1D("h_el_energy_diff", "", num_bins, 0.0, 4.0);
    h_el_energy_diff->SetStats(0);
    for(Int_t i{1}; i <= h_el_energy_diff->GetNbinsX(); ++ i)
    {
        Double_t content1{h_el_energy_original->GetBinContent(i)};
        Double_t content2{h_el_energy_reweight->GetBinContent(i)};
        Double_t error1{h_el_energy_original->GetBinError(i)};
        // note: do not use error or reweighted
        Double_t error2{0.0 * h_el_energy_reweight->GetBinError(i)};
        Double_t content{content1 - content2};
        Double_t error{std::sqrt(error1 * error1 + error2 * error2)};
        h_el_energy_diff->SetBinContent(i, content);
        h_el_energy_diff->SetBinError(i, error);
    }

    // create pull histogram
    h_el_energy_pull = new TH1D("h_el_energy_pull", "", num_bins, 0.0, 4.0);
    h_el_energy_pull->SetStats(0);
    for(Int_t i{1}; i <= h_el_energy_pull->GetNbinsX(); ++ i)
    {
        Double_t content{h_el_energy_diff->GetBinContent(i)};
        Double_t error{h_el_energy_diff->GetBinError(i)};
        if(error == 0.0)
        {
            //h_el_energy_pull->SetBinContent(i, j, 0.0);
            h_el_energy_pull->SetBinError(i, 0.0);
        }
        else
        {
            h_el_energy_pull->SetBinContent(i, content / error);
            h_el_energy_pull->SetBinError(i, 0.0);
        }
    }


    //HistogramWrapper2 hwrap_el_energy("el_energy", "x_axis_label_text", "y_axis_label_text", "E");
    //hwrap_el_energy.SetLogMode(log_mode);
    //hwrap_el_energy.Add("h_el_energy_original", h_el_energy_original);
    //hwrap_el_energy.Canvas();
    // TODO: make some nice histograms
    // X/Y AXIS: LABEL: TEXT, FONT SIZE, FONT
    // X/Y AXIS: TICK NUMBERS: FONT SIZE, FONT
    // X/Y AXIS LIMITS, RANGE, PARTICULARLY Y MAX
    // GRAPH COLORS


    ////////////////////////////////////////////////////////////////////////
    // SINGLE ELECTRON (1D) CHISQUARE CANVAS OUTPUT
    ////////////////////////////////////////////////////////////////////////

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


}