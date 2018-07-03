// self header
#include "Analysis.hpp"


// local headers
#include "aux.hpp"


// C++ headers
#include <iostream>
#include <fstream>



////////////////////////////////////////////////////////////////////////////////
// SCALE SUMMED HISTOGRAM TO MATCH AMPLITUDE
////////////////////////////////////////////////////////////////////////////////

void Analysis::SummedEnergyFit()
{

    // reset
    non_empty_bins_fit = 0;

    ////////////////////////////////////////////////////////////////////////
    // 1d histogram fit and sensitivity
    ////////////////////////////////////////////////////////////////////////

    f_el_energy_sum_original = new TF1("f_el_energy_sum_original", fit_function, 0.0, 4.0, 1 + 2 * (h_el_energy_sum_reweight->GetNbinsX() + 1));
    f_el_energy_sum_original->SetParameter(0, 1.0); // set initial amplitude parameter to 1.0
    f_el_energy_sum_original->SetNpx(100000);
    // set the "parameters" (constants)
    // these are the values from the histogram
    {
        Int_t ix{0};
        Int_t jx{0};
        Double_t par_ix{0.0};
        Double_t par_jx{0.0};

        Int_t nbx{h_el_energy_sum_reweight->GetNbinsX()};

        for(Int_t i{0}; i <= nbx; ++ i)
        {
            ix = 2 * i + 1;
            jx = 2 * i + 2;

            par_ix = h_el_energy_sum_reweight->GetXaxis()->GetBinLowEdge(i + 1);
            par_jx = h_el_energy_sum_reweight->GetBinContent(i + 1);

            //std::cout << "x=" << par_ix << " y=" << par_jx << std::endl;

            f_el_energy_sum_original->FixParameter(ix, par_ix);
            f_el_energy_sum_original->FixParameter(jx, par_jx);
        }

        ix = 2 * nbx + 3;
        ix = 2 * nbx + 4;
        
        par_ix = h_el_energy_sum_reweight->GetXaxis()->GetBinLowEdge(nbx + 1) + h_el_energy_sum_reweight->GetBinWidth(1);
        par_jx = 0.0; // content doesn't make sense here, don't use "overflow" by accident

        // final 2, marks the end point of final "bin"
        f_el_energy_sum_original->FixParameter(ix, par_ix);
        f_el_energy_sum_original->FixParameter(jx, par_jx);
        
        // TODO: am i putting the right information into the fitting function and am
        // i fitting the correct function?
        if(fit_subrange == false)
        {
            h_el_energy_sum_original->Fit("f_el_energy_sum_original", "0");
        }
        else if(fit_subrange == true)
        {
            h_el_energy_sum_original->Fit("f_el_energy_sum_original", "0", "", 2.0, 4.0);
        }

    }
    std::cout << "Amplitude parameter: " << f_el_energy_sum_original->GetParameter(0) << std::endl;
    std::cout << "                err: " << f_el_energy_sum_original->GetParError(0) << std::endl;
    //std::cout << "         chi square: " << f_el_energy_sum_original->GetChisquare() / h_el_energy_sum_original->GetNbinsX() << std::endl;
    std::cout << "         chi square: " << f_el_energy_sum_original->GetChisquare() << std::endl;
    for(Int_t i{1}; i <= h_el_energy_sum_original->GetNbinsX(); ++ i)
    {
        if(h_el_energy_sum_original->GetBinContent(i) != 0.0) ++ non_empty_bins_fit;
    }
    std::cout << " degrees of freedom: " << non_empty_bins_fit << std::endl;
    std::cout << " chi square reduced: " << f_el_energy_sum_original->GetChisquare() / (Double_t)non_empty_bins_fit << std::endl;

    // scale the red histogram using the amplitude parameter
    h_el_energy_sum_reweight->Scale(f_el_energy_sum_original->GetParameter(0));
    h_el_energy_reweight->Scale(f_el_energy_sum_original->GetParameter(0));
    h_el_energy_2d_reweight->Scale(f_el_energy_sum_original->GetParameter(0));


    ////////////////////////////////////////////////////////////////////////////
    // FIT CANVAS OUTPUT SUMMED ENERGY
    ////////////////////////////////////////////////////////////////////////////

    // print summed distribution
    Double_t max{0.0};
    Double_t min{0.0};
    const Double_t max_log_mode_2{100.0e3};
    const Double_t max_nolog_mode_2{80.0e3};
    const std::string dir_log_mode_2("el_energy_sum_log_dir");
    const std::string dir_nolog_mode_2("el_energy_sum_nolog_dir");
    std::string dir; //{std::to_string(epsilon_31)};
    if(log_mode) { min = 0.1; max = max_log_mode_2; dir = dir_log_mode_2; }
    else { min = 0.0; max = max_nolog_mode_2; dir = dir_nolog_mode_2; }
    CanvasFactorySettings settings_2("Sum Electron Energy [MeV]", "Events", min, max, log_mode);
    settings_2.SetDrawOption("E");
    CanvasFactory factory_2(settings_2);
    factory_2.Canvas("el_energy_sum", dir, h_el_energy_sum_original, "Baseline", h_el_energy_sum_reweight, "Reweighted");
    // TODO: settings should be passed to canvas in case settings should be changed?
    // or setsettings function should be provided

}


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


////////////////////////////////////////////////////////////////////////////////
// 2d histogram chisquare sensitivity
////////////////////////////////////////////////////////////////////////////////

void Analysis::SensitivityMeasurementChisquare2()
{

    // reset
    non_empty_bins_2d = 0;

    // get chi-square for single electron histograms
    if(fit_subrange == false)
    {

        for(Int_t j{1}; j <= h_el_energy_2d_original->GetNbinsY(); ++ j)
        {
            for(Int_t i{1}; i <= h_el_energy_2d_original->GetNbinsX(); ++ i)
            {
                if(h_el_energy_2d_original->GetBinContent(i, j) != 0.0) ++ non_empty_bins_2d;
            }
        }

        // get chi-square for single electron histograms
        // Note: no subrange for 2d histogram
        //if(fit_subrange == false)
        //{
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

        // get chi-square for single electron histograms
        // Note: no subrange for 2d histogram
        //if(fit_subrange == false)
        //{
        sensitivity_chisquare_2d = chi_square_test(h_el_energy_2d_reweight, h_el_energy_2d_original, 2.0, 4.0, non_empty_bins_2d);
        std::cout << "chi square of 2 electron, 2.0 MeV - 2.0 MeV: " << sensitivity_chisquare_2d << std::endl;
        std::cout << " degrees of freedom (2d): " << non_empty_bins_2d << std::endl;
        std::cout << " chi square reduced: " << sensitivity_chisquare_2d / (Double_t)non_empty_bins_2d << std::endl;
        
    }

    // NOTE TO SELF: There is no 2d fit, the 2d histogram is used to evaluate
    // the sensitivity, therefore there is only a chi-square test


    if(_batch_mode_ == false)
    {

        ////////////////////////////////////////////////////////////////////////
        // 2D ORIGINAL AND REWEIGHT CANVAS OUTPUT
        ////////////////////////////////////////////////////////////////////////

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
                }
                else
                {
                    h_el_energy_2d_pull->SetBinContent(i, j, content / error);
                }
            }
        }
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
 
}


////////////////////////////////////////////////////////////////////////////////
// 1D SINGLE ELECTRON ENERGY LOG LIKELIHOOD SENSITIVITY
////////////////////////////////////////////////////////////////////////////////

void Analysis::SensitivityMeasurementLoglikelihood1()
{

    // TODO: move histograms which are intrinsic to this function (measurement)
    // to this function - avoids "replacing existing th1" error
    // do not need these histos elsewhere

    const std::string eps_string{std::to_string(epsilon_31)};

    /*
    h_el_energy_data->Reset("ICESM");
    h_el_energy_prob->Reset("ICESM");
    h_el_energy_diff_data_orig->Reset("ICESM");
    h_el_energy_diff_data_rw->Reset("ICESM");
    h_el_energy_ll->Reset("ICESM");
    */

    delete h_ll;
    h_ll = nullptr;


    for(Int_t count{0}; count < number_of_pseudo_experiments; ++ count)
    {

        delete h_el_energy_data;
        delete h_el_energy_prob;
        delete h_el_energy_diff_data_rw;
        delete h_el_energy_diff_data_orig;
        h_el_energy_data = nullptr;
        h_el_energy_prob = nullptr;
        h_el_energy_diff_data_rw = nullptr;
        h_el_energy_diff_data_orig = nullptr;


        ////////////////////////////////////////////////////////////////////////
        // SINGLE ELECTRON ENERGY PSEUDODATA METHOD
        ////////////////////////////////////////////////////////////////////////

        // TODO: name of histograms which we iterate over
        // log likelihood method
        // create "data" histogram - poisson generated data
        std::string h_name{std::string("h_el_energy_data_") + eps_string + std::string("_") + std::to_string(count)};
        h_el_energy_data = new TH1I(h_name.c_str(), "", num_bins, 0.0, 4.0);
        for(Int_t ix{1}; ix <= h_el_energy_sum_original->GetNbinsX(); ++ ix)
        {
            // this is the input lambda value
            //Double_t content{h_el_energy_sum_original->GetBinContent(ix)};
            //Double_t lambda{h_el_energy_sum_original->GetBinContent(ix)};
            // NOTE: changed 2018-06-12: lambda was obtained from wrong distribution
            Double_t lambda{h_el_energy_original->GetBinContent(ix)};
            
            // Note: before implementing full pseudo experiments,
            // instead of drawing number randomly from a poisson with
            // mean lambda (which returns an integer),
            // instead did:
            // round the lambda value and pretend it is data
            //content = std::round(content);
            //content = std::lround(content);
            //h_el_energy_data->SetBinContent(ix, content);

            // Note: moved to pseudorandom data rather than rounding
            Int_t poisson_result{gen.Poisson(lambda)};
            h_el_energy_data->SetBinContent(ix, poisson_result);
        }

        // compute poisson likelihood for each bin
        Double_t likelihood{1.0};
        Double_t likelihood_sum{0.0};
        std::vector<Double_t> likelihood_sum_vec;
        std::string h_name_prob{std::string("h_el_energy_prob_") + eps_string + std::string("_") + std::to_string(count)}; 
        h_el_energy_prob = new TH1D(h_name_prob.c_str(), "", num_bins, 0.0, 4.0);
        for(Int_t ix{1}; ix <= h_el_energy_sum_original->GetNbinsX(); ++ ix)
        {
            //Double_t lambda{h_el_energy_sum_reweight->GetBinContent(ix)}; // reweighted
            Double_t lambda{h_el_energy_reweight->GetBinContent(ix)}; // reweighted
            // NOTE: changed 2018-06-12: lambda was obtained from wrong distribution
            // NOTE: above lambda is the output lambda, it is for the other distribution
            // rather than the input lambda in the previous code block

            //Int_t data{std::lround(h_el_energy_data->GetBinContent(ix))}; // rounded data (not -> baseline <- ?)
            Int_t data{(Int_t)std::lround(h_el_energy_data->GetBinContent(ix))}; // rounded data, is already rounded, is integer
            // Note: lround required as GetBinContent returns Double_t
            // TODO: do not need round?
            Double_t poi{TMath::Poisson(data, lambda)};
            //std::cout << "bin index = " << ix << ", poisson = " << poi << ", lambda = " << lambda << ", data = " << data << std::endl;
            likelihood *= poi;

            Double_t poi_log{std::log(poi)};
            likelihood_sum += poi_log;
            likelihood_sum_vec.push_back(poi_log);

            h_el_energy_prob->SetBinContent(ix, poi);

            if(h_el_energy_sum_original->GetBinCenter(ix) > 2.0)
            {
                //std::cout << "Bin center: " << h_el_energy_sum_original->GetBinCenter(ix) << " lambda=" << lambda << " data=" << data << " prob=" << poi << std::endl; 
            }
                
            if(count == 0)
            {
                //std::cout << "lambda=," << lambda << ",data=," << data << ",prob=," << poi << std::endl;
            }
        }
        //std::cout << "likelihood=" << likelihood << std::endl;
        //std::cout << "likelihood = " << likelihood << std::endl;
        Double_t log_likelihood{std::log(likelihood)};
        vec_ll.push_back(-2.0 * log_likelihood);

        std::cout << "LL=" << -2.0 * log_likelihood << ", LL2=" << -2.0 * likelihood_sum << ", diff=" << (-2.0 * log_likelihood) - (-2.0 * likelihood_sum) << std::endl;
        std::sort(likelihood_sum_vec.begin(), likelihood_sum_vec.end());
        Double_t likelihood_sum_vec_sum{0.0};
        for(std::vector<Double_t>::const_iterator it{likelihood_sum_vec.cbegin()}; it != likelihood_sum_vec.cend(); ++ it)
        {
            likelihood_sum_vec_sum += *it;
        }
        std::cout << "LL3=" << -2.0 * likelihood_sum_vec_sum << std::endl;

        // create difference histogram, data - reweight
        std::string h_name_diff_data_rw{std::string("h_el_energy_diff_data_rw_") + eps_string + std::string("_") + std::to_string(count)};
        h_el_energy_diff_data_rw = new TH1D(h_name_diff_data_rw.c_str(), "", num_bins, 0.0, 4.0);
        h_el_energy_diff_data_rw->SetStats(0);
        for(Int_t i{1}; i <= h_el_energy_diff_data_rw->GetNbinsX(); ++ i)
        {
            Double_t content1{h_el_energy_data->GetBinContent(i)};
            Double_t content2{h_el_energy_reweight->GetBinContent(i)};
            Double_t error1{h_el_energy_data->GetBinError(i)};
            // note: do not use error or reweighted
            Double_t error2{0.0 * h_el_energy_reweight->GetBinError(i)};
            Double_t content{content1 - content2};
            Double_t error{std::sqrt(error1 * error1 + error2 * error2)};
            h_el_energy_diff_data_rw->SetBinContent(i, content);
            h_el_energy_diff_data_rw->SetBinError(i, error);
        }

        // create difference histogram, data - original
        std::string h_name_diff_data_orig{std::string("h_el_energy_diff_data_orig_") + eps_string + std::string("_") + std::to_string(count)};
        h_el_energy_diff_data_orig = new TH1D(h_name_diff_data_orig.c_str(), "", num_bins, 0.0, 4.0);
        h_el_energy_diff_data_orig->SetStats(0);
        for(Int_t i{1}; i <= h_el_energy_diff_data_orig->GetNbinsX(); ++ i)
        {
            Double_t content1{h_el_energy_data->GetBinContent(i)};
            Double_t content2{h_el_energy_original->GetBinContent(i)};
            std::cout << content1 << "   " << content2 << std::endl;
            Double_t error1{h_el_energy_data->GetBinError(i)};
            // note: do not use error or reweighted
            Double_t error2{0.0 * h_el_energy_original->GetBinError(i)};
            Double_t content{content1 - content2};
            Double_t error{std::sqrt(error1 * error1 + error2 * error2)};
            h_el_energy_diff_data_orig->SetBinContent(i, content);
            h_el_energy_diff_data_orig->SetBinError(i, error);
        }


        ////////////////////////////////////////////////////////////////////
        // CANVAS OUTPUT
        ////////////////////////////////////////////////////////////////////

        if(false)
        {

            //const Double_t canvas_max{2.5e5};
            const Double_t canvas_max{1.0e6};
            //const Double_t canvas_min{0.0};
            const Double_t canvas_min{1.0e-1};
            const std::string canvas_dir("./pseudoexperiments1d");
            CanvasFactorySettings settings("Single Electron Energy [MeV]", "Events", canvas_min, canvas_max, false);
            settings.SetLogMode(true);
            settings.SetDrawOption("E");
            settings.SetOutputPNGOnly();
            settings.SetOutputROOT(false);
            settings.SetOutputPNG(false);
            CanvasFactory factory(settings);
            std::string c_name_data{std::string("el_energy_data/el_energy_data_") + eps_string + std::string("_") + std::to_string(count)};
            factory.Canvas(c_name_data, canvas_dir, h_el_energy_original, "Baseline", h_el_energy_reweight, "Reweighted", h_el_energy_data, "Pseudodata");

            settings.SetMax(1.1);
            settings.SetMin(1.0e-10); // 1.0e-35
            settings.SetDrawOption("hist");
            settings.SetLogMode(true);
            factory.Settings(settings);
            std::string c_name_prob{std::string("el_energy_prob/el_energy_prob") + eps_string + std::string("_") + std::to_string(count)};
            factory.Canvas(c_name_prob.c_str(), canvas_dir, h_el_energy_prob, "Probability");

            std::string c_name_diff_data_rw{canvas_dir + std::string("/el_energy_diff_data_rw/c_el_energy_diff_data_rw_") + std::to_string(count)};
            c_el_energy_diff_data_rw = new TCanvas(c_name_diff_data_rw.c_str(), "", 800, 600);
            //c_el_energy_diff_data_rw->SetRightMargin(0.12);
            //c_el_energy_diff_data_rw->SetLogz();
            h_el_energy_diff_data_rw->GetXaxis()->SetTitle("Single Electron Energy [MeV]");
            h_el_energy_diff_data_rw->GetYaxis()->SetTitle("Events (Pseudodata - Reweight)");
            h_el_energy_diff_data_rw->Draw("E");
            //c_el_energy_diff_data_rw->SaveAs((c_name_diff_data_rw + std::string(".C")).c_str());
            c_el_energy_diff_data_rw->SaveAs((c_name_diff_data_rw + std::string(".png")).c_str());
            //c_el_energy_diff_data_rw->SaveAs((c_name_diff_data_rw + std::string(".pdf")).c_str());

            std::string c_name_diff_data_orig{canvas_dir + std::string("/el_energy_diff_data_orig/c_el_energy_diff_data_orig_") + std::to_string(count)};
            c_el_energy_diff_data_orig = new TCanvas(c_name_diff_data_orig.c_str(), "", 800, 600);
            //c_el_energy_diff_data_orig->SetRightMargin(0.12);
            //c_el_energy_diff_data_orig->SetLogz();
            h_el_energy_diff_data_orig->GetXaxis()->SetTitle("Single Electron Energy [MeV]");
            h_el_energy_diff_data_orig->GetYaxis()->SetTitle("Events (Pseudodata - Baseline)");
            h_el_energy_diff_data_orig->Draw("E");
            //c_el_energy_diff_data_orig->SaveAs((c_name_diff_data_orig + std::string(".C")).c_str());
            c_el_energy_diff_data_orig->SaveAs((c_name_diff_data_orig + std::string(".png")).c_str());
            //c_el_energy_diff_data_orig->SaveAs((c_name_diff_data_orig + std::string(".pdf")).c_str());

            delete c_el_energy_diff_data_rw;
            delete c_el_energy_diff_data_orig;

        }




    }


    ////////////////////////////////////////////////////////////////////////////
    // GET DISTRIBUTION OF LL MEASUREMENTS
    ////////////////////////////////////////////////////////////////////////////

    //Double_t min{std::min_element(vec_ll.begin(), vec_ll.end())};
    //Double_t max{std::max_element(vec_ll.begin(), vec_ll.end())};
    std::pair<std::vector<Double_t>::iterator, std::vector<Double_t>::iterator> min_max_pair{std::minmax_element(vec_ll.begin(), vec_ll.end())};
    Double_t min{*min_max_pair.first};
    Double_t max{*min_max_pair.second};
    h_ll = new TH1D("h_ll", "", 100, min, max);
    //h_ll->GetXaxis()->SetTitle("Log Likelihood Value");
    h_ll->GetXaxis()->SetTitle("Chi Square Value");
    h_ll->GetYaxis()->SetTitle("Number of Pseudo Experiments");
    // fill the histogram
    for(std::vector<Double_t>::const_iterator it{vec_ll.cbegin()}; it != vec_ll.cend(); ++ it)
    {
        h_ll->Fill(*it);
    }

    c_ll = new TCanvas("c_ll", "", 800, 600);
    h_ll->Draw("E");
    c_ll->SaveAs("c_ll.png");

    delete c_ll;
        
}


////////////////////////////////////////////////////////////////////////////////
// 2D SINGLE ELECTRON ENERGY LOG LIKELIHOOD SENSITIVITY
////////////////////////////////////////////////////////////////////////////////

void Analysis::SensitivityMeasurementLoglikelihood2()
{

    
    delete h_ll_2d;
    h_ll_2d = nullptr;
    
    
    for(Int_t count{0}; count < number_of_pseudo_experiments_2d; ++ count)
    {


        delete h_el_energy_2d_data;
        delete h_el_energy_2d_prob;
        delete h_el_energy_2d_diff_data_rw;
        delete h_el_energy_2d_diff_data_orig;
        h_el_energy_2d_data = nullptr;
        h_el_energy_2d_prob = nullptr;
        h_el_energy_2d_diff_data_rw = nullptr;
        h_el_energy_2d_diff_data_orig = nullptr;

        ////////////////////////////////////////////////////////////////////////
        // INDEPENDENT SINGLE ELECTRON ENERGY
        ////////////////////////////////////////////////////////////////////////

        // log likelihood method
        // create 2d "data" histogram - poisson generated data

        // chisquare method
        // create 2d "data" histogram
        // TODO which one

        std::string h_name_2d{std::string("h_el_energy_2d_data_") + std::to_string(count)};
        h_el_energy_2d_data = new TH2I(h_name_2d.c_str(), "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
        h_el_energy_2d_data->SetStats(0);
        for(Int_t jx{1}; jx <= h_el_energy_2d_original->GetNbinsY(); ++ jx)
        {
            for(Int_t ix{1}; ix <= h_el_energy_2d_original->GetNbinsX(); ++ ix)
            {
                // this is the input lambda value
                Double_t lambda{h_el_energy_2d_original->GetBinContent(ix, jx)};
                
                // pseudorandom data
                Int_t poisson_result{gen.Poisson(lambda)};
                h_el_energy_2d_data->SetBinContent(ix, jx, poisson_result);
            }
        }

        // compute poisson likelihood for each bin
        Double_t likelihood_2d{1.0};
        std::string h_name_prob_2d{std::string("h_el_energy_2d_prob_") + std::to_string(count)}; 
        h_el_energy_2d_prob = new TH2D(h_name_prob_2d.c_str(), "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
        h_el_energy_2d_prob->SetStats(0);
        for(Int_t jx{1}; jx <= h_el_energy_2d_original->GetNbinsY(); ++ jx)
        {
            for(Int_t ix{1}; ix <= h_el_energy_2d_original->GetNbinsX(); ++ ix)
            {
                Double_t lambda{h_el_energy_2d_reweight->GetBinContent(ix, jx)}; // reweighted
                // NOTE: above lambda is the output lambda, it is for the other distribution
                // rather than the input lambda in the previous code block

                //Int_t data{std::lround(h_el_energy_data->GetBinContent(ix))}; // rounded data (not -> baseline <- ?)
                Int_t data{(Int_t)std::lround(h_el_energy_2d_data->GetBinContent(ix, jx))};
                // Note: lround required because function returns Double_t
                Double_t poi{TMath::Poisson(data, lambda)};
                likelihood_2d *= poi;

                h_el_energy_2d_prob->SetBinContent(ix, jx, poi);
            }
        }
        std::cout << "likelihood (2d) = " << likelihood_2d << std::endl;
        Double_t log_likelihood_2d{std::log(likelihood_2d)};
        vec_ll_2d.push_back(-2.0 * log_likelihood_2d);
    
        
        // create difference histogram - difference between data and reweighted
        std::string h_name_2d_diff_data_rw{std::string("h_el_energy_2d_diff_data_rw_") + std::to_string(count)};
        h_el_energy_2d_diff_data_rw = new TH2D(h_name_2d_diff_data_rw.c_str(), "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
        h_el_energy_2d_diff_data_rw->SetStats(0);
        for(Int_t j{1}; j <= h_el_energy_2d_diff_data_rw->GetNbinsY(); ++ j)
        {
            for(Int_t i{1}; i <= h_el_energy_2d_diff_data_rw->GetNbinsX(); ++ i)
            {
                Double_t content1{h_el_energy_2d_data->GetBinContent(i, j)};
                Double_t content2{h_el_energy_2d_reweight->GetBinContent(i, j)};
                Double_t error1{h_el_energy_2d_data->GetBinError(i, j)};
                // note: do not use error or reweighted
                Double_t error2{0.0 * h_el_energy_2d_reweight->GetBinError(i, j)};
                Double_t content{content1 - content2};
                Double_t error{std::sqrt(error1 * error1 + error2 * error2)};
                h_el_energy_2d_diff_data_rw->SetBinContent(i, j, content);
                h_el_energy_2d_diff_data_rw->SetBinError(i, j, error);
            }
        }

        // create difference histogram - difference between data and original
        std::string h_name_2d_diff_data_orig{std::string("h_el_energy_2d_diff_data_orig_") + std::to_string(count)};
        h_el_energy_2d_diff_data_orig = new TH2D(h_name_2d_diff_data_orig.c_str(), "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
        h_el_energy_2d_diff_data_orig->SetStats(0);
        for(Int_t j{1}; j <= h_el_energy_2d_diff_data_orig->GetNbinsY(); ++ j)
        {
            for(Int_t i{1}; i <= h_el_energy_2d_diff_data_orig->GetNbinsX(); ++ i)
            {
                Double_t content1{h_el_energy_2d_data->GetBinContent(i, j)};
                Double_t content2{h_el_energy_2d_original->GetBinContent(i, j)};
                Double_t error1{h_el_energy_2d_data->GetBinError(i, j)};
                // note: do not use error or reweighted
                Double_t error2{0.0 * h_el_energy_2d_original->GetBinError(i, j)};
                Double_t content{content1 - content2};
                Double_t error{std::sqrt(error1 * error1 + error2 * error2)};
                h_el_energy_2d_diff_data_orig->SetBinContent(i, j, content);
                h_el_energy_2d_diff_data_orig->SetBinError(i, j, error);
            }
        }

        ////////////////////////////////////////////////////////////////////
        // CANVAS OUTPUT
        ////////////////////////////////////////////////////////////////////

        if(false)
        {

            // TODO: multiple experiemnts problem
            // if statement for batch mode
            c_el_energy_2d_data = new TCanvas("c_el_energy_2d_data", "", 800, 600);
            c_el_energy_2d_data->SetRightMargin(0.17);
            //c_el_energy_2d_data->SetLogz();
            h_el_energy_2d_data->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
            h_el_energy_2d_data->GetYaxis()->SetTitle("High Energy Electron [MeV]");
            h_el_energy_2d_data->GetZaxis()->SetTitle("Events");
            h_el_energy_2d_data->GetZaxis()->SetTitleOffset(1.6);
            h_el_energy_2d_data->Draw("colz");
            c_el_energy_2d_data->SaveAs("c_el_energy_2d_data.C");
            c_el_energy_2d_data->SaveAs("c_el_energy_2d_data.png");
            c_el_energy_2d_data->SaveAs("c_el_energy_2d_data.pdf");

            c_el_energy_2d_prob = new TCanvas("c_el_energy_2d_prob", "", 800, 600);
            c_el_energy_2d_prob->SetRightMargin(0.15);
            c_el_energy_2d_prob->SetLogz();
            h_el_energy_2d_prob->SetMinimum(1.0e-6);
            h_el_energy_2d_prob->SetMaximum(1.0e0);
            h_el_energy_2d_prob->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
            h_el_energy_2d_prob->GetYaxis()->SetTitle("High Energy Electron [MeV]");
            h_el_energy_2d_prob->GetZaxis()->SetTitle("Events");
            h_el_energy_2d_prob->GetZaxis()->SetTitleOffset(1.2);
            h_el_energy_2d_prob->Draw("colz");
            c_el_energy_2d_prob->SaveAs("c_el_energy_2d_prob.C");
            c_el_energy_2d_prob->SaveAs("c_el_energy_2d_prob.png");
            c_el_energy_2d_prob->SaveAs("c_el_energy_2d_prob.pdf");

            c_el_energy_2d_diff_data_rw = new TCanvas("c_el_energy_2d_diff_data_rw", "", 800, 600);
            c_el_energy_2d_diff_data_rw->SetRightMargin(0.15);
            //c_el_energy_2d_diff_data_rw->SetLogz();
            h_el_energy_2d_diff_data_rw->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
            h_el_energy_2d_diff_data_rw->GetYaxis()->SetTitle("High Energy Electron [MeV]");
            h_el_energy_2d_diff_data_rw->GetZaxis()->SetTitle("Events");
            h_el_energy_2d_diff_data_rw->GetZaxis()->SetTitleOffset(1.2);
            h_el_energy_2d_diff_data_rw->Draw("colz");
            c_el_energy_2d_diff_data_rw->SaveAs("c_el_energy_2d_diff_data_rw.C");
            c_el_energy_2d_diff_data_rw->SaveAs("c_el_energy_2d_diff_data_rw.png");
            c_el_energy_2d_diff_data_rw->SaveAs("c_el_energy_2d_diff_data_rw.pdf");

            c_el_energy_2d_diff_data_orig = new TCanvas("c_el_energy_2d_diff_data_orig", "", 800, 600);
            c_el_energy_2d_diff_data_orig->SetRightMargin(0.15);
            //c_el_energy_2d_diff_data_orig->SetLogz();
            h_el_energy_2d_diff_data_orig->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
            h_el_energy_2d_diff_data_orig->GetYaxis()->SetTitle("High Energy Electron [MeV]");
            h_el_energy_2d_diff_data_orig->GetZaxis()->SetTitle("Events");
            h_el_energy_2d_diff_data_orig->GetZaxis()->SetTitleOffset(1.2);
            h_el_energy_2d_diff_data_orig->Draw("colz");
            c_el_energy_2d_diff_data_orig->SaveAs("c_el_energy_2d_diff_data_orig.C");
            c_el_energy_2d_diff_data_orig->SaveAs("c_el_energy_2d_diff_data_orig.png");
            c_el_energy_2d_diff_data_orig->SaveAs("c_el_energy_2d_diff_data_orig.pdf");


            delete c_el_energy_2d_data;
            delete c_el_energy_2d_prob;
            delete c_el_energy_2d_diff_data_rw;
            delete c_el_energy_2d_diff_data_orig;

        }

        ////////////////////////////////////////////////////////////////////////
        // SAVE TO ROOT FILE
        ////////////////////////////////////////////////////////////////////////

        /*
        TFile *f_el_energy_2d_loglikelihood = new TFile("f_el_energy_2d_loglikelihood.root", "RECREATE");
        h_el_energy_2d_original->Write();
        //c_el_energy_2d_original->Write(); // Note: not defined in this scope
        h_el_energy_2d_reweight->Write();
        //c_el_energy_2d_reweight->Write(); // Note: not defined in this scope
        h_el_energy_2d_diff_data_rw->Write();
        c_el_energy_2d_diff_data_rw->Write();
        h_el_energy_2d_data->Write();
        c_el_energy_2d_data->Write();
        h_el_energy_2d_prob->Write();
        c_el_energy_2d_prob->Write();
        f_el_energy_2d_loglikelihood->Close();
        //delete f_el_energy_2d;
        */

    }


    ////////////////////////////////////////////////////////////////////////////
    // GET DISTRIBUTION OF LL MEASUREMENTS
    ////////////////////////////////////////////////////////////////////////////

    std::pair<std::vector<Double_t>::iterator, std::vector<Double_t>::iterator> min_max_pair_2d{std::minmax_element(vec_ll_2d.begin(), vec_ll_2d.end())};
    Double_t min_2d{*min_max_pair_2d.first};
    Double_t max_2d{*min_max_pair_2d.second};
    h_ll_2d = new TH1D("h_ll_2d", "", 100, min_2d, max_2d);
    //h_ll_2d->GetXaxis()->SetTitle("Log Likelihood Value");
    h_ll_2d->GetXaxis()->SetTitle("Chi Square Value");
    h_ll_2d->GetYaxis()->SetTitle("Number of Pseudo Experiments");
    // fill the histogram
    for(std::vector<Double_t>::const_iterator it{vec_ll_2d.cbegin()}; it != vec_ll_2d.cend(); ++ it)
    {
        h_ll_2d->Fill(*it);
    }

    
    c_ll_2d = new TCanvas("c_ll_2d", "", 800, 600);
    h_ll_2d->Draw("E");
    c_ll_2d->SaveAs("c_ll_2d.png");
    delete c_ll_2d;

}


////////////////////////////////////////////////////////////////////////////////
// PRINT OUTPUT TO FILE
////////////////////////////////////////////////////////////////////////////////

void Analysis::PrintOutputToFile()
{


    // print entries
    //std::cout << "Number of entries in each histogram: h_el_energy_original: " << h_el_energy_original->GetEntries() << " h_el_energy_reweight: " << h_el_energy_reweight->GetEntries() << std::endl;


    {
        // add data to data output file
        std::ofstream of_data(output_filename.c_str(), std::ios::app);
        if(of_data.tellp() == 0)
        {
            of_data << "epsilon_31,chisquare (fit),degrees of freedom,chisquare (fit reduced),chisquare (sensitivity),chisquare (sensitivity reduced),-2log(l)" << std::endl;
        }
        of_data << epsilon_31 << ','
                << f_el_energy_sum_original->GetChisquare() << ','
                << non_empty_bins_fit << ','
                << f_el_energy_sum_original->GetChisquare() / (Double_t)non_empty_bins_fit << ','
                << sensitivity_chisquare << ','
                << sensitivity_chisquare / (Double_t)non_empty_bins; //<< ','
                //<< -2.0 * log_likelihood
                // TODO what goes here
        for(std::vector<Double_t>::const_iterator it{vec_ll.cbegin()}; it != vec_ll.cend(); ++ it)
        {
            of_data << ',' << *it;
        }
        of_data << std::endl;


        // new output format, do not change above - will break chi-square code
        std::string filename_chisq(output_filename.substr(0, output_filename.find(".txt")) + std::string("_chisq") + std::string(".txt"));
        std::cout << filename_chisq << std::endl;
        std::ofstream of_data_chisq(filename_chisq.c_str(), std::ios::app);
        if(of_data_chisq.tellp() == 0)
        {
            of_data_chisq << "epsilon_31,fit (summed): chisq,reduced chisq,dof,sensitivity (1d single): chisq,reduced chisq,dof,sensitivity (2d single): chisq,reduced chisq,dof" << std::endl;
        }
        of_data_chisq << epsilon_31 << ','
                      << f_el_energy_sum_original->GetChisquare() << ','
                      << f_el_energy_sum_original->GetChisquare() / (Double_t)non_empty_bins_fit << ','
                      << non_empty_bins_fit << ','
                      << sensitivity_chisquare << ','
                      << sensitivity_chisquare / (Double_t)non_empty_bins << ','
                      << non_empty_bins << ','
                      << sensitivity_chisquare_2d << ','
                      << sensitivity_chisquare_2d / (Double_t)non_empty_bins_2d << ','
                      << non_empty_bins_2d << std::endl;

    }

    {
        // 1d ll data
        std::string filename_ll(output_filename.substr(0, output_filename.find(".txt")) + std::string("_ll") + std::string(".txt"));
        std::cout << filename_ll << std::endl;
        std::ofstream of_data_ll(filename_ll.c_str(), std::ios::app);
        for(std::vector<Double_t>::const_iterator it{vec_ll.cbegin()}; it != vec_ll.cend(); ++ it)
        {
            of_data_ll << *it;
            if(it + 1 != vec_ll.cend()) of_data_ll << ',';
        }
        of_data_ll << std::endl;

    }

    {
        // 2d ll data
        std::string filename_ll_2d(output_filename.substr(0, output_filename.find(".txt")) + std::string("_ll_2d") + std::string(".txt"));
        std::cout << filename_ll_2d << std::endl;
        std::ofstream of_data_ll_2d(filename_ll_2d.c_str(), std::ios::app);
        for(std::vector<Double_t>::const_iterator it{vec_ll_2d.cbegin()}; it != vec_ll_2d.cend(); ++ it)
        {
            of_data_ll_2d << *it;
            if(it + 1 != vec_ll_2d.cend()) of_data_ll_2d << ',';
        }
        of_data_ll_2d << std::endl;
    }


}
