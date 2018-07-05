#include "Analysis.hpp"


////////////////////////////////////////////////////////////////////////////////
// 2D SINGLE ELECTRON ENERGY LOG LIKELIHOOD SENSITIVITY
////////////////////////////////////////////////////////////////////////////////

void Analysis::SensitivityMeasurementLoglikelihood2()
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
    
    delete h_ll_2d;
    h_ll_2d = nullptr;
    
    
    for(Int_t count{0}; count < number_of_pseudo_experiments_2d; ++ count)
    {
    
        std::cout << "count=" << count << std::endl;


        delete h_el_energy_2d_data;
        std::cout << "1 delete done" << std::endl;
        delete h_el_energy_2d_prob;
        std::cout << "2 delete done" << std::endl;
        delete h_el_energy_2d_diff_data_rw;
        std::cout << "3 delete done" << std::endl;
        delete h_el_energy_2d_diff_data_orig;
        std::cout << "4 delete done" << std::endl;
        h_el_energy_2d_data = nullptr;
        h_el_energy_2d_prob = nullptr;
        h_el_energy_2d_diff_data_rw = nullptr;
        h_el_energy_2d_diff_data_orig = nullptr;
        
        std::cout << "delete done" << std::endl;


        ////////////////////////////////////////////////////////////////////////
        // INDEPENDENT SINGLE ELECTRON ENERGY
        ////////////////////////////////////////////////////////////////////////

        // log likelihood method
        // create 2d "data" histogram - poisson generated data
        std::string h_name_2d{std::string("h_el_energy_2d_data_") + eps_string + std::string("_") + std::to_string(count)};
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
        std::cout << "data done" << std::endl;

        // compute poisson likelihood for each bin
        Double_t likelihood_2d{1.0};
        Double_t likelihood_sum_2d{0.0};
        std::string h_name_prob_2d{std::string("h_el_energy_2d_prob_") + eps_string + std::string("_") + std::to_string(count)}; 
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

                Double_t poi_log{std::log(poi)};
                likelihood_sum_2d += poi_log;

                h_el_energy_2d_prob->SetBinContent(ix, jx, poi);
            }
        }
        std::cout << "prob done" << std::endl;
        //std::cout << "likelihood (2d) = " << likelihood_2d << std::endl;
        Double_t log_likelihood_2d{std::log(likelihood_2d)};
        //vec_ll_2d.push_back(-2.0 * log_likelihood_2d);
        vec_ll_2d.push_back(-2.0 * likelihood_sum_2d);
    
        //std::cout << "LL=" << -2.0 * log_likelihood_2d << ", LL2=" << -2.0 * likelihood_sum_2d << ", diff=" << (-2.0 * log_likelihood_2d) - (-2.0 * likelihood_sum_2d) << std::endl;
        
        
        // create difference histogram - difference between data and reweighted
        std::string h_name_2d_diff_data_rw{std::string("h_el_energy_2d_diff_data_rw_") + eps_string + std::string("_") + std::to_string(count)};
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
        std::cout << "diff done" << std::endl;

        // create difference histogram - difference between data and original
        std::string h_name_2d_diff_data_orig{std::string("h_el_energy_2d_diff_data_orig_") + eps_string + std::string("_") + std::to_string(count)};
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
        std::cout << "other diff done" << std::endl;


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
    std::cout << "creating output" << std::endl;


    ////////////////////////////////////////////////////////////////////////////
    // GET DISTRIBUTION OF LL MEASUREMENTS
    ////////////////////////////////////////////////////////////////////////////

    std::pair<std::vector<Double_t>::iterator, std::vector<Double_t>::iterator> min_max_pair_2d{std::minmax_element(vec_ll_2d.begin(), vec_ll_2d.end())};
    Double_t min_2d{*min_max_pair_2d.first};
    Double_t max_2d{*min_max_pair_2d.second};
    std::string name_ll_2d{std::string("h_ll_2d_") + eps_string};
    Double_t center_2d{0.5 * (max_2d + min_2d)};
    Double_t width_2d{0.5 * (max_2d - min_2d)};
    width_2d *= 1.1;
    min_2d = center_2d - width_2d;
    max_2d = center_2d + width_2d;
    h_ll_2d = new TH1D(name_ll_2d.c_str(), "", 100, min_2d, max_2d);
    h_ll_2d->GetXaxis()->SetTitle("Chi Square Value");
    h_ll_2d->GetYaxis()->SetTitle("Number of Pseudo Experiments");
    // fill the histogram
    for(std::vector<Double_t>::const_iterator it{vec_ll_2d.cbegin()}; it != vec_ll_2d.cend(); ++ it)
    {
        h_ll_2d->Fill(*it);
    }
    
    std::string c_name_ll_2d{std::string("c_ll_2d_") + eps_string};
    c_ll_2d = new TCanvas(c_name_ll_2d.c_str(), "", 800, 600);
    h_ll_2d->Draw("E");
    c_ll_2d->SaveAs((c_name_ll_2d + std::string(".png")).c_str());
    delete c_ll_2d;

}