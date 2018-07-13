#include "Analysis.hpp"


////////////////////////////////////////////////////////////////////////////////
// 1D SINGLE ELECTRON ENERGY LOG LIKELIHOOD SENSITIVITY
////////////////////////////////////////////////////////////////////////////////

void SubAnalysis::SensitivityMeasurementLoglikelihood1()
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

    vec_ll.clear();


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

        // log likelihood method
        // create "data" histogram - poisson generated data
        std::string h_name{std::string("h_el_energy_data") + h_name_append + std::string("_") + eps_string + std::string("_") + std::to_string(count)};
        h_el_energy_data = new TH1I(h_name.c_str(), "", num_bins, 0.0, 4.0);
        h_el_energy_data->SetStats(0);
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
            Int_t poisson_result{gen->Poisson(lambda)};
            h_el_energy_data->SetBinContent(ix, poisson_result);
        }

        // compute poisson likelihood for each bin
        Double_t likelihood{1.0};
        Double_t likelihood_sum{0.0};
        std::vector<Double_t> likelihood_sum_vec;
        std::string h_name_prob{std::string("h_el_energy_prob") + h_name_append + std::string("_") + eps_string + std::string("_") + std::to_string(count)}; 
        h_el_energy_prob = new TH1D(h_name_prob.c_str(), "", num_bins, 0.0, 4.0);
        h_el_energy_prob->SetStats(0);
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
        }
        //std::cout << "likelihood = " << likelihood << std::endl;
        Double_t log_likelihood{std::log(likelihood)};
        vec_ll.push_back(-2.0 * log_likelihood);

        //std::cout << "LL=" << -2.0 * log_likelihood << ", LL2=" << -2.0 * likelihood_sum << ", diff=" << (-2.0 * log_likelihood) - (-2.0 * likelihood_sum) << std::endl;
        std::sort(likelihood_sum_vec.begin(), likelihood_sum_vec.end());
        Double_t likelihood_sum_vec_sum{0.0};
        for(std::vector<Double_t>::const_iterator it{likelihood_sum_vec.cbegin()}; it != likelihood_sum_vec.cend(); ++ it)
        {
            likelihood_sum_vec_sum += *it;
        }
        //std::cout << "LL3=" << -2.0 * likelihood_sum_vec_sum << std::endl;


        // create difference histogram, data - reweight
        std::string h_name_diff_data_rw{std::string("h_el_energy_diff_data_rw") + std::string("_") + h_name_append + eps_string + std::string("_") + std::to_string(count)};
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
        std::string h_name_diff_data_orig{std::string("h_el_energy_diff_data_orig") + std::string("_") + h_name_append + eps_string + std::string("_") + std::to_string(count)};
        h_el_energy_diff_data_orig = new TH1D(h_name_diff_data_orig.c_str(), "", num_bins, 0.0, 4.0);
        h_el_energy_diff_data_orig->SetStats(0);
        for(Int_t i{1}; i <= h_el_energy_diff_data_orig->GetNbinsX(); ++ i)
        {
            Double_t content1{h_el_energy_data->GetBinContent(i)};
            Double_t content2{h_el_energy_original->GetBinContent(i)};
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
    
        #if 0
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
        #endif
    
    }


    ////////////////////////////////////////////////////////////////////////////
    // GET DISTRIBUTION OF LL MEASUREMENTS
    ////////////////////////////////////////////////////////////////////////////

    //Double_t min{std::min_element(vec_ll.begin(), vec_ll.end())};
    //Double_t max{std::max_element(vec_ll.begin(), vec_ll.end())};
    std::pair<std::vector<Double_t>::iterator, std::vector<Double_t>::iterator> min_max_pair{std::minmax_element(vec_ll.begin(), vec_ll.end())};
    Double_t min{*min_max_pair.first};
    Double_t max{*min_max_pair.second};
    std::string name_ll{std::string("h_ll") + h_name_append + std::string("_") + eps_string};
    Double_t center{0.5 * (max + min)};
    Double_t width{0.5 * (max - min)};
    width *= 1.1;
    min = center - width;
    max = center + width;
    h_ll = new TH1D(name_ll.c_str(), "", 100, min, max);
    h_ll->GetXaxis()->SetTitle("Chi Square Value");
    h_ll->GetYaxis()->SetTitle("Number of Pseudo Experiments");
    // fill the histogram
    for(std::vector<Double_t>::const_iterator it{vec_ll.cbegin()}; it != vec_ll.cend(); ++ it)
    {
        h_ll->Fill(*it);
    }

    #if 0
    std::string c_name_ll{std::string("c_ll/") + std::string("c_ll_") + eps_string};
    c_ll = new TCanvas(c_name_ll.c_str(), "", 800, 600);
    h_ll->Draw("E");
    c_ll->SaveAs((c_name_ll + std::string(".png")).c_str());
    delete c_ll;
    #endif
        
}
