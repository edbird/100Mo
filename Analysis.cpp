// self header
#include "Analysis.hpp"


// local headers
#include "ReWeight.hpp"
#include "read_data.hpp"
#include "aux.hpp"


Analysis::Analysis(const std::string& filename, const std::string& output_filename)
    : epsilon_31{0.368}
    , f{nullptr}
    , t{nullptr}
    , num_bins{40}
    , filename{filename}
    , output_filename{output_filename}
    , _batch_mode_{false}
    , log_mode{false}
    , _gen_weight_enable_{false}
    , _energy_cut_enable_{false}
    , _canvas_enable_raw_data_{false}
    , _canvas_enable_decay_rate_{true}
    , _canvas_enable_single_electron_projection_{false}
    , g_el_energy_single_0{nullptr}
    , g_el_energy_single_1{nullptr}
    , g_el_energy_single_2{nullptr}
    , g_el_energy_sum_0{nullptr}
    , g_el_energy_sum_1{nullptr}
    , g_el_energy_sum_2{nullptr}
    , h_nEqNull{nullptr}
    , h_nEqTwo{nullptr}
    , psiN0{0.0}
    , psiN2{0.0}
    , sensitivity_chisquare{0.0}
    , sensitivity_chisquare_2d{0.0}
    , fit_subrange{false}
    , number_of_pseudo_experiments{1}
    , number_of_pseudo_experiments_2d{1}
    , h_data_0{nullptr}
    , h_data_1{nullptr}
    , h_data_2{nullptr}
    , h_single_electron_0{nullptr}
    , h_single_electron_1{nullptr}
    , h_single_electron_2{nullptr}
    , h_test_single_original{nullptr}
    , h_test_single_reweight{nullptr}
    , h_test_sum_original{nullptr}
    , h_test_sum_reweight{nullptr}
    , h_el_energy_original{nullptr}
    , h_el_energy_reweight{nullptr}
    , h_el_energy_sum_original{nullptr}
    , h_el_energy_sum_reweight{nullptr}
    , h_el_energy_2d_original{nullptr}
    , h_el_energy_2d_reweight{nullptr}
    , h_el_energy_2d_diff{nullptr}
    , h_el_energy_2d_pull{nullptr}
    , h_el_energy_2d_diff_data_rw{nullptr}
    , h_el_energy_2d_diff_data_orig{nullptr}
    , h_gen_weight{nullptr}
    , c_nEqNull{nullptr}
    , c_nEqTwo{nullptr}
    , c_data_0{nullptr}
    , c_data_1{nullptr}
    , c_data_2{nullptr}
    , c_single_electron{nullptr}
    , c_test_single{nullptr}
    , c_test_sum{nullptr}
    , c_gen_weight{nullptr}
    , c_el_energy_diff{nullptr}
    , c_el_energy_pull{nullptr}
    , c_el_energy_2d_original{nullptr}
    , c_el_energy_2d_reweight{nullptr}
    , c_el_energy_2d_diff{nullptr}
    , c_el_energy_2d_pull{nullptr}
    , c_ll{nullptr}
    , c_el_energy_2d_diff_data_rw{nullptr}
    , c_el_energy_2d_diff_data_orig{nullptr}
    , c_el_energy_2d_data{nullptr}
    , c_el_energy_2d_prob{nullptr}
    , c_ll_2d{nullptr}
{


    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////


}


Analysis::~Analysis()
{

}


////////////////////////////////////////////////////////////////////////////////
// ANALYSIS FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

void Analysis::SetEpsilon31(const Double_t epsilon)
{
    epsilon_31 = epsilon;
}

////////////////////////////////////////////////////////////////////////////////
// READ DATA IN FROM FILES
////////////////////////////////////////////////////////////////////////////////
void Analysis::ReadData()
{

    ////////////////////////////////////////////////////////////////////////////
    // READ DATA FOR PHASE SPACE FACTORS AND DECAY RATE DATA
    ////////////////////////////////////////////////////////////////////////////

    // read
    read_data_helper("nEqNull.dat", "nEqTwo.dat", "psiN0.txt", "psiN2.txt", data_nEqNull, data_nEqTwo, psiN0, psiN2);

    //std::cout << psiN0 << std::endl;
    //std::cout << psiN2 << std::endl;


    ////////////////////////////////////////////////////////////////////////////
    // READ DATA FOR ELECTRON ENERGY SINGLE AND SUM HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////

    // allocate data storage
    std::vector<std::vector<double>> data_electron_energy_single_0_0;
    std::vector<std::vector<double>> data_electron_energy_single_0_4;
    std::vector<std::vector<double>> data_electron_energy_single_0_8;
    std::vector<std::vector<double>> data_electron_energy_sum_0_0;
    std::vector<std::vector<double>> data_electron_energy_sum_0_4;
    std::vector<std::vector<double>> data_electron_energy_sum_0_8;

    // read
    read_data_helper_2("./data/mo100/single0.0.dat", data_electron_energy_single_0_0);
    read_data_helper_2("./data/mo100/single0.4.dat", data_electron_energy_single_0_4);
    read_data_helper_2("./data/mo100/single0.8.dat", data_electron_energy_single_0_8);
    read_data_helper_2("./data/mo100/sum0.0.dat", data_electron_energy_sum_0_0);
    read_data_helper_2("./data/mo100/sum0.4.dat", data_electron_energy_sum_0_4);
    read_data_helper_2("./data/mo100/sum0.8.dat", data_electron_energy_sum_0_8);


    ////////////////////////////////////////////////////////////////////////////
    // CONVERT INPUT DATA TO GRAPH FORMAT
    ////////////////////////////////////////////////////////////////////////////

    Int_t data_size{(Int_t)data_electron_energy_single_0_0.size()};

    if(data_electron_energy_single_0_4.size() != data_size)
        throw "error: size";
    
    if(data_electron_energy_single_0_8.size() != data_size)
        throw "error: size";
    
    if(data_electron_energy_sum_0_0.size() != data_size)
        throw "error: size";
    
    if(data_electron_energy_sum_0_4.size() != data_size)
        throw "error: size";

    if(data_electron_energy_sum_0_8.size() != data_size)
        throw "error: size";

    Double_t *data_x = new Double_t[data_size];

    Double_t *data_single_0 = new Double_t[data_size];
    Double_t *data_single_1 = new Double_t[data_size];
    Double_t *data_single_2 = new Double_t[data_size];

    Double_t *data_sum_0 = new Double_t[data_size];
    Double_t *data_sum_1 = new Double_t[data_size];
    Double_t *data_sum_2 = new Double_t[data_size];

    // SCALE THE DATA SO THAT X RUNS FROM 0.0 TO Q_BB RATHER THAN 0.0 TO 1.0
    // SCALE THE DATA SO THAT THE INTEGRAL IS UNCHANGED, IE: INTEGRAL = 1.0
    // THIS IS DONE BY SCALING THE Y VALUES IN THE INVERSE WAY TO THE SCALING
    // OF THE X VALUES

    for(std::size_t i{0}; i < data_size; ++ i)
        data_x[i] = bb_Q * data_electron_energy_single_0_0[i][0];

    for(std::size_t i{0}; i < data_size; ++ i)
        data_single_0[i] = (1.0 / bb_Q) * data_electron_energy_single_0_0[i][1];

    for(std::size_t i{0}; i < data_size; ++ i)
        data_single_1[i] = (1.0 / bb_Q) * data_electron_energy_single_0_4[i][1];
    
    for(std::size_t i{0}; i < data_size; ++ i)
        data_single_2[i] = (1.0 / bb_Q) * data_electron_energy_single_0_8[i][1];

    for(std::size_t i{0}; i < data_size; ++ i)
        data_sum_0[i] = (1.0 / bb_Q) * data_electron_energy_sum_0_0[i][1];

    for(std::size_t i{0}; i < data_size; ++ i)
        data_sum_1[i] = (1.0 / bb_Q) * data_electron_energy_sum_0_4[i][1];
    
    for(std::size_t i{0}; i < data_size; ++ i)
        data_sum_2[i] = (1.0 / bb_Q) * data_electron_energy_sum_0_8[i][1];

    g_el_energy_single_0 = new TGraph(data_size, data_x, data_single_0);
    g_el_energy_single_0->SetLineColor(2);

    g_el_energy_single_1 = new TGraph(data_size, data_x, data_single_1);
    g_el_energy_single_1->SetLineColor(3);
    
    g_el_energy_single_2 = new TGraph(data_size, data_x, data_single_2);
    g_el_energy_single_2->SetLineColor(4);

    g_el_energy_sum_0 = new TGraph(data_size, data_x, data_sum_0);
    g_el_energy_sum_0->SetLineColor(2);
    
    g_el_energy_sum_1 = new TGraph(data_size, data_x, data_sum_1);
    g_el_energy_sum_1->SetLineColor(3);
    
    g_el_energy_sum_2 = new TGraph(data_size, data_x, data_sum_2);
    g_el_energy_sum_2->SetLineColor(4);

    // compute integral to test
    // trapezium rule
    /*
    Double_t data_single_0_integral{0.0};
    for(Int_t i{0 + 1}; i < data_size - 1; ++ i)
    {
        data_single_0_integral += data_single_0[i];
    }
    data_single_0_integral *= 2.0;
    data_single_0_integral += data_single_0[0];
    data_single_0_integral += data_single_0[data_size - 1];
    data_single_0_integral *= (data_x[1] - data_x[0]) / 2.0;
    std::cout << "data_single_0_integral=" << data_single_0_integral << std::endl;
    */

    ////////////////////////////////////////////////////////////////////////////
    // TODO: apply phase space factor to data
    ////////////////////////////////////////////////////////////////////////////

    // ...
    // Note: Do not do this

    ////////////////////////////////////////////////////////////////////////////
    // CONVERT INPUT DATA TO HISTOGRAM FORMAT
    // Note: Added 2018-04-23 (After INTERMEDIATE DATA below)
    // Note: These histograms do NOT have the phase space variable included
    ////////////////////////////////////////////////////////////////////////////

    // don't plot raw data
    //TH2D *h_nEqNull = nullptr; //= new TH2D("h_nEqNull", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_nEqTwo = nullptr; //= new TH2D("h_nEqTwo", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_ratio = new TH2D("h_ratio", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //h_nEqNull->SetStats(0);
    //h_nEqTwo->SetStats(0);
    //h_ratio->SetStats(0);
    convert_data_to_histogram_format(data_nEqNull, data_nEqTwo, h_nEqNull, h_nEqTwo);
    // TODO: move code from this function including alloc. of histograms
    // into main code

    ////////////////////////////////////////////////////////////////////////////
    // CREATE INTERMEDIATE DATA
    // Format: Nx3 array, each array a different value of epsilon
    ////////////////////////////////////////////////////////////////////////////

    // create data array for "complete data"
    // (with phase space factors)
    
    // epsilon_31 = 0.0
    create_data_with_phase_space_factor(data_0, 0.0, data_nEqNull, data_nEqTwo, psiN0, psiN2);

    // epsilon_31 = 0.4
    create_data_with_phase_space_factor(data_1, 0.4, data_nEqNull, data_nEqTwo, psiN0, psiN2);
    
    // epsilon_31 = 0.8
    create_data_with_phase_space_factor(data_2, 0.8, data_nEqNull, data_nEqTwo, psiN0, psiN2);

    // create data array for ratio
    /*
    std::vector<std::vector<double>> ratio;
    ratio.resize(data_nEqNull.size());
    for(std::size_t i{0}; i < ratio.size(); ++ i)
    {
        ratio[i].resize(data_nEqNull[i].size());
    }

    // epsilon_31 = 0.8
    const double epsilon_31_2{0.8};
    for(std::size_t i{0}; i < ratio.size(); ++ i)
    {
        for(std::size_t j{0}; j < 2 /*ratio[i].size()*//*; ++ j)
        {
            ratio[i][j] = data_nEqNull[i][j];
        }
        ratio[i][2] = (psiN0 / (psiN0 + epsilon_31 * psiN2)) * ((data_nEqNull[i][2] + epsilon_31 * data_nEqTwo[i][2]) / data_nEqNull[i][2]);
    }
    std::cout << "Finished constructing ratio" << std::endl;
    */

    std::cout << "Finished constructing intermediate data" << std::endl;

    //std::cout << ratio << std::endl;


}


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
            TCanvas *c_nEqNull = new TCanvas("c_nEqNull", "", 4000, 3000);
            h_nEqNull->Draw("colz");
            c_nEqNull->SaveAs("c_nEqNull.png");
            c_nEqNull->SaveAs("c_nEqNull.pdf");
            c_nEqNull->SaveAs("c_nEqNull.C");
            delete c_nEqNull;

            TCanvas *c_nEqTwo = new TCanvas("c_nEqTwo", "", 4000, 3000);
            h_nEqTwo->Draw("colz");
            c_nEqTwo->SaveAs("c_nEqTwo.png");
            c_nEqTwo->SaveAs("c_nEqTwo.pdf");
            c_nEqTwo->SaveAs("c_nEqTwo.C");
            delete c_nEqTwo;
        }


        if(_canvas_enable_decay_rate_)
        {
            TCanvas *c_data_0 = new TCanvas("c_data_0", "", 4000, 3000);
            h_data_0->Draw("colz");
            c_data_0->SaveAs("c_data_0.png");
            c_data_0->SaveAs("c_data_0.pdf");
            c_data_0->SaveAs("c_data_0.C");
            delete c_data_0;

            TCanvas *c_data_1 = new TCanvas("c_data_1", "", 4000, 3000);
            h_data_1->Draw("colz");
            c_data_1->SaveAs("c_data_1.png");
            c_data_1->SaveAs("c_data_1.pdf");
            c_data_1->SaveAs("c_data_1.C");
            delete c_data_1;

            TCanvas *c_data_2 = new TCanvas("c_data_2", "", 4000, 3000);
            h_data_2->Draw("colz");
            c_data_2->SaveAs("c_data_2.png");
            c_data_2->SaveAs("c_data_2.pdf");
            c_data_2->SaveAs("c_data_2.C");
            delete c_data_2;
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
            TCanvas *c_single_electron = new TCanvas("c_single_electron", "", 4000, 3000);
            h_single_electron_0->SetMaximum(3.5);
            h_single_electron_0->GetXaxis()->SetTitle("Energy [MeV]");
            h_single_electron_0->GetYaxis()->SetTitle("Events");
            h_single_electron_0->Draw("hist");
            h_single_electron_1->Draw("histsame");
            h_single_electron_2->Draw("histsame");
            c_single_electron->SaveAs("c_single_electron.png");
            c_single_electron->SaveAs("c_single_electron.pdf");
            c_single_electron->SaveAs("c_single_electron.C");
            delete c_single_electron;
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

void Analysis::CanvasSingleElectronTest()
{

    if(_batch_mode_ == false)
    {
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
        //std::cout << "chi1: " << chi_square_test(g_el_energy_single_0, h_test_single_original) / h_test_single_original->GetNbinsX() << std::endl;
        //std::cout << "chi2: " << chi_square_test(g_el_energy_single_2, h_test_single_reweight) / h_test_single_reweight->GetNbinsX() << std::endl;
        
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
        //std::cout << "chi3: " << chi_square_test(g_el_energy_sum_0, h_test_sum_original) / h_test_sum_original->GetNbinsX() << std::endl;
        //std::cout << "chi4: " << chi_square_test(g_el_energy_sum_2, h_test_sum_reweight) / h_test_sum_reweight->GetNbinsX() << std::endl;

    }
}


////////////////////////////////////////////////////////////////////////////////
// CREATE HISTOGRAMS REQUIRED FOR EVENT LOOP
////////////////////////////////////////////////////////////////////////////////

void Analysis::InitEventLoop()
{

    ////////////////////////////////////////////////////////////////////////////
    // HISTOGRAMS FOR SINGLE AND SUM ENERGY HISTOGRAMS FOR BASELINE (EPS=BASE)
    // AND RESCALED HISTOGRAMS (EPS=SOME OTHER VALUE)
    ////////////////////////////////////////////////////////////////////////////
    
    // load from NEMO-3 data (MC), the reconstructed electron energy (2 in same histogram)
    h_el_energy_original = new TH1D("h_el_energy_original", "", num_bins, 0.0, 4.0);
    //h_el_energy_original->SetStats(0);
    //h_el_energy_original->SetLineColor(2);
    //h_el_energy_original->SetMarkerColor(2);
    // same as above but re-weighted
    h_el_energy_reweight = new TH1D("h_el_energy_reweight", "", num_bins, 0.0, 4.0);
    //h_el_energy_reweight->SetStats(0);
    //h_el_energy_reweight->SetLineColor(4);
    //h_el_energy_reweight->SetMarkerColor(4);

    // load from NEMO-3 data (mc), the sum of the two reconstructed electron energies
    h_el_energy_sum_original = new TH1D("h_el_energy_sum_original", "", num_bins, 0.0, 4.0);
    h_el_energy_sum_original->SetStats(0);
    h_el_energy_sum_original->SetLineColor(2);
    h_el_energy_sum_original->SetMarkerColor(2);
    // same as above but re-weighted
    h_el_energy_sum_reweight = new TH1D("h_el_energy_sum_reweight", "", num_bins, 0.0, 4.0);
    h_el_energy_sum_reweight->SetStats(0);
    h_el_energy_sum_reweight->SetLineColor(4);
    h_el_energy_sum_reweight->SetMarkerColor(4);

    // 2d fit version of h_el_energy_sum_original
    // instead of using a 1d histogram containing the sum of the 2 electron
    // energies, use a 2d histogram containing both energies
    // TODO: should I sort by energy?
    h_el_energy_2d_original = new TH2D("h_el_energy_2d_original", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_el_energy_2d_original->SetStats(0);
    
    h_el_energy_2d_reweight = new TH2D("h_el_energy_2d_reweight", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_el_energy_2d_reweight->SetStats(0);

    // chisquare method
    h_el_energy_2d_diff = new TH2D("h_el_energy_2d_diff", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_el_energy_2d_diff->SetStats(0);
    // chisquare method
    h_el_energy_2d_pull = new TH2D("h_el_energy_2d_pull", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_el_energy_2d_pull->SetStats(0);

// NOTE: defined when required in code below
    h_el_energy_2d_diff_data_rw = new TH2D("h_el_energy_2d_diff_data_rw", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_el_energy_2d_diff_data_rw->SetStats(0);
    h_el_energy_2d_diff_data_orig = new TH2D("h_el_energy_2d_diff_data_orig", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_el_energy_2d_diff_data_orig->SetStats(0);
//#if TWO_D_METHOD_LOGLIKELIHOOD
//    // pseudodata histogram
//    TH2D *h_el_energy_2d_data = new TH2D("h_el_energy_2d_data", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
//    h_el_energy_2d_data->SetStats(0);
//
//    TH2D *h_el_energy_2d_prob = new TH2D("h_el_energy_2d_prob", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
//    h_el_energy_2d_prob->SetStats(0);
//#endif

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



    // NEMO-3 data/MC read from tree
    // input file

    TH2D *h_gen_weight = new TH2D("h_gen_weight", "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
    h_gen_weight->SetStats(0);


    ////////////////////////////////////////////////////////////////////////////
    // initialize file and branch addresses
    ////////////////////////////////////////////////////////////////////////////

    // read file using program arguments
    f = new TFile(filename.c_str());
    t = (TTree*)f->Get("NewElectronNtuplizer/NewElectronNtuplizer");

    t->SetBranchAddress("nElectrons", &nElectrons);
    t->SetBranchAddress("trueT1", &trueT1);
    t->SetBranchAddress("trueT2", &trueT2);
    t->SetBranchAddress("el_energy_", el_energy_);
    // new method requires weight to be saved to tree
    gen_weight = 1.0;
    if(_gen_weight_enable_ == true)
    {
        t->SetBranchAddress("gen_weight", &gen_weight);
    }


}


////////////////////////////////////////////////////////////////////////////////
// EVENT LOOP
////////////////////////////////////////////////////////////////////////////////

void Analysis::EventLoop()
{

    std::cout << "Processing data" << std::endl;
    Long64_t prog_c{-1};
    //const Double epsilon_31{0.0};
    for(Long64_t ix{0}; ix < t->GetEntries(); ++ ix)
    {

        // TODO remove
        //if(ix == 1000) break;

        t->GetEntry(ix);

        if(nElectrons != 2) continue;

        // cut both electrons > 300 keV
        if(_energy_cut_enable_)
        {
            if(el_energy_[0] < 0.3) continue;
            if(el_energy_[1] < 0.3) continue;
        }

        Double_t T1{trueT1 / bb_Q};
        Double_t T2{trueT2 / bb_Q};

        // ReWeight = baseline 0.0, ReWeight2 = baseline = 0.382
        Double_t weight{ReWeight2(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")};

        h_el_energy_original->Fill(el_energy_[0], 1.0 * gen_weight);
        h_el_energy_original->Fill(el_energy_[1], 1.0 * gen_weight);
        
        h_el_energy_sum_original->Fill(el_energy_[0] + el_energy_[1], 1.0 * gen_weight);

        if(el_energy_[0] <= el_energy_[1])
        {
            h_el_energy_2d_original->Fill(el_energy_[0], el_energy_[1], 1.0 * gen_weight);
        }
        else
        {
            h_el_energy_2d_original->Fill(el_energy_[1], el_energy_[0], 1.0 * gen_weight);
        }

        // test histograms
        h_test_single_original->Fill(trueT1, 1.0 * gen_weight);
        h_test_single_original->Fill(trueT2, 1.0 * gen_weight);

        h_test_sum_original->Fill(trueT1 + trueT2, 1.0 * gen_weight);

        /*
        if(!(weight < 3.0 && weight > 0.001))
        {
            std::cout << "ix=" << ix << std::endl;
            std::cout << "trueT1=" << trueT1 << " trueT2=" << trueT2 << " sum=" << trueT1 + trueT2 << std::endl;

            ReWeight(T1, T2, epsilon_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true");
            
            std::cout << "weight=" << weight << std::endl;
            std::cin.get();
        }
        */

        h_el_energy_reweight->Fill(el_energy_[0], weight * gen_weight);
        h_el_energy_reweight->Fill(el_energy_[1], weight * gen_weight);

        h_el_energy_sum_reweight->Fill(el_energy_[0] + el_energy_[1], weight * gen_weight);

        if(el_energy_[0] <= el_energy_[1])
        {
            h_el_energy_2d_reweight->Fill(el_energy_[0], el_energy_[1], weight * gen_weight);
        }
        else
        {
            h_el_energy_2d_reweight->Fill(el_energy_[1], el_energy_[0], weight * gen_weight);
        }

        // test histograms
        h_test_single_reweight->Fill(trueT1, weight * gen_weight);
        h_test_single_reweight->Fill(trueT2, weight * gen_weight);

        h_test_sum_reweight->Fill(trueT1 + trueT2, weight * gen_weight);

        h_gen_weight->Fill(trueT1, trueT2, gen_weight);

        //if(ix % 10 == 9) std::cin.get();

        if((ix * 10) / t->GetEntries() > prog_c)
        {
            prog_c = (ix * 10) / t->GetEntries();
            //std::cout << 10 * (ix * 10) / t->GetEntries() << " %" << std::endl;
        }

    }

}


////////////////////////////////////////////////////////////////////////////////
// POST PROCESSING
////////////////////////////////////////////////////////////////////////////////

void Analysis::PostProcess()
{

    
    // rescale histograms pre other re-scaling
    // NOTE: to get the same >>> amplitude <<<
    // histograms must be scaled by the bin width
    // the bin width for one set is 100 keV
    // the bin width for the other set is 1 MeV / 1000 = 1 keV (?)
    // factor 1000 discrepency?
    // TODO: resolve this issue
    //h_el_energy_original->Scale(1.0 / h_el_energy_original->Integral());
    //h_el_energy_reweight->Scale(1.0 / h_el_energy_reweight->Integral());
    //h_el_energy_sum_original->Scale(1.0 / h_el_energy_sum_original->Integral());
    //h_el_energy_sum_reweight->Scale(1.0 / h_el_energy_sum_reweight->Integral());

    // NOTE: must scale all by same amount here
    // TODO: removed scaling
    /*
    Double_t scale_factor{1.0 / h_el_energy_sum_original->Integral()};
    h_el_energy_original->Scale(scale_factor);
    h_el_energy_reweight->Scale(scale_factor);
    h_el_energy_sum_original->Scale(scale_factor);
    h_el_energy_sum_reweight->Scale(scale_factor);
    */
    
    h_test_single_original->Scale((1.0 / 0.1) / h_test_single_original->Integral());
    h_test_single_reweight->Scale((1.0 / 0.1) / h_test_single_reweight->Integral());
    h_test_sum_original->Scale((1.0 / 0.1) / h_test_sum_original->Integral());
    h_test_sum_reweight->Scale((1.0 / 0.1) / h_test_sum_reweight->Integral());

    std::cout << "h_el_energy_original_integral=" << h_el_energy_original->Integral() << std::endl;

    // calculate integral manually
    Double_t h_el_energy_original_integral{0.0};
    for(Int_t i{1}; i <= h_el_energy_original->GetNbinsX(); ++ i)
    {
        h_el_energy_original_integral += h_el_energy_original->GetBinContent(i);
    }
    std::cout << "h_el_energy_original_integral=" << h_el_energy_original_integral << std::endl;

    if(_batch_mode_ == false)
    {
        TCanvas *c_gen_weight = new TCanvas("c_gen_weight", "", 800, 600);
        h_gen_weight->Draw("colz");
        c_gen_weight->SaveAs("c_gen_weight.C");
        c_gen_weight->SaveAs("c_gen_weight.png");
        c_gen_weight->SaveAs("c_gen_weight.pdf");
        delete c_gen_weight;
    }


}


////////////////////////////////////////////////////////////////////////////////
// SCALE SUMMED HISTOGRAM TO MATCH AMPLITUDE
////////////////////////////////////////////////////////////////////////////////

void Analysis::SummedEnergyFit()
{

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
    std::string dir;
    if(log_mode) { min = 0.1; max = max_log_mode_2; dir = dir_log_mode_2; }
    else { min = 0.0; max = max_nolog_mode_2; dir = dir_nolog_mode_2; }
    CanvasFactorySettings settings_2("Energy [MeV]", "Events", min, max, log_mode);
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
    for(Int_t i{1}; i <= h_el_energy_pull->GetNbinsX(); ++ i)
    {
        Double_t content{h_el_energy_diff->GetBinContent(i)};
        Double_t error{h_el_energy_diff->GetBinError(i)};
        if(error == 0.0)
        {
            //h_el_energy_pull->SetBinContent(i, j, 0.0);
        }
        else
        {
            h_el_energy_pull->SetBinContent(i, content / error);
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
        std::string dir;
        if(log_mode) { min = 0.1; max = max_log_mode; dir = dir_log_mode; }
        else { min = 0.0; max = max_nolog_mode; dir = dir_nolog_mode; }
        CanvasFactorySettings settings("Energy [MeV]", "Events", min, max, log_mode);
        settings.SetDrawOption("E");
        CanvasFactory factory(settings);
        factory.Canvas("el_energy", dir, h_el_energy_original, "Baseline", h_el_energy_reweight, "Reweighted");

        TCanvas *c_el_energy_diff = new TCanvas("c_el_energy_diff", "", 800, 600);
        //c_el_energy_diff->SetRightMargin(0.12);
        //c_el_energy_diff->SetLogz();
        h_el_energy_diff->GetXaxis()->SetTitle("Energy [MeV]");
        h_el_energy_diff->GetYaxis()->SetTitle("Events");
        h_el_energy_diff->Draw("E");
        c_el_energy_diff->SaveAs("c_el_energy_diff.C");
        c_el_energy_diff->SaveAs("c_el_energy_diff.png");
        c_el_energy_diff->SaveAs("c_el_energy_diff.pdf");
        
        TCanvas *c_el_energy_pull = new TCanvas("c_el_energy_pull", "", 800, 600);
        c_el_energy_pull->SetRightMargin(0.12);
        //c_el_energy_pull->SetLogz();
        h_el_energy_pull->GetXaxis()->SetTitle("Energy [MeV]");
        h_el_energy_pull->GetYaxis()->SetTitle("Events");
        h_el_energy_pull->Draw("E");
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


    // NOTE TO SELF: There is no 2d fit, the 2d histogram is used to evaluate
    // the sensitivity, therefore there is only a chi-square test

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
    //}
    //else if(fit_subrange == true)
    //{
    //    sensitivity_chisquare_2d = chi_square_test(h_el_energy_reweight, h_el_energy_original, 2.0, 4.0);
    //    std::cout << "chi square of single electron, 2.0 MeV - 4.0 MeV: " << sensitivity_chisquare_2d << std::endl;
    //    std::cout << " degrees of freedom: " << non_empty_bins_2d << std::endl;
    //    std::cout << " chi square reduced: " << sensitivity_chisquare_2d / (Double_t)non_empty_bins_2d << std::endl;
    //}



    if(_batch_mode_ == false)
    {

        ////////////////////////////////////////////////////////////////////////
        // 2D ORIGINAL AND REWEIGHT CANVAS OUTPUT
        ////////////////////////////////////////////////////////////////////////

        // 2d original and reweight histograms
        TCanvas *c_el_energy_2d_original = new TCanvas("c_el_energy_2d_original", "", 800, 600);
        c_el_energy_2d_original->SetRightMargin(0.15);
        //c_el_energy_2d_original->SetLogz();
        h_el_energy_2d_original->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
        h_el_energy_2d_original->GetYaxis()->SetTitle("High Energy Electron [MeV]");
        h_el_energy_2d_original->Draw("colz");
        c_el_energy_2d_original->SaveAs("c_el_energy_2d_original.C");
        c_el_energy_2d_original->SaveAs("c_el_energy_2d_original.png");
        c_el_energy_2d_original->SaveAs("c_el_energy_2d_original.pdf");

        TCanvas *c_el_energy_2d_reweight = new TCanvas("c_el_energy_2d_reweight", "", 800, 600);
        c_el_energy_2d_reweight->SetRightMargin(0.15);
        //c_el_energy_2d_reweight->SetLogz();
        h_el_energy_2d_reweight->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
        h_el_energy_2d_reweight->GetYaxis()->SetTitle("High Energy Electron [MeV]");
        h_el_energy_2d_reweight->Draw("colz");
        c_el_energy_2d_reweight->SaveAs("c_el_energy_2d_reweight.C");
        c_el_energy_2d_reweight->SaveAs("c_el_energy_2d_reweight.png");
        c_el_energy_2d_reweight->SaveAs("c_el_energy_2d_reweight.pdf");


        ////////////////////////////////////////////////////////////////////////
        // CREATE DIFFERENCE AND PULL HISTOGRAMS AND CANVAS OUTPUT
        ////////////////////////////////////////////////////////////////////////

        // create difference histogram
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
        TCanvas *c_el_energy_2d_diff = new TCanvas("c_el_energy_2d_diff", "", 800, 600);
        c_el_energy_2d_diff->SetRightMargin(0.12);
        //c_el_energy_2d_diff->SetLogz();
        h_el_energy_2d_diff->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
        h_el_energy_2d_diff->GetYaxis()->SetTitle("High Energy Electron [MeV]");
        h_el_energy_2d_diff->Draw("colz");
        c_el_energy_2d_diff->SaveAs("c_el_energy_2d_diff.C");
        c_el_energy_2d_diff->SaveAs("c_el_energy_2d_diff.png");
        c_el_energy_2d_diff->SaveAs("c_el_energy_2d_diff.pdf");
        
        // create pull histogram
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
        TCanvas *c_el_energy_2d_pull = new TCanvas("c_el_energy_2d_pull", "", 800, 600);
        c_el_energy_2d_pull->SetRightMargin(0.12);
        //c_el_energy_2d_pull->SetLogz();
        h_el_energy_2d_pull->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
        h_el_energy_2d_pull->GetYaxis()->SetTitle("High Energy Electron [MeV]");
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

    for(Int_t count{0}; count < number_of_pseudo_experiments; ++ count)
    {
        ////////////////////////////////////////////////////////////////////////
        // SINGLE ELECTRON ENERGY PSEUDODATA METHOD
        ////////////////////////////////////////////////////////////////////////

        {
            // TODO: name of histograms which we iterate over
            // log likelihood method
            // create "data" histogram - poisson generated data
            std::string h_name{std::string("h_el_energy_data_") + std::to_string(count)};
            TH1I *h_el_energy_data = new TH1I(h_name.c_str(), "", num_bins, 0.0, 4.0);
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
            std::string h_name_prob{std::string("h_el_energy_prob_") + std::to_string(count)}; 
            TH1D *h_el_energy_prob = new TH1D(h_name_prob.c_str(), "", num_bins, 0.0, 4.0);
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

                h_el_energy_prob->SetBinContent(ix, poi);
            }
            //std::cout << "likelihood = " << likelihood << std::endl;
            Double_t log_likelihood{std::log(likelihood)};
            vec_ll.push_back(-2.0 * log_likelihood);


            ////////////////////////////////////////////////////////////////////
            // CANVAS OUTPUT
            ////////////////////////////////////////////////////////////////////

            //const Double_t canvas_max{2.5e5};
            const Double_t canvas_max{1.0e6};
            //const Double_t canvas_min{0.0};
            const Double_t canvas_min{1.0e-1};
            const std::string canvas_dir(".");
            CanvasFactorySettings settings("Energy [MeV]", "Events", canvas_min, canvas_max, false);
            settings.SetLogMode(true);
            settings.SetDrawOption("E");
            CanvasFactory factory(settings);
            factory.Canvas("el_energy_data", canvas_dir, h_el_energy_original, "Baseline", h_el_energy_reweight, "Reweighted", h_el_energy_data, "Pseudodata");

            settings.SetMax(1.1);
            settings.SetMin(1.0e-5); // 1.0e-35
            settings.SetDrawOption("hist");
            settings.SetLogMode(true);
            factory.Settings(settings);
            factory.Canvas("el_energy_prob", canvas_dir, h_el_energy_prob, "Probability");
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
    TH1D *h_ll = new TH1D("h_ll", "", 100, min, max);
    h_ll->GetXaxis()->SetTitle("Log Likelihood Value");
    h_ll->GetYaxis()->SetTitle("Number of Pseudo Experiments");
    // fill the histogram
    for(std::vector<Double_t>::const_iterator it{vec_ll.cbegin()}; it != vec_ll.cend(); ++ it)
    {
        h_ll->Fill(*it);
    }

    TCanvas *c_ll = new TCanvas("c_ll", "", 800, 600);
    h_ll->Draw("E");
    c_ll->SaveAs("c_ll.png");
    delete c_ll;
        
}


////////////////////////////////////////////////////////////////////////////////
// 2D SINGLE ELECTRON ENERGY LOG LIKELIHOOD SENSITIVITY
////////////////////////////////////////////////////////////////////////////////

void Analysis::SensitivityMeasurementLoglikelihood2()
{

    for(Int_t count{0}; count < number_of_pseudo_experiments; ++ count)
    {
        ////////////////////////////////////////////////////////////////////////
        // INDEPENDENT SINGLE ELECTRON ENERGY
        ////////////////////////////////////////////////////////////////////////

        {
            // log likelihood method
            // create 2d "data" histogram - poisson generated data

            // chisquare method
            // create 2d "data" histogram
            // TODO which one

            std::string h_name_2d{std::string("h_el_energy_2d_data_") + std::to_string(count)};
            TH2I *h_el_energy_2d_data = new TH2I(h_name_2d.c_str(), "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
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
            TH2D *h_el_energy_2d_prob = new TH2D(h_name_prob_2d.c_str(), "", num_bins, 0.0, 4.0, num_bins, 0.0, 4.0);
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
                    if(poi != 1.0)
                    {
                        if(poi < 1.0e-14 || poi > 1.0)
                        {
                        //    std::cout << ix << ", " << jx << " -> " << poi << " data=" << data << " lambda=" << lambda << std::endl;
                        }
                    }
                }
            }
            std::cout << "likelihood (2d) = " << likelihood_2d << std::endl;
            Double_t log_likelihood_2d{std::log(likelihood_2d)};
            vec_ll_2d.push_back(-2.0 * log_likelihood_2d);
        
            
            // create difference histogram - difference between data and reweighted
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
            for(Int_t j{1}; j <= h_el_energy_2d_diff_data_orig->GetNbinsY(); ++ j)
            {
                for(Int_t i{1}; i <= h_el_energy_2d_diff_data_orig->GetNbinsX(); ++ i)
                {
                    Double_t content1{h_el_energy_2d_data->GetBinContent(i, j)};
                    Double_t content2{h_el_energy_2d_reweight->GetBinContent(i, j)};
                    Double_t error1{h_el_energy_2d_data->GetBinError(i, j)};
                    // note: do not use error or reweighted
                    Double_t error2{0.0 * h_el_energy_2d_reweight->GetBinError(i, j)};
                    Double_t content{content1 - content2};
                    Double_t error{std::sqrt(error1 * error1 + error2 * error2)};
                    h_el_energy_2d_diff_data_orig->SetBinContent(i, j, content);
                    h_el_energy_2d_diff_data_orig->SetBinError(i, j, error);
                }
            }

            ////////////////////////////////////////////////////////////////////
            // CANVAS OUTPUT
            ////////////////////////////////////////////////////////////////////

            // TODO: multiple experiemnts problem
            // if statement for batch mode
            TCanvas *c_el_energy_2d_data = new TCanvas("c_el_energy_2d_data", "", 800, 600);
            c_el_energy_2d_data->SetRightMargin(0.12);
            //c_el_energy_2d_data->SetLogz();
            h_el_energy_2d_data->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
            h_el_energy_2d_data->GetYaxis()->SetTitle("High Energy Electron [MeV]");
            h_el_energy_2d_data->Draw("colz");
            c_el_energy_2d_data->SaveAs("c_el_energy_2d_data.C");
            c_el_energy_2d_data->SaveAs("c_el_energy_2d_data.png");
            c_el_energy_2d_data->SaveAs("c_el_energy_2d_data.pdf");

            TCanvas *c_el_energy_2d_prob = new TCanvas("c_el_energy_2d_prob", "", 800, 600);
            c_el_energy_2d_prob->SetRightMargin(0.12);
            c_el_energy_2d_prob->SetLogz();
            h_el_energy_2d_prob->SetMinimum(1.0e-4);
            h_el_energy_2d_prob->SetMaximum(1.0e0);
            h_el_energy_2d_prob->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
            h_el_energy_2d_prob->GetYaxis()->SetTitle("High Energy Electron [MeV]");
            h_el_energy_2d_prob->Draw("colz");
            c_el_energy_2d_prob->SaveAs("c_el_energy_2d_prob.C");
            c_el_energy_2d_prob->SaveAs("c_el_energy_2d_prob.png");
            c_el_energy_2d_prob->SaveAs("c_el_energy_2d_prob.pdf");

            TCanvas *c_el_energy_2d_diff_data_rw = new TCanvas("c_el_energy_2d_diff_data_rw", "", 800, 600);
            c_el_energy_2d_diff_data_rw->SetRightMargin(0.12);
            //c_el_energy_2d_diff_data_rw->SetLogz();
            h_el_energy_2d_diff_data_rw->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
            h_el_energy_2d_diff_data_rw->GetYaxis()->SetTitle("High Energy Electron [MeV]");
            h_el_energy_2d_diff_data_rw->Draw("colz");
            c_el_energy_2d_diff_data_rw->SaveAs("c_el_energy_2d_diff_data_rw.C");
            c_el_energy_2d_diff_data_rw->SaveAs("c_el_energy_2d_diff_data_rw.png");
            c_el_energy_2d_diff_data_rw->SaveAs("c_el_energy_2d_diff_data_rw.pdf");

            TCanvas *c_el_energy_2d_diff_data_orig = new TCanvas("c_el_energy_2d_diff_data_orig", "", 800, 600);
            c_el_energy_2d_diff_data_orig->SetRightMargin(0.12);
            //c_el_energy_2d_diff_data_orig->SetLogz();
            h_el_energy_2d_diff_data_orig->GetXaxis()->SetTitle("Low Energy Electron [MeV]");
            h_el_energy_2d_diff_data_orig->GetYaxis()->SetTitle("High Energy Electron [MeV]");
            h_el_energy_2d_diff_data_orig->Draw("colz");
            c_el_energy_2d_diff_data_orig->SaveAs("c_el_energy_2d_diff_data_orig.C");
            c_el_energy_2d_diff_data_orig->SaveAs("c_el_energy_2d_diff_data_orig.png");
            c_el_energy_2d_diff_data_orig->SaveAs("c_el_energy_2d_diff_data_orig.pdf");


            ////////////////////////////////////////////////////////////////////////
            // SAVE TO ROOT FILE
            ////////////////////////////////////////////////////////////////////////

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

        }
    }


    ////////////////////////////////////////////////////////////////////////////
    // GET DISTRIBUTION OF LL MEASUREMENTS
    ////////////////////////////////////////////////////////////////////////////

    std::pair<std::vector<Double_t>::iterator, std::vector<Double_t>::iterator> min_max_pair_2d{std::minmax_element(vec_ll_2d.begin(), vec_ll_2d.end())};
    Double_t min_2d{*min_max_pair_2d.first};
    Double_t max_2d{*min_max_pair_2d.second};
    TH1D *h_ll_2d = new TH1D("h_ll_2d", "", 100, min_2d, max_2d);
    h_ll_2d->GetXaxis()->SetTitle("Log Likelihood Value");
    h_ll_2d->GetYaxis()->SetTitle("Number of Pseudo Experiments");
    // fill the histogram
    for(std::vector<Double_t>::const_iterator it{vec_ll_2d.cbegin()}; it != vec_ll_2d.cend(); ++ it)
    {
        h_ll_2d->Fill(*it);
    }

    
    TCanvas *c_ll_2d = new TCanvas("c_ll_2d", "", 800, 600);
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



////////////////////////////////////////////////////////////////////////////////
// FLAG FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

void Analysis::SetBatchMode(const bool mode)
{
    _batch_mode_ = mode;
}

void Analysis::SetLogMode(const bool mode)
{
    log_mode = mode;
}

void Analysis::SetEnergyCutEnabled(const bool enabled)
{
    _energy_cut_enable_ = enabled;
}

void Analysis::SetGenWeightEnabled(const bool enabled)
{
    _gen_weight_enable_ = enabled;
}

void Analysis::SetCanvasEnableRawData(const bool mode)
{
    _canvas_enable_raw_data_ = mode;
}

void Analysis::SetCanvasEnableDecayRate(const bool mode)
{
    _canvas_enable_decay_rate_ = mode;
}



void Analysis::SetFitSubrange(const bool flag)
{
    fit_subrange = flag;
}


void Analysis::SetNumberPseudoExperiments(const Int_t number_of_pseudo_experiments_1d, const Int_t number_of_pseudo_experiments_2d)
{
    this->number_of_pseudo_experiments = number_of_pseudo_experiments_1d;
    this->number_of_pseudo_experiments_2d = number_of_pseudo_experiments_2d;
}


////////////////////////////////////////////////////////////////////////////////
// STATIC CONSTANTS
////////////////////////////////////////////////////////////////////////////////

// Q value of decay, MeV
const Double_t Analysis::bb_Q{3.034};
