// self header
#include "Analysis.hpp"


// local headers
#include "read_data.hpp"


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