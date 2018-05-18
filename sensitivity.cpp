

#include "program_arguments.hpp"
#include "read_data.hpp"


//#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"


#include <vector>
#include <iostream>



int main(int argc, char* argv[])
{

    ////////////////////////////////////////////////////////////////////////////
    // PROCESS PROGRAM ARGUMENTS
    ////////////////////////////////////////////////////////////////////////////

    ProgramArguments pa;
    pa.Add("help", "--help", "false");
    pa.Add("input_file", "--input-file", "of_data_testing.txt");
    pa.Print();
    pa.Process(argc, argv);

    std::string arg_input_file{pa.Get("input_file")};

    std::cout << "input_file=" << arg_input_file << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    // READ DATA IN FROM FILES
    ////////////////////////////////////////////////////////////////////////////

    std::vector<std::vector<double>> data;
    read_data_helper_2(arg_input_file.c_str(), data, ',');

    //std::cout << data.size() << std::endl;
    //std::cout << data.at(0).size() << std::endl;
    //std::cout << data.at(1).size() << std::endl;

    /*
    for(std::size_t j{0}; j < data.size(); ++ j)
    {
        for(std::size_t i{0}; i < data.at(j).size(); ++ i)
        {
            std::cout << data.at(j).at(i) << ", ";
        }
        std::cout << std::endl;
    }
    */

    //std::cout << data.at(0).at(0) << " " << data.at(0).at(1) << std::endl;
    // the first line should be ignored
    
    // print output file with subset of data
    if(arg_input_file == std::string("of_data_texting.txt"))
    {
        std::vector<std::vector<double>> data_out;
        for(std::size_t j{1}; j < data.size() - 1; ++ j)
        {
            std::vector<double> temp;
            for(std::size_t i{0}; i < data.at(j).size() && i < 100; ++ i)
            {
                temp.push_back(data.at(j).at(i));
            }
            data_out.push_back(temp);
        }

        write_data_helper_2("of_data_testing_out.txt", data_out, ',');
    }


    Int_t data_size{data.size() - 2};
    std::cout << "data_size=" << data_size << std::endl;
    Double_t *data_x{new Double_t[data_size - 1]};
    Double_t *data_y{new Double_t[data_size - 1]}; // I am now confused as fuck

    // start from 1, not 0 to avoid header line -------------VVV
    // start from 1, not 0 ------VVV--- to avoid epsilon value
    for(std::size_t experiment_ix{1}; experiment_ix < data.at(1).size(); ++ experiment_ix)
    {
        //std::cout << "experiment_ix=" << experiment_ix << std::endl;
        // "epsilon" ix
        // start from 1 ------VVV--- not 0 to avoid header line
        for(std::size_t eps_ix{1}; eps_ix < data_size; ++ eps_ix)
        {
            //std::cout << "eps_ix=" << eps_ix; 
            data_x[eps_ix - 1] = data.at(eps_ix).at(0);
            data_y[eps_ix - 1] = data.at(eps_ix).at(experiment_ix);
            //std::cout << ", " << data_x[eps_ix - 1] << "," << data_y[eps_ix - 1];
            //std::cout << std::endl;

            //if(eps_ix >= 100) break;
        }

        TGraph *g_temp = new TGraph(data_size - 1, data_x, data_y);
        std::string canvas_name{std::string("sensitivity_") + std::to_string(experiment_ix)};
        TCanvas *c_temp = new TCanvas(canvas_name.c_str(), "", 800, 600);
        g_temp->Draw("ALP");
        c_temp->SaveAs((canvas_name + std::string(".png")).c_str());
        delete g_temp;

    }


    //std::cin.get();

   
    return 0;
}
