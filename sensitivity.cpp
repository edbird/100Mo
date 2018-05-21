

#include "CanvasFactory.hpp"


#include "program_arguments.hpp"
#include "read_data.hpp"


//#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"


#include <vector>
#include <iostream>



Double_t fit_quadratic(Double_t *x, Double_t *par)
{

    Double_t a{par[0]};
    Double_t b{par[1]};
    Double_t c{par[2]};

    return a * (*x) * (*x) + b * (*x) + c;

}


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
    if(arg_input_file == std::string("of_data_testing.txt"))
    {
        std::vector<std::vector<double>> data_out;
        for(std::size_t j{1}; j < data.size() - 1; ++ j)
        {
            std::vector<double> temp;
            temp.push_back(data.at(j).at(0));
            for(std::size_t i{6}; i < data.at(j).size() && i < 100; ++ i)
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

    // epsilon, minimum and width
    std::vector<Double_t> v_minimum;
    std::vector<Double_t> v_width;

    // start from 1, not 0 to avoid header line -------------VVV
    // start from 6, not 0 ------VVV--- to avoid epsilon value and other chi-square values
    for(std::size_t experiment_ix{6}; experiment_ix < data.at(1).size(); ++ experiment_ix)
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

        // find the minimum
        Double_t x0{data_x[0]};
        Double_t y0{data_y[0]};
        Int_t x0_index{0};
        for(Int_t ix{0}; ix < data_size - 1; ++ ix)
        {
            if(data_y[ix] < y0)
            {
                x0 = data_x[ix];
                y0 = data_y[ix];
                x0_index = ix;
            }
        }
        //std::cout << "Experiment " << experiment_ix << std::endl;
        //std::cout << "the minimum is at x=" << x0 << " y=" << y0 << " index is " << x0_index << std::endl;
        v_minimum.push_back(x0);
        // find the 1 sigma width
        Double_t x_high;
        Double_t x_low;
        Double_t x_high_simple;
        Double_t x_low_simple;
        for(Int_t ix{x0_index}; ix < data_size - 1; ++ ix)
        {
            if(data_y[ix] - y0 >= 1.0)
            {
                //std::cout << "found delta >= 1 sigma: delta=" << data_y[ix] - y0 << std::endl;
                x_high_simple = data_x[ix];
                Double_t y_max_1sigma{y0 + 1.0};
                Double_t y_minus{data_y[ix - 1]};
                Double_t y_plus{data_y[ix]};
                Double_t x_minus{data_x[ix - 1]};
                Double_t x_plus{data_x[ix]};
                x_high = x_minus + (x_plus - x_minus) * ((y_max_1sigma - y_minus) / (y_plus - y_minus));
                //std::cout << "interpolate: " << x_minus << " " << x_high << " " << x_plus << std::endl;
                break;
            }
        }
        for(Int_t ix{x0_index}; ix >= 0; -- ix)
        {
            if(data_y[ix] - y0 >= 1.0)
            {
                //std::cout << "found delta >= 1 sigma: delta=" << data_y[ix] - y0 << std::endl;
                x_low_simple = data_x[ix];
                Double_t y_max_1sigma{y0 + 1.0};
                Double_t y_minus{data_y[ix - 1]};
                Double_t y_plus{data_y[ix]};
                Double_t x_minus{data_x[ix - 1]};
                Double_t x_plus{data_x[ix]};
                x_low = x_minus + (x_plus - x_minus) * ((y_max_1sigma - y_minus) / (y_plus - y_minus));
                //std::cout << "interpolate: " << x_minus << " " << x_low << " " << x_plus << std::endl;
                break;
            }
        }

        // TODO interpolate

        //std::cout << "width=" << x_high_simple - x_low_simple << std::endl;
        //std::cout << "width=" << x_high - x_low << std::endl;
        v_width.push_back(x_high - x_low);
        //std::cin.get();

        


        
        // data does not fit a quadratic well
#if 0
        Double_t min{data_x[0]};
        Double_t max{data_x[data_size - 2]};
        TF1 *f_quad = new TF1("f_quad", fit_quadratic, min, max, 3);
        {
            Double_t x0{data_x[0]};
            Double_t y0{data_y[0]};
            // find minimum y
            for(Int_t ix{0}; ix < data_size - 1; ++ ix)
            {
                if(data_y[ix] < y0)
                {
                    x0 = data_x[ix];
                    y0 = data_y[ix];
                }
            }
            // assume last point is max y
            Double_t xmax{data_x[data_size - 2]};
            Double_t ymax{data_y[data_size - 2]};
            Double_t beta{y0};
            Double_t alpha{(ymax - beta) / std::pow(xmax - x0, 2.0)};
            Double_t a{alpha};
            Double_t b{-2.0 * alpha * x0};
            Double_t c{alpha * x0 * x0 + beta};
            f_quad->SetParameter(0, a);
            f_quad->SetParameter(1, b);
            f_quad->SetParameter(2, c);
        }
        TGraph *g_temp = new TGraph(data_size - 1, data_x, data_y);
        //g_temp->Fit("pol2", "M");
        g_temp->Fit("f_quad");
        std::string canvas_name{std::string("sensitivity_") + std::to_string(experiment_ix)};
        TCanvas *c_temp = new TCanvas(canvas_name.c_str(), "", 800, 600);
        //g_temp->Draw("ALP");
        g_temp->SetMarkerStyle(33);
        g_temp->Draw("AP");
        c_temp->SaveAs((canvas_name + std::string(".png")).c_str());
        delete g_temp;
        
        if(experiment_ix == 6)
        {
            std::ofstream ofs_exp_6("ofs_experiment_6.txt");
            for(Int_t i{0}; i < data_size - 1; ++ i)
            {
                ofs_exp_6 << data_x[i] << " " << data_y[i] << std::endl;
            }
            ofs_exp_6.close();
        }
#endif


    }

    TH1D *h_minimum = new TH1D("h_minimum", "", 50, 0.0, 0.8);
    TH1D *h_width = new TH1D("h_width", "", 100, 0.1, 0.15);

    std::ofstream ofs_min_width("ofs_min_width.txt");
    for(std::size_t i{0}; i < v_minimum.size(); ++ i)
    {
        ofs_min_width << v_minimum.at(i) << " " << v_width.at(i) << std::endl;
        h_minimum->Fill(v_minimum.at(i));
        h_width->Fill(v_width.at(i));
    }
    ofs_min_width.close();

    CanvasFactorySettings settings("Minimum", "Experiments", 0.0, 150000.0, false);
    CanvasFactory factory(settings);
    factory.Canvas("h_minimum", ".", h_minimum, "Minimum");
    
    CanvasFactorySettings settings2("Width", "Experiments", 0.0, 60000.0, false);
    CanvasFactory factory2(settings2);
    factory2.Canvas("h_width", ".", h_width, "Width");


    //std::cin.get();

   
    return 0;
}
