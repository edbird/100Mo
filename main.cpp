





#include <iostream>


//#include "ReWeight.hpp"
//#include "read_data.hpp"
//#include "aux.hpp"

#include "program_arguments.hpp"


#include "Analysis.hpp"

//#include "HistogramWrapper.hpp"
// TODO: make some nice histograms
// X/Y AXIS: LABEL: TEXT, FONT SIZE, FONT
// X/Y AXIS: TICK NUMBERS: FONT SIZE, FONT
// X/Y AXIS LIMITS, RANGE, PARTICULARLY Y MAX
// GRAPH COLORS




int main(int argc, char* argv[])
{
    
    ////////////////////////////////////////////////////////////////////////////
    // PROCESS PROGRAM ARGUMENTS
    ////////////////////////////////////////////////////////////////////////////

    // new program arguments processing method
    ProgramArguments pa;
    pa.Add("help", "--help", "false");
    pa.Add("filename", "--filename", "NewElectronNtuplizerExe_Int_ManDB_output.root");
    pa.Add("epsilon", "--epsilon", "0.0");
    pa.Add("batch_mode", "--batch-mode", "false");
    //pa.Add("300 keV energy cut", "--energy-cut", "false");
    pa.Add("energy_cut", "--energy-cut", "false");
    pa.Add("fit_subrange", "--fit-subrange", "false");
    pa.Add("log_mode", "--log-mode", "true");
    pa.Add("output_filename", "--output-file", "of_data.txt");
    pa.Print();
    pa.Process(argc, argv);

    // this argument is a string
    std::string help{pa.Get("help")};
    // this argument is a string
    std::string filename{pa.Get("filename")};
    // this argument is to be converted to a double
    std::string arg_epsilon_31{pa.Get("epsilon")};
    // this argument is to be converted to a bool
    std::string arg_batch_mode{pa.Get("batch_mode")};
    std::string arg_energy_cut{pa.Get("energy_cut")};
    std::string arg_fit_subrange{pa.Get("fit_subrange")};
    std::string arg_log_mode{pa.Get("log_mode")};
    std::string arg_output_filename{pa.Get("output_filename")};

    double epsilon_31{std::stod(arg_epsilon_31)};
    bool batch_mode{false};
    if(arg_batch_mode == std::string("true"))
    {
        batch_mode = true;
        std::cout << "[ INFO ] : Batch mode: " << "true" << std::endl;
    }
    else
    {
        std::cout << "[ INFO ] : Batch mode: " << "false" << std::endl;
    }

    if(arg_energy_cut == std::string("true"))
    {
        std::cout << "[ INFO ] : Energy cut: " << "true" << std::endl;
    }
    else
    {
        std::cout << "[ INFO ] : Energy cut: " << "false" << std::endl;
    }

    if(arg_fit_subrange == std::string("true"))
    {
        std::cout << "[ INFO ] : Fit subrange: " << "true" << std::endl;
    }
    else
    {
        std::cout << "[ INFO ] : Fit subrange: " << "false" << std::endl;
    }

    if(arg_log_mode == std::string("true"))
    {
        std::cout << "[ INFO ] : Log mode: " << "true" << std::endl;
    }
    else
    {
        std::cout << "[ INFO ] : Log mode: " << "false" << std::endl;
    }

    std::cout << "Writing data to file " << arg_output_filename << " append mode" << std::endl;

    // process gathered argument data
    bool gen_weight_enable{false};
    if(filename != pa.GetDefault("filename"))
    {
        gen_weight_enable = true;
    }
    // TODO: check filename exists!

    //epsilon_31 = std::stod(arg_epsilon_31);
    // todo: check valid
    if(0.0 <= epsilon_31 && epsilon_31 <= 10.0)
    {
        std::cout << "[ INFO ] : Set epsilon_31 = " << epsilon_31 << std::endl;
    }
    else
    {
        std::cout << "invalid epsilon_31 value" << std::endl;
        throw "invalid epsilon_31 value";
    }

    /*
    if(arg_batch_mode == std::string("true"))
    {
        batch_mode = true;
        std::string batch_mode_enable_string("false");
        if(batch_mode)
        {
            batch_mode_enable_string = std::string("true");
        }
        std::cout << "[ INFO ] : Batch mode: " << batch_mode_enable_string << std::endl;
    }
    */


    Analysis analysis(filename, arg_output_filename);
    analysis.SetGenWeightEnabled(gen_weight_enable);
    analysis.SetBatchMode(batch_mode);
    analysis.SetCanvasEnableRawData(false);
    analysis.SetCanvasEnableDecayRate(false);

    if(arg_fit_subrange == std::string("true"))
    {
        analysis.SetFitSubrange(true);
    }

    // enable/disable printing with log canvas
    bool log_mode{false};
    if(arg_log_mode == std::string("true"))
    {
        log_mode = true;
    }
    analysis.SetLogMode(log_mode);
    

    if(arg_energy_cut == std::string("true"))
    {
        analysis.SetEnergyCutEnabled(true);
    }


    analysis.SetEpsilon31(epsilon_31);

    // run analysis
    analysis.ReadData();
    //analysis.CanvasRawData();
    analysis.CanvasDecayRate();
    analysis.CanvasSingleElectronProjection();
    analysis.CanvasSingleElectronTest();
    analysis.InitEventLoop();
    analysis.EventLoop();
    analysis.PostProcess();
    analysis.SummedEnergyFit();
    analysis.SensitivityMeasurementChisquare1();
    analysis.SensitivityMeasurementChisquare2();
    analysis.SensitivityMeasurementLoglikelihood1();
    analysis.SensitivityMeasurementLoglikelihood2();
    analysis.PrintOutputToFile();




    ////////////////////////////////////////////////////////////////////////////
    // ELECTRON ENERGY CONVERSION
    ////////////////////////////////////////////////////////////////////////////

    // This is a test code

    /*
    TH2D *h_convert_from{nullptr};
    TH2D *h_convert_to{nullptr};

    // convert from epsilon_31 = 0.0 to epsilon_31 = 0.8 
    h_convert_from = h_data_0;
    h_convert_to = h_data_2;

    // electron energy in units of Q value
    // the sum must be less than or equal to 1.0
    // (1.0 - E1 - E2 = neutrino energy)
    Double_t electron_energy_1{0.33};
    Double_t electron_energy_2{0.33};

    // get weights
    Double_t weight_in{h_convert_from->GetBinContent(h_convert_from->FindBin(electron_energy_1, electron_energy_2))};
    Double_t weight_out{h_convert_to->GetBinContent(h_convert_to->FindBin(electron_energy_1, electron_energy_2))};
    Double_t weight_factor{weight_out / weight_in};

    std::cout << "Weight factor is " << weight_factor << std::endl;
    */




    //ReWeight(0.5, 0.5, 0.0, data_nEqNull, data_nEqTwo); 
    //ReWeight(0.45, 0.3, 0.8, h_nEqNull, h_nEqTwo, psiN0, psiN2); 


   



    return 0;

}
