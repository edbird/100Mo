





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
    pa.Add("epsilon", "--epsilon", "0.8");
    pa.Add("batch_mode", "--batch-mode", "false");
    //pa.Add("300 keV energy cut", "--energy-cut", "false");
    pa.Add("energy_cut", "--energy-cut", "false");
    pa.Add("fit_subrange", "--fit-subrange", "false");
    pa.Add("log_mode", "--log-mode", "false");
    pa.Add("output_filename", "--output-file", "of_data.txt");
    pa.Add("eps_min", "--epsilon-min", "0.0");
    pa.Add("eps_max", "--epsilon-max", "1.0");
    pa.Add("eps_num", "--epsilon-num", "-1");
    pa.Add("syst_energy_mult", "--systematic-energy-mult", "1.0");
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
    std::string arg_eps_min{pa.Get("eps_min")};
    std::string arg_eps_max{pa.Get("eps_max")};
    std::string arg_eps_num{pa.Get("eps_num")};
    std::string arg_systematic_energy_mult{pa.Get("syst_energy_mult")};

    double systematic_energy_mult{std::stod(arg_systematic_energy_mult)};
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
        std::cout << "gen_weight_enable" << std::endl;
    }
    else
    {
        std::cout << "gen_weight_enable = false" << std::endl;
    }
    // TODO: check filename exists!

    //epsilon_31 = std::stod(arg_epsilon_31);
    // todo: check valid
    /*if(0.0 <= epsilon_31 && epsilon_31 <= 10.0)
    {
        std::cout << "[ INFO ] : Set epsilon_31 = " << epsilon_31 << std::endl;
    }
    else
    {
        std::cout << "invalid epsilon_31 value" << std::endl;
        throw "invalid epsilon_31 value";
    }*/

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

    // set systematic energy multiplier
    // note: could have multiple energy calibration points, hence linear
    // calibration between these points
    //analysis.SetSystematicEnergyMultiplier(systematic_energy_mult);
    analysis.SetSystematicEnergyMultiplierHighLow(1.0 + 0.002, 1.0 - 0.002); // 0.2% energy mult
    analysis.SetSystematicEnergyMultiplierEnabled(false);
    #if 0
        analysis.SetSystematicEnergyOffset(0.003, -0.003); // 0.003 MeV = 3 keV
    #endif
    #if 1
        analysis.SetSystematicEfficiency(1.0 + 0.05, 1.0 - 0.05);
    #endif
    // Notes: to avoid problem described below, use energy mulitplier variable
    // for all systematics
    // problem: number of "histograms" to create (data sets) goes as 3^m where
    // m is the number of systematics
    // Further Notes: change of plan: going to implement something half way
    // between these 2 options
    // need to eventually use code to fix systematic parameters, therefore they
    // should all be included in the code simultaniously
    // however tracking 3^m datasets is not convenient/trivial
    // instead, track only datasets for high/low/default systematic energy
    // multiplier
    // but include all other systmatics, although there will still only be
    // 3 tracked datasets, and they will be indexed using the systematic energy
    // multiplier
    // NEW PROBLEM: when systematic energy multiplier is the same, we can't
    // track all 3 datasets!
    // SOLUTION: introduce a flag to enable/disable the systematic energy
    // multiplier

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


    // THIS CODE FOR PRODUCING THE TEST HISTOGRAMS

#if 0
    
    //analysis.SetEpsilon31(epsilon_31);
    analysis.SetEpsilon31(0.5);

    // run analysis
    analysis.ReadData();
    //analysis.CanvasRawData();
    analysis.CanvasDecayRate();
    analysis.CanvasSingleElectronProjection();
    analysis.InitSingleElectronTest();
    analysis.InitEventLoopTree();

    // event loop
    analysis.InitEventLoopHistogram();
    analysis.EventLoop();
    analysis.PostProcess();

    analysis.SummedEnergyFit();
    analysis.SensitivityMeasurementChisquare1();
    analysis.SensitivityMeasurementChisquare2();
    analysis.SensitivityMeasurementLoglikelihood1();
    analysis.SensitivityMeasurementLoglikelihood2();
    analysis.PrintOutputToFile();

    // note: assumed this can go here?
    analysis.CanvasSingleElectronTest();
    

#endif


    // THIS CODE FOR REGULAR DATA ANALYSIS RUN

#if 1

    // run analysis
    analysis.ReadData();
    //analysis.CanvasRawData();
    analysis.CanvasDecayRate();
    analysis.CanvasSingleElectronProjection();
    ///analysis.InitSingleElectronTest(); // have not integrated this into the multi run mode
    analysis.InitEventLoopTree();
//<<<<<<< HEAD
//    
//    // set number of pseudoexperiments
//    analysis.SetNumberOfPseudoexperiments(10000, 1);
//=======
    //analysis.InitEventLoopHistogram();
    analysis.SetNumberOfPseudoexperiments(100, 100);
//>>>>>>> 39030f95ed5a3920b8a999e9180245e89cd9b007

    //Double_t eps_incr{0.025};
    Double_t eps_min{std::stod(arg_eps_min)};
    Double_t eps_max{std::stod(arg_eps_max)};
    Int_t eps_num{std::stoi(arg_eps_num)};
    if(eps_num < 1)
    {
//<<<<<<< HEAD
        
        if(pa.WasProvided("epsilon") == true)
        {
            std::string arg_eps_value{pa.Get("epsilon")};
            Double_t eps_value{std::stod(arg_eps_value)};
            analysis.AddEpsilonValue(eps_value);
        }
        else
        {
            std::cout << "error" << std::endl;
        }

    }
    else if(eps_num == 1)
    {
        // TODO: should switch mode here depending on if eps_min, eps_max are
        // provided or a single eps is provided as prog arg
        std::cout << "Running: eps=" << eps_min << std::endl;
        analysis.AddEpsilonValue(eps_min);
    }
    else
    {
        Double_t eps_incr{(eps_max - eps_min) / (Double_t)(eps_num - 1)};
        for(Double_t eps{eps_min}; eps <= eps_max + 0.5 * eps_incr; eps += eps_incr)
        {
            std::cout << "Running: eps=" << eps << std::endl;
            analysis.AddEpsilonValue(eps);
        }
    }
//=======
//        std::cout << "Running: eps=" << eps << std::endl;
//        analysis.AddEpsilonValue(eps);
//    }
//
//>>>>>>> 39030f95ed5a3920b8a999e9180245e89cd9b007

    //analysis.AddEpsilonValue(0.368);
    //analysis.AddEpsilonValue(0.5);

    analysis.RunOverEpsilonVector();
    
//<<<<<<< HEAD
//    
//=======
    /* 
    analysis.SummedEnergyFit();
    analysis.SensitivityMeasurementChisquare1();
    analysis.SensitivityMeasurementChisquare2();
    analysis.SensitivityMeasurementLoglikelihood1();
    analysis.SensitivityMeasurementLoglikelihood2();
    analysis.PrintOutputToFile();
    */
//>>>>>>> 39030f95ed5a3920b8a999e9180245e89cd9b007
#endif


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
