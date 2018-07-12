#include "Analysis.hpp"


////////////////////////////////////////////////////////////////////////////////
// 1D SINGLE ELECTRON ENERGY LOG LIKELIHOOD SENSITIVITY
////////////////////////////////////////////////////////////////////////////////

void Analysis::SensitivityMeasurementLoglikelihood1()
{

    const std::string eps_string{std::to_string(epsilon_31)};

    #if 0
    std::string c_name_ll{std::string("c_ll/") + std::string("c_ll_") + eps_string};
    c_ll = new TCanvas(c_name_ll.c_str(), "", 800, 600);
    h_ll->Draw("E");
    c_ll->SaveAs((c_name_ll + std::string(".png")).c_str());
    delete c_ll;
    #endif
        
}
