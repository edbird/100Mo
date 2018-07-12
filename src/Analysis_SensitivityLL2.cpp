#include "Analysis.hpp"


////////////////////////////////////////////////////////////////////////////////
// 2D SINGLE ELECTRON ENERGY LOG LIKELIHOOD SENSITIVITY
////////////////////////////////////////////////////////////////////////////////

void Analysis::SensitivityMeasurementLoglikelihood2()
{
    
    const std::string eps_string{std::to_string(epsilon_31)};
    
    #if 0    
    std::string c_name_ll_2d{std::string("c_ll_2d/") + std::string("c_ll_2d_") + eps_string};
    c_ll_2d = new TCanvas(c_name_ll_2d.c_str(), "", 800, 600);
    h_ll_2d->Draw("E");
    c_ll_2d->SaveAs((c_name_ll_2d + std::string(".png")).c_str());
    delete c_ll_2d;
    #endif
    
}
