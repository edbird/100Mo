#include "ratio_fit_function.hpp"

// fit function for ratio histograms to calculate chisquare
Double_t ratio_chisquare_function(Double_t *x)
{

    Int_t number_of_bins{40};

    Double_t x_val{new Double_t[number_of_bins + 1]};
    Double_t y_val{new Double_t[number_of_bins + 1]};

    for(Int_t i{0}; i <= number_of_bins; ++ i)
    {
        Int_t ix{2 * i + 0};
        Int_t jx{2 * i + 1};

        x_val = par[ix];
        y_val = par[jx];
    }

    // systematic: energy scale parameter
    Double_t m_exp{par[2 * number_of_bins + 2]};
    Double_t m_sig{par[2 * number_of_bins + 3]};
    std::cout << "m_exp=" << m_exp << " m_sig=" << m_sig << std::endl;

    Double_t xx{x[0]};





    Double_t chi{0.0};
    for()
    {
        Double_t d_i{};
        Double_t m_i{};
        Double_t s_i{};

        Double_t chi_1{(d_i - m_i) / s_i};
        Double_t chi_2{(m - m_exp) / m_sig};
        
        chi_1 = chi_1 * chi_1;
        chi_2 = chi_2 * chi_2;
        
        chi += chi_1 + chi_2;
    }


}
