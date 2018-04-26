#ifndef REWEIGHT_HPP
#define REWEIGHT_HPP

#include <cmath>

// electron energies: T1, T2 (in units of Q value)
// theory parameter: epsilon
// data: G0, G2 (units of Q value)
//Double_t ReWeight(const Double_t T1, const Double_t T2, const Double_t epsilon,
//                  const std::vector<std::vector<double>>& data_nEqNull,
//                  const std::vector<std::vector<double>>& data_nEqTwo)
Double_t ReWeight(const Double_t T1, const Double_t T2, const Double_t epsilon,
                  const TH2D* const h_nEqNull,
                  const TH2D* const h_nEqTwo,
                  const Double_t psiN0, const Double_t psiN2,
                  const std::string& debug)
{

    // TODO: NO INTERPOLATION DONE YET

    /*
    std::size_t l{data_nEqNull.size()};
    Int_t sub_length{1001};

    // "primary" index
    Int_t index1{(Int_t)(Double_t)(sub_length) * T1}; index1 *= sub_length;
    Int_t index3{(Int_t)(Double_t)(sub_length) * T2};
    //Int_t index1{(Int_t)((Double_t)(sub_length * sub_length) * T1 + (Double_t)(sub_length) * T2)};
    Int_t index2{index1 + 1};
    //Int_t index3{index1 + sub_length};
    //Int_t index4{index1 + sub_length + 1};
    Int_t index4{index3 + 1};

    Double_t energy1{data_nEqNull.at(index1).at(0)};
    Double_t energy2{data_nEqNull.at(index2).at(1)};
    Double_t energy3{data_nEqNull.at(index3).at(0)};
    Double_t energy4{data_nEqNull.at(index4).at(1)};

    Double_t probability1{data_nEqNull.at(index1).at(2)};
    Double_t probability2{data_nEqNull.at(index2).at(2)};
    Double_t probability3{data_nEqNull.at(index3).at(2)};
    Double_t probability4{data_nEqNull.at(index4).at(2)};
    */

    // interpolate in 2D

    // convert from input energy T1 to "higher and lower energies" as used in
    // histogram data_nEqNull / data_nEqTwo
    
    //std::cout << "T1=" << T1 << std::endl;
    //std::cout << "T2=" << T2 << std::endl;
    Int_t bin_x{h_nEqNull->GetXaxis()->FindBin(T1)};
    Int_t bin_y{h_nEqNull->GetYaxis()->FindBin(T2)};
    //std::cout << "bin_x=" << bin_x << std::endl;
    //std::cout << "bin_y=" << bin_y << std::endl;
    //Double_t bin_x_low{h_nEqNull->GetXaxis()->GetBinLowEdge(bin_x)};
    //Double_t bin_x_high{h_nEqNull->GetXaxis()->GetBinLowEdge(bin_x) + h_nEqNull->GetXaxis()->GetBinWidth(bin_x)};
    //Double_t bin_y_low{h_nEqNull->GetXaxis()->GetBinLowEdge(bin_y)};
    //Double_t bin_y_high{h_nEqNull->GetXaxis()->GetBinLowEdge(bin_y) + h_nEqNull->GetXaxis()->GetBinWidth(bin_y)};
    //std::cout << bin_x_low << " < " << T1 << " < " << bin_x_high << std::endl;
    //std::cout << bin_y_low << " < " << T2 << " < " << bin_y_high << std::endl;

    // get the weight for this T1, T2
    // the input data is for epsilon = 0.0
    Double_t phase_1{1.0 / psiN0};
    Double_t weight_1{phase_1 * h_nEqNull->GetBinContent(bin_x, bin_y)};
    //std::cout << "weight_1=" << weight_1 << std::endl;

    // get the weight for this T1, T2
    // the input data is for epsilon = some arbitary value
    Double_t phase_2{1.0 / (psiN0 + epsilon * psiN2)};
    Double_t weight_2{phase_2 * (h_nEqNull->GetBinContent(bin_x, bin_y) + epsilon * h_nEqTwo->GetBinContent(bin_x, bin_y))};
    //std::cout << "weight_2=" << weight_2 << std::endl;

    if(debug == "true")
    {
        if(std::isnan(weight_2 / weight_1))
        {
            std::cout << "T1=" << T1 << " T2=" << T2 << std::endl;
            std::cout << "bin_x=" << bin_x << " bin_y=" << bin_y << std::endl;
            std::cout << "h_nEqNull->GetBinContent(bin_x, bin_y)=" << h_nEqNull->GetBinContent(bin_x, bin_y) << std::endl;
            std::cout << "h_nEqTwo->GetBinContent(bin_x, bin_y)=" << h_nEqTwo->GetBinContent(bin_x, bin_y) << std::endl;
            std::cout << "weight_1=" << weight_1 << " weight_2=" << weight_2 << std::endl;
        }
    }

    //std::cout << "reweight factor: " << weight_2 / weight_1 << std::endl;
    return weight_2 / weight_1;

    //std::cout << "T1=" << T1 << " T2=" << T2 << " energy1=" << energy1 << " energy2=" << energy2 << " energy3=" << energy3 << " energy4=" << energy4 << std::endl;
    //std::cout << "T1=" << T1 << " probability1=" << probability1 << " probability2=" << probability2 << " probability3=" << probability3 << " probability4=" << probability4 << std::endl;

}


Double_t chi_square_test(const TH1D* const histo_base, const TH1D* const histo_test)
{
    if(histo_base->GetNbinsX() != histo_test->GetNbinsX())
    {
        throw "invalid bin number error";
    }

    Double_t sum{0.0};
    for(Int_t ix{1}; ix <= histo_base->GetNbinsX(); ++ ix)
    {
        Double_t content_1{histo_base->GetBinContent(ix)};
        Double_t content_2{histo_test->GetBinContent(ix)};
        Double_t delta{content_1 - content_2};
        Double_t error{histo_test->GetBinError(ix)};
        if(error <= 0.0)
        {
            if(content_2 == 0.0)
            {
                continue;
            }
            else
            {
                std::cerr << "Warning: Skipping bin with index " << ix << " in chi_square_test, bin error is 0.0 but content != 0.0" << std::endl;
            }
        }
        Double_t chi{std::pow(delta / error, 2.0)};
        sum += chi;
    }
    sum /= (Double_t)histo_base->GetNbinsX();
    return sum;
}


#endif
