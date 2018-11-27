

int main()
{

    ////////////////////////////////////////////////////////////////////////////
    // input
    ////////////////////////////////////////////////////////////////////////////


    // uncertainty (u=statistical, s=systematic)
    // for epsilon paramter
    double u_stat{0.022}; // statistical
    double s_mult{0.019}; // multiplication
    double s_add{0.133}; // addition
    double s_eff{0.0}; // efficiency
    double s_nlcorr{0.115} // nonlinear optical correction

    std::vector<double> s_vec;
    s_vec.push_back(s_mult);
    s_vec.push_back(s_add);
    s_vec.push_back(s_eff);
    s_vec.push_back(s_nlcorr);



    ////////////////////////////////////////////////////////////////////////////
    // main prog
    ////////////////////////////////////////////////////////////////////////////

    double s_all{0.0};
    for(std::vector<double>::iterator it{s_vec.begin()};
        it != s_vec.end(); ++ it)
    {
        s_all += (*it) * (*it);
    }
    s_all = std::sqrt(s_all);


    return 0;
}
