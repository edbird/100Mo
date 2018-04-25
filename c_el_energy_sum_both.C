void c_el_energy_sum_both()
{
//=========Macro generated from canvas: e_el_energy_sum_both/e_el_energy_sum_both
//=========  (Wed Apr 25 11:11:09 2018) by ROOT version6.08/06
   TCanvas *e_el_energy_sum_both = new TCanvas("e_el_energy_sum_both", "e_el_energy_sum_both",0,0,800,600);
   e_el_energy_sum_both->SetHighLightColor(2);
   e_el_energy_sum_both->Range(0,0,1,1);
   e_el_energy_sum_both->SetFillColor(0);
   e_el_energy_sum_both->SetBorderMode(0);
   e_el_energy_sum_both->SetBorderSize(2);
   e_el_energy_sum_both->SetFrameBorderMode(0);
   
   TH1D *h_el_energy_sum_original__6 = new TH1D("h_el_energy_sum_original__6","",100,0,4);
   h_el_energy_sum_original__6->SetBinContent(1,1002623);
   h_el_energy_sum_original__6->SetBinError(1,2341.717);
   h_el_energy_sum_original__6->SetEntries(498501);
   h_el_energy_sum_original__6->SetStats(0);
   h_el_energy_sum_original__6->SetLineColor(2);
   h_el_energy_sum_original__6->SetMarkerColor(2);
   h_el_energy_sum_original__6->GetXaxis()->SetLabelFont(42);
   h_el_energy_sum_original__6->GetXaxis()->SetLabelSize(0.035);
   h_el_energy_sum_original__6->GetXaxis()->SetTitleSize(0.035);
   h_el_energy_sum_original__6->GetXaxis()->SetTitleFont(42);
   h_el_energy_sum_original__6->GetYaxis()->SetLabelFont(42);
   h_el_energy_sum_original__6->GetYaxis()->SetLabelSize(0.035);
   h_el_energy_sum_original__6->GetYaxis()->SetTitleSize(0.035);
   h_el_energy_sum_original__6->GetYaxis()->SetTitleFont(42);
   h_el_energy_sum_original__6->GetZaxis()->SetLabelFont(42);
   h_el_energy_sum_original__6->GetZaxis()->SetLabelSize(0.035);
   h_el_energy_sum_original__6->GetZaxis()->SetTitleSize(0.035);
   h_el_energy_sum_original__6->GetZaxis()->SetTitleFont(42);
   h_el_energy_sum_original__6->Draw("E");
   
   TH1D *h_el_energy_sum_reweight__7 = new TH1D("h_el_energy_sum_reweight__7","",100,0,4);
   h_el_energy_sum_reweight__7->SetBinContent(1,1002156);
   h_el_energy_sum_reweight__7->SetBinError(1,2300.566);
   h_el_energy_sum_reweight__7->SetEntries(498501);
   h_el_energy_sum_reweight__7->SetStats(0);
   h_el_energy_sum_reweight__7->SetLineColor(3);
   h_el_energy_sum_reweight__7->SetMarkerColor(3);
   h_el_energy_sum_reweight__7->GetXaxis()->SetLabelFont(42);
   h_el_energy_sum_reweight__7->GetXaxis()->SetLabelSize(0.035);
   h_el_energy_sum_reweight__7->GetXaxis()->SetTitleSize(0.035);
   h_el_energy_sum_reweight__7->GetXaxis()->SetTitleFont(42);
   h_el_energy_sum_reweight__7->GetYaxis()->SetLabelFont(42);
   h_el_energy_sum_reweight__7->GetYaxis()->SetLabelSize(0.035);
   h_el_energy_sum_reweight__7->GetYaxis()->SetTitleSize(0.035);
   h_el_energy_sum_reweight__7->GetYaxis()->SetTitleFont(42);
   h_el_energy_sum_reweight__7->GetZaxis()->SetLabelFont(42);
   h_el_energy_sum_reweight__7->GetZaxis()->SetLabelSize(0.035);
   h_el_energy_sum_reweight__7->GetZaxis()->SetTitleSize(0.035);
   h_el_energy_sum_reweight__7->GetZaxis()->SetTitleFont(42);
   h_el_energy_sum_reweight__7->Draw("Esame");
   e_el_energy_sum_both->Modified();
   e_el_energy_sum_both->cd();
   e_el_energy_sum_both->SetSelected(e_el_energy_sum_both);
}
