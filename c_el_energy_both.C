void c_el_energy_both()
{
//=========Macro generated from canvas: e_el_energy_both/e_el_energy_both
//=========  (Thu Apr 26 15:24:52 2018) by ROOT version6.08/06
   TCanvas *e_el_energy_both = new TCanvas("e_el_energy_both", "e_el_energy_both",0,0,800,600);
   e_el_energy_both->SetHighLightColor(2);
   e_el_energy_both->Range(0,0,1,1);
   e_el_energy_both->SetFillColor(0);
   e_el_energy_both->SetBorderMode(0);
   e_el_energy_both->SetBorderSize(2);
   e_el_energy_both->SetFrameBorderMode(0);
   
   TH1D *h_el_energy_original__4 = new TH1D("h_el_energy_original__4","",40,0,4);
   h_el_energy_original__4->SetBinContent(1,10);
   h_el_energy_original__4->SetBinError(1,0.01651513);
   h_el_energy_original__4->SetEntries(997002);
   h_el_energy_original__4->SetStats(0);
   h_el_energy_original__4->SetLineColor(2);
   h_el_energy_original__4->SetMarkerColor(2);
   h_el_energy_original__4->GetXaxis()->SetLabelFont(42);
   h_el_energy_original__4->GetXaxis()->SetLabelSize(0.035);
   h_el_energy_original__4->GetXaxis()->SetTitleSize(0.035);
   h_el_energy_original__4->GetXaxis()->SetTitleFont(42);
   h_el_energy_original__4->GetYaxis()->SetLabelFont(42);
   h_el_energy_original__4->GetYaxis()->SetLabelSize(0.035);
   h_el_energy_original__4->GetYaxis()->SetTitleSize(0.035);
   h_el_energy_original__4->GetYaxis()->SetTitleFont(42);
   h_el_energy_original__4->GetZaxis()->SetLabelFont(42);
   h_el_energy_original__4->GetZaxis()->SetLabelSize(0.035);
   h_el_energy_original__4->GetZaxis()->SetTitleSize(0.035);
   h_el_energy_original__4->GetZaxis()->SetTitleFont(42);
   h_el_energy_original__4->Draw("E");
   
   TH1D *h_el_energy_reweight__5 = new TH1D("h_el_energy_reweight__5","",40,0,4);
   h_el_energy_reweight__5->SetBinContent(1,10);
   h_el_energy_reweight__5->SetBinError(1,0.01623247);
   h_el_energy_reweight__5->SetEntries(997002);
   h_el_energy_reweight__5->SetStats(0);
   h_el_energy_reweight__5->SetLineColor(3);
   h_el_energy_reweight__5->SetMarkerColor(3);
   h_el_energy_reweight__5->GetXaxis()->SetLabelFont(42);
   h_el_energy_reweight__5->GetXaxis()->SetLabelSize(0.035);
   h_el_energy_reweight__5->GetXaxis()->SetTitleSize(0.035);
   h_el_energy_reweight__5->GetXaxis()->SetTitleFont(42);
   h_el_energy_reweight__5->GetYaxis()->SetLabelFont(42);
   h_el_energy_reweight__5->GetYaxis()->SetLabelSize(0.035);
   h_el_energy_reweight__5->GetYaxis()->SetTitleSize(0.035);
   h_el_energy_reweight__5->GetYaxis()->SetTitleFont(42);
   h_el_energy_reweight__5->GetZaxis()->SetLabelFont(42);
   h_el_energy_reweight__5->GetZaxis()->SetLabelSize(0.035);
   h_el_energy_reweight__5->GetZaxis()->SetTitleSize(0.035);
   h_el_energy_reweight__5->GetZaxis()->SetTitleFont(42);
   h_el_energy_reweight__5->Draw("Esame");
   e_el_energy_both->Modified();
   e_el_energy_both->cd();
   e_el_energy_both->SetSelected(e_el_energy_both);
}
