void c_el_energy_origional()
{
//=========Macro generated from canvas: c_el_energy_origional/c_el_energy_origional
//=========  (Mon Apr 23 16:39:25 2018) by ROOT version6.08/06
   TCanvas *c_el_energy_origional = new TCanvas("c_el_energy_origional", "c_el_energy_origional",0,0,800,600);
   c_el_energy_origional->SetHighLightColor(2);
   c_el_energy_origional->Range(0,0,1,1);
   c_el_energy_origional->SetFillColor(0);
   c_el_energy_origional->SetBorderMode(0);
   c_el_energy_origional->SetBorderSize(2);
   c_el_energy_origional->SetFrameBorderMode(0);
   
   TH1D *h_el_energy_origional__12 = new TH1D("h_el_energy_origional__12","h_el_energy_origional",100,0,4);
   h_el_energy_origional__12->SetBinContent(1,612);
   h_el_energy_origional__12->SetBinContent(2,11565);
   h_el_energy_origional__12->SetBinContent(3,23830);
   h_el_energy_origional__12->SetBinContent(4,37014);
   h_el_energy_origional__12->SetBinContent(5,50513);
   h_el_energy_origional__12->SetBinContent(6,60261);
   h_el_energy_origional__12->SetBinContent(7,68780);
   h_el_energy_origional__12->SetBinContent(8,75109);
   h_el_energy_origional__12->SetBinContent(9,79627);
   h_el_energy_origional__12->SetBinContent(10,81410);
   h_el_energy_origional__12->SetBinContent(11,82573);
   h_el_energy_origional__12->SetBinContent(12,81943);
   h_el_energy_origional__12->SetBinContent(13,80054);
   h_el_energy_origional__12->SetBinContent(14,77544);
   h_el_energy_origional__12->SetBinContent(15,74247);
   h_el_energy_origional__12->SetBinContent(16,70513);
   h_el_energy_origional__12->SetBinContent(17,66078);
   h_el_energy_origional__12->SetBinContent(18,62133);
   h_el_energy_origional__12->SetBinContent(19,58135);
   h_el_energy_origional__12->SetBinContent(20,53224);
   h_el_energy_origional__12->SetBinContent(21,49061);
   h_el_energy_origional__12->SetBinContent(22,44743);
   h_el_energy_origional__12->SetBinContent(23,40617);
   h_el_energy_origional__12->SetBinContent(24,36757);
   h_el_energy_origional__12->SetBinContent(25,32410);
   h_el_energy_origional__12->SetBinContent(26,29231);
   h_el_energy_origional__12->SetBinContent(27,25609);
   h_el_energy_origional__12->SetBinContent(28,22987);
   h_el_energy_origional__12->SetBinContent(29,20000);
   h_el_energy_origional__12->SetBinContent(30,17359);
   h_el_energy_origional__12->SetBinContent(31,15287);
   h_el_energy_origional__12->SetBinContent(32,12978);
   h_el_energy_origional__12->SetBinContent(33,11118);
   h_el_energy_origional__12->SetBinContent(34,9408);
   h_el_energy_origional__12->SetBinContent(35,7974);
   h_el_energy_origional__12->SetBinContent(36,6724);
   h_el_energy_origional__12->SetBinContent(37,5620);
   h_el_energy_origional__12->SetBinContent(38,4505);
   h_el_energy_origional__12->SetBinContent(39,3867);
   h_el_energy_origional__12->SetBinContent(40,3133);
   h_el_energy_origional__12->SetBinContent(41,2515);
   h_el_energy_origional__12->SetBinContent(42,2012);
   h_el_energy_origional__12->SetBinContent(43,1632);
   h_el_energy_origional__12->SetBinContent(44,1232);
   h_el_energy_origional__12->SetBinContent(45,962);
   h_el_energy_origional__12->SetBinContent(46,743);
   h_el_energy_origional__12->SetBinContent(47,568);
   h_el_energy_origional__12->SetBinContent(48,454);
   h_el_energy_origional__12->SetBinContent(49,318);
   h_el_energy_origional__12->SetBinContent(50,230);
   h_el_energy_origional__12->SetBinContent(51,180);
   h_el_energy_origional__12->SetBinContent(52,120);
   h_el_energy_origional__12->SetBinContent(53,89);
   h_el_energy_origional__12->SetBinContent(54,54);
   h_el_energy_origional__12->SetBinContent(55,39);
   h_el_energy_origional__12->SetBinContent(56,25);
   h_el_energy_origional__12->SetBinContent(57,13);
   h_el_energy_origional__12->SetBinContent(58,16);
   h_el_energy_origional__12->SetBinContent(59,6);
   h_el_energy_origional__12->SetBinContent(61,2);
   h_el_energy_origional__12->SetBinContent(68,1);
   h_el_energy_origional__12->SetEntries(1605764);
   h_el_energy_origional__12->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   h_el_energy_origional__12->SetLineColor(ci);
   h_el_energy_origional__12->GetXaxis()->SetLabelFont(42);
   h_el_energy_origional__12->GetXaxis()->SetLabelSize(0.035);
   h_el_energy_origional__12->GetXaxis()->SetTitleSize(0.035);
   h_el_energy_origional__12->GetXaxis()->SetTitleFont(42);
   h_el_energy_origional__12->GetYaxis()->SetLabelFont(42);
   h_el_energy_origional__12->GetYaxis()->SetLabelSize(0.035);
   h_el_energy_origional__12->GetYaxis()->SetTitleSize(0.035);
   h_el_energy_origional__12->GetYaxis()->SetTitleFont(42);
   h_el_energy_origional__12->GetZaxis()->SetLabelFont(42);
   h_el_energy_origional__12->GetZaxis()->SetLabelSize(0.035);
   h_el_energy_origional__12->GetZaxis()->SetTitleSize(0.035);
   h_el_energy_origional__12->GetZaxis()->SetTitleFont(42);
   h_el_energy_origional__12->Draw("E");
   c_el_energy_origional->Modified();
   c_el_energy_origional->cd();
   c_el_energy_origional->SetSelected(c_el_energy_origional);
}
