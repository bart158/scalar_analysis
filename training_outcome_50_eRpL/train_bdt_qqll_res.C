#ifdef __CLING__
#pragma cling optimize(0)
#endif
void train_bdt_qqll_res()
{
//=========Macro generated from canvas: canvas1/TMVA comparison BDT
//=========  (Wed Aug 14 11:36:29 2024) by ROOT version 6.32.02
   TCanvas *canvas1 = new TCanvas("canvas1", "TMVA comparison BDT",200,81,600,468);
   gStyle->SetOptStat(0);
   canvas1->Range(-0.5480649,-1.323864,0.3697165,9.708336);
   canvas1->SetFillColor(0);
   canvas1->SetBorderMode(0);
   canvas1->SetBorderSize(2);
   canvas1->SetTickx(1);
   canvas1->SetTicky(1);
   canvas1->SetRightMargin(0.05);
   canvas1->SetBottomMargin(0.12);
   canvas1->SetFrameBorderMode(0);
   canvas1->SetFrameBorderMode(0);
   
   TH2F *frameBDT__37 = new TH2F("frameBDT__37","TMVA overtraining check for classifier: BDT",500,-0.4516979,0.3238274,500,0,8.605116);
   frameBDT__37->SetStats(0);
   frameBDT__37->GetXaxis()->SetTitle("BDT response");
   frameBDT__37->GetXaxis()->SetLabelOffset(0.012);
   frameBDT__37->GetXaxis()->SetLabelSize(0.04);
   frameBDT__37->GetXaxis()->SetTitleSize(0.045);
   frameBDT__37->GetXaxis()->SetTitleOffset(1.25);
   frameBDT__37->GetYaxis()->SetTitle("(1/N) dN^{ }/^{ }dx");
   frameBDT__37->GetYaxis()->SetLabelOffset(0.012);
   frameBDT__37->GetYaxis()->SetLabelSize(0.04);
   frameBDT__37->GetYaxis()->SetTitleSize(0.045);
   frameBDT__37->GetYaxis()->SetTitleOffset(1.2);
   frameBDT__37->GetZaxis()->SetLabelSize(0.04);
   frameBDT__37->GetZaxis()->SetTitleSize(0.04);
   frameBDT__37->GetZaxis()->SetTitleOffset(1);
   frameBDT__37->Draw("");
   
   TLegend *leg = new TLegend(0.105,0.78,0.505,0.9,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1);
   TLegendEntry *entry=leg->AddEntry("MVA_BDT_S","Signal (test sample)","F");
   entry->SetFillColor(38);
   entry->SetFillStyle(1001);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ee");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("MVA_BDT_B","Background (test sample)","F");
   entry->SetFillColor(2);
   entry->SetFillStyle(3554);
   entry->SetLineColor(2);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   leg->Draw();
   
   TH1D *MVA_BDT_S__38 = new TH1D("MVA_BDT_S__38","TMVA overtraining check for classifier: BDT",40,-0.4516979,0.3238274);
   MVA_BDT_S__38->SetBinContent(13,0.0002343031);
   MVA_BDT_S__38->SetBinContent(14,0.002976872);
   MVA_BDT_S__38->SetBinContent(15,0.002976872);
   MVA_BDT_S__38->SetBinContent(16,0.01324424);
   MVA_BDT_S__38->SetBinContent(17,0.008462011);
   MVA_BDT_S__38->SetBinContent(18,0.03863027);
   MVA_BDT_S__38->SetBinContent(19,0.05006916);
   MVA_BDT_S__38->SetBinContent(20,0.09826389);
   MVA_BDT_S__38->SetBinContent(21,0.1620459);
   MVA_BDT_S__38->SetBinContent(22,0.191346);
   MVA_BDT_S__38->SetBinContent(23,0.3522204);
   MVA_BDT_S__38->SetBinContent(24,0.565465);
   MVA_BDT_S__38->SetBinContent(25,0.6841403);
   MVA_BDT_S__38->SetBinContent(26,1.008466);
   MVA_BDT_S__38->SetBinContent(27,1.484102);
   MVA_BDT_S__38->SetBinContent(28,2.05894);
   MVA_BDT_S__38->SetBinContent(29,2.744182);
   MVA_BDT_S__38->SetBinContent(30,3.809014);
   MVA_BDT_S__38->SetBinContent(31,4.762203);
   MVA_BDT_S__38->SetBinContent(32,5.809477);
   MVA_BDT_S__38->SetBinContent(33,6.177894);
   MVA_BDT_S__38->SetBinContent(34,6.591427);
   MVA_BDT_S__38->SetBinContent(35,6.216221);
   MVA_BDT_S__38->SetBinContent(36,4.639298);
   MVA_BDT_S__38->SetBinContent(37,2.764026);
   MVA_BDT_S__38->SetBinContent(38,1.021173);
   MVA_BDT_S__38->SetBinContent(39,0.2812441);
   MVA_BDT_S__38->SetBinContent(40,0.04020133);
   MVA_BDT_S__38->SetBinError(13,0.0002343031);
   MVA_BDT_S__38->SetBinError(14,0.002530058);
   MVA_BDT_S__38->SetBinError(15,0.002530058);
   MVA_BDT_S__38->SetBinError(16,0.005623317);
   MVA_BDT_S__38->SetBinError(17,0.004369644);
   MVA_BDT_S__38->SetBinError(18,0.009428842);
   MVA_BDT_S__38->SetBinError(19,0.0106957);
   MVA_BDT_S__38->SetBinError(20,0.01511148);
   MVA_BDT_S__38->SetBinError(21,0.01935165);
   MVA_BDT_S__38->SetBinError(22,0.02093773);
   MVA_BDT_S__38->SetBinError(23,0.02850615);
   MVA_BDT_S__38->SetBinError(24,0.03639457);
   MVA_BDT_S__38->SetBinError(25,0.03968394);
   MVA_BDT_S__38->SetBinError(26,0.04820625);
   MVA_BDT_S__38->SetBinError(27,0.05849809);
   MVA_BDT_S__38->SetBinError(28,0.06891144);
   MVA_BDT_S__38->SetBinError(29,0.07955856);
   MVA_BDT_S__38->SetBinError(30,0.09406983);
   MVA_BDT_S__38->SetBinError(31,0.1049911);
   MVA_BDT_S__38->SetBinError(32,0.1161837);
   MVA_BDT_S__38->SetBinError(33,0.1193114);
   MVA_BDT_S__38->SetBinError(34,0.1235387);
   MVA_BDT_S__38->SetBinError(35,0.1197069);
   MVA_BDT_S__38->SetBinError(36,0.1035126);
   MVA_BDT_S__38->SetBinError(37,0.08023017);
   MVA_BDT_S__38->SetBinError(38,0.04859056);
   MVA_BDT_S__38->SetBinError(39,0.02567259);
   MVA_BDT_S__38->SetBinError(40,0.009745505);
   MVA_BDT_S__38->SetEntries(37685);
   MVA_BDT_S__38->SetFillColor(38);

   ci = TColor::GetColor("#0000ee");
   MVA_BDT_S__38->SetLineColor(ci);
   MVA_BDT_S__38->GetXaxis()->SetLabelFont(42);
   MVA_BDT_S__38->GetXaxis()->SetTitleOffset(1);
   MVA_BDT_S__38->GetXaxis()->SetTitleFont(42);
   MVA_BDT_S__38->GetYaxis()->SetLabelFont(42);
   MVA_BDT_S__38->GetYaxis()->SetTitleFont(42);
   MVA_BDT_S__38->GetZaxis()->SetLabelFont(42);
   MVA_BDT_S__38->GetZaxis()->SetTitleOffset(1);
   MVA_BDT_S__38->GetZaxis()->SetTitleFont(42);
   MVA_BDT_S__38->Draw("samehist");
   
   TH1D *MVA_BDT_B__39 = new TH1D("MVA_BDT_B__39","MVA_BDT_B",40,-0.4516979,0.3238274);
   MVA_BDT_B__39->SetBinContent(1,0.04884955);
   MVA_BDT_B__39->SetBinContent(2,0.1314867);
   MVA_BDT_B__39->SetBinContent(3,0.201091);
   MVA_BDT_B__39->SetBinContent(4,0.5966221);
   MVA_BDT_B__39->SetBinContent(5,1.0766);
   MVA_BDT_B__39->SetBinContent(6,1.310159);
   MVA_BDT_B__39->SetBinContent(7,2.013202);
   MVA_BDT_B__39->SetBinContent(8,2.561011);
   MVA_BDT_B__39->SetBinContent(9,3.113303);
   MVA_BDT_B__39->SetBinContent(10,3.434253);
   MVA_BDT_B__39->SetBinContent(11,3.684496);
   MVA_BDT_B__39->SetBinContent(12,4.065674);
   MVA_BDT_B__39->SetBinContent(13,4.18896);
   MVA_BDT_B__39->SetBinContent(14,4.174456);
   MVA_BDT_B__39->SetBinContent(15,3.889024);
   MVA_BDT_B__39->SetBinContent(16,3.64155);
   MVA_BDT_B__39->SetBinContent(17,3.140305);
   MVA_BDT_B__39->SetBinContent(18,2.504143);
   MVA_BDT_B__39->SetBinContent(19,2.136082);
   MVA_BDT_B__39->SetBinContent(20,1.506474);
   MVA_BDT_B__39->SetBinContent(21,1.215367);
   MVA_BDT_B__39->SetBinContent(22,0.8480747);
   MVA_BDT_B__39->SetBinContent(23,0.6123981);
   MVA_BDT_B__39->SetBinContent(24,0.453518);
   MVA_BDT_B__39->SetBinContent(25,0.2837101);
   MVA_BDT_B__39->SetBinContent(26,0.2227721);
   MVA_BDT_B__39->SetBinContent(27,0.1352761);
   MVA_BDT_B__39->SetBinContent(28,0.1323661);
   MVA_BDT_B__39->SetBinContent(29,0.05021623);
   MVA_BDT_B__39->SetBinContent(30,0.06883308);
   MVA_BDT_B__39->SetBinContent(31,0.06688416);
   MVA_BDT_B__39->SetBinContent(32,0.02173858);
   MVA_BDT_B__39->SetBinContent(33,0.01970785);
   MVA_BDT_B__39->SetBinContent(34,0.01970785);
   MVA_BDT_B__39->SetBinContent(35,0.005037029);
   MVA_BDT_B__39->SetBinContent(36,0.00459676);
   MVA_BDT_B__39->SetBinError(1,0.01331739);
   MVA_BDT_B__39->SetBinError(2,0.02126649);
   MVA_BDT_B__39->SetBinError(3,0.02507636);
   MVA_BDT_B__39->SetBinError(4,0.04131238);
   MVA_BDT_B__39->SetBinError(5,0.05326098);
   MVA_BDT_B__39->SetBinError(6,0.05935876);
   MVA_BDT_B__39->SetBinError(7,0.07775906);
   MVA_BDT_B__39->SetBinError(8,0.08742566);
   MVA_BDT_B__39->SetBinError(9,0.0979686);
   MVA_BDT_B__39->SetBinError(10,0.1049622);
   MVA_BDT_B__39->SetBinError(11,0.1108167);
   MVA_BDT_B__39->SetBinError(12,0.117538);
   MVA_BDT_B__39->SetBinError(13,0.1198301);
   MVA_BDT_B__39->SetBinError(14,0.1206119);
   MVA_BDT_B__39->SetBinError(15,0.1175027);
   MVA_BDT_B__39->SetBinError(16,0.1141215);
   MVA_BDT_B__39->SetBinError(17,0.1066792);
   MVA_BDT_B__39->SetBinError(18,0.09541153);
   MVA_BDT_B__39->SetBinError(19,0.08861321);
   MVA_BDT_B__39->SetBinError(20,0.07411174);
   MVA_BDT_B__39->SetBinError(21,0.0670478);
   MVA_BDT_B__39->SetBinError(22,0.05610103);
   MVA_BDT_B__39->SetBinError(23,0.04831218);
   MVA_BDT_B__39->SetBinError(24,0.04125718);
   MVA_BDT_B__39->SetBinError(25,0.03253247);
   MVA_BDT_B__39->SetBinError(26,0.02898586);
   MVA_BDT_B__39->SetBinError(27,0.02328768);
   MVA_BDT_B__39->SetBinError(28,0.02285904);
   MVA_BDT_B__39->SetBinError(29,0.0133229);
   MVA_BDT_B__39->SetBinError(30,0.01624817);
   MVA_BDT_B__39->SetBinError(31,0.01675409);
   MVA_BDT_B__39->SetBinError(32,0.009299945);
   MVA_BDT_B__39->SetBinError(33,0.009225092);
   MVA_BDT_B__39->SetBinError(34,0.009225092);
   MVA_BDT_B__39->SetBinError(35,0.004617796);
   MVA_BDT_B__39->SetBinError(36,0.00459676);
   MVA_BDT_B__39->SetEntries(41017);
   MVA_BDT_B__39->SetFillColor(2);
   MVA_BDT_B__39->SetFillStyle(3554);
   MVA_BDT_B__39->SetLineColor(2);
   MVA_BDT_B__39->GetXaxis()->SetLabelFont(42);
   MVA_BDT_B__39->GetXaxis()->SetTitleOffset(1);
   MVA_BDT_B__39->GetXaxis()->SetTitleFont(42);
   MVA_BDT_B__39->GetYaxis()->SetLabelFont(42);
   MVA_BDT_B__39->GetYaxis()->SetTitleFont(42);
   MVA_BDT_B__39->GetZaxis()->SetLabelFont(42);
   MVA_BDT_B__39->GetZaxis()->SetTitleOffset(1);
   MVA_BDT_B__39->GetZaxis()->SetTitleFont(42);
   MVA_BDT_B__39->Draw("samehist");
   
   leg = new TLegend(0.53,0.78,0.95,0.9,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1);
   entry=leg->AddEntry("MVA_BDT_Train_S","Signal (training sample)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ee");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(0.7);
   entry->SetTextFont(62);
   entry=leg->AddEntry("MVA_BDT_Train_B","Background (training sample)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(0.7);
   entry->SetTextFont(62);
   leg->Draw();
   
   TH1D *MVA_BDT_Train_S__40 = new TH1D("MVA_BDT_Train_S__40","MVA_BDT_Train_S",40,-0.4516979,0.3238274);
   MVA_BDT_Train_S__40->SetBinContent(13,0.0002338619);
   MVA_BDT_Train_S__40->SetBinContent(14,0.007510628);
   MVA_BDT_Train_S__40->SetBinContent(15,0.002737405);
   MVA_BDT_Train_S__40->SetBinContent(16,0.005007086);
   MVA_BDT_Train_S__40->SetBinContent(17,0.01141734);
   MVA_BDT_Train_S__40->SetBinContent(18,0.01986342);
   MVA_BDT_Train_S__40->SetBinContent(19,0.02493944);
   MVA_BDT_Train_S__40->SetBinContent(20,0.06012692);
   MVA_BDT_Train_S__40->SetBinContent(21,0.09007344);
   MVA_BDT_Train_S__40->SetBinContent(22,0.1412176);
   MVA_BDT_Train_S__40->SetBinContent(23,0.2892024);
   MVA_BDT_Train_S__40->SetBinContent(24,0.487327);
   MVA_BDT_Train_S__40->SetBinContent(25,0.6680646);
   MVA_BDT_Train_S__40->SetBinContent(26,1.155281);
   MVA_BDT_Train_S__40->SetBinContent(27,1.634285);
   MVA_BDT_Train_S__40->SetBinContent(28,2.202633);
   MVA_BDT_Train_S__40->SetBinContent(29,2.874413);
   MVA_BDT_Train_S__40->SetBinContent(30,3.891378);
   MVA_BDT_Train_S__40->SetBinContent(31,4.681925);
   MVA_BDT_Train_S__40->SetBinContent(32,5.56718);
   MVA_BDT_Train_S__40->SetBinContent(33,6.286029);
   MVA_BDT_Train_S__40->SetBinContent(34,6.61932);
   MVA_BDT_Train_S__40->SetBinContent(35,6.169354);
   MVA_BDT_Train_S__40->SetBinContent(36,4.761186);
   MVA_BDT_Train_S__40->SetBinContent(37,2.679383);
   MVA_BDT_Train_S__40->SetBinContent(38,0.9583255);
   MVA_BDT_Train_S__40->SetBinContent(39,0.2623231);
   MVA_BDT_Train_S__40->SetBinContent(40,0.02720912);
   MVA_BDT_Train_S__40->SetBinContent(41,0.002737405);
   MVA_BDT_Train_S__40->SetBinError(13,0.0002338619);
   MVA_BDT_Train_S__40->SetBinError(14,0.004336263);
   MVA_BDT_Train_S__40->SetBinError(15,0.002514442);
   MVA_BDT_Train_S__40->SetBinError(16,0.003540544);
   MVA_BDT_Train_S__40->SetBinError(17,0.005039747);
   MVA_BDT_Train_S__40->SetBinError(18,0.006664908);
   MVA_BDT_Train_S__40->SetBinError(19,0.007161727);
   MVA_BDT_Train_S__40->SetBinError(20,0.01130072);
   MVA_BDT_Train_S__40->SetBinError(21,0.01383951);
   MVA_BDT_Train_S__40->SetBinError(22,0.01748637);
   MVA_BDT_Train_S__40->SetBinError(23,0.025552);
   MVA_BDT_Train_S__40->SetBinError(24,0.03346228);
   MVA_BDT_Train_S__40->SetBinError(25,0.03913231);
   MVA_BDT_Train_S__40->SetBinError(26,0.05176338);
   MVA_BDT_Train_S__40->SetBinError(27,0.06171377);
   MVA_BDT_Train_S__40->SetBinError(28,0.07149858);
   MVA_BDT_Train_S__40->SetBinError(29,0.08140766);
   MVA_BDT_Train_S__40->SetBinError(30,0.09493642);
   MVA_BDT_Train_S__40->SetBinError(31,0.1040606);
   MVA_BDT_Train_S__40->SetBinError(32,0.1134003);
   MVA_BDT_Train_S__40->SetBinError(33,0.12072);
   MVA_BDT_Train_S__40->SetBinError(34,0.1235514);
   MVA_BDT_Train_S__40->SetBinError(35,0.1192328);
   MVA_BDT_Train_S__40->SetBinError(36,0.1048023);
   MVA_BDT_Train_S__40->SetBinError(37,0.07870962);
   MVA_BDT_Train_S__40->SetBinError(38,0.04704148);
   MVA_BDT_Train_S__40->SetBinError(39,0.02463418);
   MVA_BDT_Train_S__40->SetBinError(40,0.007583097);
   MVA_BDT_Train_S__40->SetBinError(41,0.002514442);
   MVA_BDT_Train_S__40->SetEntries(37685);

   ci = TColor::GetColor("#0000ee");
   MVA_BDT_Train_S__40->SetLineColor(ci);

   ci = TColor::GetColor("#0000ee");
   MVA_BDT_Train_S__40->SetMarkerColor(ci);
   MVA_BDT_Train_S__40->SetMarkerStyle(20);
   MVA_BDT_Train_S__40->SetMarkerSize(0.7);
   MVA_BDT_Train_S__40->GetXaxis()->SetLabelFont(42);
   MVA_BDT_Train_S__40->GetXaxis()->SetTitleOffset(1);
   MVA_BDT_Train_S__40->GetXaxis()->SetTitleFont(42);
   MVA_BDT_Train_S__40->GetYaxis()->SetLabelFont(42);
   MVA_BDT_Train_S__40->GetYaxis()->SetTitleFont(42);
   MVA_BDT_Train_S__40->GetZaxis()->SetLabelFont(42);
   MVA_BDT_Train_S__40->GetZaxis()->SetTitleOffset(1);
   MVA_BDT_Train_S__40->GetZaxis()->SetTitleFont(42);
   MVA_BDT_Train_S__40->Draw("e1same");
   
   TH1D *MVA_BDT_Train_B__41 = new TH1D("MVA_BDT_Train_B__41","MVA_BDT_Train_B",40,-0.4516979,0.3238274);
   MVA_BDT_Train_B__41->SetBinContent(0,0.002782228);
   MVA_BDT_Train_B__41->SetBinContent(1,0.02239738);
   MVA_BDT_Train_B__41->SetBinContent(2,0.1167294);
   MVA_BDT_Train_B__41->SetBinContent(3,0.2210872);
   MVA_BDT_Train_B__41->SetBinContent(4,0.5616639);
   MVA_BDT_Train_B__41->SetBinContent(5,1.047687);
   MVA_BDT_Train_B__41->SetBinContent(6,1.480341);
   MVA_BDT_Train_B__41->SetBinContent(7,2.019054);
   MVA_BDT_Train_B__41->SetBinContent(8,2.617655);
   MVA_BDT_Train_B__41->SetBinContent(9,3.121321);
   MVA_BDT_Train_B__41->SetBinContent(10,3.573569);
   MVA_BDT_Train_B__41->SetBinContent(11,3.799289);
   MVA_BDT_Train_B__41->SetBinContent(12,3.950323);
   MVA_BDT_Train_B__41->SetBinContent(13,4.118946);
   MVA_BDT_Train_B__41->SetBinContent(14,4.075423);
   MVA_BDT_Train_B__41->SetBinContent(15,3.714994);
   MVA_BDT_Train_B__41->SetBinContent(16,3.460546);
   MVA_BDT_Train_B__41->SetBinContent(17,3.243587);
   MVA_BDT_Train_B__41->SetBinContent(18,2.567118);
   MVA_BDT_Train_B__41->SetBinContent(19,2.21514);
   MVA_BDT_Train_B__41->SetBinContent(20,1.693184);
   MVA_BDT_Train_B__41->SetBinContent(21,1.211849);
   MVA_BDT_Train_B__41->SetBinContent(22,1.004345);
   MVA_BDT_Train_B__41->SetBinContent(23,0.6679799);
   MVA_BDT_Train_B__41->SetBinContent(24,0.4448779);
   MVA_BDT_Train_B__41->SetBinContent(25,0.207179);
   MVA_BDT_Train_B__41->SetBinContent(26,0.1767243);
   MVA_BDT_Train_B__41->SetBinContent(27,0.1137972);
   MVA_BDT_Train_B__41->SetBinContent(28,0.03389823);
   MVA_BDT_Train_B__41->SetBinContent(29,0.03070127);
   MVA_BDT_Train_B__41->SetBinContent(30,0.01247141);
   MVA_BDT_Train_B__41->SetBinContent(31,0.005557301);
   MVA_BDT_Train_B__41->SetBinContent(32,0.02942108);
   MVA_BDT_Train_B__41->SetBinContent(33,0.01282779);
   MVA_BDT_Train_B__41->SetBinContent(34,0.006258419);
   MVA_BDT_Train_B__41->SetBinError(0,0.001488388);
   MVA_BDT_Train_B__41->SetBinError(1,0.00735972);
   MVA_BDT_Train_B__41->SetBinError(2,0.01995509);
   MVA_BDT_Train_B__41->SetBinError(3,0.0280682);
   MVA_BDT_Train_B__41->SetBinError(4,0.040011);
   MVA_BDT_Train_B__41->SetBinError(5,0.0521596);
   MVA_BDT_Train_B__41->SetBinError(6,0.0649908);
   MVA_BDT_Train_B__41->SetBinError(7,0.07791749);
   MVA_BDT_Train_B__41->SetBinError(8,0.08874482);
   MVA_BDT_Train_B__41->SetBinError(9,0.09774518);
   MVA_BDT_Train_B__41->SetBinError(10,0.1073388);
   MVA_BDT_Train_B__41->SetBinError(11,0.1138363);
   MVA_BDT_Train_B__41->SetBinError(12,0.1156305);
   MVA_BDT_Train_B__41->SetBinError(13,0.1200407);
   MVA_BDT_Train_B__41->SetBinError(14,0.1197569);
   MVA_BDT_Train_B__41->SetBinError(15,0.1153106);
   MVA_BDT_Train_B__41->SetBinError(16,0.1108542);
   MVA_BDT_Train_B__41->SetBinError(17,0.1085328);
   MVA_BDT_Train_B__41->SetBinError(18,0.09820173);
   MVA_BDT_Train_B__41->SetBinError(19,0.09137316);
   MVA_BDT_Train_B__41->SetBinError(20,0.08011801);
   MVA_BDT_Train_B__41->SetBinError(21,0.06769672);
   MVA_BDT_Train_B__41->SetBinError(22,0.06268502);
   MVA_BDT_Train_B__41->SetBinError(23,0.05015552);
   MVA_BDT_Train_B__41->SetBinError(24,0.04058786);
   MVA_BDT_Train_B__41->SetBinError(25,0.02515722);
   MVA_BDT_Train_B__41->SetBinError(26,0.02565002);
   MVA_BDT_Train_B__41->SetBinError(27,0.02066917);
   MVA_BDT_Train_B__41->SetBinError(28,0.01068916);
   MVA_BDT_Train_B__41->SetBinError(29,0.009829841);
   MVA_BDT_Train_B__41->SetBinError(30,0.005223831);
   MVA_BDT_Train_B__41->SetBinError(31,0.004706491);
   MVA_BDT_Train_B__41->SetBinError(32,0.01147322);
   MVA_BDT_Train_B__41->SetBinError(33,0.006751434);
   MVA_BDT_Train_B__41->SetBinError(34,0.004823796);
   MVA_BDT_Train_B__41->SetEntries(41017);
   MVA_BDT_Train_B__41->SetLineColor(2);
   MVA_BDT_Train_B__41->SetMarkerColor(2);
   MVA_BDT_Train_B__41->SetMarkerStyle(20);
   MVA_BDT_Train_B__41->SetMarkerSize(0.7);
   MVA_BDT_Train_B__41->GetXaxis()->SetLabelFont(42);
   MVA_BDT_Train_B__41->GetXaxis()->SetTitleOffset(1);
   MVA_BDT_Train_B__41->GetXaxis()->SetTitleFont(42);
   MVA_BDT_Train_B__41->GetYaxis()->SetLabelFont(42);
   MVA_BDT_Train_B__41->GetYaxis()->SetTitleFont(42);
   MVA_BDT_Train_B__41->GetZaxis()->SetLabelFont(42);
   MVA_BDT_Train_B__41->GetZaxis()->SetTitleOffset(1);
   MVA_BDT_Train_B__41->GetZaxis()->SetTitleFont(42);
   MVA_BDT_Train_B__41->Draw("e1same");
   TText *text = new TText(0.12,0.74,"Kolmogorov-Smirnov test: signal (background) probability = 0.347 (0.464)");
   text->SetNDC();
   text->SetTextSize(0.032);
   text->Draw();
   
   TH2F *frameBDT__42 = new TH2F("frameBDT__42","TMVA overtraining check for classifier: BDT",500,-0.4516979,0.3238274,500,0,8.605116);
   frameBDT__42->SetStats(0);
   frameBDT__42->GetXaxis()->SetTitle("BDT response");
   frameBDT__42->GetXaxis()->SetLabelOffset(0.012);
   frameBDT__42->GetXaxis()->SetLabelSize(0.04);
   frameBDT__42->GetXaxis()->SetTitleSize(0.045);
   frameBDT__42->GetXaxis()->SetTitleOffset(1.25);
   frameBDT__42->GetYaxis()->SetTitle("(1/N) dN^{ }/^{ }dx");
   frameBDT__42->GetYaxis()->SetLabelOffset(0.012);
   frameBDT__42->GetYaxis()->SetLabelSize(0.04);
   frameBDT__42->GetYaxis()->SetTitleSize(0.045);
   frameBDT__42->GetYaxis()->SetTitleOffset(1.2);
   frameBDT__42->GetZaxis()->SetLabelSize(0.04);
   frameBDT__42->GetZaxis()->SetTitleSize(0.04);
   frameBDT__42->GetZaxis()->SetTitleOffset(1);
   frameBDT__42->Draw("sameaxis");
   text = new TText(0.975,0.115,"U/O-flow (S,B): (0.0, 0.0)% / (0.0, 0.0)%");
   text->SetNDC();
   text->SetTextSize(0.03);
   text->SetTextAngle(90);
   text->Draw();
   
   TPaveText *pt = new TPaveText(0.01,0.9382432,0.71,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(1);
   pt->SetFillColor(0);
   TText *pt_LaTex = pt->AddText("TMVA overtraining check for classifier: BDT");
   pt->Draw();
   canvas1->Modified();
   canvas1->SetSelected(canvas1);
}
