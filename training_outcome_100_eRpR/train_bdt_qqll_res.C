#ifdef __CLING__
#pragma cling optimize(0)
#endif
void train_bdt_qqll_res()
{
//=========Macro generated from canvas: canvas1/TMVA comparison BDT
//=========  (Wed Aug 14 11:46:30 2024) by ROOT version 6.32.02
   TCanvas *canvas1 = new TCanvas("canvas1", "TMVA comparison BDT",200,81,600,468);
   gStyle->SetOptStat(0);
   canvas1->Range(-0.6750363,-1.424388,0.3764648,10.44551);
   canvas1->SetFillColor(0);
   canvas1->SetBorderMode(0);
   canvas1->SetBorderSize(2);
   canvas1->SetTickx(1);
   canvas1->SetTicky(1);
   canvas1->SetRightMargin(0.05);
   canvas1->SetBottomMargin(0.12);
   canvas1->SetFrameBorderMode(0);
   canvas1->SetFrameBorderMode(0);
   
   TH2F *frameBDT__37 = new TH2F("frameBDT__37","TMVA overtraining check for classifier: BDT",500,-0.5646287,0.3238898,500,0,9.258521);
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
   
   TH1D *MVA_BDT_S__38 = new TH1D("MVA_BDT_S__38","TMVA overtraining check for classifier: BDT",40,-0.5646287,0.3238897);
   MVA_BDT_S__38->SetBinContent(17,0.004678393);
   MVA_BDT_S__38->SetBinContent(18,0.002339197);
   MVA_BDT_S__38->SetBinContent(19,0.02474794);
   MVA_BDT_S__38->SetBinContent(20,0.02474794);
   MVA_BDT_S__38->SetBinContent(21,0.05575918);
   MVA_BDT_S__38->SetBinContent(22,0.1147645);
   MVA_BDT_S__38->SetBinContent(23,0.1787446);
   MVA_BDT_S__38->SetBinContent(24,0.3762791);
   MVA_BDT_S__38->SetBinContent(25,0.6612531);
   MVA_BDT_S__38->SetBinContent(26,1.212802);
   MVA_BDT_S__38->SetBinContent(27,2.027814);
   MVA_BDT_S__38->SetBinContent(28,3.076822);
   MVA_BDT_S__38->SetBinContent(29,4.257679);
   MVA_BDT_S__38->SetBinContent(30,5.736251);
   MVA_BDT_S__38->SetBinContent(31,6.822525);
   MVA_BDT_S__38->SetBinContent(32,6.844535);
   MVA_BDT_S__38->SetBinContent(33,5.868194);
   MVA_BDT_S__38->SetBinContent(34,3.973324);
   MVA_BDT_S__38->SetBinContent(35,2.25063);
   MVA_BDT_S__38->SetBinContent(36,0.9524142);
   MVA_BDT_S__38->SetBinContent(37,0.4395812);
   MVA_BDT_S__38->SetBinContent(38,0.07960021);
   MVA_BDT_S__38->SetBinContent(39,0.02633284);
   MVA_BDT_S__38->SetBinContent(40,0.006941282);
   MVA_BDT_S__38->SetBinError(17,0.003308123);
   MVA_BDT_S__38->SetBinError(18,0.002339197);
   MVA_BDT_S__38->SetBinError(19,0.006951721);
   MVA_BDT_S__38->SetBinError(20,0.006951721);
   MVA_BDT_S__38->SetBinError(21,0.01060657);
   MVA_BDT_S__38->SetBinError(22,0.0154457);
   MVA_BDT_S__38->SetBinError(23,0.01880326);
   MVA_BDT_S__38->SetBinError(24,0.02747116);
   MVA_BDT_S__38->SetBinError(25,0.03594667);
   MVA_BDT_S__38->SetBinError(26,0.04873935);
   MVA_BDT_S__38->SetBinError(27,0.06264622);
   MVA_BDT_S__38->SetBinError(28,0.07761779);
   MVA_BDT_S__38->SetBinError(29,0.09064509);
   MVA_BDT_S__38->SetBinError(30,0.1060416);
   MVA_BDT_S__38->SetBinError(31,0.1157054);
   MVA_BDT_S__38->SetBinError(32,0.1154886);
   MVA_BDT_S__38->SetBinError(33,0.10684);
   MVA_BDT_S__38->SetBinError(34,0.08801269);
   MVA_BDT_S__38->SetBinError(35,0.066335);
   MVA_BDT_S__38->SetBinError(36,0.04291553);
   MVA_BDT_S__38->SetBinError(37,0.0295441);
   MVA_BDT_S__38->SetBinError(38,0.01235828);
   MVA_BDT_S__38->SetBinError(39,0.007295841);
   MVA_BDT_S__38->SetBinError(40,0.003556764);
   MVA_BDT_S__38->SetEntries(29128);
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
   
   TH1D *MVA_BDT_B__39 = new TH1D("MVA_BDT_B__39","MVA_BDT_B",40,-0.5646287,0.3238897);
   MVA_BDT_B__39->SetBinContent(1,0.003612027);
   MVA_BDT_B__39->SetBinContent(2,0.04010981);
   MVA_BDT_B__39->SetBinContent(3,0.08820691);
   MVA_BDT_B__39->SetBinContent(4,0.2096719);
   MVA_BDT_B__39->SetBinContent(5,0.410757);
   MVA_BDT_B__39->SetBinContent(6,0.8780825);
   MVA_BDT_B__39->SetBinContent(7,1.342482);
   MVA_BDT_B__39->SetBinContent(8,1.878122);
   MVA_BDT_B__39->SetBinContent(9,2.383557);
   MVA_BDT_B__39->SetBinContent(10,2.985948);
   MVA_BDT_B__39->SetBinContent(11,3.367164);
   MVA_BDT_B__39->SetBinContent(12,3.192424);
   MVA_BDT_B__39->SetBinContent(13,3.029155);
   MVA_BDT_B__39->SetBinContent(14,2.568829);
   MVA_BDT_B__39->SetBinContent(15,2.465264);
   MVA_BDT_B__39->SetBinContent(16,2.271878);
   MVA_BDT_B__39->SetBinContent(17,2.176605);
   MVA_BDT_B__39->SetBinContent(18,2.038371);
   MVA_BDT_B__39->SetBinContent(19,2.081606);
   MVA_BDT_B__39->SetBinContent(20,1.822763);
   MVA_BDT_B__39->SetBinContent(21,1.611138);
   MVA_BDT_B__39->SetBinContent(22,1.569903);
   MVA_BDT_B__39->SetBinContent(23,1.249653);
   MVA_BDT_B__39->SetBinContent(24,1.074848);
   MVA_BDT_B__39->SetBinContent(25,0.9670412);
   MVA_BDT_B__39->SetBinContent(26,0.8462566);
   MVA_BDT_B__39->SetBinContent(27,0.7168967);
   MVA_BDT_B__39->SetBinContent(28,0.6095253);
   MVA_BDT_B__39->SetBinContent(29,0.4166336);
   MVA_BDT_B__39->SetBinContent(30,0.2892279);
   MVA_BDT_B__39->SetBinContent(31,0.2186168);
   MVA_BDT_B__39->SetBinContent(32,0.1194129);
   MVA_BDT_B__39->SetBinContent(33,0.05900697);
   MVA_BDT_B__39->SetBinContent(34,0.0224225);
   MVA_BDT_B__39->SetBinContent(35,0.009126739);
   MVA_BDT_B__39->SetBinContent(36,0.003557668);
   MVA_BDT_B__39->SetBinContent(38,0.0008834068);
   MVA_BDT_B__39->SetBinError(1,0.001725613);
   MVA_BDT_B__39->SetBinError(2,0.008042783);
   MVA_BDT_B__39->SetBinError(3,0.01200801);
   MVA_BDT_B__39->SetBinError(4,0.01875868);
   MVA_BDT_B__39->SetBinError(5,0.0262978);
   MVA_BDT_B__39->SetBinError(6,0.03754018);
   MVA_BDT_B__39->SetBinError(7,0.04588369);
   MVA_BDT_B__39->SetBinError(8,0.0560169);
   MVA_BDT_B__39->SetBinError(9,0.0646028);
   MVA_BDT_B__39->SetBinError(10,0.07439863);
   MVA_BDT_B__39->SetBinError(11,0.07999493);
   MVA_BDT_B__39->SetBinError(12,0.07772543);
   MVA_BDT_B__39->SetBinError(13,0.07731967);
   MVA_BDT_B__39->SetBinError(14,0.07126086);
   MVA_BDT_B__39->SetBinError(15,0.07055279);
   MVA_BDT_B__39->SetBinError(16,0.06760488);
   MVA_BDT_B__39->SetBinError(17,0.06574404);
   MVA_BDT_B__39->SetBinError(18,0.06387796);
   MVA_BDT_B__39->SetBinError(19,0.0651914);
   MVA_BDT_B__39->SetBinError(20,0.06043738);
   MVA_BDT_B__39->SetBinError(21,0.05658671);
   MVA_BDT_B__39->SetBinError(22,0.05600781);
   MVA_BDT_B__39->SetBinError(23,0.04958458);
   MVA_BDT_B__39->SetBinError(24,0.04591238);
   MVA_BDT_B__39->SetBinError(25,0.04366134);
   MVA_BDT_B__39->SetBinError(26,0.04084942);
   MVA_BDT_B__39->SetBinError(27,0.03777993);
   MVA_BDT_B__39->SetBinError(28,0.03569664);
   MVA_BDT_B__39->SetBinError(29,0.028988);
   MVA_BDT_B__39->SetBinError(30,0.02369292);
   MVA_BDT_B__39->SetBinError(31,0.02122865);
   MVA_BDT_B__39->SetBinError(32,0.015902);
   MVA_BDT_B__39->SetBinError(33,0.009864834);
   MVA_BDT_B__39->SetBinError(34,0.006227445);
   MVA_BDT_B__39->SetBinError(35,0.004462051);
   MVA_BDT_B__39->SetBinError(36,0.002816395);
   MVA_BDT_B__39->SetBinError(38,0.0008834068);
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
   
   TH1D *MVA_BDT_Train_S__40 = new TH1D("MVA_BDT_Train_S__40","MVA_BDT_Train_S",40,-0.5646287,0.3238897);
   MVA_BDT_Train_S__40->SetBinContent(18,0.00459101);
   MVA_BDT_Train_S__40->SetBinContent(19,0.007677058);
   MVA_BDT_Train_S__40->SetBinContent(20,0.02002125);
   MVA_BDT_Train_S__40->SetBinContent(21,0.02702195);
   MVA_BDT_Train_S__40->SetBinContent(22,0.06488313);
   MVA_BDT_Train_S__40->SetBinContent(23,0.1317192);
   MVA_BDT_Train_S__40->SetBinContent(24,0.3089005);
   MVA_BDT_Train_S__40->SetBinContent(25,0.5801009);
   MVA_BDT_Train_S__40->SetBinContent(26,1.200549);
   MVA_BDT_Train_S__40->SetBinContent(27,1.988911);
   MVA_BDT_Train_S__40->SetBinContent(28,3.159056);
   MVA_BDT_Train_S__40->SetBinContent(29,4.437052);
   MVA_BDT_Train_S__40->SetBinContent(30,5.873891);
   MVA_BDT_Train_S__40->SetBinContent(31,6.732252);
   MVA_BDT_Train_S__40->SetBinContent(32,7.12194);
   MVA_BDT_Train_S__40->SetBinContent(33,5.754042);
   MVA_BDT_Train_S__40->SetBinContent(34,3.944587);
   MVA_BDT_Train_S__40->SetBinContent(35,2.144136);
   MVA_BDT_Train_S__40->SetBinContent(36,1.00108);
   MVA_BDT_Train_S__40->SetBinContent(37,0.3799645);
   MVA_BDT_Train_S__40->SetBinContent(38,0.1102691);
   MVA_BDT_Train_S__40->SetBinContent(39,0.02461226);
   MVA_BDT_Train_S__40->SetBinContent(40,0.001504962);
   MVA_BDT_Train_S__40->SetBinContent(41,0.002333567);
   MVA_BDT_Train_S__40->SetBinError(18,0.002672867);
   MVA_BDT_Train_S__40->SetBinError(19,0.003627117);
   MVA_BDT_Train_S__40->SetBinError(20,0.006099429);
   MVA_BDT_Train_S__40->SetBinError(21,0.007317078);
   MVA_BDT_Train_S__40->SetBinError(22,0.01140149);
   MVA_BDT_Train_S__40->SetBinError(23,0.01570665);
   MVA_BDT_Train_S__40->SetBinError(24,0.02406785);
   MVA_BDT_Train_S__40->SetBinError(25,0.03356284);
   MVA_BDT_Train_S__40->SetBinError(26,0.04839755);
   MVA_BDT_Train_S__40->SetBinError(27,0.06193388);
   MVA_BDT_Train_S__40->SetBinError(28,0.07836129);
   MVA_BDT_Train_S__40->SetBinError(29,0.09279651);
   MVA_BDT_Train_S__40->SetBinError(30,0.1072221);
   MVA_BDT_Train_S__40->SetBinError(31,0.1145274);
   MVA_BDT_Train_S__40->SetBinError(32,0.1181191);
   MVA_BDT_Train_S__40->SetBinError(33,0.1057244);
   MVA_BDT_Train_S__40->SetBinError(34,0.08803357);
   MVA_BDT_Train_S__40->SetBinError(35,0.06494812);
   MVA_BDT_Train_S__40->SetBinError(36,0.04397677);
   MVA_BDT_Train_S__40->SetBinError(37,0.02753509);
   MVA_BDT_Train_S__40->SetBinError(38,0.01456401);
   MVA_BDT_Train_S__40->SetBinError(39,0.006659373);
   MVA_BDT_Train_S__40->SetBinError(40,0.001064169);
   MVA_BDT_Train_S__40->SetBinError(41,0.002333567);
   MVA_BDT_Train_S__40->SetEntries(29128);

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
   
   TH1D *MVA_BDT_Train_B__41 = new TH1D("MVA_BDT_Train_B__41","MVA_BDT_Train_B",40,-0.5646287,0.3238897);
   MVA_BDT_Train_B__41->SetBinContent(0,0.002019816);
   MVA_BDT_Train_B__41->SetBinContent(1,0.009981358);
   MVA_BDT_Train_B__41->SetBinContent(2,0.02418282);
   MVA_BDT_Train_B__41->SetBinContent(3,0.07749132);
   MVA_BDT_Train_B__41->SetBinContent(4,0.1949901);
   MVA_BDT_Train_B__41->SetBinContent(5,0.4067462);
   MVA_BDT_Train_B__41->SetBinContent(6,0.8642038);
   MVA_BDT_Train_B__41->SetBinContent(7,1.453929);
   MVA_BDT_Train_B__41->SetBinContent(8,1.888195);
   MVA_BDT_Train_B__41->SetBinContent(9,2.537549);
   MVA_BDT_Train_B__41->SetBinContent(10,3.027668);
   MVA_BDT_Train_B__41->SetBinContent(11,3.239716);
   MVA_BDT_Train_B__41->SetBinContent(12,3.227477);
   MVA_BDT_Train_B__41->SetBinContent(13,2.92646);
   MVA_BDT_Train_B__41->SetBinContent(14,2.620656);
   MVA_BDT_Train_B__41->SetBinContent(15,2.386503);
   MVA_BDT_Train_B__41->SetBinContent(16,2.283301);
   MVA_BDT_Train_B__41->SetBinContent(17,2.22387);
   MVA_BDT_Train_B__41->SetBinContent(18,2.060358);
   MVA_BDT_Train_B__41->SetBinContent(19,1.947023);
   MVA_BDT_Train_B__41->SetBinContent(20,1.769255);
   MVA_BDT_Train_B__41->SetBinContent(21,1.597729);
   MVA_BDT_Train_B__41->SetBinContent(22,1.494872);
   MVA_BDT_Train_B__41->SetBinContent(23,1.371318);
   MVA_BDT_Train_B__41->SetBinContent(24,1.130973);
   MVA_BDT_Train_B__41->SetBinContent(25,0.9970335);
   MVA_BDT_Train_B__41->SetBinContent(26,1.007677);
   MVA_BDT_Train_B__41->SetBinContent(27,0.7667653);
   MVA_BDT_Train_B__41->SetBinContent(28,0.6458724);
   MVA_BDT_Train_B__41->SetBinContent(29,0.3754202);
   MVA_BDT_Train_B__41->SetBinContent(30,0.2272537);
   MVA_BDT_Train_B__41->SetBinContent(31,0.1218073);
   MVA_BDT_Train_B__41->SetBinContent(32,0.07148675);
   MVA_BDT_Train_B__41->SetBinContent(33,0.02774681);
   MVA_BDT_Train_B__41->SetBinContent(34,0.0114755);
   MVA_BDT_Train_B__41->SetBinContent(35,0.001774203);
   MVA_BDT_Train_B__41->SetBinError(0,0.002019816);
   MVA_BDT_Train_B__41->SetBinError(1,0.003729037);
   MVA_BDT_Train_B__41->SetBinError(2,0.006157236);
   MVA_BDT_Train_B__41->SetBinError(3,0.01056334);
   MVA_BDT_Train_B__41->SetBinError(4,0.01741486);
   MVA_BDT_Train_B__41->SetBinError(5,0.02605737);
   MVA_BDT_Train_B__41->SetBinError(6,0.03747436);
   MVA_BDT_Train_B__41->SetBinError(7,0.04938287);
   MVA_BDT_Train_B__41->SetBinError(8,0.056497);
   MVA_BDT_Train_B__41->SetBinError(9,0.06750544);
   MVA_BDT_Train_B__41->SetBinError(10,0.07431287);
   MVA_BDT_Train_B__41->SetBinError(11,0.07827737);
   MVA_BDT_Train_B__41->SetBinError(12,0.07892992);
   MVA_BDT_Train_B__41->SetBinError(13,0.07619544);
   MVA_BDT_Train_B__41->SetBinError(14,0.07229426);
   MVA_BDT_Train_B__41->SetBinError(15,0.06959621);
   MVA_BDT_Train_B__41->SetBinError(16,0.06691647);
   MVA_BDT_Train_B__41->SetBinError(17,0.06684109);
   MVA_BDT_Train_B__41->SetBinError(18,0.0639024);
   MVA_BDT_Train_B__41->SetBinError(19,0.0627826);
   MVA_BDT_Train_B__41->SetBinError(20,0.05905265);
   MVA_BDT_Train_B__41->SetBinError(21,0.05650949);
   MVA_BDT_Train_B__41->SetBinError(22,0.05460865);
   MVA_BDT_Train_B__41->SetBinError(23,0.05288578);
   MVA_BDT_Train_B__41->SetBinError(24,0.04786311);
   MVA_BDT_Train_B__41->SetBinError(25,0.04496678);
   MVA_BDT_Train_B__41->SetBinError(26,0.045912);
   MVA_BDT_Train_B__41->SetBinError(27,0.03991564);
   MVA_BDT_Train_B__41->SetBinError(28,0.03624795);
   MVA_BDT_Train_B__41->SetBinError(29,0.02694508);
   MVA_BDT_Train_B__41->SetBinError(30,0.02052781);
   MVA_BDT_Train_B__41->SetBinError(31,0.01504114);
   MVA_BDT_Train_B__41->SetBinError(32,0.01232149);
   MVA_BDT_Train_B__41->SetBinError(33,0.007639167);
   MVA_BDT_Train_B__41->SetBinError(34,0.004426782);
   MVA_BDT_Train_B__41->SetBinError(35,0.001254551);
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
   TText *text = new TText(0.12,0.74,"Kolmogorov-Smirnov test: signal (background) probability = 0.206 (0.453)");
   text->SetNDC();
   text->SetTextSize(0.032);
   text->Draw();
   
   TH2F *frameBDT__42 = new TH2F("frameBDT__42","TMVA overtraining check for classifier: BDT",500,-0.5646287,0.3238898,500,0,9.258521);
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
