#ifdef __CLING__
#pragma cling optimize(0)
#endif
void train_bdt_qqll_res()
{
//=========Macro generated from canvas: canvas1/TMVA comparison BDT
//=========  (Sun Jun 23 22:07:18 2024) by ROOT version 6.32.00
   TCanvas *canvas1 = new TCanvas("canvas1", "TMVA comparison BDT",200,81,600,468);
   gStyle->SetOptStat(0);
   canvas1->Range(-0.6206534,-1.08173,0.4142477,7.932683);
   canvas1->SetFillColor(0);
   canvas1->SetBorderMode(0);
   canvas1->SetBorderSize(2);
   canvas1->SetTickx(1);
   canvas1->SetTicky(1);
   canvas1->SetRightMargin(0.05);
   canvas1->SetBottomMargin(0.12);
   canvas1->SetFrameBorderMode(0);
   canvas1->SetFrameBorderMode(0);
   
   TH2F *frameBDT__37 = new TH2F("frameBDT__37","TMVA overtraining check for classifier: BDT",500,-0.5119888,0.3625026,500,0,7.031242);
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
   
   TH1D *MVA_BDT_S__38 = new TH1D("MVA_BDT_S__38","TMVA overtraining check for classifier: BDT",40,-0.5119888,0.3625026);
   MVA_BDT_S__38->SetBinContent(15,0.003965284);
   MVA_BDT_S__38->SetBinContent(16,0.003965284);
   MVA_BDT_S__38->SetBinContent(17,0.01453938);
   MVA_BDT_S__38->SetBinContent(18,0.01982642);
   MVA_BDT_S__38->SetBinContent(19,0.03568756);
   MVA_BDT_S__38->SetBinContent(20,0.09384506);
   MVA_BDT_S__38->SetBinContent(21,0.1744725);
   MVA_BDT_S__38->SetBinContent(22,0.260387);
   MVA_BDT_S__38->SetBinContent(23,0.4110678);
   MVA_BDT_S__38->SetBinContent(24,0.6410543);
   MVA_BDT_S__38->SetBinContent(25,0.9781035);
   MVA_BDT_S__38->SetBinContent(26,1.477729);
   MVA_BDT_S__38->SetBinContent(27,1.912589);
   MVA_BDT_S__38->SetBinContent(28,2.549678);
   MVA_BDT_S__38->SetBinContent(29,3.255498);
   MVA_BDT_S__38->SetBinContent(30,4.212454);
   MVA_BDT_S__38->SetBinContent(31,4.953962);
   MVA_BDT_S__38->SetBinContent(32,5.408648);
   MVA_BDT_S__38->SetBinContent(33,5.268541);
   MVA_BDT_S__38->SetBinContent(34,4.579903);
   MVA_BDT_S__38->SetBinContent(35,3.804029);
   MVA_BDT_S__38->SetBinContent(36,2.848396);
   MVA_BDT_S__38->SetBinContent(37,1.732829);
   MVA_BDT_S__38->SetBinContent(38,0.8406403);
   MVA_BDT_S__38->SetBinContent(39,0.2260212);
   MVA_BDT_S__38->SetBinContent(40,0.03304404);
   MVA_BDT_S__38->SetBinError(15,0.002289358);
   MVA_BDT_S__38->SetBinError(16,0.002289358);
   MVA_BDT_S__38->SetBinError(17,0.004383787);
   MVA_BDT_S__38->SetBinError(18,0.00511916);
   MVA_BDT_S__38->SetBinError(19,0.006868074);
   MVA_BDT_S__38->SetBinError(20,0.01113736);
   MVA_BDT_S__38->SetBinError(21,0.01518588);
   MVA_BDT_S__38->SetBinError(22,0.01855181);
   MVA_BDT_S__38->SetBinError(23,0.02330952);
   MVA_BDT_S__38->SetBinError(24,0.02910878);
   MVA_BDT_S__38->SetBinError(25,0.0359558);
   MVA_BDT_S__38->SetBinError(26,0.04419509);
   MVA_BDT_S__38->SetBinError(27,0.05027908);
   MVA_BDT_S__38->SetBinError(28,0.05805227);
   MVA_BDT_S__38->SetBinError(29,0.0655972);
   MVA_BDT_S__38->SetBinError(30,0.07461809);
   MVA_BDT_S__38->SetBinError(31,0.08091944);
   MVA_BDT_S__38->SetBinError(32,0.08455142);
   MVA_BDT_S__38->SetBinError(33,0.08344911);
   MVA_BDT_S__38->SetBinError(34,0.0778045);
   MVA_BDT_S__38->SetBinError(35,0.07090853);
   MVA_BDT_S__38->SetBinError(36,0.06135878);
   MVA_BDT_S__38->SetBinError(37,0.04785799);
   MVA_BDT_S__38->SetBinError(38,0.03333356);
   MVA_BDT_S__38->SetBinError(39,0.01728427);
   MVA_BDT_S__38->SetBinError(40,0.006608807);
   MVA_BDT_S__38->SetEntries(34606);
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
   
   TH1D *MVA_BDT_B__39 = new TH1D("MVA_BDT_B__39","MVA_BDT_B",40,-0.5119888,0.3625026);
   MVA_BDT_B__39->SetBinContent(1,0.01855556);
   MVA_BDT_B__39->SetBinContent(2,0.07497438);
   MVA_BDT_B__39->SetBinContent(3,0.1498565);
   MVA_BDT_B__39->SetBinContent(4,0.2548997);
   MVA_BDT_B__39->SetBinContent(5,0.49662);
   MVA_BDT_B__39->SetBinContent(6,1.052004);
   MVA_BDT_B__39->SetBinContent(7,1.685776);
   MVA_BDT_B__39->SetBinContent(8,2.394806);
   MVA_BDT_B__39->SetBinContent(9,2.760369);
   MVA_BDT_B__39->SetBinContent(10,2.796048);
   MVA_BDT_B__39->SetBinContent(11,2.899339);
   MVA_BDT_B__39->SetBinContent(12,2.917945);
   MVA_BDT_B__39->SetBinContent(13,3.038367);
   MVA_BDT_B__39->SetBinContent(14,3.335707);
   MVA_BDT_B__39->SetBinContent(15,3.464921);
   MVA_BDT_B__39->SetBinContent(16,3.215286);
   MVA_BDT_B__39->SetBinContent(17,2.863844);
   MVA_BDT_B__39->SetBinContent(18,2.549851);
   MVA_BDT_B__39->SetBinContent(19,2.239447);
   MVA_BDT_B__39->SetBinContent(20,2.086757);
   MVA_BDT_B__39->SetBinContent(21,1.550434);
   MVA_BDT_B__39->SetBinContent(22,1.225372);
   MVA_BDT_B__39->SetBinContent(23,0.9413076);
   MVA_BDT_B__39->SetBinContent(24,0.6061106);
   MVA_BDT_B__39->SetBinContent(25,0.4220033);
   MVA_BDT_B__39->SetBinContent(26,0.2401561);
   MVA_BDT_B__39->SetBinContent(27,0.1686678);
   MVA_BDT_B__39->SetBinContent(28,0.1014604);
   MVA_BDT_B__39->SetBinContent(29,0.06258549);
   MVA_BDT_B__39->SetBinContent(30,0.03862377);
   MVA_BDT_B__39->SetBinContent(31,0.02856802);
   MVA_BDT_B__39->SetBinContent(32,0.02593666);
   MVA_BDT_B__39->SetBinContent(33,0.009268743);
   MVA_BDT_B__39->SetBinContent(34,0.01203996);
   MVA_BDT_B__39->SetBinContent(35,0.007410744);
   MVA_BDT_B__39->SetBinContent(36,0.005556394);
   MVA_BDT_B__39->SetBinError(1,0.005958032);
   MVA_BDT_B__39->SetBinError(2,0.01164854);
   MVA_BDT_B__39->SetBinError(3,0.01642226);
   MVA_BDT_B__39->SetBinError(4,0.02038195);
   MVA_BDT_B__39->SetBinError(5,0.02795819);
   MVA_BDT_B__39->SetBinError(6,0.03923636);
   MVA_BDT_B__39->SetBinError(7,0.04913642);
   MVA_BDT_B__39->SetBinError(8,0.05845139);
   MVA_BDT_B__39->SetBinError(9,0.06452286);
   MVA_BDT_B__39->SetBinError(10,0.06681413);
   MVA_BDT_B__39->SetBinError(11,0.06916389);
   MVA_BDT_B__39->SetBinError(12,0.07008569);
   MVA_BDT_B__39->SetBinError(13,0.07215162);
   MVA_BDT_B__39->SetBinError(14,0.07725043);
   MVA_BDT_B__39->SetBinError(15,0.07956629);
   MVA_BDT_B__39->SetBinError(16,0.07687699);
   MVA_BDT_B__39->SetBinError(17,0.07319286);
   MVA_BDT_B__39->SetBinError(18,0.0694717);
   MVA_BDT_B__39->SetBinError(19,0.06547029);
   MVA_BDT_B__39->SetBinError(20,0.06352128);
   MVA_BDT_B__39->SetBinError(21,0.05457084);
   MVA_BDT_B__39->SetBinError(22,0.04897769);
   MVA_BDT_B__39->SetBinError(23,0.04312987);
   MVA_BDT_B__39->SetBinError(24,0.03430935);
   MVA_BDT_B__39->SetBinError(25,0.02832481);
   MVA_BDT_B__39->SetBinError(26,0.02150207);
   MVA_BDT_B__39->SetBinError(27,0.01827825);
   MVA_BDT_B__39->SetBinError(28,0.01376274);
   MVA_BDT_B__39->SetBinError(29,0.01155762);
   MVA_BDT_B__39->SetBinError(30,0.008645652);
   MVA_BDT_B__39->SetBinError(31,0.007084706);
   MVA_BDT_B__39->SetBinError(32,0.006972072);
   MVA_BDT_B__39->SetBinError(33,0.004002318);
   MVA_BDT_B__39->SetBinError(34,0.00504476);
   MVA_BDT_B__39->SetBinError(35,0.004278595);
   MVA_BDT_B__39->SetBinError(36,0.0032961);
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
   
   TH1D *MVA_BDT_Train_S__40 = new TH1D("MVA_BDT_Train_S__40","MVA_BDT_Train_S",40,-0.5119888,0.3625026);
   MVA_BDT_Train_S__40->SetBinContent(12,0.001321838);
   MVA_BDT_Train_S__40->SetBinContent(13,0.001321838);
   MVA_BDT_Train_S__40->SetBinContent(15,0.003965513);
   MVA_BDT_Train_S__40->SetBinContent(16,0.006609189);
   MVA_BDT_Train_S__40->SetBinContent(17,0.005287351);
   MVA_BDT_Train_S__40->SetBinContent(18,0.009252865);
   MVA_BDT_Train_S__40->SetBinContent(19,0.02511492);
   MVA_BDT_Train_S__40->SetBinContent(20,0.05022984);
   MVA_BDT_Train_S__40->SetBinContent(21,0.09120681);
   MVA_BDT_Train_S__40->SetBinContent(22,0.2656894);
   MVA_BDT_Train_S__40->SetBinContent(23,0.4084479);
   MVA_BDT_Train_S__40->SetBinContent(24,0.7045396);
   MVA_BDT_Train_S__40->SetBinContent(25,0.9636198);
   MVA_BDT_Train_S__40->SetBinContent(26,1.364137);
   MVA_BDT_Train_S__40->SetBinContent(27,1.915343);
   MVA_BDT_Train_S__40->SetBinContent(28,2.483733);
   MVA_BDT_Train_S__40->SetBinContent(29,3.420916);
   MVA_BDT_Train_S__40->SetBinContent(30,4.211375);
   MVA_BDT_Train_S__40->SetBinContent(31,4.888156);
   MVA_BDT_Train_S__40->SetBinContent(32,5.346834);
   MVA_BDT_Train_S__40->SetBinContent(33,5.215972);
   MVA_BDT_Train_S__40->SetBinContent(34,4.677984);
   MVA_BDT_Train_S__40->SetBinContent(35,3.908674);
   MVA_BDT_Train_S__40->SetBinContent(36,2.956951);
   MVA_BDT_Train_S__40->SetBinContent(37,1.701205);
   MVA_BDT_Train_S__40->SetBinContent(38,0.8406889);
   MVA_BDT_Train_S__40->SetBinContent(39,0.2458618);
   MVA_BDT_Train_S__40->SetBinContent(40,0.02643676);
   MVA_BDT_Train_S__40->SetBinContent(41,0.002643676);
   MVA_BDT_Train_S__40->SetBinError(12,0.001321838);
   MVA_BDT_Train_S__40->SetBinError(13,0.001321838);
   MVA_BDT_Train_S__40->SetBinError(15,0.00228949);
   MVA_BDT_Train_S__40->SetBinError(16,0.002955719);
   MVA_BDT_Train_S__40->SetBinError(17,0.002643676);
   MVA_BDT_Train_S__40->SetBinError(18,0.003497254);
   MVA_BDT_Train_S__40->SetBinError(19,0.005761758);
   MVA_BDT_Train_S__40->SetBinError(20,0.008148356);
   MVA_BDT_Train_S__40->SetBinError(21,0.01098001);
   MVA_BDT_Train_S__40->SetBinError(22,0.01874029);
   MVA_BDT_Train_S__40->SetBinError(23,0.02323579);
   MVA_BDT_Train_S__40->SetBinError(24,0.030517);
   MVA_BDT_Train_S__40->SetBinError(25,0.03568962);
   MVA_BDT_Train_S__40->SetBinError(26,0.04246372);
   MVA_BDT_Train_S__40->SetBinError(27,0.05031673);
   MVA_BDT_Train_S__40->SetBinError(28,0.05729828);
   MVA_BDT_Train_S__40->SetBinError(29,0.06724505);
   MVA_BDT_Train_S__40->SetBinError(30,0.07461069);
   MVA_BDT_Train_S__40->SetBinError(31,0.08038252);
   MVA_BDT_Train_S__40->SetBinError(32,0.0840693);
   MVA_BDT_Train_S__40->SetBinError(33,0.08303415);
   MVA_BDT_Train_S__40->SetBinError(34,0.07863546);
   MVA_BDT_Train_S__40->SetBinError(35,0.0718793);
   MVA_BDT_Train_S__40->SetBinError(36,0.06251888);
   MVA_BDT_Train_S__40->SetBinError(37,0.04742064);
   MVA_BDT_Train_S__40->SetBinError(38,0.03333548);
   MVA_BDT_Train_S__40->SetBinError(39,0.01802746);
   MVA_BDT_Train_S__40->SetBinError(40,0.005911439);
   MVA_BDT_Train_S__40->SetBinError(41,0.001869361);
   MVA_BDT_Train_S__40->SetEntries(34606);

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
   
   TH1D *MVA_BDT_Train_B__41 = new TH1D("MVA_BDT_Train_B__41","MVA_BDT_Train_B",40,-0.5119888,0.3625026);
   MVA_BDT_Train_B__41->SetBinContent(0,0.002820012);
   MVA_BDT_Train_B__41->SetBinContent(1,0.009398751);
   MVA_BDT_Train_B__41->SetBinContent(2,0.08855577);
   MVA_BDT_Train_B__41->SetBinContent(3,0.1439029);
   MVA_BDT_Train_B__41->SetBinContent(4,0.2675447);
   MVA_BDT_Train_B__41->SetBinContent(5,0.5047205);
   MVA_BDT_Train_B__41->SetBinContent(6,1.013496);
   MVA_BDT_Train_B__41->SetBinContent(7,1.749193);
   MVA_BDT_Train_B__41->SetBinContent(8,2.457943);
   MVA_BDT_Train_B__41->SetBinContent(9,2.662315);
   MVA_BDT_Train_B__41->SetBinContent(10,2.743908);
   MVA_BDT_Train_B__41->SetBinContent(11,2.88272);
   MVA_BDT_Train_B__41->SetBinContent(12,2.887491);
   MVA_BDT_Train_B__41->SetBinContent(13,3.160359);
   MVA_BDT_Train_B__41->SetBinContent(14,3.256145);
   MVA_BDT_Train_B__41->SetBinContent(15,3.343357);
   MVA_BDT_Train_B__41->SetBinContent(16,3.191425);
   MVA_BDT_Train_B__41->SetBinContent(17,2.991313);
   MVA_BDT_Train_B__41->SetBinContent(18,2.6313);
   MVA_BDT_Train_B__41->SetBinContent(19,2.345323);
   MVA_BDT_Train_B__41->SetBinContent(20,1.968734);
   MVA_BDT_Train_B__41->SetBinContent(21,1.651257);
   MVA_BDT_Train_B__41->SetBinContent(22,1.331959);
   MVA_BDT_Train_B__41->SetBinContent(23,0.95082);
   MVA_BDT_Train_B__41->SetBinContent(24,0.6215906);
   MVA_BDT_Train_B__41->SetBinContent(25,0.4106033);
   MVA_BDT_Train_B__41->SetBinContent(26,0.2061968);
   MVA_BDT_Train_B__41->SetBinContent(27,0.1126824);
   MVA_BDT_Train_B__41->SetBinContent(28,0.06504652);
   MVA_BDT_Train_B__41->SetBinContent(29,0.03640554);
   MVA_BDT_Train_B__41->SetBinContent(30,0.01383413);
   MVA_BDT_Train_B__41->SetBinContent(31,0.02512801);
   MVA_BDT_Train_B__41->SetBinContent(32,0.006546397);
   MVA_BDT_Train_B__41->SetBinContent(33,0.002494186);
   MVA_BDT_Train_B__41->SetBinContent(34,0.005610237);
   MVA_BDT_Train_B__41->SetBinContent(35,0.001558026);
   MVA_BDT_Train_B__41->SetBinError(0,0.00179548);
   MVA_BDT_Train_B__41->SetBinError(1,0.003669876);
   MVA_BDT_Train_B__41->SetBinError(2,0.01299877);
   MVA_BDT_Train_B__41->SetBinError(3,0.01572252);
   MVA_BDT_Train_B__41->SetBinError(4,0.02239525);
   MVA_BDT_Train_B__41->SetBinError(5,0.02904196);
   MVA_BDT_Train_B__41->SetBinError(6,0.03761474);
   MVA_BDT_Train_B__41->SetBinError(7,0.04968733);
   MVA_BDT_Train_B__41->SetBinError(8,0.0601662);
   MVA_BDT_Train_B__41->SetBinError(9,0.06357086);
   MVA_BDT_Train_B__41->SetBinError(10,0.06560148);
   MVA_BDT_Train_B__41->SetBinError(11,0.06983408);
   MVA_BDT_Train_B__41->SetBinError(12,0.06992938);
   MVA_BDT_Train_B__41->SetBinError(13,0.07461871);
   MVA_BDT_Train_B__41->SetBinError(14,0.07640927);
   MVA_BDT_Train_B__41->SetBinError(15,0.0780074);
   MVA_BDT_Train_B__41->SetBinError(16,0.07708573);
   MVA_BDT_Train_B__41->SetBinError(17,0.07472688);
   MVA_BDT_Train_B__41->SetBinError(18,0.07103723);
   MVA_BDT_Train_B__41->SetBinError(19,0.06727815);
   MVA_BDT_Train_B__41->SetBinError(20,0.06150171);
   MVA_BDT_Train_B__41->SetBinError(21,0.05711703);
   MVA_BDT_Train_B__41->SetBinError(22,0.05135922);
   MVA_BDT_Train_B__41->SetBinError(23,0.04321978);
   MVA_BDT_Train_B__41->SetBinError(24,0.03492544);
   MVA_BDT_Train_B__41->SetBinError(25,0.02838996);
   MVA_BDT_Train_B__41->SetBinError(26,0.01944262);
   MVA_BDT_Train_B__41->SetBinError(27,0.01377365);
   MVA_BDT_Train_B__41->SetBinError(28,0.01144152);
   MVA_BDT_Train_B__41->SetBinError(29,0.008648187);
   MVA_BDT_Train_B__41->SetBinError(30,0.00488519);
   MVA_BDT_Train_B__41->SetBinError(31,0.006883098);
   MVA_BDT_Train_B__41->SetBinError(32,0.003856082);
   MVA_BDT_Train_B__41->SetBinError(33,0.002494186);
   MVA_BDT_Train_B__41->SetBinError(34,0.00332804);
   MVA_BDT_Train_B__41->SetBinError(35,0.001558026);
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
   TText *text = new TText(0.12,0.74,"Kolmogorov-Smirnov test: signal (background) probability = 0.255 (0.652)");
   text->SetNDC();
   text->SetTextSize(0.032);
   text->Draw();
   
   TH2F *frameBDT__42 = new TH2F("frameBDT__42","TMVA overtraining check for classifier: BDT",500,-0.5119888,0.3625026,500,0,7.031242);
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
