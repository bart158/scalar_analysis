#ifdef __CLING__
#pragma cling optimize(0)
#endif
void train_bdt_qqll_res()
{
//=========Macro generated from canvas: canvas1/TMVA comparison BDT
//=========  (Wed Aug 14 11:48:10 2024) by ROOT version 6.32.02
   TCanvas *canvas1 = new TCanvas("canvas1", "TMVA comparison BDT",200,81,600,468);
   gStyle->SetOptStat(0);
   canvas1->Range(-0.6227902,-1.288735,0.3838944,9.450727);
   canvas1->SetFillColor(0);
   canvas1->SetBorderMode(0);
   canvas1->SetBorderSize(2);
   canvas1->SetTickx(1);
   canvas1->SetTicky(1);
   canvas1->SetRightMargin(0.05);
   canvas1->SetBottomMargin(0.12);
   canvas1->SetFrameBorderMode(0);
   canvas1->SetFrameBorderMode(0);
   
   TH2F *frameBDT__37 = new TH2F("frameBDT__37","TMVA overtraining check for classifier: BDT",500,-0.5170884,0.3335602,500,0,8.376781);
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
   
   TH1D *MVA_BDT_S__38 = new TH1D("MVA_BDT_S__38","TMVA overtraining check for classifier: BDT",40,-0.5170884,0.3335602);
   MVA_BDT_S__38->SetBinContent(14,0.0001392563);
   MVA_BDT_S__38->SetBinContent(17,0.00782401);
   MVA_BDT_S__38->SetBinContent(18,0.003912005);
   MVA_BDT_S__38->SetBinContent(19,0.00782401);
   MVA_BDT_S__38->SetBinContent(20,0.02221872);
   MVA_BDT_S__38->SetBinContent(21,0.07434093);
   MVA_BDT_S__38->SetBinContent(22,0.09571129);
   MVA_BDT_S__38->SetBinContent(23,0.2112739);
   MVA_BDT_S__38->SetBinContent(24,0.3839849);
   MVA_BDT_S__38->SetBinContent(25,0.50765);
   MVA_BDT_S__38->SetBinContent(26,0.7721474);
   MVA_BDT_S__38->SetBinContent(27,1.304814);
   MVA_BDT_S__38->SetBinContent(28,1.90858);
   MVA_BDT_S__38->SetBinContent(29,2.784706);
   MVA_BDT_S__38->SetBinContent(30,3.838582);
   MVA_BDT_S__38->SetBinContent(31,5.363456);
   MVA_BDT_S__38->SetBinContent(32,6.304124);
   MVA_BDT_S__38->SetBinContent(33,6.443677);
   MVA_BDT_S__38->SetBinContent(34,6.016851);
   MVA_BDT_S__38->SetBinContent(35,5.115151);
   MVA_BDT_S__38->SetBinContent(36,3.377851);
   MVA_BDT_S__38->SetBinContent(37,1.708777);
   MVA_BDT_S__38->SetBinContent(38,0.6280867);
   MVA_BDT_S__38->SetBinContent(39,0.1328818);
   MVA_BDT_S__38->SetBinContent(40,0.008381036);
   MVA_BDT_S__38->SetBinError(14,0.0001392563);
   MVA_BDT_S__38->SetBinError(17,0.005146077);
   MVA_BDT_S__38->SetBinError(18,0.003638826);
   MVA_BDT_S__38->SetBinError(19,0.005146077);
   MVA_BDT_S__38->SetBinError(20,0.008903471);
   MVA_BDT_S__38->SetBinError(21,0.01625663);
   MVA_BDT_S__38->SetBinError(22,0.01818613);
   MVA_BDT_S__38->SetBinError(23,0.02721053);
   MVA_BDT_S__38->SetBinError(24,0.03672183);
   MVA_BDT_S__38->SetBinError(25,0.04224561);
   MVA_BDT_S__38->SetBinError(26,0.05193812);
   MVA_BDT_S__38->SetBinError(27,0.06763592);
   MVA_BDT_S__38->SetBinError(28,0.08171539);
   MVA_BDT_S__38->SetBinError(29,0.09878079);
   MVA_BDT_S__38->SetBinError(30,0.1160692);
   MVA_BDT_S__38->SetBinError(31,0.1371168);
   MVA_BDT_S__38->SetBinError(32,0.1486384);
   MVA_BDT_S__38->SetBinError(33,0.1501482);
   MVA_BDT_S__38->SetBinError(34,0.1452624);
   MVA_BDT_S__38->SetBinError(35,0.1338973);
   MVA_BDT_S__38->SetBinError(36,0.1087238);
   MVA_BDT_S__38->SetBinError(37,0.07739104);
   MVA_BDT_S__38->SetBinError(38,0.04685134);
   MVA_BDT_S__38->SetBinError(39,0.02151452);
   MVA_BDT_S__38->SetBinError(40,0.005153608);
   MVA_BDT_S__38->SetEntries(24999);
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
   
   TH1D *MVA_BDT_B__39 = new TH1D("MVA_BDT_B__39","MVA_BDT_B",40,-0.5170884,0.3335602);
   MVA_BDT_B__39->SetBinContent(1,0.01010482);
   MVA_BDT_B__39->SetBinContent(2,0.03700404);
   MVA_BDT_B__39->SetBinContent(3,0.1649723);
   MVA_BDT_B__39->SetBinContent(4,0.4843541);
   MVA_BDT_B__39->SetBinContent(5,1.008115);
   MVA_BDT_B__39->SetBinContent(6,1.772296);
   MVA_BDT_B__39->SetBinContent(7,2.57467);
   MVA_BDT_B__39->SetBinContent(8,3.709718);
   MVA_BDT_B__39->SetBinContent(9,4.220522);
   MVA_BDT_B__39->SetBinContent(10,4.840048);
   MVA_BDT_B__39->SetBinContent(11,4.577473);
   MVA_BDT_B__39->SetBinContent(12,4.609056);
   MVA_BDT_B__39->SetBinContent(13,4.055317);
   MVA_BDT_B__39->SetBinContent(14,3.160951);
   MVA_BDT_B__39->SetBinContent(15,2.349441);
   MVA_BDT_B__39->SetBinContent(16,1.842421);
   MVA_BDT_B__39->SetBinContent(17,1.299241);
   MVA_BDT_B__39->SetBinContent(18,1.071436);
   MVA_BDT_B__39->SetBinContent(19,0.8955494);
   MVA_BDT_B__39->SetBinContent(20,0.6887886);
   MVA_BDT_B__39->SetBinContent(21,0.6909352);
   MVA_BDT_B__39->SetBinContent(22,0.6297974);
   MVA_BDT_B__39->SetBinContent(23,0.4951676);
   MVA_BDT_B__39->SetBinContent(24,0.4602613);
   MVA_BDT_B__39->SetBinContent(25,0.3835994);
   MVA_BDT_B__39->SetBinContent(26,0.2959287);
   MVA_BDT_B__39->SetBinContent(27,0.2083703);
   MVA_BDT_B__39->SetBinContent(28,0.1800549);
   MVA_BDT_B__39->SetBinContent(29,0.1213501);
   MVA_BDT_B__39->SetBinContent(30,0.05145433);
   MVA_BDT_B__39->SetBinContent(31,0.03137717);
   MVA_BDT_B__39->SetBinContent(32,0.03446994);
   MVA_BDT_B__39->SetBinContent(33,0.03001909);
   MVA_BDT_B__39->SetBinContent(34,0.01745497);
   MVA_BDT_B__39->SetBinContent(35,0.01273043);
   MVA_BDT_B__39->SetBinContent(36,0.004191263);
   MVA_BDT_B__39->SetBinContent(38,0.004309077);
   MVA_BDT_B__39->SetBinError(1,0.006017488);
   MVA_BDT_B__39->SetBinError(2,0.01121346);
   MVA_BDT_B__39->SetBinError(3,0.02373766);
   MVA_BDT_B__39->SetBinError(4,0.04103031);
   MVA_BDT_B__39->SetBinError(5,0.05738536);
   MVA_BDT_B__39->SetBinError(6,0.07466877);
   MVA_BDT_B__39->SetBinError(7,0.09129259);
   MVA_BDT_B__39->SetBinError(8,0.1113548);
   MVA_BDT_B__39->SetBinError(9,0.1177831);
   MVA_BDT_B__39->SetBinError(10,0.1282312);
   MVA_BDT_B__39->SetBinError(11,0.1255546);
   MVA_BDT_B__39->SetBinError(12,0.128533);
   MVA_BDT_B__39->SetBinError(13,0.1210309);
   MVA_BDT_B__39->SetBinError(14,0.1070741);
   MVA_BDT_B__39->SetBinError(15,0.09218727);
   MVA_BDT_B__39->SetBinError(16,0.08178134);
   MVA_BDT_B__39->SetBinError(17,0.06761328);
   MVA_BDT_B__39->SetBinError(18,0.06177375);
   MVA_BDT_B__39->SetBinError(19,0.05593575);
   MVA_BDT_B__39->SetBinError(20,0.04881138);
   MVA_BDT_B__39->SetBinError(21,0.04954128);
   MVA_BDT_B__39->SetBinError(22,0.04721721);
   MVA_BDT_B__39->SetBinError(23,0.04192769);
   MVA_BDT_B__39->SetBinError(24,0.04088181);
   MVA_BDT_B__39->SetBinError(25,0.03790933);
   MVA_BDT_B__39->SetBinError(26,0.03378893);
   MVA_BDT_B__39->SetBinError(27,0.02803232);
   MVA_BDT_B__39->SetBinError(28,0.02638428);
   MVA_BDT_B__39->SetBinError(29,0.02185955);
   MVA_BDT_B__39->SetBinError(30,0.0139809);
   MVA_BDT_B__39->SetBinError(31,0.01112846);
   MVA_BDT_B__39->SetBinError(32,0.01186089);
   MVA_BDT_B__39->SetBinError(33,0.01109732);
   MVA_BDT_B__39->SetBinError(34,0.008402852);
   MVA_BDT_B__39->SetBinError(35,0.00726117);
   MVA_BDT_B__39->SetBinError(36,0.004191263);
   MVA_BDT_B__39->SetBinError(38,0.004192918);
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
   
   TH1D *MVA_BDT_Train_S__40 = new TH1D("MVA_BDT_Train_S__40","MVA_BDT_Train_S",40,-0.5170884,0.3335602);
   MVA_BDT_Train_S__40->SetBinContent(12,0.0001395475);
   MVA_BDT_Train_S__40->SetBinContent(14,0.00364109);
   MVA_BDT_Train_S__40->SetBinContent(16,0.0001395475);
   MVA_BDT_Train_S__40->SetBinContent(17,0.00364109);
   MVA_BDT_Train_S__40->SetBinContent(18,0.0001395475);
   MVA_BDT_Train_S__40->SetBinContent(19,0.004059732);
   MVA_BDT_Train_S__40->SetBinContent(20,0.02296292);
   MVA_BDT_Train_S__40->SetBinContent(21,0.03374664);
   MVA_BDT_Train_S__40->SetBinContent(22,0.03583985);
   MVA_BDT_Train_S__40->SetBinContent(23,0.08483574);
   MVA_BDT_Train_S__40->SetBinContent(24,0.130622);
   MVA_BDT_Train_S__40->SetBinContent(25,0.2761004);
   MVA_BDT_Train_S__40->SetBinContent(26,0.8034488);
   MVA_BDT_Train_S__40->SetBinContent(27,1.318187);
   MVA_BDT_Train_S__40->SetBinContent(28,2.076839);
   MVA_BDT_Train_S__40->SetBinContent(29,2.954352);
   MVA_BDT_Train_S__40->SetBinContent(30,4.202673);
   MVA_BDT_Train_S__40->SetBinContent(31,5.092811);
   MVA_BDT_Train_S__40->SetBinContent(32,6.100096);
   MVA_BDT_Train_S__40->SetBinContent(33,6.353843);
   MVA_BDT_Train_S__40->SetBinContent(34,6.094819);
   MVA_BDT_Train_S__40->SetBinContent(35,5.196702);
   MVA_BDT_Train_S__40->SetBinContent(36,3.525102);
   MVA_BDT_Train_S__40->SetBinContent(37,1.903069);
   MVA_BDT_Train_S__40->SetBinContent(38,0.6130215);
   MVA_BDT_Train_S__40->SetBinContent(39,0.1692914);
   MVA_BDT_Train_S__40->SetBinContent(40,0.02282337);
   MVA_BDT_Train_S__40->SetBinContent(41,0.003780638);
   MVA_BDT_Train_S__40->SetBinError(12,0.0001395475);
   MVA_BDT_Train_S__40->SetBinError(14,0.00364109);
   MVA_BDT_Train_S__40->SetBinError(16,0.0001395475);
   MVA_BDT_Train_S__40->SetBinError(17,0.00364109);
   MVA_BDT_Train_S__40->SetBinError(18,0.0001395475);
   MVA_BDT_Train_S__40->SetBinError(19,0.003649104);
   MVA_BDT_Train_S__40->SetBinError(20,0.008927542);
   MVA_BDT_Train_S__40->SetBinError(21,0.01092951);
   MVA_BDT_Train_S__40->SetBinError(22,0.01094286);
   MVA_BDT_Train_S__40->SetBinError(23,0.01672055);
   MVA_BDT_Train_S__40->SetBinError(24,0.02095135);
   MVA_BDT_Train_S__40->SetBinError(25,0.03072033);
   MVA_BDT_Train_S__40->SetBinError(26,0.05305656);
   MVA_BDT_Train_S__40->SetBinError(27,0.06806983);
   MVA_BDT_Train_S__40->SetBinError(28,0.08545177);
   MVA_BDT_Train_S__40->SetBinError(29,0.1018311);
   MVA_BDT_Train_S__40->SetBinError(30,0.1215065);
   MVA_BDT_Train_S__40->SetBinError(31,0.1337794);
   MVA_BDT_Train_S__40->SetBinError(32,0.1462552);
   MVA_BDT_Train_S__40->SetBinError(33,0.149351);
   MVA_BDT_Train_S__40->SetBinError(34,0.1463398);
   MVA_BDT_Train_S__40->SetBinError(35,0.1351135);
   MVA_BDT_Train_S__40->SetBinError(36,0.1112975);
   MVA_BDT_Train_S__40->SetBinError(37,0.08180026);
   MVA_BDT_Train_S__40->SetBinError(38,0.04637837);
   MVA_BDT_Train_S__40->SetBinError(39,0.02444072);
   MVA_BDT_Train_S__40->SetBinError(40,0.008926452);
   MVA_BDT_Train_S__40->SetBinError(41,0.003643763);
   MVA_BDT_Train_S__40->SetEntries(24999);

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
   
   TH1D *MVA_BDT_Train_B__41 = new TH1D("MVA_BDT_Train_B__41","MVA_BDT_Train_B",40,-0.5170884,0.3335602);
   MVA_BDT_Train_B__41->SetBinContent(1,0.01374733);
   MVA_BDT_Train_B__41->SetBinContent(2,0.0555066);
   MVA_BDT_Train_B__41->SetBinContent(3,0.1498337);
   MVA_BDT_Train_B__41->SetBinContent(4,0.5011759);
   MVA_BDT_Train_B__41->SetBinContent(5,0.9704697);
   MVA_BDT_Train_B__41->SetBinContent(6,1.795533);
   MVA_BDT_Train_B__41->SetBinContent(7,2.62447);
   MVA_BDT_Train_B__41->SetBinContent(8,3.577081);
   MVA_BDT_Train_B__41->SetBinContent(9,4.386559);
   MVA_BDT_Train_B__41->SetBinContent(10,4.947107);
   MVA_BDT_Train_B__41->SetBinContent(11,4.702975);
   MVA_BDT_Train_B__41->SetBinContent(12,4.425132);
   MVA_BDT_Train_B__41->SetBinContent(13,4.106198);
   MVA_BDT_Train_B__41->SetBinContent(14,3.054002);
   MVA_BDT_Train_B__41->SetBinContent(15,2.298567);
   MVA_BDT_Train_B__41->SetBinContent(16,1.843186);
   MVA_BDT_Train_B__41->SetBinContent(17,1.307135);
   MVA_BDT_Train_B__41->SetBinContent(18,1.032011);
   MVA_BDT_Train_B__41->SetBinContent(19,0.907488);
   MVA_BDT_Train_B__41->SetBinContent(20,0.7575961);
   MVA_BDT_Train_B__41->SetBinContent(21,0.8249943);
   MVA_BDT_Train_B__41->SetBinContent(22,0.7615692);
   MVA_BDT_Train_B__41->SetBinContent(23,0.7122702);
   MVA_BDT_Train_B__41->SetBinContent(24,0.4950793);
   MVA_BDT_Train_B__41->SetBinContent(25,0.2490534);
   MVA_BDT_Train_B__41->SetBinContent(26,0.1427845);
   MVA_BDT_Train_B__41->SetBinContent(27,0.1023099);
   MVA_BDT_Train_B__41->SetBinContent(28,0.09717781);
   MVA_BDT_Train_B__41->SetBinContent(29,0.05542196);
   MVA_BDT_Train_B__41->SetBinContent(30,0.0569615);
   MVA_BDT_Train_B__41->SetBinContent(31,0.01825374);
   MVA_BDT_Train_B__41->SetBinContent(32,0.008480889);
   MVA_BDT_Train_B__41->SetBinContent(33,0.02528709);
   MVA_BDT_Train_B__41->SetBinContent(34,0.005041999);
   MVA_BDT_Train_B__41->SetBinContent(35,0.008325317);
   MVA_BDT_Train_B__41->SetBinContent(36,0.004162659);
   MVA_BDT_Train_B__41->SetBinError(1,0.007261968);
   MVA_BDT_Train_B__41->SetBinError(2,0.01402526);
   MVA_BDT_Train_B__41->SetBinError(3,0.02241609);
   MVA_BDT_Train_B__41->SetBinError(4,0.04143951);
   MVA_BDT_Train_B__41->SetBinError(5,0.05616077);
   MVA_BDT_Train_B__41->SetBinError(6,0.07578992);
   MVA_BDT_Train_B__41->SetBinError(7,0.09205802);
   MVA_BDT_Train_B__41->SetBinError(8,0.10871);
   MVA_BDT_Train_B__41->SetBinError(9,0.1204293);
   MVA_BDT_Train_B__41->SetBinError(10,0.1299773);
   MVA_BDT_Train_B__41->SetBinError(11,0.1275124);
   MVA_BDT_Train_B__41->SetBinError(12,0.125324);
   MVA_BDT_Train_B__41->SetBinError(13,0.1218066);
   MVA_BDT_Train_B__41->SetBinError(14,0.1046739);
   MVA_BDT_Train_B__41->SetBinError(15,0.09084182);
   MVA_BDT_Train_B__41->SetBinError(16,0.08127576);
   MVA_BDT_Train_B__41->SetBinError(17,0.06806992);
   MVA_BDT_Train_B__41->SetBinError(18,0.05975332);
   MVA_BDT_Train_B__41->SetBinError(19,0.05626968);
   MVA_BDT_Train_B__41->SetBinError(20,0.05052824);
   MVA_BDT_Train_B__41->SetBinError(21,0.05466054);
   MVA_BDT_Train_B__41->SetBinError(22,0.05239287);
   MVA_BDT_Train_B__41->SetBinError(23,0.05112843);
   MVA_BDT_Train_B__41->SetBinError(24,0.04259703);
   MVA_BDT_Train_B__41->SetBinError(25,0.02987057);
   MVA_BDT_Train_B__41->SetBinError(26,0.02232991);
   MVA_BDT_Train_B__41->SetBinError(27,0.01923217);
   MVA_BDT_Train_B__41->SetBinError(28,0.0191587);
   MVA_BDT_Train_B__41->SetBinError(29,0.01451815);
   MVA_BDT_Train_B__41->SetBinError(30,0.01506547);
   MVA_BDT_Train_B__41->SetBinError(31,0.008368346);
   MVA_BDT_Train_B__41->SetBinError(32,0.005888944);
   MVA_BDT_Train_B__41->SetBinError(33,0.01019876);
   MVA_BDT_Train_B__41->SetBinError(34,0.004207016);
   MVA_BDT_Train_B__41->SetBinError(35,0.005886888);
   MVA_BDT_Train_B__41->SetBinError(36,0.004162659);
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
   TText *text = new TText(0.12,0.74,"Kolmogorov-Smirnov test: signal (background) probability = 0.028 (0.099)");
   text->SetNDC();
   text->SetTextSize(0.032);
   text->Draw();
   
   TH2F *frameBDT__42 = new TH2F("frameBDT__42","TMVA overtraining check for classifier: BDT",500,-0.5170884,0.3335602,500,0,8.376781);
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
