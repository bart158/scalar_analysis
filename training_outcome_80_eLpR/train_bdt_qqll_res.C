#ifdef __CLING__
#pragma cling optimize(0)
#endif
void train_bdt_qqll_res()
{
//=========Macro generated from canvas: canvas1/TMVA comparison BDT
//=========  (Wed Aug 14 11:40:27 2024) by ROOT version 6.32.02
   TCanvas *canvas1 = new TCanvas("canvas1", "TMVA comparison BDT",200,81,600,468);
   gStyle->SetOptStat(0);
   canvas1->Range(-0.5745338,-1.096854,0.3871879,8.043597);
   canvas1->SetFillColor(0);
   canvas1->SetBorderMode(0);
   canvas1->SetBorderSize(2);
   canvas1->SetTickx(1);
   canvas1->SetTicky(1);
   canvas1->SetRightMargin(0.05);
   canvas1->SetBottomMargin(0.12);
   canvas1->SetFrameBorderMode(0);
   canvas1->SetFrameBorderMode(0);
   
   TH2F *frameBDT__37 = new TH2F("frameBDT__37","TMVA overtraining check for classifier: BDT",500,-0.473553,0.3391018,500,0,7.129551);
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
   
   TH1D *MVA_BDT_S__38 = new TH1D("MVA_BDT_S__38","TMVA overtraining check for classifier: BDT",40,-0.473553,0.3391018);
   MVA_BDT_S__38->SetBinContent(13,0.005878114);
   MVA_BDT_S__38->SetBinContent(14,0.003104971);
   MVA_BDT_S__38->SetBinContent(15,0.008983085);
   MVA_BDT_S__38->SetBinContent(16,0.02084992);
   MVA_BDT_S__38->SetBinContent(17,0.04103619);
   MVA_BDT_S__38->SetBinContent(18,0.08273604);
   MVA_BDT_S__38->SetBinContent(19,0.1095747);
   MVA_BDT_S__38->SetBinContent(20,0.1591436);
   MVA_BDT_S__38->SetBinContent(21,0.2554056);
   MVA_BDT_S__38->SetBinContent(22,0.5327752);
   MVA_BDT_S__38->SetBinContent(23,0.8140239);
   MVA_BDT_S__38->SetBinContent(24,1.206396);
   MVA_BDT_S__38->SetBinContent(25,1.658868);
   MVA_BDT_S__38->SetBinContent(26,2.120553);
   MVA_BDT_S__38->SetBinContent(27,3.014651);
   MVA_BDT_S__38->SetBinContent(28,3.993491);
   MVA_BDT_S__38->SetBinContent(29,4.871581);
   MVA_BDT_S__38->SetBinContent(30,5.192855);
   MVA_BDT_S__38->SetBinContent(31,5.089625);
   MVA_BDT_S__38->SetBinContent(32,4.736646);
   MVA_BDT_S__38->SetBinContent(33,4.119894);
   MVA_BDT_S__38->SetBinContent(34,3.774444);
   MVA_BDT_S__38->SetBinContent(35,3.20121);
   MVA_BDT_S__38->SetBinContent(36,2.236907);
   MVA_BDT_S__38->SetBinContent(37,1.280497);
   MVA_BDT_S__38->SetBinContent(38,0.5345528);
   MVA_BDT_S__38->SetBinContent(39,0.140735);
   MVA_BDT_S__38->SetBinContent(40,0.01497181);
   MVA_BDT_S__38->SetBinError(13,0.004079742);
   MVA_BDT_S__38->SetBinError(14,0.002887992);
   MVA_BDT_S__38->SetBinError(15,0.004998479);
   MVA_BDT_S__38->SetBinError(16,0.007634502);
   MVA_BDT_S__38->SetBinError(17,0.01079342);
   MVA_BDT_S__38->SetBinError(18,0.0152666);
   MVA_BDT_S__38->SetBinError(19,0.01755025);
   MVA_BDT_S__38->SetBinError(20,0.02101064);
   MVA_BDT_S__38->SetBinError(21,0.02660828);
   MVA_BDT_S__38->SetBinError(22,0.03850203);
   MVA_BDT_S__38->SetBinError(23,0.04759453);
   MVA_BDT_S__38->SetBinError(24,0.05793318);
   MVA_BDT_S__38->SetBinError(25,0.06786649);
   MVA_BDT_S__38->SetBinError(26,0.07673639);
   MVA_BDT_S__38->SetBinError(27,0.0915341);
   MVA_BDT_S__38->SetBinError(28,0.1054409);
   MVA_BDT_S__38->SetBinError(29,0.1163016);
   MVA_BDT_S__38->SetBinError(30,0.1200387);
   MVA_BDT_S__38->SetBinError(31,0.11892);
   MVA_BDT_S__38->SetBinError(32,0.1148184);
   MVA_BDT_S__38->SetBinError(33,0.106975);
   MVA_BDT_S__38->SetBinError(34,0.1024392);
   MVA_BDT_S__38->SetBinError(35,0.09444011);
   MVA_BDT_S__38->SetBinError(36,0.07892543);
   MVA_BDT_S__38->SetBinError(37,0.0598362);
   MVA_BDT_S__38->SetBinError(38,0.03860829);
   MVA_BDT_S__38->SetBinError(39,0.01978455);
   MVA_BDT_S__38->SetBinError(40,0.006453009);
   MVA_BDT_S__38->SetEntries(33002);
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
   
   TH1D *MVA_BDT_B__39 = new TH1D("MVA_BDT_B__39","MVA_BDT_B",40,-0.473553,0.3391018);
   MVA_BDT_B__39->SetBinContent(1,0.01760837);
   MVA_BDT_B__39->SetBinContent(2,0.05607249);
   MVA_BDT_B__39->SetBinContent(3,0.1785294);
   MVA_BDT_B__39->SetBinContent(4,0.4855192);
   MVA_BDT_B__39->SetBinContent(5,0.9252818);
   MVA_BDT_B__39->SetBinContent(6,1.226872);
   MVA_BDT_B__39->SetBinContent(7,1.702692);
   MVA_BDT_B__39->SetBinContent(8,2.121974);
   MVA_BDT_B__39->SetBinContent(9,2.189736);
   MVA_BDT_B__39->SetBinContent(10,2.700824);
   MVA_BDT_B__39->SetBinContent(11,2.831723);
   MVA_BDT_B__39->SetBinContent(12,3.189099);
   MVA_BDT_B__39->SetBinContent(13,3.325815);
   MVA_BDT_B__39->SetBinContent(14,3.432868);
   MVA_BDT_B__39->SetBinContent(15,3.26539);
   MVA_BDT_B__39->SetBinContent(16,3.145772);
   MVA_BDT_B__39->SetBinContent(17,2.89221);
   MVA_BDT_B__39->SetBinContent(18,2.636088);
   MVA_BDT_B__39->SetBinContent(19,2.388822);
   MVA_BDT_B__39->SetBinContent(20,2.313053);
   MVA_BDT_B__39->SetBinContent(21,1.923521);
   MVA_BDT_B__39->SetBinContent(22,1.491684);
   MVA_BDT_B__39->SetBinContent(23,1.402287);
   MVA_BDT_B__39->SetBinContent(24,0.9901297);
   MVA_BDT_B__39->SetBinContent(25,0.7941233);
   MVA_BDT_B__39->SetBinContent(26,0.5396401);
   MVA_BDT_B__39->SetBinContent(27,0.3866471);
   MVA_BDT_B__39->SetBinContent(28,0.2151085);
   MVA_BDT_B__39->SetBinContent(29,0.17494);
   MVA_BDT_B__39->SetBinContent(30,0.1016175);
   MVA_BDT_B__39->SetBinContent(31,0.09686467);
   MVA_BDT_B__39->SetBinContent(32,0.04532609);
   MVA_BDT_B__39->SetBinContent(33,0.02356306);
   MVA_BDT_B__39->SetBinContent(34,0.004745779);
   MVA_BDT_B__39->SetBinContent(35,0.004910808);
   MVA_BDT_B__39->SetBinContent(36,0.0001650288);
   MVA_BDT_B__39->SetBinContent(37,0.0001650288);
   MVA_BDT_B__39->SetBinError(1,0.007826573);
   MVA_BDT_B__39->SetBinError(2,0.01226073);
   MVA_BDT_B__39->SetBinError(3,0.0210893);
   MVA_BDT_B__39->SetBinError(4,0.03454864);
   MVA_BDT_B__39->SetBinError(5,0.04899443);
   MVA_BDT_B__39->SetBinError(6,0.05717489);
   MVA_BDT_B__39->SetBinError(7,0.07172306);
   MVA_BDT_B__39->SetBinError(8,0.08344942);
   MVA_BDT_B__39->SetBinError(9,0.0849195);
   MVA_BDT_B__39->SetBinError(10,0.09676645);
   MVA_BDT_B__39->SetBinError(11,0.09909659);
   MVA_BDT_B__39->SetBinError(12,0.106908);
   MVA_BDT_B__39->SetBinError(13,0.1097732);
   MVA_BDT_B__39->SetBinError(14,0.1125577);
   MVA_BDT_B__39->SetBinError(15,0.1105803);
   MVA_BDT_B__39->SetBinError(16,0.1093768);
   MVA_BDT_B__39->SetBinError(17,0.1061463);
   MVA_BDT_B__39->SetBinError(18,0.1011884);
   MVA_BDT_B__39->SetBinError(19,0.09700973);
   MVA_BDT_B__39->SetBinError(20,0.09594551);
   MVA_BDT_B__39->SetBinError(21,0.0877669);
   MVA_BDT_B__39->SetBinError(22,0.07730016);
   MVA_BDT_B__39->SetBinError(23,0.07586063);
   MVA_BDT_B__39->SetBinError(24,0.06337957);
   MVA_BDT_B__39->SetBinError(25,0.05718834);
   MVA_BDT_B__39->SetBinError(26,0.04732516);
   MVA_BDT_B__39->SetBinError(27,0.03951303);
   MVA_BDT_B__39->SetBinError(28,0.02948832);
   MVA_BDT_B__39->SetBinError(29,0.02668868);
   MVA_BDT_B__39->SetBinError(30,0.0203317);
   MVA_BDT_B__39->SetBinError(31,0.02026463);
   MVA_BDT_B__39->SetBinError(32,0.01397057);
   MVA_BDT_B__39->SetBinError(33,0.009903911);
   MVA_BDT_B__39->SetBinError(34,0.004421885);
   MVA_BDT_B__39->SetBinError(35,0.004424963);
   MVA_BDT_B__39->SetBinError(36,0.0001650288);
   MVA_BDT_B__39->SetBinError(37,0.0001650288);
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
   
   TH1D *MVA_BDT_Train_S__40 = new TH1D("MVA_BDT_Train_S__40","MVA_BDT_Train_S",40,-0.473553,0.3391018);
   MVA_BDT_Train_S__40->SetBinContent(10,0.0001103106);
   MVA_BDT_Train_S__40->SetBinContent(12,0.0002206212);
   MVA_BDT_Train_S__40->SetBinContent(13,0.00575195);
   MVA_BDT_Train_S__40->SetBinContent(14,0.002986286);
   MVA_BDT_Train_S__40->SetBinContent(15,0.008958857);
   MVA_BDT_Train_S__40->SetBinContent(16,0.01791771);
   MVA_BDT_Train_S__40->SetBinContent(17,0.01548298);
   MVA_BDT_Train_S__40->SetBinContent(18,0.04755995);
   MVA_BDT_Train_S__40->SetBinContent(19,0.02067548);
   MVA_BDT_Train_S__40->SetBinContent(20,0.07033133);
   MVA_BDT_Train_S__40->SetBinContent(21,0.191335);
   MVA_BDT_Train_S__40->SetBinContent(22,0.4133129);
   MVA_BDT_Train_S__40->SetBinContent(23,0.7645916);
   MVA_BDT_Train_S__40->SetBinContent(24,1.213984);
   MVA_BDT_Train_S__40->SetBinContent(25,1.701734);
   MVA_BDT_Train_S__40->SetBinContent(26,2.294909);
   MVA_BDT_Train_S__40->SetBinContent(27,2.85708);
   MVA_BDT_Train_S__40->SetBinContent(28,4.169982);
   MVA_BDT_Train_S__40->SetBinContent(29,4.883034);
   MVA_BDT_Train_S__40->SetBinContent(30,5.48427);
   MVA_BDT_Train_S__40->SetBinContent(31,5.158308);
   MVA_BDT_Train_S__40->SetBinContent(32,4.550319);
   MVA_BDT_Train_S__40->SetBinContent(33,4.211214);
   MVA_BDT_Train_S__40->SetBinContent(34,3.840931);
   MVA_BDT_Train_S__40->SetBinContent(35,3.224866);
   MVA_BDT_Train_S__40->SetBinContent(36,2.268576);
   MVA_BDT_Train_S__40->SetBinContent(37,1.190976);
   MVA_BDT_Train_S__40->SetBinContent(38,0.4731569);
   MVA_BDT_Train_S__40->SetBinContent(39,0.121114);
   MVA_BDT_Train_S__40->SetBinContent(40,0.01769709);
   MVA_BDT_Train_S__40->SetBinError(10,0.0001103106);
   MVA_BDT_Train_S__40->SetBinError(12,0.0001560027);
   MVA_BDT_Train_S__40->SetBinError(13,0.004067243);
   MVA_BDT_Train_S__40->SetBinError(14,0.00287809);
   MVA_BDT_Train_S__40->SetBinError(15,0.004984998);
   MVA_BDT_Train_S__40->SetBinError(16,0.007049852);
   MVA_BDT_Train_S__40->SetBinError(17,0.00644033);
   MVA_BDT_Train_S__40->SetBinError(18,0.0115113);
   MVA_BDT_Train_S__40->SetBinError(19,0.007071394);
   MVA_BDT_Train_S__40->SetBinError(20,0.01351835);
   MVA_BDT_Train_S__40->SetBinError(21,0.02285185);
   MVA_BDT_Train_S__40->SetBinError(22,0.03369404);
   MVA_BDT_Train_S__40->SetBinError(23,0.04596312);
   MVA_BDT_Train_S__40->SetBinError(24,0.05806192);
   MVA_BDT_Train_S__40->SetBinError(25,0.06865519);
   MVA_BDT_Train_S__40->SetBinError(26,0.07981081);
   MVA_BDT_Train_S__40->SetBinError(27,0.08890031);
   MVA_BDT_Train_S__40->SetBinError(28,0.1076458);
   MVA_BDT_Train_S__40->SetBinError(29,0.1164445);
   MVA_BDT_Train_S__40->SetBinError(30,0.1233872);
   MVA_BDT_Train_S__40->SetBinError(31,0.1196058);
   MVA_BDT_Train_S__40->SetBinError(32,0.1123204);
   MVA_BDT_Train_S__40->SetBinError(33,0.1080725);
   MVA_BDT_Train_S__40->SetBinError(34,0.1032881);
   MVA_BDT_Train_S__40->SetBinError(35,0.09462553);
   MVA_BDT_Train_S__40->SetBinError(36,0.07929263);
   MVA_BDT_Train_S__40->SetBinError(37,0.05748927);
   MVA_BDT_Train_S__40->SetBinError(38,0.03617905);
   MVA_BDT_Train_S__40->SetBinError(39,0.0184248);
   MVA_BDT_Train_S__40->SetBinError(40,0.007048125);
   MVA_BDT_Train_S__40->SetEntries(33002);

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
   
   TH1D *MVA_BDT_Train_B__41 = new TH1D("MVA_BDT_Train_B__41","MVA_BDT_Train_B",40,-0.473553,0.3391018);
   MVA_BDT_Train_B__41->SetBinContent(0,0.001894695);
   MVA_BDT_Train_B__41->SetBinContent(1,0.02194833);
   MVA_BDT_Train_B__41->SetBinContent(2,0.05614989);
   MVA_BDT_Train_B__41->SetBinContent(3,0.1568975);
   MVA_BDT_Train_B__41->SetBinContent(4,0.3766474);
   MVA_BDT_Train_B__41->SetBinContent(5,0.9070312);
   MVA_BDT_Train_B__41->SetBinContent(6,1.168169);
   MVA_BDT_Train_B__41->SetBinContent(7,1.540916);
   MVA_BDT_Train_B__41->SetBinContent(8,2.003499);
   MVA_BDT_Train_B__41->SetBinContent(9,2.377952);
   MVA_BDT_Train_B__41->SetBinContent(10,2.655995);
   MVA_BDT_Train_B__41->SetBinContent(11,2.865097);
   MVA_BDT_Train_B__41->SetBinContent(12,3.180198);
   MVA_BDT_Train_B__41->SetBinContent(13,3.406238);
   MVA_BDT_Train_B__41->SetBinContent(14,3.466207);
   MVA_BDT_Train_B__41->SetBinContent(15,3.471429);
   MVA_BDT_Train_B__41->SetBinContent(16,3.228504);
   MVA_BDT_Train_B__41->SetBinContent(17,3.001974);
   MVA_BDT_Train_B__41->SetBinContent(18,2.790586);
   MVA_BDT_Train_B__41->SetBinContent(19,2.423904);
   MVA_BDT_Train_B__41->SetBinContent(20,2.14187);
   MVA_BDT_Train_B__41->SetBinContent(21,1.977315);
   MVA_BDT_Train_B__41->SetBinContent(22,1.883228);
   MVA_BDT_Train_B__41->SetBinContent(23,1.457057);
   MVA_BDT_Train_B__41->SetBinContent(24,1.015184);
   MVA_BDT_Train_B__41->SetBinContent(25,0.6656407);
   MVA_BDT_Train_B__41->SetBinContent(26,0.4324806);
   MVA_BDT_Train_B__41->SetBinContent(27,0.195086);
   MVA_BDT_Train_B__41->SetBinContent(28,0.1385338);
   MVA_BDT_Train_B__41->SetBinContent(29,0.08382234);
   MVA_BDT_Train_B__41->SetBinContent(30,0.07199682);
   MVA_BDT_Train_B__41->SetBinContent(31,0.03239657);
   MVA_BDT_Train_B__41->SetBinContent(32,0.01363629);
   MVA_BDT_Train_B__41->SetBinContent(33,0.004976933);
   MVA_BDT_Train_B__41->SetBinContent(35,0.004329681);
   MVA_BDT_Train_B__41->SetBinContent(36,0.004329681);
   MVA_BDT_Train_B__41->SetBinContent(40,0.0001618132);
   MVA_BDT_Train_B__41->SetBinError(0,0.001030846);
   MVA_BDT_Train_B__41->SetBinError(1,0.008816417);
   MVA_BDT_Train_B__41->SetBinError(2,0.01202529);
   MVA_BDT_Train_B__41->SetBinError(3,0.01790375);
   MVA_BDT_Train_B__41->SetBinError(4,0.02832101);
   MVA_BDT_Train_B__41->SetBinError(5,0.04757731);
   MVA_BDT_Train_B__41->SetBinError(6,0.05423762);
   MVA_BDT_Train_B__41->SetBinError(7,0.06662486);
   MVA_BDT_Train_B__41->SetBinError(8,0.08007864);
   MVA_BDT_Train_B__41->SetBinError(9,0.08957313);
   MVA_BDT_Train_B__41->SetBinError(10,0.09469928);
   MVA_BDT_Train_B__41->SetBinError(11,0.09971857);
   MVA_BDT_Train_B__41->SetBinError(12,0.1056612);
   MVA_BDT_Train_B__41->SetBinError(13,0.1102371);
   MVA_BDT_Train_B__41->SetBinError(14,0.1123727);
   MVA_BDT_Train_B__41->SetBinError(15,0.1135101);
   MVA_BDT_Train_B__41->SetBinError(16,0.1100397);
   MVA_BDT_Train_B__41->SetBinError(17,0.1074736);
   MVA_BDT_Train_B__41->SetBinError(18,0.1040536);
   MVA_BDT_Train_B__41->SetBinError(19,0.09739735);
   MVA_BDT_Train_B__41->SetBinError(20,0.0913188);
   MVA_BDT_Train_B__41->SetBinError(21,0.0882996);
   MVA_BDT_Train_B__41->SetBinError(22,0.08705154);
   MVA_BDT_Train_B__41->SetBinError(23,0.07682765);
   MVA_BDT_Train_B__41->SetBinError(24,0.06357315);
   MVA_BDT_Train_B__41->SetBinError(25,0.05141216);
   MVA_BDT_Train_B__41->SetBinError(26,0.04135078);
   MVA_BDT_Train_B__41->SetBinError(27,0.02700367);
   MVA_BDT_Train_B__41->SetBinError(28,0.02312142);
   MVA_BDT_Train_B__41->SetBinError(29,0.01799302);
   MVA_BDT_Train_B__41->SetBinError(30,0.01686763);
   MVA_BDT_Train_B__41->SetBinError(31,0.01148042);
   MVA_BDT_Train_B__41->SetBinError(32,0.007506206);
   MVA_BDT_Train_B__41->SetBinError(33,0.004341759);
   MVA_BDT_Train_B__41->SetBinError(35,0.004329681);
   MVA_BDT_Train_B__41->SetBinError(36,0.004329681);
   MVA_BDT_Train_B__41->SetBinError(40,0.0001618132);
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
   TText *text = new TText(0.12,0.74,"Kolmogorov-Smirnov test: signal (background) probability = 0.202 (0.051)");
   text->SetNDC();
   text->SetTextSize(0.032);
   text->Draw();
   
   TH2F *frameBDT__42 = new TH2F("frameBDT__42","TMVA overtraining check for classifier: BDT",500,-0.473553,0.3391018,500,0,7.129551);
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
