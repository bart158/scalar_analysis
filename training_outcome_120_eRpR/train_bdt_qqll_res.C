#ifdef __CLING__
#pragma cling optimize(0)
#endif
void train_bdt_qqll_res()
{
//=========Macro generated from canvas: canvas1/TMVA comparison BDT
//=========  (Mon Jul  1 17:59:41 2024) by ROOT version 6.32.02
   TCanvas *canvas1 = new TCanvas("canvas1", "TMVA comparison BDT",200,81,600,468);
   gStyle->SetOptStat(0);
   canvas1->Range(-0.5983755,-1.354657,0.397432,9.934148);
   canvas1->SetFillColor(0);
   canvas1->SetBorderMode(0);
   canvas1->SetBorderSize(2);
   canvas1->SetTickx(1);
   canvas1->SetTicky(1);
   canvas1->SetRightMargin(0.05);
   canvas1->SetBottomMargin(0.12);
   canvas1->SetFrameBorderMode(0);
   canvas1->SetFrameBorderMode(0);
   
   TH2F *frameBDT__37 = new TH2F("frameBDT__37","TMVA overtraining check for classifier: BDT",500,-0.4938157,0.3476416,500,0,8.805267);
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
   
   TH1D *MVA_BDT_S__38 = new TH1D("MVA_BDT_S__38","TMVA overtraining check for classifier: BDT",40,-0.4938157,0.3476416);
   MVA_BDT_S__38->SetBinContent(11,0.001901539);
   MVA_BDT_S__38->SetBinContent(16,0.003803079);
   MVA_BDT_S__38->SetBinContent(17,0.003803079);
   MVA_BDT_S__38->SetBinContent(18,0.01140924);
   MVA_BDT_S__38->SetBinContent(19,0.02662155);
   MVA_BDT_S__38->SetBinContent(20,0.0437354);
   MVA_BDT_S__38->SetBinContent(21,0.08747081);
   MVA_BDT_S__38->SetBinContent(22,0.1464185);
   MVA_BDT_S__38->SetBinContent(23,0.2262832);
   MVA_BDT_S__38->SetBinContent(24,0.4145356);
   MVA_BDT_S__38->SetBinContent(25,0.6465234);
   MVA_BDT_S__38->SetBinContent(26,0.965982);
   MVA_BDT_S__38->SetBinContent(27,1.61821);
   MVA_BDT_S__38->SetBinContent(28,2.241915);
   MVA_BDT_S__38->SetBinContent(29,3.295368);
   MVA_BDT_S__38->SetBinContent(30,4.622642);
   MVA_BDT_S__38->SetBinContent(31,6.012667);
   MVA_BDT_S__38->SetBinContent(32,6.39868);
   MVA_BDT_S__38->SetBinContent(33,6.773283);
   MVA_BDT_S__38->SetBinContent(34,5.81871);
   MVA_BDT_S__38->SetBinContent(35,4.352623);
   MVA_BDT_S__38->SetBinContent(36,2.375023);
   MVA_BDT_S__38->SetBinContent(37,1.081976);
   MVA_BDT_S__38->SetBinContent(38,0.3080494);
   MVA_BDT_S__38->SetBinContent(39,0.05514464);
   MVA_BDT_S__38->SetBinContent(40,0.003803079);
   MVA_BDT_S__38->SetBinError(11,0.001901539);
   MVA_BDT_S__38->SetBinError(16,0.002689183);
   MVA_BDT_S__38->SetBinError(17,0.002689183);
   MVA_BDT_S__38->SetBinError(18,0.004657801);
   MVA_BDT_S__38->SetBinError(19,0.007114908);
   MVA_BDT_S__38->SetBinError(20,0.009119462);
   MVA_BDT_S__38->SetBinError(21,0.01289687);
   MVA_BDT_S__38->SetBinError(22,0.01668594);
   MVA_BDT_S__38->SetBinError(23,0.02074334);
   MVA_BDT_S__38->SetBinError(24,0.02807589);
   MVA_BDT_S__38->SetBinError(25,0.03506265);
   MVA_BDT_S__38->SetBinError(26,0.04285852);
   MVA_BDT_S__38->SetBinError(27,0.05547152);
   MVA_BDT_S__38->SetBinError(28,0.06529234);
   MVA_BDT_S__38->SetBinError(29,0.07915978);
   MVA_BDT_S__38->SetBinError(30,0.09375572);
   MVA_BDT_S__38->SetBinError(31,0.1069267);
   MVA_BDT_S__38->SetBinError(32,0.1103057);
   MVA_BDT_S__38->SetBinError(33,0.1134886);
   MVA_BDT_S__38->SetBinError(34,0.105188);
   MVA_BDT_S__38->SetBinError(35,0.09097628);
   MVA_BDT_S__38->SetBinError(36,0.06720267);
   MVA_BDT_S__38->SetBinError(37,0.04535879);
   MVA_BDT_S__38->SetBinError(38,0.02420264);
   MVA_BDT_S__38->SetBinError(39,0.0102401);
   MVA_BDT_S__38->SetBinError(40,0.002689183);
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
   
   TH1D *MVA_BDT_B__39 = new TH1D("MVA_BDT_B__39","MVA_BDT_B",40,-0.4938157,0.3476416);
   MVA_BDT_B__39->SetBinContent(1,0.04359294);
   MVA_BDT_B__39->SetBinContent(2,0.1028906);
   MVA_BDT_B__39->SetBinContent(3,0.4860431);
   MVA_BDT_B__39->SetBinContent(4,1.146327);
   MVA_BDT_B__39->SetBinContent(5,2.318017);
   MVA_BDT_B__39->SetBinContent(6,3.402951);
   MVA_BDT_B__39->SetBinContent(7,4.398908);
   MVA_BDT_B__39->SetBinContent(8,4.748905);
   MVA_BDT_B__39->SetBinContent(9,4.832204);
   MVA_BDT_B__39->SetBinContent(10,4.33607);
   MVA_BDT_B__39->SetBinContent(11,3.932633);
   MVA_BDT_B__39->SetBinContent(12,3.528318);
   MVA_BDT_B__39->SetBinContent(13,3.004553);
   MVA_BDT_B__39->SetBinContent(14,2.528323);
   MVA_BDT_B__39->SetBinContent(15,1.883959);
   MVA_BDT_B__39->SetBinContent(16,1.430086);
   MVA_BDT_B__39->SetBinContent(17,1.032814);
   MVA_BDT_B__39->SetBinContent(18,0.8768255);
   MVA_BDT_B__39->SetBinContent(19,0.5919893);
   MVA_BDT_B__39->SetBinContent(20,0.5349211);
   MVA_BDT_B__39->SetBinContent(21,0.4739923);
   MVA_BDT_B__39->SetBinContent(22,0.4382745);
   MVA_BDT_B__39->SetBinContent(23,0.3460158);
   MVA_BDT_B__39->SetBinContent(24,0.2853909);
   MVA_BDT_B__39->SetBinContent(25,0.245466);
   MVA_BDT_B__39->SetBinContent(26,0.1938505);
   MVA_BDT_B__39->SetBinContent(27,0.139552);
   MVA_BDT_B__39->SetBinContent(28,0.08606272);
   MVA_BDT_B__39->SetBinContent(29,0.07592146);
   MVA_BDT_B__39->SetBinContent(30,0.03809256);
   MVA_BDT_B__39->SetBinContent(31,0.0193218);
   MVA_BDT_B__39->SetBinContent(32,0.01993974);
   MVA_BDT_B__39->SetBinContent(33,0.005721291);
   MVA_BDT_B__39->SetBinContent(34,0.002797014);
   MVA_BDT_B__39->SetBinContent(35,0.002797014);
   MVA_BDT_B__39->SetBinContent(36,0.000932338);
   MVA_BDT_B__39->SetBinContent(37,0.002122813);
   MVA_BDT_B__39->SetBinError(1,0.009588628);
   MVA_BDT_B__39->SetBinError(2,0.01457173);
   MVA_BDT_B__39->SetBinError(3,0.03061658);
   MVA_BDT_B__39->SetBinError(4,0.04553173);
   MVA_BDT_B__39->SetBinError(5,0.06526842);
   MVA_BDT_B__39->SetBinError(6,0.08039406);
   MVA_BDT_B__39->SetBinError(7,0.0927031);
   MVA_BDT_B__39->SetBinError(8,0.09753741);
   MVA_BDT_B__39->SetBinError(9,0.09983954);
   MVA_BDT_B__39->SetBinError(10,0.0948383);
   MVA_BDT_B__39->SetBinError(11,0.09143114);
   MVA_BDT_B__39->SetBinError(12,0.08623027);
   MVA_BDT_B__39->SetBinError(13,0.07975847);
   MVA_BDT_B__39->SetBinError(14,0.07346026);
   MVA_BDT_B__39->SetBinError(15,0.0632422);
   MVA_BDT_B__39->SetBinError(16,0.05498999);
   MVA_BDT_B__39->SetBinError(17,0.04530433);
   MVA_BDT_B__39->SetBinError(18,0.04263513);
   MVA_BDT_B__39->SetBinError(19,0.03388875);
   MVA_BDT_B__39->SetBinError(20,0.03241127);
   MVA_BDT_B__39->SetBinError(21,0.03077337);
   MVA_BDT_B__39->SetBinError(22,0.030139);
   MVA_BDT_B__39->SetBinError(23,0.02635349);
   MVA_BDT_B__39->SetBinError(24,0.02370991);
   MVA_BDT_B__39->SetBinError(25,0.02313617);
   MVA_BDT_B__39->SetBinError(26,0.01996051);
   MVA_BDT_B__39->SetBinError(27,0.01667354);
   MVA_BDT_B__39->SetBinError(28,0.0135796);
   MVA_BDT_B__39->SetBinError(29,0.01346828);
   MVA_BDT_B__39->SetBinError(30,0.008783296);
   MVA_BDT_B__39->SetBinError(31,0.00645342);
   MVA_BDT_B__39->SetBinError(32,0.006715616);
   MVA_BDT_B__39->SetBinError(33,0.002287296);
   MVA_BDT_B__39->SetBinError(34,0.001614857);
   MVA_BDT_B__39->SetBinError(35,0.001614857);
   MVA_BDT_B__39->SetBinError(36,0.000932338);
   MVA_BDT_B__39->SetBinError(37,0.002122813);
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
   
   TH1D *MVA_BDT_Train_S__40 = new TH1D("MVA_BDT_Train_S__40","MVA_BDT_Train_S",40,-0.4938157,0.3476416);
   MVA_BDT_Train_S__40->SetBinContent(12,0.001901539);
   MVA_BDT_Train_S__40->SetBinContent(14,0.001901539);
   MVA_BDT_Train_S__40->SetBinContent(17,0.003803079);
   MVA_BDT_Train_S__40->SetBinContent(18,0.001901539);
   MVA_BDT_Train_S__40->SetBinContent(19,0.009507696);
   MVA_BDT_Train_S__40->SetBinContent(20,0.03232617);
   MVA_BDT_Train_S__40->SetBinContent(21,0.03803079);
   MVA_BDT_Train_S__40->SetBinContent(22,0.05894772);
   MVA_BDT_Train_S__40->SetBinContent(23,0.1255016);
   MVA_BDT_Train_S__40->SetBinContent(24,0.3270648);
   MVA_BDT_Train_S__40->SetBinContent(25,0.6370157);
   MVA_BDT_Train_S__40->SetBinContent(26,1.042044);
   MVA_BDT_Train_S__40->SetBinContent(27,1.665748);
   MVA_BDT_Train_S__40->SetBinContent(28,2.498623);
   MVA_BDT_Train_S__40->SetBinContent(29,3.588205);
   MVA_BDT_Train_S__40->SetBinContent(30,4.615036);
   MVA_BDT_Train_S__40->SetBinContent(31,5.664685);
   MVA_BDT_Train_S__40->SetBinContent(32,6.480446);
   MVA_BDT_Train_S__40->SetBinContent(33,6.617357);
   MVA_BDT_Train_S__40->SetBinContent(34,5.698913);
   MVA_BDT_Train_S__40->SetBinContent(35,4.377343);
   MVA_BDT_Train_S__40->SetBinContent(36,2.688777);
   MVA_BDT_Train_S__40->SetBinContent(37,1.032536);
   MVA_BDT_Train_S__40->SetBinContent(38,0.2814278);
   MVA_BDT_Train_S__40->SetBinContent(39,0.04183386);
   MVA_BDT_Train_S__40->SetBinContent(40,0.005704618);
   MVA_BDT_Train_S__40->SetBinError(12,0.001901539);
   MVA_BDT_Train_S__40->SetBinError(14,0.001901539);
   MVA_BDT_Train_S__40->SetBinError(17,0.002689183);
   MVA_BDT_Train_S__40->SetBinError(18,0.001901539);
   MVA_BDT_Train_S__40->SetBinError(19,0.004251971);
   MVA_BDT_Train_S__40->SetBinError(20,0.007840247);
   MVA_BDT_Train_S__40->SetBinError(21,0.008503942);
   MVA_BDT_Train_S__40->SetBinError(22,0.01058732);
   MVA_BDT_Train_S__40->SetBinError(23,0.01544818);
   MVA_BDT_Train_S__40->SetBinError(24,0.02493845);
   MVA_BDT_Train_S__40->SetBinError(25,0.03480388);
   MVA_BDT_Train_S__40->SetBinError(26,0.04451389);
   MVA_BDT_Train_S__40->SetBinError(27,0.05628042);
   MVA_BDT_Train_S__40->SetBinError(28,0.06892916);
   MVA_BDT_Train_S__40->SetBinError(29,0.08260213);
   MVA_BDT_Train_S__40->SetBinError(30,0.09367856);
   MVA_BDT_Train_S__40->SetBinError(31,0.1037864);
   MVA_BDT_Train_S__40->SetBinError(32,0.1110082);
   MVA_BDT_Train_S__40->SetBinError(33,0.1121747);
   MVA_BDT_Train_S__40->SetBinError(34,0.1040995);
   MVA_BDT_Train_S__40->SetBinError(35,0.09123426);
   MVA_BDT_Train_S__40->SetBinError(36,0.07150395);
   MVA_BDT_Train_S__40->SetBinError(37,0.04431035);
   MVA_BDT_Train_S__40->SetBinError(38,0.02313322);
   MVA_BDT_Train_S__40->SetBinError(39,0.00891901);
   MVA_BDT_Train_S__40->SetBinError(40,0.003293563);
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
   
   TH1D *MVA_BDT_Train_B__41 = new TH1D("MVA_BDT_Train_B__41","MVA_BDT_Train_B",40,-0.4938157,0.3476416);
   MVA_BDT_Train_B__41->SetBinContent(0,0.001731086);
   MVA_BDT_Train_B__41->SetBinContent(1,0.02616344);
   MVA_BDT_Train_B__41->SetBinContent(2,0.1624936);
   MVA_BDT_Train_B__41->SetBinContent(3,0.438474);
   MVA_BDT_Train_B__41->SetBinContent(4,1.1807);
   MVA_BDT_Train_B__41->SetBinContent(5,2.146379);
   MVA_BDT_Train_B__41->SetBinContent(6,3.432173);
   MVA_BDT_Train_B__41->SetBinContent(7,4.483528);
   MVA_BDT_Train_B__41->SetBinContent(8,4.77767);
   MVA_BDT_Train_B__41->SetBinContent(9,4.662783);
   MVA_BDT_Train_B__41->SetBinContent(10,4.446892);
   MVA_BDT_Train_B__41->SetBinContent(11,4.071035);
   MVA_BDT_Train_B__41->SetBinContent(12,3.542018);
   MVA_BDT_Train_B__41->SetBinContent(13,3.089572);
   MVA_BDT_Train_B__41->SetBinContent(14,2.457647);
   MVA_BDT_Train_B__41->SetBinContent(15,1.857204);
   MVA_BDT_Train_B__41->SetBinContent(16,1.36627);
   MVA_BDT_Train_B__41->SetBinContent(17,0.9984657);
   MVA_BDT_Train_B__41->SetBinContent(18,0.8527149);
   MVA_BDT_Train_B__41->SetBinContent(19,0.6602607);
   MVA_BDT_Train_B__41->SetBinContent(20,0.6079587);
   MVA_BDT_Train_B__41->SetBinContent(21,0.5754013);
   MVA_BDT_Train_B__41->SetBinContent(22,0.4908689);
   MVA_BDT_Train_B__41->SetBinContent(23,0.3995427);
   MVA_BDT_Train_B__41->SetBinContent(24,0.3159066);
   MVA_BDT_Train_B__41->SetBinContent(25,0.1738347);
   MVA_BDT_Train_B__41->SetBinContent(26,0.1286291);
   MVA_BDT_Train_B__41->SetBinContent(27,0.07087031);
   MVA_BDT_Train_B__41->SetBinContent(28,0.05315588);
   MVA_BDT_Train_B__41->SetBinContent(29,0.02156124);
   MVA_BDT_Train_B__41->SetBinContent(30,0.02464785);
   MVA_BDT_Train_B__41->SetBinContent(31,0.003876679);
   MVA_BDT_Train_B__41->SetBinContent(32,0.009422892);
   MVA_BDT_Train_B__41->SetBinContent(33,0.002811565);
   MVA_BDT_Train_B__41->SetBinContent(34,0.004711446);
   MVA_BDT_Train_B__41->SetBinContent(37,0.0009371884);
   MVA_BDT_Train_B__41->SetBinError(0,0.001228249);
   MVA_BDT_Train_B__41->SetBinError(1,0.007063833);
   MVA_BDT_Train_B__41->SetBinError(2,0.01817133);
   MVA_BDT_Train_B__41->SetBinError(3,0.02886155);
   MVA_BDT_Train_B__41->SetBinError(4,0.04676161);
   MVA_BDT_Train_B__41->SetBinError(5,0.06220294);
   MVA_BDT_Train_B__41->SetBinError(6,0.0807796);
   MVA_BDT_Train_B__41->SetBinError(7,0.09396655);
   MVA_BDT_Train_B__41->SetBinError(8,0.09785422);
   MVA_BDT_Train_B__41->SetBinError(9,0.09729815);
   MVA_BDT_Train_B__41->SetBinError(10,0.09624241);
   MVA_BDT_Train_B__41->SetBinError(11,0.09313497);
   MVA_BDT_Train_B__41->SetBinError(12,0.0867097);
   MVA_BDT_Train_B__41->SetBinError(13,0.08102375);
   MVA_BDT_Train_B__41->SetBinError(14,0.07232736);
   MVA_BDT_Train_B__41->SetBinError(15,0.06302567);
   MVA_BDT_Train_B__41->SetBinError(16,0.0532328);
   MVA_BDT_Train_B__41->SetBinError(17,0.04496392);
   MVA_BDT_Train_B__41->SetBinError(18,0.04161512);
   MVA_BDT_Train_B__41->SetBinError(19,0.03636241);
   MVA_BDT_Train_B__41->SetBinError(20,0.03472311);
   MVA_BDT_Train_B__41->SetBinError(21,0.03425925);
   MVA_BDT_Train_B__41->SetBinError(22,0.03173469);
   MVA_BDT_Train_B__41->SetBinError(23,0.02834031);
   MVA_BDT_Train_B__41->SetBinError(24,0.02544164);
   MVA_BDT_Train_B__41->SetBinError(25,0.01735688);
   MVA_BDT_Train_B__41->SetBinError(26,0.01571431);
   MVA_BDT_Train_B__41->SetBinError(27,0.0124741);
   MVA_BDT_Train_B__41->SetBinError(28,0.01022063);
   MVA_BDT_Train_B__41->SetBinError(29,0.006200241);
   MVA_BDT_Train_B__41->SetBinError(30,0.007438865);
   MVA_BDT_Train_B__41->SetBinError(31,0.001878737);
   MVA_BDT_Train_B__41->SetBinError(32,0.004428455);
   MVA_BDT_Train_B__41->SetBinError(33,0.001623258);
   MVA_BDT_Train_B__41->SetBinError(34,0.003131391);
   MVA_BDT_Train_B__41->SetBinError(37,0.0009371884);
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
   TText *text = new TText(0.12,0.74,"Kolmogorov-Smirnov test: signal (background) probability = 0.201 (0.354)");
   text->SetNDC();
   text->SetTextSize(0.032);
   text->Draw();
   
   TH2F *frameBDT__42 = new TH2F("frameBDT__42","TMVA overtraining check for classifier: BDT",500,-0.4938157,0.3476416,500,0,8.805267);
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
