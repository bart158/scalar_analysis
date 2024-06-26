#ifdef __CLING__
#pragma cling optimize(0)
#endif
void train_bdt_qqll_res()
{
//=========Macro generated from canvas: canvas1/TMVA comparison BDT
//=========  (Mon Jul  1 17:46:50 2024) by ROOT version 6.32.02
   TCanvas *canvas1 = new TCanvas("canvas1", "TMVA comparison BDT",200,81,600,468);
   gStyle->SetOptStat(0);
   canvas1->Range(-0.6048397,-1.162767,0.4101682,8.526956);
   canvas1->SetFillColor(0);
   canvas1->SetBorderMode(0);
   canvas1->SetBorderSize(2);
   canvas1->SetTickx(1);
   canvas1->SetTicky(1);
   canvas1->SetRightMargin(0.05);
   canvas1->SetBottomMargin(0.12);
   canvas1->SetFrameBorderMode(0);
   canvas1->SetFrameBorderMode(0);
   
   TH2F *frameBDT__37 = new TH2F("frameBDT__37","TMVA overtraining check for classifier: BDT",500,-0.4982639,0.3594178,500,0,7.557984);
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
   
   TH1D *MVA_BDT_S__38 = new TH1D("MVA_BDT_S__38","TMVA overtraining check for classifier: BDT",40,-0.4982639,0.3594178);
   MVA_BDT_S__38->SetBinContent(15,0.004043);
   MVA_BDT_S__38->SetBinContent(16,0.002695333);
   MVA_BDT_S__38->SetBinContent(17,0.006738333);
   MVA_BDT_S__38->SetBinContent(18,0.03099633);
   MVA_BDT_S__38->SetBinContent(19,0.056602);
   MVA_BDT_S__38->SetBinContent(20,0.09837966);
   MVA_BDT_S__38->SetBinContent(21,0.1684583);
   MVA_BDT_S__38->SetBinContent(22,0.2654903);
   MVA_BDT_S__38->SetBinContent(23,0.464945);
   MVA_BDT_S__38->SetBinContent(24,0.6980913);
   MVA_BDT_S__38->SetBinContent(25,1.030965);
   MVA_BDT_S__38->SetBinContent(26,1.497258);
   MVA_BDT_S__38->SetBinContent(27,1.999937);
   MVA_BDT_S__38->SetBinContent(28,2.706115);
   MVA_BDT_S__38->SetBinContent(29,3.533582);
   MVA_BDT_S__38->SetBinContent(30,4.369135);
   MVA_BDT_S__38->SetBinContent(31,5.555082);
   MVA_BDT_S__38->SetBinContent(32,5.813834);
   MVA_BDT_S__38->SetBinContent(33,5.52139);
   MVA_BDT_S__38->SetBinContent(34,4.77074);
   MVA_BDT_S__38->SetBinContent(35,3.852979);
   MVA_BDT_S__38->SetBinContent(36,2.421757);
   MVA_BDT_S__38->SetBinContent(37,1.231767);
   MVA_BDT_S__38->SetBinContent(38,0.4352963);
   MVA_BDT_S__38->SetBinContent(39,0.09029366);
   MVA_BDT_S__38->SetBinContent(40,0.01078133);
   MVA_BDT_S__38->SetBinError(15,0.002334227);
   MVA_BDT_S__38->SetBinError(16,0.001905888);
   MVA_BDT_S__38->SetBinError(17,0.003013474);
   MVA_BDT_S__38->SetBinError(18,0.006463182);
   MVA_BDT_S__38->SetBinError(19,0.008733878);
   MVA_BDT_S__38->SetBinError(20,0.01151447);
   MVA_BDT_S__38->SetBinError(21,0.01506737);
   MVA_BDT_S__38->SetBinError(22,0.0189154);
   MVA_BDT_S__38->SetBinError(23,0.0250318);
   MVA_BDT_S__38->SetBinError(24,0.03067237);
   MVA_BDT_S__38->SetBinError(25,0.03727462);
   MVA_BDT_S__38->SetBinError(26,0.04491998);
   MVA_BDT_S__38->SetBinError(27,0.05191579);
   MVA_BDT_S__38->SetBinError(28,0.0603899);
   MVA_BDT_S__38->SetBinError(29,0.0690079);
   MVA_BDT_S__38->SetBinError(30,0.0767342);
   MVA_BDT_S__38->SetBinError(31,0.08652398);
   MVA_BDT_S__38->SetBinError(32,0.08851616);
   MVA_BDT_S__38->SetBinError(33,0.08626119);
   MVA_BDT_S__38->SetBinError(34,0.08018333);
   MVA_BDT_S__38->SetBinError(35,0.07205922);
   MVA_BDT_S__38->SetBinError(36,0.05712899);
   MVA_BDT_S__38->SetBinError(37,0.04074324);
   MVA_BDT_S__38->SetBinError(38,0.02422054);
   MVA_BDT_S__38->SetBinError(39,0.01103113);
   MVA_BDT_S__38->SetBinError(40,0.003811777);
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
   
   TH1D *MVA_BDT_B__39 = new TH1D("MVA_BDT_B__39","MVA_BDT_B",40,-0.4982639,0.3594178);
   MVA_BDT_B__39->SetBinContent(1,0.007224541);
   MVA_BDT_B__39->SetBinContent(2,0.06632255);
   MVA_BDT_B__39->SetBinContent(3,0.1464214);
   MVA_BDT_B__39->SetBinContent(4,0.3769706);
   MVA_BDT_B__39->SetBinContent(5,0.6926446);
   MVA_BDT_B__39->SetBinContent(6,1.251588);
   MVA_BDT_B__39->SetBinContent(7,2.095278);
   MVA_BDT_B__39->SetBinContent(8,2.450819);
   MVA_BDT_B__39->SetBinContent(9,2.634547);
   MVA_BDT_B__39->SetBinContent(10,2.784931);
   MVA_BDT_B__39->SetBinContent(11,2.866111);
   MVA_BDT_B__39->SetBinContent(12,3.088182);
   MVA_BDT_B__39->SetBinContent(13,3.08449);
   MVA_BDT_B__39->SetBinContent(14,3.257522);
   MVA_BDT_B__39->SetBinContent(15,3.17674);
   MVA_BDT_B__39->SetBinContent(16,3.271092);
   MVA_BDT_B__39->SetBinContent(17,2.829067);
   MVA_BDT_B__39->SetBinContent(18,2.593738);
   MVA_BDT_B__39->SetBinContent(19,2.227156);
   MVA_BDT_B__39->SetBinContent(20,1.990085);
   MVA_BDT_B__39->SetBinContent(21,1.648794);
   MVA_BDT_B__39->SetBinContent(22,1.228532);
   MVA_BDT_B__39->SetBinContent(23,0.9504291);
   MVA_BDT_B__39->SetBinContent(24,0.6284967);
   MVA_BDT_B__39->SetBinContent(25,0.4865343);
   MVA_BDT_B__39->SetBinContent(26,0.2699742);
   MVA_BDT_B__39->SetBinContent(27,0.2017005);
   MVA_BDT_B__39->SetBinContent(28,0.1029482);
   MVA_BDT_B__39->SetBinContent(29,0.07421874);
   MVA_BDT_B__39->SetBinContent(30,0.05074746);
   MVA_BDT_B__39->SetBinContent(31,0.02744097);
   MVA_BDT_B__39->SetBinContent(32,0.03466707);
   MVA_BDT_B__39->SetBinContent(33,0.02706577);
   MVA_BDT_B__39->SetBinContent(34,0.009914946);
   MVA_BDT_B__39->SetBinContent(35,0.0007969703);
   MVA_BDT_B__39->SetBinContent(36,0.004160502);
   MVA_BDT_B__39->SetBinError(1,0.004344389);
   MVA_BDT_B__39->SetBinError(2,0.01481734);
   MVA_BDT_B__39->SetBinError(3,0.02083517);
   MVA_BDT_B__39->SetBinError(4,0.03522844);
   MVA_BDT_B__39->SetBinError(5,0.04323663);
   MVA_BDT_B__39->SetBinError(6,0.05453128);
   MVA_BDT_B__39->SetBinError(7,0.06991948);
   MVA_BDT_B__39->SetBinError(8,0.07854749);
   MVA_BDT_B__39->SetBinError(9,0.08416898);
   MVA_BDT_B__39->SetBinError(10,0.08997157);
   MVA_BDT_B__39->SetBinError(11,0.09120576);
   MVA_BDT_B__39->SetBinError(12,0.09796059);
   MVA_BDT_B__39->SetBinError(13,0.09645583);
   MVA_BDT_B__39->SetBinError(14,0.1012604);
   MVA_BDT_B__39->SetBinError(15,0.1005411);
   MVA_BDT_B__39->SetBinError(16,0.1031983);
   MVA_BDT_B__39->SetBinError(17,0.09610487);
   MVA_BDT_B__39->SetBinError(18,0.09315485);
   MVA_BDT_B__39->SetBinError(19,0.08632371);
   MVA_BDT_B__39->SetBinError(20,0.08200617);
   MVA_BDT_B__39->SetBinError(21,0.07478888);
   MVA_BDT_B__39->SetBinError(22,0.06491339);
   MVA_BDT_B__39->SetBinError(23,0.05764826);
   MVA_BDT_B__39->SetBinError(24,0.04731366);
   MVA_BDT_B__39->SetBinError(25,0.04162019);
   MVA_BDT_B__39->SetBinError(26,0.03082156);
   MVA_BDT_B__39->SetBinError(27,0.02668396);
   MVA_BDT_B__39->SetBinError(28,0.0189506);
   MVA_BDT_B__39->SetBinError(29,0.0163111);
   MVA_BDT_B__39->SetBinError(30,0.01335657);
   MVA_BDT_B__39->SetBinError(31,0.0102581);
   MVA_BDT_B__39->SetBinError(32,0.01178939);
   MVA_BDT_B__39->SetBinError(33,0.01025467);
   MVA_BDT_B__39->SetBinError(34,0.005937569);
   MVA_BDT_B__39->SetBinError(35,0.0005635431);
   MVA_BDT_B__39->SetBinError(36,0.004160502);
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
   
   TH1D *MVA_BDT_Train_S__40 = new TH1D("MVA_BDT_Train_S__40","MVA_BDT_Train_S",40,-0.4982639,0.3594178);
   MVA_BDT_Train_S__40->SetBinContent(12,0.001347667);
   MVA_BDT_Train_S__40->SetBinContent(14,0.004043);
   MVA_BDT_Train_S__40->SetBinContent(15,0.002695333);
   MVA_BDT_Train_S__40->SetBinContent(16,0.004043);
   MVA_BDT_Train_S__40->SetBinContent(17,0.004043);
   MVA_BDT_Train_S__40->SetBinContent(18,0.016172);
   MVA_BDT_Train_S__40->SetBinContent(19,0.02156267);
   MVA_BDT_Train_S__40->SetBinContent(20,0.05121133);
   MVA_BDT_Train_S__40->SetBinContent(21,0.117247);
   MVA_BDT_Train_S__40->SetBinContent(22,0.258752);
   MVA_BDT_Train_S__40->SetBinContent(23,0.4662927);
   MVA_BDT_Train_S__40->SetBinContent(24,0.7412167);
   MVA_BDT_Train_S__40->SetBinContent(25,0.9918826);
   MVA_BDT_Train_S__40->SetBinContent(26,1.411007);
   MVA_BDT_Train_S__40->SetBinContent(27,1.989156);
   MVA_BDT_Train_S__40->SetBinContent(28,2.664337);
   MVA_BDT_Train_S__40->SetBinContent(29,3.506629);
   MVA_BDT_Train_S__40->SetBinContent(30,4.68988);
   MVA_BDT_Train_S__40->SetBinContent(31,5.293635);
   MVA_BDT_Train_S__40->SetBinContent(32,5.750494);
   MVA_BDT_Train_S__40->SetBinContent(33,5.608989);
   MVA_BDT_Train_S__40->SetBinContent(34,4.951327);
   MVA_BDT_Train_S__40->SetBinContent(35,3.710126);
   MVA_BDT_Train_S__40->SetBinContent(36,2.516094);
   MVA_BDT_Train_S__40->SetBinContent(37,1.234463);
   MVA_BDT_Train_S__40->SetBinContent(38,0.505375);
   MVA_BDT_Train_S__40->SetBinContent(39,0.109161);
   MVA_BDT_Train_S__40->SetBinContent(40,0.016172);
   MVA_BDT_Train_S__40->SetBinError(12,0.001347667);
   MVA_BDT_Train_S__40->SetBinError(14,0.002334227);
   MVA_BDT_Train_S__40->SetBinError(15,0.001905888);
   MVA_BDT_Train_S__40->SetBinError(16,0.002334227);
   MVA_BDT_Train_S__40->SetBinError(17,0.002334227);
   MVA_BDT_Train_S__40->SetBinError(18,0.004668454);
   MVA_BDT_Train_S__40->SetBinError(19,0.005390667);
   MVA_BDT_Train_S__40->SetBinError(20,0.008307575);
   MVA_BDT_Train_S__40->SetBinError(21,0.0125702);
   MVA_BDT_Train_S__40->SetBinError(22,0.01867382);
   MVA_BDT_Train_S__40->SetBinError(23,0.02506805);
   MVA_BDT_Train_S__40->SetBinError(24,0.03160558);
   MVA_BDT_Train_S__40->SetBinError(25,0.03656128);
   MVA_BDT_Train_S__40->SetBinError(26,0.04360696);
   MVA_BDT_Train_S__40->SetBinError(27,0.05177566);
   MVA_BDT_Train_S__40->SetBinError(28,0.05992193);
   MVA_BDT_Train_S__40->SetBinError(29,0.06874421);
   MVA_BDT_Train_S__40->SetBinError(30,0.07950091);
   MVA_BDT_Train_S__40->SetBinError(31,0.08446333);
   MVA_BDT_Train_S__40->SetBinError(32,0.08803265);
   MVA_BDT_Train_S__40->SetBinError(33,0.08694278);
   MVA_BDT_Train_S__40->SetBinError(34,0.08168683);
   MVA_BDT_Train_S__40->SetBinError(35,0.07071077);
   MVA_BDT_Train_S__40->SetBinError(36,0.05823105);
   MVA_BDT_Train_S__40->SetBinError(37,0.04078779);
   MVA_BDT_Train_S__40->SetBinError(38,0.02609745);
   MVA_BDT_Train_S__40->SetBinError(39,0.012129);
   MVA_BDT_Train_S__40->SetBinError(40,0.004668454);
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
   
   TH1D *MVA_BDT_Train_B__41 = new TH1D("MVA_BDT_Train_B__41","MVA_BDT_Train_B",40,-0.4982639,0.3594178);
   MVA_BDT_Train_B__41->SetBinContent(1,0.0189574);
   MVA_BDT_Train_B__41->SetBinContent(2,0.05801508);
   MVA_BDT_Train_B__41->SetBinContent(3,0.1522929);
   MVA_BDT_Train_B__41->SetBinContent(4,0.305161);
   MVA_BDT_Train_B__41->SetBinContent(5,0.6166843);
   MVA_BDT_Train_B__41->SetBinContent(6,1.292731);
   MVA_BDT_Train_B__41->SetBinContent(7,2.091939);
   MVA_BDT_Train_B__41->SetBinContent(8,2.589301);
   MVA_BDT_Train_B__41->SetBinContent(9,2.679255);
   MVA_BDT_Train_B__41->SetBinContent(10,2.738801);
   MVA_BDT_Train_B__41->SetBinContent(11,2.953807);
   MVA_BDT_Train_B__41->SetBinContent(12,2.946466);
   MVA_BDT_Train_B__41->SetBinContent(13,3.022924);
   MVA_BDT_Train_B__41->SetBinContent(14,3.055721);
   MVA_BDT_Train_B__41->SetBinContent(15,3.373556);
   MVA_BDT_Train_B__41->SetBinContent(16,3.216527);
   MVA_BDT_Train_B__41->SetBinContent(17,2.872885);
   MVA_BDT_Train_B__41->SetBinContent(18,2.673545);
   MVA_BDT_Train_B__41->SetBinContent(19,2.317302);
   MVA_BDT_Train_B__41->SetBinContent(20,2.244265);
   MVA_BDT_Train_B__41->SetBinContent(21,1.607259);
   MVA_BDT_Train_B__41->SetBinContent(22,1.37285);
   MVA_BDT_Train_B__41->SetBinContent(23,1.067029);
   MVA_BDT_Train_B__41->SetBinContent(24,0.5798852);
   MVA_BDT_Train_B__41->SetBinContent(25,0.3254967);
   MVA_BDT_Train_B__41->SetBinContent(26,0.1655172);
   MVA_BDT_Train_B__41->SetBinContent(27,0.1237048);
   MVA_BDT_Train_B__41->SetBinContent(28,0.06519608);
   MVA_BDT_Train_B__41->SetBinContent(29,0.0356247);
   MVA_BDT_Train_B__41->SetBinContent(30,0.03248007);
   MVA_BDT_Train_B__41->SetBinContent(31,0.01106987);
   MVA_BDT_Train_B__41->SetBinContent(32,0.01304175);
   MVA_BDT_Train_B__41->SetBinContent(33,0.009232484);
   MVA_BDT_Train_B__41->SetBinContent(34,0.004616242);
   MVA_BDT_Train_B__41->SetBinContent(35,0.004212752);
   MVA_BDT_Train_B__41->SetBinError(1,0.007657117);
   MVA_BDT_Train_B__41->SetBinError(2,0.01374546);
   MVA_BDT_Train_B__41->SetBinError(3,0.02245802);
   MVA_BDT_Train_B__41->SetBinError(4,0.03040393);
   MVA_BDT_Train_B__41->SetBinError(5,0.04040558);
   MVA_BDT_Train_B__41->SetBinError(6,0.05582314);
   MVA_BDT_Train_B__41->SetBinError(7,0.07177301);
   MVA_BDT_Train_B__41->SetBinError(8,0.08162376);
   MVA_BDT_Train_B__41->SetBinError(9,0.08634468);
   MVA_BDT_Train_B__41->SetBinError(10,0.08828472);
   MVA_BDT_Train_B__41->SetBinError(11,0.09393174);
   MVA_BDT_Train_B__41->SetBinError(12,0.09361675);
   MVA_BDT_Train_B__41->SetBinError(13,0.09592672);
   MVA_BDT_Train_B__41->SetBinError(14,0.09718921);
   MVA_BDT_Train_B__41->SetBinError(15,0.1042334);
   MVA_BDT_Train_B__41->SetBinError(16,0.1027007);
   MVA_BDT_Train_B__41->SetBinError(17,0.09753371);
   MVA_BDT_Train_B__41->SetBinError(18,0.09454994);
   MVA_BDT_Train_B__41->SetBinError(19,0.08830272);
   MVA_BDT_Train_B__41->SetBinError(20,0.08879869);
   MVA_BDT_Train_B__41->SetBinError(21,0.07396515);
   MVA_BDT_Train_B__41->SetBinError(22,0.06986852);
   MVA_BDT_Train_B__41->SetBinError(23,0.06208935);
   MVA_BDT_Train_B__41->SetBinError(24,0.04444536);
   MVA_BDT_Train_B__41->SetBinError(25,0.03283713);
   MVA_BDT_Train_B__41->SetBinError(26,0.02162965);
   MVA_BDT_Train_B__41->SetBinError(27,0.02062842);
   MVA_BDT_Train_B__41->SetBinError(28,0.01486435);
   MVA_BDT_Train_B__41->SetBinError(29,0.0112652);
   MVA_BDT_Train_B__41->SetBinError(30,0.01119447);
   MVA_BDT_Train_B__41->SetBinError(31,0.006124391);
   MVA_BDT_Train_B__41->SetBinError(32,0.007307849);
   MVA_BDT_Train_B__41->SetBinError(33,0.005984996);
   MVA_BDT_Train_B__41->SetBinError(34,0.004232031);
   MVA_BDT_Train_B__41->SetBinError(35,0.004212752);
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
   TText *text = new TText(0.12,0.74,"Kolmogorov-Smirnov test: signal (background) probability = 0.203 (0.144)");
   text->SetNDC();
   text->SetTextSize(0.032);
   text->Draw();
   
   TH2F *frameBDT__42 = new TH2F("frameBDT__42","TMVA overtraining check for classifier: BDT",500,-0.4982639,0.3594178,500,0,7.557984);
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
