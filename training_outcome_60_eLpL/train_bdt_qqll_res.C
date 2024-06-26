#ifdef __CLING__
#pragma cling optimize(0)
#endif
void train_bdt_qqll_res()
{
//=========Macro generated from canvas: canvas1/TMVA comparison BDT
//=========  (Mon Jul  1 17:49:02 2024) by ROOT version 6.32.02
   TCanvas *canvas1 = new TCanvas("canvas1", "TMVA comparison BDT",200,81,600,468);
   gStyle->SetOptStat(0);
   canvas1->Range(-0.5917432,-1.1261,0.4359098,8.258068);
   canvas1->SetFillColor(0);
   canvas1->SetBorderMode(0);
   canvas1->SetBorderSize(2);
   canvas1->SetTickx(1);
   canvas1->SetTicky(1);
   canvas1->SetRightMargin(0.05);
   canvas1->SetBottomMargin(0.12);
   canvas1->SetFrameBorderMode(0);
   canvas1->SetFrameBorderMode(0);
   
   TH2F *frameBDT__37 = new TH2F("frameBDT__37","TMVA overtraining check for classifier: BDT",500,-0.4838397,0.3845271,500,0,7.319651);
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
   
   TH1D *MVA_BDT_S__38 = new TH1D("MVA_BDT_S__38","TMVA overtraining check for classifier: BDT",40,-0.4838397,0.3845271);
   MVA_BDT_S__38->SetBinContent(11,0.001263294);
   MVA_BDT_S__38->SetBinContent(12,0.001263294);
   MVA_BDT_S__38->SetBinContent(13,0.007579764);
   MVA_BDT_S__38->SetBinContent(14,0.003789882);
   MVA_BDT_S__38->SetBinContent(15,0.007579764);
   MVA_BDT_S__38->SetBinContent(16,0.021476);
   MVA_BDT_S__38->SetBinContent(17,0.02779247);
   MVA_BDT_S__38->SetBinContent(18,0.05053176);
   MVA_BDT_S__38->SetBinContent(19,0.0619014);
   MVA_BDT_S__38->SetBinContent(20,0.1225395);
   MVA_BDT_S__38->SetBinContent(21,0.1617016);
   MVA_BDT_S__38->SetBinContent(22,0.288031);
   MVA_BDT_S__38->SetBinContent(23,0.4497326);
   MVA_BDT_S__38->SetBinContent(24,0.696075);
   MVA_BDT_S__38->SetBinContent(25,0.9449439);
   MVA_BDT_S__38->SetBinContent(26,1.358041);
   MVA_BDT_S__38->SetBinContent(27,1.855779);
   MVA_BDT_S__38->SetBinContent(28,2.582173);
   MVA_BDT_S__38->SetBinContent(29,3.53217);
   MVA_BDT_S__38->SetBinContent(30,4.419002);
   MVA_BDT_S__38->SetBinContent(31,5.194665);
   MVA_BDT_S__38->SetBinContent(32,5.630501);
   MVA_BDT_S__38->SetBinContent(33,5.5686);
   MVA_BDT_S__38->SetBinContent(34,4.887684);
   MVA_BDT_S__38->SetBinContent(35,3.887155);
   MVA_BDT_S__38->SetBinContent(36,2.463423);
   MVA_BDT_S__38->SetBinContent(37,1.29614);
   MVA_BDT_S__38->SetBinContent(38,0.4358364);
   MVA_BDT_S__38->SetBinContent(39,0.09601034);
   MVA_BDT_S__38->SetBinContent(40,0.01010635);
   MVA_BDT_S__38->SetBinError(11,0.001263294);
   MVA_BDT_S__38->SetBinError(12,0.001263294);
   MVA_BDT_S__38->SetBinError(13,0.003094426);
   MVA_BDT_S__38->SetBinError(14,0.002188089);
   MVA_BDT_S__38->SetBinError(15,0.003094426);
   MVA_BDT_S__38->SetBinError(16,0.005208694);
   MVA_BDT_S__38->SetBinError(17,0.005925374);
   MVA_BDT_S__38->SetBinError(18,0.007989772);
   MVA_BDT_S__38->SetBinError(19,0.008843058);
   MVA_BDT_S__38->SetBinError(20,0.012442);
   MVA_BDT_S__38->SetBinError(21,0.01429254);
   MVA_BDT_S__38->SetBinError(22,0.01907532);
   MVA_BDT_S__38->SetBinError(23,0.02383578);
   MVA_BDT_S__38->SetBinError(24,0.02965379);
   MVA_BDT_S__38->SetBinError(25,0.03455057);
   MVA_BDT_S__38->SetBinError(26,0.04141986);
   MVA_BDT_S__38->SetBinError(27,0.04841894);
   MVA_BDT_S__38->SetBinError(28,0.0571143);
   MVA_BDT_S__38->SetBinError(29,0.06679947);
   MVA_BDT_S__38->SetBinError(30,0.07471612);
   MVA_BDT_S__38->SetBinError(31,0.08100857);
   MVA_BDT_S__38->SetBinError(32,0.08433847);
   MVA_BDT_S__38->SetBinError(33,0.08387359);
   MVA_BDT_S__38->SetBinError(34,0.07857851);
   MVA_BDT_S__38->SetBinError(35,0.07007582);
   MVA_BDT_S__38->SetBinError(36,0.05578555);
   MVA_BDT_S__38->SetBinError(37,0.04046487);
   MVA_BDT_S__38->SetBinError(38,0.02346464);
   MVA_BDT_S__38->SetBinError(39,0.01101314);
   MVA_BDT_S__38->SetBinError(40,0.003573135);
   MVA_BDT_S__38->SetEntries(36463);
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
   
   TH1D *MVA_BDT_B__39 = new TH1D("MVA_BDT_B__39","MVA_BDT_B",40,-0.4838397,0.3845271);
   MVA_BDT_B__39->SetBinContent(1,0.01234248);
   MVA_BDT_B__39->SetBinContent(2,0.05849179);
   MVA_BDT_B__39->SetBinContent(3,0.1539983);
   MVA_BDT_B__39->SetBinContent(4,0.5133953);
   MVA_BDT_B__39->SetBinContent(5,1.223615);
   MVA_BDT_B__39->SetBinContent(6,1.855745);
   MVA_BDT_B__39->SetBinContent(7,2.423863);
   MVA_BDT_B__39->SetBinContent(8,2.778893);
   MVA_BDT_B__39->SetBinContent(9,3.078536);
   MVA_BDT_B__39->SetBinContent(10,3.365229);
   MVA_BDT_B__39->SetBinContent(11,3.452015);
   MVA_BDT_B__39->SetBinContent(12,3.509764);
   MVA_BDT_B__39->SetBinContent(13,3.81026);
   MVA_BDT_B__39->SetBinContent(14,3.409582);
   MVA_BDT_B__39->SetBinContent(15,3.305989);
   MVA_BDT_B__39->SetBinContent(16,2.798696);
   MVA_BDT_B__39->SetBinContent(17,2.542593);
   MVA_BDT_B__39->SetBinContent(18,1.97061);
   MVA_BDT_B__39->SetBinContent(19,1.693027);
   MVA_BDT_B__39->SetBinContent(20,1.295084);
   MVA_BDT_B__39->SetBinContent(21,0.8963651);
   MVA_BDT_B__39->SetBinContent(22,0.6271354);
   MVA_BDT_B__39->SetBinContent(23,0.4351398);
   MVA_BDT_B__39->SetBinContent(24,0.315845);
   MVA_BDT_B__39->SetBinContent(25,0.1747154);
   MVA_BDT_B__39->SetBinContent(26,0.1238445);
   MVA_BDT_B__39->SetBinContent(27,0.06750003);
   MVA_BDT_B__39->SetBinContent(28,0.05933599);
   MVA_BDT_B__39->SetBinContent(29,0.03582334);
   MVA_BDT_B__39->SetBinContent(30,0.04639089);
   MVA_BDT_B__39->SetBinContent(31,0.01166198);
   MVA_BDT_B__39->SetBinContent(32,0.0117021);
   MVA_BDT_B__39->SetBinContent(33,0.001118568);
   MVA_BDT_B__39->SetBinContent(34,0.004763729);
   MVA_BDT_B__39->SetBinContent(35,0.0004174913);
   MVA_BDT_B__39->SetBinError(1,0.005324127);
   MVA_BDT_B__39->SetBinError(2,0.01226178);
   MVA_BDT_B__39->SetBinError(3,0.01739057);
   MVA_BDT_B__39->SetBinError(4,0.02987309);
   MVA_BDT_B__39->SetBinError(5,0.04791968);
   MVA_BDT_B__39->SetBinError(6,0.06040655);
   MVA_BDT_B__39->SetBinError(7,0.07184029);
   MVA_BDT_B__39->SetBinError(8,0.07903145);
   MVA_BDT_B__39->SetBinError(9,0.08425792);
   MVA_BDT_B__39->SetBinError(10,0.0897197);
   MVA_BDT_B__39->SetBinError(11,0.09129548);
   MVA_BDT_B__39->SetBinError(12,0.09317929);
   MVA_BDT_B__39->SetBinError(13,0.09859871);
   MVA_BDT_B__39->SetBinError(14,0.09318785);
   MVA_BDT_B__39->SetBinError(15,0.09273471);
   MVA_BDT_B__39->SetBinError(16,0.08484314);
   MVA_BDT_B__39->SetBinError(17,0.08182527);
   MVA_BDT_B__39->SetBinError(18,0.07123595);
   MVA_BDT_B__39->SetBinError(19,0.06679743);
   MVA_BDT_B__39->SetBinError(20,0.05845408);
   MVA_BDT_B__39->SetBinError(21,0.04818043);
   MVA_BDT_B__39->SetBinError(22,0.04030333);
   MVA_BDT_B__39->SetBinError(23,0.03334504);
   MVA_BDT_B__39->SetBinError(24,0.02869058);
   MVA_BDT_B__39->SetBinError(25,0.02130777);
   MVA_BDT_B__39->SetBinError(26,0.01788875);
   MVA_BDT_B__39->SetBinError(27,0.01268838);
   MVA_BDT_B__39->SetBinError(28,0.01291213);
   MVA_BDT_B__39->SetBinError(29,0.01002695);
   MVA_BDT_B__39->SetBinError(30,0.01149875);
   MVA_BDT_B__39->SetBinError(31,0.004972812);
   MVA_BDT_B__39->SetBinError(32,0.004858737);
   MVA_BDT_B__39->SetBinError(33,0.0006549951);
   MVA_BDT_B__39->SetBinError(34,0.003579889);
   MVA_BDT_B__39->SetBinError(35,0.0004174913);
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
   
   TH1D *MVA_BDT_Train_S__40 = new TH1D("MVA_BDT_Train_S__40","MVA_BDT_Train_S",40,-0.4838397,0.3845271);
   MVA_BDT_Train_S__40->SetBinContent(13,0.001263294);
   MVA_BDT_Train_S__40->SetBinContent(14,0.001263294);
   MVA_BDT_Train_S__40->SetBinContent(15,0.00631647);
   MVA_BDT_Train_S__40->SetBinContent(16,0.007579764);
   MVA_BDT_Train_S__40->SetBinContent(17,0.01515953);
   MVA_BDT_Train_S__40->SetBinContent(18,0.01642282);
   MVA_BDT_Train_S__40->SetBinContent(19,0.03284564);
   MVA_BDT_Train_S__40->SetBinContent(20,0.07074446);
   MVA_BDT_Train_S__40->SetBinContent(21,0.1313826);
   MVA_BDT_Train_S__40->SetBinContent(22,0.3410894);
   MVA_BDT_Train_S__40->SetBinContent(23,0.4535225);
   MVA_BDT_Train_S__40->SetBinContent(24,0.7061813);
   MVA_BDT_Train_S__40->SetBinContent(25,0.9613667);
   MVA_BDT_Train_S__40->SetBinContent(26,1.408573);
   MVA_BDT_Train_S__40->SetBinContent(27,2.013691);
   MVA_BDT_Train_S__40->SetBinContent(28,2.717345);
   MVA_BDT_Train_S__40->SetBinContent(29,3.533433);
   MVA_BDT_Train_S__40->SetBinContent(30,4.369734);
   MVA_BDT_Train_S__40->SetBinContent(31,5.176979);
   MVA_BDT_Train_S__40->SetBinContent(32,5.44606);
   MVA_BDT_Train_S__40->SetBinContent(33,5.483959);
   MVA_BDT_Train_S__40->SetBinContent(34,4.719666);
   MVA_BDT_Train_S__40->SetBinContent(35,3.882102);
   MVA_BDT_Train_S__40->SetBinContent(36,2.66934);
   MVA_BDT_Train_S__40->SetBinContent(37,1.39973);
   MVA_BDT_Train_S__40->SetBinContent(38,0.4093072);
   MVA_BDT_Train_S__40->SetBinContent(39,0.07958752);
   MVA_BDT_Train_S__40->SetBinContent(40,0.008843058);
   MVA_BDT_Train_S__40->SetBinError(13,0.001263294);
   MVA_BDT_Train_S__40->SetBinError(14,0.001263294);
   MVA_BDT_Train_S__40->SetBinError(15,0.002824811);
   MVA_BDT_Train_S__40->SetBinError(16,0.003094426);
   MVA_BDT_Train_S__40->SetBinError(17,0.004376179);
   MVA_BDT_Train_S__40->SetBinError(18,0.004554871);
   MVA_BDT_Train_S__40->SetBinError(19,0.006441561);
   MVA_BDT_Train_S__40->SetBinError(20,0.009453626);
   MVA_BDT_Train_S__40->SetBinError(21,0.01288312);
   MVA_BDT_Train_S__40->SetBinError(22,0.02075804);
   MVA_BDT_Train_S__40->SetBinError(23,0.023936);
   MVA_BDT_Train_S__40->SetBinError(24,0.02986829);
   MVA_BDT_Train_S__40->SetBinError(25,0.03484952);
   MVA_BDT_Train_S__40->SetBinError(26,0.04218343);
   MVA_BDT_Train_S__40->SetBinError(27,0.05043692);
   MVA_BDT_Train_S__40->SetBinError(28,0.05859015);
   MVA_BDT_Train_S__40->SetBinError(29,0.06681141);
   MVA_BDT_Train_S__40->SetBinError(30,0.07429844);
   MVA_BDT_Train_S__40->SetBinError(31,0.08087055);
   MVA_BDT_Train_S__40->SetBinError(32,0.08294561);
   MVA_BDT_Train_S__40->SetBinError(33,0.08323372);
   MVA_BDT_Train_S__40->SetBinError(34,0.0772161);
   MVA_BDT_Train_S__40->SetBinError(35,0.07003025);
   MVA_BDT_Train_S__40->SetBinError(36,0.05807031);
   MVA_BDT_Train_S__40->SetBinError(37,0.0420508);
   MVA_BDT_Train_S__40->SetBinError(38,0.02273929);
   MVA_BDT_Train_S__40->SetBinError(39,0.01002708);
   MVA_BDT_Train_S__40->SetBinError(40,0.003342362);
   MVA_BDT_Train_S__40->SetEntries(36463);

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
   
   TH1D *MVA_BDT_Train_B__41 = new TH1D("MVA_BDT_Train_B__41","MVA_BDT_Train_B",40,-0.4838397,0.3845271);
   MVA_BDT_Train_B__41->SetBinContent(1,0.00376343);
   MVA_BDT_Train_B__41->SetBinContent(2,0.04188155);
   MVA_BDT_Train_B__41->SetBinContent(3,0.1697006);
   MVA_BDT_Train_B__41->SetBinContent(4,0.5427116);
   MVA_BDT_Train_B__41->SetBinContent(5,1.249639);
   MVA_BDT_Train_B__41->SetBinContent(6,1.877947);
   MVA_BDT_Train_B__41->SetBinContent(7,2.34648);
   MVA_BDT_Train_B__41->SetBinContent(8,2.894382);
   MVA_BDT_Train_B__41->SetBinContent(9,3.109834);
   MVA_BDT_Train_B__41->SetBinContent(10,3.482942);
   MVA_BDT_Train_B__41->SetBinContent(11,3.472862);
   MVA_BDT_Train_B__41->SetBinContent(12,3.586901);
   MVA_BDT_Train_B__41->SetBinContent(13,3.55422);
   MVA_BDT_Train_B__41->SetBinContent(14,3.347381);
   MVA_BDT_Train_B__41->SetBinContent(15,3.099288);
   MVA_BDT_Train_B__41->SetBinContent(16,3.078274);
   MVA_BDT_Train_B__41->SetBinContent(17,2.597494);
   MVA_BDT_Train_B__41->SetBinContent(18,1.940319);
   MVA_BDT_Train_B__41->SetBinContent(19,1.833573);
   MVA_BDT_Train_B__41->SetBinContent(20,1.375956);
   MVA_BDT_Train_B__41->SetBinContent(21,0.9343806);
   MVA_BDT_Train_B__41->SetBinContent(22,0.6494203);
   MVA_BDT_Train_B__41->SetBinContent(23,0.416775);
   MVA_BDT_Train_B__41->SetBinContent(24,0.1851867);
   MVA_BDT_Train_B__41->SetBinContent(25,0.08678587);
   MVA_BDT_Train_B__41->SetBinContent(26,0.05579173);
   MVA_BDT_Train_B__41->SetBinContent(27,0.03983022);
   MVA_BDT_Train_B__41->SetBinContent(28,0.0294645);
   MVA_BDT_Train_B__41->SetBinContent(29,0.01853798);
   MVA_BDT_Train_B__41->SetBinContent(30,0.01927912);
   MVA_BDT_Train_B__41->SetBinContent(31,0.01876088);
   MVA_BDT_Train_B__41->SetBinContent(32,0.0005116852);
   MVA_BDT_Train_B__41->SetBinContent(33,0.002374419);
   MVA_BDT_Train_B__41->SetBinContent(35,0.0004200772);
   MVA_BDT_Train_B__41->SetBinContent(36,0.0004200772);
   MVA_BDT_Train_B__41->SetBinError(1,0.002215756);
   MVA_BDT_Train_B__41->SetBinError(2,0.009971545);
   MVA_BDT_Train_B__41->SetBinError(3,0.01804026);
   MVA_BDT_Train_B__41->SetBinError(4,0.03141008);
   MVA_BDT_Train_B__41->SetBinError(5,0.04838586);
   MVA_BDT_Train_B__41->SetBinError(6,0.05972764);
   MVA_BDT_Train_B__41->SetBinError(7,0.07054843);
   MVA_BDT_Train_B__41->SetBinError(8,0.08136626);
   MVA_BDT_Train_B__41->SetBinError(9,0.08524488);
   MVA_BDT_Train_B__41->SetBinError(10,0.09143494);
   MVA_BDT_Train_B__41->SetBinError(11,0.09178284);
   MVA_BDT_Train_B__41->SetBinError(12,0.09438351);
   MVA_BDT_Train_B__41->SetBinError(13,0.09454436);
   MVA_BDT_Train_B__41->SetBinError(14,0.09226011);
   MVA_BDT_Train_B__41->SetBinError(15,0.08983101);
   MVA_BDT_Train_B__41->SetBinError(16,0.09066588);
   MVA_BDT_Train_B__41->SetBinError(17,0.08321611);
   MVA_BDT_Train_B__41->SetBinError(18,0.07102151);
   MVA_BDT_Train_B__41->SetBinError(19,0.07048223);
   MVA_BDT_Train_B__41->SetBinError(20,0.06111828);
   MVA_BDT_Train_B__41->SetBinError(21,0.04919236);
   MVA_BDT_Train_B__41->SetBinError(22,0.04127149);
   MVA_BDT_Train_B__41->SetBinError(23,0.03251874);
   MVA_BDT_Train_B__41->SetBinError(24,0.02025247);
   MVA_BDT_Train_B__41->SetBinError(25,0.01342685);
   MVA_BDT_Train_B__41->SetBinError(26,0.01037618);
   MVA_BDT_Train_B__41->SetBinError(27,0.009945727);
   MVA_BDT_Train_B__41->SetBinError(28,0.00863655);
   MVA_BDT_Train_B__41->SetBinError(29,0.006814833);
   MVA_BDT_Train_B__41->SetBinError(30,0.007518225);
   MVA_BDT_Train_B__41->SetBinError(31,0.007471261);
   MVA_BDT_Train_B__41->SetBinError(32,0.0004299498);
   MVA_BDT_Train_B__41->SetBinError(33,0.001645265);
   MVA_BDT_Train_B__41->SetBinError(35,0.0004200772);
   MVA_BDT_Train_B__41->SetBinError(36,0.0004200772);
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
   TText *text = new TText(0.12,0.74,"Kolmogorov-Smirnov test: signal (background) probability = 0.291 (0.244)");
   text->SetNDC();
   text->SetTextSize(0.032);
   text->Draw();
   
   TH2F *frameBDT__42 = new TH2F("frameBDT__42","TMVA overtraining check for classifier: BDT",500,-0.4838397,0.3845271,500,0,7.319651);
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
