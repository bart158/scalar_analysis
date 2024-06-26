#ifdef __CLING__
#pragma cling optimize(0)
#endif
void train_bdt_qqll_res()
{
//=========Macro generated from canvas: canvas1/TMVA comparison BDT
//=========  (Mon Jul  1 17:40:15 2024) by ROOT version 6.32.02
   TCanvas *canvas1 = new TCanvas("canvas1", "TMVA comparison BDT",200,81,600,468);
   gStyle->SetOptStat(0);
   canvas1->Range(-0.6226496,-1.450859,0.3621557,10.63963);
   canvas1->SetFillColor(0);
   canvas1->SetBorderMode(0);
   canvas1->SetBorderSize(2);
   canvas1->SetTickx(1);
   canvas1->SetTicky(1);
   canvas1->SetRightMargin(0.05);
   canvas1->SetBottomMargin(0.12);
   canvas1->SetFrameBorderMode(0);
   canvas1->SetFrameBorderMode(0);
   
   TH2F *frameBDT__37 = new TH2F("frameBDT__37","TMVA overtraining check for classifier: BDT",500,-0.519245,0.3129154,500,0,9.430585);
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
   
   TH1D *MVA_BDT_S__38 = new TH1D("MVA_BDT_S__38","TMVA overtraining check for classifier: BDT",40,-0.519245,0.3129154);
   MVA_BDT_S__38->SetBinContent(9,0.0002092104);
   MVA_BDT_S__38->SetBinContent(13,0.0002092104);
   MVA_BDT_S__38->SetBinContent(15,0.0002092104);
   MVA_BDT_S__38->SetBinContent(16,0.004900339);
   MVA_BDT_S__38->SetBinContent(17,0.0102191);
   MVA_BDT_S__38->SetBinContent(18,0.01266927);
   MVA_BDT_S__38->SetBinContent(19,0.01694198);
   MVA_BDT_S__38->SetBinContent(20,0.03247984);
   MVA_BDT_S__38->SetBinContent(21,0.07042728);
   MVA_BDT_S__38->SetBinContent(22,0.06636379);
   MVA_BDT_S__38->SetBinContent(23,0.1457549);
   MVA_BDT_S__38->SetBinContent(24,0.2187531);
   MVA_BDT_S__38->SetBinContent(25,0.3762519);
   MVA_BDT_S__38->SetBinContent(26,0.5696061);
   MVA_BDT_S__38->SetBinContent(27,0.9626189);
   MVA_BDT_S__38->SetBinContent(28,1.556843);
   MVA_BDT_S__38->SetBinContent(29,2.257354);
   MVA_BDT_S__38->SetBinContent(30,3.422497);
   MVA_BDT_S__38->SetBinContent(31,4.996547);
   MVA_BDT_S__38->SetBinContent(32,5.829738);
   MVA_BDT_S__38->SetBinContent(33,7.190499);
   MVA_BDT_S__38->SetBinContent(34,7.254296);
   MVA_BDT_S__38->SetBinContent(35,5.936021);
   MVA_BDT_S__38->SetBinContent(36,4.140517);
   MVA_BDT_S__38->SetBinContent(37,2.041832);
   MVA_BDT_S__38->SetBinContent(38,0.7474506);
   MVA_BDT_S__38->SetBinContent(39,0.1851668);
   MVA_BDT_S__38->SetBinContent(40,0.02127504);
   MVA_BDT_S__38->SetBinError(9,0.0002092104);
   MVA_BDT_S__38->SetBinError(13,0.0002092104);
   MVA_BDT_S__38->SetBinError(15,0.0002092104);
   MVA_BDT_S__38->SetBinError(16,0.003182976);
   MVA_BDT_S__38->SetBinError(17,0.00451112);
   MVA_BDT_S__38->SetBinError(18,0.005041416);
   MVA_BDT_S__38->SetBinError(19,0.005951126);
   MVA_BDT_S__38->SetBinError(20,0.008123114);
   MVA_BDT_S__38->SetBinError(21,0.01211499);
   MVA_BDT_S__38->SetBinError(22,0.01169687);
   MVA_BDT_S__38->SetBinError(23,0.01742634);
   MVA_BDT_S__38->SetBinError(24,0.02112901);
   MVA_BDT_S__38->SetBinError(25,0.02768084);
   MVA_BDT_S__38->SetBinError(26,0.0341522);
   MVA_BDT_S__38->SetBinError(27,0.0445156);
   MVA_BDT_S__38->SetBinError(28,0.05675789);
   MVA_BDT_S__38->SetBinError(29,0.06823877);
   MVA_BDT_S__38->SetBinError(30,0.084204);
   MVA_BDT_S__38->SetBinError(31,0.1018556);
   MVA_BDT_S__38->SetBinError(32,0.1095089);
   MVA_BDT_S__38->SetBinError(33,0.12197);
   MVA_BDT_S__38->SetBinError(34,0.1224159);
   MVA_BDT_S__38->SetBinError(35,0.110541);
   MVA_BDT_S__38->SetBinError(36,0.09254984);
   MVA_BDT_S__38->SetBinError(37,0.0647148);
   MVA_BDT_S__38->SetBinError(38,0.03936512);
   MVA_BDT_S__38->SetBinError(39,0.0196156);
   MVA_BDT_S__38->SetBinError(40,0.006393394);
   MVA_BDT_S__38->SetEntries(39509);
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
   
   TH1D *MVA_BDT_B__39 = new TH1D("MVA_BDT_B__39","MVA_BDT_B",40,-0.519245,0.3129154);
   MVA_BDT_B__39->SetBinContent(1,0.008758894);
   MVA_BDT_B__39->SetBinContent(2,0.01512685);
   MVA_BDT_B__39->SetBinContent(3,0.02001091);
   MVA_BDT_B__39->SetBinContent(4,0.08337004);
   MVA_BDT_B__39->SetBinContent(5,0.1485388);
   MVA_BDT_B__39->SetBinContent(6,0.3660914);
   MVA_BDT_B__39->SetBinContent(7,0.6502555);
   MVA_BDT_B__39->SetBinContent(8,1.217648);
   MVA_BDT_B__39->SetBinContent(9,1.852001);
   MVA_BDT_B__39->SetBinContent(10,2.572796);
   MVA_BDT_B__39->SetBinContent(11,3.277275);
   MVA_BDT_B__39->SetBinContent(12,3.822311);
   MVA_BDT_B__39->SetBinContent(13,4.037774);
   MVA_BDT_B__39->SetBinContent(14,4.222814);
   MVA_BDT_B__39->SetBinContent(15,4.111372);
   MVA_BDT_B__39->SetBinContent(16,3.797098);
   MVA_BDT_B__39->SetBinContent(17,3.566804);
   MVA_BDT_B__39->SetBinContent(18,3.104512);
   MVA_BDT_B__39->SetBinContent(19,2.619362);
   MVA_BDT_B__39->SetBinContent(20,2.103884);
   MVA_BDT_B__39->SetBinContent(21,1.538654);
   MVA_BDT_B__39->SetBinContent(22,1.238001);
   MVA_BDT_B__39->SetBinContent(23,1.000631);
   MVA_BDT_B__39->SetBinContent(24,0.7687016);
   MVA_BDT_B__39->SetBinContent(25,0.580912);
   MVA_BDT_B__39->SetBinContent(26,0.4880288);
   MVA_BDT_B__39->SetBinContent(27,0.3254275);
   MVA_BDT_B__39->SetBinContent(28,0.1550524);
   MVA_BDT_B__39->SetBinContent(29,0.1374818);
   MVA_BDT_B__39->SetBinContent(30,0.08838106);
   MVA_BDT_B__39->SetBinContent(31,0.04989834);
   MVA_BDT_B__39->SetBinContent(32,0.04929857);
   MVA_BDT_B__39->SetBinContent(33,0.02058235);
   MVA_BDT_B__39->SetBinContent(34,0.01941162);
   MVA_BDT_B__39->SetBinContent(35,0.00938619);
   MVA_BDT_B__39->SetBinError(1,0.006059996);
   MVA_BDT_B__39->SetBinError(2,0.007546014);
   MVA_BDT_B__39->SetBinError(3,0.007743416);
   MVA_BDT_B__39->SetBinError(4,0.01648767);
   MVA_BDT_B__39->SetBinError(5,0.02075138);
   MVA_BDT_B__39->SetBinError(6,0.03415265);
   MVA_BDT_B__39->SetBinError(7,0.04448755);
   MVA_BDT_B__39->SetBinError(8,0.06000472);
   MVA_BDT_B__39->SetBinError(9,0.07593451);
   MVA_BDT_B__39->SetBinError(10,0.08927099);
   MVA_BDT_B__39->SetBinError(11,0.1012237);
   MVA_BDT_B__39->SetBinError(12,0.1093679);
   MVA_BDT_B__39->SetBinError(13,0.1125281);
   MVA_BDT_B__39->SetBinError(14,0.1163432);
   MVA_BDT_B__39->SetBinError(15,0.115223);
   MVA_BDT_B__39->SetBinError(16,0.1106246);
   MVA_BDT_B__39->SetBinError(17,0.1067291);
   MVA_BDT_B__39->SetBinError(18,0.09938381);
   MVA_BDT_B__39->SetBinError(19,0.09088659);
   MVA_BDT_B__39->SetBinError(20,0.08105485);
   MVA_BDT_B__39->SetBinError(21,0.06819317);
   MVA_BDT_B__39->SetBinError(22,0.06192789);
   MVA_BDT_B__39->SetBinError(23,0.05729984);
   MVA_BDT_B__39->SetBinError(24,0.04961853);
   MVA_BDT_B__39->SetBinError(25,0.04310604);
   MVA_BDT_B__39->SetBinError(26,0.04141048);
   MVA_BDT_B__39->SetBinError(27,0.03408716);
   MVA_BDT_B__39->SetBinError(28,0.02284538);
   MVA_BDT_B__39->SetBinError(29,0.0215432);
   MVA_BDT_B__39->SetBinError(30,0.01793908);
   MVA_BDT_B__39->SetBinError(31,0.01254902);
   MVA_BDT_B__39->SetBinError(32,0.01312813);
   MVA_BDT_B__39->SetBinError(33,0.008640734);
   MVA_BDT_B__39->SetBinError(34,0.008616986);
   MVA_BDT_B__39->SetBinError(35,0.006084636);
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
   
   TH1D *MVA_BDT_Train_S__40 = new TH1D("MVA_BDT_Train_S__40","MVA_BDT_Train_S",40,-0.519245,0.3129154);
   MVA_BDT_Train_S__40->SetBinContent(16,0.002238518);
   MVA_BDT_Train_S__40->SetBinContent(17,0.004477037);
   MVA_BDT_Train_S__40->SetBinContent(18,0.004477037);
   MVA_BDT_Train_S__40->SetBinContent(19,0.01041695);
   MVA_BDT_Train_S__40->SetBinContent(20,0.01895306);
   MVA_BDT_Train_S__40->SetBinContent(21,0.03181751);
   MVA_BDT_Train_S__40->SetBinContent(22,0.06208375);
   MVA_BDT_Train_S__40->SetBinContent(23,0.0884396);
   MVA_BDT_Train_S__40->SetBinContent(24,0.1521952);
   MVA_BDT_Train_S__40->SetBinContent(25,0.3566198);
   MVA_BDT_Train_S__40->SetBinContent(26,0.6463493);
   MVA_BDT_Train_S__40->SetBinContent(27,1.009146);
   MVA_BDT_Train_S__40->SetBinContent(28,1.555835);
   MVA_BDT_Train_S__40->SetBinContent(29,2.351273);
   MVA_BDT_Train_S__40->SetBinContent(30,3.358984);
   MVA_BDT_Train_S__40->SetBinContent(31,4.707232);
   MVA_BDT_Train_S__40->SetBinContent(32,6.176777);
   MVA_BDT_Train_S__40->SetBinContent(33,7.104947);
   MVA_BDT_Train_S__40->SetBinContent(34,7.254007);
   MVA_BDT_Train_S__40->SetBinContent(35,6.074376);
   MVA_BDT_Train_S__40->SetBinContent(36,4.092045);
   MVA_BDT_Train_S__40->SetBinContent(37,2.109682);
   MVA_BDT_Train_S__40->SetBinContent(38,0.6779297);
   MVA_BDT_Train_S__40->SetBinContent(39,0.1823087);
   MVA_BDT_Train_S__40->SetBinContent(40,0.03504066);
   MVA_BDT_Train_S__40->SetBinContent(41,0.004477037);
   MVA_BDT_Train_S__40->SetBinError(16,0.002238518);
   MVA_BDT_Train_S__40->SetBinError(17,0.003165743);
   MVA_BDT_Train_S__40->SetBinError(18,0.003165743);
   MVA_BDT_Train_S__40->SetBinError(19,0.00451105);
   MVA_BDT_Train_S__40->SetBinError(20,0.006348708);
   MVA_BDT_Train_S__40->SetBinError(21,0.008106189);
   MVA_BDT_Train_S__40->SetBinError(22,0.01104582);
   MVA_BDT_Train_S__40->SetBinError(23,0.01315103);
   MVA_BDT_Train_S__40->SetBinError(24,0.01718456);
   MVA_BDT_Train_S__40->SetBinError(25,0.02699521);
   MVA_BDT_Train_S__40->SetBinError(26,0.03653037);
   MVA_BDT_Train_S__40->SetBinError(27,0.04568606);
   MVA_BDT_Train_S__40->SetBinError(28,0.05665726);
   MVA_BDT_Train_S__40->SetBinError(29,0.06966206);
   MVA_BDT_Train_S__40->SetBinError(30,0.08338674);
   MVA_BDT_Train_S__40->SetBinError(31,0.0985681);
   MVA_BDT_Train_S__40->SetBinError(32,0.11291);
   MVA_BDT_Train_S__40->SetBinError(33,0.1211158);
   MVA_BDT_Train_S__40->SetBinError(34,0.1223262);
   MVA_BDT_Train_S__40->SetBinError(35,0.1117633);
   MVA_BDT_Train_S__40->SetBinError(36,0.09183216);
   MVA_BDT_Train_S__40->SetBinError(37,0.06628296);
   MVA_BDT_Train_S__40->SetBinError(38,0.03723573);
   MVA_BDT_Train_S__40->SetBinError(39,0.01946371);
   MVA_BDT_Train_S__40->SetBinError(40,0.008687358);
   MVA_BDT_Train_S__40->SetBinError(41,0.003165743);
   MVA_BDT_Train_S__40->SetEntries(39509);

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
   
   TH1D *MVA_BDT_Train_B__41 = new TH1D("MVA_BDT_Train_B__41","MVA_BDT_Train_B",40,-0.519245,0.3129154);
   MVA_BDT_Train_B__41->SetBinContent(1,0.0001960227);
   MVA_BDT_Train_B__41->SetBinContent(2,0.01058853);
   MVA_BDT_Train_B__41->SetBinContent(3,0.00826692);
   MVA_BDT_Train_B__41->SetBinContent(4,0.0555104);
   MVA_BDT_Train_B__41->SetBinContent(5,0.1221465);
   MVA_BDT_Train_B__41->SetBinContent(6,0.3473633);
   MVA_BDT_Train_B__41->SetBinContent(7,0.7227802);
   MVA_BDT_Train_B__41->SetBinContent(8,1.211794);
   MVA_BDT_Train_B__41->SetBinContent(9,1.897898);
   MVA_BDT_Train_B__41->SetBinContent(10,2.677468);
   MVA_BDT_Train_B__41->SetBinContent(11,3.166384);
   MVA_BDT_Train_B__41->SetBinContent(12,3.556172);
   MVA_BDT_Train_B__41->SetBinContent(13,4.087859);
   MVA_BDT_Train_B__41->SetBinContent(14,4.165871);
   MVA_BDT_Train_B__41->SetBinContent(15,4.123839);
   MVA_BDT_Train_B__41->SetBinContent(16,4.007436);
   MVA_BDT_Train_B__41->SetBinContent(17,3.326564);
   MVA_BDT_Train_B__41->SetBinContent(18,3.175127);
   MVA_BDT_Train_B__41->SetBinContent(19,2.673962);
   MVA_BDT_Train_B__41->SetBinContent(20,2.256964);
   MVA_BDT_Train_B__41->SetBinContent(21,1.654034);
   MVA_BDT_Train_B__41->SetBinContent(22,1.444713);
   MVA_BDT_Train_B__41->SetBinContent(23,1.063076);
   MVA_BDT_Train_B__41->SetBinContent(24,0.7961654);
   MVA_BDT_Train_B__41->SetBinContent(25,0.5471968);
   MVA_BDT_Train_B__41->SetBinContent(26,0.3823618);
   MVA_BDT_Train_B__41->SetBinContent(27,0.2239597);
   MVA_BDT_Train_B__41->SetBinContent(28,0.1140172);
   MVA_BDT_Train_B__41->SetBinContent(29,0.09117783);
   MVA_BDT_Train_B__41->SetBinContent(30,0.06089109);
   MVA_BDT_Train_B__41->SetBinContent(31,0.03083451);
   MVA_BDT_Train_B__41->SetBinContent(32,0.04971112);
   MVA_BDT_Train_B__41->SetBinContent(33,0.01281231);
   MVA_BDT_Train_B__41->SetBinContent(34,0.0006488559);
   MVA_BDT_Train_B__41->SetBinContent(35,0.001861531);
   MVA_BDT_Train_B__41->SetBinError(1,0.0001960227);
   MVA_BDT_Train_B__41->SetBinError(2,0.00625705);
   MVA_BDT_Train_B__41->SetBinError(3,0.002499574);
   MVA_BDT_Train_B__41->SetBinError(4,0.0123332);
   MVA_BDT_Train_B__41->SetBinError(5,0.01886144);
   MVA_BDT_Train_B__41->SetBinError(6,0.03256097);
   MVA_BDT_Train_B__41->SetBinError(7,0.04809985);
   MVA_BDT_Train_B__41->SetBinError(8,0.06065404);
   MVA_BDT_Train_B__41->SetBinError(9,0.0767122);
   MVA_BDT_Train_B__41->SetBinError(10,0.09066964);
   MVA_BDT_Train_B__41->SetBinError(11,0.09927161);
   MVA_BDT_Train_B__41->SetBinError(12,0.1052647);
   MVA_BDT_Train_B__41->SetBinError(13,0.1149482);
   MVA_BDT_Train_B__41->SetBinError(14,0.1157548);
   MVA_BDT_Train_B__41->SetBinError(15,0.1155451);
   MVA_BDT_Train_B__41->SetBinError(16,0.1137695);
   MVA_BDT_Train_B__41->SetBinError(17,0.1033202);
   MVA_BDT_Train_B__41->SetBinError(18,0.100539);
   MVA_BDT_Train_B__41->SetBinError(19,0.09257018);
   MVA_BDT_Train_B__41->SetBinError(20,0.08550475);
   MVA_BDT_Train_B__41->SetBinError(21,0.0718634);
   MVA_BDT_Train_B__41->SetBinError(22,0.06879163);
   MVA_BDT_Train_B__41->SetBinError(23,0.05888041);
   MVA_BDT_Train_B__41->SetBinError(24,0.05105937);
   MVA_BDT_Train_B__41->SetBinError(25,0.04172374);
   MVA_BDT_Train_B__41->SetBinError(26,0.03582521);
   MVA_BDT_Train_B__41->SetBinError(27,0.0267082);
   MVA_BDT_Train_B__41->SetBinError(28,0.01821769);
   MVA_BDT_Train_B__41->SetBinError(29,0.01730038);
   MVA_BDT_Train_B__41->SetBinError(30,0.01417152);
   MVA_BDT_Train_B__41->SetBinError(31,0.01007665);
   MVA_BDT_Train_B__41->SetBinError(32,0.01392665);
   MVA_BDT_Train_B__41->SetBinError(33,0.006360604);
   MVA_BDT_Train_B__41->SetBinError(34,0.0004616533);
   MVA_BDT_Train_B__41->SetBinError(35,0.0008555139);
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
   TText *text = new TText(0.12,0.74,"Kolmogorov-Smirnov test: signal (background) probability = 0.268 (0.391)");
   text->SetNDC();
   text->SetTextSize(0.032);
   text->Draw();
   
   TH2F *frameBDT__42 = new TH2F("frameBDT__42","TMVA overtraining check for classifier: BDT",500,-0.519245,0.3129154,500,0,9.430585);
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
