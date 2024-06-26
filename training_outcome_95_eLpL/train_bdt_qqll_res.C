#ifdef __CLING__
#pragma cling optimize(0)
#endif
void train_bdt_qqll_res()
{
//=========Macro generated from canvas: canvas1/TMVA comparison BDT
//=========  (Mon Jul  1 17:54:54 2024) by ROOT version 6.32.02
   TCanvas *canvas1 = new TCanvas("canvas1", "TMVA comparison BDT",200,81,600,468);
   gStyle->SetOptStat(0);
   canvas1->Range(-0.6118295,-1.744272,0.3120082,12.79133);
   canvas1->SetFillColor(0);
   canvas1->SetBorderMode(0);
   canvas1->SetBorderSize(2);
   canvas1->SetTickx(1);
   canvas1->SetTicky(1);
   canvas1->SetRightMargin(0.05);
   canvas1->SetBottomMargin(0.12);
   canvas1->SetFrameBorderMode(0);
   canvas1->SetFrameBorderMode(0);
   
   TH2F *frameBDT__37 = new TH2F("frameBDT__37","TMVA overtraining check for classifier: BDT",500,-0.5148265,0.2658163,500,0,11.33777);
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
   
   TH1D *MVA_BDT_S__38 = new TH1D("MVA_BDT_S__38","TMVA overtraining check for classifier: BDT",40,-0.5148265,0.2658163);
   MVA_BDT_S__38->SetBinContent(10,0.003033742);
   MVA_BDT_S__38->SetBinContent(13,0.003033742);
   MVA_BDT_S__38->SetBinContent(16,0.00160713);
   MVA_BDT_S__38->SetBinContent(17,0.006469265);
   MVA_BDT_S__38->SetBinContent(18,0.01414388);
   MVA_BDT_S__38->SetBinContent(19,0.0351588);
   MVA_BDT_S__38->SetBinContent(20,0.06951404);
   MVA_BDT_S__38->SetBinContent(21,0.1154219);
   MVA_BDT_S__38->SetBinContent(22,0.1578536);
   MVA_BDT_S__38->SetBinContent(23,0.2473577);
   MVA_BDT_S__38->SetBinContent(24,0.2950533);
   MVA_BDT_S__38->SetBinContent(25,0.5405826);
   MVA_BDT_S__38->SetBinContent(26,0.9781402);
   MVA_BDT_S__38->SetBinContent(27,1.637901);
   MVA_BDT_S__38->SetBinContent(28,2.834488);
   MVA_BDT_S__38->SetBinContent(29,4.540138);
   MVA_BDT_S__38->SetBinContent(30,6.604323);
   MVA_BDT_S__38->SetBinContent(31,8.050877);
   MVA_BDT_S__38->SetBinContent(32,8.367108);
   MVA_BDT_S__38->SetBinContent(33,7.398471);
   MVA_BDT_S__38->SetBinContent(34,4.943736);
   MVA_BDT_S__38->SetBinContent(35,2.727381);
   MVA_BDT_S__38->SetBinContent(36,1.168719);
   MVA_BDT_S__38->SetBinContent(37,0.4000872);
   MVA_BDT_S__38->SetBinContent(38,0.07759044);
   MVA_BDT_S__38->SetBinContent(39,0.01820245);
   MVA_BDT_S__38->SetBinContent(40,0.003435524);
   MVA_BDT_S__38->SetBinError(10,0.003033742);
   MVA_BDT_S__38->SetBinError(13,0.003033742);
   MVA_BDT_S__38->SetBinError(16,0.0008035648);
   MVA_BDT_S__38->SetBinError(17,0.00430913);
   MVA_BDT_S__38->SetBinError(18,0.006133636);
   MVA_BDT_S__38->SetBinError(19,0.009693968);
   MVA_BDT_S__38->SetBinError(20,0.01369756);
   MVA_BDT_S__38->SetBinError(21,0.018051);
   MVA_BDT_S__38->SetBinError(22,0.02094524);
   MVA_BDT_S__38->SetBinError(23,0.02596693);
   MVA_BDT_S__38->SetBinError(24,0.02837657);
   MVA_BDT_S__38->SetBinError(25,0.0383509);
   MVA_BDT_S__38->SetBinError(26,0.05173605);
   MVA_BDT_S__38->SetBinError(27,0.06671435);
   MVA_BDT_S__38->SetBinError(28,0.08785344);
   MVA_BDT_S__38->SetBinError(29,0.1115875);
   MVA_BDT_S__38->SetBinError(30,0.1342238);
   MVA_BDT_S__38->SetBinError(31,0.1479153);
   MVA_BDT_S__38->SetBinError(32,0.1510115);
   MVA_BDT_S__38->SetBinError(33,0.1419705);
   MVA_BDT_S__38->SetBinError(34,0.1157774);
   MVA_BDT_S__38->SetBinError(35,0.08571943);
   MVA_BDT_S__38->SetBinError(36,0.05586131);
   MVA_BDT_S__38->SetBinError(37,0.03272617);
   MVA_BDT_S__38->SetBinError(38,0.01438185);
   MVA_BDT_S__38->SetBinError(39,0.007431119);
   MVA_BDT_S__38->SetBinError(40,0.003060231);
   MVA_BDT_S__38->SetEntries(30024);
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
   
   TH1D *MVA_BDT_B__39 = new TH1D("MVA_BDT_B__39","MVA_BDT_B",40,-0.5148265,0.2658163);
   MVA_BDT_B__39->SetBinContent(1,0.00759651);
   MVA_BDT_B__39->SetBinContent(2,0.03175649);
   MVA_BDT_B__39->SetBinContent(3,0.09278035);
   MVA_BDT_B__39->SetBinContent(4,0.2931268);
   MVA_BDT_B__39->SetBinContent(5,0.6768314);
   MVA_BDT_B__39->SetBinContent(6,1.189052);
   MVA_BDT_B__39->SetBinContent(7,1.666551);
   MVA_BDT_B__39->SetBinContent(8,1.924449);
   MVA_BDT_B__39->SetBinContent(9,2.230564);
   MVA_BDT_B__39->SetBinContent(10,2.399662);
   MVA_BDT_B__39->SetBinContent(11,2.541689);
   MVA_BDT_B__39->SetBinContent(12,2.387379);
   MVA_BDT_B__39->SetBinContent(13,2.221335);
   MVA_BDT_B__39->SetBinContent(14,2.265918);
   MVA_BDT_B__39->SetBinContent(15,2.338154);
   MVA_BDT_B__39->SetBinContent(16,2.317704);
   MVA_BDT_B__39->SetBinContent(17,2.377502);
   MVA_BDT_B__39->SetBinContent(18,2.644171);
   MVA_BDT_B__39->SetBinContent(19,2.692616);
   MVA_BDT_B__39->SetBinContent(20,2.799291);
   MVA_BDT_B__39->SetBinContent(21,2.297945);
   MVA_BDT_B__39->SetBinContent(22,2.005409);
   MVA_BDT_B__39->SetBinContent(23,1.779454);
   MVA_BDT_B__39->SetBinContent(24,1.43212);
   MVA_BDT_B__39->SetBinContent(25,1.292507);
   MVA_BDT_B__39->SetBinContent(26,1.142049);
   MVA_BDT_B__39->SetBinContent(27,1.038815);
   MVA_BDT_B__39->SetBinContent(28,1.138762);
   MVA_BDT_B__39->SetBinContent(29,1.044943);
   MVA_BDT_B__39->SetBinContent(30,1.005247);
   MVA_BDT_B__39->SetBinContent(31,0.8173136);
   MVA_BDT_B__39->SetBinContent(32,0.5298321);
   MVA_BDT_B__39->SetBinContent(33,0.3109597);
   MVA_BDT_B__39->SetBinContent(34,0.1784277);
   MVA_BDT_B__39->SetBinContent(35,0.07548898);
   MVA_BDT_B__39->SetBinContent(36,0.04112325);
   MVA_BDT_B__39->SetBinContent(37,0.01130007);
   MVA_BDT_B__39->SetBinError(1,0.004215413);
   MVA_BDT_B__39->SetBinError(2,0.007497198);
   MVA_BDT_B__39->SetBinError(3,0.01187819);
   MVA_BDT_B__39->SetBinError(4,0.02116892);
   MVA_BDT_B__39->SetBinError(5,0.03363985);
   MVA_BDT_B__39->SetBinError(6,0.04490185);
   MVA_BDT_B__39->SetBinError(7,0.056285);
   MVA_BDT_B__39->SetBinError(8,0.06286332);
   MVA_BDT_B__39->SetBinError(9,0.07002392);
   MVA_BDT_B__39->SetBinError(10,0.07388435);
   MVA_BDT_B__39->SetBinError(11,0.07916583);
   MVA_BDT_B__39->SetBinError(12,0.0771588);
   MVA_BDT_B__39->SetBinError(13,0.07490551);
   MVA_BDT_B__39->SetBinError(14,0.07823252);
   MVA_BDT_B__39->SetBinError(15,0.08003091);
   MVA_BDT_B__39->SetBinError(16,0.07945631);
   MVA_BDT_B__39->SetBinError(17,0.08176874);
   MVA_BDT_B__39->SetBinError(18,0.08843795);
   MVA_BDT_B__39->SetBinError(19,0.08895164);
   MVA_BDT_B__39->SetBinError(20,0.09225783);
   MVA_BDT_B__39->SetBinError(21,0.082816);
   MVA_BDT_B__39->SetBinError(22,0.07789666);
   MVA_BDT_B__39->SetBinError(23,0.0741691);
   MVA_BDT_B__39->SetBinError(24,0.06586439);
   MVA_BDT_B__39->SetBinError(25,0.06256081);
   MVA_BDT_B__39->SetBinError(26,0.05910435);
   MVA_BDT_B__39->SetBinError(27,0.05615114);
   MVA_BDT_B__39->SetBinError(28,0.05935776);
   MVA_BDT_B__39->SetBinError(29,0.05740085);
   MVA_BDT_B__39->SetBinError(30,0.05694681);
   MVA_BDT_B__39->SetBinError(31,0.05163835);
   MVA_BDT_B__39->SetBinError(32,0.04162006);
   MVA_BDT_B__39->SetBinError(33,0.03184933);
   MVA_BDT_B__39->SetBinError(34,0.02424073);
   MVA_BDT_B__39->SetBinError(35,0.01619809);
   MVA_BDT_B__39->SetBinError(36,0.01200515);
   MVA_BDT_B__39->SetBinError(37,0.006272629);
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
   
   TH1D *MVA_BDT_Train_S__40 = new TH1D("MVA_BDT_Train_S__40","MVA_BDT_Train_S",40,-0.5148265,0.2658163);
   MVA_BDT_Train_S__40->SetBinContent(17,0.003422865);
   MVA_BDT_Train_S__40->SetBinContent(18,0.006445428);
   MVA_BDT_Train_S__40->SetBinContent(19,0.02195855);
   MVA_BDT_Train_S__40->SetBinContent(20,0.0348088);
   MVA_BDT_Train_S__40->SetBinContent(21,0.03641001);
   MVA_BDT_Train_S__40->SetBinContent(22,0.08655208);
   MVA_BDT_Train_S__40->SetBinContent(23,0.1555895);
   MVA_BDT_Train_S__40->SetBinContent(24,0.2564538);
   MVA_BDT_Train_S__40->SetBinContent(25,0.4868474);
   MVA_BDT_Train_S__40->SetBinContent(26,0.798247);
   MVA_BDT_Train_S__40->SetBinContent(27,1.472812);
   MVA_BDT_Train_S__40->SetBinContent(28,2.666249);
   MVA_BDT_Train_S__40->SetBinContent(29,4.316691);
   MVA_BDT_Train_S__40->SetBinContent(30,6.329197);
   MVA_BDT_Train_S__40->SetBinContent(31,8.302931);
   MVA_BDT_Train_S__40->SetBinContent(32,8.721362);
   MVA_BDT_Train_S__40->SetBinContent(33,7.664259);
   MVA_BDT_Train_S__40->SetBinContent(34,5.154933);
   MVA_BDT_Train_S__40->SetBinContent(35,2.90397);
   MVA_BDT_Train_S__40->SetBinContent(36,1.318905);
   MVA_BDT_Train_S__40->SetBinContent(37,0.3785167);
   MVA_BDT_Train_S__40->SetBinContent(38,0.109172);
   MVA_BDT_Train_S__40->SetBinContent(39,0.0110692);
   MVA_BDT_Train_S__40->SetBinContent(40,0.003022563);
   MVA_BDT_Train_S__40->SetBinError(17,0.003048955);
   MVA_BDT_Train_S__40->SetBinError(18,0.004293253);
   MVA_BDT_Train_S__40->SetBinError(19,0.008016963);
   MVA_BDT_Train_S__40->SetBinError(20,0.009234045);
   MVA_BDT_Train_S__40->SetBinError(21,0.009268686);
   MVA_BDT_Train_S__40->SetBinError(22,0.01499566);
   MVA_BDT_Train_S__40->SetBinError(23,0.02007732);
   MVA_BDT_Train_S__40->SetBinError(24,0.02594855);
   MVA_BDT_Train_S__40->SetBinError(25,0.03633691);
   MVA_BDT_Train_S__40->SetBinError(26,0.04620154);
   MVA_BDT_Train_S__40->SetBinError(27,0.06316488);
   MVA_BDT_Train_S__40->SetBinError(28,0.08472433);
   MVA_BDT_Train_S__40->SetBinError(29,0.1080144);
   MVA_BDT_Train_S__40->SetBinError(30,0.1309846);
   MVA_BDT_Train_S__40->SetBinError(31,0.1502525);
   MVA_BDT_Train_S__40->SetBinError(32,0.1540329);
   MVA_BDT_Train_S__40->SetBinError(33,0.144628);
   MVA_BDT_Train_S__40->SetBinError(34,0.1181541);
   MVA_BDT_Train_S__40->SetBinError(35,0.08874539);
   MVA_BDT_Train_S__40->SetBinError(36,0.0601588);
   MVA_BDT_Train_S__40->SetBinError(37,0.03199024);
   MVA_BDT_Train_S__40->SetBinError(38,0.01769701);
   MVA_BDT_Train_S__40->SetBinError(39,0.005311202);
   MVA_BDT_Train_S__40->SetBinError(40,0.003022563);
   MVA_BDT_Train_S__40->SetEntries(30024);

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
   
   TH1D *MVA_BDT_Train_B__41 = new TH1D("MVA_BDT_Train_B__41","MVA_BDT_Train_B",40,-0.5148265,0.2658163);
   MVA_BDT_Train_B__41->SetBinContent(0,0.001702635);
   MVA_BDT_Train_B__41->SetBinContent(1,0.005655927);
   MVA_BDT_Train_B__41->SetBinContent(2,0.02229555);
   MVA_BDT_Train_B__41->SetBinContent(3,0.0974957);
   MVA_BDT_Train_B__41->SetBinContent(4,0.2599269);
   MVA_BDT_Train_B__41->SetBinContent(5,0.6523859);
   MVA_BDT_Train_B__41->SetBinContent(6,1.187915);
   MVA_BDT_Train_B__41->SetBinContent(7,1.615416);
   MVA_BDT_Train_B__41->SetBinContent(8,1.871548);
   MVA_BDT_Train_B__41->SetBinContent(9,2.266301);
   MVA_BDT_Train_B__41->SetBinContent(10,2.495942);
   MVA_BDT_Train_B__41->SetBinContent(11,2.437532);
   MVA_BDT_Train_B__41->SetBinContent(12,2.353471);
   MVA_BDT_Train_B__41->SetBinContent(13,2.320585);
   MVA_BDT_Train_B__41->SetBinContent(14,2.226484);
   MVA_BDT_Train_B__41->SetBinContent(15,2.288754);
   MVA_BDT_Train_B__41->SetBinContent(16,2.230167);
   MVA_BDT_Train_B__41->SetBinContent(17,2.583837);
   MVA_BDT_Train_B__41->SetBinContent(18,2.587826);
   MVA_BDT_Train_B__41->SetBinContent(19,2.720834);
   MVA_BDT_Train_B__41->SetBinContent(20,2.586463);
   MVA_BDT_Train_B__41->SetBinContent(21,2.557789);
   MVA_BDT_Train_B__41->SetBinContent(22,2.11901);
   MVA_BDT_Train_B__41->SetBinContent(23,1.772269);
   MVA_BDT_Train_B__41->SetBinContent(24,1.575336);
   MVA_BDT_Train_B__41->SetBinContent(25,1.304033);
   MVA_BDT_Train_B__41->SetBinContent(26,1.320577);
   MVA_BDT_Train_B__41->SetBinContent(27,1.201693);
   MVA_BDT_Train_B__41->SetBinContent(28,1.210159);
   MVA_BDT_Train_B__41->SetBinContent(29,1.21584);
   MVA_BDT_Train_B__41->SetBinContent(30,0.8849378);
   MVA_BDT_Train_B__41->SetBinContent(31,0.5376621);
   MVA_BDT_Train_B__41->SetBinContent(32,0.359537);
   MVA_BDT_Train_B__41->SetBinContent(33,0.2092623);
   MVA_BDT_Train_B__41->SetBinContent(34,0.0930372);
   MVA_BDT_Train_B__41->SetBinContent(35,0.04113197);
   MVA_BDT_Train_B__41->SetBinContent(36,0.0190154);
   MVA_BDT_Train_B__41->SetBinContent(37,0.007699397);
   MVA_BDT_Train_B__41->SetBinError(0,0.001702635);
   MVA_BDT_Train_B__41->SetBinError(1,0.002552008);
   MVA_BDT_Train_B__41->SetBinError(2,0.005862425);
   MVA_BDT_Train_B__41->SetBinError(3,0.01312006);
   MVA_BDT_Train_B__41->SetBinError(4,0.01957553);
   MVA_BDT_Train_B__41->SetBinError(5,0.03228912);
   MVA_BDT_Train_B__41->SetBinError(6,0.04498424);
   MVA_BDT_Train_B__41->SetBinError(7,0.05462982);
   MVA_BDT_Train_B__41->SetBinError(8,0.06105645);
   MVA_BDT_Train_B__41->SetBinError(9,0.07110007);
   MVA_BDT_Train_B__41->SetBinError(10,0.07679791);
   MVA_BDT_Train_B__41->SetBinError(11,0.07655279);
   MVA_BDT_Train_B__41->SetBinError(12,0.07658346);
   MVA_BDT_Train_B__41->SetBinError(13,0.0776748);
   MVA_BDT_Train_B__41->SetBinError(14,0.07663529);
   MVA_BDT_Train_B__41->SetBinError(15,0.07894439);
   MVA_BDT_Train_B__41->SetBinError(16,0.07829946);
   MVA_BDT_Train_B__41->SetBinError(17,0.08622973);
   MVA_BDT_Train_B__41->SetBinError(18,0.08701716);
   MVA_BDT_Train_B__41->SetBinError(19,0.0896096);
   MVA_BDT_Train_B__41->SetBinError(20,0.08814174);
   MVA_BDT_Train_B__41->SetBinError(21,0.0883109);
   MVA_BDT_Train_B__41->SetBinError(22,0.08034161);
   MVA_BDT_Train_B__41->SetBinError(23,0.07331576);
   MVA_BDT_Train_B__41->SetBinError(24,0.06985873);
   MVA_BDT_Train_B__41->SetBinError(25,0.06323182);
   MVA_BDT_Train_B__41->SetBinError(26,0.06433048);
   MVA_BDT_Train_B__41->SetBinError(27,0.06149365);
   MVA_BDT_Train_B__41->SetBinError(28,0.06191583);
   MVA_BDT_Train_B__41->SetBinError(29,0.06272584);
   MVA_BDT_Train_B__41->SetBinError(30,0.053272);
   MVA_BDT_Train_B__41->SetBinError(31,0.04075529);
   MVA_BDT_Train_B__41->SetBinError(32,0.03326192);
   MVA_BDT_Train_B__41->SetBinError(33,0.02555904);
   MVA_BDT_Train_B__41->SetBinError(34,0.01720762);
   MVA_BDT_Train_B__41->SetBinError(35,0.01162837);
   MVA_BDT_Train_B__41->SetBinError(36,0.008113814);
   MVA_BDT_Train_B__41->SetBinError(37,0.00513586);
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
   TText *text = new TText(0.12,0.74,"Kolmogorov-Smirnov test: signal (background) probability =     0 (0.009)");
   text->SetNDC();
   text->SetTextSize(0.032);
   text->Draw();
   
   TH2F *frameBDT__42 = new TH2F("frameBDT__42","TMVA overtraining check for classifier: BDT",500,-0.5148265,0.2658163,500,0,11.33777);
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
