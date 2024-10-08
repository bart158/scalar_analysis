#ifdef __CLING__
#pragma cling optimize(0)
#endif
void train_bdt_qqll_res()
{
//=========Macro generated from canvas: canvas1/TMVA comparison BDT
//=========  (Wed Aug 14 12:55:21 2024) by ROOT version 6.32.02
   TCanvas *canvas1 = new TCanvas("canvas1", "TMVA comparison BDT",200,81,600,468);
   gStyle->SetOptStat(0);
   canvas1->Range(-0.6629098,-1.389935,0.3566272,10.19285);
   canvas1->SetFillColor(0);
   canvas1->SetBorderMode(0);
   canvas1->SetBorderSize(2);
   canvas1->SetTickx(1);
   canvas1->SetTicky(1);
   canvas1->SetRightMargin(0.05);
   canvas1->SetBottomMargin(0.12);
   canvas1->SetFrameBorderMode(0);
   canvas1->SetFrameBorderMode(0);
   
   TH2F *frameBDT__37 = new TH2F("frameBDT__37","TMVA overtraining check for classifier: BDT",500,-0.5558584,0.3056503,500,0,9.034575);
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
   
   TH1D *MVA_BDT_S__38 = new TH1D("MVA_BDT_S__38","TMVA overtraining check for classifier: BDT",40,-0.5558584,0.3056503);
   MVA_BDT_S__38->SetBinContent(18,0.00443055);
   MVA_BDT_S__38->SetBinContent(19,0.02549562);
   MVA_BDT_S__38->SetBinContent(20,0.03311312);
   MVA_BDT_S__38->SetBinContent(21,0.04710453);
   MVA_BDT_S__38->SetBinContent(22,0.1101438);
   MVA_BDT_S__38->SetBinContent(23,0.1810326);
   MVA_BDT_S__38->SetBinContent(24,0.3333864);
   MVA_BDT_S__38->SetBinContent(25,0.6093279);
   MVA_BDT_S__38->SetBinContent(26,1.018037);
   MVA_BDT_S__38->SetBinContent(27,1.684263);
   MVA_BDT_S__38->SetBinContent(28,2.755942);
   MVA_BDT_S__38->SetBinContent(29,3.752058);
   MVA_BDT_S__38->SetBinContent(30,4.936593);
   MVA_BDT_S__38->SetBinContent(31,6.155477);
   MVA_BDT_S__38->SetBinContent(32,6.773202);
   MVA_BDT_S__38->SetBinContent(33,6.466642);
   MVA_BDT_S__38->SetBinContent(34,5.158592);
   MVA_BDT_S__38->SetBinContent(35,3.387886);
   MVA_BDT_S__38->SetBinContent(36,1.800624);
   MVA_BDT_S__38->SetBinContent(37,0.8500644);
   MVA_BDT_S__38->SetBinContent(38,0.2339672);
   MVA_BDT_S__38->SetBinContent(39,0.08853491);
   MVA_BDT_S__38->SetBinContent(40,0.02425202);
   MVA_BDT_S__38->SetBinError(18,0.002621011);
   MVA_BDT_S__38->SetBinError(19,0.006525742);
   MVA_BDT_S__38->SetBinError(20,0.007401228);
   MVA_BDT_S__38->SetBinError(21,0.008810072);
   MVA_BDT_S__38->SetBinError(22,0.0134852);
   MVA_BDT_S__38->SetBinError(23,0.01708115);
   MVA_BDT_S__38->SetBinError(24,0.02360902);
   MVA_BDT_S__38->SetBinError(25,0.03172722);
   MVA_BDT_S__38->SetBinError(26,0.04126196);
   MVA_BDT_S__38->SetBinError(27,0.0528795);
   MVA_BDT_S__38->SetBinError(28,0.06806977);
   MVA_BDT_S__38->SetBinError(29,0.0793422);
   MVA_BDT_S__38->SetBinError(30,0.09092493);
   MVA_BDT_S__38->SetBinError(31,0.1012537);
   MVA_BDT_S__38->SetBinError(32,0.1062341);
   MVA_BDT_S__38->SetBinError(33,0.1040591);
   MVA_BDT_S__38->SetBinError(34,0.09295877);
   MVA_BDT_S__38->SetBinError(35,0.07524252);
   MVA_BDT_S__38->SetBinError(36,0.0546427);
   MVA_BDT_S__38->SetBinError(37,0.0379197);
   MVA_BDT_S__38->SetBinError(38,0.01951187);
   MVA_BDT_S__38->SetBinError(39,0.01211687);
   MVA_BDT_S__38->SetBinError(40,0.006406152);
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
   
   TH1D *MVA_BDT_B__39 = new TH1D("MVA_BDT_B__39","MVA_BDT_B",40,-0.5558584,0.3056503);
   MVA_BDT_B__39->SetBinContent(1,0.007722718);
   MVA_BDT_B__39->SetBinContent(2,0.02865164);
   MVA_BDT_B__39->SetBinContent(3,0.08770755);
   MVA_BDT_B__39->SetBinContent(4,0.3204677);
   MVA_BDT_B__39->SetBinContent(5,0.6196595);
   MVA_BDT_B__39->SetBinContent(6,0.9683939);
   MVA_BDT_B__39->SetBinContent(7,1.402289);
   MVA_BDT_B__39->SetBinContent(8,2.015844);
   MVA_BDT_B__39->SetBinContent(9,2.381575);
   MVA_BDT_B__39->SetBinContent(10,2.539589);
   MVA_BDT_B__39->SetBinContent(11,2.600839);
   MVA_BDT_B__39->SetBinContent(12,2.630827);
   MVA_BDT_B__39->SetBinContent(13,2.588105);
   MVA_BDT_B__39->SetBinContent(14,2.526185);
   MVA_BDT_B__39->SetBinContent(15,2.571438);
   MVA_BDT_B__39->SetBinContent(16,2.565228);
   MVA_BDT_B__39->SetBinContent(17,2.487867);
   MVA_BDT_B__39->SetBinContent(18,2.28826);
   MVA_BDT_B__39->SetBinContent(19,2.232181);
   MVA_BDT_B__39->SetBinContent(20,2.112962);
   MVA_BDT_B__39->SetBinContent(21,1.83675);
   MVA_BDT_B__39->SetBinContent(22,1.619452);
   MVA_BDT_B__39->SetBinContent(23,1.442914);
   MVA_BDT_B__39->SetBinContent(24,1.229256);
   MVA_BDT_B__39->SetBinContent(25,1.11557);
   MVA_BDT_B__39->SetBinContent(26,0.984658);
   MVA_BDT_B__39->SetBinContent(27,0.8150151);
   MVA_BDT_B__39->SetBinContent(28,0.7310995);
   MVA_BDT_B__39->SetBinContent(29,0.525314);
   MVA_BDT_B__39->SetBinContent(30,0.3993418);
   MVA_BDT_B__39->SetBinContent(31,0.2999742);
   MVA_BDT_B__39->SetBinContent(32,0.2123265);
   MVA_BDT_B__39->SetBinContent(33,0.1463108);
   MVA_BDT_B__39->SetBinContent(34,0.04726408);
   MVA_BDT_B__39->SetBinContent(35,0.02929792);
   MVA_BDT_B__39->SetBinContent(36,0.01479782);
   MVA_BDT_B__39->SetBinContent(38,0.005038423);
   MVA_BDT_B__39->SetBinError(1,0.003286457);
   MVA_BDT_B__39->SetBinError(2,0.00421488);
   MVA_BDT_B__39->SetBinError(3,0.008898472);
   MVA_BDT_B__39->SetBinError(4,0.01834812);
   MVA_BDT_B__39->SetBinError(5,0.02554215);
   MVA_BDT_B__39->SetBinError(6,0.03450825);
   MVA_BDT_B__39->SetBinError(7,0.04280115);
   MVA_BDT_B__39->SetBinError(8,0.05370958);
   MVA_BDT_B__39->SetBinError(9,0.06026952);
   MVA_BDT_B__39->SetBinError(10,0.06358395);
   MVA_BDT_B__39->SetBinError(11,0.0648555);
   MVA_BDT_B__39->SetBinError(12,0.06679496);
   MVA_BDT_B__39->SetBinError(13,0.06750667);
   MVA_BDT_B__39->SetBinError(14,0.06682839);
   MVA_BDT_B__39->SetBinError(15,0.06890358);
   MVA_BDT_B__39->SetBinError(16,0.06937604);
   MVA_BDT_B__39->SetBinError(17,0.06829699);
   MVA_BDT_B__39->SetBinError(18,0.06617152);
   MVA_BDT_B__39->SetBinError(19,0.06553281);
   MVA_BDT_B__39->SetBinError(20,0.06460508);
   MVA_BDT_B__39->SetBinError(21,0.05990665);
   MVA_BDT_B__39->SetBinError(22,0.05647135);
   MVA_BDT_B__39->SetBinError(23,0.05365131);
   MVA_BDT_B__39->SetBinError(24,0.04991002);
   MVA_BDT_B__39->SetBinError(25,0.04739266);
   MVA_BDT_B__39->SetBinError(26,0.04460624);
   MVA_BDT_B__39->SetBinError(27,0.04044815);
   MVA_BDT_B__39->SetBinError(28,0.03893336);
   MVA_BDT_B__39->SetBinError(29,0.03310866);
   MVA_BDT_B__39->SetBinError(30,0.02939027);
   MVA_BDT_B__39->SetBinError(31,0.02511844);
   MVA_BDT_B__39->SetBinError(32,0.02168809);
   MVA_BDT_B__39->SetBinError(33,0.0181112);
   MVA_BDT_B__39->SetBinError(34,0.01026343);
   MVA_BDT_B__39->SetBinError(35,0.007997964);
   MVA_BDT_B__39->SetBinError(36,0.00572843);
   MVA_BDT_B__39->SetBinError(38,0.003562703);
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
   
   TH1D *MVA_BDT_Train_S__40 = new TH1D("MVA_BDT_Train_S__40","MVA_BDT_Train_S",40,-0.5558584,0.3056503);
   MVA_BDT_Train_S__40->SetBinContent(16,0.001945623);
   MVA_BDT_Train_S__40->SetBinContent(17,0.001945623);
   MVA_BDT_Train_S__40->SetBinContent(18,0.007626391);
   MVA_BDT_Train_S__40->SetBinContent(19,0.01151764);
   MVA_BDT_Train_S__40->SetBinContent(20,0.01914403);
   MVA_BDT_Train_S__40->SetBinContent(21,0.03260729);
   MVA_BDT_Train_S__40->SetBinContent(22,0.06085502);
   MVA_BDT_Train_S__40->SetBinContent(23,0.1620998);
   MVA_BDT_Train_S__40->SetBinContent(24,0.2465346);
   MVA_BDT_Train_S__40->SetBinContent(25,0.5126853);
   MVA_BDT_Train_S__40->SetBinContent(26,0.976494);
   MVA_BDT_Train_S__40->SetBinContent(27,1.714091);
   MVA_BDT_Train_S__40->SetBinContent(28,2.782577);
   MVA_BDT_Train_S__40->SetBinContent(29,3.9187);
   MVA_BDT_Train_S__40->SetBinContent(30,5.154499);
   MVA_BDT_Train_S__40->SetBinContent(31,6.307097);
   MVA_BDT_Train_S__40->SetBinContent(32,6.949674);
   MVA_BDT_Train_S__40->SetBinContent(33,6.42803);
   MVA_BDT_Train_S__40->SetBinContent(34,5.001876);
   MVA_BDT_Train_S__40->SetBinContent(35,3.248492);
   MVA_BDT_Train_S__40->SetBinContent(36,1.694246);
   MVA_BDT_Train_S__40->SetBinContent(37,0.7828034);
   MVA_BDT_Train_S__40->SetBinContent(38,0.2828808);
   MVA_BDT_Train_S__40->SetBinContent(39,0.1109729);
   MVA_BDT_Train_S__40->SetBinContent(40,0.02077745);
   MVA_BDT_Train_S__40->SetBinContent(41,0.005136294);
   MVA_BDT_Train_S__40->SetBinError(16,0.001945623);
   MVA_BDT_Train_S__40->SetBinError(17,0.001945623);
   MVA_BDT_Train_S__40->SetBinError(18,0.003495902);
   MVA_BDT_Train_S__40->SetBinError(19,0.004448846);
   MVA_BDT_Train_S__40->SetBinError(20,0.005658053);
   MVA_BDT_Train_S__40->SetBinError(21,0.007455953);
   MVA_BDT_Train_S__40->SetBinError(22,0.009739107);
   MVA_BDT_Train_S__40->SetBinError(23,0.01613795);
   MVA_BDT_Train_S__40->SetBinError(24,0.01985404);
   MVA_BDT_Train_S__40->SetBinError(25,0.02917112);
   MVA_BDT_Train_S__40->SetBinError(26,0.03989969);
   MVA_BDT_Train_S__40->SetBinError(27,0.05354857);
   MVA_BDT_Train_S__40->SetBinError(28,0.06825305);
   MVA_BDT_Train_S__40->SetBinError(29,0.08134411);
   MVA_BDT_Train_S__40->SetBinError(30,0.09310211);
   MVA_BDT_Train_S__40->SetBinError(31,0.1025744);
   MVA_BDT_Train_S__40->SetBinError(32,0.1078848);
   MVA_BDT_Train_S__40->SetBinError(33,0.1035298);
   MVA_BDT_Train_S__40->SetBinError(34,0.09123571);
   MVA_BDT_Train_S__40->SetBinError(35,0.07367516);
   MVA_BDT_Train_S__40->SetBinError(36,0.05322781);
   MVA_BDT_Train_S__40->SetBinError(37,0.03613655);
   MVA_BDT_Train_S__40->SetBinError(38,0.02182288);
   MVA_BDT_Train_S__40->SetBinError(39,0.01358346);
   MVA_BDT_Train_S__40->SetBinError(40,0.005473398);
   MVA_BDT_Train_S__40->SetBinError(41,0.003020106);
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
   
   TH1D *MVA_BDT_Train_B__41 = new TH1D("MVA_BDT_Train_B__41","MVA_BDT_Train_B",40,-0.5558584,0.3056503);
   MVA_BDT_Train_B__41->SetBinContent(0,0.0006387634);
   MVA_BDT_Train_B__41->SetBinContent(1,0.00472446);
   MVA_BDT_Train_B__41->SetBinContent(2,0.02661879);
   MVA_BDT_Train_B__41->SetBinContent(3,0.1054721);
   MVA_BDT_Train_B__41->SetBinContent(4,0.2767597);
   MVA_BDT_Train_B__41->SetBinContent(5,0.6456193);
   MVA_BDT_Train_B__41->SetBinContent(6,0.9286394);
   MVA_BDT_Train_B__41->SetBinContent(7,1.41624);
   MVA_BDT_Train_B__41->SetBinContent(8,1.99057);
   MVA_BDT_Train_B__41->SetBinContent(9,2.404519);
   MVA_BDT_Train_B__41->SetBinContent(10,2.57326);
   MVA_BDT_Train_B__41->SetBinContent(11,2.625754);
   MVA_BDT_Train_B__41->SetBinContent(12,2.54789);
   MVA_BDT_Train_B__41->SetBinContent(13,2.591736);
   MVA_BDT_Train_B__41->SetBinContent(14,2.673663);
   MVA_BDT_Train_B__41->SetBinContent(15,2.609357);
   MVA_BDT_Train_B__41->SetBinContent(16,2.68213);
   MVA_BDT_Train_B__41->SetBinContent(17,2.416289);
   MVA_BDT_Train_B__41->SetBinContent(18,2.444834);
   MVA_BDT_Train_B__41->SetBinContent(19,2.201238);
   MVA_BDT_Train_B__41->SetBinContent(20,1.981564);
   MVA_BDT_Train_B__41->SetBinContent(21,1.796836);
   MVA_BDT_Train_B__41->SetBinContent(22,1.565333);
   MVA_BDT_Train_B__41->SetBinContent(23,1.466753);
   MVA_BDT_Train_B__41->SetBinContent(24,1.243557);
   MVA_BDT_Train_B__41->SetBinContent(25,1.106641);
   MVA_BDT_Train_B__41->SetBinContent(26,0.9818876);
   MVA_BDT_Train_B__41->SetBinContent(27,0.9735031);
   MVA_BDT_Train_B__41->SetBinContent(28,0.7937259);
   MVA_BDT_Train_B__41->SetBinContent(29,0.5802447);
   MVA_BDT_Train_B__41->SetBinContent(30,0.3550696);
   MVA_BDT_Train_B__41->SetBinContent(31,0.1976356);
   MVA_BDT_Train_B__41->SetBinContent(32,0.1127218);
   MVA_BDT_Train_B__41->SetBinContent(33,0.06867944);
   MVA_BDT_Train_B__41->SetBinContent(34,0.02401505);
   MVA_BDT_Train_B__41->SetBinContent(35,0.01669281);
   MVA_BDT_Train_B__41->SetBinError(0,0.0006387634);
   MVA_BDT_Train_B__41->SetBinError(1,0.00145663);
   MVA_BDT_Train_B__41->SetBinError(2,0.005073163);
   MVA_BDT_Train_B__41->SetBinError(3,0.009739759);
   MVA_BDT_Train_B__41->SetBinError(4,0.01623033);
   MVA_BDT_Train_B__41->SetBinError(5,0.02504005);
   MVA_BDT_Train_B__41->SetBinError(6,0.03170754);
   MVA_BDT_Train_B__41->SetBinError(7,0.04230194);
   MVA_BDT_Train_B__41->SetBinError(8,0.05293245);
   MVA_BDT_Train_B__41->SetBinError(9,0.06016774);
   MVA_BDT_Train_B__41->SetBinError(10,0.06391081);
   MVA_BDT_Train_B__41->SetBinError(11,0.06579148);
   MVA_BDT_Train_B__41->SetBinError(12,0.06536877);
   MVA_BDT_Train_B__41->SetBinError(13,0.06703646);
   MVA_BDT_Train_B__41->SetBinError(14,0.07049321);
   MVA_BDT_Train_B__41->SetBinError(15,0.06934257);
   MVA_BDT_Train_B__41->SetBinError(16,0.07134423);
   MVA_BDT_Train_B__41->SetBinError(17,0.06770886);
   MVA_BDT_Train_B__41->SetBinError(18,0.06920773);
   MVA_BDT_Train_B__41->SetBinError(19,0.06526719);
   MVA_BDT_Train_B__41->SetBinError(20,0.06189201);
   MVA_BDT_Train_B__41->SetBinError(21,0.05954647);
   MVA_BDT_Train_B__41->SetBinError(22,0.05587512);
   MVA_BDT_Train_B__41->SetBinError(23,0.05403813);
   MVA_BDT_Train_B__41->SetBinError(24,0.049958);
   MVA_BDT_Train_B__41->SetBinError(25,0.04720687);
   MVA_BDT_Train_B__41->SetBinError(26,0.04428891);
   MVA_BDT_Train_B__41->SetBinError(27,0.04505302);
   MVA_BDT_Train_B__41->SetBinError(28,0.04066815);
   MVA_BDT_Train_B__41->SetBinError(29,0.0349095);
   MVA_BDT_Train_B__41->SetBinError(30,0.02715529);
   MVA_BDT_Train_B__41->SetBinError(31,0.02013147);
   MVA_BDT_Train_B__41->SetBinError(32,0.0154887);
   MVA_BDT_Train_B__41->SetBinError(33,0.01231224);
   MVA_BDT_Train_B__41->SetBinError(34,0.007144376);
   MVA_BDT_Train_B__41->SetBinError(35,0.006369758);
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
   TText *text = new TText(0.12,0.74,"Kolmogorov-Smirnov test: signal (background) probability = 0.071 (0.233)");
   text->SetNDC();
   text->SetTextSize(0.032);
   text->Draw();
   
   TH2F *frameBDT__42 = new TH2F("frameBDT__42","TMVA overtraining check for classifier: BDT",500,-0.5558584,0.3056503,500,0,9.034575);
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
