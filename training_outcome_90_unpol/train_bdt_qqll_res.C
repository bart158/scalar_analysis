#ifdef __CLING__
#pragma cling optimize(0)
#endif
void train_bdt_qqll_res()
{
//=========Macro generated from canvas: canvas1/TMVA comparison BDT
//=========  (Sun Jun 23 22:08:58 2024) by ROOT version 6.32.00
   TCanvas *canvas1 = new TCanvas("canvas1", "TMVA comparison BDT",200,81,600,468);
   gStyle->SetOptStat(0);
   canvas1->Range(-0.6329625,-1.71071,0.2913666,12.5452);
   canvas1->SetFillColor(0);
   canvas1->SetBorderMode(0);
   canvas1->SetBorderSize(2);
   canvas1->SetTickx(1);
   canvas1->SetTicky(1);
   canvas1->SetRightMargin(0.05);
   canvas1->SetBottomMargin(0.12);
   canvas1->SetFrameBorderMode(0);
   canvas1->SetFrameBorderMode(0);
   
   TH2F *frameBDT__37 = new TH2F("frameBDT__37","TMVA overtraining check for classifier: BDT",500,-0.535908,0.2451501,500,0,11.11961);
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
   
   TH1D *MVA_BDT_S__38 = new TH1D("MVA_BDT_S__38","TMVA overtraining check for classifier: BDT",40,-0.535908,0.2451501);
   MVA_BDT_S__38->SetBinContent(13,0.003294155);
   MVA_BDT_S__38->SetBinContent(14,0.003294155);
   MVA_BDT_S__38->SetBinContent(15,0.004941233);
   MVA_BDT_S__38->SetBinContent(16,0.001647078);
   MVA_BDT_S__38->SetBinContent(17,0.009882465);
   MVA_BDT_S__38->SetBinContent(18,0.01976493);
   MVA_BDT_S__38->SetBinContent(19,0.0296474);
   MVA_BDT_S__38->SetBinContent(20,0.04611817);
   MVA_BDT_S__38->SetBinContent(21,0.05764771);
   MVA_BDT_S__38->SetBinContent(22,0.09058926);
   MVA_BDT_S__38->SetBinContent(23,0.1614136);
   MVA_BDT_S__38->SetBinContent(24,0.2684736);
   MVA_BDT_S__38->SetBinContent(25,0.3689454);
   MVA_BDT_S__38->SetBinContent(26,0.6291836);
   MVA_BDT_S__38->SetBinContent(27,1.031071);
   MVA_BDT_S__38->SetBinContent(28,1.648725);
   MVA_BDT_S__38->SetBinContent(29,2.643559);
   MVA_BDT_S__38->SetBinContent(30,3.979339);
   MVA_BDT_S__38->SetBinContent(31,5.784536);
   MVA_BDT_S__38->SetBinContent(32,7.59138);
   MVA_BDT_S__38->SetBinContent(33,8.546685);
   MVA_BDT_S__38->SetBinContent(34,7.372319);
   MVA_BDT_S__38->SetBinContent(35,5.255824);
   MVA_BDT_S__38->SetBinContent(36,3.252978);
   MVA_BDT_S__38->SetBinContent(37,1.686607);
   MVA_BDT_S__38->SetBinContent(38,0.5863596);
   MVA_BDT_S__38->SetBinContent(39,0.1218837);
   MVA_BDT_S__38->SetBinContent(40,0.01647078);
   MVA_BDT_S__38->SetBinError(13,0.002329319);
   MVA_BDT_S__38->SetBinError(14,0.002329319);
   MVA_BDT_S__38->SetBinError(15,0.002852822);
   MVA_BDT_S__38->SetBinError(16,0.001647078);
   MVA_BDT_S__38->SetBinError(17,0.004034499);
   MVA_BDT_S__38->SetBinError(18,0.005705644);
   MVA_BDT_S__38->SetBinError(19,0.006987958);
   MVA_BDT_S__38->SetBinError(20,0.008715515);
   MVA_BDT_S__38->SetBinError(21,0.009744242);
   MVA_BDT_S__38->SetBinError(22,0.01221505);
   MVA_BDT_S__38->SetBinError(23,0.01630524);
   MVA_BDT_S__38->SetBinError(24,0.02102848);
   MVA_BDT_S__38->SetBinError(25,0.0246512);
   MVA_BDT_S__38->SetBinError(26,0.03219183);
   MVA_BDT_S__38->SetBinError(27,0.04120987);
   MVA_BDT_S__38->SetBinError(28,0.0521112);
   MVA_BDT_S__38->SetBinError(29,0.06598596);
   MVA_BDT_S__38->SetBinError(30,0.08095851);
   MVA_BDT_S__38->SetBinError(31,0.09760932);
   MVA_BDT_S__38->SetBinError(32,0.1118195);
   MVA_BDT_S__38->SetBinError(33,0.1186468);
   MVA_BDT_S__38->SetBinError(34,0.1101943);
   MVA_BDT_S__38->SetBinError(35,0.09304166);
   MVA_BDT_S__38->SetBinError(36,0.07319773);
   MVA_BDT_S__38->SetBinError(37,0.05270648);
   MVA_BDT_S__38->SetBinError(38,0.031077);
   MVA_BDT_S__38->SetBinError(39,0.0141687);
   MVA_BDT_S__38->SetBinError(40,0.005208516);
   MVA_BDT_S__38->SetEntries(31093);
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
   
   TH1D *MVA_BDT_B__39 = new TH1D("MVA_BDT_B__39","MVA_BDT_B",40,-0.535908,0.2451501);
   MVA_BDT_B__39->SetBinContent(1,0.01571716);
   MVA_BDT_B__39->SetBinContent(2,0.1032144);
   MVA_BDT_B__39->SetBinContent(3,0.3820152);
   MVA_BDT_B__39->SetBinContent(4,0.6212423);
   MVA_BDT_B__39->SetBinContent(5,0.8172814);
   MVA_BDT_B__39->SetBinContent(6,1.027749);
   MVA_BDT_B__39->SetBinContent(7,1.298402);
   MVA_BDT_B__39->SetBinContent(8,1.455238);
   MVA_BDT_B__39->SetBinContent(9,1.614573);
   MVA_BDT_B__39->SetBinContent(10,1.664322);
   MVA_BDT_B__39->SetBinContent(11,1.761765);
   MVA_BDT_B__39->SetBinContent(12,1.841674);
   MVA_BDT_B__39->SetBinContent(13,2.075016);
   MVA_BDT_B__39->SetBinContent(14,2.019233);
   MVA_BDT_B__39->SetBinContent(15,2.22963);
   MVA_BDT_B__39->SetBinContent(16,2.464749);
   MVA_BDT_B__39->SetBinContent(17,2.482625);
   MVA_BDT_B__39->SetBinContent(18,2.65714);
   MVA_BDT_B__39->SetBinContent(19,2.631005);
   MVA_BDT_B__39->SetBinContent(20,2.563315);
   MVA_BDT_B__39->SetBinContent(21,2.298004);
   MVA_BDT_B__39->SetBinContent(22,2.069907);
   MVA_BDT_B__39->SetBinContent(23,1.88774);
   MVA_BDT_B__39->SetBinContent(24,1.665628);
   MVA_BDT_B__39->SetBinContent(25,1.523602);
   MVA_BDT_B__39->SetBinContent(26,1.374948);
   MVA_BDT_B__39->SetBinContent(27,1.274982);
   MVA_BDT_B__39->SetBinContent(28,1.281489);
   MVA_BDT_B__39->SetBinContent(29,1.204842);
   MVA_BDT_B__39->SetBinContent(30,1.224411);
   MVA_BDT_B__39->SetBinContent(31,1.034808);
   MVA_BDT_B__39->SetBinContent(32,0.9793983);
   MVA_BDT_B__39->SetBinContent(33,0.7501529);
   MVA_BDT_B__39->SetBinContent(34,0.4522488);
   MVA_BDT_B__39->SetBinContent(35,0.2882522);
   MVA_BDT_B__39->SetBinContent(36,0.1035948);
   MVA_BDT_B__39->SetBinContent(37,0.05736047);
   MVA_BDT_B__39->SetBinContent(38,0.009741759);
   MVA_BDT_B__39->SetBinContent(39,0.005563045);
   MVA_BDT_B__39->SetBinError(1,0.003744254);
   MVA_BDT_B__39->SetBinError(2,0.009492228);
   MVA_BDT_B__39->SetBinError(3,0.01710976);
   MVA_BDT_B__39->SetBinError(4,0.0237428);
   MVA_BDT_B__39->SetBinError(5,0.02757223);
   MVA_BDT_B__39->SetBinError(6,0.03468588);
   MVA_BDT_B__39->SetBinError(7,0.04141663);
   MVA_BDT_B__39->SetBinError(8,0.0441);
   MVA_BDT_B__39->SetBinError(9,0.0489976);
   MVA_BDT_B__39->SetBinError(10,0.05106873);
   MVA_BDT_B__39->SetBinError(11,0.05345592);
   MVA_BDT_B__39->SetBinError(12,0.0553204);
   MVA_BDT_B__39->SetBinError(13,0.0616409);
   MVA_BDT_B__39->SetBinError(14,0.06116723);
   MVA_BDT_B__39->SetBinError(15,0.06590878);
   MVA_BDT_B__39->SetBinError(16,0.07070753);
   MVA_BDT_B__39->SetBinError(17,0.07118074);
   MVA_BDT_B__39->SetBinError(18,0.07469887);
   MVA_BDT_B__39->SetBinError(19,0.07509198);
   MVA_BDT_B__39->SetBinError(20,0.07446348);
   MVA_BDT_B__39->SetBinError(21,0.07087416);
   MVA_BDT_B__39->SetBinError(22,0.06747732);
   MVA_BDT_B__39->SetBinError(23,0.06525408);
   MVA_BDT_B__39->SetBinError(24,0.06027818);
   MVA_BDT_B__39->SetBinError(25,0.05803388);
   MVA_BDT_B__39->SetBinError(26,0.05560311);
   MVA_BDT_B__39->SetBinError(27,0.05273364);
   MVA_BDT_B__39->SetBinError(28,0.0542753);
   MVA_BDT_B__39->SetBinError(29,0.0529992);
   MVA_BDT_B__39->SetBinError(30,0.05333319);
   MVA_BDT_B__39->SetBinError(31,0.04904601);
   MVA_BDT_B__39->SetBinError(32,0.04867324);
   MVA_BDT_B__39->SetBinError(33,0.04267326);
   MVA_BDT_B__39->SetBinError(34,0.03350156);
   MVA_BDT_B__39->SetBinError(35,0.02642796);
   MVA_BDT_B__39->SetBinError(36,0.01598436);
   MVA_BDT_B__39->SetBinError(37,0.01196761);
   MVA_BDT_B__39->SetBinError(38,0.004691143);
   MVA_BDT_B__39->SetBinError(39,0.003933667);
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
   
   TH1D *MVA_BDT_Train_S__40 = new TH1D("MVA_BDT_Train_S__40","MVA_BDT_Train_S",40,-0.535908,0.2451501);
   MVA_BDT_Train_S__40->SetBinContent(15,0.00164713);
   MVA_BDT_Train_S__40->SetBinContent(16,0.00164713);
   MVA_BDT_Train_S__40->SetBinContent(17,0.006588522);
   MVA_BDT_Train_S__40->SetBinContent(18,0.009882783);
   MVA_BDT_Train_S__40->SetBinContent(19,0.0164713);
   MVA_BDT_Train_S__40->SetBinContent(20,0.02470696);
   MVA_BDT_Train_S__40->SetBinContent(21,0.05270818);
   MVA_BDT_Train_S__40->SetBinContent(22,0.09059218);
   MVA_BDT_Train_S__40->SetBinContent(23,0.1581245);
   MVA_BDT_Train_S__40->SetBinContent(24,0.2256569);
   MVA_BDT_Train_S__40->SetBinContent(25,0.3129548);
   MVA_BDT_Train_S__40->SetBinContent(26,0.6110854);
   MVA_BDT_Train_S__40->SetBinContent(27,0.9487472);
   MVA_BDT_Train_S__40->SetBinContent(28,1.528537);
   MVA_BDT_Train_S__40->SetBinContent(29,2.585995);
   MVA_BDT_Train_S__40->SetBinContent(30,3.864168);
   MVA_BDT_Train_S__40->SetBinContent(31,5.970848);
   MVA_BDT_Train_S__40->SetBinContent(32,7.734925);
   MVA_BDT_Train_S__40->SetBinContent(33,8.553549);
   MVA_BDT_Train_S__40->SetBinContent(34,7.535622);
   MVA_BDT_Train_S__40->SetBinContent(35,5.344938);
   MVA_BDT_Train_S__40->SetBinContent(36,3.175668);
   MVA_BDT_Train_S__40->SetBinContent(37,1.717957);
   MVA_BDT_Train_S__40->SetBinContent(38,0.5814371);
   MVA_BDT_Train_S__40->SetBinContent(39,0.1433004);
   MVA_BDT_Train_S__40->SetBinContent(40,0.01482417);
   MVA_BDT_Train_S__40->SetBinContent(41,0.00164713);
   MVA_BDT_Train_S__40->SetBinError(15,0.00164713);
   MVA_BDT_Train_S__40->SetBinError(16,0.00164713);
   MVA_BDT_Train_S__40->SetBinError(17,0.003294261);
   MVA_BDT_Train_S__40->SetBinError(18,0.004034629);
   MVA_BDT_Train_S__40->SetBinError(19,0.005208684);
   MVA_BDT_Train_S__40->SetBinError(20,0.006379309);
   MVA_BDT_Train_S__40->SetBinError(21,0.009317577);
   MVA_BDT_Train_S__40->SetBinError(22,0.01221545);
   MVA_BDT_Train_S__40->SetBinError(23,0.01613852);
   MVA_BDT_Train_S__40->SetBinError(24,0.01927917);
   MVA_BDT_Train_S__40->SetBinError(25,0.02270413);
   MVA_BDT_Train_S__40->SetBinError(26,0.03172597);
   MVA_BDT_Train_S__40->SetBinError(27,0.03953113);
   MVA_BDT_Train_S__40->SetBinError(28,0.05017669);
   MVA_BDT_Train_S__40->SetBinError(29,0.06526462);
   MVA_BDT_Train_S__40->SetBinError(30,0.07977963);
   MVA_BDT_Train_S__40->SetBinError(31,0.09917039);
   MVA_BDT_Train_S__40->SetBinError(32,0.1128735);
   MVA_BDT_Train_S__40->SetBinError(33,0.1186963);
   MVA_BDT_Train_S__40->SetBinError(34,0.1114098);
   MVA_BDT_Train_S__40->SetBinError(35,0.09382863);
   MVA_BDT_Train_S__40->SetBinError(36,0.07232385);
   MVA_BDT_Train_S__40->SetBinError(37,0.05319492);
   MVA_BDT_Train_S__40->SetBinError(38,0.03094677);
   MVA_BDT_Train_S__40->SetBinError(39,0.01536341);
   MVA_BDT_Train_S__40->SetBinError(40,0.004941391);
   MVA_BDT_Train_S__40->SetBinError(41,0.00164713);
   MVA_BDT_Train_S__40->SetEntries(31093);

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
   
   TH1D *MVA_BDT_Train_B__41 = new TH1D("MVA_BDT_Train_B__41","MVA_BDT_Train_B",40,-0.535908,0.2451501);
   MVA_BDT_Train_B__41->SetBinContent(1,0.01295625);
   MVA_BDT_Train_B__41->SetBinContent(2,0.1087193);
   MVA_BDT_Train_B__41->SetBinContent(3,0.3737194);
   MVA_BDT_Train_B__41->SetBinContent(4,0.6072083);
   MVA_BDT_Train_B__41->SetBinContent(5,0.8595265);
   MVA_BDT_Train_B__41->SetBinContent(6,1.028214);
   MVA_BDT_Train_B__41->SetBinContent(7,1.262467);
   MVA_BDT_Train_B__41->SetBinContent(8,1.434142);
   MVA_BDT_Train_B__41->SetBinContent(9,1.56821);
   MVA_BDT_Train_B__41->SetBinContent(10,1.657428);
   MVA_BDT_Train_B__41->SetBinContent(11,1.823407);
   MVA_BDT_Train_B__41->SetBinContent(12,1.920058);
   MVA_BDT_Train_B__41->SetBinContent(13,2.01096);
   MVA_BDT_Train_B__41->SetBinContent(14,2.077003);
   MVA_BDT_Train_B__41->SetBinContent(15,2.279518);
   MVA_BDT_Train_B__41->SetBinContent(16,2.414414);
   MVA_BDT_Train_B__41->SetBinContent(17,2.500244);
   MVA_BDT_Train_B__41->SetBinContent(18,2.648378);
   MVA_BDT_Train_B__41->SetBinContent(19,2.697442);
   MVA_BDT_Train_B__41->SetBinContent(20,2.315224);
   MVA_BDT_Train_B__41->SetBinContent(21,2.247099);
   MVA_BDT_Train_B__41->SetBinContent(22,2.091561);
   MVA_BDT_Train_B__41->SetBinContent(23,1.812647);
   MVA_BDT_Train_B__41->SetBinContent(24,1.725387);
   MVA_BDT_Train_B__41->SetBinContent(25,1.562876);
   MVA_BDT_Train_B__41->SetBinContent(26,1.403244);
   MVA_BDT_Train_B__41->SetBinContent(27,1.393239);
   MVA_BDT_Train_B__41->SetBinContent(28,1.315694);
   MVA_BDT_Train_B__41->SetBinContent(29,1.398021);
   MVA_BDT_Train_B__41->SetBinContent(30,1.379099);
   MVA_BDT_Train_B__41->SetBinContent(31,1.269832);
   MVA_BDT_Train_B__41->SetBinContent(32,0.9272863);
   MVA_BDT_Train_B__41->SetBinContent(33,0.6018265);
   MVA_BDT_Train_B__41->SetBinContent(34,0.2703);
   MVA_BDT_Train_B__41->SetBinContent(35,0.1482592);
   MVA_BDT_Train_B__41->SetBinContent(36,0.04615058);
   MVA_BDT_Train_B__41->SetBinContent(37,0.01526691);
   MVA_BDT_Train_B__41->SetBinContent(38,0.005552963);
   MVA_BDT_Train_B__41->SetBinError(1,0.003628029);
   MVA_BDT_Train_B__41->SetBinError(2,0.01013207);
   MVA_BDT_Train_B__41->SetBinError(3,0.0171098);
   MVA_BDT_Train_B__41->SetBinError(4,0.02334579);
   MVA_BDT_Train_B__41->SetBinError(5,0.02914429);
   MVA_BDT_Train_B__41->SetBinError(6,0.03516786);
   MVA_BDT_Train_B__41->SetBinError(7,0.03920294);
   MVA_BDT_Train_B__41->SetBinError(8,0.04400711);
   MVA_BDT_Train_B__41->SetBinError(9,0.04763456);
   MVA_BDT_Train_B__41->SetBinError(10,0.05100977);
   MVA_BDT_Train_B__41->SetBinError(11,0.05506671);
   MVA_BDT_Train_B__41->SetBinError(12,0.05763704);
   MVA_BDT_Train_B__41->SetBinError(13,0.05917649);
   MVA_BDT_Train_B__41->SetBinError(14,0.06212364);
   MVA_BDT_Train_B__41->SetBinError(15,0.06708899);
   MVA_BDT_Train_B__41->SetBinError(16,0.06945442);
   MVA_BDT_Train_B__41->SetBinError(17,0.07144599);
   MVA_BDT_Train_B__41->SetBinError(18,0.07468784);
   MVA_BDT_Train_B__41->SetBinError(19,0.07633376);
   MVA_BDT_Train_B__41->SetBinError(20,0.07067136);
   MVA_BDT_Train_B__41->SetBinError(21,0.06990414);
   MVA_BDT_Train_B__41->SetBinError(22,0.06815494);
   MVA_BDT_Train_B__41->SetBinError(23,0.06311285);
   MVA_BDT_Train_B__41->SetBinError(24,0.06237643);
   MVA_BDT_Train_B__41->SetBinError(25,0.05897928);
   MVA_BDT_Train_B__41->SetBinError(26,0.0560761);
   MVA_BDT_Train_B__41->SetBinError(27,0.05625028);
   MVA_BDT_Train_B__41->SetBinError(28,0.0548267);
   MVA_BDT_Train_B__41->SetBinError(29,0.05682617);
   MVA_BDT_Train_B__41->SetBinError(30,0.05735803);
   MVA_BDT_Train_B__41->SetBinError(31,0.05512296);
   MVA_BDT_Train_B__41->SetBinError(32,0.04707624);
   MVA_BDT_Train_B__41->SetBinError(33,0.03785197);
   MVA_BDT_Train_B__41->SetBinError(34,0.02553727);
   MVA_BDT_Train_B__41->SetBinError(35,0.01890811);
   MVA_BDT_Train_B__41->SetBinError(36,0.01057684);
   MVA_BDT_Train_B__41->SetBinError(37,0.005929472);
   MVA_BDT_Train_B__41->SetBinError(38,0.003926537);
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
   TText *text = new TText(0.12,0.74,"Kolmogorov-Smirnov test: signal (background) probability = 0.017 (0.015)");
   text->SetNDC();
   text->SetTextSize(0.032);
   text->Draw();
   
   TH2F *frameBDT__42 = new TH2F("frameBDT__42","TMVA overtraining check for classifier: BDT",500,-0.535908,0.2451501,500,0,11.11961);
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
