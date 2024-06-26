#ifdef __CLING__
#pragma cling optimize(0)
#endif
void train_bdt_qqll_res()
{
//=========Macro generated from canvas: canvas1/TMVA comparison BDT
//=========  (Mon Jul  1 17:41:49 2024) by ROOT version 6.32.02
   TCanvas *canvas1 = new TCanvas("canvas1", "TMVA comparison BDT",200,81,600,468);
   gStyle->SetOptStat(0);
   canvas1->Range(-0.6026714,-1.442664,0.4141938,10.57954);
   canvas1->SetFillColor(0);
   canvas1->SetBorderMode(0);
   canvas1->SetBorderSize(2);
   canvas1->SetTickx(1);
   canvas1->SetTicky(1);
   canvas1->SetRightMargin(0.05);
   canvas1->SetBottomMargin(0.12);
   canvas1->SetFrameBorderMode(0);
   canvas1->SetFrameBorderMode(0);
   
   TH2F *frameBDT__37 = new TH2F("frameBDT__37","TMVA overtraining check for classifier: BDT",500,-0.4959005,0.3633506,500,0,9.377316);
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
   
   TH1D *MVA_BDT_S__38 = new TH1D("MVA_BDT_S__38","TMVA overtraining check for classifier: BDT",40,-0.4959005,0.3633506);
   MVA_BDT_S__38->SetBinContent(13,0.002397248);
   MVA_BDT_S__38->SetBinContent(14,0.003595873);
   MVA_BDT_S__38->SetBinContent(15,0.003595873);
   MVA_BDT_S__38->SetBinContent(16,0.01198624);
   MVA_BDT_S__38->SetBinContent(17,0.01318487);
   MVA_BDT_S__38->SetBinContent(18,0.02517111);
   MVA_BDT_S__38->SetBinContent(19,0.03595873);
   MVA_BDT_S__38->SetBinContent(20,0.07311608);
   MVA_BDT_S__38->SetBinContent(21,0.1330473);
   MVA_BDT_S__38->SetBinContent(22,0.1797936);
   MVA_BDT_S__38->SetBinContent(23,0.3440051);
   MVA_BDT_S__38->SetBinContent(24,0.5190043);
   MVA_BDT_S__38->SetBinContent(25,0.7958865);
   MVA_BDT_S__38->SetBinContent(26,1.28972);
   MVA_BDT_S__38->SetBinContent(27,1.998107);
   MVA_BDT_S__38->SetBinContent(28,3.103238);
   MVA_BDT_S__38->SetBinContent(29,4.638676);
   MVA_BDT_S__38->SetBinContent(30,6.164524);
   MVA_BDT_S__38->SetBinContent(31,7.026335);
   MVA_BDT_S__38->SetBinContent(32,6.930445);
   MVA_BDT_S__38->SetBinContent(33,5.842094);
   MVA_BDT_S__38->SetBinContent(34,3.968645);
   MVA_BDT_S__38->SetBinContent(35,2.16951);
   MVA_BDT_S__38->SetBinContent(36,0.9109544);
   MVA_BDT_S__38->SetBinContent(37,0.2744849);
   MVA_BDT_S__38->SetBinContent(38,0.08150645);
   MVA_BDT_S__38->SetBinContent(39,0.01198624);
   MVA_BDT_S__38->SetBinContent(40,0.001198624);
   MVA_BDT_S__38->SetBinError(13,0.001695111);
   MVA_BDT_S__38->SetBinError(14,0.002076078);
   MVA_BDT_S__38->SetBinError(15,0.002076078);
   MVA_BDT_S__38->SetBinError(16,0.003790383);
   MVA_BDT_S__38->SetBinError(17,0.003975387);
   MVA_BDT_S__38->SetBinError(18,0.005492786);
   MVA_BDT_S__38->SetBinError(19,0.006565135);
   MVA_BDT_S__38->SetBinError(20,0.009361554);
   MVA_BDT_S__38->SetBinError(21,0.01262829);
   MVA_BDT_S__38->SetBinError(22,0.01468009);
   MVA_BDT_S__38->SetBinError(23,0.02030598);
   MVA_BDT_S__38->SetBinError(24,0.02494175);
   MVA_BDT_S__38->SetBinError(25,0.03088638);
   MVA_BDT_S__38->SetBinError(26,0.0393178);
   MVA_BDT_S__38->SetBinError(27,0.04893852);
   MVA_BDT_S__38->SetBinError(28,0.06098866);
   MVA_BDT_S__38->SetBinError(29,0.0745656);
   MVA_BDT_S__38->SetBinError(30,0.08595899);
   MVA_BDT_S__38->SetBinError(31,0.0917711);
   MVA_BDT_S__38->SetBinError(32,0.09114274);
   MVA_BDT_S__38->SetBinError(33,0.0836808);
   MVA_BDT_S__38->SetBinError(34,0.06897038);
   MVA_BDT_S__38->SetBinError(35,0.05099438);
   MVA_BDT_S__38->SetBinError(36,0.03304379);
   MVA_BDT_S__38->SetBinError(37,0.01813848);
   MVA_BDT_S__38->SetBinError(38,0.009884108);
   MVA_BDT_S__38->SetBinError(39,0.003790383);
   MVA_BDT_S__38->SetBinError(40,0.001198624);
   MVA_BDT_S__38->SetEntries(38838);
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
   
   TH1D *MVA_BDT_B__39 = new TH1D("MVA_BDT_B__39","MVA_BDT_B",40,-0.4959005,0.3633506);
   MVA_BDT_B__39->SetBinContent(1,0.008556292);
   MVA_BDT_B__39->SetBinContent(2,0.0682671);
   MVA_BDT_B__39->SetBinContent(3,0.1661568);
   MVA_BDT_B__39->SetBinContent(4,0.4335211);
   MVA_BDT_B__39->SetBinContent(5,0.6869377);
   MVA_BDT_B__39->SetBinContent(6,1.048568);
   MVA_BDT_B__39->SetBinContent(7,1.678988);
   MVA_BDT_B__39->SetBinContent(8,2.267638);
   MVA_BDT_B__39->SetBinContent(9,2.550026);
   MVA_BDT_B__39->SetBinContent(10,3.126909);
   MVA_BDT_B__39->SetBinContent(11,3.365222);
   MVA_BDT_B__39->SetBinContent(12,3.568357);
   MVA_BDT_B__39->SetBinContent(13,3.86711);
   MVA_BDT_B__39->SetBinContent(14,3.801104);
   MVA_BDT_B__39->SetBinContent(15,3.444323);
   MVA_BDT_B__39->SetBinContent(16,3.298218);
   MVA_BDT_B__39->SetBinContent(17,2.974827);
   MVA_BDT_B__39->SetBinContent(18,2.498552);
   MVA_BDT_B__39->SetBinContent(19,2.1412);
   MVA_BDT_B__39->SetBinContent(20,1.552248);
   MVA_BDT_B__39->SetBinContent(21,1.224808);
   MVA_BDT_B__39->SetBinContent(22,0.8693336);
   MVA_BDT_B__39->SetBinContent(23,0.6040299);
   MVA_BDT_B__39->SetBinContent(24,0.3958737);
   MVA_BDT_B__39->SetBinContent(25,0.2958216);
   MVA_BDT_B__39->SetBinContent(26,0.2450825);
   MVA_BDT_B__39->SetBinContent(27,0.1059341);
   MVA_BDT_B__39->SetBinContent(28,0.09473036);
   MVA_BDT_B__39->SetBinContent(29,0.06266842);
   MVA_BDT_B__39->SetBinContent(30,0.05236685);
   MVA_BDT_B__39->SetBinContent(31,0.02993572);
   MVA_BDT_B__39->SetBinContent(32,0.01854752);
   MVA_BDT_B__39->SetBinContent(33,0.00499083);
   MVA_BDT_B__39->SetBinContent(34,0.00131457);
   MVA_BDT_B__39->SetBinError(1,0.004584184);
   MVA_BDT_B__39->SetBinError(2,0.01439724);
   MVA_BDT_B__39->SetBinError(3,0.02236917);
   MVA_BDT_B__39->SetBinError(4,0.03570748);
   MVA_BDT_B__39->SetBinError(5,0.04252037);
   MVA_BDT_B__39->SetBinError(6,0.05198458);
   MVA_BDT_B__39->SetBinError(7,0.06949487);
   MVA_BDT_B__39->SetBinError(8,0.08096884);
   MVA_BDT_B__39->SetBinError(9,0.0851374);
   MVA_BDT_B__39->SetBinError(10,0.09673151);
   MVA_BDT_B__39->SetBinError(11,0.1006744);
   MVA_BDT_B__39->SetBinError(12,0.1048705);
   MVA_BDT_B__39->SetBinError(13,0.1100238);
   MVA_BDT_B__39->SetBinError(14,0.1092683);
   MVA_BDT_B__39->SetBinError(15,0.1042054);
   MVA_BDT_B__39->SetBinError(16,0.1023415);
   MVA_BDT_B__39->SetBinError(17,0.09734029);
   MVA_BDT_B__39->SetBinError(18,0.08844271);
   MVA_BDT_B__39->SetBinError(19,0.08332072);
   MVA_BDT_B__39->SetBinError(20,0.07061268);
   MVA_BDT_B__39->SetBinError(21,0.06377811);
   MVA_BDT_B__39->SetBinError(22,0.054006);
   MVA_BDT_B__39->SetBinError(23,0.04457807);
   MVA_BDT_B__39->SetBinError(24,0.03641352);
   MVA_BDT_B__39->SetBinError(25,0.0313075);
   MVA_BDT_B__39->SetBinError(26,0.02879373);
   MVA_BDT_B__39->SetBinError(27,0.01792327);
   MVA_BDT_B__39->SetBinError(28,0.01775922);
   MVA_BDT_B__39->SetBinError(29,0.01526105);
   MVA_BDT_B__39->SetBinError(30,0.01403248);
   MVA_BDT_B__39->SetBinError(31,0.01039483);
   MVA_BDT_B__39->SetBinError(32,0.008417449);
   MVA_BDT_B__39->SetBinError(33,0.004226746);
   MVA_BDT_B__39->SetBinError(34,0.0009976179);
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
   
   TH1D *MVA_BDT_Train_S__40 = new TH1D("MVA_BDT_Train_S__40","MVA_BDT_Train_S",40,-0.4959005,0.3633506);
   MVA_BDT_Train_S__40->SetBinContent(15,0.001198624);
   MVA_BDT_Train_S__40->SetBinContent(16,0.002397248);
   MVA_BDT_Train_S__40->SetBinContent(17,0.008390369);
   MVA_BDT_Train_S__40->SetBinContent(18,0.007191745);
   MVA_BDT_Train_S__40->SetBinContent(19,0.02636973);
   MVA_BDT_Train_S__40->SetBinContent(20,0.03116423);
   MVA_BDT_Train_S__40->SetBinContent(21,0.06112983);
   MVA_BDT_Train_S__40->SetBinContent(22,0.1833895);
   MVA_BDT_Train_S__40->SetBinContent(23,0.2948616);
   MVA_BDT_Train_S__40->SetBinContent(24,0.5142098);
   MVA_BDT_Train_S__40->SetBinContent(25,0.8534204);
   MVA_BDT_Train_S__40->SetBinContent(26,1.32448);
   MVA_BDT_Train_S__40->SetBinContent(27,2.017285);
   MVA_BDT_Train_S__40->SetBinContent(28,3.056492);
   MVA_BDT_Train_S__40->SetBinContent(29,4.659052);
   MVA_BDT_Train_S__40->SetBinContent(30,5.94997);
   MVA_BDT_Train_S__40->SetBinContent(31,7.107841);
   MVA_BDT_Train_S__40->SetBinContent(32,7.21332);
   MVA_BDT_Train_S__40->SetBinContent(33,5.89843);
   MVA_BDT_Train_S__40->SetBinContent(34,3.926693);
   MVA_BDT_Train_S__40->SetBinContent(35,2.068825);
   MVA_BDT_Train_S__40->SetBinContent(36,0.9672897);
   MVA_BDT_Train_S__40->SetBinContent(37,0.299656);
   MVA_BDT_Train_S__40->SetBinContent(38,0.06472571);
   MVA_BDT_Train_S__40->SetBinContent(39,0.01438349);
   MVA_BDT_Train_S__40->SetBinError(15,0.001198624);
   MVA_BDT_Train_S__40->SetBinError(16,0.001695111);
   MVA_BDT_Train_S__40->SetBinError(17,0.003171262);
   MVA_BDT_Train_S__40->SetBinError(18,0.002936018);
   MVA_BDT_Train_S__40->SetBinError(19,0.005622046);
   MVA_BDT_Train_S__40->SetBinError(20,0.006111808);
   MVA_BDT_Train_S__40->SetBinError(21,0.008559889);
   MVA_BDT_Train_S__40->SetBinError(22,0.01482616);
   MVA_BDT_Train_S__40->SetBinError(23,0.01879969);
   MVA_BDT_Train_S__40->SetBinError(24,0.02482628);
   MVA_BDT_Train_S__40->SetBinError(25,0.03198328);
   MVA_BDT_Train_S__40->SetBinError(26,0.03984411);
   MVA_BDT_Train_S__40->SetBinError(27,0.04917282);
   MVA_BDT_Train_S__40->SetBinError(28,0.06052755);
   MVA_BDT_Train_S__40->SetBinError(29,0.0747292);
   MVA_BDT_Train_S__40->SetBinError(30,0.08444986);
   MVA_BDT_Train_S__40->SetBinError(31,0.09230185);
   MVA_BDT_Train_S__40->SetBinError(32,0.09298419);
   MVA_BDT_Train_S__40->SetBinError(33,0.08408329);
   MVA_BDT_Train_S__40->SetBinError(34,0.06860488);
   MVA_BDT_Train_S__40->SetBinError(35,0.04979703);
   MVA_BDT_Train_S__40->SetBinError(36,0.03405021);
   MVA_BDT_Train_S__40->SetBinError(37,0.01895191);
   MVA_BDT_Train_S__40->SetBinError(38,0.008808053);
   MVA_BDT_Train_S__40->SetBinError(39,0.004152156);
   MVA_BDT_Train_S__40->SetEntries(38838);

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
   
   TH1D *MVA_BDT_Train_B__41 = new TH1D("MVA_BDT_Train_B__41","MVA_BDT_Train_B",40,-0.4959005,0.3633506);
   MVA_BDT_Train_B__41->SetBinContent(0,0.004169556);
   MVA_BDT_Train_B__41->SetBinContent(1,0.01659921);
   MVA_BDT_Train_B__41->SetBinContent(2,0.05470587);
   MVA_BDT_Train_B__41->SetBinContent(3,0.1796091);
   MVA_BDT_Train_B__41->SetBinContent(4,0.3466279);
   MVA_BDT_Train_B__41->SetBinContent(5,0.7390933);
   MVA_BDT_Train_B__41->SetBinContent(6,1.094744);
   MVA_BDT_Train_B__41->SetBinContent(7,1.562034);
   MVA_BDT_Train_B__41->SetBinContent(8,2.322943);
   MVA_BDT_Train_B__41->SetBinContent(9,2.636255);
   MVA_BDT_Train_B__41->SetBinContent(10,3.065811);
   MVA_BDT_Train_B__41->SetBinContent(11,3.446632);
   MVA_BDT_Train_B__41->SetBinContent(12,3.579071);
   MVA_BDT_Train_B__41->SetBinContent(13,3.817181);
   MVA_BDT_Train_B__41->SetBinContent(14,3.830427);
   MVA_BDT_Train_B__41->SetBinContent(15,3.696683);
   MVA_BDT_Train_B__41->SetBinContent(16,3.442643);
   MVA_BDT_Train_B__41->SetBinContent(17,3.062849);
   MVA_BDT_Train_B__41->SetBinContent(18,2.359835);
   MVA_BDT_Train_B__41->SetBinContent(19,2.02556);
   MVA_BDT_Train_B__41->SetBinContent(20,1.466901);
   MVA_BDT_Train_B__41->SetBinContent(21,1.180642);
   MVA_BDT_Train_B__41->SetBinContent(22,0.9357426);
   MVA_BDT_Train_B__41->SetBinContent(23,0.6384249);
   MVA_BDT_Train_B__41->SetBinContent(24,0.4168142);
   MVA_BDT_Train_B__41->SetBinContent(25,0.2632438);
   MVA_BDT_Train_B__41->SetBinContent(26,0.1352898);
   MVA_BDT_Train_B__41->SetBinContent(27,0.08051409);
   MVA_BDT_Train_B__41->SetBinContent(28,0.0740723);
   MVA_BDT_Train_B__41->SetBinContent(29,0.03047174);
   MVA_BDT_Train_B__41->SetBinContent(30,0.02446349);
   MVA_BDT_Train_B__41->SetBinContent(31,0.01238922);
   MVA_BDT_Train_B__41->SetBinContent(32,0.00496826);
   MVA_BDT_Train_B__41->SetBinContent(33,0.008926472);
   MVA_BDT_Train_B__41->SetBinError(0,0.004169556);
   MVA_BDT_Train_B__41->SetBinError(1,0.007420305);
   MVA_BDT_Train_B__41->SetBinError(2,0.01344055);
   MVA_BDT_Train_B__41->SetBinError(3,0.02314292);
   MVA_BDT_Train_B__41->SetBinError(4,0.03049999);
   MVA_BDT_Train_B__41->SetBinError(5,0.0441109);
   MVA_BDT_Train_B__41->SetBinError(6,0.05317282);
   MVA_BDT_Train_B__41->SetBinError(7,0.06575311);
   MVA_BDT_Train_B__41->SetBinError(8,0.08128563);
   MVA_BDT_Train_B__41->SetBinError(9,0.08678834);
   MVA_BDT_Train_B__41->SetBinError(10,0.09561285);
   MVA_BDT_Train_B__41->SetBinError(11,0.1023526);
   MVA_BDT_Train_B__41->SetBinError(12,0.1041878);
   MVA_BDT_Train_B__41->SetBinError(13,0.1094913);
   MVA_BDT_Train_B__41->SetBinError(14,0.1102248);
   MVA_BDT_Train_B__41->SetBinError(15,0.1081824);
   MVA_BDT_Train_B__41->SetBinError(16,0.1038978);
   MVA_BDT_Train_B__41->SetBinError(17,0.09901383);
   MVA_BDT_Train_B__41->SetBinError(18,0.08587253);
   MVA_BDT_Train_B__41->SetBinError(19,0.07996872);
   MVA_BDT_Train_B__41->SetBinError(20,0.06748432);
   MVA_BDT_Train_B__41->SetBinError(21,0.06149478);
   MVA_BDT_Train_B__41->SetBinError(22,0.05621938);
   MVA_BDT_Train_B__41->SetBinError(23,0.04623955);
   MVA_BDT_Train_B__41->SetBinError(24,0.03692963);
   MVA_BDT_Train_B__41->SetBinError(25,0.02940505);
   MVA_BDT_Train_B__41->SetBinError(26,0.02087061);
   MVA_BDT_Train_B__41->SetBinError(27,0.01549174);
   MVA_BDT_Train_B__41->SetBinError(28,0.01590187);
   MVA_BDT_Train_B__41->SetBinError(29,0.01040269);
   MVA_BDT_Train_B__41->SetBinError(30,0.009410833);
   MVA_BDT_Train_B__41->SetBinError(31,0.006069005);
   MVA_BDT_Train_B__41->SetBinError(32,0.004207631);
   MVA_BDT_Train_B__41->SetBinError(33,0.00591314);
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
   TText *text = new TText(0.12,0.74,"Kolmogorov-Smirnov test: signal (background) probability = 0.093 (0.127)");
   text->SetNDC();
   text->SetTextSize(0.032);
   text->Draw();
   
   TH2F *frameBDT__42 = new TH2F("frameBDT__42","TMVA overtraining check for classifier: BDT",500,-0.4959005,0.3633506,500,0,9.377316);
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
