#include <RtypesCore.h>
#include <TCanvas.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"

// qqqq = 0; qqtt = 1; qqll = 2; qqvv = 3; qqlv = 4; qqtv = 5; ttll = 6; tttt = 7; qq = 8
TH1F* bg_template = new TH1F("bg_temp", "Background template", 100, -.6, .4);
TH1F* sig_template = new TH1F("sig_temp", "Signal template", 100, -.6, .4);

double bg_temp(Double_t x){
    return bg_template->GetBinContent(bg_template->FindBin(x));
}

double sig_temp(Double_t x){
    return sig_template->GetBinContent(sig_template->FindBin(x));
}

double fit_func(double *x, double *par){
    //std::cout << "fit_func called\n";
    return par[0] * bg_temp(x[0]) + par[1] * sig_temp(x[0]);
}

void find_limit(){
    std::string procId[9] = {"qqqq", "qqtt", "qqll", "qqvv", "qqlv", "qqtv", "ttll", "tttt", "qq"};
    TFile* bdt_output_tree = new TFile("train_bdt_qqll.root");
    TTree* bdt_tree = (TTree*)bdt_output_tree->Get("dataset/TrainTree");

    Float_t BDTresp;
    Float_t Iproc;
    bdt_tree->SetBranchAddress("BDT", &BDTresp);
    bdt_tree->SetBranchAddress("Iproc", &Iproc);
    Long64_t genEntries = bdt_tree->GetEntries();
    TH1F* bdt_hists[10];
    bdt_hists[9] = new TH1F("bdt_dist","BDT response distribution",100,-.6,.4);
    for (int i = 0; i <  sizeof(procId) / sizeof(std::string); i++){
        TString name = "bdt_dist " + procId[i];
        TString title = "BDT response distribution for " + procId[i];
        bdt_hists[i] = new TH1F(name,title,100,-.6,.4);
    }
    for(Int_t entry = 0; entry < genEntries; ++entry){
        bdt_tree->GetEntry(entry);
        Int_t indexproc = std::round(Iproc);
        if (indexproc >= 0){
            bdt_hists[indexproc]->Fill(BDTresp);
        }
        else{
            bdt_hists[9]->Fill(BDTresp);
        }
        //bdt_dist->Fill(BDTresp);
    }
    //std::cout<< sizeof(bdt_hists) / sizeof(TH1F*) << "\n";
    /*
    TCanvas* c = new TCanvas;
    c->SetLogy();
    for(int i = 0; i < sizeof(bdt_hists) / sizeof(TH1F*); i++){
        bdt_hists[i]->GetYaxis()->SetRangeUser(1, 5000);
        bdt_hists[i]->Draw("SAME");
    }
    */
    TRandom3* rnd = new TRandom3(0);
    Int_t randbin;
    Double_t alpha = 0.001;
    TH1F* gen_bg = new TH1F("gen_bg", "Generated background", 100, -.6, .4);
    TH1F* gen_sig = new TH1F("gen_sig", "Generated signal", 100, -.6, .4);

    for(int i = 0; i < 100; i++){
        randbin = rnd->Poisson(alpha * bdt_hists[9]->GetBinContent(i+1));
            for(int k = 0; k < randbin; k++){
                gen_sig->Fill(0.01*i - 0.599);
            }
    }

    for(int i = 0; i < sizeof(bdt_hists) / sizeof(TH1F*) - 1; i++){
        //if(bdt_hists[i]->Integral()>0) gen_bg->FillRandom(bdt_hists[i], bdt_hists[i]->Integral());
        bg_template->Add(bdt_hists[i]);
        for(int j = 0; j < 100; j++){
            randbin = rnd->Poisson(bdt_hists[i]->GetBinContent(j+1));
            for(int k = 0; k < randbin; k++){
                gen_bg->Fill(0.01*j - 0.599);
            }
        }
    }

    *sig_template = *bdt_hists[9];

    *bg_template = *bg_template * (1/bg_template->Integral());
    *sig_template = *sig_template * (1/sig_template->Integral());

    std::cout << gen_bg->Integral()<< "\n";
    std::cout << bdt_hists[9]->Integral() << "\n";

    std::cout << gen_sig->Integral() << "\n";
    TH1F* data = new TH1F("data", "Combined generated data", 100, -0.6, 0.4);
    data->Add(gen_bg, gen_sig);


    auto func = new TF1("fitting_func",&fit_func,-0.6,0.4,2);

    func->SetParNames("Background","Signal");

    //hist_to_fit.at(0)->GetListOfFunctions()->ls();
    auto fitresult = data->Fit(func, "S");
    std::cout << fitresult << "\n";

    /*
    TCanvas* c1 = new TCanvas;
    gen_bg->SetLineColor(kRed);
    gen_bg->Draw();
    gen_sig->Draw("SAME");
    */
    //bdt_dist->Draw();
}