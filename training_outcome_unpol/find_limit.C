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

// qqqq = 0; qqtt = 1; qqll = 2; qqvv = 3; qqlv = 4; qqtv = 5; ttll = 6; tttt = 7; qq = 8

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
    bdt_hists[9] = new TH1F("bdt_dist","BDT response distribution",100,-1,1);
    for (int i = 0; i <  sizeof(procId) / sizeof(std::string); i++){
        TString name = "bdt_dist " + procId[i];
        TString title = "BDT response distribution for " + procId[i];
        bdt_hists[i] = new TH1F(name,title,100,-1,1);
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
    std::cout<< sizeof(bdt_hists) / sizeof(TH1F*) << "\n";
    TCanvas* c = new TCanvas;
    c->SetLogy();
    for(int i = 0; i < sizeof(bdt_hists) / sizeof(TH1F*); i++){
        bdt_hists[i]->GetYaxis()->SetRangeUser(1, 5000);
        bdt_hists[i]->Draw("SAME");
    }
    //bdt_dist->Draw();
}