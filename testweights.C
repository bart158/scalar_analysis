#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include <RtypesCore.h>

struct event_head {
    Int_t Ievt;
    Int_t Iproc;
    Int_t Ipol;
    Float_t Ms;
    Float_t Cs;
    Float_t Lgen; // Lgen
};

void testweights(){
    
    Double_t weights[5][4] = {
    {0.585, 0.035, 0.315, 0.065},
    {0.035, 0.585, 0.065, 0.315},
    {0.315, 0.065, 0.585, 0.035},
    {0.065, 0.315, 0.035, 0.585},
    {0.25, 0.25, 0.25, 0.25}};

    Double_t Lexp[5][4];

    for(int i = 0; i < 5; i++){
      for(int j = 0; j < 4; j++){
        if(i < 2){
          Lexp[i][j] = 900 * weights[i][j];
        }
        else if(i < 4){
          Lexp[i][j] = 100 * weights[i][j];
        }
        else{
          Lexp[i][j] = 2000 * weights[i][j];

        }
      }
    }

    TFile* pretrainfile = new TFile("trees_for_training/bbll_sig_40_eLpR_new.root");
    TTree* pretraintree = (TTree*)pretrainfile->FindObjectAny("tree");

    Long64_t entries_pretrain = pretraintree->GetEntries();
    Float_t preLgen;
    event_head head;
    pretraintree->SetBranchAddress("header", &head);
    pretraintree->GetEntry(0);
    
    preLgen = head.Lgen;

    TFile* posttrainfile = new TFile("training_outcome_40_eLpR/train_bdt_qqll.root");
    TTree* posttraintree = (TTree*)posttrainfile->FindObjectAny("TrainTree");
    TTree* posttesttree = (TTree*)posttrainfile->FindObjectAny("TestTree");

    Float_t bdt;
    Float_t weight;
    Float_t Iproc;
    Float_t Ipol;
    Float_t postLgen;
    Float_t postLgen_s;
    posttraintree->SetBranchAddress("BDT", &bdt);
    posttraintree->SetBranchAddress("weight", &weight);
    posttraintree->SetBranchAddress("Iproc", &Iproc);
    posttraintree->SetBranchAddress("Ipol", &Ipol);
    posttraintree->SetBranchAddress("Lgen", &postLgen);
    TH1F* bdt_hist = new TH1F("h_bdt", "BDT resp weighted", 100,-.6,.4);
    Long64_t entries_posttrain = posttraintree->GetEntries();
    for(Long64_t ievt = 0; ievt < entries_posttrain; ievt++){
        posttraintree->GetEntry(ievt);
        if((int)Iproc != -1 || (int)Ipol != 0) continue;
        postLgen_s = postLgen;
        bdt_hist->Fill(bdt, Lexp[0][0] / preLgen);
    }

    posttesttree->SetBranchAddress("BDT", &bdt);
    posttesttree->SetBranchAddress("weight", &weight);
    posttesttree->SetBranchAddress("Iproc", &Iproc);
    posttesttree->SetBranchAddress("Ipol", &Ipol);
    Long64_t entries_posttest = posttesttree->GetEntries();
    for(Long64_t ievt = 0; ievt < entries_posttest; ievt++){
        posttesttree->GetEntry(ievt);
        if((int)Iproc != -1 || (int)Ipol != 0) continue;
        bdt_hist->Fill(bdt, Lexp[0][0] / preLgen);
    }

    std::cout << "weighted preLgen: " << preLgen << "\n";
    std::cout << "weighted postLgen_s: " <<postLgen_s << "\n";
    std::cout << "weighted histogram integral: " << bdt_hist->Integral() << "\n";
    std::cout << "weighted histogram entries: " << bdt_hist->GetEntries() << "\n";
    std::cout << "expected number of events beamsetting=0: " << entries_pretrain * Lexp[0][0] / preLgen << "\n";
     std::cout << "expected number of events postLgen beamsetting=0: " << entries_pretrain * Lexp[0][0] / postLgen_s << "\n";
/*     std::cout << "expected number of events beamsetting=1: " << entries_pretrain * Lexp[1][0] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=2: " << entries_pretrain * Lexp[2][0] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=3: " << entries_pretrain * Lexp[3][0] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=4: " << entries_pretrain * Lexp[4][0] / preLgen << "\n";

    std::cout << "expected number of events beamsetting=0: " << entries_pretrain * Lexp[0][1] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=1: " << entries_pretrain * Lexp[1][1] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=2: " << entries_pretrain * Lexp[2][1] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=3: " << entries_pretrain * Lexp[3][1] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=4: " << entries_pretrain * Lexp[4][1] / preLgen << "\n";

    std::cout << "expected number of events beamsetting=0: " << entries_pretrain * Lexp[0][2] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=1: " << entries_pretrain * Lexp[1][2] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=2: " << entries_pretrain * Lexp[2][2] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=3: " << entries_pretrain * Lexp[3][2] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=4: " << entries_pretrain * Lexp[4][2] / preLgen << "\n";

    std::cout << "expected number of events beamsetting=0: " << entries_pretrain * Lexp[0][3] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=1: " << entries_pretrain * Lexp[1][3] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=2: " << entries_pretrain * Lexp[2][3] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=3: " << entries_pretrain * Lexp[3][3] / preLgen << "\n";
    std::cout << "expected number of events beamsetting=4: " << entries_pretrain * Lexp[4][3] / preLgen << "\n"; */
    std::cout << "expected number of entries: " << entries_pretrain << "\n";

    std::cout << "Lgen from train output: " << postLgen << "\n";
}