#include <string.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TChain.h"

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
    return bg_temp(x[0]) + par[0] * sig_temp(x[0]);
}

void get_testtree_entries(TTree *bdt_tree, TH1F* bdt_hists[], int beamsetting){
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

    Float_t BDTresp;
    Float_t Iproc;
    Float_t Ipol;
    Float_t Lgen;

    bdt_tree->SetBranchAddress("BDT", &BDTresp);
    bdt_tree->SetBranchAddress("Iproc", &Iproc);
    bdt_tree->SetBranchAddress("Ipol", &Ipol);
    bdt_tree->SetBranchAddress("Lgen", &Lgen);

    Long64_t genEntries = bdt_tree->GetEntries();

    for(Int_t entry = 0; entry < genEntries; ++entry){
        bdt_tree->GetEntry(entry);
        Int_t indexproc = std::round(Iproc);
        if (indexproc >= 0){
            bdt_hists[indexproc]->Fill(BDTresp, 2*Lexp[beamsetting][(int)std::round(Ipol)]/Lgen);
        }
        else{
            bdt_hists[10]->Fill(BDTresp, 2*Lexp[beamsetting][(int)std::round(Ipol)]/Lgen);
        }
        //bdt_dist->Fill(BDTresp);
    }
}

void sum_low_bins(TH1F* sig_hist, TH1F* bg_hist, std::vector<float>* sig_mod, std::vector<float>* bg_mod){
    for(int i = 100; i > 0; i--){
        //std::cout << sig_comb->GetBinContent(i) << " " << bg_comb->GetBinContent(i) << "\n";
        if(sig_mod->empty()){
            sig_mod->push_back(sig_hist->GetBinContent(i));
            bg_mod->push_back(bg_hist->GetBinContent(i));
            continue;
        }
        if(bg_hist->GetBinContent(i+1) < 300){
            sig_mod->back() += sig_hist->GetBinContent(i);
            bg_mod->back() += bg_hist->GetBinContent(i);
        }
        else{
            sig_mod->push_back(sig_hist->GetBinContent(i));
            bg_mod->push_back(bg_hist->GetBinContent(i));
        }
    }
}
void calc_lim_v2_with_bmst(std::string mass_point = "80", int bmst = 4){
    std::string procId[10] = {"qqqq", "qqtt", "qqll", "qqvv", "qqlv", "qqtv", "ttll", "tttt", "qq", "qqllvv"};
    TChain* bdt_tree = new TChain;
    std::string path0 = "../training_outcome_";
    path0 += mass_point;
    TString rfile0 = path0;
    rfile0 += "_eLpR/train_bdt_qqll.root";
    path0 += "_eLpR/train_bdt_qqll.root?#dataset/TestTree";
    char* char_path0 = new char[path0.length() + 1];
    strcpy(char_path0, path0.c_str());

    std::string path1 = "../training_outcome_";
    path1 += mass_point;
    path1 += "_eRpL/train_bdt_qqll.root?#dataset/TestTree";
    TString rfile1 = path1;
    rfile1 += "_eRpL/train_bdt_qqll.root";
    char* char_path1 = new char[path1.length() + 1];
    strcpy(char_path1, path1.c_str());

    std::string path2 = "../training_outcome_";
    path2 += mass_point;
    TString rfile2 = path2;
    rfile2 += "_eLpL/train_bdt_qqll.root";
    path2 += "_eLpL/train_bdt_qqll.root?#dataset/TestTree";
    char* char_path2 = new char[path2.length() + 1];
    strcpy(char_path2, path2.c_str());

    std::string path3 = "../training_outcome_";
    path3 += mass_point;
    TString rfile3 = path3;
    rfile3 += "_eRpR/train_bdt_qqll.root";
    path3 += "_eRpR/train_bdt_qqll.root?#dataset/TestTree";
    char* char_path3 = new char[path3.length() + 1];
    strcpy(char_path3, path3.c_str());

    std::string path4 = "../training_outcome_";
    path4 += mass_point;
    TString rfile4 = path4;
    rfile4 += "_unpol/train_bdt_qqll.root";
    path4 += "_unpol/train_bdt_qqll.root?#dataset/TestTree";
    char* char_path4 = new char[path4.length() + 1];
    strcpy(char_path4, path4.c_str());

    bdt_tree->Add(char_path0);
    bdt_tree->Add(char_path1);
    bdt_tree->Add(char_path2);
    bdt_tree->Add(char_path3);

    TFile* fileLR = new TFile(rfile0);
    TFile* fileRL = new TFile(rfile1);
    TFile* fileLL = new TFile(rfile2);
    TFile* fileRR = new TFile(rfile3);
    TFile* fileunpol = new TFile(rfile4);
    TTree* treeLR = (TTree*)fileLR->Get("dataset/TestTree");
    TTree* treeRL = (TTree*)fileRL->Get("dataset/TestTree");
    TTree* treeLL = (TTree*)fileLL->Get("dataset/TestTree");
    TTree* treeRR = (TTree*)fileRR->Get("dataset/TestTree");
    TTree* treeunpol = (TTree*)fileunpol->Get("dataset/TestTree");




    TH1F* bdt_hists[11];
    bdt_hists[10] = new TH1F("bdt_dist","BDT response distribution",100,-.6,.4);
    for (int i = 0; i <  sizeof(procId) / sizeof(std::string); i++){
        TString name = "bdt_dist " + procId[i];
        TString title = "BDT response distribution for " + procId[i];
        bdt_hists[i] = new TH1F(name,title,100,-.6,.4);
    }

    if(bmst == 0){
        get_testtree_entries(treeLR, bdt_hists, bmst);
    }
    else if(bmst == 1){
        get_testtree_entries(treeRL, bdt_hists, bmst);
    }
    else if (bmst == 2) {
        get_testtree_entries(treeLL, bdt_hists, bmst);
    }
    else if (bmst == 3) {
        get_testtree_entries(treeRR, bdt_hists, bmst);
    }
    else if (bmst == 4){
        get_testtree_entries(treeunpol, bdt_hists, bmst);
    }
    else if(bmst == 5){
        get_testtree_entries(treeLR, bdt_hists, 0);
        get_testtree_entries(treeRL, bdt_hists, 1);
        get_testtree_entries(treeLL, bdt_hists, 2);
        get_testtree_entries(treeRR, bdt_hists, 3);
    }

    for(int i = 0; i < sizeof(bdt_hists) / sizeof(TH1F*) - 1; i++){
        //if(bdt_hists[i]->Integral()>0) gen_bg->FillRandom(bdt_hists[i], bdt_hists[i]->Integral());
        bg_template->Add(bdt_hists[i]);
    }
    *sig_template = *bdt_hists[10];


    TH1F* sig_comb;
    TH1F* bg_comb;
    sig_comb = sig_template;
    bg_comb = bg_template;

    
    std::vector<float> sig_mod_bin;
    std::vector<float> bg_mod_bin;
    
    sum_low_bins(sig_comb, bg_comb, &sig_mod_bin, &bg_mod_bin);

    std::cout << "-------------------\n";
    double sec_derivative_sum = 0;
    for(int i = 0; i < sig_mod_bin.size(); i++){
        //std::cout << sig_mod_bin.at(i) << " " << bg_mod_bin.at(i) << "\n";
        sec_derivative_sum += ((double)sig_mod_bin.at(i) * (double)sig_mod_bin.at(i))/(double)bg_mod_bin.at(i);
    }

    double sigma_alpha = sqrt(1.0/sec_derivative_sum);
    double alpha_lim = 1.64 * sigma_alpha;
    std::cout << "*---------------------------*\n";
    std::cout << "95\% CL on alpha combined: " << alpha_lim << "\n";
    std::cout << "*---------------------------*\n";
    std::string pol_sig{""};
    std::string outfilename = "alpha_val_";
    switch (bmst) {
        case 0:
            outfilename+="eLpR_";
            break;
        case 1:
            outfilename+="eRpL_";
            break;
        case 2:
            outfilename+="eLpL_";
            break;
        case 3:
            outfilename+="eRpR_";
            break;
        case 4:
            outfilename+="unpol_";
            break;
        case 5:
            outfilename+="comb_";
            break;
    }
    outfilename += (std::string)mass_point;
    outfilename += ".txt";
    std::ofstream outfile;
    outfile.open(outfilename.c_str());
    outfile << alpha_lim;
    outfile.close();

    TCanvas c;
    bg_comb->Draw();
    sig_comb->Draw("same");
    c.SetLogy(1);
    c.SaveAs("output_hists_bdt.png");
}
