#include <RtypesCore.h>
#include <TCanvas.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>

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

TH1F* sig_summed;
TH1F* bg_summed;

double bg_sum_temp(Double_t x){
    return bg_summed->GetBinContent(bg_summed->FindBin(x));
}

double sig_sum_temp(Double_t x){
    return sig_summed->GetBinContent(sig_summed->FindBin(x));
}

double fit_sum_func(double *x, double *par){
    //std::cout << "fit_func called\n";
    return bg_sum_temp(x[0]) + par[0] * sig_sum_temp(x[0]);
}

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

void get_traintree_entries(TTree *bdt_tree, TH1F* bdt_hists[], int beamsetting){
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
            bdt_hists[indexproc]->Fill(BDTresp, Lexp[beamsetting][(int)std::round(Ipol)]/Lgen);
        }
        else{
            bdt_hists[10]->Fill(BDTresp, Lexp[beamsetting][(int)std::round(Ipol)]/Lgen);
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
void find_limit_with_loop3_verify(std::string mass_point = "80"){
    std::string procId[10] = {"qqqq", "qqtt", "qqll", "qqvv", "qqlv", "qqtv", "ttll", "tttt", "qq", "qqllvv"};
    TChain* bdt_tree = new TChain;
    std::string path0 = "../training_outcome_";
    path0 += mass_point;
    TString rfile0 = path0;
    rfile0 += "_eLpR/train_bdt_qqll.root";
    path0 += "_eLpR/train_bdt_qqll.root?#dataset/TrainTree";
    char* char_path0 = new char[path0.length() + 1];
    strcpy(char_path0, path0.c_str());

    std::string path1 = "../training_outcome_";
    path1 += mass_point;
    path1 += "_eRpL/train_bdt_qqll.root?#dataset/TrainTree";
    TString rfile1 = path1;
    rfile1 += "_eRpL/train_bdt_qqll.root";
    char* char_path1 = new char[path1.length() + 1];
    strcpy(char_path1, path1.c_str());

    std::string path2 = "../training_outcome_";
    path2 += mass_point;
    TString rfile2 = path2;
    rfile2 += "_eLpL/train_bdt_qqll.root";
    path2 += "_eLpL/train_bdt_qqll.root?#dataset/TrainTree";
    char* char_path2 = new char[path2.length() + 1];
    strcpy(char_path2, path2.c_str());

    std::string path3 = "../training_outcome_";
    path3 += mass_point;
    TString rfile3 = path3;
    rfile3 += "_eRpR/train_bdt_qqll.root";
    path3 += "_eRpR/train_bdt_qqll.root?#dataset/TrainTree";
    char* char_path3 = new char[path3.length() + 1];
    strcpy(char_path3, path3.c_str());

    bdt_tree->Add(char_path0);
    bdt_tree->Add(char_path1);
    bdt_tree->Add(char_path2);
    bdt_tree->Add(char_path3);

    TFile* fileLR = new TFile(rfile0);
    TFile* fileRL = new TFile(rfile1);
    TFile* fileLL = new TFile(rfile2);
    TFile* fileRR = new TFile(rfile3);
    TTree* treeLR = (TTree*)fileLR->Get("dataset/TrainTree");
    TTree* treeRL = (TTree*)fileRL->Get("dataset/TrainTree");
    TTree* treeLL = (TTree*)fileLL->Get("dataset/TrainTree");
    TTree* treeRR = (TTree*)fileRR->Get("dataset/TrainTree");




    TH1F* bdt_hists[11];
    bdt_hists[10] = new TH1F("bdt_dist","BDT response distribution",100,-.6,.4);
    for (int i = 0; i <  sizeof(procId) / sizeof(std::string); i++){
        TString name = "bdt_dist " + procId[i];
        TString title = "BDT response distribution for " + procId[i];
        bdt_hists[i] = new TH1F(name,title,100,-.6,.4);
    }

    get_traintree_entries(treeLR, bdt_hists, 0);
    get_traintree_entries(treeRL, bdt_hists, 1);
    get_traintree_entries(treeLL, bdt_hists, 2);
    get_traintree_entries(treeRR, bdt_hists, 3);

    for(int i = 0; i < sizeof(bdt_hists) / sizeof(TH1F*) - 1; i++){
        //if(bdt_hists[i]->Integral()>0) gen_bg->FillRandom(bdt_hists[i], bdt_hists[i]->Integral());
        bg_template->Add(bdt_hists[i]);
    }
    *sig_template = *bdt_hists[10];
    //std::cout<< sizeof(bdt_hists) / sizeof(TH1F*) << "\n";
    /*
    TCanvas* c = new TCanvas;
    c->SetLogy();
    for(int i = 0; i < sizeof(bdt_hists) / sizeof(TH1F*); i++){
        bdt_hists[i]->GetYaxis()->SetRangeUser(1, 5000);
        bdt_hists[i]->Draw("SAME");
    }
    */
    for(int i = 0; i < bg_template->GetNbinsX()+1; i++){
        std::cout << sig_template->GetBinContent(i) << " ";
    }
    std::cout << "\n";

    std::vector<float> sig_bin_vec, bg_bin_vec;

    sum_low_bins(sig_template, bg_template, &sig_bin_vec, &bg_bin_vec);
    double sum_sig_vec = 0;
    for(int i = 0; i < sig_bin_vec.size(); i++){
        sum_sig_vec += sig_bin_vec.at(i);
    }
    std::cout << sig_template->GetEntries() << " " << sig_template->Integral() << " " << sum_sig_vec << "\n";
    if(sig_bin_vec.size() != bg_bin_vec.size()){
        std::cout << "signal and background have different binning!!!\n";
        return;
    }
    sig_summed = new TH1F("sig_sum", "Signal summed temp", sig_bin_vec.size(), 0, sig_bin_vec.size());
    bg_summed = new TH1F("bg_sum", "Background summed temp", bg_bin_vec.size(), 0, bg_bin_vec.size());

    for(int i = 0; i < sig_bin_vec.size(); i++){
        sig_summed->SetBinContent(i+1, sig_bin_vec.at(i));
        bg_summed->SetBinContent(i+1, bg_bin_vec.at(i));
    }
    TCanvas* c_summed = new TCanvas;
    sig_summed->Draw();
    bg_summed->SetLineColor(kRed);
    bg_summed->Draw("SAME");
    c_summed->SetLogy();
    c_summed->SaveAs("summ_hists.png");
    Double_t mean_sig = 0;
    Double_t mean_sig_err = 0;
    Double_t alpha_50_new_m = 0.00261137;
    Double_t alpha_test = 0.000445;
    Double_t alpha = 0.00454707;
    std::ifstream get_alpha_txt;
    std::string alpha_path = "../new_calc_lim_out/alpha_val_comb_" + mass_point + ".txt";
    get_alpha_txt.open(alpha_path.c_str());
    std::string alpha_filein;
    get_alpha_txt >> alpha_filein;
    alpha = std::stod(alpha_filein);
    Int_t steps = 10000;
    TH1F* sig_dist = new TH1F("hsig", "Distribution of fitted signal parameter", 50, 0, 0);
    Double_t me_ratio = 0;

        mean_sig = 0;
        mean_sig_err = 0;
        std::cout << alpha <<"\n";
        sig_dist->Reset("ICES");
        for(int m = 0; m < steps; m++){
            
            TRandom3* rnd = new TRandom3(0);
            Int_t randbin;
            
            TH1F* gen_bg = new TH1F("gen_bg", "Generated background", bg_bin_vec.size(), 0, bg_bin_vec.size());
            TH1F* gen_sig = new TH1F("gen_sig", "Generated signal", sig_bin_vec.size(), 0, sig_bin_vec.size());
            for(int i = 0; i < sig_bin_vec.size(); i++){
                gen_sig->SetBinContent(i+1, rnd->Poisson(alpha * sig_bin_vec.at(i)));
                gen_bg->SetBinContent(i+1, rnd->Poisson(bg_bin_vec.at(i)));
            }
            

            //*bg_template = *bg_template * (1/bg_template->Integral());
            //*sig_template = *sig_template * (1/sig_template->Integral());

            TH1F* data = new TH1F("data", "Combined generated data", bg_bin_vec.size(), 0, bg_bin_vec.size());
            data->Add(gen_bg, gen_sig);


            auto func = new TF1("fitting_func",&fit_sum_func, 0, bg_bin_vec.size(),1);

            func->SetParNames("Signal");

            //hist_to_fit.at(0)->GetListOfFunctions()->ls();
            auto fitresult = data->Fit(func, "SQL");
            mean_sig += func->GetParameter(0) / steps;
            mean_sig_err += func->GetParError(0) / steps;
            sig_dist->Fill(func->GetParameter(0));
            delete rnd;
            delete gen_bg;
            delete gen_sig;
            delete data;
            delete func;
        }
        me_ratio = mean_sig/mean_sig_err;


    std::cout << "final alpha: " << alpha << "\n";
    std::cout << "mean signal: " << mean_sig << " mean error: " << mean_sig_err << "ratio: " << mean_sig/mean_sig_err << "\n";
    std::cout << "integral from 0: " << 1 - sig_dist->Integral(1,sig_dist->FindBin(0) + 1)/sig_dist->Integral() << "\n";

    TCanvas* csig = new TCanvas;
    sig_dist->Draw();
    csig->SaveAs("signal_dist_new.png");

    std::ofstream file;

    file.open("alpha_val_new.txt");
    file << alpha;
    file.close();
    /*
    TCanvas* c1 = new TCanvas;
    gen_bg->SetLineColor(kRed);
    gen_bg->Draw();
    gen_sig->Draw("SAME");
    */
    //bdt_dist->Draw();
}