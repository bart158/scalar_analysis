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

void find_limit_with_loop2(std::string mass_point){
    std::string procId[10] = {"qqqq", "qqtt", "qqll", "qqvv", "qqlv", "qqtv", "ttll", "tttt", "qq", "qqllvv"};
    TChain* bdt_tree = new TChain;
    std::string path0 = "../training_outcome_";
    path0 += mass_point;
    path0 += "_eLpR/train_bdt_qqll.root?#dataset/TrainTree";
    char* char_path0 = new char[path0.length() + 1];
    strcpy(char_path0, path0.c_str());

    std::string path1 = "../training_outcome_";
    path1 += mass_point;
    path1 += "_eRpL/train_bdt_qqll.root?#dataset/TrainTree";
    char* char_path1 = new char[path1.length() + 1];
    strcpy(char_path1, path1.c_str());

    std::string path2 = "../training_outcome_";
    path2 += mass_point;
    path2 += "_eLpL/train_bdt_qqll.root?#dataset/TrainTree";
    char* char_path2 = new char[path2.length() + 1];
    strcpy(char_path2, path2.c_str());

    std::string path3 = "../training_outcome_";
    path3 += mass_point;
    path3 += "_eRpR/train_bdt_qqll.root?#dataset/TrainTree";
    char* char_path3 = new char[path3.length() + 1];
    strcpy(char_path3, path3.c_str());

    bdt_tree->Add(char_path0);
    bdt_tree->Add(char_path1);
    bdt_tree->Add(char_path2);
    bdt_tree->Add(char_path3);


    Float_t BDTresp;
    Float_t Iproc;
    bdt_tree->SetBranchAddress("BDT", &BDTresp);
    bdt_tree->SetBranchAddress("Iproc", &Iproc);
    Long64_t genEntries = bdt_tree->GetEntries();
    TH1F* bdt_hists[11];
    bdt_hists[10] = new TH1F("bdt_dist","BDT response distribution",100,-.6,.4);
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
            bdt_hists[10]->Fill(BDTresp);
        }
        //bdt_dist->Fill(BDTresp);
    }

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
    int ibin = 0;
    int first_bg_bins_sum = 1;
    int first_sig_bins_sum = 1;
    double sec_der_sum = 0;
    int bg_bins[100];
    int sig_bins[100];
    while(first_bg_bins_sum < 10){
        first_bg_bins_sum+=bg_template->GetBinContent(ibin);
        first_sig_bins_sum+=sig_template->GetBinContent(ibin);
        ibin++;
    }
    sec_der_sum+=first_sig_bins_sum*first_sig_bins_sum/(double)first_bg_bins_sum;
    int last_bg_bins_sum = 0;
    int last_sig_bins_sum = 0;
    int jbin = sig_template->GetNbinsX();
    while(last_bg_bins_sum < 10){
        last_bg_bins_sum+=bg_template->GetBinContent(jbin);
        last_sig_bins_sum+=sig_template->GetBinContent(jbin);
        jbin--;
    }
    sec_der_sum+=last_sig_bins_sum*last_sig_bins_sum/(double)last_bg_bins_sum;
    std::cout << first_sig_bins_sum <<" ";
    double sig;
    for(int i = ibin; i < jbin+1; i++){
        sig = sig_template->GetBinContent(i);
        sec_der_sum+=sig*sig/bg_template->GetBinContent(i);
        bg_bins[i-ibin] = bg_template->GetBinContent(i);
        sig_bins[i-ibin] = sig;
        std::cout << sig_bins[i-ibin] <<" ";
    }
    std::cout << last_sig_bins_sum <<"\n";
    double alpha_org = 1.64 * TMath::Sqrt(1/sec_der_sum);
    std::cout << "alpha from the first method: " << alpha_org << "\n";

    Double_t mean_sig = 0;
    Double_t mean_sig_err = 0;
    Double_t alpha = alpha_org;
    Int_t steps = 10000;
    TH1F* sig_dist = new TH1F("hsig", "Distribution of fitted signal parameter", 50, -.002, .003);
    Double_t me_ratio = 0;
    Double_t a = 0;
    Double_t b = 0.005;
    while(me_ratio < 1.60 || me_ratio > 1.68){
        mean_sig = 0;
        mean_sig_err = 0;
        alpha = (a+b)/2;
        std::cout << alpha <<"\n";
        sig_dist->Reset("ICES");
        for(int m = 0; m < steps; m++){
            
            TRandom3* rnd = new TRandom3(0);
            Int_t randbin;
            
            TH1F* gen_bg = new TH1F("gen_bg", "Generated background", 100, -.6, .4);
            TH1F* gen_sig = new TH1F("gen_sig", "Generated signal", 100, -.6, .4);

            for(int i = 0; i < 100; i++){
                randbin = rnd->Poisson(alpha * bdt_hists[10]->GetBinContent(i+1));
                    for(int k = 0; k < randbin; k++){
                        gen_sig->Fill(0.01*i - 0.599);
                    }
            }

            for(int i = 0; i < sizeof(bdt_hists) / sizeof(TH1F*) - 1; i++){
                //if(bdt_hists[i]->Integral()>0) gen_bg->FillRandom(bdt_hists[i], bdt_hists[i]->Integral());
                for(int j = 0; j < 100; j++){
                    randbin = rnd->Poisson(bdt_hists[i]->GetBinContent(j+1));
                    for(int k = 0; k < randbin; k++){
                        gen_bg->Fill(0.01*j - 0.599);
                    }
                }
            }
            

            //*bg_template = *bg_template * (1/bg_template->Integral());
            //*sig_template = *sig_template * (1/sig_template->Integral());

            TH1F* data = new TH1F("data", "Combined generated data", 100, -0.6, 0.4);
            data->Add(gen_bg, gen_sig);


            auto func = new TF1("fitting_func",&fit_func,-0.6,0.4,2);

            func->SetParNames("Background","Signal");

            //hist_to_fit.at(0)->GetListOfFunctions()->ls();
            auto fitresult = data->Fit(func, "SQ");
            mean_sig += func->GetParameter(1) / steps;
            mean_sig_err += func->GetParError(1) / steps;
            sig_dist->Fill(func->GetParameter(1));
            delete rnd;
            delete gen_bg;
            delete gen_sig;
            delete data;
            delete func;
        }
        me_ratio = mean_sig/mean_sig_err;
        if(me_ratio < 1.6){
            a = alpha;
        }
        else if(me_ratio > 1.68){
            b = alpha;
        }
    }

    std::cout << "final alpha: " << alpha << "\n";
    std::cout << "mean signal: " << mean_sig << " mean error: " << mean_sig_err << "ratio: " << mean_sig/mean_sig_err << "\n";
    std::cout << "integral from 0: " << 1 - sig_dist->Integral(1,sig_dist->FindBin(0) + 1)/sig_dist->Integral() << "\n";

    TCanvas* csig = new TCanvas;
    sig_dist->Draw();
    csig->SaveAs("signal_dist.png");

    std::ofstream file;

    file.open("alpha_val.txt");
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