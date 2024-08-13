#include <string.h>
#include <iostream>
#include <math.h>
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"


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
void calc_lim(){
    std::string pathLR = "training_output_eLpR/classification_output.root";
    std::string pathRL = "training_output_eRpL/classification_output.root";
    std::string pathLL = "training_output_eLpL/classification_output.root";
    std::string pathRR = "training_output_eRpR/classification_output.root";

    TFile* fileLR = new TFile(pathLR.c_str(), "OPEN");
    TFile* fileRL = new TFile(pathRL.c_str(), "OPEN");
    TFile* fileLL = new TFile(pathLL.c_str(), "OPEN");
    TFile* fileRR = new TFile(pathRR.c_str(), "OPEN");

    std::vector<TH1F*> sig_hists;
    std::vector<TH1F*> bg_hists;

    sig_hists.push_back((TH1F*)fileLR->FindObjectAny("BDTSignalResponse"));
    sig_hists.push_back((TH1F*)fileRL->FindObjectAny("BDTSignalResponse"));
    sig_hists.push_back((TH1F*)fileLL->FindObjectAny("BDTSignalResponse"));
    sig_hists.push_back((TH1F*)fileRR->FindObjectAny("BDTSignalResponse"));

    bg_hists.push_back((TH1F*)fileLR->FindObjectAny("BDTBackgroundResponse"));
    bg_hists.push_back((TH1F*)fileRL->FindObjectAny("BDTBackgroundResponse"));
    bg_hists.push_back((TH1F*)fileLL->FindObjectAny("BDTBackgroundResponse"));
    bg_hists.push_back((TH1F*)fileRR->FindObjectAny("BDTBackgroundResponse"));

    std::vector<std::vector<float>> sig_mod_pol(4);
    std::vector<std::vector<float>> bg_mod_pol(4);
    std::vector<double> sec_der_pol{0, 0, 0, 0};
    std::vector<double> alpha_lim_pol(4);
    for(int i = 0; i < sig_hists.size(); i++){
        sum_low_bins(sig_hists.at(i), bg_hists.at(i), &(sig_mod_pol.at(i)), &(bg_mod_pol.at(i)));
        for(int j = 0; j < sig_mod_pol.at(i).size(); j++){
            if(sig_mod_pol.at(i).at(j) == 0 && bg_mod_pol.at(i).at(j) == 0) continue;
            sec_der_pol.at(i) += ((double)sig_mod_pol.at(i).at(j) * (double)sig_mod_pol.at(i).at(j))/(double)bg_mod_pol.at(i).at(j);
        }
        alpha_lim_pol.at(i) = 1.64 * sqrt(1.0/sec_der_pol.at(i));
    }

    TH1F* sig_comb = new TH1F("h_sig_comb", "Combined signal response for all polarizations", 100, -0.6, 0.4);
    TH1F* bg_comb = new TH1F("h_bg_comb", "Combined background response for all polarizations", 100, -0.6, 0.4);

    for(int i = 0; i < 4; i++){
        *sig_comb = *sig_comb + *(sig_hists.at(i));
        *bg_comb = *bg_comb + *(bg_hists.at(i));
    }
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
    for(int i = 0; i < alpha_lim_pol.size(); i++){
        if(i == 0) {pol_sig = "eLpR";}
        else if(i == 1) {pol_sig = "eRpL";}
        else if(i == 2) {pol_sig = "eLpL";}
        else {pol_sig = "eRpR";}
        std::cout << "*---------------------------*\n";
        std::cout << "95\% CL on alpha for " << pol_sig << ": " << alpha_lim_pol.at(i) << "\n";
        std::cout << "*---------------------------*\n";
    }
    TCanvas c;
    bg_comb->Draw();
    sig_comb->Draw("same");
    c.SetLogy(1);
    c.SaveAs("output_hists_bdt.png");
}
