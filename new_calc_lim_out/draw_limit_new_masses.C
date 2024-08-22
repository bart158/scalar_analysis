#include "TH1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"

#include <Rtypes.h>
#include <RtypesCore.h>
#include <TMVA/Factory.h>
#include <TPad.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <istream>
#include <string>

void draw_limit_new_masses(){
    Double_t masses[12] = {20, 30, 40, 50, 60, 70, 80, 90, 95, 100, 110, 120};
    Double_t alpha_val_unpol[12] = {0.00036, 0.000625, 0.000575, 0.00043, 0.000557, 0.00054, 0.00077, 0.00255, 0.00225, 0.00155, 0.00088, 0.00077};
    Double_t alpha_val_pol[5][12];
    Double_t alpha_val_comb[12];
    std::string polarizations[5] = {"eLpR", "eRpL", "eLpL", "eRpR", "unpol"};
    
    std::ifstream file;
    std::string line;
    for(int i = 0; i < 5; i++){
        for(int j = 0; j < 12; j++){
            std::string path;
            path = "alpha_val_" + polarizations[i] + "_" + std::to_string((int)masses[j]) + ".txt";
            std::cout << path << "\n";
            file.open(path);
            getline(file, line);
            std::cout << line << "\n";
            alpha_val_pol[i][j] = std::stod(line);
            file.close();
        }
    }
    for(int j = 0; j < 12; j++){
        std::string path;
        path = "alpha_val_comb_" + std::to_string((int)masses[j]) + ".txt";
        std::cout << path << "\n";
        file.open(path);
        getline(file, line);
        std::cout << line << "\n";
        alpha_val_comb[j] = std::stod(line);
        file.close();
    }
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Alpha 95\% CL;mass [GeV];alpha val");
    TGraph *g = new TGraph(12,masses,alpha_val_pol[4]);
    TGraph *g0 = new TGraph(12,masses,alpha_val_pol[0]);
    TGraph *g1 = new TGraph(12,masses,alpha_val_pol[1]);
    TGraph *g2 = new TGraph(12,masses,alpha_val_pol[2]);
    TGraph *g3 = new TGraph(12,masses,alpha_val_pol[3]);
    TGraph *g4 = new TGraph(12, masses, alpha_val_comb);

    g->SetTitle("95\% limit on S->bbll cross section");
    g->SetMarkerSize(1);
    g->SetMarkerStyle(8);
    g->GetXaxis()->SetTitle("Mass [GeV]");
    g->GetYaxis()->SetTitle("95\% limit on sigma(S->bbll) [fb]");

    g0->SetMarkerSize(1);
    g0->SetMarkerStyle(8);
    g1->SetMarkerSize(1);
    g1->SetMarkerStyle(8);
    g2->SetMarkerSize(1);
    g2->SetMarkerStyle(8);
    g3->SetMarkerSize(1);
    g3->SetMarkerStyle(8);
    g4->SetMarkerSize(1);
    g4->SetMarkerStyle(8);

    g0->SetLineColor(kRed);
    g1->SetLineColor(kBlue);
    g2->SetLineColor(kGreen);
    g3->SetLineColor(kOrange);
    g4->SetLineColor(kViolet);

    TLegend* legend = new TLegend(0.1, 0.7, 0.3, 0.9);
    legend->AddEntry(g, "unpolarized");
    legend->AddEntry(g0, "LR");
    legend->AddEntry(g1, "RL");
    legend->AddEntry(g2, "LL");
    legend->AddEntry(g3, "RR");
    legend->AddEntry(g4, "comb.");
    TCanvas* c = new TCanvas;
    c->SetLogy(true);
    
    mg->Add(g);
    mg->Add(g0);
    mg->Add(g1);
    mg->Add(g2);
    mg->Add(g3);
    mg->Add(g4);

    mg->Draw("APL");
    legend->Draw();

    //c->DrawClone();
    c->SaveAs("limit_graph_new_masses_v2_org.png");

}

/*# Mass  bbvv_eLpR   bbvv_eRpL   bbll_eLpR   bbll_eRpL   bbqq_eLpR   bbqq_eRpL  
    30     335.85      215.26      169.49      108.62      1145.5      733.07  
    50     277.91      178.17      140.24      89.824      945.93       605.9  
    80     208.55      133.78      105.08      67.365      710.23      454.45  
    95      173.6      111.02       87.47      56.143      590.77      378.02  
   110     137.67      88.288      69.407        44.5      468.84      299.66  
*/