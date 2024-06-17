#include "TH1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"

#include <RtypesCore.h>
#include <TPad.h>
#include <iostream>

void draw_limit(){
    Double_t masses[5] = {30, 50, 80, 95, 110};
    Double_t alpha_val[5] = {0.00063, 0.00043, 0.00077, 0.00225, 0.00088};
    Double_t gen_cros[2][5] = {{169.49, 140.24, 105.08, 87.47, 69.407}, {108.62, 89.824, 67.365, 56.143, 44.5}};
    Double_t avg_cros_lim[5];

    for (int i = 0; i < 5; i++){
        avg_cros_lim[i] = (gen_cros[0][i] + gen_cros[1][i]) * alpha_val[i]/2;
    }
    TGraph *g = new TGraph(5,masses,avg_cros_lim);
    g->SetTitle("95\% limit on S->bbll cross section");
    g->SetMarkerSize(1);
    g->SetMarkerStyle(8);
    g->GetXaxis()->SetTitle("Mass [GeV]");
    g->GetYaxis()->SetTitle("95\% limit on sigma(S->bbll) [fb]");

    TLegend* legend = new TLegend(0.1, 0.8, 0.4, 0.9);
    legend->AddEntry(g, "unpolarized");
    TCanvas* c = new TCanvas;
    
    g->Draw("APL");
    legend->Draw();

    //c->DrawClone();
    c->SaveAs("limit_graph.png");

}

/*# Mass  bbvv_eLpR   bbvv_eRpL   bbll_eLpR   bbll_eRpL   bbqq_eLpR   bbqq_eRpL  
    30     335.85      215.26      169.49      108.62      1145.5      733.07  
    50     277.91      178.17      140.24      89.824      945.93       605.9  
    80     208.55      133.78      105.08      67.365      710.23      454.45  
    95      173.6      111.02       87.47      56.143      590.77      378.02  
   110     137.67      88.288      69.407        44.5      468.84      299.66  
*/