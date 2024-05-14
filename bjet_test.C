/*

Check consistency of event reconstruction by comparing total energy of
all energy flow objects with energy of objects reconstructed at the top
level (jets, isolated electrons, muons and photons).

root -l 'rec_check.C++("events/bbll_sig_50_eLpR.root","Jet_N2")'

*/

#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TGraph.h"
#include "TLegend.h"

#include <iostream>
#include <fstream>    
#include <iomanip>
using namespace std;


//#ifdef __CLING__
R__LOAD_LIBRARY(delphes/libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
//#else
//class ExRootTreeReader;
//class ExRootResult;
//#endif


void bjet_test(int Bbit = 2, const char *genFile="bbvv_sig_110_eRpL.root",
	       const char *jetName="Jet_N2",
	       int nbin=300, double emin=-30., double emax=270.)
{
  // (re)book control histrogram

  if( gDirectory->FindObject("heflow")!=NULL) gDirectory->Delete("heflow");
  TH1F *heflow = new TH1F("heflow","Total energy of energy flow objects",nbin, emin, emax);

  if( gDirectory->FindObject("hetot")!=NULL) gDirectory->Delete("hetot");
  TH1F *hetot = new TH1F("hetot","Total energy of final state objects",nbin, emin, emax);

    
  if( gDirectory->FindObject("hedif")!=NULL) gDirectory->Delete("hedif");
  TH1F *hedif = new TH1F("hedif","Energy 'lost' in reconstruction",nbin, emin, -emin);

  if( gDirectory->FindObject("hbdif")!=NULL) gDirectory->Delete("hbdif");
  TH1F *hbdif = new TH1F("hbdif","Energy difference between bquark and reconstructed jet",nbin, -120, 120);

  if( gDirectory->FindObject("hjcnt")!=NULL) gDirectory->Delete("hjcnt");
  TH1I *hjcnt = new TH1I("hjcnt","Number of electrons in a jet",20, 0, 20);

    if( gDirectory->FindObject("hbdifnolep")!=NULL) gDirectory->Delete("hbdifnolep");
  TH1F *hbdifnolep = new TH1F("hbdifnolep","Energy difference between bquark and reconstructed jet not containg electrons",nbin, -120, 120);

  if( gDirectory->FindObject("hbdiflep")!=NULL) gDirectory->Delete("hbdiflep");
  TH1F *hbdiflep = new TH1F("hbdiflep","Energy difference between bquark and reconstructed jet containing electrons",nbin, -120, 120);


  gSystem->Load("delphes/libDelphes");
  
  // Create chain of input root trees
  
  TChain chain("Delphes");
  chain.Add(genFile);

  // Create object of class ExRootTreeReader

  ExRootTreeReader *genReader = new ExRootTreeReader(&chain);
  Long64_t genEntries = genReader->GetEntries();

  if(genEntries <=0)
    {
      cout << "No events in input file " << genFile << " ?" << endl;
      return;
    }
  
  cout << "Input chain contains " << genEntries << " events" << endl;


  // Get pointers to branches used in this analysis
  
  
  // Final reconstruction results

  TClonesArray *branchJet = genReader->UseBranch(jetName);
  Jet *jet;
  Jet *jet1;
  Jet *jet2;
  
  TClonesArray *branchParticle = genReader->UseBranch("Particle");
  GenParticle *particle;

  GenParticle *bquark;
  GenParticle *antibquark;

  Double_t energyDiffAntiB = 0;
  Double_t energyDiffB = 0;

  Double_t quarkE[200000];
  Double_t jetE[200000];
  int _nEvt=0;
  
  // Loop over tree events
  
  for(Int_t entry = 0; entry < genEntries; ++entry)
  {
    bool debug =  entry<10 || entry%5000==0;
    genReader->ReadEntry(entry);
    TLorentzVector EFlowSum;
    EFlowSum.SetPxPyPzE(0.,0.,0.,0.);
    TLorentzVector ERecSum;
    ERecSum.SetPxPyPzE(0.,0.,0.,0.);

    // Jets
    
    int nJet = branchJet->GetEntriesFast();
    int nBtag = 0;
    for(Int_t ijet=0; ijet<nJet; ijet++)
      {
      jet = (Jet*) branchJet->At(ijet);

      //TLorentzVector v = jet->P4() ;
      //if(debug) cout << "\n";
      if((jet->BTag & Bbit)&& nBtag==1){
        jet1 = jet;
        nBtag++;
      }
      if((jet->BTag & Bbit)&&nBtag==0){
        jet2 = jet;
        nBtag++;
      }
      //if(debug) cout << entry << " -> " << v.E() << "(j)\n"  ;
      //ERecSum+=v;
      }
    


    if(nBtag != 2) continue;

    int nParticle = branchParticle->GetEntriesFast();
    for(Int_t iparticle = 0; iparticle<nParticle; iparticle++){
      particle = (GenParticle *)branchParticle->At(iparticle);
      if(particle->PID == 5) bquark = particle;
      if(particle->PID == -5) antibquark = particle;
      /*if(entry == 1){
        std::cout << particle->PID << "\n";
      }*/
    }

    TLorentzVector v1 = jet1->P4();
    TLorentzVector v2 = jet2->P4();
    TVector3 jet1p(v1.Px(), v1.Py(), v1.Pz());
    TVector3 jet2p(v2.Px(), v2.Py(), v2.Pz());
    TVector3 quark1p(bquark->Px, bquark->Py, bquark->Pz);
    TVector3 quark2p(antibquark->Px, antibquark->Py, antibquark->Pz);

    Double_t anglediff1 = jet1p.Angle(quark1p) * jet1p.Angle(quark1p) + jet2p.Angle(quark2p) * jet2p.Angle(quark2p);
    Double_t anglediff2 = jet1p.Angle(quark2p) * jet1p.Angle(quark2p) + jet2p.Angle(quark1p) * jet2p.Angle(quark1p);

    if(anglediff2 < anglediff1){
      jet = jet1;
      jet1 = jet2;
      jet2 = jet;
      v1 = jet1->P4();
      v2 = jet2->P4();
    }
    Int_t jet1_electron_count = 0;
    Int_t jet2_electron_count = 0;
    TRefArray jet1particles = jet1->Particles;
    TRefArray jet2particles = jet2->Particles;
    //cout << "Jet1 particles:\n";
    for(int i = 0; i < jet1particles.GetSize(); i++){
      GenParticle *a;
      a=(GenParticle *)jet1particles[i];
      if(a != nullptr && (a->PID == 11 || a->PID == -11)) jet1_electron_count++;
    }
    for(int i = 0; i < jet2particles.GetSize(); i++){
      GenParticle *a;
      a=(GenParticle *)jet2particles[i];
      if(a != nullptr && (a->PID == 11 || a->PID == -11)) jet2_electron_count++;
    }
  
    energyDiffB = v1.E() - bquark->E;
    energyDiffAntiB = v2.E() - antibquark->E;

    quarkE[2*entry] = bquark->E;
    quarkE[2*entry + 1] = antibquark->E;
    
    jetE[2*entry] = v1.E();
    jetE[2*entry+1] = v2.E();

    if(debug){
      cout << entry << " Quark energ: " << quarkE[2*entry+1] << " Jet energ: " << jetE[2*entry+1] << "\n";
    }
    hbdif->Fill(energyDiffAntiB);
    hbdif->Fill(energyDiffB);
    //heflow->Fill(EFlowSum.E());
    //hetot->Fill(ERecSum.E());
    //hedif->Fill(EFlowSum.E()-ERecSum.E());
    hjcnt->Fill(jet1_electron_count);
    hjcnt->Fill(jet2_electron_count);
    if(jet1_electron_count == 0)
    {
      hbdifnolep->Fill(energyDiffB);
    }
    else
    {
      hbdiflep->Fill(energyDiffB);
    }
    if(jet2_electron_count == 0)
    {
      hbdifnolep->Fill(energyDiffAntiB);
    }
    else
    {
      hbdiflep->Fill(energyDiffAntiB);
    }
    _nEvt++;
  }
  TCanvas *c1 = new TCanvas;
  TGraph* g = new TGraph(200000, quarkE, jetE);
  TAxis* xaxis = g->GetXaxis();
  xaxis->SetLimits(0., 120.);
  g->SetTitle("Correlation between quark energy and reconstructed jet energy;B Quark Energy; B Jet Energy");
  g->SetMarkerStyle(1);
  g->Draw("ap");
  c1->SaveAs("energ_corr_graph.png");
  
  TCanvas *c2 = new TCanvas;
  hbdif->Draw();
  c2->SaveAs("Bjet_reconstruct_energy_diff.png");

  TCanvas *c3 = new TCanvas;
  hjcnt->Draw();
  c3->SaveAs("number_of_electrons_per_jet.png");

  TCanvas *c4 = new TCanvas;
  auto legend = new TLegend(0.1,0.8,0.3,0.9);
  legend->AddEntry(hbdiflep, "With electrons", "l");
  legend->AddEntry(hbdifnolep, "Without electrons","l");
  
  hbdiflep->Scale(1./hbdiflep->Integral());
  hbdifnolep->Scale(1./hbdifnolep->Integral());
  
  hbdifnolep->SetLineColor(2);
  hbdifnolep->Draw("hist");
  hbdiflep->Draw("hist same");
  legend->Draw();
  c4->SaveAs("energdiff_lep_nolep.png");


  // Open canvas
/*
   gStyle->SetPadBottomMargin(0.10);
   gStyle->SetPadTopMargin(0.08);
   gStyle->SetPadRightMargin(0.05);
   gStyle->SetPadLeftMargin(0.15);
   
  
  TCanvas  *ch1 = (TCanvas *) gROOT->FindObject("ch1");
  if(!ch1)  ch1 = new TCanvas("ch1","Delphes histograms",20,20,900,450);

  ch1->Clear();
  ch1->Divide(2,1);
 
  ch1->cd(1);
  heflow->SetLineWidth(2);
  heflow->SetLineColor(3);
  heflow->GetXaxis()->SetTitle("E [GeV]");
  heflow->GetYaxis()->SetTitle("Entries");
  heflow->SetMinimum(0.7);
  heflow->Draw();
  
  hetot->SetLineColor(4);
  hetot->Draw("same");
 
  ch1->cd(2);
  hedif->SetLineColor(2);
  hedif->GetXaxis()->SetTitle("E_{Eflow} - E_{rec} [GeV]");
  hedif->GetYaxis()->SetTitle("Entries");
  hedif->SetMinimum(0.7);
  
  hedif->Draw();
*/
}
