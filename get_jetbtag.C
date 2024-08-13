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

#include <TFile.h>
#include <iostream>
#include <fstream>    
#include <iomanip>
#include <string>
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


void get_jetbtag(int Bbit = 2, const char *genFile="bbll_sig_80_eRpL.root",
	       bool isBg = false, const char *jetName="Jet_N2",
	       int nbin=300, double emin=-30., double emax=270.)
{
  // (re)book control histrogram




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

  TFile* jetbtagfile = new TFile("jetbtag.root", "RECREATE");
  
  TTree* eventtree
  
  // Final reconstruction results

  TClonesArray *branchJet = genReader->UseBranch(jetName);
  Jet *jet;
  Jet *jet1;
  Jet *jet2;


  TClonesArray *branchParticle = genReader->UseBranch("Particle");
  GenParticle *particle;

  GenParticle *bquark;
  GenParticle *antibquark;

  int _nEvt=0;

  int n_phi_pi = 0;
  
  // Loop over tree events
  
  for(Int_t entry = 0; entry < genEntries; ++entry)
  {
    bool debug =  entry<10 || entry%5000==0;
    genReader->ReadEntry(entry);

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

    //muons
    int nMuons = branchMuon->GetEntriesFast();
    int nElectrons = branchElectron->GetEntriesFast();

    if(!(nMuons == 2 || nElectrons == 2)) continue;

    Double_t pllx, plly;

    TLorentzVector lep1P4, lep2P4;


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

    _nEvt++;
  }

  string path = "bbll_80_eRpL";

  if(isBg) path = "bg_qqll_eRpL";

  //const char* cpath = path.c_str();


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
