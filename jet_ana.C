/*

Check consistency of event reconstruction by comparing total energy of
all energy flow objects with energy of objects reconstructed at the top
level (jets, isolated electrons, muons and photons).
More information about jets (generator level particles and EFlow level objects)
is included in printout.

root -l 'jet_ana.C++("events/bbll_sig_50_eLpR.root","Jet_N2")'

*/

#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TClonesArray.h"

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


void jet_ana(const char *genFile="bbvv_sig_80_eRpL.root",
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
  
  // Generator level

  TClonesArray *branchParticle = genReader->UseBranch("Particle");
  GenParticle *particle;

  // Energy flow objects
  
  TClonesArray *branchCharged = genReader->UseBranch("EFlowTrack");
  Track *charged;
  
  TClonesArray *branchNeutral = genReader->UseBranch("EFlowNeutralHadron");
  Tower *neutral;
  
  TClonesArray *branchEFPhoton = genReader->UseBranch("EFlowPhoton");
  Tower *efphoton;
  
  // Final reconstruction results

  TClonesArray *branchJet = genReader->UseBranch(jetName);
  Jet *jet;
  
  TClonesArray *branchElectron = genReader->UseBranch("Electron");
  Electron *electron;

  TClonesArray *branchMuon = genReader->UseBranch("Muon");
  Muon *muon;

  TClonesArray *branchPhoton = genReader->UseBranch("Photon");
  Photon *photon;


  int _nEvt=0;
  
  // Loop over tree events
  
  for(Int_t entry = 0; entry < genEntries; ++entry)
  {
    bool debug =  entry<10 || entry%5000==0;
    
    // Load selected branches with data from specified event
    
    genReader->ReadEntry(entry);

    // Generator level particles
    
    int nParticle = branchParticle->GetEntriesFast();


    // Energy flow objects
    
    int nCh = branchCharged->GetEntriesFast();
    int nNeu = branchNeutral->GetEntriesFast();
    int nPh = branchEFPhoton->GetEntriesFast();

   // Fill energy sum for EFlow objects
   
    TLorentzVector EFlowSum;
    EFlowSum.SetPxPyPzE(0.,0.,0.,0.);

    for(Int_t ineu=0; ineu<nNeu; ineu++)
      {
      neutral = (Tower*) branchNeutral->At(ineu);

      TLorentzVector v = neutral->P4() ;
      
      EFlowSum+=v;
      }

    for(Int_t iphot=0; iphot<nPh; iphot++)
      {
      efphoton = (Tower*) branchEFPhoton->At(iphot);

      TLorentzVector v = efphoton->P4() ;
      
      EFlowSum+=v;
      }

    for(Int_t ich=0; ich<nCh; ich++)
      {
      charged = (Track*) branchCharged->At(ich);

      TLorentzVector v = charged->P4() ;
      
      EFlowSum+=v;
      }


    if(debug) cout << "Event " << _nEvt << endl << "EFlow object energy sum: " << EFlowSum.E() << " GeV    "
		   << "SIM: " << nParticle << " particles,  EFlow: " 
		   << nCh << " charged, " << nNeu << " neutral, " << nPh << " photons" << endl;
    

    // Final reconstruction objects

    TLorentzVector ERecSum;
    ERecSum.SetPxPyPzE(0.,0.,0.,0.);

    // Jets
    
    int nJet = branchJet->GetEntriesFast();

    for(Int_t ijet=0; ijet<nJet; ijet++)
      {
      jet = (Jet*) branchJet->At(ijet);

      TLorentzVector v = jet->P4() ;
      
      ERecSum+=v;

      if(debug) cout << " -> " << v.E() << " GeV jet " << endl ;

      TRefArray JetPart =  jet->Particles;     // references to generated particles
      TRefArray JetConst = jet->Constituents;  // references to constituents

      int nPart  = JetPart.GetEntries();
      int nConst = JetConst.GetEntries();

      // Generator level particles in a jet
      
      if(debug) cout << "    " << nPart << " generated particles inside jet " << endl;

      double Egen = 0;
      
      for(Int_t ip=0; ip<nPart; ip++)
	{
	  particle = (GenParticle *) JetPart.At(ip);
	  if(particle == NULL) cerr << "GenParticle reference missing ?! for event " << entry << " jet " << ijet << endl;

	  Egen += particle->E;
	  
	  if(debug) cout << "       PID: " << particle->PID << "  E = " << particle->E << " GeV " << endl;
	}

      if(debug) cout << "          Total energy: " << Egen << " GeV " << endl;
      

      // Reconstructed objects inside jet
      
      if(debug) cout << "    " << nConst << " reconstructed objects inside jet " << endl;

      double Erec = 0;
      
      for(Int_t ic=0; ic<nConst; ic++)
	{
	  // Need to check the type
	  
	  const char *cname = JetConst.At(ic)->ClassName();

	  int unknown = 1;
	  
	  if( strcmp(cname,"Tower") == 0)
	    {
            unknown = 0;
	    Tower *ctow = (Tower *) JetConst.At(ic);

	    Erec += ctow->E;
	    
  	    if(debug) cout <<  "       CAL tower with E = " << ctow->E << endl;
	    }

	  
	  if( strcmp(cname,"Track") == 0)
	    {
            unknown = 0;
	    Track *ctrk = (Track *) JetConst.At(ic);

	    // Energy is not stored for track (only Pt and Eta), need to use 4-vector

	    TLorentzVector ctp4 = ctrk->P4();

	    Erec += ctp4.E();
	    
  	    if(debug) cout <<  "       "<< ctrk->PID << " track with E = " << ctp4.E() << endl;
	    }

	  if(unknown)
	    cerr << "Constituent type unknown: " << cname << " for event " << entry << " jet " << ijet << endl;
	  
	}
      
      if(debug) cout << "          Total energy: " << Erec << " GeV " << endl;
      }
    
    // Muons

    int nMuon = branchMuon->GetEntriesFast();

    for(Int_t imuon=0; imuon<nMuon; imuon++)
      {
      muon = (Muon*) branchMuon->At(imuon);

      TLorentzVector v = muon->P4() ;
      
      if(debug) cout << " -> " << v.E() << " GeV isolated muon" << endl ;
      
      ERecSum+=v;
      }
    
    // Electrons

    int nElectron = branchElectron->GetEntriesFast();

    for(Int_t ielectron=0; ielectron<nElectron; ielectron++)
      {
      electron = (Electron*) branchElectron->At(ielectron);

      TLorentzVector v = electron->P4() ;
      
      if(debug) cout << " -> " << v.E() << " GeV isolated electron " << endl ;
      
      ERecSum+=v;
      }

    // Photons

    int nPhoton = branchPhoton->GetEntriesFast();

    for(Int_t iphoton=0; iphoton<nPhoton; iphoton++)
      {
      photon = (Photon*) branchPhoton->At(iphoton);

      TLorentzVector v = photon->P4() ;
      
      if(debug) cout << " -> " << v.E() << " GeV isolated photon " << endl;
      
      ERecSum+=v;
      }

    // What is left

    if(debug) cout << "   ==>>  EFlow - Rec = "
		   << EFlowSum.E()-ERecSum.E() << " GeV" << endl ;
    
    heflow->Fill(EFlowSum.E());
    hetot->Fill(ERecSum.E());
    hedif->Fill(EFlowSum.E()-ERecSum.E());
        
    _nEvt++;
  }
  
  // Open canvas

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

}
