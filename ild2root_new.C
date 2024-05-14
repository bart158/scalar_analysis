/*

Extracts information from Delphes tree to text file

root -l 'ild2root_new.C++("data/sm_lep90_pythia6_analysis.root","export/analysis.root")'


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


void ild2root_new(const char *genFile="bbll_sig_80_eRpL.root",const char *outFile="export/analysis.root")
{
  // Open canvas

  TCanvas  *ch1 = (TCanvas *) gROOT->FindObject("ch1");
  if(!ch1)  ch1 = new TCanvas("ch1","Delphes histograms",20,20,600,450);
  ch1->cd();
   
  // (re)book control histrogram

  int nbin=480;
  double mmin=0.;
  double mmax=120.;

  if( gDirectory->FindObject("hetot")!=NULL) gDirectory->Delete("hetot");
  TH1F *hetot = new TH1F("hetot","Total event energy",nbin, mmin, mmax);

    
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
  

  TClonesArray *branchCharged = genReader->UseBranch("EFlowTrack");
  Track *charged;
  
  TClonesArray *branchNeutral = genReader->UseBranch("EFlowNeutralHadron");
  Tower *neutral;
  
  TClonesArray *branchPhoton = genReader->UseBranch("EFlowPhoton");
  Tower *photon;

  

  TClonesArray *branchJet = genReader->UseBranch("Jet");
  Jet *jet;

  TClonesArray *branchElectron = genReader->UseBranch("Electron");
  Electron *electron;

  TClonesArray *branchMuon = genReader->UseBranch("Muon");
  Muon *muon;



  // New root file
  
  TFile* outputFile;
  TTree* _tree;
  
  outputFile = NULL;
  outputFile = new TFile(outFile, "CREATE" );

  if(!outputFile->IsOpen())
    {
    cout << "Error opening output root file " << outFile << " !?!" << endl;
    return;
    }
   
  // Create TTree now
   
   _tree = new TTree("P4", "Root tree for P4 analysis" );

  // Numbers

   int eventNum[4];
    
  // Global four-momentum
  
  TLorentzVector *EFlowSum = new TLorentzVector();
  EFlowSum->SetPxPyPzE(0,0,0,0);


  // Highest Pt objects
  
  Int_t Nmax = 3;
  int P4type[3];

  TLorentzVector *P4obj = new TLorentzVector[Nmax];

  for(Int_t iobj=0;iobj<Nmax;iobj++)
      {
      P4type[iobj]=0;
      P4obj[iobj].SetPxPyPzE(0,0,0,0);
      }


  _tree->Branch("eventNum", &eventNum, "nEvt/I:nElectron:nMuon:nJet");

  _tree->Branch("ObjectID", &P4type, "ID1/I:ID2:ID3");

  _tree->Branch("Total","TLorentzVector",&EFlowSum,16000,0);

  _tree->Branch("P1", "TLorentzVector", &P4obj[0]);  
  _tree->Branch("P2", "TLorentzVector", &P4obj[1]);  
  _tree->Branch("P3", "TLorentzVector", &P4obj[2]);  
          
  
  int _nEvt=0;
  
  // Loop over tree events
  
  for(Int_t entry = 0; entry < genEntries; ++entry)
  {
    bool debug =  entry<10 || entry%100==0;
    
    // Load selected branches with data from specified event
    
    genReader->ReadEntry(entry);

    int nCh = branchCharged->GetEntriesFast();
    int nNeu = branchNeutral->GetEntriesFast();
    int nPh = branchPhoton->GetEntriesFast();

   // Fill energy sum
   
    EFlowSum->SetPxPyPzE(0.,0.,0.,0.);

    for(Int_t ineu=0; ineu<nNeu; ineu++)
      {
      neutral = (Tower*) branchNeutral->At(ineu);

      TLorentzVector v = neutral->P4() ;
      
      *EFlowSum+=v;
      }

    for(Int_t iphot=0; iphot<nPh; iphot++)
      {
      photon = (Tower*) branchPhoton->At(iphot);

      TLorentzVector v = photon->P4() ;
      
      *EFlowSum+=v;
      }

    for(Int_t ich=0; ich<nCh; ich++)
      {
      charged = (Track*) branchCharged->At(ich);

      TLorentzVector v = charged->P4() ;
      
      *EFlowSum+=v;
      }

    hetot->Fill(EFlowSum->E());

    if(debug)
      {
	if(debug) cout << entry << " " << nCh
		       << "+" << nNeu
		       << "+" << nPh
		       << " -> ";

      cout << "E= " << EFlowSum->E() << " GeV  ";
      cout << "Pt= " << EFlowSum->Pt() << " GeV  ";
      cout << "Pz= " << EFlowSum->Pz() << " GeV  ";
      cout << "M= " << EFlowSum->M() << " GeV  ";
      cout << endl;
      }

    // Fill objects

    Int_t Nobj = 0;
    
    int nElectron = branchElectron->GetEntriesFast();
    int nMuon = branchMuon->GetEntriesFast();
    int nJet = branchJet->GetEntriesFast();

    eventNum[0]= _nEvt;
    eventNum[1]= nElectron;
    eventNum[2]= nMuon;
    eventNum[3]= nJet;
    
    for(Int_t iobj=0;iobj<Nmax;iobj++)
      {
	P4type[iobj]=0;
         P4obj[iobj].SetPxPyPzE(0,0,0,0);
      }

    for(Int_t ijet=0; ijet<nJet && ijet<Nmax; ijet++)
      {
	jet = (Jet *) branchJet->At(ijet);
	
	P4type[ijet]=3;
	P4obj[ijet]=jet->P4();
      }


    Int_t iele = 0;
    
    for(Int_t iobj=0; iele<nElectron && iobj<Nmax; iobj++)
      {
	electron = (Electron *) branchElectron->At(iele);

	if(electron->P4().Pt() < P4obj[iobj].Pt()) continue;

	// Shift possible contents
	
	for(Int_t iob2=Nmax-1; iob2>iobj;iob2--)
	  {
	  P4type[iob2]=P4type[iob2-1];
  	  P4obj[iob2] =P4obj[iob2-1];
	  }  
	 
	P4type[iobj]=1;
	P4obj[iobj]=electron->P4();
	iele++;
      }




    Int_t imu = 0;
    
    for(Int_t iobj=0; imu<nMuon && iobj<Nmax; iobj++)
      {
	muon = (Muon *) branchMuon->At(imu);

	if(muon->P4().Pt() < P4obj[iobj].Pt()) continue;

	// Shift possible contents
	
	for(Int_t iob2=Nmax-1; iob2>iobj;iob2--)
	  {
	  P4type[iob2]=P4type[iob2-1];
  	  P4obj[iob2] =P4obj[iob2-1];
	  }  
	 
	P4type[iobj]=2;
	P4obj[iobj]=muon->P4();
	imu++;
      }

    _tree->Fill();

    _nEvt++;
  }
  
 _tree->Write();
 outputFile->Close();

 ch1->Clear();

 hetot->SetLineColor(2);
 hetot->GetXaxis()->SetTitle("E [GeV]");
 hetot->GetYaxis()->SetTitle("Entries");
 hetot->SetMinimum(0.7);
  
 hetot->Draw();

}
