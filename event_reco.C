/*


root -l 'event_reco.C++("bbll_sig_80_eRpL.root","reco/bbll_sig_80_eRpL_reco.root","reco/bbll_sig_80_eRpL_reco",0,-1,80.,359.465771,0.032352,15,7,7,7,7,100,200.)'


*/

#include "CLICdpStyle.C"

#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"

#include <iostream>
#include <fstream>    
#include <iomanip>   

#include "vector"

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


void event_reco(const char *genFile="export/analysis.root",
		const char *outFile="reco/event_reco.root",
		string plname="reco/event_reco",
		int Iproc=0, int Ipol=-1, double Ms=95., double Cs=359.465771, double w=0.032352,
	       int imask = 15, int emask = 7, int mmask = 7 , int tmask=7, int amask = 7,
	       int nbin=100, double mmax=200.)
{
  // Open output file

  TFile *fout = new TFile(outFile,"CREATE");

  if( ! fout->IsOpen() )
    {
     cout << "Error creating output root file !?" << endl;
     exit(2);
    }

  // Define tree

  TTree *reco = new TTree("reco","Reconstructed events");

  // Tree header

   struct event_head {
      Int_t Ievt;
      Int_t Iproc;
      Int_t Ipol;
      Float_t Ms;
      Float_t Cs;
      Float_t w;
      };

   event_head myheader;

   reco->Branch("header",&myheader,"Ievt/I:Iproc/I:Ipol/I:Ms/F:Cs/F:w/F");

   struct event_reco {
      Int_t Nel;
      Int_t Nmu;
      Int_t Nta;
      Int_t Nisr;
      Int_t Nph;
     Float_t Mjj;
     Float_t Mtt;
     Float_t Mcorr;
     Float_t Mrec;
     Float_t Eisr;
     Float_t Eph;
     Float_t En1;
     Float_t En2;
     Float_t Eraw;
     Float_t Etot;
     Float_t y23;
     Float_t y34;
     Float_t y45;
     Int_t  QT[2];
     Int_t  BJ[2];
      };

   event_reco myevent;

   reco->Branch("event",&myevent,"Nel/I:Nmu/I:Nta/I:Nisr/I:Nph/I:Mjj/F:Mtt/F:Mcorr/F:Mrec/F:Eisr/F:Eph/F:En1/F:En2/F:Eraw/F:Etot/F:y23/F:y34/F:y45/F:Q1/I:Q2/I:B1/I:B2/I" );

    // Four-vectors
   
   TLorentzVector h1p4 ;
   reco->Branch("H","TLorentzVector",&h1p4);

   TLorentzVector zp4 ;
   reco->Branch("Z","TLorentzVector",&zp4);

   TLorentzVector hzp4 ;
   reco->Branch("HZ","TLorentzVector",&hzp4);
   
   TLorentzVector h1p4corr2;
   reco->Branch("Hc","TLorentzVector",&h1p4corr2);

   TLorentzVector hzp4corr2;
   reco->Branch("HZc","TLorentzVector",&hzp4corr2);
   
   TLorentzVector p4rec0;
   reco->Branch("Rec0","TLorentzVector",&p4rec0);

   TLorentzVector p4rec1;
   reco->Branch("Rec1","TLorentzVector",&p4rec1);

   TLorentzVector p4rec2;
   reco->Branch("Rec2","TLorentzVector",&p4rec2);

   TLorentzVector p4isr;
   reco->Branch("ISR","TLorentzVector",&p4isr);

   TLorentzVector p4ph;
   reco->Branch("Ph","TLorentzVector",&p4ph);

   TLorentzVector myrec[4];
   reco->Branch("T1","TLorentzVector",myrec);
   reco->Branch("T2","TLorentzVector",myrec+1);
   reco->Branch("J1","TLorentzVector",myrec+2);
   reco->Branch("J2","TLorentzVector",myrec+3);

   TLorentzVector T1star;
   reco->Branch("T1s","TLorentzVector",&T1star);
   
  // Setup graphics output

  CLICdpStyle();
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
 
  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(1);

  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);

  // Load Delphe library

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

  TClonesArray *branchParticle = genReader->UseBranch("Particle");
  GenParticle *particle;
  
  TClonesArray *branchJet4 = genReader->UseBranch("Jet_N4");
  TClonesArray *branchJet3 = genReader->UseBranch("Jet_N3");
  TClonesArray *branchJet2 = genReader->UseBranch("Jet_N2");
  TClonesArray *branchJet;
  Jet *jet;

  TClonesArray *branchElectron = genReader->UseBranch("Electron");
  Electron *electron;
  
  TClonesArray *branchMuon = genReader->UseBranch("Muon");
  Muon *muon;
  
  TClonesArray *branchPhoton = genReader->UseBranch("Photon");
  TClonesArray *branchBCal = genReader->UseBranch("BCalPhoton");
  Photon *photon;




  
  TH2D *hmtt = new TH2D("hmtt","Invariant Mass 2-D", nbin/2, 0., mmax, nbin/2, 0., mmax);
  hmtt->GetXaxis()->SetTitle("m_{#tau#tau} [GeV]");
  hmtt->GetYaxis()->SetTitle("m_{jj} [GeV]");

  TH2D *hmcr1 = new TH2D("hmcr1","Recoil vs reconstructed mass", nbin/2, 0., mmax, nbin/2, 0., mmax);
  hmcr1->GetXaxis()->SetTitle("m_{#tau#tau}^{corr} [GeV]");
  hmcr1->GetYaxis()->SetTitle("m_{recoil} [GeV]");

  TH2D *hmcr2 = new TH2D("hmcr2","Recoil vs reconstructed mass", nbin/2, 0., mmax, nbin/2, 0., mmax);
  hmcr2->GetXaxis()->SetTitle("m_{#tau#tau}^{corr} [GeV]");
  hmcr2->GetYaxis()->SetTitle("m_{recoil} [GeV]");

  TH1D *hmt1 = new TH1D("hmt1","Invariant Mass  - tau pair", nbin, 0., mmax);
  hmt1->GetXaxis()->SetTitle("m_{#tau#tau} [GeV]");
  hmt1->GetYaxis()->SetTitle("Events");

  TH1D *hmt2 = new TH1D("hmt2","Invariant Mass  - jet pair", nbin, 0., mmax);
  hmt2->GetXaxis()->SetTitle("m_{jj} [GeV]");
  hmt2->GetYaxis()->SetTitle("Events");

  TH1D *hmrec = new TH1D("hmrec","Recoil Mass", nbin, 0., mmax);
  hmrec->GetXaxis()->SetTitle("m_{recoil} [GeV]");
  hmrec->GetYaxis()->SetTitle("Events");

  TH1D *hmc1 = new TH1D("hmc1","Corrected tau pair invariant mass", nbin, 0., mmax);
  hmc1->GetXaxis()->SetTitle("m_{#tau#tau}^{corr} [GeV]");
  hmc1->GetYaxis()->SetTitle("Events");

  TH1D *hmc2 = new TH1D("hmc2","Corrected tau pair invariant mass with E_{#nu} cut", nbin, 0., mmax);
  hmc2->GetXaxis()->SetTitle("m_{#tau#tau}^{corr} [GeV]");
  hmc2->GetYaxis()->SetTitle("Events");

  TH1D *henu = new TH1D("henu","Neutrino energy from correction", nbin, -mmax, mmax);
  henu->GetXaxis()->SetTitle("E_{#nu} [GeV]");
  henu->GetYaxis()->SetTitle("Events");

  TH1D *hetot = new TH1D("hetot","Total energy", nbin, 0, 2*mmax);
  hetot->GetXaxis()->SetTitle("E_{tot} [GeV]");
  hetot->GetYaxis()->SetTitle("Events");

  TH1D *hetc1 = new TH1D("hetc1","Total energy with correction", nbin, 0, 2*mmax);
  hetc1->GetXaxis()->SetTitle("E_{tot}^{corr} [GeV]");
  hetc1->GetYaxis()->SetTitle("Events");

  TH1D *hetc2 = new TH1D("hetc2","Total energy with correction and  E_{#nu} cut", nbin, 0, 2*mmax);
  hetc2->GetXaxis()->SetTitle("E_{tot}^{corr} [GeV]");
  hetc2->GetYaxis()->SetTitle("Events");

  // Loop over tree events

  int nEvt=0;

  for(Int_t entry = 0; entry < genEntries; ++entry)
  {
    // Load selected branches with data from specified event
    
    // bool debug =  entry<10 || entry%1000==0;
    
    genReader->ReadEntry(entry);

    int nParticle = branchParticle->GetEntriesFast();

    int nMuon = branchMuon->GetEntriesFast();
    int nElectron = branchElectron->GetEntriesFast();
    int nPhoton = branchPhoton->GetEntriesFast();
    int nBcal = branchBCal->GetEntriesFast();


    // Photons

    int nIsr = 0;
    int nPh = 0;
    
    p4isr.SetPxPyPzE(0.,0.,0.,0.);
    p4ph.SetPxPyPzE(0.,0.,0.,0.);

    // Loop over photons in BCAL
    
    for(int iPh=0; iPh<nBcal; iPh++)
      {
      photon = (Photon*) branchBCal->At(iPh);

      // All BCAL photons are considered ISR photons
      
      nIsr++;
      p4isr+=photon->P4();
      }
    
    // Loop over photons in calorimeter
    
    for(int iPh=0; iPh<nPhoton; iPh++)
      {
      photon = (Photon*) branchPhoton->At(iPh);

      // ISR photon candidate selection
      
      if(abs(photon->Eta) > 2. && photon->E > 10.)
	{
	  nIsr++;
	  p4isr+=photon->P4();
	}
      else
	{
	  nPh++;
	  p4ph+= photon->P4();
	}

      }
    
    // Select only given ISR photon tag configuration

    if ( ( 1<<nIsr & imask)  == 0) continue;
    
    // Electrons
    
    int nEle = 0;
    TLorentzVector myelectron[4];
    int myqel[4];

    int Qtot = 0;
    
    for(int iEle=0; iEle<nElectron; iEle++)
      {
      electron = (Electron*) branchElectron->At(iEle);

      // Electron candidate selection
      
      if(abs(electron->Eta) > 2. || electron->PT < 5.) continue;

      myelectron[nEle] = electron->P4();
      myqel[nEle]=electron->Charge;
      Qtot += electron->Charge;
      nEle++;

      if(nEle==4)break;
      }
    
    // Select only given electron tag configuration

    if ( ( 1<<nEle & emask)  == 0) continue;

    
    // Muons
    
    int nMu = 0;
    TLorentzVector mymuon[4];
    int myqmu[4];
    
    for(int iMu=0; iMu<nMuon; iMu++)
      {
      muon = (Muon*) branchMuon->At(iMu);

      // Muon candidate selection
      
      if(abs(muon->Eta) > 2. || muon->PT < 5.) continue;

      mymuon[nMu] = muon->P4();
      myqmu[nMu]=muon->Charge;
      Qtot += muon->Charge;
      nMu++;

      if(nMu==4)break;
      }
    
    // Select only given muon tag configuration

    if ( ( 1<<nMu & mmask)  == 0) continue;


    // At most two leptons allowed

    if ( nEle + nMu > 2 ) continue;

    // Oposit sign required for two leptons (not relevant in Delphes)
    
    if ( nEle + nMu == 2 && Qtot != 0) continue;
    
    // Jets
    
    int nJet = 4 - nMu - nEle;

    if (nJet < 2) continue;

    if (nJet == 2) branchJet = branchJet2;
    if (nJet == 3) branchJet = branchJet3;
    if (nJet == 4) branchJet = branchJet4;
    
    if( branchJet->GetEntriesFast() != nJet) continue;
    
    TLorentzVector myjet[4];
    TLorentzVector jtmp;
    int mytag[4];
    int mybtag[4];
    int nTau=0;
    int itmp;
    
    for(int ijet=0;ijet<nJet;ijet++)
      {
      Jet* j1= (Jet*) branchJet->At(ijet);

       myjet[ijet] = j1->P4();
       mytag[ijet] = j1->TauTag;
       mybtag[ijet] = j1->BTag;
       nTau+=j1->TauTag;

       if(ijet==0)
	 {
	   myevent.y23 = j1->ExclYmerge23; 
	   myevent.y34 = j1->ExclYmerge34;
	   myevent.y45 = j1->ExclYmerge45;
	 }
      }

    // Select only given tag configuration

    if ( ( 1<<nTau & tmask)  == 0) continue;

    
    // Check combined tag number

    int ntag = nEle + nMu + nTau;

    if ( ( 1<<ntag & amask)  == 0) continue;

    
    // Count selected events only
    
    nEvt++;
    bool debug = ( nEvt<=5 || nEvt%1000==0 ) && (plname.length() == 0); 
    
    
    // Sort jets. Those with tau tag go first. For same tag value, lower masses go first
    
    for(int ijet=1;ijet<nJet;)
      {
	if(mytag[ijet] > mytag[ijet-1] || (mytag[ijet] == mytag[ijet-1] && myjet[ijet].M() < myjet[ijet-1].M()) )
	 {
	  jtmp = myjet[ijet];
	  myjet[ijet]=myjet[ijet-1];
	  myjet[ijet-1]=jtmp;

	  itmp = mytag[ijet];
	  mytag[ijet]=mytag[ijet-1];
	  mytag[ijet-1]=itmp;

	  itmp = mybtag[ijet];
	  mybtag[ijet]=mybtag[ijet-1];
	  mybtag[ijet-1]=itmp;

	  ijet = 1;
	 }
       else
	 ijet++;

       }

    if(debug)
      {
	cout << "Event " << nEvt << " (" << entry << ")" << endl;
      
	cout << nIsr << " ISR photons with E = " << p4isr.E() << endl;
      
	cout << nPh  << " central/soft photons with E = " << p4ph.E() << endl;
      
      for(int iEle=0;iEle<nEle;iEle++)
	cout << "Electron " << iEle << " PT = " << myelectron[iEle].Pt() << " Eta = " << myelectron[iEle].Eta() << endl;
      
      for(int iMu=0;iMu<nMu;iMu++)
	cout << "Muon " << iMu << " PT = " << mymuon[iMu].Pt() << " Eta = " << mymuon[iMu].Eta() << endl;
      
      for(int ijet=0;ijet<nJet;ijet++)
	  cout << "Jet " << ijet << " Tag = " << mytag[ijet] << " M = " << myjet[ijet].M() << endl;
      }


     // Combine leptons and jets
    
     // TLorentzVector myrec[4];

     myevent.QT[0] =  myevent.QT[1] = 0;
       
     for(int iEle=0;iEle<nEle;iEle++)
       {
       myrec[iEle]=myelectron[iEle];
       myevent.QT[iEle]=myqel[iEle];
       }
     
     for(int iMu=0;iMu<nMu;iMu++)
       {
       myrec[nEle+iMu] = mymuon[iMu];
       myevent.QT[nEle+iMu]=myqmu[iMu];
       }
       
     for(int ijet=0;ijet<nJet;ijet++)
       myrec[nEle+nMu+ijet]=myjet[ijet];
     

     myevent.BJ[0] = mybtag[nJet-2];
     myevent.BJ[1] = mybtag[nJet-1];
     
     // First two reconstructed objects should be tau decay products...
     
     
     h1p4 = myrec[0]+myrec[1];
     zp4 = myrec[2]+myrec[3];

     hzp4 = h1p4 + zp4;

     // Include measured photons in total sum

     hzp4 += p4isr + p4ph;
    
   	
     if(debug) cout << " m_tt = " << h1p4.M() << " GeV" 
		      << ", m_jj = " << zp4.M() << " GeV" 
                      << endl;

    
     hmt1->Fill(h1p4.M());
     hmt2->Fill(zp4.M());
     hmtt->Fill(h1p4.M(),zp4.M());
     hetot->Fill(hzp4.E());

     
     // Correcting for neutrinos from tau decays

     TVector3 nt1 = myrec[0].Vect().Unit();
     TVector3 nt2 = myrec[1].Vect().Unit();
 
     double wg = nt1.X()*nt2.Y() - nt2.X()*nt1.Y();
     
     double En1 = (nt2.X()*hzp4.Py() - nt2.Y()*hzp4.Px())/wg;
     double En2 = (nt1.Y()*hzp4.Px() - nt1.X()*hzp4.Py())/wg;

     TVector3 Pn1 = En1 * nt1;
     TVector3 Pn2 = En2 * nt2;

     TLorentzVector n1p4(Pn1,En1);
     TLorentzVector n2p4(Pn2,En2);

     if(En1>0 && En1<250) myrec[0] += n1p4;
     if(En2>0 && En2<250) myrec[1] += n2p4;
     
     TLorentzVector hzp4corr = hzp4 + n1p4 + n2p4;

     TLorentzVector h1p4corr = h1p4 + n1p4 + n2p4;

     h1p4corr2 = myrec[0] + myrec[1];

     hzp4corr2 = h1p4corr2 + zp4;

     // Recoil mass

     p4rec0.SetPxPyPzE(0.,0.,0.,250.);

     // Subtract reconstructed Z
     
     p4rec0 -= zp4;

     // Subtract ISR photons

     p4rec1 = p4rec0 - p4isr;

     // Subtract also soft photons;

     p4rec2 = p4rec1 - p4ph;
     
     if(debug){
       cout << "Estimated corrections: En1 = " << En1 << "     En2 = " << En2 << endl;
       cout << "   tau jets corrected invariant mass = " << h1p4corr2.M() << " GeV"  << endl;
       cout << "   recoil mass = " << p4rec1.M() << " GeV"  << endl;
     }

     henu->Fill(En1);
     henu->Fill(En2);
     
     hmc1->Fill(h1p4corr.M());
     hmc2->Fill(h1p4corr2.M());

     hetc1->Fill(hzp4corr.E());
     hetc2->Fill(hzp4corr2.E());

     hmrec->Fill(p4rec1.M());
     
     hmcr1->Fill(h1p4corr.M(),p4rec1.M());
     hmcr2->Fill(h1p4corr2.M(),p4rec1.M());

     // Fill output tree

     myheader.Ievt  = nEvt;
     myheader.Iproc = Iproc;
     myheader.Ipol = Ipol;
     myheader.Ms = Ms;
     myheader.Cs = Cs;
     myheader.w = w;

     
     myevent.Nel = nEle;
     myevent.Nmu = nMu;
     myevent.Nta = nTau;
     myevent.Nph = nPh;
     myevent.Nisr = nIsr;
     
     myevent.Mjj = zp4.M();
     myevent.Mtt = h1p4.M();
     myevent.Mcorr = h1p4corr2.M();
     myevent.Mrec = p4rec2.M();
     myevent.Eisr = p4isr.E();
     myevent.Eph = p4ph.E();
     myevent.En1 = En1;
     myevent.En2 = En2;
     myevent.Eraw = hzp4.E();
     myevent.Etot = hzp4corr2.E();

     // Calculate tau parameters in H frame
     
     TVector3 Hboost;
     Hboost = h1p4corr2.BoostVector();
     
    // Reference frame defined by pair momentum direction (new Z direction)
	
    TVector3 zlab(0., 0., 1.);

    TVector3 zpair(1.,0.,0.);
    if(Hboost.Mag()>0) zpair=Hboost.Unit();
	
    // Vector normal to the scattering plane - Y axis

    TVector3 ypair = zlab.Cross(zpair);
    if(ypair.Mag()==0.)
	  ypair.SetXYZ(0.,1.,0.);
    else
	  ypair.SetMag(1.);

    // X axis from the vector product

    TVector3 xpair = ypair.Cross(zpair);

    // Tau momentum in the Higgs reference frame (Z along Higgs direction)
    
    T1star.SetPxPyPzE(xpair*myrec[0].Vect(),ypair*myrec[0].Vect(),zpair*myrec[0].Vect(),myrec[0].E());

     // Boost to the pair rest frame

     TVector3 bpair(0.,0.,-Hboost.Mag());
	
     T1star.Boost(bpair);

     reco->Fill();
     }

  cout << nEvt << " events accepted out of " << genEntries
       << " (" << 100.*nEvt/genEntries << ")" <<  endl;

  cout << "Input chain contains " << w*genEntries << " weighted events" << endl;
  cout << w*nEvt << " weighted events accepted out of " << genEntries << endl;

   // Open canvas

  TCanvas  *ch0 = (TCanvas *) gROOT->FindObject("ch0");
  if(!ch0)  ch0 = new TCanvas("ch0","ch0: 2D raw mass plot",20,20,800,600);
  ch0->cd();

  hmtt->Draw("COL");

  string pdfname, pngname;
  
  if(plname.length() > 0)
    {
    pdfname = plname +"_rawmass_2d.pdf";
    pngname = plname +"_rawmass_2d.png";
    ch0->Print(pdfname.c_str());
    ch0->Print(pngname.c_str());
    }
  
  TCanvas  *ch1 = (TCanvas *) gROOT->FindObject("ch1");
  if(!ch1)  ch1 = new TCanvas("ch1","ch1: 2D mass plot",550,20,800,600);
  ch1->cd();

  hmcr2->Draw("COL");

  if(plname.length() > 0)
    {
    pdfname = plname +"_mass_2d.pdf";
    pngname = plname +"_mass_2d.png";
    ch1->Print(pdfname.c_str());
    ch1->Print(pngname.c_str());
    }
  
  // gStyle->SetOptLogy(1);

  TCanvas  *ch2 = (TCanvas *) gROOT->FindObject("ch2");
  if(!ch2)  ch2 = new TCanvas("ch2","ch2: reconstructed scalar mass",1100,20,800,600);
  ch2->cd();

  hmc2->SetLineWidth(2);
  hmc2->SetLineStyle(1);
  hmc2->SetLineColor(4);
  hmc2->Draw();

  hmt2->SetLineWidth(2);
  hmt2->SetLineStyle(7);
  hmt2->SetLineColor(4);
  // hmt2->Draw("same");
  
  hmt1->SetLineWidth(2);
  hmt1->SetLineStyle(1);
  hmt1->SetLineColor(2);
  hmt1->Draw("same");

  hmrec->SetLineWidth(2);
  hmrec->SetLineStyle(1);
  hmrec->SetLineColor(6);
  //  hmrec->Draw("same");

  hmc1->SetLineWidth(2);
  hmc1->SetLineStyle(1);
  hmc1->SetLineColor(3);
  //  hmc1->Draw("same");
  // hmc1->Fit("gaus","","same",mfit-wfit,mfit+wfit);

  if(plname.length() > 0)
    {
    pdfname = plname +"_mass.pdf";
    pngname = plname +"_mass.png";
    ch2->Print(pdfname.c_str());
    ch2->Print(pngname.c_str());
    }
  


  
  // gStyle->SetOptLogy(1);

  TCanvas  *ch3 = (TCanvas *) gROOT->FindObject("ch3");
  if(!ch3)  ch3 = new TCanvas("ch3","ch3: nutrino energy from fit",20,320,800,600);
  ch3->cd();

  henu->SetLineWidth(2);
  henu->SetLineStyle(1);
  henu->SetLineColor(4);
  henu->Draw();
  
  
  if(plname.length() > 0)
    {
    pdfname = plname +"_enu.pdf";
    pngname = plname +"_enu.png";
    ch3->Print(pdfname.c_str());
    ch3->Print(pngname.c_str());
    }


  
  TCanvas  *ch4 = (TCanvas *) gROOT->FindObject("ch4");
  if(!ch4)  ch4 = new TCanvas("ch4","ch4: total measured energy",620,320,800,600);
  ch4->cd();

  hetc1->SetLineWidth(2);
  hetc1->SetLineStyle(1);
  hetc1->SetLineColor(3);
  // hetc1->Draw();
  
  hetc2->SetLineWidth(2);
  hetc2->SetLineStyle(1);
  hetc2->SetLineColor(4);
  hetc2->Draw();
  
  hetot->SetLineWidth(2);
  hetot->SetLineStyle(1);
  hetot->SetLineColor(2);
  hetot->Draw("same");


  if(plname.length() > 0)
    {
    pdfname = plname +"_etot.pdf";
    pngname = plname +"_etot.png";
    ch4->Print(pdfname.c_str());
    ch4->Print(pngname.c_str());
    }
  
  
  TCanvas  *ch5 = (TCanvas *) gROOT->FindObject("ch5");
  if(!ch5)  ch5 = new TCanvas("ch5","ch5: measured Z mass",1120,320,800,600);
  ch5->cd();

  hmt2->SetLineWidth(2);
  hmt2->SetLineStyle(1);
  hmt2->SetLineColor(4);
  hmt2->Draw();

  if(plname.length() > 0)
    {
    pdfname = plname +"_mz.pdf";
    pngname = plname +"_mz.png";
    ch5->Print(pdfname.c_str());
    ch5->Print(pngname.c_str());
    }

  
  TCanvas  *ch6 = (TCanvas *) gROOT->FindObject("ch6");
  if(!ch6)  ch6 = new TCanvas("ch6","ch6: measured recoil mass",1120,620,800,600);
  ch6->cd();

  hmrec->SetLineWidth(2);
  hmrec->SetLineStyle(1);
  hmrec->SetLineColor(4);
  hmrec->Draw();

  if(plname.length() > 0)
    {
    pdfname = plname +"_mrec.pdf";
    pngname = plname +"_mrec.png";
    ch6->Print(pdfname.c_str());
    ch6->Print(pngname.c_str());
    }
  

    // Condenced output
    
    cout <<  "EFF  " << imask << " " << emask << " " <<  mmask  << " " <<  tmask   << " " <<  amask  << "   "
         << nEvt << " " << w*nEvt << " " << 100.*nEvt/genEntries << endl;


    // Close output file

  
  fout->Write();
  
  fout->Close();


}
