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


void recoil_recon(int Bbit = 2, const char *genFile="bbll_sig_80_eRpL.root",
	       bool isBg = false, const char *jetName="Jet_N2",
	       int nbin=300, double emin=-30., double emax=270.)
{
  // (re)book control histrogram


  if( gDirectory->FindObject("hbdif")!=NULL) gDirectory->Delete("hbdif");
  TH1F *hbdif = new TH1F("hbdif","Energy difference between bquark and measured jet",nbin, -120, 120);

  if( gDirectory->FindObject("hrecondif")!=NULL) gDirectory->Delete("hrecondif");
  TH1F *hrecondif = new TH1F("hrecondif","Energy difference between bquark and reconstructed jet",nbin, -120, 120);

  if( gDirectory->FindObject("hreconmh")!=NULL) gDirectory->Delete("hreconmh");
  TH1F *hreconmh = new TH1F("hreconmh","Higgs mass from reconstructed jet", nbin, 0, 300);

  if( gDirectory->FindObject("hmeasurmh")!=NULL) gDirectory->Delete("hmeasurmh");
  TH1F *hmeasurmh = new TH1F("hmeasurmh","Higgs mass from measured jet", nbin, 0, 300);

  if( gDirectory->FindObject("hrecoilmh")!=NULL) gDirectory->Delete("hrecoilmh");
  TH1F *hrecoilmh = new TH1F("hrecoilmh","Higgs mass from recoil mass", nbin, 0, 300);

  if( gDirectory->FindObject("hfrommz")!=NULL) gDirectory->Delete("hfrommz");
  TH1F *hfrommz = new TH1F("hfrommz","Higgs mass from #sqrt{s} and M_Z", nbin, 0, 300);

  if( gDirectory->FindObject("hllminv")!=NULL) gDirectory->Delete("hllminv");
  TH1F *hllminv = new TH1F("hllminv","Invariant mass of the lepton pair", nbin, 0, 300);

  if( gDirectory->FindObject("hllcostheta")!=NULL) gDirectory->Delete("hllcostheta");
  TH1F *hllcostheta = new TH1F("hllcostheta","#cos(#theta) of leptons", nbin/2, -1, 1);

  if( gDirectory->FindObject("hjjcostheta")!=NULL) gDirectory->Delete("hjjcostheta");
  TH1F *hjjcostheta = new TH1F("hjjcostheta","#cos(#theta) of jets", nbin/2, -1, 1);

  if( gDirectory->FindObject("hjetpt")!=NULL) gDirectory->Delete("hjetpt");
  TH1F *hjetpt = new TH1F("hjetpt","Jet pt", nbin, 0, 150);
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

  TClonesArray *branchMuon = genReader->UseBranch("Muon");
  Muon *muon1;
  Muon *muon2;
  
  TClonesArray *branchElectron = genReader->UseBranch("Electron");
  Electron *electron1;
  Electron *electron2;

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

    //muons
    int nMuons = branchMuon->GetEntriesFast();
    int nElectrons = branchElectron->GetEntriesFast();

    if(!(nMuons == 2 || nElectrons == 2)) continue;

    Double_t pllx, plly;

    TLorentzVector lep1P4, lep2P4;

    if(nElectrons != 2){
      muon1 = (Muon *)branchMuon->At(0);
      muon2 = (Muon *)branchMuon->At(1);
      if(muon1->P4().E()+muon2->P4().E() < 100)continue;
      lep1P4 = muon1->P4();
      lep2P4 = muon2->P4();
      pllx = muon1->P4().Px() + muon2->P4().Px();
      plly = muon1->P4().Py() + muon2->P4().Py();
    }
    else if(nMuons != 2){
      electron1 = (Electron *)branchElectron->At(0);
      electron2 = (Electron *)branchElectron->At(1);
      lep1P4 = electron1->P4();
      lep2P4 = electron2->P4();
      pllx = electron1->P4().Px() + electron2->P4().Px();
      plly = electron1->P4().Py() + electron2->P4().Py();
    }
    else{
      muon1 = (Muon *)branchMuon->At(0);
      muon2 = (Muon *)branchMuon->At(1);
      electron1 = (Electron *)branchElectron->At(0);
      electron2 = (Electron *)branchElectron->At(1);
      if((muon1->P4().E() + muon2->P4().E()) > (electron1->P4().E() + electron2->P4().E())){
        if(muon1->P4().E()+muon2->P4().E() < 100)continue;
        pllx = muon1->P4().Px() + muon2->P4().Px();
        plly = muon1->P4().Py() + muon2->P4().Py();
        lep1P4 = muon1->P4();
        lep2P4 = muon2->P4();
      }
      else{
        pllx = electron1->P4().Px() + electron2->P4().Px();
        plly = electron1->P4().Py() + electron2->P4().Py();
        lep1P4 = electron1->P4();
        lep2P4 = electron2->P4();
      }
    }

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

    jet1p.SetXYZ(v1.Px(), v1.Py(), v1.Pz());
    jet2p.SetXYZ(v2.Px(), v2.Py(), v2.Pz());

    Double_t sqrt_s, alpha, px, py, pt, phi1, phi2, phi12, phi, theta1, theta2, p1, p2;
    sqrt_s = 250;
    alpha = 0;

    px = -pllx + sqrt_s * TMath::Sin(alpha/2);
    py = -plly;

    pt = TMath::Sqrt(px*px + py*py);

    phi = TMath::ATan2(py,px);

    phi1 = jet1p.Phi();
    theta1 = jet1p.Theta();

    phi2 = jet2p.Phi();
    theta2 = jet2p.Theta();

    phi12 = phi1 - phi2;

    p1 = pt*TMath::Sin(phi - phi2)/(TMath::Sin(phi12)*TMath::Sin(theta1));
    p2 = pt*TMath::Sin(phi1 - phi)/(TMath::Sin(phi12)*TMath::Sin(theta2));

    Double_t E1, E2;
    Double_t m1 = v1.M();
    Double_t m2 = v2.M();

    E1 = TMath::Sqrt(p1*p1 + m1*m1);
    E2 = TMath::Sqrt(p2*p2 + m2*m2);

    Double_t reconmH = TMath::Sqrt((1+p2/p1)*m1*m1+(1+p1/p2)*m2*m2+2*pt*pt*TMath::Sin(phi-phi2)*TMath::Sin(phi1-phi)*(1-TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Cos(phi12)-TMath::Cos(theta1)*TMath::Cos(theta2))/(TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Sin(phi12)*TMath::Sin(phi12)));
    //Double_t measuredmH = TMath::Sqrt((1+jet2p.Mag()/jet1p.Mag())*m1*m1+(1+jet1p.Mag()/jet2p.Mag())*m2*m2+2*pt*pt*TMath::Sin(phi-phi2)*TMath::Sin(phi1-phi)*(1-TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Cos(phi12)-TMath::Cos(theta1)*TMath::Cos(theta2))/(TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Sin(phi12)*TMath::Sin(phi12)));
    Double_t measuredmH = (jet1->P4()+jet2->P4()).M();
    Double_t jet1pmag = jet1p.Mag();
    Double_t jet2pmag = jet2p.Mag();
    TLorentzVector jet1corrP4(jet1p.X()*p1/jet1pmag, jet1p.Y()*p1/jet1pmag, jet1p.Z()*p1/jet1pmag, E1);
    TLorentzVector jet2corrP4(jet2p.X()*p2/jet2pmag, jet2p.Y()*p2/jet2pmag, jet2p.Z()*p2/jet2pmag, E2);
    TLorentzVector scalarcorrP4 = jet1corrP4 + jet2corrP4;
    TLorentzVector ZP4 = lep1P4+lep2P4;
    TLorentzVector reconHP4( -ZP4.Px(), -ZP4.Py(), -ZP4.Pz(), 250 - ZP4.E());
    Double_t recoilmass = reconHP4.M();

    TVector3* ZCM_boostvec = new TVector3(-ZP4.Px()/ZP4.E(), -ZP4.Py()/ZP4.E(), -ZP4.Pz()/ZP4.E());
    TVector3* SCM_boostvec = new TVector3(-scalarcorrP4.Px()/scalarcorrP4.E(), -scalarcorrP4.Py()/scalarcorrP4.E(), -scalarcorrP4.Pz()/scalarcorrP4.E());
    TLorentzVector lep1P4_ZCM = lep1P4;
    TLorentzVector lep2P4_ZCM = lep2P4;
    lep1P4_ZCM.Boost(*ZCM_boostvec);
    lep2P4_ZCM.Boost(*ZCM_boostvec);

    TLorentzVector jet1corrP4_SCM= jet1corrP4;
    TLorentzVector jet2corrP4_SCM= jet2corrP4;
    jet1corrP4_SCM.Boost(*SCM_boostvec);
    jet2corrP4_SCM.Boost(*SCM_boostvec);

    hllcostheta->Fill(lep1P4_ZCM.CosTheta());
    hllcostheta->Fill(lep2P4_ZCM.CosTheta());
    hjjcostheta->Fill(jet1corrP4_SCM.CosTheta());
    hjjcostheta->Fill(jet2corrP4_SCM.CosTheta());
    if(reconmH > 300){
      std::cout << jet1->P4().P() << " " << jet2->P4().P() << " " << p1 << " " << p2 << " " << reconmH << " phi: " << phi <<" phi12: " << phi12 << "\n";
      n_phi_pi++;
    }
    
    hrecondif->Fill(bquark->E - E1);
    hrecondif->Fill(antibquark->E - E2);

    hbdif->Fill(bquark->E - v1.E());
    hbdif->Fill(antibquark->E - v2.E());

    hreconmh->Fill(reconmH);
    hmeasurmh->Fill(measuredmH);
    hrecoilmh->Fill(recoilmass);
    hfrommz->Fill(TMath::Sqrt(sqrt_s*sqrt_s - 2*ZP4.E()*sqrt_s + 91.19*91.19));
    //std::cout << ZP4.E() << "\n";

    hllminv->Fill(ZP4.M());

    hjetpt->Fill(jet1corrP4.Pt());
    hjetpt->Fill(jet2corrP4.Pt());
    _nEvt++;
  }

  string path = "bbll_80_eRpL";

  if(isBg) path = "bg_qqll_eRpL";

  //const char* cpath = path.c_str();

  std::cout << "Number of events where mH > 300 GeV: " << n_phi_pi << "\n";
  TCanvas *c1 = new TCanvas;
  auto legend = new TLegend(0.1,0.7,0.45,0.9);
  legend->SetTextSize(0.03);
  hrecondif->SetLineColor(2);
  legend->AddEntry(hbdif, "Jet momentum measurement", "l");
  legend->AddEntry(hrecondif, "Jet momentum reconstruction","l");
  hrecondif->Draw();
  hbdif->Draw("SAME");
  legend->Draw();
  
  c1->SaveAs((path + "/" + "jet_energ_diff_msr_vs_rcnstr.png").c_str());

  TCanvas *c2 = new TCanvas;
  auto legend2 = new TLegend(0.5,0.4,0.9,0.7);
  legend2->SetTextSize(0.03);
  hreconmh->SetLineColor(2);
  hrecoilmh->SetLineColor(3);
  hfrommz->SetLineColor(6);
  legend2->AddEntry(hreconmh, "#splitline{Reconstructed Higgs mass }{from reconstructed jet momentum}", "l");
  legend2->AddEntry(hmeasurmh, "#splitline{Reconstructed Higgs mass }{from measured jet momentum}", "l");
  legend2->AddEntry(hrecoilmh, "#splitline{Reconstructed Higgs mass }{from recoil mass}", "l");
  legend2->AddEntry(hfrommz, "#splitline{Reconstructed Higgs mass }{using M_Z}", "l");
  hrecoilmh->Draw();
  hreconmh->Draw("SAME");
  hmeasurmh->Draw("SAME");
  hfrommz->Draw("SAME");
  legend2->Draw();
  c2->SaveAs((path + "/" + "reconstructed_higgs_mass.png").c_str());

  TCanvas *c3 = new TCanvas;
  hllminv->Draw();
  c3->SaveAs((path + "/" + "leplepinvmass.png").c_str());

  TCanvas *c4 = new TCanvas;
  hllcostheta->Draw();
  c4->SaveAs((path + "/" + "lepcostheta.png").c_str());

  TCanvas *c5 = new TCanvas;
  hjjcostheta->Draw();
  c5->SaveAs((path + "/" + "jetcostheta.png").c_str());

  TCanvas *c6 = new TCanvas;
  hjetpt->Draw();
  c6->SaveAs((path + "/" + "jetpt.png").c_str());


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
