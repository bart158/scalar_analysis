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

#include <RtypesCore.h>
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
/*
# Proc     eLpL        eRpR        eLpR        eRpL     
  qqqq       -           -         69.9234     103.269  
  qqtt       -           -         337.801     632.644  
  qqll     383.667     384.444     97.2588     155.698  
  qqvv       -           -          93.925     254.328  
  qqlv     849.756     849.384     102.891     576.273  
  qqtv       -           -          106.57      1154.4  
  ttll     688.559     688.318     487.072     556.402  
  tttt       -           -         8356.19     13483.2  
    qq       -           -         15.7433     14.4222 
qqllvv       -           -         1110.82     1879.71  
*/


/*
# Mass  bbvv_eLpR   bbvv_eRpL   bbll_eLpR   bbll_eRpL   bbqq_eLpR   bbqq_eRpL  
    20       -           -         520.55       812.6        -           -     
    30     297.75      464.54      590.02      920.67      87.295      136.41  
    40       -           -          652.3      1014.4        -           -     
    50     359.83      561.27      713.05      1113.3      105.72      165.04  
    60       -           -         779.52      1217.8        -           -     
    70       -           -         856.09      1335.9        -           -     
    80      479.5       747.5       951.7      1484.5       140.8      220.04  
    90       -           -         1071.6      1670.1        -           -     
    95     576.05      900.74      1143.2      1781.2      169.27      264.54  
   100       -           -         1225.9      1915.7        -           -     
   110      726.4      1132.7      1440.8      2247.2      213.29      333.71  
   120       -           -         1753.8      2737.8        -           -     
*/
// eLpR = 0; eRpL = 1; eLpL = 2; eRpR = 3
// qqqq = 0; qqtt = 1; qqll = 2; qqvv = 3; qqlv = 4; qqtv = 5; ttll = 6; tttt = 7; qq = 8; qqllvv = 9
Double_t Lgen_bg_arr[4][10] = {
    {69.9234, 337.801, 97.2588, 93.925, 102.891, 106.57, 487.072, 8356.19, 15.7433, 1110.82},
    {103.269, 632.644, 155.698, 254.328, 576.273, 1154.4, 556.402, 13483.2, 14.4222, 1879.71},
    {0, 0, 383.667, 0, 849.756, 0, 688.559, 0, 0, 0},
    {0, 0, 384.444, 0, 849.384, 0, 688.318, 0, 0, 0}
};
// eLpR = 0; eRpL = 1; eLpL = 2; eRpR = 3
// 30 = 0; 50 = 1, 80 = 2, 95 = 3, 110 = 4, 20 = 5, 40 = 6, 60 = 7, 70 = 8, 90 = 9, 100 = 10, 120 = 11
Double_t Lgen_sig_bbll_arr[2][12] = {
    {590.02, 713.05, 951.7, 1143.2, 1440.8, 520.55, 652.3, 779.52, 856.09, 1071.6, 1225.9, 1753.8},
    {920.67, 1113.3, 1484.5, 1781.2, 2247.2, 812.6, 1014.4, 1217.8, 1335.9, 1670.1, 1915.7, 2737.8}
};

void make_new_ttree(const char *genFile="SM_bg_eRpL/qqllvv_bg_eRpL.root",
		const char *outFile="trees_for_training/qqllvv_bg_eRpL_new.root",
		string plname="trees_for_training/qqllvv_bg_eRpL_new",
		int Iproc=9, int Ipol=1, int Ms=0, double Cs=359.465771, double w=0.032352,
	       int imask = 15, int emask = 7, int mmask = 7 , int tmask=7, int amask = 7,
	       int bmask = 7, int nbin=100, double mmax=200., int Bbit = 2){



    TFile *fout = new TFile(outFile, "RECREATE");

    if( ! fout->IsOpen() )
    {
        cout << "Error creating output root file !?" << endl;
        exit(2);
    }

    TTree *newtree = new TTree("tree", "New tree for TMVA");

    struct event_head {
        Int_t Ievt;
        Int_t Iproc;
        Int_t Ipol;
        Float_t Ms;
        Float_t Cs;
        Float_t Lgen; // Lgen
    };
    
    event_head myheader;

    newtree->Branch("header",&myheader,"Ievt/I:Iproc/I:Ipol/I:Ms/F:Cs/F:Lgen/F");

    struct event_reco{
        Int_t Nel;
        Int_t Nmu;
        Int_t Nisr;
        Int_t Nph;
        Int_t jet1btag;
        Int_t jet2btag;
        Float_t Mjj;
        Float_t Mcorr;
        Float_t Mrec;
        Float_t Etot;
        Float_t y23;
        Float_t y34;
        Float_t y45;
        Float_t costhetajetcms;
        Float_t jetpt;
    };

    event_reco myevent;
    newtree->Branch("event",&myevent,"Nel/I:Nmu/I:Nisr/I:Nph/I:jet1btag/I:jet2btag/I:Mjj/F:Mcorr/F:Mrec/F:Etot/F:y23/F:y34/F:y45/F:costhetajetcms/F:jetpt/F" );
    TLorentzVector h1p4;
    newtree->Branch("H", "TLorentzVector", &h1p4);
    TLorentzVector zp4;
    newtree->Branch("Z", "TLorentzVector", &zp4);
    TLorentzVector hzp4;
    newtree->Branch("HZ", "TLorentzVector", &hzp4);
    TLorentzVector myrec[4];
    newtree->Branch("LEP1","TLorentzVector",myrec);
    newtree->Branch("LEP2","TLorentzVector",myrec+1);
    newtree->Branch("J1","TLorentzVector",myrec+2);
    newtree->Branch("J2","TLorentzVector",myrec+3);
    TLorentzVector p4isr;
    newtree->Branch("ISR","TLorentzVector",&p4isr);

    TLorentzVector p4ph;
    newtree->Branch("Ph","TLorentzVector",&p4ph);

    TLorentzVector recoilhp4;
    newtree->Branch("recoilHP4", "TLorentzVector", &recoilhp4);
    TLorentzVector hp4corr;
    newtree->Branch("Hcorr", "TLorentzVector", &hp4corr);
    TLorentzVector hzp4corr;
    newtree->Branch("HZcorr", "TLorentzVector", &hzp4corr);

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
    Double_t Lgen;
    int nEvt = 0;
    for(Int_t entry = 0; entry < genEntries; ++entry){
        // eLpR = 0; eRpL = 1; eLpL = 2; eRpR = 3
        //std::cout << Iproc << "\n";
        if(Iproc >= 0){
            Lgen = Lgen_bg_arr[Ipol][Iproc];
        }
        else{
            //std::cout << "else" << "\n";
            Lgen = Lgen_sig_bbll_arr[Ipol][Ms];
        }
        myheader.Ipol = Ipol;
        myheader.Lgen = Lgen;
        myheader.Iproc = Iproc;
        genReader->ReadEntry(entry);
        //std::cout << "testtttt\n";

        int nParticle = branchParticle->GetEntriesFast();
        int nMuon = branchMuon->GetEntriesFast();
        int nElectron = branchElectron->GetEntriesFast();
        int nPhoton = branchPhoton->GetEntriesFast();
        int nBcal = branchBCal->GetEntriesFast();

        int nEle = 0;
        TLorentzVector myelectron[4];
        int myqel[4];

        int Qtot = 0;

        // Photons

        int nIsr = 0;
        int nPh = 0;
        
        p4isr.SetPxPyPzE(0.,0.,0.,0.);
        p4ph.SetPxPyPzE(0.,0.,0.,0.);

        // Loop over photons in BCAL
        
        for(int iPh=0; iPh<nBcal; iPh++){
            photon = (Photon*) branchBCal->At(iPh);
            // All BCAL photons are considered ISR photons
            nIsr++;
            p4isr+=photon->P4();
        }
        // Loop over photons in calorimeter    
        for(int iPh=0; iPh<nPhoton; iPh++){
            photon = (Photon*) branchPhoton->At(iPh);

            // ISR photon candidate selection
            
            if(abs(photon->Eta) > 2. && photon->E > 10.){
                nIsr++;
                p4isr+=photon->P4();
            }
            else{
                nPh++;
                p4ph+= photon->P4();
            }
        }
        if ( ( 1<<nIsr & imask)  == 0) continue;
        //Zapytać o selekcję isr
        for(int iEle = 0; iEle<nElectron; iEle++){
            electron = (Electron*) branchElectron->At(iEle);

            if(abs(electron->Eta) > 2. || electron->PT < 5.) continue;

            myelectron[nEle] = electron->P4();
            myqel[nEle] = electron->Charge;
            Qtot += electron->Charge;
            nEle++;

            if(nEle==4) break;
        }
        //if ( ( 1<<nEle & emask)  == 0) continue;

        myevent.Nel = nEle;

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
        
        myevent.Nmu = nMu;
        //if ( ( 1<<nMu & mmask)  == 0) continue;


        // At most two leptons allowed

        if ( nEle + nMu > 2 ) continue;
        
        //Only pairs of leptons allowed
        if ( nEle != 2 && nMu != 2) continue;
        
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
        int nBquark = 0;
        int itmp;
        
        for(int ijet=0;ijet<nJet;ijet++)
        {
            Jet* j1= (Jet*) branchJet->At(ijet);

            myjet[ijet] = j1->P4();
            mytag[ijet] = j1->TauTag;
            mybtag[ijet] = j1->BTag;
            nTau+=j1->TauTag;
            nBquark+=j1->BTag;

            if(ijet==0)
            {
                myevent.y23 = j1->ExclYmerge23; 
                myevent.y34 = j1->ExclYmerge34;
                myevent.y45 = j1->ExclYmerge45;
            }
        }

        // Select only given tag configuration

        //if ( ( 1<<nBquark & bmask)  == 0) continue;

        
        // Check combined tag number

        int ntag = nEle + nMu + nBquark;

        //if ( ( 1<<ntag & amask)  == 0) continue;

        
        // Count selected events only
        
        nEvt++;
        bool debug = ( nEvt<=5 || nEvt%1000==0 ) && (plname.length() == 0); 
        
        
        // Sort jets. Those with tau tag go first. For same tag value, lower masses go first
        
        for(int ijet=1;ijet<nJet;)
        {
            if(mybtag[ijet] > mybtag[ijet-1] || (mybtag[ijet] == mybtag[ijet-1] && myjet[ijet].M() < myjet[ijet-1].M()) )
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

        if (nElectron == 2){
            for(int iEle = 0; iEle < nEle; iEle++){
                myrec[iEle] = myelectron[iEle];
            }
        }
        else{
            for(int iMu = 0; iMu < nMu; iMu++){
                myrec[iMu] = mymuon[iMu];
            }
            
            if(myrec[0].E()+myrec[1].E() < 30){
                std::cout << "r";
                continue;
            }
        }
        for(int ijet = 0; ijet<nJet; ijet++){
            myrec[nEle+nMu+ijet] = myjet[ijet];
        }
        
        myevent.Nisr = nIsr;
        myevent.Nph = nPh;
        myevent.jet1btag = (mybtag[0] & Bbit)/Bbit;
        myevent.jet2btag = (mybtag[1] & Bbit)/Bbit;
        
        zp4 = myrec[0] + myrec[1];
        h1p4 = myrec[2] + myrec[3];

        myevent.Mjj = h1p4.M();

        hzp4 = h1p4 + zp4;
        

        TLorentzVector temprecoilhp4(0, 0, 0, 250);
        temprecoilhp4 = temprecoilhp4 - zp4;
        recoilhp4 = temprecoilhp4;
        myevent.Mrec = recoilhp4.M();

        Double_t sqrt_s = 250;
        Double_t alpha = 0;
        Double_t pllx = myrec[0].Px() + myrec[1].Px();
        Double_t plly = myrec[0].Py() + myrec[1].Py();

        Double_t px, py, pt, phi1, phi2, phi12, phi, theta1, theta2, p1, p2;
        px = -pllx + sqrt_s * TMath::Sin(alpha/2);
        py = -plly;

        pt = TMath::Sqrt(px*px + py*py);

        phi = TMath::ATan2(py,px);

        phi1 = myrec[2].Vect().Phi();
        theta1 = myrec[2].Vect().Theta();

        phi2 = myrec[3].Vect().Phi();
        theta2 = myrec[3].Vect().Theta();

        phi12 = phi1 - phi2;

        p1 = pt*TMath::Sin(phi - phi2)/(TMath::Sin(phi12)*TMath::Sin(theta1));
        p2 = pt*TMath::Sin(phi1 - phi)/(TMath::Sin(phi12)*TMath::Sin(theta2));

        Double_t E1, E2;
        Double_t m1 = myrec[2].M();
        Double_t m2 = myrec[3].M();

        E1 = TMath::Sqrt(p1*p1 + m1*m1);
        E2 = TMath::Sqrt(p2*p2 + m2*m2);

        Double_t jet1pmag = myrec[2].P();
        Double_t jet2pmag = myrec[3].P();

        myrec[2] = TLorentzVector(myrec[2].X() * p1 / jet1pmag, myrec[2].Y() * p1 / jet1pmag, myrec[2].Z() * p1 / jet1pmag, E1);
        myrec[3] = TLorentzVector(myrec[3].X() * p2 / jet2pmag, myrec[3].Y() * p2 / jet2pmag, myrec[3].Z() * p2 / jet2pmag, E2);
        hp4corr = myrec[2] + myrec[3];
        hzp4corr = zp4 + hp4corr;
        if (TMath::Abs(phi12) > 3) continue;
        
        myevent.Mcorr = TMath::Sqrt((1+p2/p1)*m1*m1+(1+p1/p2)*m2*m2+2*pt*pt*TMath::Sin(phi-phi2)*TMath::Sin(phi1-phi)*(1-TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Cos(phi12)-TMath::Cos(theta1)*TMath::Cos(theta2))/(TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Sin(phi12)*TMath::Sin(phi12)));
        myevent.Etot = zp4.E() + E1 + E2;
        if(myevent.Etot > 500) continue;
        if(myevent.Mcorr != myevent.Mcorr) continue;
        //std::cout << myevent.Nmu;
        /*
        if(myevent.Mcorr > 300){
            std::cout << "phi1: " << phi1 << " phi2: " << phi12 << " mcorr: " << myevent.Mcorr << "\n";
        }*/
        if(myevent.Mcorr < 0) continue;
        TVector3* SCM_boostvec = new TVector3(-h1p4.Px()/h1p4.E(), -h1p4.Py()/h1p4.E(), -h1p4.Pz()/h1p4.E());
        TLorentzVector jet1P4_SCM= myrec[2];
        TLorentzVector jet2P4_SCM= myrec[3];
        jet1P4_SCM.Boost(*SCM_boostvec);
        jet2P4_SCM.Boost(*SCM_boostvec);
        myevent.costhetajetcms = jet1P4_SCM.CosTheta();
        if(myrec[2].Pt() > 250) continue;
        myevent.jetpt = myrec[2].Pt();
        //jetpt_branch = myrec[2].Pt();
        //if(entry%500 == 0) std::cout << myrec[2].Pt() << " ";
        
        newtree->Fill();

    }
    std::cout <<  "\n";
    fout->Write();
    fout->Close();
}