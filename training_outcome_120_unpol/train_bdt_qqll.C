/**********************************************************************************
 *                                                                                *
 * Based on TMVAClassification.C - TMVA tutororial example                        *
 * General routine for training and testing BDT method                            *
 * A.F.Zarnecki                                                                   *
 *                                                                                *
 **********************************************************************************/

#include <RtypesCore.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/mvaeffs.h"
#include "TMVA/variables.h"
#include "TMVA/mvas.h"

using namespace std;
string treename = "tree";         // tree name (should be the same in all input files)


            void train_bdt_qqll(string variables = "Nel Nmu Nisr Nph jet1btag jet2btag Mjj Mcorr Mrec Etot log10(y23) log10(y34) log10(y45) costhetajetcms jetpt  Z.M() Z.Pt() Z.E()", // input BDT variable names
			   string spectator = "Iproc Ipol Lgen", // spectator variable names (for output tree)
			   double nsgen= 100.0,     // expected signal events
			   double nbgen=1000.0,     // expected background events
			   string weight =  "w",     // event weigth formula
			   string outfile = "train_bdt_qqll", // output file name root
			   string Ntrees = "1000",  // number of trees
			   string transf = "I",     // input variable transformation options
			   string bdtpar = "UseNvars=6:PruneMethod=CostComplexity:PruneStrength=-1",
         int beamsetting = 4) // BDT options
{
   //---------------------------------------------------------------
   // This loads the library
   
   TMVA::Tools::Instance();

   // to get access to the GUI and all tmva macros
   //    TString tmva_dir(TString(gRootDir) + "/tmva");
   //   if(gSystem->Getenv("TMVASYS"))
   //      tmva_dir = TString(gSystem->Getenv("TMVASYS"));
   //   gROOT->SetMacroPath(tmva_dir + "/test/:" + gROOT->GetMacroPath() );
   //
   //  gROOT->SetMacroPath(TString("test/:") + gROOT->GetMacroPath() );
   //  gROOT->ProcessLine(".L TMVAGui.C");
   //  gROOT->ProcessLine(".L efficiencies.C");

   
   // Default MVA methods to be trained + tested
   
   std::map<std::string,int> Use;

 // --- Boosted Decision Trees

   Use["BDT"]             = 1; // uses Adaptive Boost
   
 /*  
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
 */

 // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;


   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   
   string outroot = outfile + ".root";
   
   TFile* outputFile = TFile::Open( outroot.c_str(), "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is 
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string


   /*  V1
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                 "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );
   */


   string factory_opt = "!V:!Silent:Color:DrawProgressBar:Transformations=" + transf + ":AnalysisType=Classification";

   TMVA::Factory *factory = new TMVA::Factory( outfile.c_str(), outputFile, factory_opt.c_str() );

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");


   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   

   TString varlist(variables);
   TObjArray* vararray = varlist.Tokenize(" ");
   int nVar = vararray->GetEntries();

   for(int iVar=0; iVar<nVar; iVar++)
     {
       string varName(((TObjString *)vararray->At(iVar))->String());
       cout << "Adding variable: " << varName << endl;
       dataloader->AddVariable( varName.c_str(), "", "", 'F' ); 
     }
   
   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   //dataloader->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
   //dataloader->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );

   TString speclist(spectator);
   TObjArray* specarray = speclist.Tokenize(" ");
   int nSpec = specarray->GetEntries();

   for(int iSpec=0; iSpec<nSpec; iSpec++)
     {
       string specName(((TObjString *)specarray->At(iSpec))->String());
       cout << "Adding spectator variable: " << specName << endl;
       dataloader->AddSpectator( specName.c_str(), "", "", 'F' ); 
     }
   
  
   // --- Register the training and test trees
 
   // global event weights per tree (see below for setting event-wise weights)


   Double_t weights[5][4] = {
    {0.585, 0.035, 0.315, 0.065},
    {0.035, 0.585, 0.065, 0.315},
    {0.315, 0.065, 0.585, 0.035},
    {0.065, 0.315, 0.035, 0.585},
    {0.25, 0.25, 0.25, 0.25}};

    Double_t Lexp[5][4];

    for(int i = 0; i < 5; i++){
      for(int j = 0; j < 4; j++){
        if(i < 2){
          Lexp[i][j] = 900 * weights[i][j];
        }
        else if(i < 4){
          Lexp[i][j] = 100 * weights[i][j];
        }
        else{
          Lexp[i][j] = 2000 * weights[i][j];

        }
      }
    }
   
   // Define signal and background trees
        struct event_head {
        Int_t Ievt;
        Int_t Iproc;
        Int_t Ipol;
        Float_t Ms;
        Float_t Cs;
        Float_t Lgen; // Lgen
    };
    event_head header;
    TString path = "../trees_for_training/";
    TFile* sigfile = new TFile(path + "bbll_sig_120_eRpL_new.root");
    TFile* sigfile2 = new TFile(path + "bbll_sig_120_eLpR_new.root");
        
    if (!sigfile)
          {
           cout << "==> Abort, could not open signal event file " << endl;
           return;
          }
  

   TTree *sigev = (TTree *)(sigfile->FindObjectAny(treename.c_str()));
   TTree *sigev2 = (TTree *)(sigfile2->FindObjectAny(treename.c_str()));
   Int_t Ipol;
   Double_t Lgen;
   sigev->SetBranchAddress("header", &header);
   
   sigev->GetEntry(0);
   Ipol = header.Ipol;
   Lgen = header.Lgen;
   dataloader->AddSignalTree( sigev,  Lexp[beamsetting][Ipol]/Lgen,   TMVA::Types::kMaxTreeType );

   sigev2->SetBranchAddress("header", &header);
   
   sigev2->GetEntry(0);
   Ipol = header.Ipol;
   Lgen = header.Lgen;
   dataloader->AddSignalTree( sigev2,  Lexp[beamsetting][Ipol]/Lgen,   TMVA::Types::kMaxTreeType );

  TFile* bg_files[26];
  TTree* bgevt[26];
  TString bg_filenames[26] = {
    "qq_bg_eRpL_new.root", "qqll_bg_eRpL_new.root", "qqlv_bg_eRpL_new.root", "qqqq_bg_eRpL_new.root",
    "qqtt_bg_eRpL_new.root", "qqtv_bg_eRpL_new.root", "qqvv_bg_eRpL_new.root",
    "ttll_bg_eRpL_new.root", "tttt_bg_eRpL_new.root",
    "qq_bg_eLpR_new.root", "qqll_bg_eLpR_new.root", "qqlv_bg_eLpR_new.root", "qqqq_bg_eLpR_new.root",
    "qqtt_bg_eLpR_new.root", "qqtv_bg_eLpR_new.root", "qqvv_bg_eLpR_new.root",
    "ttll_bg_eLpR_new.root", "tttt_bg_eLpR_new.root", "qqll_bg_eLpL_new.root", "qqlv_bg_eLpL_new.root", "ttll_bg_eLpL_new.root",
    "qqll_bg_eRpR_new.root", "qqlv_bg_eRpR_new.root", "ttll_bg_eRpR_new.root", "qqllvv_bg_eRpL_new.root", "qqllvv_bg_eLpR_new.root"  };

  for(int i = 0; i < 26; i++){
    bg_files[i] = new TFile(path + bg_filenames[i]);
    if (!bg_files[i])
      {
        cout << "==> Abort, could not open background event file " << endl;
        return;
      }
    bgevt[i] = (TTree *)(bg_files[i]->FindObjectAny(treename.c_str()));
    if(bgevt[i]->GetEntries() == 0) continue;
    bgevt[i]->SetBranchAddress("header", &header);
    
    bgevt[i]->GetEntry(0);
    Ipol = header.Ipol;
    Lgen = header.Lgen;
    //std::cout << Lexp[beamsetting][Ipol]/Lgen << "\n";
    dataloader->AddBackgroundTree( bgevt[i], Lexp[beamsetting][Ipol]/Lgen,  TMVA::Types::kMaxTreeType  );
    
  }


	 
   // From User Guide
   
   //dataloader->SetWeightExpression( "Ipol" );

   // Apply additional cuts on the signal and background samples (can be different)
   
   TCut mycuts =  "";
   
   TCut mycutb =  "";  // Should in principle be empty - all events contribute
   
   // Tell the factory how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
   //    dataloader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   // To also specify the number of testing events, use:
   //    dataloader->PrepareTrainingAndTestTree( mycut,
   //                          "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
  
   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // ---- Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable
  


   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT",
                 "!H:!V:NTrees=" + Ntrees + ":MinNodeSize=3%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseRandomisedTrees=True:SeparationType=GiniIndex:nCuts=100" + bdtpar );


   // For an example of the category classifier usage, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","GA");

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;
   delete dataloader;

   // Launch variable distributions and significance estimate
   
   TCanvas  *canvas1 = (TCanvas *) gROOT->FindObject("canvas1");
   if(canvas1)  gROOT->Delete("canvas1");
   
   TCanvas  *canvas2 = (TCanvas *) gROOT->FindObject("canvas2");
   if(canvas2)  gROOT->Delete("canvas2");

   
   TString fin=outroot.c_str();
   TString dataset="dataset";
   Bool_t useTMVAStyle = kFALSE;

   gStyle->SetOptLogy(0);
  
   TMVA::variables(dataset,fin,"InputVariables_Id","TMVA Input Variables", kFALSE, useTMVAStyle);
   
   canvas1 = (TCanvas *) gROOT->FindObject("canvas1");

   string outpdf = outfile + "_var.pdf";
   string outpng = outfile + "_var.png";
   string outc   = outfile + "_var.C";

   canvas1->Print(outpdf.c_str());
   canvas1->Print(outpng.c_str());
   canvas1->Print(outc.c_str());

   canvas2 = (TCanvas *) gROOT->FindObject("canvas2");

   if(canvas2)
     {
      outpdf = outfile + "_var2.pdf";
      outpng = outfile + "_var2.png";
      outc   = outfile + "_var2.C";

      canvas2->Print(outpdf.c_str());
      canvas2->Print(outpng.c_str());
      canvas2->Print(outc.c_str());

      gROOT->Delete("canvas2");
     }
   
   gROOT->Delete("canvas1");

   TMVA::HistType htype =  TMVA::kCompareType;
   
   TMVA::mvas(dataset,fin,htype,useTMVAStyle);

   canvas1 = (TCanvas *) gROOT->FindObject("canvas1");

   outpdf = outfile + "_res.pdf";
   outpng = outfile + "_res.png";
   outc   = outfile + "_res.C";

   canvas1->Print(outpdf.c_str());
   canvas1->Print(outpng.c_str());
   canvas1->Print(outc.c_str());

 
   TMVA::TMVAGlob::Initialize( useTMVAStyle );
   
   TString formula="S/sqrt(S+B)";


 // Number of events expected in each tree

   std::cout << "Number of expected signal events for significance calculation: " << nsgen << std::endl;

   std::cout << "Number of expected background events for significance calculation: " << nbgen << std::endl;
   
   TMVA::StatDialogMVAEffs* gGui = new TMVA::StatDialogMVAEffs(dataset,gClient->GetRoot(), nsgen, nbgen);

   TFile* file = TMVA::TMVAGlob::OpenFile( fin );
   
   gGui->ReadHistograms(file);
   gGui->SetFormula(formula);
   gGui->UpdateSignificanceHists();
   gGui->DrawHistograms();
   
   //   gGui->RaiseDialog();   

   canvas1 = (TCanvas *) gROOT->FindObject("canvas1");

   outpdf = outfile + "_sig.pdf";
   outpng = outfile + "_sig.png";
   outc   = outfile + "_sig.C";

   canvas1->Print(outpdf.c_str());
   canvas1->Print(outpng.c_str());
   canvas1->Print(outc.c_str());
  

   // Launch the GUI for the root macros
   
    if (!gROOT->IsBatch()) TMVA::TMVAGui( outroot.c_str() );

}
