#include "TFile.h"
#include "TTree.h"
#include <RtypesCore.h>
#include <fstream>
#include <iostream>

void get_preselec_eff(std::string prefile, std::string postfile){
    std::cout << "inprefile: " << prefile << "\n";
    std::cout << "inpostfile: " << postfile << "\n";
    TFile* presel = new TFile(prefile.c_str(), "OPEN");
    TFile* postsel = new TFile(postfile.c_str(), "OPEN");

    TTree* pretree = (TTree*)presel->FindObjectAny("Delphes");
    TTree* posttree = (TTree*)postsel->FindObjectAny("tree");

    Long64_t nr_before = pretree->GetEntries();
    Long64_t nr_after = posttree->GetEntries();

    Double_t eff = (Double_t)nr_after/nr_before;

    std::ofstream outfile;

    std::string outfilename = prefile.substr(11, prefile.size()-5);
    outfilename = "efficiencies/" + outfilename + ".txt";

    outfile.open(outfilename);

    outfile << eff << "\n" << nr_before << "\n" << nr_after;
    outfile.close();
}