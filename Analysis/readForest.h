//!
//! Raghav Kunnawalkam Elayavalli
//! Nov 3rd 2016 
//!
//! Macro to read forest files and produce root histograms for the radial moment analysis.
//! Basically we need to plot the first and second radial moment for bachward vs forward
//! and for different pT and Ntrk bins 
//!
//!


// C++, C, etc.
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
// ROOT
#include <TSystem.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
// Fastjet 
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h" 


const int minArgs = 1;

const std::string filelist = "filelists/forests.txt";
const int startFile = 1;
const int endFile = 2;
const std::string jetalgo = "ak";
const int radius = 4;
const std::string jettype = "PF";
const bool printDebug = true;
const std::string outputName = "MC_variables.root";

int readForest_MC(std::string inFilelist = filelist,
		  int startfile = startFile,
		  int endfile = endFile,
		  std::string jetAlgo = jetalgo,
		  int radius = radius,
		  std::string jetType = jettype,
		  bool debugMode = printDebug,       
		  std::string outfile = outputName);
int readForest_Data(std::string inFilelist = filelist,
		    int startfile = startFile,
		    int endfile = endFile,
		    std::string jetAlgo = jetalgo,
		    int radius = radius,
		    std::string jetType = jettype,
		    bool debugMode = printDebug,       
		    std::string outfile = outputName);

const int readForestsArgCount = 8+minArgs;

double deltaphi(double phi1, double phi2){
  double pi=TMath::Pi();
  double dphi=TMath::Abs(phi1-phi2);
  if(dphi > pi)dphi -=2.*pi;
  return TMath::Abs(dphi);
}

double deltaR(double eta1, double phi1, double eta2, double phi2){
  double delR =TMath::Sqrt((eta1-eta2)*(eta1-eta2) + deltaphi(phi1, phi2)*deltaphi(phi1, phi2));
  return delR;
}

double getPX(double pt, double phi){
  double px = pt * TMath::Cos(phi);
  return px;
}
double getPY(double pt, double phi){
  double py = pt * TMath::Sin(phi);
  return py;
}
double getPZ(double pt, double eta){
  double pz = pt * TMath::SinH(eta);
  return pz;
}

bool isinMatchCollection(std::vector<unsigned> coll, unsigned location){
  bool decision = false;
  for(unsigned ic = 0; ic<coll.size(); ++ic){
    if(location == coll.at(ic))
      decision = true;
  }
  return decision;
}

std::vector<unsigned> matchedJets(std::vector<fastjet::PseudoJet> recoJets, std::vector<fastjet::PseudoJet> genJets){

  std::vector<unsigned> matGenJet;
  /* matGenJet.resize(recoJets.size()); */
  for(unsigned rj = 0; rj<recoJets.size(); ++rj){
    for(unsigned gj = 0; gj<genJets.size(); ++gj){
      if(deltaR(recoJets.at(rj).eta(), recoJets.at(rj).phi(), genJets.at(gj).eta(), genJets.at(gj).phi()) < 0.3){
	if(matGenJet.size()>0 && isinMatchCollection(matGenJet, gj))
	  continue;
	else
	  matGenJet.push_back(gj);
      }
    }
  }
  return matGenJet;
}

std::string dir_MC[] = {
  //! jet tree
  Form("%s%d%sJetAnalyzer/t", jetalgo.c_str(), radius, jettype.c_str()),
  //! event selection 
  "hltanalysis/HltTree", "skimanalysis/HltTree", "hiEvtAnalyzer/HiTree",
  //! objects for clustering/constituents 
  "pfcandAnalyzer/pfTree", "ppTrack/trackTree"
  //! MC related trees 
  , "HiGenParticleAna/hi", "runAnalyzer/run" 
};
const int _ntrees_MC = sizeof(dir_MC)/sizeof(std::string);

std::string dir_Data[] = {
  //! jet tree
  Form("%s%d%sJetAnalyzer/t", jetalgo.c_str(), radius, jettype.c_str()),
  //! event selection 
  "hltanalysis/HltTree", "skimanalysis/HltTree", "hiEvtAnalyzer/HiTree",
  //! objects for clustering/constituents 
  "pfcandAnalyzer/pfTree", "ppTrack/trackTree"
};
const int _ntrees_Data = sizeof(dir_Data)/sizeof(std::string);
