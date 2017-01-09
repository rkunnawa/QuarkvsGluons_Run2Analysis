#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

void TMVAClassification( ) { 
  
  // Create a new root output file.
  TString outfileName = "QGTest_Discriminant.root";
  TFile*  outputFile = TFile::Open(outfileName,"RECREATE");

  TMVA::Factory *factory = new TMVA::Factory( "TMVA", outputFile, "!V:!Silent");

  factory->AddVariable( "jtnCands", "Number of Constituents", "#", 'I' );
  factory->AddVariable( "jtnChCands", "Number of Charged Constituents", "#", 'I' );
  factory->AddVariable( "jtnNeCands", "Number of Neutral Constituents", "#", 'I' );
  factory->AddVariable( "jtm", "Mass", "GeV/c2", 'F' );
  factory->AddVariable( "jtMByPt", "Mass over pT", "1/c", 'F' );
  factory->AddVariable( "jtRMSCand", "RMS of constituents pT", "#", 'F' );
  factory->AddVariable( "jtAxis1", "Major axis", "#", 'F' );
  factory->AddVariable( "jtAxis2", "Minor axis", "#", 'F' );
  factory->AddVariable( "jtSigma", "Sigma of constituents pT", "#", 'F' );
  factory->AddVariable( "jtR", "distribution of constituents", "#", 'F' );
  factory->AddVariable( "jtpTD", "pTD: sum square const pT/ jet pT", "#", 'F' );
  factory->AddVariable( "jtpull", "pull", "#", 'F' );
  factory->AddVariable( "jtrm0p5", "radial moment (0.5)", "#", 'F' );
  factory->AddVariable( "jtrm1", "first radial moment", "#", 'F' );
  factory->AddVariable( "jtrm2", "second radial moment", "#", 'F' );
  factory->AddVariable( "jtrm3", "third radial moment", "#", 'F' );
  
  // read training and test data
  // load the signal and background event samples from ROOT trees

  TFile *fin = TFile::Open("HiForestAOD.root"); 
   
  TTree *inputTree     = (TTree*)fin->Get("ak4PFJetAnalyzer/t");

  TCut sig = "subid==0 && refdrjt<0.4 && refparton_flavor>=-6 && refparton_flavor<=6 && refpt>80 && fabs(jteta)<2.0";
  TCut bkg = "subid==0 && refdrjt<0.4 && refparton_flavor==21 && refpt>80 && fabs(jteta)<2.0";

  factory->SetInputTrees(inputTree, sig, bkg);

  // tell the factory to use all remaining events in the trees after training for testing:
  factory->PrepareTrainingAndTestTree( "","", "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

  // factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:" );
  // factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood", "H:!V:" );
  factory->BookMethod( TMVA::Types::kBDT, "BDTDefault", "H:!V");
  // factory->BookMethod( TMVA::Types::kMLP, "MLPDefault", "H:!V:" );

  // Train MVAs using the set of training events
  factory->TrainAllMethods();

  // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();    

  // --------------------------------------------------------------
   
  // Save the output
  outputFile->Close();

  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) {
    if (gROOT->GetVersionInt() > 60402) TMVA::TMVAGui( outfileName ); 
    //     else                                {gROOT->ProcessLine(".L TMVAGui.C"); TMVAGui(outfileName);}
  } 
  delete factory;
}
