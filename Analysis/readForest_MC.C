//!
//! Raghav Kunnawalkam Elayavalli
//! Nov 3rd 2016 
//!
//! Macro to read forest files and produce root histograms for the radial moment analysis.
//! Basically we need to plot the first and second radial moment for bachward vs forward
//! and for different pT and Ntrk bins 
//! this takes either MC or Data.
//! for MC there needs to be comparison with gen moments and possbly for unfolding
//! for MC there should also be quark vs gluon 
//! general analysis structure:
//! need to do corrections on the fly! 
//!

#include "readForest.h"

using namespace std;
using namespace fastjet;

int readForest_MC(std::string fileList,
		  int startfile,
		  int endfile,
		  std::string jetAlgo,
		  int radius,
		  std::string jetType,
		  bool debugMode,
		  std::string outfile){
  
  TStopwatch timer;
  timer.Start();
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;

  std::string dataset = "MC";
  
  std::cout<<"reading filelist "<< fileList << std::endl;
  std::cout<<"Running on " << dataset << std::endl;
  std::cout<<"reading files #'s "<< startfile << " to " << endfile<<std::endl;
  std::cout<<"jet Algo = " << jetAlgo;
  std::cout<<", radius = " << radius;
  std::cout<<", jetType = " << jetType;
  std::cout<<"debugMode = "<<debugMode<<std::endl;
  
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  //! declare the output tree for clean events
  //! List of variables:
  //! centrality
  //! jet gen parton flavor 
  //! jet pT
  //! eta
  //! phi
  //! mass
  //! radial moment - 1
  //! radial moment - 2
  //! radial moment - 0.5
  //! number of charged constituents  

  TTree *jetTree = new TTree("jetTree","QG Discriminant Tree");
  jetTree->Branch("vz",&vz,"vz/F");
  jetTree->Branch("hibin",&hibin,"hibin/I");
  jetTree->Branch("pthat",&pthat,"pthat/F");
  jetTree->Branch("weight",&weight,"weight/F");
  jetTree->Branch("njets",&njets,"njets/I");
  jetTree->Branch("jetpt",jetpt,"jetpt[njets]/F");
  jetTree->Branch("jeteta",jeteta,"jeteta[njets]/F");
  jetTree->Branch("jetphi",jetphi,"jetphi[njets]/F");
  jetTree->Branch("jetmass",jetmass,"jetmass[njets]/F");
  jetTree->Branch("jetflav",jetflav,"jetflav[njets]/I");
  jetTree->Branch("nPFCand",nPFCand,"nPFCand[njets]/I");  
  jetTree->Branch("nCharged",nCharged,"nCharged[njets]/I");  
  jetTree->Branch("nNeutral",nNeutral,"nNeutral[njets]/I");  
  jetTree->Branch("nPFMult",nPFMult,"nPFMult[njets]/I");  
  jetTree->Branch("jetMByPt",jetMByPt,"jetMByPt[njets]/F");  
  jetTree->Branch("RMSCand",RMSCand,"RMSCand[njets]/F");  
  jetTree->Branch("Axis1",Axis1,"Axis1[njets]/F");  
  jetTree->Branch("Axis2",Axis2,"Axis2[njets]/F");  
  jetTree->Branch("Sigma",Sigma,"Sigma[njets]/F");  
  jetTree->Branch("R",R,"R[njets]/F");  
  jetTree->Branch("pTD",pTD,"pTD[njets]/F");  
  jetTree->Branch("pull",pull,"pull[njets]/F");  
  
  //! setup the analysis framework
  std::ifstream instr_Forest(fileList.c_str(),std::ifstream::in);
  std::string filename_Forest;
  
  if(debugMode)std::cout<<"reading from "<<startfile<<" to "<<endfile<<std::endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }


  TChain * jetTree[_ntrees_MC];

  for(int t = 0;t<_ntrees_MC;t++){
    jetTree[t] = new TChain(std::string(dir_MC[t]).data());
  }//!tree loop ends
  
  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;

    for(int t = 0; t<_ntrees_MC; ++t){
      jetTree[t]->Add(filename_Forest.c_str());    
      if(debugMode)std::cout << "Tree loaded  " << std::string(dir_MC[t]).data() << std::endl;
      if(debugMode)std::cout << "Entries : " << jetTree[t]->GetEntries() << std::endl;
    }
    std::cout<<"Total number of events loaded in HiForest = "<<jetTree[2]->GetEntries()<<std::endl;

  }

  //! friend the trees
  for(int t = 1; t<_ntrees_MC; ++t)
    jetTree[0]->AddFriend(jetTree[t]);

  //! Variables for HiForest
  
  //! jet info
  int nref_F;
  float pthat_F;
  int   subid_F[1000];
  float refdrjt_F[1000];
  int refparton_F[1000];
  float refpt_F[1000];
  float geneta_F[1000];
  float pt_F[1000];
  float eta_F[1000];
  float phi_F[1000];
  float rawpt_F[1000];
  float jtpu_F[1000];
  float chMax_F[1000];
  float chSum_F[1000];
  int chN_F[1000];
  // float chHardSum_F[1000];
  float trkMax_F[1000];
  float trkSum_F[1000];
  // float trkHardSum_F[1000];
  float phMax_F[1000];
  float phSum_F[1000];
  // float phHardSum_F[1000];
  float neSum_F[1000];
  float neMax_F[1000];
  int neN_F[1000];
  float eMax_F[1000];
  float eSum_F[1000];
  float muMax_F[1000];
  float muSum_F[1000];
  jetTree[0]->SetBranchAddress("nref",&nref_F);
  jetTree[0]->SetBranchAddress("pthat",&pthat_F);
  jetTree[0]->SetBranchAddress("subid",subid_F);
  jetTree[0]->SetBranchAddress("refdrjt",refdrjt_F);
  jetTree[0]->SetBranchAddress("refparton_flavor",refparton_F);
  jetTree[0]->SetBranchAddress("refpt",refpt_F);
  jetTree[0]->SetBranchAddress("refeta",geneta_F);
  jetTree[0]->SetBranchAddress("jtpt",pt_F);
  jetTree[0]->SetBranchAddress("jteta",eta_F);
  jetTree[0]->SetBranchAddress("jtphi",phi_F);
  jetTree[0]->SetBranchAddress("rawpt",rawpt_F);
  jetTree[0]->SetBranchAddress("jtpu",jtpu_F);
  jetTree[0]->SetBranchAddress("chargedMax",&chMax_F);
  jetTree[0]->SetBranchAddress("chargedSum",&chSum_F);
  jetTree[0]->SetBranchAddress("chargedN",&chN_F);
  // jetTree[0]->SetBranchAddress("chargedHardSum",&chSum_F);
  jetTree[0]->SetBranchAddress("trackMax",&trkMax_F);
  jetTree[0]->SetBranchAddress("trackSum",&trkSum_F);
  // jetTree[0]->SetBranchAddress("trackHardSum",&trkSum_F);
  jetTree[0]->SetBranchAddress("photonMax",&phMax_F);
  jetTree[0]->SetBranchAddress("photonSum",&phSum_F);
  // jetTree[0]->SetBranchAddress("photonHardSum",&phSum_F);
  jetTree[0]->SetBranchAddress("neutralMax",&neMax_F);
  jetTree[0]->SetBranchAddress("neutralSum",&neSum_F);
  jetTree[0]->SetBranchAddress("neutralN",&neN_F);
  jetTree[0]->SetBranchAddress("eSum",eSum_F);
  jetTree[0]->SetBranchAddress("eMax",eMax_F);
  jetTree[0]->SetBranchAddress("muSum",muSum_F);
  jetTree[0]->SetBranchAddress("muMax",muMax_F);

  //! hlt analysis 
  int minbiasBit, minbiasBit_p;
  int jet40, jet40_p;
  int jet60, jet60_p;
  int jet80, jet80_p;
  int jet100, jet100_p;
  jetTree[1]->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_v1", &minbiasBit);
  jetTree[1]->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_v1_Prescl", &minbiasBit_p);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet40_Eta5p1_v2",&jet40);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet40_Eta5p1_v2_Prescl",&jet40_p);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v2",&jet60);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v2_Prescl",&jet60_p);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v2",&jet80);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v2_Prescl",&jet80_p);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet100_Eta5p1_v2",&jet100);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet100_Eta5p1_v2_Prescl",&jet100_p);
  
  //! noise filters 
  int pBeamScrapingFilter_F;
  int pPAprimaryVertexFilter_F;
  int pHBHENoiseFilter_F;
  jetTree[2]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_F);
  jetTree[2]->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter_F);
  jetTree[2]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_F);

  //! event information 
  ULong64_t evt_F;
  UInt_t run_F;
  UInt_t lumi_F;
  float vz_F;
  int hiNpix_F;
  int hiNpixelTracks_F;
  int hiBin_F;
  float hiHF_F;
  int hiNtracks_F;
  int hiNtracksPtCut_F;
  int hiNtracksEtaCut_F;
  int hiNtracksEtaPtCut_F;
  jetTree[3]->SetBranchAddress("evt",&evt_F);
  jetTree[3]->SetBranchAddress("run",&run_F);
  jetTree[3]->SetBranchAddress("lumi",&lumi_F);
  jetTree[3]->SetBranchAddress("vz",&vz_F);
  jetTree[3]->SetBranchAddress("hiBin",&hiBin_F);
  jetTree[3]->SetBranchAddress("hiHF", &hiHF_F);
  jetTree[3]->SetBranchAddress("hiNpix",&hiNpix_F);
  jetTree[3]->SetBranchAddress("hiNpixelTracks",&hiNpixelTracks_F);
  jetTree[3]->SetBranchAddress("hiNtracks",&hiNtracks_F);
  jetTree[3]->SetBranchAddress("hiNtracksPtCut",&hiNtracksPtCut_F);
  jetTree[3]->SetBranchAddress("hiNtracksEtaCut",&hiNtracksEtaCut_F);
  jetTree[3]->SetBranchAddress("hiNtracksEtaPtCut",&hiNtracksEtaPtCut_F);

  //! Particle flow candidates: vectors
  int nPart_F;
  std::vector<int> *pfID_F = 0;
  std::vector<float> *pfEnergy_F = 0;
  std::vector<float> *pfPt_F = 0;
  std::vector<float> *pfEta_F = 0;
  std::vector<float> *pfPhi_F = 0;
  jetTree[4]->SetBranchAddress("nPFpart", &nPart_F);
  jetTree[4]->SetBranchAddress("pfId", &pfID_F);
  jetTree[4]->SetBranchAddress("pfEnergy", &pfEnergy_F);
  jetTree[4]->SetBranchAddress("pfPt", &pfPt_F);
  jetTree[4]->SetBranchAddress("pfEta", &pfEta_F);
  jetTree[4]->SetBranchAddress("pfPhi", &pfPhi_F);
  
  //! Track trees
  int nTrk_F;
  int trkPt_F[1000];
  int trkEta_F[1000];
  int trkPhi_F[1000];
  //!possibly need to add other high quality cuts
  jetTree[5]->SetBranchAddress("nTrk", &nTrk_F);
  jetTree[5]->SetBranchAddress("trkPt", &trkPt_F);
  jetTree[5]->SetBranchAddress("trkEta", &trkEta_F);
  jetTree[5]->SetBranchAddress("trkPhi", &trkPhi_F);
  
  //! genparticles Tree 
  int mult_F;
  std::vector<float> *genp_pt_F = 0;
  std::vector<float> *genp_eta_F = 0;
  std::vector<float> *genp_phi_F = 0;
  std::vector<int> *genp_pdgid_F = 0;
  jetTree[6]->SetBranchAddress("mult", &mult_F);
  jetTree[6]->SetBranchAddress("pt", &genp_pt_F);
  jetTree[6]->SetBranchAddress("eta", &genp_eta_F);
  jetTree[6]->SetBranchAddress("phi", &genp_phi_F);
  jetTree[6]->SetBranchAddress("pdg", &genp_pdgid_F);

  //! runAnalyzer
  float xsec_F;
  jetTree[7]->SetBranchAddress("xsec", &xsec_F);

  //! run the analysis: 
  if(debugMode) std::cout<<"Running through all the events now"<<std::endl;
  Long64_t nentries = jetTree[0]->GetEntries();
  if(debugMode) nentries = 1000;
  
  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)std::cout<<nEvt<<"/"<<nentries<<std::endl;
    if(debugMode){
      std::cout<<"*******************************************"<<std::endl;
      std::cout<<"Start of Event: "<<"  nEvt = "<<nEvt<<std::endl;
    }
    for(int t = 0; t<_ntrees_MC; ++t)
      jetTree[t]->GetEntry(nEvt);

    double weight = xsec_F;

    //! get the particles from the ParticleFlow collection and make a jet collection from them
    vector<fastjet::PseudoJet> pfObjects;
    for(unsigned ip = 0; ip < pfPt_F->size(); ++ip){
      fastjet::PseudoJet particle(getPX(pfPt_F->at(ip), pfPhi_F->at(ip)),
				  getPY(pfPt_F->at(ip), pfPhi_F->at(ip)),
				  getPZ(pfPt_F->at(ip), pfEta_F->at(ip)),
				  pfEnergy_F->at(ip));
      pfObjects.push_back(particle);
    }

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, radius*0.1);
    fastjet::ClusterSequence cs(pfObjects, jet_def);
    vector<fastjet::PseudoJet> recoJets = fastjet::sorted_by_pt(cs.inclusive_jets(ptCut));

    for (unsigned ij = 0; ij < recoJets.size(); ++ij){
      fastjet::PseudoJet jet = recoJets.at(ij); 

      if(debugMode){
	std::cout<<"Ungroomed Jet "<<std::endl;
	std::cout<<"    raw jet pT = "<<jet.pt()<<", eta = "<<jet.eta()<<", phi = "<<jet.phi()<<std::endl;
      }
      
      //! get the JEC
      vector<JetCorrectorParameters> vpar;   
      FactorizedJetCorrector *JEC = new FactorizedJetCorrector(vpar);
      std::string L2Name;
      std::string L3Name;
      std::string L2L3Res; 
      // L2Name="Summer16_25nsV5_MC_L2Relative_AK4PF.txt";
      // L3Name="Summer16_25nsV5_MC_L3Absolute_AK4PF.txt";
      // L2L3Res="Summer16_25nsV5_MC_L2L3Residual_AK4PF.txt";
      L2Name="Spring16_25nsV8_MC_L2Relative_AK4PF.txt";
      L3Name="Spring16_25nsV8_MC_L3Absolute_AK4PF.txt";
      L2L3Res="Spring16_25nsV8_MC_L2L3Residual_AK4PF.txt";
      // cout << "Using .txt files to update JECs..." << endl;
      // cout << "L2: "<< L2Name << endl;
      // cout << "L3: "<< L3Name << endl;
      JetCorrectorParameters *parl2 = new JetCorrectorParameters(L2Name.c_str());
      JetCorrectorParameters *parl3 = new JetCorrectorParameters(L3Name.c_str());
      JetCorrectorParameters *parl2l3Res = new JetCorrectorParameters(L2L3Res.c_str());
      vpar.push_back(*parl2);
      vpar.push_back(*parl3);
      vpar.push_back(*parl2l3Res);
      JEC = new FactorizedJetCorrector(vpar);
      JEC->setJetEta(jet.eta());
      JEC->setJetPt(jet.pt());
      float jetcorr = JEC->getCorrection();
      float jetpt = jetcorr * jet.pt();

      double radMom1 = 0.0;
      double radMom2 = 0.0;
	
      for(size_t b = 0; b < _nbeta; ++b){
	double rm = 0.0;
	//! loop over the candidates of the jet and include them in the grid
	for ( unsigned ic = 0; ic<jet.constituents().size(); ++ic ){
	  fastjet::PseudoJet c = jet.constituents().at(ic);
	  double delR = deltaR(jet.eta(), jet.phi_std(), c.eta(), c.phi_std());
	  rm+=(c.pt()*pow(delR,_betaValues[b]))/jetpt;
	}
	if(b==0)
	  radMom1 = rm;
	else
	  radMom2 = rm;
      }
		
      
    }
    
    //! get the genparticles and make the genjet collections
    // vector<fastjet::PseudoJet> genParticles;    
    
  }//! end of event loop
  
  //! write to output file
  TFile *fout = new TFile(Form("%s", outfile.c_str()),"RECREATE");
  fout->cd();
  fout->Close();

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  return 0;
    
}


int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=9 && argc!=1){
    std::cout<<"for tests on default inputs, do..." <<std::endl;
    std::cout<<"./readForest_MC.exe";
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readForest_MC.exe ";
    std::cout<<"<inputFileList> <startFile> <endFile> ";
    std::cout<<"<jetAlgo> <jetRadius> <jetType> <debugMode> ";
    std::cout<<"<outputFilename> ";
    std::cout<<std::endl<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl;
    return rStatus;
  }
  
  rStatus=1;
  if(argc==1)
    rStatus = readForest_MC();
  else{
    std::string inputFileList=argv[1];
    int startfile= atoi(argv[2]);
    int endfile= atoi(argv[3]);  
    std::string jetAlgo=argv[4];
    int jetRadius= atoi(argv[5]);
    std::string jetType=argv[6];
    bool debug=(bool)atoi(argv[7]);
    std::string outputFileName=argv[8];
    
    rStatus = readForest_MC(inputFileList,
			    startfile,
			    endfile,
			    jetAlgo,
			    jetRadius,
			    jetType,
			    debug,
			    outputFileName);
  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}

