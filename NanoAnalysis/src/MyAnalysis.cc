#define MyAnalysis_cxx
#include "MyAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "PU_reWeighting.h"
#include "sumOfWeights.h"
#include "sumOfWeightsSignal.h"
#include "lepton_candidate.h"
#include "jet_candidate.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TDirectory.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <time.h>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>
#include "RoccoR.h"
#include "BTagCalibrationStandalone.h"
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
//#include "CondFormats/Serialization/interface/Archive.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "GEScaleSyst.h"
#include "Utils.h"
#include "correction.h"
#include "WCPoint.h"
#include "WCFit.h"
#include "TH1EFT.h"
#endif

using namespace correction;

float topPtPowheg(float pt){
  return (0.973 - (0.000134 * pt) + (0.103 * exp(pt * (-0.0118))));
}

float topPtMGLO(float x){
  return (0.688 -  0.0000174*x + 0.824*exp(-0.0000253*x)/(pow(x,0.2185)));
}

void MyAnalysis::Loop(TString fname, TString data, TString dataset ,string year, TString run, float xs, float lumi, float Nevent, int iseft, int nRuns)
//void MyAnalysis::Loop(TString fname, TString sname, TString data, TString dataset ,TString year, TString run, float xs, float lumi, float Nevent)
{

  TH1EFT  *crossSection = new TH1EFT("crossSection","crossSection",1,0,1);
  TH2F  btagEff_b_H;
  TH2F  btagEff_c_H;
  TH2F  btagEff_udsg_H;
  TH2F  sf_triggeree_H;
  TH2F  sf_triggeremu_H;
  TH2F  sf_triggermumu_H;
  std::string rochesterFile;
  std::string btagFile;
  GEScaleSyst *GE = new GEScaleSyst();

  RoccoR  rc;
  if(year == "2016preVFP")    rochesterFile = "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/RoccoR2016aUL.txt";
  if(year == "2016postVFP")    rochesterFile = "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/RoccoR2016bUL.txt";
  if(year == "2017")    rochesterFile = "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/RoccoR2017UL.txt";
  if(year == "2018")    rochesterFile = "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/RoccoR2018UL.txt";
  rc.init(rochesterFile);

  if(data == "mc"){
    TFile *f_btagEff_Map = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/btagEff.root");
    btagEff_b_H = *(TH2F*)f_btagEff_Map->Get((year + "_h2_BTaggingEff_b").c_str());
    btagEff_c_H = *(TH2F*)f_btagEff_Map->Get((year +"_h2_BTaggingEff_c").c_str());
    btagEff_udsg_H = *(TH2F*)f_btagEff_Map->Get((year +"_h2_BTaggingEff_udsg").c_str());

    TFile *f_trigger = new TFile(("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/TriggerSF_" + year + "_UL.root").c_str());
    sf_triggeree_H = *(TH2F*)f_trigger->Get("h2D_SF_ee_lepABpt_FullError");
    sf_triggeremu_H = *(TH2F*)f_trigger->Get("h2D_SF_emu_lepABpt_FullError");
    sf_triggermumu_H = *(TH2F*)f_trigger->Get("h2D_SF_mumu_lepABpt_FullError");

    f_btagEff_Map->Close();
    f_trigger->Close();
  }

   string eleSF="";
   string muSF="";
   string bSF="";

   if(year == "2016preVFP"){
     eleSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/EGM/2017_UL/electron.json";
     muSF= "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/MUO/2017_UL/muon_Z.json";
     bSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/BTV/2017_UL/btagging.json";
   }

   if(year == "2016postVFP"){
     eleSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/EGM/2017_UL/electron.json";
     muSF= "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/MUO/2017_UL/muon_Z.json";
     bSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/BTV/2017_UL/btagging.json";
   }

   if(year == "2017"){
     eleSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/EGM/2017_UL/electron.json";
     muSF= "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/MUO/2017_UL/muon_Z.json";
     bSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/BTV/2017_UL/btagging.json";
   }

   if(year == "2018"){
     eleSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/EGM/2017_UL/electron.json";
     muSF= "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/MUO/2017_UL/muon_Z.json";
     bSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/BTV/2017_UL/btagging.json";
   }

   auto csetFileEleSF = CorrectionSet::from_file(eleSF);
   auto csetEleIdReco = csetFileEleSF->at("UL-Electron-ID-SF");

   auto csetFileMuSF = CorrectionSet::from_file(muSF);
   auto csetMuHighPtId = csetFileMuSF->at("NUM_HighPtID_DEN_TrackerMuons");
   auto csetMuHighPtIso = csetFileMuSF->at("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut");

   auto csetFilebSF = CorrectionSet::from_file(bSF);
   auto csetLightJetSF = csetFilebSF->at("deepCSV_incl");
   auto csetBcJetSF = csetFilebSF->at("deepCSV_mujets");

   TRandom3 Tr;

   TMVA::Tools::Instance();
   TMVA::Reader *readerMVA = new TMVA::Reader( "!Color:!Silent" );
   Float_t leading_pt, jet_leading_pt, deltaR_ll, MET, n_jet;

   readerMVA->AddVariable ("leading_pt", &leading_pt);
   readerMVA->AddVariable ("jet_leading_pt", &jet_leading_pt);
   readerMVA->AddVariable ("deltaR_ll", &deltaR_ll);
   readerMVA->AddVariable ("MET", &MET);
   readerMVA->AddVariable ("n_jet", &n_jet);
   readerMVA->BookMVA( "BDT_1b_tc_BDT", "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/TMVA_BDT_1b_BDT.weights.xml");

   Double_t ptBins[11] = {30., 40., 60., 80., 100., 150., 200., 300., 400., 500., 1000.};
   Double_t etaBins [4]= {0., 0.6, 1.2, 2.4};
   TH2D *h2_BTaggingEff_Denom_b    = new TH2D("h2_BTaggingEff_Denom_b"   , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Denom_c    = new TH2D("h2_BTaggingEff_Denom_c"   , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Denom_udsg = new TH2D("h2_BTaggingEff_Denom_udsg", ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Num_b      = new TH2D("h2_BTaggingEff_Num_b"     , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Num_c      = new TH2D("h2_BTaggingEff_Num_c"     , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Num_udsg   = new TH2D("h2_BTaggingEff_Num_udsg"  , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);


  typedef vector<std::shared_ptr<TH1EFT>> Dim1;
  typedef vector<Dim1> Dim2;
  typedef vector<Dim2> Dim3;
  typedef vector<Dim3> Dim4;

  std::vector<TString> regions{"ll","llOffZ","llB1", "llBg1", };
  std::vector<TString> channels{"ee", "emu", "mumu"};
  std::vector<TString> vars   {"lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","llM","llPt","llDr","llDphi","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx", "llMZw","BDT","muPt","Electron_pt[l]","muEta","eleEta"};
  std::vector<int>    nbins   {60      ,20       ,25       ,25      ,20       ,25       ,30   ,20    ,25    ,15      ,20      ,20       ,25       ,10    ,6      ,30   ,20      ,70    ,80      ,100, 70,70,20,20};
  std::vector<float> lowEdge  {0       ,-3       ,-4       ,0       ,-3       ,-4       ,0    ,0     ,0     ,0       ,0       ,-3       ,-4       ,0     ,0      ,0    ,-4      ,0     ,70      ,-0.4, 0,0,-3,-3};
  std::vector<float> highEdge {1500     ,3        ,4        ,1000     ,3        ,4        ,500  ,200   ,7     ,4       ,300     ,3        ,4        ,10    ,6      ,210  ,4       ,70    ,110,     0.6 ,1500,1500,3,3};


  Dim3 Hists(channels.size(),Dim2(regions.size(),Dim1(vars.size())));
  std::stringstream name;
  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int l=0;l<vars.size();++l){
        name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l];
        std::shared_ptr<TH1EFT> h_test(new TH1EFT((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]));
        h_test->StatOverflows(kTRUE);
        h_test->Sumw2(kTRUE);
        Hists[i][k][l] = h_test;
        name.str("");
      }
    }
  }

   std::vector<string> wc_names_lst={};
   std::vector<string> wc_names_lst_BNV={"cT", "cS"};
   std::vector<string> wc_names_lst_FCNC={"ctlS", "cte", "ctl", "ctlT", "ctZ", "cpt", "cpQM", "ctA", "cQe", "ctG", "cQlM"};
   if (fname.Contains("BNV")) wc_names_lst = wc_names_lst_BNV;
   if (fname.Contains("FCNC")) wc_names_lst = wc_names_lst_FCNC;

//  std::vector<TString> sys{"eleRecoSf", "eleIDSf", "muIdSf", "muIsoSf", "bcTagSF", "udsgTagSF","pu", "prefiring", "trigSF","jes", "jer","unclusMET","muonScale","electronScale","muonRes" };
//  Dim4 HistsSysUp(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(sys.size()))));
//  Dim4 HistsSysDown(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(sys.size()))));
//
//  for (int i=0;i<channels.size();++i){
//    for (int k=0;k<regions.size();++k){
//      for (int l=0;l<vars.size();++l){
//        for (int n=0;n<sys.size();++n){
//          name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l]<<"_"<<sys[n]<<"_Up";
//          std::shared_ptr<TH1F> h_test(new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]));
//          h_test->StatOverflows(kTRUE);
//          h_test->Sumw2(kTRUE);
//          HistsSysUp[i][k][l][n] = h_test;
//          name.str("");
//          name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l]<<"_"<<sys[n]<<"_Down";
//          std::shared_ptr<TH1F> h_test2(new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]));
//          h_test2->StatOverflows(kTRUE);
//          h_test2->Sumw2(kTRUE);
//          HistsSysDown[i][k][l][n] = h_test2;
//          name.str("");
//        }
//      }
//    }
//  }

//  Dim4 HistsJECUp(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(sysJecNames.size()))));
//  Dim4 HistsJECDown(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(sysJecNames.size()))));
//  for (int i=0;i<channels.size();++i){
//    for (int k=0;k<regions.size();++k){
//      for (int l=0;l<vars.size();++l){
//        for (int n=0;n<sysJecNames.size();++n){
//          name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l]<<"_"<<sysJecNames[n]<<"_Up";
//          std::shared_ptr<TH1F> h_test(new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]));
//          h_test->StatOverflows(kTRUE);
//          h_test->Sumw2(kTRUE);
//          HistsJECUp[i][k][l][n] = h_test;
//          name.str("");
//          name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l]<<"_"<<sysJecNames[n]<<"_Down";
//          std::shared_ptr<TH1F> h_test2(new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]));
//          h_test2->StatOverflows(kTRUE);
//          h_test2->Sumw2(kTRUE);
//          HistsJECDown[i][k][l][n] = h_test2;
//          name.str("");
//        }
//      }
//    }
//  }


//  int reweightSizeQscalePDF = 150;
//  int reweightSizePS = 14;
//  Dim4 HistsSysReweightsQscalePDF(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(reweightSizeQscalePDF))));
//  Dim4 HistsSysReweightsPS(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(reweightSizePS))));
//  sumOfWeights SW;
//  sumOfWeightsSignal SWS;
//  std::vector<float> SLW;
//  std::vector<float> SGW;
//
//  if (fname.Contains("TTTo2L2Nu")){
//     SLW = SW.LHEWeight(year);
//     SGW = SW.GenWeight(year);
//  }
//  if (fname.Contains("LFV")){
//     SLW = SWS.LHEWeightSignal(sname);
//     SGW = SWS.GenWeightSignal(sname);
//  }
//
//  if (fname.Contains("TTTo2L2Nu") || fname.Contains("LFV")){
//    for (int i=0;i<channels.size();++i){
//      for (int k=0;k<regions.size();++k){
//        for (int l=0;l<vars.size();++l){
//          for (int n=0;n<reweightSizeQscalePDF;++n){
//            name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l]<<"_QscalePDF_"<<n;
//            std::shared_ptr<TH1F> h_test(new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]));
//            h_test->StatOverflows(kTRUE);
//            h_test->Sumw2(kTRUE);
//            HistsSysReweightsQscalePDF[i][k][l][n] = h_test;
//            name.str("");
//          }
//          for (int n=0;n<reweightSizePS;++n){
//            name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l]<<"_PS_"<<n;
//            std::shared_ptr<TH1F> h_test(new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]));
//            h_test->StatOverflows(kTRUE);
//            h_test->Sumw2(kTRUE);
//            HistsSysReweightsPS[i][k][l][n] = h_test;
//            name.str("");
//          }
//        }
//      }
//    }
//  }

  TFile file_out ("ANoutput.root","RECREATE");
  TTree tree_out("analysis","main analysis") ;
  std::vector<lepton_candidate*> *selectedLeptons;
  std::vector<lepton_candidate*> *selectedLeptonsMuScaleUp;
  std::vector<lepton_candidate*> *selectedLeptonsMuScaleDown;
  std::vector<lepton_candidate*> *selectedLeptonsEleScaleUp;
  std::vector<lepton_candidate*> *selectedLeptonsEleScaleDown;
  std::vector<lepton_candidate*> *selectedLeptonsMuResUp;
  std::vector<lepton_candidate*> *selectedLeptonsMuResDown;
  std::vector<jet_candidate*> *selectedJets;
  std::vector<jet_candidate*> *selectedJetsJerUp;
  std::vector<jet_candidate*> *selectedJetsJerDown;
  std::vector<jet_candidate*> *selectedJetsJesUp;
  std::vector<jet_candidate*> *selectedJetsJesDown;
  std::vector<jet_candidate*> *JECJetsUp;
  std::vector<jet_candidate*> *JECJetsDown;
  std::vector<std::vector<jet_candidate*>> *JECsysUp;
  std::vector<std::vector<jet_candidate*>> *JECsysDown;
  std::vector<int> *JECsysNbtagUp;
  std::vector<int> *JECsysNbtagDown;
  std::vector<float> *JECsysMVAUp;
  std::vector<float> *JECsysMVADown;
  std::vector<float> *JECsysMETUp;
  std::vector<float> *JECsysMETDown;
  WCFit *eft_fit;


  TLorentzVector wp, wm, b, ab, top, atop;
  std::vector<float> nominalWeights;
//  nominalWeights.assign(sys.size(), 1);
//  std::vector<float> sysUpWeights;
//  sysUpWeights.assign(sys.size(), 1);
//  std::vector<float> sysDownWeights;
//  sysDownWeights.assign(sys.size(), 1);
  bool triggerPassEE;
  bool triggerPassEMu;
  bool triggerPassMuMu;
  bool metFilterPass;
  bool ifTopPt=false;
  int ch;
  float sf_Ele_Reco;
  float sf_Ele_ID;
  float sf_Mu_ID;
  float sf_Mu_ISO;
  float sf_Trigger;
  float weight_PU;
  float weight_Lumi;
  float weight_lep;
  float weight_lepB;
  float weight_prefiring;
  float weight_topPtPowheg;
  float weight_topPtMGLO;
  float MVAoutputJerUp;
  float MVAoutputJerDown;
  float MVAoutputJesUp;
  float MVAoutputJesDown;
  float MVAoutputMuScaleUp;
  float MVAoutputMuScaleDown;
  float MVAoutputEleScaleUp;
  float MVAoutputEleScaleDown;
  float MVAoutputMuResUp;
  float MVAoutputMuResDown;
  float MVAoutputUnclusMETUp;
  float MVAoutputUnclusMETDown;
  float metUnclusMETUp;
  float metUnclusMETDown;
  float MVAoutput;
  double P_bjet_data;
  double P_bjet_mc;
  int nAccept=0;
  int nbjetGen;
  int nbjet;
  int nbjetJesUp;
  int nbjetJesDown;
  int nbjetJerUp;
  int nbjetJerDown;
  float JECMETUpx;
  float JECMETUpy;
  float JECMETDownx;
  float JECMETDowny;
  float pt_res;
  double muPtSFRochester;
  int R;

  if (fname.Contains("TTTo2L2Nu") || fname.Contains("TTsys") || fname.Contains("LFVTt")) ifTopPt=true;

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  Long64_t ntr = fChain->GetEntries ();

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    displayProgress(jentry, ntr) ;

    triggerPassEE = false;
    triggerPassEMu = false;
    triggerPassMuMu = false;
    metFilterPass = false;
    ch =10;
    sf_Ele_Reco =1;
    sf_Ele_ID =1;
    sf_Mu_ID =1;
    sf_Mu_ISO =1;
    sf_Trigger =1;
    weight_PU =1;
    weight_Lumi =1;
    weight_lep =1;
    weight_lepB =1;
    weight_prefiring =1;
    weight_topPtPowheg =1;
    weight_topPtMGLO =1;
    P_bjet_data =1;
    P_bjet_mc =1;
    MVAoutput=0;
    nbjetGen=0;
    nbjet=0;
    nbjetJerUp=0;
    nbjetJerDown=0;
    nbjetJesUp=0;
    nbjetJesDown=0;
    metUnclusMETUp=0;
    metUnclusMETDown=0;
    

//    for (int n=0;n<sys.size();++n){
//      nominalWeights[n] =1;
//      sysUpWeights[n] =1;
//      sysDownWeights[n] =1;
//    }
//MET filters

    if (iseft) {
      eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, 1.0/nRuns);
//      for (UInt_t i=0;i<nWCnames;++i){
//         char ch[4];
////         for(int j = 0; j < 4; j++) ch[j] = (WCnames[i] >> (4-1-j)*8) & 0xFF;
//         cout<< " - "<<WCnames[i] <<endl;
////         for(int n=0 ; n<4 ; ++n) cout << ch[n];
////         cout<<" : "<<endl;;
//      }
      
    }
    else {
      eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0); 
    }
    crossSection->Fill(0.5, *eft_fit);

    delete eft_fit; 
//You should add Flag_BadPFMuonDzFilter to the MET filter list but since it is not available in v8, lets remove it for now.
   if(year == "2017" || year == "2018"){
     if ( Flag_goodVertices  &&  Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter &&  Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_eeBadScFilter && Flag_ecalBadCalibFilter) metFilterPass = true;
   }
   else{
     if ( Flag_goodVertices  &&  Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter &&  Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_eeBadScFilter) metFilterPass = true;
   }

//trigger MC
      if(data == "mc" && year == "2016"){
        if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele27_WPTight_Gsf ) triggerPassEE =true;
        if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele27_WPTight_Gsf || HLT_IsoMu24 || HLT_IsoTkMu24) triggerPassEMu =true;
        if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL || HLT_IsoMu24 || HLT_IsoTkMu24) triggerPassMuMu =true;
      }

      if(data == "mc" && year == "2017"){
        if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele35_WPTight_Gsf) triggerPassEE =true;
        if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele35_WPTight_Gsf || HLT_IsoMu27) triggerPassEMu =true;
        if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_IsoMu27) triggerPassMuMu =true;
      }

      if(data == "mc" && year == "2018"){
        if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele32_WPTight_Gsf) triggerPassEE =true;
        if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele32_WPTight_Gsf || HLT_IsoMu24) triggerPassEMu =true;
        if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_IsoMu24) triggerPassMuMu =true;
      }

//trigger DATA
    if(data == "data"){
      if(year == "2016"){
        if(run == "H"){
          if(dataset=="MuonEG"){
            if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) triggerPassEMu =true;
          }
          if(dataset=="SingleElectron"){
            if(!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && HLT_Ele27_WPTight_Gsf) triggerPassEE =true;
            if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && HLT_Ele27_WPTight_Gsf) triggerPassEMu =true;
          }
          if(dataset=="SingleMuon"){
            if(!(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ) && (HLT_IsoMu24 || HLT_IsoTkMu24)) triggerPassMuMu =true;
            if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele27_WPTight_Gsf) && (HLT_IsoMu24 || HLT_IsoTkMu24)) triggerPassEMu =true;
          }
          if(dataset=="DoubleEG"){
            if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) triggerPassEE =true;
          }
          if(dataset=="DoubleMuon"){
            if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ) triggerPassMuMu =true;
          }
        }
        if(run != "H"){
          if(dataset=="MuonEG"){
            if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL) triggerPassEMu =true;
          }
          if(dataset=="SingleElectron"){
            if(!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && HLT_Ele27_WPTight_Gsf) triggerPassEE =true;
            if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL) && HLT_Ele27_WPTight_Gsf) triggerPassEMu =true;
          }
          if(dataset=="SingleMuon"){
            if(!(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL) && (HLT_IsoMu24 || HLT_IsoTkMu24)) triggerPassMuMu =true;
            if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele27_WPTight_Gsf) && (HLT_IsoMu24 || HLT_IsoTkMu24)) triggerPassEMu =true;
          }
          if(dataset=="DoubleEG"){
            if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) triggerPassEE =true;
          }
          if(dataset=="DoubleMuon"){
            if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ) triggerPassMuMu =true;
          }
        }
      }
      if(year == "2017"){
        if(dataset=="MuonEG"){
          if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL ) triggerPassEMu =true;
        }
        if(dataset=="SingleElectron"){
          if(!HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL && HLT_Ele35_WPTight_Gsf) triggerPassEE =true;
          if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL) && HLT_Ele35_WPTight_Gsf) triggerPassEMu =true;
        }
        if(dataset=="SingleMuon"){
           if(!HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 && HLT_IsoMu27) triggerPassMuMu =true;
          if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele35_WPTight_Gsf) && HLT_IsoMu27) triggerPassEMu =true;
        }
        if(dataset=="DoubleEG"){
          if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL) triggerPassEE =true;
        }
        if(dataset=="DoubleMuon"){
          if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8)triggerPassMuMu =true;
        }
      }
      if(year == "2018"){
        if(dataset=="MuonEG"){
          if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL) triggerPassEMu =true;
        }
        if(dataset=="EGamma"){
          if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele32_WPTight_Gsf) triggerPassEE =true;
          if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL) && HLT_Ele32_WPTight_Gsf) triggerPassEMu =true;
      }
      if(dataset=="SingleMuon"){
        if(!HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 && HLT_IsoMu24) triggerPassMuMu =true;
        if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele32_WPTight_Gsf) && HLT_IsoMu24) triggerPassEMu =true;
      }
      if(dataset=="DoubleMuon"){
        if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8) triggerPassMuMu =true;
      }
    }
  }

    if(!(triggerPassEE || triggerPassEMu || triggerPassMuMu)) continue;
    if(!metFilterPass) continue;

  selectedLeptons = new std::vector<lepton_candidate*>();
  selectedLeptonsMuScaleUp = new std::vector<lepton_candidate*>();
  selectedLeptonsMuScaleDown = new std::vector<lepton_candidate*>();
  selectedLeptonsEleScaleUp = new std::vector<lepton_candidate*>();
  selectedLeptonsEleScaleDown = new std::vector<lepton_candidate*>();
  selectedLeptonsMuResUp = new std::vector<lepton_candidate*>();
  selectedLeptonsMuResDown = new std::vector<lepton_candidate*>();

// electron
    for (int l=0;l<nElectron;l++){
      if(abs(Electron_eta[l]) > 2.4 || (abs(Electron_eta[l])> 1.4442 && (abs(Electron_eta[l])< 1.566))) continue;
      if(Electron_cutBased[l] < 4) continue;
      if (data == "mc"){
        if((Electron_pt[l] + 0.004*Electron_pt[l]) > 20 && abs(Electron_eta[l])< 1.4442) selectedLeptonsEleScaleUp->push_back(new lepton_candidate(Electron_pt[l]+ 0.004*Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
        if((Electron_pt[l] + 0.008*Electron_pt[l]) > 20 && abs(Electron_eta[l])> 1.4442) selectedLeptonsEleScaleUp->push_back(new lepton_candidate(Electron_pt[l]+ 0.008*Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
        if((Electron_pt[l] - 0.004*Electron_pt[l]) > 20 && abs(Electron_eta[l])< 1.4442) selectedLeptonsEleScaleDown->push_back(new lepton_candidate(Electron_pt[l]- 0.004*Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
        if((Electron_pt[l] - 0.008*Electron_pt[l]) > 20 && abs(Electron_eta[l])> 1.4442) selectedLeptonsEleScaleDown->push_back(new lepton_candidate(Electron_pt[l]- 0.008*Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
      }
      if(Electron_pt[l] <20) continue;
      selectedLeptons->push_back(new lepton_candidate(Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
      if (data == "mc"){
        selectedLeptonsMuScaleUp->push_back(new lepton_candidate(Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
        selectedLeptonsMuScaleDown->push_back(new lepton_candidate(Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
        selectedLeptonsMuResUp->push_back(new lepton_candidate(Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
        selectedLeptonsMuResDown->push_back(new lepton_candidate(Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
        sf_Ele_Reco = sf_Ele_Reco * csetEleIdReco->evaluate({year, "sf", "RecoAbove20", Electron_eta[l],Electron_pt[l]}); 
//        nominalWeights[0] = nominalWeights[0] * cset_2017->evaluate({"2017", "sf", "RecoAbove20", Electron_eta[l],Electron_pt[l]});
//        sysUpWeights[0] = sysUpWeights[0] * cset_2017->evaluate({"2017", "sfup", "RecoAbove20", Electron_eta[l],Electron_pt[l]});
//        sysDownWeights[0] = sysDownWeights[0] * cset_2017->evaluate({"2017", "sfdown", "RecoAbove20", Electron_eta[l],Electron_pt[l]});
//
        sf_Ele_ID = sf_Ele_ID * csetEleIdReco->evaluate({year, "sf", "Tight", Electron_eta[l],Electron_pt[l]});
//        nominalWeights[1] = nominalWeights[1] * cset_2017->evaluate({"2017", "sf", "Tight", Electron_eta[l],Electron_pt[l]});
//        sysUpWeights[1] = sysUpWeights[1] * cset_2017->evaluate({"2017", "sfup", "Tight", Electron_eta[l],Electron_pt[l]});
//        sysDownWeights[1] = sysDownWeights[1] * cset_2017->evaluate({"2017", "sfdown", "Tight", Electron_eta[l],Electron_pt[l]});
      }
    }

// Muon selection
    for (int l=0;l<nMuon;l++){
      if(abs(Muon_eta[l]) > 2.4) continue;
      if(Muon_highPtId[l]<2) continue;
      if(Muon_tkIsoId[l]<1) continue;
      if(Muon_tunepRelPt[l]==1){
//      if(Muon_pt[l] == (*mu_it_pt)[l]){
        if(data == "data") {
          muPtSFRochester = rc.kScaleDT(Muon_charge[l], Muon_pt[l],Muon_eta[l],Muon_phi[l], 0, 0);
          if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
        }
        if (data == "mc"){
          if (Muon_genPartIdx[l]>=0) muPtSFRochester = rc.kSpreadMC(Muon_charge[l], Muon_pt[l],Muon_eta[l],Muon_phi[l], GenPart_pt[Muon_genPartIdx[l]],0, 0);
          if (Muon_genPartIdx[l]<0) muPtSFRochester = rc.kSmearMC(Muon_charge[l], Muon_pt[l],Muon_eta[l],Muon_phi[l], Muon_nTrackerLayers[l] , gRandom->Rndm(),0, 0);
          if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
//          if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptonsEleScaleUp->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
//          if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptonsEleScaleDown->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
//          if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptonsMuResUp->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
//          if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptonsMuResDown->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
//          if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptonsMuScaleUp->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
//          if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptonsMuScaleDown->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
        }
      }
      if(Muon_tunepRelPt[l]!=1){
//      if(Muon_pt[l] != (*mu_it_pt)[l]){
        pt_res = 0.02;
        if(abs(Muon_eta[l]) < 1.2) pt_res = 0.01;
        if (data == "data" && Muon_tunepRelPt[l]*Muon_pt[l]>20) selectedLeptons->push_back(new lepton_candidate(Muon_tunepRelPt[l]*Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
        if(data == "mc"){
          if (Muon_tunepRelPt[l]*Muon_pt[l]>20) selectedLeptons->push_back(new lepton_candidate(Muon_tunepRelPt[l]*Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
//          if (GE->GEScaleCorrPt(1600, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]) > 20) selectedLeptons->push_back(new lepton_candidate(GE->GEScaleCorrPt(1600, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]),Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
//          if (GE->GEScaleCorrPt(1600, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]) > 20) selectedLeptonsEleScaleUp->push_back(new lepton_candidate(GE->GEScaleCorrPt(1600, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]),Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
//          if (GE->GEScaleCorrPt(1600, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]) > 20) selectedLeptonsEleScaleDown->push_back(new lepton_candidate(GE->GEScaleCorrPt(1600, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]),Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
//
//          if (GE->GEScaleCorrPt(1600, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l])*(Tr.Gaus(1, pt_res)) > 20) selectedLeptonsMuResUp->push_back(new lepton_candidate((Tr.Gaus(1, pt_res))*GE->GEScaleCorrPt(1600, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]),Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
//          if (GE->GEScaleCorrPt(1600, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l])*(Tr.Gaus(1, pt_res)) > 20) selectedLeptonsMuResDown->push_back(new lepton_candidate((Tr.Gaus(1, pt_res))*GE->GEScaleCorrPt(1600, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]),Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
//
//          if (GE->GEScaleCorrPt(1601, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]) > 20) selectedLeptonsMuScaleUp->push_back(new lepton_candidate(GE->GEScaleCorrPt(1601, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]),Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
//          if (GE->GEScaleCorrPt(1602, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]) > 20) selectedLeptonsMuScaleDown->push_back(new lepton_candidate(GE->GEScaleCorrPt(1602, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]),Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
        }
      }
      if(Muon_pt[l] <20) continue;
      if (data == "mc") {
        sf_Mu_ID = sf_Mu_ID * csetMuHighPtId->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_pt[l], "sf"});
//        nominalWeights[2] = nominalWeights[2] * csetMuHighPtId->evaluate({"2017", abs(Muon_eta[l]),  Muon_pt[l], "sf"});
//        sysUpWeights[2] = sysUpWeights[2] * csetMuHighPtId->evaluate({"2017", abs(Muon_eta[l]),  Muon_pt[l], "systup"});
//        sysDownWeights[2] = sysDownWeights[2] * csetMuHighPtId->evaluate({"2017", abs(Muon_eta[l]),  Muon_pt[l], "systdown"});

        sf_Mu_ISO = sf_Mu_ISO * csetMuHighPtIso->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_pt[l], "sf"});
//        nominalWeights[3] = nominalWeights[3] * csetMuHighPtIso->evaluate({"2017", abs(Muon_eta[l]),  Muon_pt[l], "sf"});
//        sysUpWeights[3] = sysUpWeights[3] * csetMuHighPtIso->evaluate({"2017", abs(Muon_eta[l]),  Muon_pt[l], "systup"});
//        sysDownWeights[3] = sysDownWeights[3] * csetMuHighPtIso->evaluate({"2017", abs(Muon_eta[l]),  Muon_pt[l], "systdown"});
      }
    }
    sort(selectedLeptons->begin(), selectedLeptons->end(), ComparePtLep);
//    sort(selectedLeptonsMuScaleUp->begin(), selectedLeptonsMuScaleUp->end(), ComparePtLep);
//    sort(selectedLeptonsMuScaleDown->begin(), selectedLeptonsMuScaleDown->end(), ComparePtLep);
//    sort(selectedLeptonsMuResDown->begin(), selectedLeptonsMuResDown->end(), ComparePtLep);
//    sort(selectedLeptonsMuResUp->begin(), selectedLeptonsMuResUp->end(), ComparePtLep);
//    sort(selectedLeptonsEleScaleUp->begin(), selectedLeptonsEleScaleUp->end(), ComparePtLep);
//    sort(selectedLeptonsEleScaleDown->begin(), selectedLeptonsEleScaleDown->end(), ComparePtLep);
    if(selectedLeptons->size()!=2 ||
      ((*selectedLeptons)[0]->pt_ <25) ||
      ((*selectedLeptons)[0]->charge_ * (*selectedLeptons)[1]->charge_ == 1) ||
      ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()<20) {
      for (int l=0;l<selectedLeptons->size();l++){
        delete (*selectedLeptons)[l];
      }
      for (int l=0;l<selectedLeptonsMuScaleUp->size();l++){
        delete (*selectedLeptonsMuScaleUp)[l];
      }
      for (int l=0;l<selectedLeptonsMuScaleDown->size();l++){
        delete (*selectedLeptonsMuScaleDown)[l];
      }
      for (int l=0;l<selectedLeptonsEleScaleUp->size();l++){
        delete (*selectedLeptonsEleScaleUp)[l];
      }
      for (int l=0;l<selectedLeptonsEleScaleDown->size();l++){
        delete (*selectedLeptonsEleScaleDown)[l];
      }
      for (int l=0;l<selectedLeptonsMuResUp->size();l++){
        delete (*selectedLeptonsMuResUp)[l];
      }
      for (int l=0;l<selectedLeptonsMuResDown->size();l++){
        delete (*selectedLeptonsMuResDown)[l];
      }
      selectedLeptons->clear();
      selectedLeptons->shrink_to_fit();
      delete selectedLeptons;
      selectedLeptonsMuScaleUp->clear();
      selectedLeptonsMuScaleUp->shrink_to_fit();
      delete selectedLeptonsMuScaleUp;
      selectedLeptonsMuScaleDown->clear();
      selectedLeptonsMuScaleDown->shrink_to_fit();
      delete selectedLeptonsMuScaleDown;
      selectedLeptonsMuResUp->clear();
      selectedLeptonsMuResUp->shrink_to_fit();
      delete selectedLeptonsMuResUp;
      selectedLeptonsMuResDown->clear();
      selectedLeptonsMuResDown->shrink_to_fit();
      delete selectedLeptonsMuResDown;
      selectedLeptonsEleScaleUp->clear();
      selectedLeptonsEleScaleUp->shrink_to_fit();
      delete selectedLeptonsEleScaleUp;
      selectedLeptonsEleScaleDown->clear();
      selectedLeptonsEleScaleDown->shrink_to_fit();
      delete selectedLeptonsEleScaleDown;
      continue;
    }
//remove HEM effected region from 2018 samples
    if (year == "2018" && (*selectedLeptons)[0]->lep_ == 1 && (*selectedLeptons)[0]->eta_ < -1.4 && (*selectedLeptons)[0]->phi_< -0.8 && (*selectedLeptons)[0]->phi_ > -1.6) continue;
    if (year == "2018" && (*selectedLeptons)[1]->lep_ == 1 && (*selectedLeptons)[1]->eta_ < -1.4 && (*selectedLeptons)[1]->phi_< -0.8 && (*selectedLeptons)[1]->phi_ > -1.6) continue;

//categorize dilepton channels
    if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 2) ch = 0;
    if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 11) ch = 1;
    if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 20) ch = 2;
    if(ch ==0 && !triggerPassEE) continue;
    if(ch ==1 && !triggerPassEMu) continue;
    if(ch ==2 && !triggerPassMuMu) continue;
//jets

    selectedJets = new std::vector<jet_candidate*>();
    selectedJetsJerUp = new std::vector<jet_candidate*>();
    selectedJetsJerDown = new std::vector<jet_candidate*>();
    selectedJetsJesUp = new std::vector<jet_candidate*>();
    selectedJetsJesDown = new std::vector<jet_candidate*>();
    bool jetlepfail;
    for (int l=0;l<nJet;l++){
      if(Jet_jetId[l]<6) continue;
      jetlepfail = false;
      for (int i=0;i<selectedLeptons->size();i++){
        if(deltaR((*selectedLeptons)[i]->eta_,(*selectedLeptons)[i]->phi_,Jet_eta[l],Jet_phi[l]) < 0.4 ) jetlepfail=true;
      }
      if(jetlepfail) continue;
      if(data == "mc" && abs(Jet_eta[l]) < 2.4){
        if(Jet_pt[l] >30) selectedJets->push_back(new jet_candidate(Jet_pt[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepB[l], year,Jet_partonFlavour[l]));
//        if((*jet_SmearedJetResUp_pt)[l] >30) selectedJetsJerUp->push_back(new jet_candidate((*jet_SmearedJetResUp_pt)[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepB[l], year,Jet_partonFlavour[l]));
//        if((*jet_SmearedJetResDown_pt)[l] >30) selectedJetsJerDown->push_back(new jet_candidate((*jet_SmearedJetResDown_pt)[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepB[l], year,Jet_partonFlavour[l]));
//        if((*jet_SmearedJetEnUp_pt)[l] >30) selectedJetsJesUp->push_back(new jet_candidate((*jet_SmearedJetEnUp_pt)[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepB[l], year,Jet_partonFlavour[l]));
//        if((*jet_SmearedJetEnDown_pt)[l] >30) selectedJetsJesDown->push_back(new jet_candidate((*jet_SmearedJetEnDown_pt)[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepB[l], year,Jet_partonFlavour[l]));
      }
      if(data == "data" && Jet_pt[l] >30 && abs(Jet_eta[l]) < 2.4){
        selectedJets->push_back(new jet_candidate(Jet_pt[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepB[l],year,0));
//        selectedJetsJerUp->push_back(new jet_candidate(Jet_pt[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepB[l],year,0));
//        selectedJetsJerDown->push_back(new jet_candidate(Jet_pt[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepB[l],year,0));
//        selectedJetsJesUp->push_back(new jet_candidate(Jet_pt[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepB[l],year,0));
//        selectedJetsJesDown->push_back(new jet_candidate(Jet_pt[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepB[l],year,0));
      }
    }


    sort(selectedJets->begin(), selectedJets->end(), ComparePtJet);
    sort(selectedJetsJerUp->begin(), selectedJetsJerUp->end(), ComparePtJet);
    sort(selectedJetsJerDown->begin(), selectedJetsJerDown->end(), ComparePtJet);
    sort(selectedJetsJesUp->begin(), selectedJetsJesUp->end(), ComparePtJet);
    sort(selectedJetsJesDown->begin(), selectedJetsJesDown->end(), ComparePtJet);

    for (int l=0;l<selectedJetsJerUp->size();l++){
      if((*selectedJetsJerUp)[l]->btag_) nbjetJerUp++;
    }
    for (int l=0;l<selectedJetsJerDown->size();l++){
      if((*selectedJetsJerDown)[l]->btag_) nbjetJerDown++;
    }
    for (int l=0;l<selectedJetsJesUp->size();l++){
      if((*selectedJetsJesUp)[l]->btag_) nbjetJesUp++;
    }
    for (int l=0;l<selectedJetsJesDown->size();l++){
      if((*selectedJetsJesDown)[l]->btag_) nbjetJesDown++;
    }


//    JECsysUp = new std::vector<std::vector<jet_candidate*>>();
//    JECsysDown = new std::vector<std::vector<jet_candidate*>>();
//    JECsysMETUp = new std::vector<float>();
//    JECsysMETDown = new std::vector<float>();
//    JECsysMVAUp = new std::vector<float>();
//    JECsysMVADown = new std::vector<float>();
//
//    for (int n=0;n<sysJecNames.size();++n){
//      JECJetsUp= new std::vector<jet_candidate*>();
//      JECJetsDown= new std::vector<jet_candidate*>();
//      double sup = 0;
//      double sdw = 0;
//      JECMETUpx =   0;
//      JECMETUpy =   0;
//      JECMETDownx = 0;
//      JECMETDowny = 0;
//      for (int l=0;l<jet_pt->size();l++){
//        if(data == "data" || abs(Jet_eta[l]) > 2.4) continue;
//        if(year == "2016" && !(*jet_isJetIDTightLepVeto_2016)[l]) continue;
//        if(year == "2017" && !(*jet_isJetIDLepVeto_2017)[l]) continue;
//        if(year == "2018" && !(*jet_isJetIDLepVeto_2018)[l]) continue;
//        jetlepfail = false;
//        for (int i=0;i<selectedLeptons->size();i++){
//          if(deltaR((*selectedLeptons)[i]->eta_,(*selectedLeptons)[i]->phi_,Jet_eta[l],Jet_phi[l]) < 0.4 ) jetlepfail=true;
//        }
//        if(jetlepfail) continue;
//        JetCorrectionUncertainty *unc = vsrc[n];
//        unc->setJetPt(Jet_pt[l]);
//        unc->setJetEta(Jet_eta[l]);
//        sup = unc->getUncertainty(true);
//        if ((1+sup)*Jet_pt[l]>30) {
//          JECJetsUp->push_back(new jet_candidate((1+sup)*Jet_pt[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepB[l], year,Jet_partonFlavour[l]));
//          JECMETUpx = JECMETUpx + (sup * Jet_pt[l] * cos(Jet_phi[l]));
//          JECMETUpy = JECMETUpy + (sup * Jet_pt[l] * sin(Jet_phi[l]));
//        }
//        unc->setJetPt(Jet_pt[l]);
//        unc->setJetEta(Jet_eta[l]);
//        sdw = unc->getUncertainty(false);
//        if ((1-sdw)*Jet_pt[l]>30){
//          JECJetsDown->push_back(new jet_candidate((1-sdw)*Jet_pt[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepB[l], year,Jet_partonFlavour[l]));
//          JECMETDownx = JECMETDownx - (sdw * Jet_pt[l] * cos(Jet_phi[l]));
//          JECMETDowny = JECMETDowny - (sdw * Jet_pt[l] * sin(Jet_phi[l]));
//        }
//      }
//      sort(JECJetsUp->begin(), JECJetsUp->end(), ComparePtJet);
//      sort(JECJetsDown->begin(), JECJetsDown->end(), ComparePtJet);
//      JECsysUp->push_back(*JECJetsUp);
//      JECsysDown->push_back(*JECJetsDown);
//      JECsysMETUp->push_back(sqrt(pow(MET_pt * cos(MET_phi) - JECMETUpx,2)+pow(MET_pt * sin(MET_phi) - JECMETUpy,2)));
//      JECsysMETDown->push_back(sqrt(pow(MET_pt * cos(MET_phi) - JECMETDownx,2)+pow(MET_pt * sin(MET_phi) - JECMETDowny,2)));
//  }

    leading_pt = (*selectedLeptons)[0]->pt_;
    jet_leading_pt = 0;
    if (selectedJets->size()>0) jet_leading_pt = (*selectedJets)[0]->pt_;
    deltaR_ll = deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_);
    MET = MET_pt;
    n_jet = selectedJets->size();
    MVAoutput = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");



// Btag SF
    for (int l=0;l<selectedJets->size();l++){
      if((*selectedJets)[l]->btag_) nbjet++;
      if(data == "data") continue;
      if( abs((*selectedJets)[l]->flavor_) == 5){
        nbjetGen++;
        h2_BTaggingEff_Denom_b->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
        if( (*selectedJets)[l]->btag_ ) {
          h2_BTaggingEff_Num_b->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false);
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * csetBcJetSF->evaluate({"central", "M", 5, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
//          nominalWeights[4] = nominalWeights[4] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
//          sysUpWeights[4] = sysUpWeights[4] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("up", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
//          sysDownWeights[4] = sysDownWeights[4] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("down", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);

//         nominalWeights[5] = nominalWeights[5] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
//          sysUpWeights[5] = sysUpWeights[5] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
//          sysDownWeights[5] = sysDownWeights[5] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * csetBcJetSF->evaluate({"central", "M", 5, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_})));
//          nominalWeights[4] = nominalWeights[4]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
//          sysUpWeights[4] = sysUpWeights[4]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("up", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
//          sysDownWeights[4] = sysDownWeights[4]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("down", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));

//          nominalWeights[5] = nominalWeights[5]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
//          sysUpWeights[5] = sysUpWeights[5]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
//          sysDownWeights[5] = sysDownWeights[5]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }
      }
      if( abs((*selectedJets)[l]->flavor_) == 4){
        h2_BTaggingEff_Denom_c->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
        if( (*selectedJets)[l]->btag_) {
          h2_BTaggingEff_Num_c->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false);
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * csetBcJetSF->evaluate({"central", "M", 4, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
//          nominalWeights[4] = nominalWeights[4] * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
//          sysUpWeights[4] = sysUpWeights[4] * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("up", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
//          sysDownWeights[4] = sysDownWeights[4] * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("down", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);

//          nominalWeights[5] = nominalWeights[5] * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
//          sysUpWeights[5] = sysUpWeights[5] * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
//          sysDownWeights[5] = sysDownWeights[5] * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * csetBcJetSF->evaluate({"central", "M", 4, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_})));
//          nominalWeights[4] = nominalWeights[4]* (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
//          sysUpWeights[4] = sysUpWeights[4]* (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("up", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
//          sysDownWeights[4] = sysDownWeights[4]* (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("down", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));

//          nominalWeights[5] = nominalWeights[5]* (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
//          sysUpWeights[5] = sysUpWeights[5]* (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
//          sysDownWeights[5] = sysDownWeights[5]* (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }
      }

      if( abs((*selectedJets)[l]->flavor_) != 4 && abs((*selectedJets)[l]->flavor_) != 5){
        h2_BTaggingEff_Denom_udsg->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
        if( (*selectedJets)[l]->btag_) {
          h2_BTaggingEff_Num_udsg->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false);
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * csetLightJetSF->evaluate({"central", "M", 0, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
//          nominalWeights[4] = nominalWeights[4]* scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
//          sysUpWeights[4] = sysUpWeights[4]* scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
//          sysDownWeights[4] = sysDownWeights[4]* scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);

//          nominalWeights[5] = nominalWeights[5] * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
//          sysUpWeights[5] = sysUpWeights[5] * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("up", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
//          sysDownWeights[5] = sysDownWeights[5] * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("down", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * csetLightJetSF->evaluate({"central", "M", 0, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_})));
//          nominalWeights[4] = nominalWeights[4]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
//          sysUpWeights[4] = sysUpWeights[4]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
//          sysDownWeights[4] = sysDownWeights[4]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));

//          nominalWeights[5] = nominalWeights[5]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
//          sysUpWeights[5] = sysUpWeights[5]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("up", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
//          sysDownWeights[5] = sysDownWeights[5]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("down", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }
      }
    }


//    if (data == "mc"){
//      jet_leading_pt = 0;
//      if (selectedJetsJesUp->size()>0) jet_leading_pt = (*selectedJetsJesUp)[0]->pt_;
//      MET = MET_T1SmearJetEnUp_Pt;
//      n_jet = selectedJetsJesUp->size();
//      MVAoutputJesUp = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");
//      jet_leading_pt = 0;
//      if (selectedJetsJesDown->size()>0) jet_leading_pt = (*selectedJetsJesDown)[0]->pt_;
//      MET = MET_T1SmearJetEnDown_Pt;
//      n_jet = selectedJetsJesDown->size();
//      MVAoutputJesDown = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");
//  }


    if (data == "mc") weight_Lumi = (1000*xs*lumi)/Nevent;

    if (data == "mc" && (year == "2016" || year == "2017")) {
      weight_prefiring = PrefireWeight;
//      nominalWeights[7] = PrefireWeight;
//      sysUpWeights[7] = PrefireWeightup;
//      sysDownWeights[7] = PrefireWeightdown;
    }

    if (data == "mc" && ch==0) {
      sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"central",false, true);
//      nominalWeights[8] = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"central",false, true);
//      sysUpWeights[8] = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"up",false, true);
//      sysDownWeights[8] = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"down",false, true);
    }
    if (data == "mc" && ch==1) {
      sf_Trigger = scale_factor(&sf_triggeremu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"central",false, true);
//      nominalWeights[8] = scale_factor(&sf_triggeremu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"central",false, true);
//      sysUpWeights[8] = scale_factor(&sf_triggeremu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"up",false, true);
//      sysDownWeights[8] = scale_factor(&sf_triggeremu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"down",false, true);
    }
    if (data == "mc" && ch==2) {
      sf_Trigger = scale_factor(&sf_triggermumu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"central",false, true);
//      nominalWeights[8] = scale_factor(&sf_triggermumu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"central",false, true);
//      sysUpWeights[8] = scale_factor(&sf_triggermumu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"up",false, true);
//      sysDownWeights[8] = scale_factor(&sf_triggermumu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"down",false, true);
    }


    if (data == "mc" && ifTopPt) {
      for (int l=0;l<nGenPart;l++){
        if(GenPart_status[l]<30 && GenPart_status[l]>20 && GenPart_pdgId[l]==24) wp.SetPtEtaPhiM(GenPart_pt[l], GenPart_eta[l], GenPart_phi[l], GenPart_mass[l]) ;
        if(GenPart_status[l]<30 && GenPart_status[l]>20 && GenPart_pdgId[l]==-24) wm.SetPtEtaPhiM(GenPart_pt[l], GenPart_eta[l], GenPart_phi[l], GenPart_mass[l]) ;
        if(GenPart_status[l]<30 && GenPart_status[l]>20 && GenPart_pdgId[l]==5) b.SetPtEtaPhiM(GenPart_pt[l], GenPart_eta[l], GenPart_phi[l], GenPart_mass[l]) ;
        if(GenPart_status[l]<30 && GenPart_status[l]>20 && GenPart_pdgId[l]==-5) ab.SetPtEtaPhiM(GenPart_pt[l], GenPart_eta[l], GenPart_phi[l], GenPart_mass[l]) ;
        if(GenPart_status[l]<30 && GenPart_status[l]>20 && GenPart_pdgId[l]==6) top.SetPtEtaPhiM(GenPart_pt[l], GenPart_eta[l], GenPart_phi[l], GenPart_mass[l]) ;
        if(GenPart_status[l]<30 && GenPart_status[l]>20 && GenPart_pdgId[l]==-6) atop.SetPtEtaPhiM(GenPart_pt[l], GenPart_eta[l], GenPart_phi[l], GenPart_mass[l]) ;
      }
    weight_topPtPowheg = sqrt(topPtPowheg((wp + b).Pt()) * topPtPowheg((wm + ab).Pt()));
    weight_topPtMGLO = sqrt(topPtMGLO((atop).Pt()) * topPtMGLO((top).Pt()));
    }

    if (fname.Contains("LFVTt")) weight_topPtPowheg = weight_topPtMGLO;
//    if (data == "mc") weight_lep = sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO * sf_Trigger * weight_PU * signnum_typical(LHEWeight_originalXWGTUP) * weight_prefiring * weight_topPtPowheg;
//    if (data == "mc") weight_lepB = sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO * sf_Trigger * weight_PU * weight_Lumi * mc_w_sign *  weight_prefiring * weight_topPtPowheg * (P_bjet_data/P_bjet_mc);
    if (data == "mc") weight_lep = weight_Lumi * signnum_typical(LHEWeight_originalXWGTUP) * puWeight * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO * sf_Trigger;
    if (data == "mc") weight_lepB = weight_Lumi * signnum_typical(LHEWeight_originalXWGTUP) * puWeight * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger;
//if(isnan(weight_lepB) || isinf(weight_lepB)) cout<<weight_topPtMGLO<<"  "<<weight_topPtPowheg<<endl;
if(isnan(weight_lepB) || isinf(weight_lepB)) cout<<P_bjet_data<<"  "<<P_bjet_mc<<endl;
    nAccept++;

    Hists[ch][0][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    Hists[ch][0][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    Hists[ch][0][2]->Fill((*selectedLeptons)[0]->phi_,weight_lep);
    Hists[ch][0][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
    Hists[ch][0][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
    Hists[ch][0][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    Hists[ch][0][6]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep);
    Hists[ch][0][7]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lep);
    Hists[ch][0][8]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lep);
    Hists[ch][0][9]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lep);
    if(selectedJets->size()>0) Hists[ch][0][10]->Fill((*selectedJets)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][0][11]->Fill((*selectedJets)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][0][12]->Fill((*selectedJets)[0]->phi_,weight_lep);
    Hists[ch][0][13]->Fill(selectedJets->size(),weight_lep);
    Hists[ch][0][14]->Fill(nbjet,weight_lepB);
    Hists[ch][0][15]->Fill(MET_pt,weight_lep);
    Hists[ch][0][16]->Fill(MET_phi,weight_lep);
    Hists[ch][0][17]->Fill(PV_npvs,weight_lep);
    Hists[ch][0][18]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep);
    Hists[ch][0][19]->Fill(MVAoutput,weight_lep);
    if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][0][20]->Fill((*selectedLeptons)[0]->pt_  ,weight_lep);
    if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][0][22]->Fill((*selectedLeptons)[0]->eta_ ,weight_lep);
    if ((*selectedLeptons)[0]->lep_ == 1) Hists[ch][0][21]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep);
    if ((*selectedLeptons)[0]->lep_ == 1) Hists[ch][0][23]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep);
    if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][0][20]->Fill((*selectedLeptons)[1]->pt_  ,weight_lep);
    if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][0][22]->Fill((*selectedLeptons)[1]->eta_ ,weight_lep);
    if ((*selectedLeptons)[1]->lep_ == 1) Hists[ch][0][21]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep);
    if ((*selectedLeptons)[1]->lep_ == 1) Hists[ch][0][23]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep);


  if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ != 11 && ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()<106 && ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()>76) continue;
  Hists[ch][1][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
  Hists[ch][1][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
  Hists[ch][1][2]->Fill((*selectedLeptons)[0]->phi_,weight_lep);
  Hists[ch][1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
  Hists[ch][1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
  Hists[ch][1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
  Hists[ch][1][6]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep);
  Hists[ch][1][7]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lep);
  Hists[ch][1][8]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lep);
  Hists[ch][1][9]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lep);
  if(selectedJets->size()>0) Hists[ch][1][10]->Fill((*selectedJets)[0]->pt_,weight_lep);
  if(selectedJets->size()>0) Hists[ch][1][11]->Fill((*selectedJets)[0]->eta_,weight_lep);
  if(selectedJets->size()>0) Hists[ch][1][12]->Fill((*selectedJets)[0]->phi_,weight_lep);
  Hists[ch][1][13]->Fill(selectedJets->size(),weight_lep);
  Hists[ch][1][14]->Fill(nbjet,weight_lepB);
  Hists[ch][1][15]->Fill(MET_pt,weight_lep);
  Hists[ch][1][16]->Fill(MET_phi,weight_lep);
  Hists[ch][1][17]->Fill(PV_npvs,weight_lep);
  Hists[ch][1][18]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep);
  Hists[ch][1][19]->Fill(MVAoutput,weight_lep);
  if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][1][20]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep);
  if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][1][22]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep);
  if ((*selectedLeptons)[0]->lep_ == 1)  Hists[ch][1][21]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep);
  if ((*selectedLeptons)[0]->lep_ == 1)  Hists[ch][1][23]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep);
  if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][1][20]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep);
  if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][1][22]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep);
  if ((*selectedLeptons)[1]->lep_ == 1)  Hists[ch][1][21]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep);
  if ((*selectedLeptons)[1]->lep_ == 1)  Hists[ch][1][23]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep);



  if(nbjet==1){
    Hists[ch][2][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    Hists[ch][2][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    Hists[ch][2][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    Hists[ch][2][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    Hists[ch][2][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    Hists[ch][2][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    Hists[ch][2][6]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    Hists[ch][2][7]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][2][8]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB);
    Hists[ch][2][9]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB);
    Hists[ch][2][10]->Fill((*selectedJets)[0]->pt_,weight_lepB);
    Hists[ch][2][11]->Fill((*selectedJets)[0]->eta_,weight_lepB);
    Hists[ch][2][12]->Fill((*selectedJets)[0]->phi_,weight_lepB);
    Hists[ch][2][13]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][2][14]->Fill(nbjet,weight_lepB);
    Hists[ch][2][15]->Fill(MET_pt,weight_lepB);
    Hists[ch][2][16]->Fill(MET_phi,weight_lepB);
    Hists[ch][2][17]->Fill(PV_npvs,weight_lepB);
    Hists[ch][2][18]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    Hists[ch][2][19]->Fill(MVAoutput,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][2][20]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][2][22]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  Hists[ch][2][21]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  Hists[ch][2][23]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][2][20]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][2][22]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  Hists[ch][2][21]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  Hists[ch][2][23]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
  }
  if(nbjet>1){
    Hists[ch][3][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    Hists[ch][3][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    Hists[ch][3][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    Hists[ch][3][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    Hists[ch][3][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    Hists[ch][3][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    Hists[ch][3][6]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    Hists[ch][3][7]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][3][8]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB);
    Hists[ch][3][9]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB);
    Hists[ch][3][10]->Fill((*selectedJets)[0]->pt_,weight_lepB);
    Hists[ch][3][11]->Fill((*selectedJets)[0]->eta_,weight_lepB);
    Hists[ch][3][12]->Fill((*selectedJets)[0]->phi_,weight_lepB);
    Hists[ch][3][13]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][3][14]->Fill(nbjet,weight_lepB);
    Hists[ch][3][15]->Fill(MET_pt,weight_lepB);
    Hists[ch][3][16]->Fill(MET_phi,weight_lepB);
    Hists[ch][3][17]->Fill(PV_npvs,weight_lepB);
    Hists[ch][3][18]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    Hists[ch][3][19]->Fill(MVAoutput,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][3][20]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][3][22]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  Hists[ch][3][21]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  Hists[ch][3][23]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][3][20]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][3][22]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  Hists[ch][3][21]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  Hists[ch][3][23]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
  }

    for (int l=0;l<selectedLeptons->size();l++){
      delete (*selectedLeptons)[l];
    }
    for (int l=0;l<selectedLeptonsMuScaleUp->size();l++){
      delete (*selectedLeptonsMuScaleUp)[l];
    }
    for (int l=0;l<selectedLeptonsMuScaleDown->size();l++){
      delete (*selectedLeptonsMuScaleDown)[l];
    }
    for (int l=0;l<selectedLeptonsEleScaleUp->size();l++){
      delete (*selectedLeptonsEleScaleUp)[l];
    }
    for (int l=0;l<selectedLeptonsEleScaleDown->size();l++){
      delete (*selectedLeptonsEleScaleDown)[l];
    }
    for (int l=0;l<selectedLeptonsMuResUp->size();l++){
      delete (*selectedLeptonsMuResUp)[l];
    }
    for (int l=0;l<selectedLeptonsMuResDown->size();l++){
      delete (*selectedLeptonsMuResDown)[l];
    }
    selectedLeptons->clear();
    selectedLeptons->shrink_to_fit();
    delete selectedLeptons;
    selectedLeptonsMuScaleUp->clear();
    selectedLeptonsMuScaleUp->shrink_to_fit();
    delete selectedLeptonsMuScaleUp;
    selectedLeptonsMuScaleDown->clear();
    selectedLeptonsMuScaleDown->shrink_to_fit();
    delete selectedLeptonsMuScaleDown;
    selectedLeptonsMuResUp->clear();
    selectedLeptonsMuResUp->shrink_to_fit();
    delete selectedLeptonsMuResUp;
    selectedLeptonsMuResDown->clear();
    selectedLeptonsMuResDown->shrink_to_fit();
    delete selectedLeptonsMuResDown;
    selectedLeptonsEleScaleUp->clear();
    selectedLeptonsEleScaleUp->shrink_to_fit();
    delete selectedLeptonsEleScaleUp;
    selectedLeptonsEleScaleDown->clear();
    selectedLeptonsEleScaleDown->shrink_to_fit();
    delete selectedLeptonsEleScaleDown;

    for (int l=0;l<selectedJets->size();l++){
      delete (*selectedJets)[l];
    }
    for (int l=0;l<selectedJetsJerUp->size();l++){
      delete (*selectedJetsJerUp)[l];
    }
    for (int l=0;l<selectedJetsJerDown->size();l++){
      delete (*selectedJetsJerDown)[l];
    }
    for (int l=0;l<selectedJetsJesUp->size();l++){
      delete (*selectedJetsJesUp)[l];
    }
    for (int l=0;l<selectedJetsJesDown->size();l++){
      delete (*selectedJetsJesDown)[l];
    }

    selectedJets->clear();
    selectedJets->shrink_to_fit();
    delete selectedJets;
    selectedJetsJerUp->clear();
    selectedJetsJerUp->shrink_to_fit();
    delete selectedJetsJerUp;
    selectedJetsJerDown->clear();
    selectedJetsJerDown->shrink_to_fit();
    delete selectedJetsJerDown;
    selectedJetsJesUp->clear();
    selectedJetsJesUp->shrink_to_fit();
    delete selectedJetsJesUp;
    selectedJetsJesDown->clear();
    selectedJetsJesDown->shrink_to_fit();
    delete selectedJetsJesDown;

//    for (int l=0;l<JECsysUp->size();l++){
//      for (int n=0;n<(*JECsysUp)[l].size();n++){
//        delete (*JECsysUp)[l][n];
//      }
//    }
//    for (int l=0;l<JECsysDown->size();l++){
//      for (int n=0;n<(*JECsysDown)[l].size();n++){
//        delete (*JECsysDown)[l][n];
//      }
//    }
//
//    JECsysUp->clear();
//    JECsysUp->shrink_to_fit();
//    JECsysDown->clear();
//    JECsysDown->shrink_to_fit();
//    JECsysNbtagUp->clear();
//    JECsysNbtagDown->clear();
//    JECsysMETUp->clear();
//    JECsysMETDown->clear();
//    JECsysMVAUp->clear();
//    JECsysMVADown->clear();
//    JECsysMETUp->shrink_to_fit();
//    JECsysMETDown->shrink_to_fit();
//    JECsysMVAUp->shrink_to_fit();
//    JECsysMVADown->shrink_to_fit();
//    delete JECsysUp;
//    delete JECsysDown;
//    delete JECsysNbtagUp;
//    delete JECsysNbtagDown;
//    delete JECsysMETUp;
//    delete JECsysMETDown;
//    delete JECsysMVAUp;
//    delete JECsysMVADown;


   }

  cout<<"from "<<ntr<<" events, "<<nAccept<<" events are accepted"<<endl;
  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int l=0;l<vars.size();++l){
        Hists[i][k][l]  ->Write("",TObject::kOverwrite);
//        for (int n=0;n<sys.size();++n){
//          HistsSysUp[i][k][l][n]->Write("",TObject::kOverwrite);
//          HistsSysDown[i][k][l][n]->Write("",TObject::kOverwrite);
//        }
      }
    }
  }

   h2_BTaggingEff_Denom_b   ->Write("",TObject::kOverwrite);
   h2_BTaggingEff_Denom_c   ->Write("",TObject::kOverwrite);
   h2_BTaggingEff_Denom_udsg->Write("",TObject::kOverwrite);
   h2_BTaggingEff_Num_b     ->Write("",TObject::kOverwrite);
   h2_BTaggingEff_Num_c     ->Write("",TObject::kOverwrite);
   h2_BTaggingEff_Num_udsg  ->Write("",TObject::kOverwrite);
   crossSection             ->Write("",TObject::kOverwrite);
  file_out.cd("");

//  if (fname.Contains("TTTo2L2Nu") || fname.Contains("LFV")){
//    file_out.mkdir("reweightingSys");
//    file_out.cd("reweightingSys/");
//    for (int i=0;i<channels.size();++i){
//      for (int k=0;k<regions.size();++k){
//        for (int l=0;l<vars.size();++l){
//          if( channels[i] != "emu" || k<2) continue;
//          for (int n=0;n<reweightSizeQscalePDF;++n){
//            HistsSysReweightsQscalePDF[i][k][l][n]->Write("",TObject::kOverwrite);
//          }
//          for (int n=0;n<reweightSizePS;++n){
//            HistsSysReweightsPS[i][k][l][n]->Write("",TObject::kOverwrite);
//          }
//        }
//      }
//    }
//  }

  file_out.cd("");
  file_out.Close() ;
cout<<"Cleaning the memory"<<endl;

Hists.clear();
cout<<"1"<<endl;
//HistsSysUp.clear();
//cout<<"1"<<endl;
//HistsSysDown.clear();
//cout<<"1"<<endl;
//HistsJECUp.clear();
//cout<<"1"<<endl;
//HistsJECDown.clear();
//cout<<"1"<<endl;
//HistsSysReweightsQscalePDF.clear();
//cout<<"1"<<endl;
//HistsSysReweightsPS.clear();

  cout<<"Job is finished"<<endl;
}

