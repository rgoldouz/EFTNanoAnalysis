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
#include "CondFormats/Serialization/interface/Archive.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
//#include "Archive.h"
//#include "JetCorrectorParameters.h"
//#include "JetCorrectionUncertainty.h"
#include "GEScaleSyst.h"
#include "Utils.h"
#include "correction.h"
#include "WCPoint.h"
#include "WCFit.h"
#include "TH1EFT.h"
#include <map>
#endif

using namespace correction;

float topPtPowheg(float pt){
  return (0.973 - (0.000134 * pt) + (0.103 * exp(pt * (-0.0118))));
}

float topPtMGLO(float x){
  return (0.688 -  0.0000174*x + 0.824*exp(-0.0000253*x)/(pow(x,0.2185)));
}

void MyAnalysis::Loop(TString fname, TString data, TString dataset ,string year, TString run, float xs, float lumi, float Nevent, int iseft, int nRuns){
  TH1EFT  *crossSection = new TH1EFT("crossSection","crossSection",1,0,1);
  TH2F  btagEff_b_H;
  TH2F  btagEff_c_H;
  TH2F  btagEff_udsg_H;
  TH2F  sf_triggeree_H;
  TH2F  sf_triggeremu_H;
  TH2F  sf_triggermumu_H;
  TH2F  jetVetoMaps_H;
  std::string rochesterFile;
  std::string btagFile;
  GEScaleSyst *GE = new GEScaleSyst();
  PU wPU;
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
  string jetSF="";
  if(year == "2016preVFP"){
    eleSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/EGM/2016preVFP_UL/electron.json.gz";
    muSF= "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/MUO/2016preVFP_UL/muon_Z.json.gz";
    bSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/BTV/2016preVFP_UL/btagging.json.gz";
    jetSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/JME/2016preVFP_UL/UL16preVFP_jmar.json.gz";
    TFile *Map2016preVFP = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/hotjets-UL16.root");
    jetVetoMaps_H = *(TH2F*)Map2016preVFP->Get("h2hot_ul16_plus_hbm2_hbp12_qie11");
    Map2016preVFP->Close();
  }
  if(year == "2016postVFP"){
    eleSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/EGM/2016postVFP_UL/electron.json.gz";
    muSF= "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/MUO/2016postVFP_UL/muon_Z.json.gz";
    bSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/BTV/2016postVFP_UL/btagging.json.gz";
    jetSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/JME/2016postVFP_UL/UL16postVFP_jmar.json.gz";
    TFile *Map2016postVFP = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/hotjets-UL16.root");
    jetVetoMaps_H = *(TH2F*)Map2016postVFP->Get("h2hot_ul16_plus_hbm2_hbp12_qie11");
    Map2016postVFP->Close();
  }
  if(year == "2017"){
    eleSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/EGM/2017_UL/electron.json.gz";
    muSF= "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/MUO/2017_UL/muon_Z.json.gz";
    bSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/BTV/2017_UL/btagging.json.gz";
    jetSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/JME/2017_UL/UL17_jmar.json.gz";
    TFile *Map2017 = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/hotjets-UL17_v2.root");
    jetVetoMaps_H = *(TH2F*)Map2017->Get("h2hot_ul17_plus_hep17_plus_hbpw89");
    Map2017->Close();
  }
  if(year == "2018"){
    eleSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/EGM/2018_UL/electron.json.gz";
    muSF= "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/MUO/2018_UL/muon_Z.json.gz";
    bSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/BTV/2018_UL/btagging.json.gz";
    jetSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/JME/2018_UL/UL18_jmar.json.gz";
    TFile *Map2018 = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/hotjets-UL18.root");
    jetVetoMaps_H = *(TH2F*)Map2018->Get("h2hot_ul18_plus_hem1516_and_hbp2m1");
    Map2018->Close();
  }

  auto csetFileEleSF = CorrectionSet::from_file(eleSF);
  auto csetEleIdReco = csetFileEleSF->at("UL-Electron-ID-SF");

  auto csetFileMuSF = CorrectionSet::from_file(muSF);
  auto csetMuTightId = csetFileMuSF->at("NUM_TightID_DEN_TrackerMuons");
  auto csetMuTightRelIso = csetFileMuSF->at("NUM_TightRelIso_DEN_TightIDandIPCut");

  auto csetFilebSF = CorrectionSet::from_file(bSF);
  auto csetLightJetSF = csetFilebSF->at("deepJet_incl");
//  auto csetBcJetSF = csetFilebSF->at("deepCSV_mujets");
  auto csetBcJetSF = csetFilebSF->at("deepJet_comb");

  auto csetFileJetSF = CorrectionSet::from_file(jetSF);
  auto csetJetPuID = csetFileJetSF->at("PUJetID_eff");
  TRandom3 Tr;

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

  const std::map<TString, std::vector<float>> vars =
  {
    {"lep1Pt",                         {0,      60,   0,  1500}},
    {"lep1Eta",                        {1,      20,   -3, 3   }},
    {"lep1Phi",                        {2,      25,   -4, 4   }},
    {"lep2Pt",                         {3,      25,   0,  1000}},
    {"lep2Eta",                        {4,      20,   -3, 3   }},
    {"lep2Phi",                        {5,      25,   -4, 4   }},
    {"llM",                            {6,      30,    0, 500 }},
    {"llPt",                           {7,      20,    0, 200 }},
    {"llDr",                           {8,      25,    0, 7   }},
    {"llDphi",                         {9,      15,    0, 4   }},
    {"jet1Pt",                         {10,     20,    0, 300 }},
    {"jet1Eta",                        {11,     20,    -3, 3  }},
    {"jet1Phi",                        {12,     25,    -4, 4  }},
    {"njet",                           {13,     10,    0, 10  }},
    {"nbjet",                          {14,     6,     0, 6   }},
    {"Met",                            {15,     30,    0, 210 }},
    {"MetPhi",                         {16,     20,    -4, 4  }},
    {"nVtx",                           {17,     70,    0, 70  }},
    {"llMZw",                          {18,     80,    70, 110}},
    {"BDT",                            {19,     100,    -1, 1  }},
    {"topMass",                        {20,     25,    0, 500 }},
    {"topL1Dphi",                      {21,     15,    0, 4   }},
    {"topL1Dr",                        {22,     25,    0, 7   }},
    {"topL1DptOsumPt",                 {23,     20,    0, 1   }},
    {"topPt",                          {24,     25,    0, 500 }},
  };


  Dim3 Hists(channels.size(),Dim2(regions.size(),Dim1(vars.size())));
  std::stringstream name;
  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for( auto it = vars.cbegin() ; it != vars.cend() ; ++it ){
        name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first;
        std::shared_ptr<TH1EFT> h_test(new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3)));
        h_test->StatOverflows(kTRUE);
        h_test->Sumw2(kTRUE);
        Hists[i][k][it->second.at(0)] = h_test;
        name.str("");
      }
    }
  }

  std::vector<string> wc_names_lst={};
  std::vector<string> wc_names_lst_BNV={"cT", "cS"};
  std::vector<string> wc_names_lst_FCNC={"ctlS", "cte", "ctl", "ctlT", "ctZ", "cpt", "cpQM", "ctA", "cQe", "ctG", "cQlM"};
  if (fname.Contains("BNV")) wc_names_lst = wc_names_lst_BNV;
  if (fname.Contains("FCNC")) wc_names_lst = wc_names_lst_FCNC;

  std::vector<TString> sys{"eleRecoSf", "eleIDSf", "JetPuID", "muIdIsoSf", "bcTagSf", "LTagSf","pu", "prefiring", "trigSF","jes", "jer","muonScale","electronScale","muonRes", "unclusMET", "bcTagSfUnCorr", "LTagSfUnCorr" };
  Dim4 HistsSysUp(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(sys.size()))));
  Dim4 HistsSysDown(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(sys.size()))));

  if(data == "mc"){
    for (int i=0;i<channels.size();++i){
      for (int k=2;k<regions.size();++k){
        for( auto it = vars.cbegin() ; it != vars.cend() ; ++it ){
          for (int n=0;n<sys.size();++n){
            name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_"<<sys[n]<<"_Up";
            std::shared_ptr<TH1EFT> h_test(new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3)));
            h_test->StatOverflows(kTRUE);
            h_test->Sumw2(kTRUE);
            HistsSysUp[i][k][it->second.at(0)][n] = h_test;
            name.str("");
            name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_"<<sys[n]<<"_Down";
            std::shared_ptr<TH1EFT> h_test2(new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3)));
            h_test2->StatOverflows(kTRUE);
            h_test2->Sumw2(kTRUE);
            HistsSysDown[i][k][it->second.at(0)][n] = h_test2;
            name.str("");
          }
        }
      }
    }
  } 
  std::string JECFile;
  if(year == "2016preVFP")    JECFile = "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/Summer19UL16APV_V7_MC/Summer19UL16APV_V7_MC_UncertaintySources_AK4PFchs.txt";
  if(year == "2016postVFP")   JECFile = "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/Summer19UL16_V7_MC/Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.txt";
  if(year == "2017")          JECFile = "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/Summer19UL17_V5_MC/Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.txt";
  if(year == "2018")          JECFile = "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/Summer19UL18_V5_MC/Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.txt";

  std::vector<TString> sysJecNames{"AbsoluteMPFBias","AbsoluteScale","AbsoluteStat","FlavorQCD","Fragmentation","PileUpDataMC","PileUpPtBB","PileUpPtEC1","PileUpPtEC2","PileUpPtHF","PileUpPtRef","RelativeFSR","RelativePtBB","RelativePtEC1","RelativePtEC2","RelativePtHF","RelativeBal","RelativeSample","RelativeStatEC","RelativeStatFSR","RelativeStatHF","SinglePionECAL","SinglePionHCAL","TimePtEta"};
  const int nsrc = 24;
  const char* srcnames[nsrc] = {"AbsoluteMPFBias","AbsoluteScale","AbsoluteStat","FlavorQCD","Fragmentation","PileUpDataMC","PileUpPtBB","PileUpPtEC1","PileUpPtEC2","PileUpPtHF","PileUpPtRef","RelativeFSR","RelativePtBB","RelativePtEC1","RelativePtEC2","RelativePtHF","RelativeBal","RelativeSample","RelativeStatEC","RelativeStatFSR","RelativeStatHF","SinglePionECAL","SinglePionHCAL","TimePtEta"};
  std::vector<JetCorrectionUncertainty*> vsrc(nsrc);
  for (int isrc = 0; isrc < nsrc; isrc++) {
    JetCorrectorParameters *p = new JetCorrectorParameters(JECFile, srcnames[isrc]);
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
    vsrc[isrc] = unc;
  }
 
  JetCorrectorParameters *pJes = new JetCorrectorParameters(JECFile, "Total");
  JetCorrectionUncertainty *uncJes = new JetCorrectionUncertainty(*pJes);

  JetCorrectorParameters *pJer = new JetCorrectorParameters(JECFile, "RelativeJEREC1");
  JetCorrectionUncertainty *uncJer = new JetCorrectionUncertainty(*pJer);

  Dim4 HistsJECUp(channels.size(),Dim3(2,Dim2(1,Dim1(sysJecNames.size()))));
  Dim4 HistsJECDown(channels.size(),Dim3(2,Dim2(1,Dim1(sysJecNames.size()))));
  if(data == "mc"){
    for (int i=0;i<channels.size();++i){
      for (int k=2;k<regions.size();++k){
        for( auto it = vars.cbegin() ; it != vars.cend() ; ++it ){
          if(it->first !="BDT") continue;
          for (int n=0;n<sysJecNames.size();++n){
            name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_"<<sysJecNames[n]<<"_Up";
            std::shared_ptr<TH1EFT> h_test(new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3)));
            h_test->StatOverflows(kTRUE);
            h_test->Sumw2(kTRUE);
            HistsJECUp[i][k-2][0][n] = h_test;
            name.str("");
            name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_"<<sysJecNames[n]<<"_Down";
            std::shared_ptr<TH1EFT> h_test2(new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3)));
            h_test2->StatOverflows(kTRUE);
            h_test2->Sumw2(kTRUE);
            HistsJECDown[i][k-2][0][n] = h_test2;
            name.str("");
          }
        }
      }
    }
  }

// [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; [4] is renscfact=1d0 facscfact=1d0 ; [5] is renscfact=1d0 facscfact=2d0 ; [6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0 *
  int nScale = 9;
// LHA IDs NNPDF31_nnlo_hessian_pdfas 306000 - 306102*
  int nPdf = 100;
//[0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2[2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5;*
  int nPS = 4;
  Dim4 HistsSysReweightsQscale(channels.size(),Dim3(2,Dim2(1,Dim1(nScale))));
  Dim4 HistsSysReweightsPDF(channels.size(),Dim3(2,Dim2(1,Dim1(nPdf))));
  Dim4 HistsSysReweightsPS(channels.size(),Dim3(2,Dim2(1,Dim1(nPS))));

  if(data == "mc"){
    if ((fname.Contains("TTTo2L2Nu") && !fname.Contains("sys")) || fname.Contains("BNV")){
      for (int i=0;i<channels.size();++i){
        for (int k=2;k<regions.size();++k){
        for( auto it = vars.cbegin() ; it != vars.cend() ; ++it ){
            if(it->first !="BDT") continue;
            for (int n=0;n<nScale;++n){
              name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_Qscale_"<<n;
              std::shared_ptr<TH1EFT> h_test(new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3)));
              h_test->StatOverflows(kTRUE);
              h_test->Sumw2(kTRUE);
              HistsSysReweightsQscale[i][k-2][0][n] = h_test;
              name.str("");
            }
            for (int n=0;n<nPdf;++n){
              name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_PDF_"<<n;
              std::shared_ptr<TH1EFT> h_test(new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3)));
              h_test->StatOverflows(kTRUE);
              h_test->Sumw2(kTRUE);
              HistsSysReweightsPDF[i][k-2][0][n] = h_test;
              name.str("");
            }
            for (int n=0;n<nPS;++n){
              name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_PS_"<<n;
              std::shared_ptr<TH1EFT> h_test(new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3)));
              h_test->StatOverflows(kTRUE);
              h_test->Sumw2(kTRUE);
              HistsSysReweightsPS[i][k-2][0][n] = h_test;
              name.str("");
            }
          }
        }
      }
    }
  }

  TFile file_out ("ANoutput.root","RECREATE");
  float lep1Pt_;
  float lep1Eta_;
  float lep2Pt_;
  float lep2Eta_;
  float llM_;
  float llPt_;
  float llDr_;
  float llDphi_;
  float jet1Pt_;
  float jet1Eta_;
  float topMass_;
  float topL1Dphi_;
  float topL1Dr_;
  float topL1DptOsumPt_;
  float topPt_;
  float weight_;

  TTree tree_out("BNV","Top BNV analysis") ;
  tree_out.Branch("lep1Pt"      , &lep1Pt_ , "lep1Pt/F" ) ;
  tree_out.Branch("lep1Eta"      , &lep1Eta_ , "lep1Eta/F" ) ;
  tree_out.Branch("lep2Pt"      , &lep2Pt_ , "lep2Pt/F" ) ;
  tree_out.Branch("lep2Eta"      , &lep2Eta_ , "lep2Eta/F" ) ;
  tree_out.Branch("llM"      , &llM_ , "llM/F" ) ;
  tree_out.Branch("llPt"      , &llPt_ , "llPt/F" ) ;
  tree_out.Branch("llDr"      , &llDr_ , "llDr/F" ) ;
  tree_out.Branch("llDphi"      , &llDphi_ , "llDphi/F" ) ;
  tree_out.Branch("jet1Pt"      , &jet1Pt_ , "jet1Pt/F" ) ;
  tree_out.Branch("jet1Eta"      , &jet1Eta_ , "jet1Eta/F"  ) ;
  tree_out.Branch("topMass"      , &topMass_ , "topMass/F" ) ;
  tree_out.Branch("topL1Dphi"      , &topL1Dphi_ , "topL1Dphi/F" ) ;
  tree_out.Branch("topL1Dr"      , &topL1Dr_ , "topL1Dr/F" ) ;
  tree_out.Branch("topL1DptOsumPt"      , &topL1DptOsumPt_ , "topL1DptOsumPt/F" ) ;
  tree_out.Branch("topPt"      , &topPt_ , "topPt/F" ) ;
  tree_out.Branch("weight"      , &weight_ , "weight/F" ) ;

  TMVA::Tools::Instance();
  TMVA::Reader *readerMVA = new TMVA::Reader( "!Color:!Silent" );
  Float_t MVA_lep1Pt;
  Float_t MVA_lep2Pt;
  Float_t MVA_llM;
  Float_t MVA_llPt;
  Float_t MVA_llDr;
  Float_t MVA_llDphi;
  Float_t MVA_topL1Dphi;
  Float_t MVA_topL1Dr;
  Float_t MVA_topL1DptOsumPt;
  Float_t MVA_topPt;
  readerMVA->AddVariable ("lep1Pt"      , &MVA_lep1Pt) ;
  readerMVA->AddVariable ("lep2Pt"      , &MVA_lep2Pt) ;
  readerMVA->AddVariable ("llM"      , &MVA_llM) ;
  readerMVA->AddVariable ("llPt"      , &MVA_llPt) ;
  readerMVA->AddVariable ("llDr"      , &MVA_llDr) ;
  readerMVA->AddVariable ("llDphi"      , &MVA_llDphi) ;
  readerMVA->AddVariable ("topL1Dphi"      , &MVA_topL1Dphi) ;
  readerMVA->AddVariable ("topL1Dr"      , &MVA_topL1Dr) ;
  readerMVA->AddVariable ("topL1DptOsumPt"      , &MVA_topL1Dr) ;
  readerMVA->AddVariable ("topPt"      , &MVA_topPt) ;
  readerMVA->BookMVA( "BDTG", "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/TMVAClassification_BDTG.weights.xml");

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
  std::vector<float> *JECsysMETUp;
  std::vector<float> *JECsysMETDown;
  std::vector<float> *JECsysMETPhiUp;
  std::vector<float> *JECsysMETPhiDown;
  float MetJetsPtJesUp;
  float MetJetsPtJesDown;
  float MetJetsPtJerUp;
  float MetJetsPtJerDown;
  float MetJetsPhiJesUp;
  float MetJetsPhiJesDown;
  float MetJetsPhiJerUp;
  float MetJetsPhiJerDown;
  float MetXJetsJesUp;
  float MetXJetsJesDown;
  float MetXJetsJerUp;
  float MetXJetsJerDown;
  float MetYJetsJesUp;
  float MetYJetsJesDown;
  float MetYJetsJerUp;
  float MetYJetsJerDown;
  WCFit *eft_fit;


  TLorentzVector wp, wm, b, ab, top, atop;
  std::vector<float> nominalWeights;
  TLorentzVector recoTop, recoBjet, recoW, recoNu, recoL1, recoL2;
  nominalWeights.assign(sys.size(), 1);
  std::vector<float> sysUpWeights;
  sysUpWeights.assign(sys.size(), 1);
  std::vector<float> sysDownWeights;
  sysDownWeights.assign(sys.size(), 1);
  bool leptonPass;
  bool triggerPass;
  bool DyPass;
  bool triggerPassEE;
  bool triggerPassEMu;
  bool triggerPassMuMu;
  bool metFilterPass;
  bool ifTopPt=false;
  float ttKFactor;
  int ch;
  float sf_Ele_Reco;
  float sf_Ele_ID;
  float sf_Mu_ID;
  float sf_Mu_ISO;
  float sf_Trigger;
  float sf_JetPuId;
  float weight_PU;
  float weight_Lumi;
  float weight_lep;
  float weight_lepB;
  float weight_EFT;
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
  float metUnclusMETPhiUp;
  float metUnclusMETPhiDown;
  float MVAoutput;
  double P_bjet_data;
  double P_bjet_mc;
  int nAccept=0;
  float sumPuWeight=0;
  float sumPreFireWeight=0;
  float sumWeightMuID=0;
  float sumWeightMuIso=0;
  float sumWeighttTrigger=0;
  float sumWeightEleID=0;
  float sumWeightEleReco=0;
  float sumWeightBtag=0;
  float correctMuonPt=0;
  float correctMuonPtUp=0;
  float correctMuonPtDown=0;
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
  double sup = 0;
  double sdw = 0;
  bool jetlepfail;
  float BJetSF;
  float CJetSF;
  float LJetSF;
  float BJetSF_UpCorr;
  float CJetSF_UpCorr;
  float LJetSF_UpCorr;
  float BJetSF_UpUnCorr;
  float CJetSF_UpUnCorr;
  float LJetSF_UpUnCorr;
  float BJetSF_DownCorr;
  float CJetSF_DownCorr;
  float LJetSF_DownCorr;
  float BJetSF_DownUnCorr;
  float CJetSF_DownUnCorr;
  float LJetSF_DownUnCorr;
  float BJetEff;
  float CJetEff;
  float LJetEff;
  float SmearedMuonPt;

  if (fname.Contains("TTTo2L2Nu") || fname.Contains("sys") || fname.Contains("TTBNV")) ifTopPt=true;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  Long64_t ntr = fChain->GetEntries ();
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
//  for (Long64_t jentry=0; jentry<1000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    displayProgress(jentry, ntr) ;

    triggerPassEE = false;
    triggerPassEMu = false;
    triggerPassMuMu = false;
    metFilterPass = false;
    leptonPass = false;
    triggerPass = false;
    DyPass = false;
    ch =10;
    ttKFactor=1;
    sf_Ele_Reco =1;
    sf_Ele_ID =1;
    sf_Mu_ID =1;
    sf_Mu_ISO =1;
    sf_Trigger =1;
    sf_JetPuId =1;
    muPtSFRochester=1;
    weight_PU =1;
    weight_Lumi =1;
    weight_lep =1;
    weight_lepB =1;
    weight_EFT =1;
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
    metUnclusMETPhiUp=0;
    metUnclusMETPhiDown=0;

    BJetSF=1;
    CJetSF=1;
    LJetSF=1;
    BJetSF_UpCorr=1;
    CJetSF_UpCorr=1;
    LJetSF_UpCorr=1;
    BJetSF_UpUnCorr=1;
    CJetSF_UpUnCorr=1;
    LJetSF_UpUnCorr=1;
    BJetSF_DownCorr=1;
    CJetSF_DownCorr=1;
    LJetSF_DownCorr=1;
    BJetSF_DownUnCorr=1;
    CJetSF_DownUnCorr=1;
    LJetSF_DownUnCorr=1;
    BJetEff=1;
    CJetEff=1;
    LJetEff=1;
    
    for (int n=0;n<sys.size();++n){
      nominalWeights[n] =1;
      sysUpWeights[n] =1;
      sysDownWeights[n] =1;
    }
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
    crossSection->Fill(0.5, 1,*eft_fit);

    delete eft_fit; 
//You should add Flag_BadPFMuonDzFilter to the MET filter list but since it is not available in v8, lets remove it for now.
    if(year == "2017" || year == "2018"){
      if ( Flag_goodVertices  &&  Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter &&  Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_eeBadScFilter && Flag_ecalBadCalibFilter) metFilterPass = true;
    }
    else{
      if ( Flag_goodVertices  &&  Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter &&  Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_eeBadScFilter) metFilterPass = true;
    }

//trigger MC
    if(data == "mc" && (year == "2016preVFP" || year == "2016postVFP")){
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
      if(year == "2016preVFP" || year == "2016postVFP"){
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
      if(!Electron_convVeto[l]) continue;
      //remove HEM affected region
      if (year == "2018" && Electron_eta[l] < -1.3 && Electron_phi[l] < -0.87 && Electron_phi[l] > -1.57) continue;
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
        nominalWeights[0] = nominalWeights[0] * csetEleIdReco->evaluate({year, "sf", "RecoAbove20", Electron_eta[l],Electron_pt[l]});
        sysUpWeights[0] = sysUpWeights[0] * csetEleIdReco->evaluate({year, "sfup", "RecoAbove20", Electron_eta[l],Electron_pt[l]});
        sysDownWeights[0] = sysDownWeights[0] * csetEleIdReco->evaluate({year, "sfdown", "RecoAbove20", Electron_eta[l],Electron_pt[l]});

        sf_Ele_ID = sf_Ele_ID * csetEleIdReco->evaluate({year, "sf", "Tight", Electron_eta[l],Electron_pt[l]});
        nominalWeights[1] = nominalWeights[1] * csetEleIdReco->evaluate({year, "sf", "Tight", Electron_eta[l],Electron_pt[l]});
        sysUpWeights[1] = sysUpWeights[1] * csetEleIdReco->evaluate({year, "sfup", "Tight", Electron_eta[l],Electron_pt[l]});
        sysDownWeights[1] = sysDownWeights[1] * csetEleIdReco->evaluate({year, "sfdown", "Tight", Electron_eta[l],Electron_pt[l]});
      }
    }
// Muon selection
    for (int l=0;l<nMuon;l++){
      if(abs(Muon_eta[l]) > 2.4) continue;
//      if(Muon_highPtId[l]<2) continue;
//      if(Muon_tkIsoId[l]<1) continue;
      if(!Muon_tightId[l] || Muon_pfIsoId[l]<4) continue;
//      if(Muon_tunepRelPt[l]==1){
      if(data == "data") {
        muPtSFRochester = rc.kScaleDT(Muon_charge[l], Muon_pt[l],Muon_eta[l],Muon_phi[l], 0, 0);
        if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
      }
      if (data == "mc"){
        if (Muon_genPartIdx[l]>=0) muPtSFRochester = rc.kSpreadMC(Muon_charge[l], Muon_pt[l],Muon_eta[l],Muon_phi[l], GenPart_pt[Muon_genPartIdx[l]],0, 0);
        if (Muon_genPartIdx[l]<0) muPtSFRochester = rc.kSmearMC(Muon_charge[l], Muon_pt[l],Muon_eta[l],Muon_phi[l], Muon_nTrackerLayers[l] , gRandom->Rndm(),0, 0);
        if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
        if(Muon_pt[l]<200){
          if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptonsEleScaleUp->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptonsEleScaleDown->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptonsMuResUp->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptonsMuResDown->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptonsMuScaleUp->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          if(muPtSFRochester * Muon_pt[l] > 20) selectedLeptonsMuScaleDown->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
        }
        else{
          pt_res = 0.02;
          SmearedMuonPt = Tr.Gaus(0, pt_res);
          correctMuonPt = GE->GEScaleCorrPt(1600, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]);
          correctMuonPtUp = GE->GEScaleCorrPt(1601, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]);
          correctMuonPtDown = GE->GEScaleCorrPt(1602, Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]);
          selectedLeptonsEleScaleUp->push_back(new lepton_candidate(Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          selectedLeptonsEleScaleDown->push_back(new lepton_candidate(Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          selectedLeptonsMuResUp->push_back(new lepton_candidate((1+abs(SmearedMuonPt))*Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          selectedLeptonsMuResDown->push_back(new lepton_candidate((1-abs(SmearedMuonPt))*Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          selectedLeptonsMuScaleUp->push_back(new lepton_candidate(Muon_pt[l] + abs((correctMuonPtUp-correctMuonPt)) ,Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          selectedLeptonsMuScaleDown->push_back(new lepton_candidate(Muon_pt[l] - abs((correctMuonPtDown-correctMuonPt)),Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
        }
      }
      if(Muon_pt[l] <20) continue;
      if (data == "mc") {
        sf_Mu_ID = sf_Mu_ID * csetMuTightId->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_pt[l], "sf"});
//        nominalWeights[2] = nominalWeights[2] * csetMuTightId->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_pt[l], "sf"});
//        sysUpWeights[2] = sysUpWeights[2] * csetMuTightId->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_pt[l], "systup"});
//        sysDownWeights[2] = sysDownWeights[2] * csetMuTightId->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_pt[l], "systdown"});

        sf_Mu_ISO = sf_Mu_ISO * csetMuTightRelIso->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_pt[l], "sf"});
        nominalWeights[3] = nominalWeights[3] * csetMuTightRelIso->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_pt[l], "sf"}) * csetMuTightId->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_pt[l], "sf"});
        sysUpWeights[3] = sysUpWeights[3] * csetMuTightRelIso->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_pt[l], "systup"}) * csetMuTightId->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_pt[l], "systup"});
        sysDownWeights[3] = sysDownWeights[3] * csetMuTightRelIso->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_pt[l], "systdown"}) * csetMuTightId->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_pt[l], "systdown"});
      }
    }
    sort(selectedLeptons->begin(), selectedLeptons->end(), ComparePtLep);
    sort(selectedLeptonsMuScaleUp->begin(), selectedLeptonsMuScaleUp->end(), ComparePtLep);
    sort(selectedLeptonsMuScaleDown->begin(), selectedLeptonsMuScaleDown->end(), ComparePtLep);
    sort(selectedLeptonsMuResDown->begin(), selectedLeptonsMuResDown->end(), ComparePtLep);
    sort(selectedLeptonsMuResUp->begin(), selectedLeptonsMuResUp->end(), ComparePtLep);
    sort(selectedLeptonsEleScaleUp->begin(), selectedLeptonsEleScaleUp->end(), ComparePtLep);
    sort(selectedLeptonsEleScaleDown->begin(), selectedLeptonsEleScaleDown->end(), ComparePtLep);

//jets
    selectedJets = new std::vector<jet_candidate*>();
    selectedJetsJerUp = new std::vector<jet_candidate*>();
    selectedJetsJerDown = new std::vector<jet_candidate*>();
    selectedJetsJesUp = new std::vector<jet_candidate*>();
    selectedJetsJesDown = new std::vector<jet_candidate*>();
    MetJetsPtJesUp=0;
    MetJetsPtJesDown=0;
    MetJetsPtJerUp=0;
    MetJetsPtJerDown=0;
    MetJetsPhiJesUp=0;
    MetJetsPhiJesDown=0;
    MetJetsPhiJerUp=0;
    MetJetsPhiJerDown=0;
    MetXJetsJesUp=0;
    MetXJetsJesDown=0;
    MetXJetsJerUp=0;
    MetXJetsJerDown=0;
    MetYJetsJesUp=0;
    MetYJetsJesDown=0;
    MetYJetsJerUp=0;
    MetYJetsJerDown=0;
    for (int l=0;l<nJet;l++){
      if(Jet_jetId[l]<6) continue;
      if(Jet_puId[l]<1 && Jet_pt_nom[l]<50) continue;
//Jet Veto map
      if (jetVetoMaps_H.GetBinContent(jetVetoMaps_H.GetXaxis()->FindBin(Jet_eta[l]),jetVetoMaps_H.GetYaxis()->FindBin(Jet_phi[l]))>0) continue;
      jetlepfail = false;
      for (int i=0;i<selectedLeptons->size();i++){
        if(deltaR((*selectedLeptons)[i]->eta_,(*selectedLeptons)[i]->phi_,Jet_eta[l],Jet_phi[l]) < 0.4 ) jetlepfail=true;
      }
      if(jetlepfail) continue;
      if(data == "mc" && abs(Jet_eta[l]) < 2.4){
        if(Jet_pt_nom[l] >30) {
          selectedJets->push_back(new jet_candidate(Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l], year,Jet_partonFlavour[l]));
          if(Jet_pt_nom[l] <50){
            sf_JetPuId = sf_JetPuId * csetJetPuID->evaluate({Jet_eta[l],Jet_pt_nom[l],"nom","L"});
            nominalWeights[2] = nominalWeights[2] * csetJetPuID->evaluate({Jet_eta[l],Jet_pt_nom[l],"nom","L"});
            sysUpWeights[2] = sysUpWeights[2] * csetJetPuID->evaluate({Jet_eta[l],Jet_pt_nom[l],"up","L"});
            sysDownWeights[2] = sysDownWeights[2] * csetJetPuID->evaluate({Jet_eta[l],Jet_pt_nom[l],"down","L"});
          }
        }
        uncJes->setJetPt(Jet_pt_nom[l]);
        uncJes->setJetEta(Jet_eta[l]);
        sup = uncJes->getUncertainty(true);
        if ((1+sup)*Jet_pt_nom[l]>30) {
          selectedJetsJesUp->push_back(new jet_candidate((1+sup)*Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l], year,Jet_partonFlavour[l]));
          MetXJetsJesUp = MetXJetsJesUp + (sup * Jet_pt_nom[l] * cos(Jet_phi[l]));
          MetYJetsJesUp = MetYJetsJesUp + (sup * Jet_pt_nom[l] * sin(Jet_phi[l]));
        }

        uncJes->setJetPt(Jet_pt_nom[l]);
        uncJes->setJetEta(Jet_eta[l]);
        sdw = uncJes->getUncertainty(false);
        if ((1-sdw)*Jet_pt_nom[l]>30) {
          selectedJetsJesDown->push_back(new jet_candidate((1-sdw)*Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l], year,Jet_partonFlavour[l]));
          MetXJetsJesDown = MetXJetsJesDown - (sdw * Jet_pt_nom[l] * cos(Jet_phi[l]));
          MetYJetsJesDown = MetYJetsJesDown - (sdw * Jet_pt_nom[l] * sin(Jet_phi[l]));
        }

        uncJer->setJetPt(Jet_pt_nom[l]);
        uncJer->setJetEta(Jet_eta[l]);
        sup = uncJer->getUncertainty(true);
        if ((1+sup)*Jet_pt_nom[l]>30) {
          selectedJetsJerUp->push_back(new jet_candidate((1+sup)*Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l], year,Jet_partonFlavour[l]));
          MetXJetsJerUp = MetXJetsJerUp + (sup * Jet_pt_nom[l] * cos(Jet_phi[l]));
          MetYJetsJerUp = MetYJetsJerUp + (sup * Jet_pt_nom[l] * sin(Jet_phi[l]));
        }

        uncJer->setJetPt(Jet_pt_nom[l]);
        uncJer->setJetEta(Jet_eta[l]);
        sdw = uncJer->getUncertainty(false);
        if ((1-sdw)*Jet_pt_nom[l]>30) {
          selectedJetsJerDown->push_back(new jet_candidate((1-sdw)*Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l], year,Jet_partonFlavour[l]));
          MetXJetsJerDown = MetXJetsJerDown - (sdw * Jet_pt_nom[l] * cos(Jet_phi[l]));
          MetYJetsJerDown = MetYJetsJerDown - (sdw * Jet_pt_nom[l] * sin(Jet_phi[l]));
        }
      }
      if(data == "data" && Jet_pt_nom[l] >30 && abs(Jet_eta[l]) < 2.4){
        selectedJets->push_back(new jet_candidate(Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l],year,0));
        selectedJetsJerUp->push_back(new jet_candidate(Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l],year,0));
        selectedJetsJerDown->push_back(new jet_candidate(Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l],year,0));
        selectedJetsJesUp->push_back(new jet_candidate(Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l],year,0));
        selectedJetsJesDown->push_back(new jet_candidate(Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l],year,0));
      }
    }

    MetJetsPtJesUp = sqrt(pow(MET_T1_pt * cos(MET_T1_phi) - MetXJetsJesUp,2)+pow(MET_T1_pt * sin(MET_T1_phi) - MetYJetsJesUp,2));
    MetJetsPhiJesUp = atan2 ((MET_T1_pt * sin(MET_T1_phi) - MetYJetsJesUp),(MET_T1_pt * cos(MET_T1_phi) - MetXJetsJesUp));
    MetJetsPtJesDown = sqrt(pow(MET_T1_pt * cos(MET_T1_phi) - MetXJetsJesDown,2)+pow(MET_T1_pt * sin(MET_T1_phi) - MetYJetsJesDown,2));
    MetJetsPhiJesDown = atan2 ((MET_T1_pt * sin(MET_T1_phi) - MetYJetsJesDown),(MET_T1_pt * cos(MET_T1_phi) - MetXJetsJesDown));
    MetJetsPtJerUp = sqrt(pow(MET_T1_pt * cos(MET_T1_phi) - MetXJetsJerUp,2)+pow(MET_T1_pt * sin(MET_T1_phi) - MetYJetsJerUp,2));
    MetJetsPhiJerUp = atan2 ((MET_T1_pt * sin(MET_T1_phi) - MetYJetsJerUp),(MET_T1_pt * cos(MET_T1_phi) - MetXJetsJerUp));
    MetJetsPtJerDown = sqrt(pow(MET_T1_pt * cos(MET_T1_phi) - MetXJetsJerDown,2)+pow(MET_T1_pt * sin(MET_T1_phi) - MetYJetsJerDown,2));
    MetJetsPhiJerDown = atan2 ((MET_T1_pt * sin(MET_T1_phi) - MetYJetsJerDown),(MET_T1_pt * cos(MET_T1_phi) - MetXJetsJerDown));

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
    JECsysUp = new std::vector<std::vector<jet_candidate*>>();
    JECsysDown = new std::vector<std::vector<jet_candidate*>>();
    JECsysMETUp = new std::vector<float>();
    JECsysMETDown = new std::vector<float>();
    JECsysMETPhiUp = new std::vector<float>();
    JECsysMETPhiDown = new std::vector<float>();
    JECsysNbtagUp = new std::vector<int>();
    JECsysNbtagDown = new std::vector<int>();
    if(data == "mc"){
      for (int n=0;n<sysJecNames.size();++n){
        JECJetsUp= new std::vector<jet_candidate*>();
        JECJetsDown= new std::vector<jet_candidate*>();
        JECMETUpx =   0;
        JECMETUpy =   0;
        JECMETDownx = 0;
        JECMETDowny = 0;
        for (int l=0;l<nJet;l++){
          if(abs(Jet_eta[l]) > 2.4) continue;
          if(Jet_puId[l]<1 && Jet_pt_nom[l]<50) continue;
          if(Jet_jetId[l]<6) continue;
          if (jetVetoMaps_H.GetBinContent(jetVetoMaps_H.GetXaxis()->FindBin(Jet_eta[l]),jetVetoMaps_H.GetYaxis()->FindBin(Jet_phi[l]))>0) continue;
          jetlepfail = false;
          for (int i=0;i<selectedLeptons->size();i++){
            if(deltaR((*selectedLeptons)[i]->eta_,(*selectedLeptons)[i]->phi_,Jet_eta[l],Jet_phi[l]) < 0.4 ) jetlepfail=true;
          }
          if(jetlepfail) continue;
          JetCorrectionUncertainty *unc = vsrc[n];
          unc->setJetPt(Jet_pt_nom[l]);
          unc->setJetEta(Jet_eta[l]);
          sup = unc->getUncertainty(true);
          if ((1+sup)*Jet_pt_nom[l]>30) {
            JECJetsUp->push_back(new jet_candidate((1+sup)*Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l], year,Jet_partonFlavour[l]));
            JECMETUpx = JECMETUpx + (sup * Jet_pt_nom[l] * cos(Jet_phi[l]));
            JECMETUpy = JECMETUpy + (sup * Jet_pt_nom[l] * sin(Jet_phi[l]));
          }
          unc->setJetPt(Jet_pt_nom[l]);
          unc->setJetEta(Jet_eta[l]);
          sdw = unc->getUncertainty(false);
          if ((1-sdw)*Jet_pt_nom[l]>30){
            JECJetsDown->push_back(new jet_candidate((1-sdw)*Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l], year,Jet_partonFlavour[l]));
            JECMETDownx = JECMETDownx - (sdw * Jet_pt_nom[l] * cos(Jet_phi[l]));
            JECMETDowny = JECMETDowny - (sdw * Jet_pt_nom[l] * sin(Jet_phi[l]));
          }
        }
        sort(JECJetsUp->begin(), JECJetsUp->end(), ComparePtJet);
        sort(JECJetsDown->begin(), JECJetsDown->end(), ComparePtJet);
        JECsysUp->push_back(*JECJetsUp);
        JECsysDown->push_back(*JECJetsDown);
        JECsysMETUp->push_back(sqrt(pow(MET_T1_pt * cos(MET_T1_phi) - JECMETUpx,2)+pow(MET_T1_pt * sin(MET_T1_phi) - JECMETUpy,2)));
        JECsysMETDown->push_back(sqrt(pow(MET_T1_pt * cos(MET_T1_phi) - JECMETDownx,2)+pow(MET_T1_pt * sin(MET_T1_phi) - JECMETDowny,2)));
        JECsysMETPhiUp->push_back(atan2 ((MET_T1_pt * sin(MET_T1_phi) - JECMETUpx),(MET_T1_pt * cos(MET_T1_phi) - JECMETUpy)));
        JECsysMETPhiDown->push_back(atan2 ((MET_T1_pt * sin(MET_T1_phi) - JECMETDownx),(MET_T1_pt * cos(MET_T1_phi) - JECMETDowny)));
      }
    
      for (int n=0;n<sysJecNames.size();++n){
        int JECnbUp = 0;
        int JECnbDown = 0;
        for (int i=0;i<(*JECsysUp)[n].size();i++){
          if((*JECsysUp)[n][i]->btag_) JECnbUp++;
        }
        for (int i=0;i<(*JECsysDown)[n].size();i++){
          if((*JECsysDown)[n][i]->btag_) JECnbDown++;
        }
        JECsysNbtagUp->push_back(JECnbUp);
        JECsysNbtagDown->push_back(JECnbDown);
      }
    }
// Btag SF
    for (int l=0;l<selectedJets->size();l++){
      if((*selectedJets)[l]->btag_) nbjet++;
      if(data == "data") continue;
      BJetSF=csetBcJetSF->evaluate({"central", "M", 5, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      CJetSF=csetBcJetSF->evaluate({"central", "M", 4, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      LJetSF=csetLightJetSF->evaluate({"central", "M", 0, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      BJetSF_UpCorr=csetBcJetSF->evaluate({"up_correlated", "M", 5, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      CJetSF_UpCorr=csetBcJetSF->evaluate({"up_correlated", "M", 4, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      LJetSF_UpCorr=csetLightJetSF->evaluate({"up_correlated", "M", 0, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      BJetSF_UpUnCorr=csetBcJetSF->evaluate({"up_uncorrelated", "M", 5, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      CJetSF_UpUnCorr=csetBcJetSF->evaluate({"up_uncorrelated", "M", 4, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      LJetSF_UpUnCorr=csetLightJetSF->evaluate({"up_uncorrelated", "M", 0, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      BJetSF_DownCorr=csetBcJetSF->evaluate({"down_correlated", "M", 5, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      CJetSF_DownCorr=csetBcJetSF->evaluate({"down_correlated", "M", 4, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      LJetSF_DownCorr=csetLightJetSF->evaluate({"down_correlated", "M", 0, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      BJetSF_DownUnCorr=csetBcJetSF->evaluate({"down_uncorrelated", "M", 5, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      CJetSF_DownUnCorr=csetBcJetSF->evaluate({"down_uncorrelated", "M", 4, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      LJetSF_DownUnCorr=csetLightJetSF->evaluate({"down_uncorrelated", "M", 0, abs((*selectedJets)[l]->eta_),(*selectedJets)[l]->pt_});
      BJetEff=scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false);
      CJetEff=scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false);
      LJetEff=scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false);

//b-quark
      if( abs((*selectedJets)[l]->flavor_) == 5){
        nbjetGen++;
        h2_BTaggingEff_Denom_b->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
        if( (*selectedJets)[l]->btag_ ) {
          h2_BTaggingEff_Num_b->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
          P_bjet_mc = P_bjet_mc * BJetEff;
          P_bjet_data = P_bjet_data * BJetEff * BJetSF;
          nominalWeights[4] = nominalWeights[4] * BJetEff * BJetSF;
          sysUpWeights[4] = sysUpWeights[4] * BJetEff * BJetSF_UpCorr;
          sysDownWeights[4] = sysDownWeights[4] * BJetEff * BJetSF_DownCorr;
          nominalWeights[5] = nominalWeights[5] * BJetEff * BJetSF;
          sysUpWeights[5] = sysUpWeights[5] * BJetEff * BJetSF;
          sysDownWeights[5] = sysDownWeights[5] * BJetEff * BJetSF;

          nominalWeights[15] = nominalWeights[15] * BJetEff * BJetSF;
          sysUpWeights[15] = sysUpWeights[15] * BJetEff * BJetSF_UpUnCorr;
          sysDownWeights[15] = sysDownWeights[15] * BJetEff * BJetSF_DownUnCorr;
          nominalWeights[16] = nominalWeights[16] * BJetEff * BJetSF;
          sysUpWeights[16] = sysUpWeights[16] * BJetEff * BJetSF;
          sysDownWeights[16] = sysDownWeights[16] * BJetEff * BJetSF;
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - BJetEff);
          P_bjet_data = P_bjet_data * (1- (BJetEff * BJetSF));
          nominalWeights[4] = nominalWeights[4]* (1- (BJetEff * BJetSF));
          sysUpWeights[4] = sysUpWeights[4]* (1- (BJetEff * BJetSF_UpCorr));
          sysDownWeights[4] = sysDownWeights[4]* (1- (BJetEff * BJetSF_DownCorr));
          nominalWeights[5] = nominalWeights[5]* (1- (BJetEff * BJetSF));
          sysUpWeights[5] = sysUpWeights[5]* (1- (BJetEff * BJetSF));
          sysDownWeights[5] = sysDownWeights[5]* (1- (BJetEff * BJetSF));

          nominalWeights[15] = nominalWeights[15]* (1- (BJetEff * BJetSF));
          sysUpWeights[15] = sysUpWeights[15]* (1- (BJetEff * BJetSF_UpUnCorr));
          sysDownWeights[15] = sysDownWeights[15]* (1- (BJetEff * BJetSF_DownUnCorr));
          nominalWeights[16] = nominalWeights[16]* (1- (BJetEff * BJetSF));
          sysUpWeights[16] = sysUpWeights[16]* (1- (BJetEff * BJetSF));
          sysDownWeights[16] = sysDownWeights[16]* (1- (BJetEff * BJetSF));
        }
      }
//c-quark
      if( abs((*selectedJets)[l]->flavor_) == 4){
        h2_BTaggingEff_Denom_c->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
        if( (*selectedJets)[l]->btag_) {
          h2_BTaggingEff_Num_c->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
          P_bjet_mc = P_bjet_mc * CJetEff;
          P_bjet_data = P_bjet_data * CJetEff * CJetSF;
          nominalWeights[4] = nominalWeights[4] * CJetEff * CJetSF;
          sysUpWeights[4] = sysUpWeights[4] * CJetEff * CJetSF_UpCorr;
          sysDownWeights[4] = sysDownWeights[4] * CJetEff * CJetSF_DownCorr;
          nominalWeights[5] = nominalWeights[5] * CJetEff * CJetSF;
          sysUpWeights[5] = sysUpWeights[5] * CJetEff * CJetSF;
          sysDownWeights[5] = sysDownWeights[5] * CJetEff * CJetSF;

          nominalWeights[15] = nominalWeights[15] * CJetEff * CJetSF;
          sysUpWeights[15] = sysUpWeights[15] * CJetEff * CJetSF_UpUnCorr;
          sysDownWeights[15] = sysDownWeights[15] * CJetEff * CJetSF_DownUnCorr;
          nominalWeights[16] = nominalWeights[16] * CJetEff * CJetSF;
          sysUpWeights[16] = sysUpWeights[16] * CJetEff * CJetSF;
          sysDownWeights[16] = sysDownWeights[16] * CJetEff * CJetSF;
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - CJetEff);
          P_bjet_data = P_bjet_data * (1- (CJetEff * CJetSF));
          nominalWeights[4] = nominalWeights[4]* (1- (CJetEff * CJetSF));
          sysUpWeights[4] = sysUpWeights[4]* (1- (CJetEff * CJetSF_UpCorr));
          sysDownWeights[4] = sysDownWeights[4]* (1- (CJetEff * CJetSF_DownCorr));
          nominalWeights[5] = nominalWeights[5]* (1- (CJetEff * CJetSF));
          sysUpWeights[5] = sysUpWeights[5]* (1- (CJetEff * CJetSF));
          sysDownWeights[5] = sysDownWeights[5]* (1- (CJetEff * CJetSF));

          nominalWeights[15] = nominalWeights[15]* (1- (CJetEff * CJetSF));
          sysUpWeights[15] = sysUpWeights[15]* (1- (CJetEff * CJetSF_UpUnCorr));
          sysDownWeights[15] = sysDownWeights[15]* (1- (CJetEff * CJetSF_DownUnCorr));
          nominalWeights[16] = nominalWeights[16]* (1- (CJetEff * CJetSF));
          sysUpWeights[16] = sysUpWeights[16]* (1- (CJetEff * CJetSF));
          sysDownWeights[16] = sysDownWeights[16]* (1- (CJetEff * CJetSF));
        }
      }
//light-quark
      if( abs((*selectedJets)[l]->flavor_) != 4 && abs((*selectedJets)[l]->flavor_) != 5){
        h2_BTaggingEff_Denom_udsg->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
        if( (*selectedJets)[l]->btag_) {
          h2_BTaggingEff_Num_udsg->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
          P_bjet_mc = P_bjet_mc * LJetEff;
          P_bjet_data = P_bjet_data * LJetEff * LJetSF;
          nominalWeights[4] = nominalWeights[4]* LJetEff * LJetSF;
          sysUpWeights[4] = sysUpWeights[4]* LJetEff * LJetSF;
          sysDownWeights[4] = sysDownWeights[4]* LJetEff * LJetSF;
          nominalWeights[5] = nominalWeights[5] * LJetEff * LJetSF;
          sysUpWeights[5] = sysUpWeights[5] * LJetEff * LJetSF_UpCorr;
          sysDownWeights[5] = sysDownWeights[5] * LJetEff * LJetSF_DownCorr;

          nominalWeights[15] = nominalWeights[15]* LJetEff * LJetSF;
          sysUpWeights[15] = sysUpWeights[15]* LJetEff * LJetSF;
          sysDownWeights[15] = sysDownWeights[15]* LJetEff * LJetSF;
          nominalWeights[16] = nominalWeights[16] * LJetEff * LJetSF;
          sysUpWeights[16] = sysUpWeights[16] * LJetEff * LJetSF_UpUnCorr;
          sysDownWeights[16] = sysDownWeights[16] * LJetEff * LJetSF_DownUnCorr;
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - LJetEff);
          P_bjet_data = P_bjet_data * (1- (LJetEff * LJetSF));
          nominalWeights[4] = nominalWeights[4]* (1- (LJetEff * LJetSF));
          sysUpWeights[4] = sysUpWeights[4]* (1- (LJetEff * LJetSF));
          sysDownWeights[4] = sysDownWeights[4]* (1- (LJetEff * LJetSF));
          nominalWeights[5] = nominalWeights[5]* (1- (LJetEff * LJetSF));
          sysUpWeights[5] = sysUpWeights[5]* (1- (LJetEff * LJetSF_UpCorr));
          sysDownWeights[5] = sysDownWeights[5]* (1- (LJetEff * LJetSF_DownCorr));

          nominalWeights[15] = nominalWeights[15]* (1- (LJetEff * LJetSF));
          sysUpWeights[15] = sysUpWeights[15]* (1- (LJetEff * LJetSF));
          sysDownWeights[15] = sysDownWeights[15]* (1- (LJetEff * LJetSF));
          nominalWeights[16] = nominalWeights[16]* (1- (LJetEff * LJetSF));
          sysUpWeights[16] = sysUpWeights[16]* (1- (LJetEff * LJetSF_UpUnCorr));
          sysDownWeights[16] = sysDownWeights[16]* (1- (LJetEff * LJetSF_DownUnCorr));
        }
      }
    }

//MET Unc 

    metUnclusMETUp = sqrt(pow(MET_T1_pt * cos(MET_T1_phi) + MET_MetUnclustEnUpDeltaX,2)+pow(MET_T1_pt * sin(MET_T1_phi) + MET_MetUnclustEnUpDeltaY,2));
    metUnclusMETDown = sqrt(pow(MET_T1_pt * cos(MET_T1_phi) - MET_MetUnclustEnUpDeltaX,2)+pow(MET_T1_pt * sin(MET_T1_phi) - MET_MetUnclustEnUpDeltaY,2));
    metUnclusMETPhiUp = atan2 ((MET_T1_pt * sin(MET_T1_phi) + MET_MetUnclustEnUpDeltaY),(MET_T1_pt * cos(MET_T1_phi) + MET_MetUnclustEnUpDeltaX));
    metUnclusMETPhiDown = atan2 ((MET_T1_pt * sin(MET_T1_phi) - MET_MetUnclustEnUpDeltaY),(MET_T1_pt * cos(MET_T1_phi) - MET_MetUnclustEnUpDeltaX));
//PU reweighting
    if (data == "mc" && year == "2016preVFP") {
      weight_PU = wPU.PU_2016preVFP(int(Pileup_nTrueInt),"nominal");
      nominalWeights[6] = wPU.PU_2016preVFP(int(Pileup_nTrueInt),"nominal");
      sysUpWeights[6] = wPU.PU_2016preVFP(int(Pileup_nTrueInt),"up");
      sysDownWeights[6] = wPU.PU_2016preVFP(int(Pileup_nTrueInt),"down");
    }
    if (data == "mc" && year == "2016postVFP") {
      weight_PU = wPU.PU_2016postVFP(int(Pileup_nTrueInt),"nominal");
      nominalWeights[6] = wPU.PU_2016postVFP(int(Pileup_nTrueInt),"nominal");
      sysUpWeights[6] = wPU.PU_2016postVFP(int(Pileup_nTrueInt),"up");
      sysDownWeights[6] = wPU.PU_2016postVFP(int(Pileup_nTrueInt),"down");
    }
    if (data == "mc" && year == "2017") {
      weight_PU = wPU.PU_2017(int(Pileup_nTrueInt),"nominal");
      nominalWeights[6] = wPU.PU_2017(int(Pileup_nTrueInt),"nominal");
      sysUpWeights[6] = wPU.PU_2017(int(Pileup_nTrueInt),"up");
      sysDownWeights[6] = wPU.PU_2017(int(Pileup_nTrueInt),"down");
    }
    if (data == "mc" && year == "2018") {
      weight_PU = wPU.PU_2018(int(Pileup_nTrueInt),"nominal");
      nominalWeights[6] = wPU.PU_2018(int(Pileup_nTrueInt),"nominal");
      sysUpWeights[6] = wPU.PU_2018(int(Pileup_nTrueInt),"up");
      sysDownWeights[6] = wPU.PU_2018(int(Pileup_nTrueInt),"down");
    }


    if (data == "mc") weight_Lumi = (1000*xs*lumi)/Nevent;
    if (data == "mc"){
        weight_prefiring = L1PreFiringWeight_Nom;
        nominalWeights[7] = L1PreFiringWeight_Nom;
        sysUpWeights[7] = L1PreFiringWeight_Up;
        sysDownWeights[7] = L1PreFiringWeight_Dn;
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

    if (fname.Contains("TTBNV")) {
      weight_topPtPowheg = weight_topPtMGLO;
      ttKFactor=831.7/445.0;
    }

//remove HEM effected region from 2018 samples
//    if (year == "2018" && (*selectedLeptons)[0]->lep_ == 1 && (*selectedLeptons)[0]->eta_ < -1.4 && (*selectedLeptons)[0]->phi_< -0.8 && (*selectedLeptons)[0]->phi_ > -1.6) continue;
//    if (year == "2018" && (*selectedLeptons)[1]->lep_ == 1 && (*selectedLeptons)[1]->eta_ < -1.4 && (*selectedLeptons)[1]->phi_< -0.8 && (*selectedLeptons)[1]->phi_ > -1.6) continue;

//Check if analysis cuts are passed and then categorize dilepton channels
    if(selectedLeptons->size() ==2){ 
      if (((*selectedLeptons)[0]->pt_ > 25) &&
      ((*selectedLeptons)[0]->charge_ * (*selectedLeptons)[1]->charge_ == -1) &&
      ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()>20) leptonPass=true; 
    
      if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 2) ch = 0;
      if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 11) ch = 1;
      if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 20) ch = 2;
      if(ch ==0 && triggerPassEE) triggerPass=true;
      if(ch ==1 && triggerPassEMu) triggerPass=true;
      if(ch ==2 && triggerPassMuMu) triggerPass=true;
      DyPass = true;
      if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ != 11 && ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()<106 ) DyPass = false;

      if (data == "mc" && ch==0) {
        sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"central",false, true);
        nominalWeights[8] = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"central",false, true);
        sysUpWeights[8] = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"up",false, true);
        sysDownWeights[8] = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"down",false, true);
      }
      if (data == "mc" && ch==1) {
        sf_Trigger = scale_factor(&sf_triggeremu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"central",false, true);
        nominalWeights[8] = scale_factor(&sf_triggeremu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"central",false, true);
        sysUpWeights[8] = scale_factor(&sf_triggeremu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"up",false, true);
        sysDownWeights[8] = scale_factor(&sf_triggeremu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"down",false, true);
      }
      if (data == "mc" && ch==2) {
        sf_Trigger = scale_factor(&sf_triggermumu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"central",false, true);
        nominalWeights[8] = scale_factor(&sf_triggermumu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"central",false, true);
        sysUpWeights[8] = scale_factor(&sf_triggermumu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"up",false, true);
        sysDownWeights[8] = scale_factor(&sf_triggermumu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"down",false, true);
      }
    }
//Fill histograms
    if(leptonPass && triggerPass){
      if(ch==1 && nbjet>1){
      nAccept++;
      sumPuWeight+=weight_PU;
      sumPreFireWeight+=weight_prefiring;
      sumWeightMuID+=sf_Mu_ID;
      sumWeightMuIso+=sf_Mu_ISO;
      sumWeighttTrigger+=sf_Trigger;
      sumWeightEleID+=sf_Ele_ID;
      sumWeightEleReco+=sf_Ele_Reco;
      sumWeightBtag+=(P_bjet_data/P_bjet_mc);
      }
      if (data == "mc") weight_lep  = weight_Lumi * signnum_typical(LHEWeight_originalXWGTUP) * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowheg * ttKFactor * sf_JetPuId;
      if (data == "mc") weight_lepB = weight_lep * (P_bjet_data/P_bjet_mc);
      if (data == "mc") weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowheg * ttKFactor * sf_JetPuId;
      if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
      else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
//if(nbjet>1)   cout<<nbjet<<" = "<<P_bjet_data<<"/"<<P_bjet_mc<<"  = '"<<(P_bjet_data/P_bjet_mc)<<endl;
//cout<<"sf_Ele_Reco:"<<sf_Ele_Reco<<" sf_Ele_ID:"<<sf_Ele_ID<<" sf_Mu_ID:"<<sf_Mu_ID<<" sf_Mu_ISO:"<<sf_Mu_ISO<<" sf_Trigger:"<<sf_Trigger<<" weight_PU:"<<weight_PU<<" weight_prefiring:"<<weight_prefiring<< " weight_topPtPowheg:"<<weight_topPtPowheg<< "P_bjet_data:"<<P_bjet_data<<" P_bjet_mc:"<<P_bjet_mc<<endl;
      Hists[ch][0][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep, *eft_fit);
      Hists[ch][0][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep, *eft_fit);
      Hists[ch][0][2]->Fill((*selectedLeptons)[0]->phi_,weight_lep, *eft_fit);
      Hists[ch][0][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep, *eft_fit);
      Hists[ch][0][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep, *eft_fit);
      Hists[ch][0][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep, *eft_fit);
      Hists[ch][0][6]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep, *eft_fit);
      Hists[ch][0][7]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lep, *eft_fit);
      Hists[ch][0][8]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lep, *eft_fit);
      Hists[ch][0][9]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_)),weight_lep, *eft_fit);
      if(selectedJets->size()>0) Hists[ch][0][10]->Fill((*selectedJets)[0]->pt_,weight_lep, *eft_fit);
      if(selectedJets->size()>0) Hists[ch][0][11]->Fill((*selectedJets)[0]->eta_,weight_lep, *eft_fit);
      if(selectedJets->size()>0) Hists[ch][0][12]->Fill((*selectedJets)[0]->phi_,weight_lep, *eft_fit);
      Hists[ch][0][13]->Fill(selectedJets->size(),weight_lep, *eft_fit);
      Hists[ch][0][14]->Fill(nbjet,weight_lepB, *eft_fit);
      Hists[ch][0][15]->Fill(MET_T1_pt,weight_lep, *eft_fit);
      Hists[ch][0][16]->Fill(MET_T1_phi,weight_lep, *eft_fit);
      Hists[ch][0][17]->Fill(PV_npvs,weight_lep, *eft_fit);
      Hists[ch][0][18]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep, *eft_fit);
      Hists[ch][0][19]->Fill(MVAoutput,weight_lep, *eft_fit);
  
      if (DyPass){
        Hists[ch][1][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep, *eft_fit);
        Hists[ch][1][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep, *eft_fit);
        Hists[ch][1][2]->Fill((*selectedLeptons)[0]->phi_,weight_lep, *eft_fit);
        Hists[ch][1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep, *eft_fit);
        Hists[ch][1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep, *eft_fit);
        Hists[ch][1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep, *eft_fit);
        Hists[ch][1][6]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep, *eft_fit);
        Hists[ch][1][7]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lep, *eft_fit);
        Hists[ch][1][8]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lep, *eft_fit);
        Hists[ch][1][9]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_)),weight_lep, *eft_fit);
        if(selectedJets->size()>0) Hists[ch][1][10]->Fill((*selectedJets)[0]->pt_,weight_lep, *eft_fit);
        if(selectedJets->size()>0) Hists[ch][1][11]->Fill((*selectedJets)[0]->eta_,weight_lep, *eft_fit);
        if(selectedJets->size()>0) Hists[ch][1][12]->Fill((*selectedJets)[0]->phi_,weight_lep, *eft_fit);
        Hists[ch][1][13]->Fill(selectedJets->size(),weight_lep, *eft_fit);
        Hists[ch][1][14]->Fill(nbjet,weight_lepB, *eft_fit);
        Hists[ch][1][15]->Fill(MET_T1_pt,weight_lep, *eft_fit);
        Hists[ch][1][16]->Fill(MET_T1_phi,weight_lep, *eft_fit);
        Hists[ch][1][17]->Fill(PV_npvs,weight_lep, *eft_fit);
        Hists[ch][1][18]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep, *eft_fit);
        Hists[ch][1][19]->Fill(MVAoutput,weight_lep, *eft_fit);
      }
      if(nbjet==1 && DyPass){
        recoL2 = (*selectedLeptons)[1]->p4_;
        for (int l=0;l<selectedJets->size();l++){
          if((*selectedJets)[l]->btag_){
            recoBjet = (*selectedJets)[l]->p4_;
            break;
          }
        }
        recoNu = Wneutrino(MET_T1_pt, MET_T1_phi, (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
        recoW = recoNu + recoL2;
        recoTop = recoW + recoBjet;
        MVA_lep1Pt=(*selectedLeptons)[0]->pt_;
        MVA_lep2Pt=(*selectedLeptons)[1]->pt_;
        MVA_llM=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M();
        MVA_llPt=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt();
        MVA_llDr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_);
        MVA_llDphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_));
        MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
        MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
        MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
        MVA_topPt=recoTop.Pt();
        MVAoutput = readerMVA->EvaluateMVA( "BDTG");
        Hists[ch][2][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB, *eft_fit);
        Hists[ch][2][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB, *eft_fit);
        Hists[ch][2][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB, *eft_fit);
        Hists[ch][2][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB, *eft_fit);
        Hists[ch][2][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB, *eft_fit);
        Hists[ch][2][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB, *eft_fit);
        Hists[ch][2][6]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
        Hists[ch][2][7]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB, *eft_fit);
        Hists[ch][2][8]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB, *eft_fit);
        Hists[ch][2][9]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB, *eft_fit);
        Hists[ch][2][10]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
        Hists[ch][2][11]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
        Hists[ch][2][12]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
        Hists[ch][2][13]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
        Hists[ch][2][14]->Fill(nbjet,weight_lepB, *eft_fit);
        Hists[ch][2][15]->Fill(MET_T1_pt,weight_lepB, *eft_fit);
        Hists[ch][2][16]->Fill(MET_T1_phi,weight_lepB, *eft_fit);
        Hists[ch][2][17]->Fill(PV_npvs,weight_lepB, *eft_fit);
        Hists[ch][2][18]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
        Hists[ch][2][19]->Fill(MVAoutput,weight_lepB, *eft_fit);
        Hists[ch][2][20]->Fill(recoTop.M(),weight_lepB, *eft_fit);
        Hists[ch][2][21]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
        Hists[ch][2][22]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
        Hists[ch][2][23]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
        Hists[ch][2][24]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
    
        lep1Pt_=(*selectedLeptons)[0]->pt_;
        lep1Eta_=(*selectedLeptons)[0]->eta_;
        lep2Pt_=(*selectedLeptons)[1]->pt_;
        lep2Eta_=(*selectedLeptons)[1]->eta_;
        llM_=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M();
        llPt_=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt();
        llDr_=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_);
        llDphi_=abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_));
        jet1Pt_=(*selectedJets)[0]->pt_;
        jet1Eta_=(*selectedJets)[0]->eta_;
        topMass_=recoTop.M();
        topL1Dphi_=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
        topL1Dr_=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
        topL1DptOsumPt_=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
        topPt_=recoTop.Pt();
        weight_=weight_lepB;
        tree_out.Fill() ;
      }
      if(nbjet>1 && DyPass){
        Hists[ch][3][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB, *eft_fit);
        Hists[ch][3][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB, *eft_fit);
        Hists[ch][3][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB, *eft_fit);
        Hists[ch][3][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB, *eft_fit);
        Hists[ch][3][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB, *eft_fit);
        Hists[ch][3][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB, *eft_fit);
        Hists[ch][3][6]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
        Hists[ch][3][7]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB, *eft_fit);
        Hists[ch][3][8]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB, *eft_fit);
        Hists[ch][3][9]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_)),weight_lepB, *eft_fit);
        Hists[ch][3][10]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
        Hists[ch][3][11]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
        Hists[ch][3][12]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
        Hists[ch][3][13]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
        Hists[ch][3][14]->Fill(nbjet,weight_lepB, *eft_fit);
        Hists[ch][3][15]->Fill(MET_T1_pt,weight_lepB, *eft_fit);
        Hists[ch][3][16]->Fill(MET_T1_phi,weight_lepB, *eft_fit);
        Hists[ch][3][17]->Fill(PV_npvs,weight_lepB, *eft_fit);
        Hists[ch][3][18]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
        Hists[ch][3][19]->Fill(MVAoutput,weight_lepB, *eft_fit);
      }
      delete eft_fit;
 
  //sys Hists only for regions with b-tagged jets
      if (data == "mc") {
        if(nbjet==1 && DyPass){
          for (int n=0;n<sys.size();++n){
            if (n>8 && n<15) continue;
            if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* (sysUpWeights[n]/nominalWeights[n]));
            else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);    
            HistsSysUp[ch][2][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]) , *eft_fit);
            HistsSysUp[ch][2][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][10][n]->Fill((*selectedJets)[0]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][11][n]->Fill((*selectedJets)[0]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][12][n]->Fill((*selectedJets)[0]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][13][n]->Fill(selectedJets->size(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][14][n]->Fill(nbjet,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][15][n]->Fill(MET_T1_pt,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][16][n]->Fill(MET_T1_phi,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][17][n]->Fill(PV_npvs,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][19][n]->Fill(MVAoutput,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][20][n]->Fill(recoTop.M(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][21][n]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][22][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][23][n]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][2][24][n]->Fill(recoTop.Pt(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            delete eft_fit;
            if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* (sysDownWeights[n]/nominalWeights[n]));
            else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
            HistsSysDown[ch][2][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]) , *eft_fit);
            HistsSysDown[ch][2][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][10][n]->Fill((*selectedJets)[0]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][11][n]->Fill((*selectedJets)[0]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][12][n]->Fill((*selectedJets)[0]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][13][n]->Fill(selectedJets->size(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][14][n]->Fill(nbjet,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][15][n]->Fill(MET_T1_pt,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][16][n]->Fill(MET_T1_phi,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][17][n]->Fill(PV_npvs,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][19][n]->Fill(MVAoutput,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][20][n]->Fill(recoTop.M(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][21][n]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][22][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][2][23][n]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            delete eft_fit;
          }
          if ((fname.Contains("TTTo2L2Nu")&& !fname.Contains("sys")) || fname.Contains("BNV")){
            for (int n=0;n<nScale;++n){
              if(isnan(LHEScaleWeight[n]) || isinf(LHEScaleWeight[n])) continue;
              if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* LHEScaleWeight[n]);
              else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
              HistsSysReweightsQscale[ch][0][0][n]->Fill(MVAoutput,weight_lepB * LHEScaleWeight[n], *eft_fit);
              delete eft_fit;
            }
            for (int n=0;n<nPdf;++n){
              if(isnan(LHEPdfWeight[n]) || isinf(LHEPdfWeight[n])) continue;
              if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* LHEPdfWeight[n]);
              else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
              HistsSysReweightsPDF[ch][0][0][n]->Fill(MVAoutput,weight_lepB * LHEPdfWeight[n], *eft_fit);
              delete eft_fit;
            }
            for (int n=0;n<nPS;++n){
              if(isnan(PSWeight[n]) || isinf(PSWeight[n])) continue;
              if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* PSWeight[n]);
              else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
              HistsSysReweightsPS[ch][0][0][n]->Fill(MVAoutput,weight_lepB * PSWeight[n], *eft_fit);
              delete eft_fit;
            }
          }
          //Unclustered MET unc.
          if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
          else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
          recoL2 = (*selectedLeptons)[1]->p4_;
          for (int l=0;l<selectedJets->size();l++){
            if((*selectedJets)[l]->btag_){
              recoBjet = (*selectedJets)[l]->p4_;
              break;
            }
          }
          recoNu = Wneutrino(metUnclusMETUp, metUnclusMETPhiUp, (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
          recoW = recoNu + recoL2;
          recoTop = recoW + recoBjet;
          MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
          MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
          MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
          MVA_topPt=recoTop.Pt();
          MVAoutput = readerMVA->EvaluateMVA( "BDTG");
          HistsSysUp[ch][2][0][14]->Fill((*selectedLeptons)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][1][14]->Fill((*selectedLeptons)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][2][14]->Fill((*selectedLeptons)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][3][14]->Fill((*selectedLeptons)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][4][14]->Fill((*selectedLeptons)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][5][14]->Fill((*selectedLeptons)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][6][14]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][7][14]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][8][14]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][9][14]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][10][14]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][11][14]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][12][14]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][13][14]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][14][14]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][15][14]->Fill(metUnclusMETUp,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][16][14]->Fill(metUnclusMETPhiUp,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][17][14]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][18][14]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][19][14]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][20][14]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][21][14]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][22][14]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][23][14]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][24][14]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);

          recoNu = Wneutrino(metUnclusMETDown, metUnclusMETPhiDown, (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
          recoW = recoNu + recoL2;
          recoTop = recoW + recoBjet;
          MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
          MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
          MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
          MVA_topPt=recoTop.Pt();
          MVAoutput = readerMVA->EvaluateMVA( "BDTG");
          HistsSysDown[ch][2][0][14]->Fill((*selectedLeptons)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][1][14]->Fill((*selectedLeptons)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][2][14]->Fill((*selectedLeptons)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][3][14]->Fill((*selectedLeptons)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][4][14]->Fill((*selectedLeptons)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][5][14]->Fill((*selectedLeptons)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][6][14]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][7][14]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][8][14]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][9][14]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][10][14]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][11][14]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][12][14]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][13][14]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][14][14]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][15][14]->Fill(metUnclusMETDown,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][16][14]->Fill(metUnclusMETPhiDown,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][17][14]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][18][14]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][19][14]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][20][14]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][21][14]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][22][14]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][23][14]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][24][14]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
          delete eft_fit;
        }
        if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
        else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);

        if(nbjetJesUp ==1 && DyPass) {
          recoL2 = (*selectedLeptons)[1]->p4_;
          MVA_lep1Pt=(*selectedLeptons)[0]->pt_;
          MVA_lep2Pt=(*selectedLeptons)[1]->pt_;
          MVA_llM=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M();
          MVA_llPt=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt();
          MVA_llDr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_);
          MVA_llDphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_));
          recoNu = Wneutrino(MetJetsPtJesUp, MetJetsPhiJesUp, (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
          recoW = recoNu + recoL2;
          for (int l=0;l<selectedJetsJesUp->size();l++){
            if((*selectedJetsJesUp)[l]->btag_){
              recoBjet = (*selectedJetsJesUp)[l]->p4_;
              break;
            }
          }
          recoTop = recoW + recoBjet;
          MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
          MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
          MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
          MVA_topPt=recoTop.Pt();
          MVAoutput = readerMVA->EvaluateMVA( "BDTG");
          HistsSysUp[ch][2][0][9]->Fill((*selectedLeptons)[0]->pt_,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][1][9]->Fill((*selectedLeptons)[0]->eta_,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][2][9]->Fill((*selectedLeptons)[0]->phi_,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][3][9]->Fill((*selectedLeptons)[1]->pt_,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][4][9]->Fill((*selectedLeptons)[1]->eta_,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][5][9]->Fill((*selectedLeptons)[1]->phi_,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][6][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][7][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][8][9]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][9][9]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][10][9]->Fill((*selectedJetsJesUp)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][11][9]->Fill((*selectedJetsJesUp)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][12][9]->Fill((*selectedJetsJesUp)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][13][9]->Fill(selectedJetsJesUp->size(),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][14][9]->Fill(nbjetJesUp,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][15][9]->Fill(MetJetsPtJesUp,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][16][9]->Fill(MetJetsPhiJesUp,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][17][9]->Fill(PV_npvs,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][18][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][19][9]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][20][9]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][21][9]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][22][9]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][23][9]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][24][9]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        if(nbjetJesDown==1 && DyPass){
          recoL2 = (*selectedLeptons)[1]->p4_;
          MVA_lep1Pt=(*selectedLeptons)[0]->pt_;
          MVA_lep2Pt=(*selectedLeptons)[1]->pt_;
          MVA_llM=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M();
          MVA_llPt=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt();
          MVA_llDr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_);
          MVA_llDphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_)); 
          recoNu = Wneutrino(MetJetsPtJesDown, MetJetsPhiJesDown, (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
          recoW = recoNu + recoL2;
          for (int l=0;l<selectedJetsJesDown->size();l++){
            if((*selectedJetsJesDown)[l]->btag_){
              recoBjet = (*selectedJetsJesDown)[l]->p4_;
              break;
            }
          }
          recoTop = recoW + recoBjet;
          MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
          MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
          MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
          MVA_topPt=recoTop.Pt();
          MVAoutput = readerMVA->EvaluateMVA( "BDTG"); 
          HistsSysDown[ch][2][0][9]->Fill((*selectedLeptons)[0]->pt_,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][1][9]->Fill((*selectedLeptons)[0]->eta_,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][2][9]->Fill((*selectedLeptons)[0]->phi_,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][3][9]->Fill((*selectedLeptons)[1]->pt_,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][4][9]->Fill((*selectedLeptons)[1]->eta_,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][5][9]->Fill((*selectedLeptons)[1]->phi_,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][6][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][7][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][8][9]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][9][9]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][10][9]->Fill((*selectedJetsJesDown)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][11][9]->Fill((*selectedJetsJesDown)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][12][9]->Fill((*selectedJetsJesDown)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][13][9]->Fill(selectedJetsJesDown->size(),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][14][9]->Fill(nbjetJesDown,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][15][9]->Fill(MetJetsPtJesDown,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][16][9]->Fill(MetJetsPhiJesDown,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][17][9]->Fill(PV_npvs,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][18][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][19][9]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][20][9]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][21][9]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][22][9]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][23][9]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][24][9]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        if(nbjetJerUp ==1 && DyPass) {
          recoL2 = (*selectedLeptons)[1]->p4_;
          MVA_lep1Pt=(*selectedLeptons)[0]->pt_;
          MVA_lep2Pt=(*selectedLeptons)[1]->pt_;
          MVA_llM=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M();
          MVA_llPt=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt();
          MVA_llDr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_);
          MVA_llDphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_));
          recoNu = Wneutrino(MetJetsPtJerUp, MetJetsPhiJerUp, (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
          recoW = recoNu + recoL2;
          for (int l=0;l<selectedJetsJerUp->size();l++){
            if((*selectedJetsJerUp)[l]->btag_){
              recoBjet = (*selectedJetsJerUp)[l]->p4_;
              break;
            }
          }
          recoTop = recoW + recoBjet;
          MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
          MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
          MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
          MVA_topPt=recoTop.Pt();
          MVAoutput = readerMVA->EvaluateMVA( "BDTG");
          HistsSysUp[ch][2][0][10]->Fill((*selectedLeptons)[0]->pt_,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][1][10]->Fill((*selectedLeptons)[0]->eta_,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][2][10]->Fill((*selectedLeptons)[0]->phi_,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][3][10]->Fill((*selectedLeptons)[1]->pt_,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][4][10]->Fill((*selectedLeptons)[1]->eta_,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][5][10]->Fill((*selectedLeptons)[1]->phi_,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][6][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][7][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][8][10]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][9][10]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][10][10]->Fill((*selectedJetsJerUp)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][11][10]->Fill((*selectedJetsJerUp)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][12][10]->Fill((*selectedJetsJerUp)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][13][10]->Fill(selectedJetsJerUp->size(),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][14][10]->Fill(nbjetJerUp,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][15][10]->Fill(MetJetsPtJerUp,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][16][10]->Fill(MetJetsPhiJerUp,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][17][10]->Fill(PV_npvs,weight_lepB , *eft_fit);
          HistsSysUp[ch][2][18][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][19][10]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][20][10]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][21][10]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][22][10]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][23][10]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB , *eft_fit);
          HistsSysUp[ch][2][24][10]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        if(nbjetJerDown==1 && DyPass){
          recoL2 = (*selectedLeptons)[1]->p4_;
          MVA_lep1Pt=(*selectedLeptons)[0]->pt_;
          MVA_lep2Pt=(*selectedLeptons)[1]->pt_;
          MVA_llM=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M();
          MVA_llPt=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt();
          MVA_llDr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_);
          MVA_llDphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_));
          recoNu = Wneutrino(MetJetsPtJerDown, MetJetsPhiJerDown, (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
          recoW = recoNu + recoL2;
          for (int l=0;l<selectedJetsJerDown->size();l++){
            if((*selectedJetsJerDown)[l]->btag_){
              recoBjet = (*selectedJetsJerDown)[l]->p4_;
              break;
            }
          }
          recoTop = recoW + recoBjet;
          MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
          MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
          MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
          MVA_topPt=recoTop.Pt();
          MVAoutput = readerMVA->EvaluateMVA( "BDTG");
          HistsSysDown[ch][2][0][10]->Fill((*selectedLeptons)[0]->pt_,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][1][10]->Fill((*selectedLeptons)[0]->eta_,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][2][10]->Fill((*selectedLeptons)[0]->phi_,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][3][10]->Fill((*selectedLeptons)[1]->pt_,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][4][10]->Fill((*selectedLeptons)[1]->eta_,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][5][10]->Fill((*selectedLeptons)[1]->phi_,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][6][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][7][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][8][10]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][9][10]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][10][10]->Fill((*selectedJetsJerDown)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][11][10]->Fill((*selectedJetsJerDown)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][12][10]->Fill((*selectedJetsJerDown)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][13][10]->Fill(selectedJetsJerDown->size(),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][14][10]->Fill(nbjetJerDown,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][15][10]->Fill(MetJetsPtJerDown,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][16][10]->Fill(MetJetsPhiJerDown,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][17][10]->Fill(PV_npvs,weight_lepB , *eft_fit);
          HistsSysDown[ch][2][18][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][19][10]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][20][10]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][21][10]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][22][10]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][23][10]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB , *eft_fit);
          HistsSysDown[ch][2][24][10]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
    
        for (int n=0;n<sysJecNames.size();++n){
          recoL2 = (*selectedLeptons)[1]->p4_;
          MVA_lep1Pt=(*selectedLeptons)[0]->pt_;
          MVA_lep2Pt=(*selectedLeptons)[1]->pt_;
          MVA_llM=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M();
          MVA_llPt=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt();
          MVA_llDr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_);
          MVA_llDphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_));
          recoNu = Wneutrino((*JECsysMETUp)[n], (*JECsysMETPhiUp)[n], (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
          recoW = recoNu + recoL2;
          if ((*JECsysNbtagUp)[n]==1 && DyPass) {
            for (int l=0;l<(*JECsysUp)[n].size();l++){
              if((*JECsysUp)[n][l]->btag_){
                recoBjet = (*JECsysUp)[n][l]->p4_;
                break;
              }
            }
            recoTop = recoW + recoBjet;
            MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
            MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
            MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
            MVA_topPt=recoTop.Pt();
            MVAoutput = readerMVA->EvaluateMVA( "BDTG");
            HistsJECUp[ch][0][0][n]->Fill(MVAoutput,weight_lepB, *eft_fit);
          }

          if ((*JECsysNbtagDown)[n]==1 && DyPass) {
          recoL2 = (*selectedLeptons)[1]->p4_;
          MVA_lep1Pt=(*selectedLeptons)[0]->pt_;
          MVA_lep2Pt=(*selectedLeptons)[1]->pt_;
          MVA_llM=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M();
          MVA_llPt=((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt();
          MVA_llDr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_);
          MVA_llDphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_));
            recoNu = Wneutrino((*JECsysMETDown)[n], (*JECsysMETPhiDown)[n], (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
            recoW = recoNu + recoL2;
            for (int l=0;l<(*JECsysDown)[n].size();l++){
              if((*JECsysDown)[n][l]->btag_){
                recoBjet = (*JECsysDown)[n][l]->p4_;
                break;
              }
            }
            recoTop = recoW + recoBjet;
            MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
            MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
            MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
            MVA_topPt=recoTop.Pt();
            MVAoutput = readerMVA->EvaluateMVA( "BDTG");
            HistsJECDown[ch][0][0][n]->Fill(MVAoutput,weight_lepB, *eft_fit);
          }
        }
        delete eft_fit;    
   
//New Region 
        if(nbjet>1 && DyPass){
          for (int n=0;n<sys.size();++n){
            if (n>8 && n<15) continue;
            if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* (sysUpWeights[n]/nominalWeights[n]));
            else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);    
            HistsSysUp[ch][3][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]) , *eft_fit);
            HistsSysUp[ch][3][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][10][n]->Fill((*selectedJets)[0]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][11][n]->Fill((*selectedJets)[0]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][12][n]->Fill((*selectedJets)[0]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][13][n]->Fill(selectedJets->size(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][14][n]->Fill(nbjet,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][15][n]->Fill(MET_T1_pt,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][16][n]->Fill(MET_T1_phi,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][17][n]->Fill(PV_npvs,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][19][n]->Fill(MVAoutput,weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][20][n]->Fill(recoTop.M(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][21][n]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][22][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][23][n]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysUp[ch][3][24][n]->Fill(recoTop.Pt(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]), *eft_fit);
            delete eft_fit;
            if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* (sysDownWeights[n]/nominalWeights[n]));
            else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
            HistsSysDown[ch][3][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]) , *eft_fit);
            HistsSysDown[ch][3][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][10][n]->Fill((*selectedJets)[0]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][11][n]->Fill((*selectedJets)[0]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][12][n]->Fill((*selectedJets)[0]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][13][n]->Fill(selectedJets->size(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][14][n]->Fill(nbjet,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][15][n]->Fill(MET_T1_pt,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][16][n]->Fill(MET_T1_phi,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][17][n]->Fill(PV_npvs,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][19][n]->Fill(MVAoutput,weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][20][n]->Fill(recoTop.M(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][21][n]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][22][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            HistsSysDown[ch][3][23][n]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB * (sysDownWeights[n]/nominalWeights[n]), *eft_fit);
            delete eft_fit;
          }
          if ((fname.Contains("TTTo2L2Nu")&& !fname.Contains("sys")) || fname.Contains("BNV")){
            for (int n=0;n<nScale;++n){
              if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* LHEScaleWeight[n]);
              else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
              HistsSysReweightsQscale[ch][1][0][n]->Fill(MVAoutput,weight_lepB * LHEScaleWeight[n], *eft_fit);
              delete eft_fit;
            }
            for (int n=0;n<nPdf;++n){
              if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* LHEPdfWeight[n]);
              else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
              HistsSysReweightsPDF[ch][1][0][n]->Fill(MVAoutput,weight_lepB * LHEPdfWeight[n], *eft_fit);
              delete eft_fit;
            }
            for (int n=0;n<nPS;++n){
              if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* PSWeight[n]);
              else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
              HistsSysReweightsPS[ch][1][0][n]->Fill(MVAoutput,weight_lepB * PSWeight[n], *eft_fit);
              delete eft_fit;
            }
          }
          // Unclustered MET
          if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
          else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
          HistsSysUp[ch][3][0][14]->Fill((*selectedLeptons)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][1][14]->Fill((*selectedLeptons)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][2][14]->Fill((*selectedLeptons)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][3][14]->Fill((*selectedLeptons)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][4][14]->Fill((*selectedLeptons)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][5][14]->Fill((*selectedLeptons)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][6][14]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][7][14]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][8][14]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][9][14]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][10][14]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][11][14]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][12][14]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][13][14]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][14][14]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][15][14]->Fill(metUnclusMETUp,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][16][14]->Fill(metUnclusMETPhiUp,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][17][14]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][18][14]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][19][14]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][20][14]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][21][14]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][22][14]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][23][14]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][24][14]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);

          HistsSysDown[ch][3][0][14]->Fill((*selectedLeptons)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][1][14]->Fill((*selectedLeptons)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][2][14]->Fill((*selectedLeptons)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][3][14]->Fill((*selectedLeptons)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][4][14]->Fill((*selectedLeptons)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][5][14]->Fill((*selectedLeptons)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][6][14]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][7][14]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][8][14]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][9][14]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][10][14]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][11][14]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][12][14]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][13][14]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][14][14]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][15][14]->Fill(metUnclusMETDown,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][16][14]->Fill(metUnclusMETPhiDown,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][17][14]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][18][14]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][19][14]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][20][14]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][21][14]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][22][14]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][23][14]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][24][14]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
          delete eft_fit;
        }
 
    
        if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
        else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
//Jes unc
        if(nbjetJesUp >1 && DyPass) {
          HistsSysUp[ch][3][0][9]->Fill((*selectedLeptons)[0]->pt_,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][1][9]->Fill((*selectedLeptons)[0]->eta_,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][2][9]->Fill((*selectedLeptons)[0]->phi_,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][3][9]->Fill((*selectedLeptons)[1]->pt_,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][4][9]->Fill((*selectedLeptons)[1]->eta_,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][5][9]->Fill((*selectedLeptons)[1]->phi_,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][6][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][7][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][8][9]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][9][9]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][10][9]->Fill((*selectedJetsJesUp)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][11][9]->Fill((*selectedJetsJesUp)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][12][9]->Fill((*selectedJetsJesUp)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][13][9]->Fill(selectedJetsJesUp->size(),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][14][9]->Fill(nbjetJesUp,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][15][9]->Fill(MetJetsPtJesUp,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][16][9]->Fill(MetJetsPhiJesUp,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][17][9]->Fill(PV_npvs,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][18][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][19][9]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][20][9]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][21][9]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][22][9]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][23][9]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][24][9]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        if(nbjetJesDown>1 && DyPass){ 
          HistsSysDown[ch][3][0][9]->Fill((*selectedLeptons)[0]->pt_,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][1][9]->Fill((*selectedLeptons)[0]->eta_,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][2][9]->Fill((*selectedLeptons)[0]->phi_,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][3][9]->Fill((*selectedLeptons)[1]->pt_,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][4][9]->Fill((*selectedLeptons)[1]->eta_,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][5][9]->Fill((*selectedLeptons)[1]->phi_,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][6][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][7][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][8][9]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][9][9]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][10][9]->Fill((*selectedJetsJesDown)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][11][9]->Fill((*selectedJetsJesDown)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][12][9]->Fill((*selectedJetsJesDown)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][13][9]->Fill(selectedJetsJesDown->size(),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][14][9]->Fill(nbjetJesDown,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][15][9]->Fill(MetJetsPtJesDown,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][16][9]->Fill(MetJetsPhiJesDown,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][17][9]->Fill(PV_npvs,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][18][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][19][9]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][20][9]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][21][9]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][22][9]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][23][9]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][24][9]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        if(nbjetJerUp >1 && DyPass) {
          HistsSysUp[ch][3][0][10]->Fill((*selectedLeptons)[0]->pt_,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][1][10]->Fill((*selectedLeptons)[0]->eta_,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][2][10]->Fill((*selectedLeptons)[0]->phi_,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][3][10]->Fill((*selectedLeptons)[1]->pt_,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][4][10]->Fill((*selectedLeptons)[1]->eta_,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][5][10]->Fill((*selectedLeptons)[1]->phi_,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][6][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][7][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][8][10]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][9][10]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][10][10]->Fill((*selectedJetsJerUp)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][11][10]->Fill((*selectedJetsJerUp)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][12][10]->Fill((*selectedJetsJerUp)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][13][10]->Fill(selectedJetsJerUp->size(),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][14][10]->Fill(nbjetJerUp,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][15][10]->Fill(MetJetsPtJerUp,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][16][10]->Fill(MetJetsPhiJerUp,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][17][10]->Fill(PV_npvs,weight_lepB , *eft_fit);
          HistsSysUp[ch][3][18][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][19][10]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][20][10]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][21][10]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][22][10]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][23][10]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB , *eft_fit);
          HistsSysUp[ch][3][24][10]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        if(nbjetJerDown>1 && DyPass){
          HistsSysDown[ch][3][0][10]->Fill((*selectedLeptons)[0]->pt_,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][1][10]->Fill((*selectedLeptons)[0]->eta_,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][2][10]->Fill((*selectedLeptons)[0]->phi_,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][3][10]->Fill((*selectedLeptons)[1]->pt_,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][4][10]->Fill((*selectedLeptons)[1]->eta_,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][5][10]->Fill((*selectedLeptons)[1]->phi_,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][6][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][7][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][8][10]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][9][10]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][10][10]->Fill((*selectedJetsJerDown)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][11][10]->Fill((*selectedJetsJerDown)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][12][10]->Fill((*selectedJetsJerDown)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][13][10]->Fill(selectedJetsJerDown->size(),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][14][10]->Fill(nbjetJerDown,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][15][10]->Fill(MetJetsPtJerDown,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][16][10]->Fill(MetJetsPhiJerDown,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][17][10]->Fill(PV_npvs,weight_lepB , *eft_fit);
          HistsSysDown[ch][3][18][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][19][10]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][20][10]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][21][10]->Fill(abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][22][10]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][23][10]->Fill((abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()),weight_lepB , *eft_fit);
          HistsSysDown[ch][3][24][10]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
    
        for (int n=0;n<sysJecNames.size();++n){
          if ((*JECsysNbtagUp)[n]>1 && DyPass) {
            HistsJECUp[ch][1][0][n]->Fill(MVAoutput,weight_lepB, *eft_fit);
          }
          if ((*JECsysNbtagDown)[n]>1 && DyPass) {
            HistsJECDown[ch][1][0][n]->Fill(MVAoutput,weight_lepB, *eft_fit);
          }
        }
        delete eft_fit;
      }
    }

// Lepton scale resolution uncertainties
    if (data == "mc"){
      leptonPass=false;
      triggerPass=false;
      DyPass = true;
      ch=10;
      if(selectedLeptonsMuScaleUp->size() ==2){
        if (((*selectedLeptonsMuScaleUp)[0]->pt_ > 25) &&
        ((*selectedLeptonsMuScaleUp)[0]->charge_ * (*selectedLeptonsMuScaleUp)[1]->charge_ == -1) &&
        ((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).M()>20) leptonPass=true; 
        if ((*selectedLeptonsMuScaleUp)[0]->lep_ + (*selectedLeptonsMuScaleUp)[1]->lep_ == 2) ch = 0;
        if ((*selectedLeptonsMuScaleUp)[0]->lep_ + (*selectedLeptonsMuScaleUp)[1]->lep_ == 11) ch = 1;
        if ((*selectedLeptonsMuScaleUp)[0]->lep_ + (*selectedLeptonsMuScaleUp)[1]->lep_ == 20) ch = 2;
        if(ch ==0 && triggerPassEE) triggerPass=true;
        if(ch ==1 && triggerPassEMu) triggerPass=true;
        if(ch ==2 && triggerPassMuMu) triggerPass=true;
        if ((*selectedLeptonsMuScaleUp)[0]->lep_ + (*selectedLeptonsMuScaleUp)[1]->lep_ != 11 && ((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).M()<106 ) DyPass = false;
        if (data == "mc" && ch==0) sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptonsMuScaleUp)[0]->pt_, (*selectedLeptonsMuScaleUp)[1]->pt_,"central",false, true);
        if (data == "mc" && ch==1) sf_Trigger = scale_factor(&sf_triggeremu_H, (*selectedLeptonsMuScaleUp)[0]->pt_, (*selectedLeptonsMuScaleUp)[1]->pt_,"central",false, true);
        if (data == "mc" && ch==2) sf_Trigger = scale_factor(&sf_triggermumu_H, (*selectedLeptonsMuScaleUp)[0]->pt_, (*selectedLeptonsMuScaleUp)[1]->pt_,"central",false, true);
      }
      if(leptonPass && triggerPass){
        weight_lepB = weight_Lumi * signnum_typical(LHEWeight_originalXWGTUP) * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowheg * ttKFactor * sf_JetPuId;
        weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowheg * ttKFactor * sf_JetPuId;
        if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
        else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
        if(nbjet==1 && DyPass){
          recoL2 = (*selectedLeptonsMuScaleUp)[1]->p4_;
          for (int l=0;l<selectedJets->size();l++){
            if((*selectedJets)[l]->btag_){
              recoBjet = (*selectedJets)[l]->p4_;
              break;
            }
          }
          recoNu = Wneutrino(MET_T1_pt, MET_T1_phi, (*selectedLeptonsMuScaleUp)[1]->pt_, (*selectedLeptonsMuScaleUp)[1]->eta_, (*selectedLeptonsMuScaleUp)[1]->phi_);
          recoW = recoNu + recoL2;
          recoTop = recoW + recoBjet;
          MVA_lep1Pt=(*selectedLeptonsMuScaleUp)[0]->pt_;
          MVA_lep2Pt=(*selectedLeptonsMuScaleUp)[1]->pt_;
          MVA_llM=((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).M();
          MVA_llPt=((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).Pt();
          MVA_llDr=deltaR((*selectedLeptonsMuScaleUp)[0]->eta_,(*selectedLeptonsMuScaleUp)[0]->phi_,(*selectedLeptonsMuScaleUp)[1]->eta_,(*selectedLeptonsMuScaleUp)[1]->phi_);
          MVA_llDphi=abs(deltaPhi((*selectedLeptonsMuScaleUp)[0]->phi_,(*selectedLeptonsMuScaleUp)[1]->phi_));
          MVA_topL1Dphi=abs(deltaPhi((*selectedLeptonsMuScaleUp)[0]->phi_,recoTop.Phi()));
          MVA_topL1Dr=deltaR((*selectedLeptonsMuScaleUp)[0]->eta_,(*selectedLeptonsMuScaleUp)[0]->phi_,recoTop.Eta(),recoTop.Phi());
          MVA_topL1DptOsumPt=(abs((*selectedLeptonsMuScaleUp)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuScaleUp)[0]->pt_ + recoTop.Pt());
          MVA_topPt=recoTop.Pt();
          MVAoutput = readerMVA->EvaluateMVA( "BDTG");
  
          HistsSysUp[ch][2][0][11]->Fill((*selectedLeptonsMuScaleUp)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][1][11]->Fill((*selectedLeptonsMuScaleUp)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][2][11]->Fill((*selectedLeptonsMuScaleUp)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][3][11]->Fill((*selectedLeptonsMuScaleUp)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][4][11]->Fill((*selectedLeptonsMuScaleUp)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][5][11]->Fill((*selectedLeptonsMuScaleUp)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][6][11]->Fill(((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][7][11]->Fill(((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][8][11]->Fill(deltaR((*selectedLeptonsMuScaleUp)[0]->eta_,(*selectedLeptonsMuScaleUp)[0]->phi_,(*selectedLeptonsMuScaleUp)[1]->eta_,(*selectedLeptonsMuScaleUp)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][9][11]->Fill(deltaPhi((*selectedLeptonsMuScaleUp)[0]->phi_,(*selectedLeptonsMuScaleUp)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][10][11]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][11][11]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][12][11]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][13][11]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][14][11]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][15][11]->Fill(MET_T1_pt,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][16][11]->Fill(MET_T1_phi,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][17][11]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][18][11]->Fill(((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][19][11]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][20][11]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][21][11]->Fill(abs(deltaPhi((*selectedLeptonsMuScaleUp)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][22][11]->Fill(deltaR((*selectedLeptonsMuScaleUp)[0]->eta_,(*selectedLeptonsMuScaleUp)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][23][11]->Fill((abs((*selectedLeptonsMuScaleUp)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuScaleUp)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][24][11]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        if(nbjet>1 && DyPass){
          HistsSysUp[ch][3][0][11]->Fill((*selectedLeptonsMuScaleUp)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][1][11]->Fill((*selectedLeptonsMuScaleUp)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][2][11]->Fill((*selectedLeptonsMuScaleUp)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][3][11]->Fill((*selectedLeptonsMuScaleUp)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][4][11]->Fill((*selectedLeptonsMuScaleUp)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][5][11]->Fill((*selectedLeptonsMuScaleUp)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][6][11]->Fill(((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][7][11]->Fill(((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][8][11]->Fill(deltaR((*selectedLeptonsMuScaleUp)[0]->eta_,(*selectedLeptonsMuScaleUp)[0]->phi_,(*selectedLeptonsMuScaleUp)[1]->eta_,(*selectedLeptonsMuScaleUp)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][9][11]->Fill(deltaPhi((*selectedLeptonsMuScaleUp)[0]->phi_,(*selectedLeptonsMuScaleUp)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][10][11]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][11][11]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][12][11]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][13][11]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][14][11]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][15][11]->Fill(MET_T1_pt,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][16][11]->Fill(MET_T1_phi,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][17][11]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][18][11]->Fill(((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][19][11]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][20][11]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][21][11]->Fill(abs(deltaPhi((*selectedLeptonsMuScaleUp)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][22][11]->Fill(deltaR((*selectedLeptonsMuScaleUp)[0]->eta_,(*selectedLeptonsMuScaleUp)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][23][11]->Fill((abs((*selectedLeptonsMuScaleUp)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuScaleUp)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][24][11]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        delete eft_fit;
      } 

      leptonPass=false;
      triggerPass=false;
      DyPass = true;
      ch=10;
      if(selectedLeptonsMuScaleDown->size() ==2){
        if (((*selectedLeptonsMuScaleDown)[0]->pt_ > 25) &&
        ((*selectedLeptonsMuScaleDown)[0]->charge_ * (*selectedLeptonsMuScaleDown)[1]->charge_ == -1) &&
        ((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).M()>20) leptonPass=true; 
        if ((*selectedLeptonsMuScaleDown)[0]->lep_ + (*selectedLeptonsMuScaleDown)[1]->lep_ == 2) ch = 0;
        if ((*selectedLeptonsMuScaleDown)[0]->lep_ + (*selectedLeptonsMuScaleDown)[1]->lep_ == 11) ch = 1;
        if ((*selectedLeptonsMuScaleDown)[0]->lep_ + (*selectedLeptonsMuScaleDown)[1]->lep_ == 20) ch = 2;
        if(ch ==0 && triggerPassEE) triggerPass=true;
        if(ch ==1 && triggerPassEMu) triggerPass=true;
        if(ch ==2 && triggerPassMuMu) triggerPass=true;
        if ((*selectedLeptonsMuScaleDown)[0]->lep_ + (*selectedLeptonsMuScaleDown)[1]->lep_ != 11 && ((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).M()<106 ) DyPass = false;
        if (data == "mc" && ch==0) sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptonsMuScaleDown)[0]->pt_, (*selectedLeptonsMuScaleDown)[1]->pt_,"central",false, true);
        if (data == "mc" && ch==1) sf_Trigger = scale_factor(&sf_triggeremu_H, (*selectedLeptonsMuScaleDown)[0]->pt_, (*selectedLeptonsMuScaleDown)[1]->pt_,"central",false, true);
        if (data == "mc" && ch==2) sf_Trigger = scale_factor(&sf_triggermumu_H, (*selectedLeptonsMuScaleDown)[0]->pt_, (*selectedLeptonsMuScaleDown)[1]->pt_,"central",false, true);
      }
      if(leptonPass && triggerPass){
        weight_lepB = weight_Lumi * signnum_typical(LHEWeight_originalXWGTUP) * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowheg * ttKFactor * sf_JetPuId;
        weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowheg * ttKFactor * sf_JetPuId;
        if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
        else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
        if(nbjet==1 && DyPass){
          recoL2 = (*selectedLeptonsMuScaleDown)[1]->p4_;
          for (int l=0;l<selectedJets->size();l++){
            if((*selectedJets)[l]->btag_){
              recoBjet = (*selectedJets)[l]->p4_;
              break;
            }
          }
          recoNu = Wneutrino(MET_T1_pt, MET_T1_phi, (*selectedLeptonsMuScaleDown)[1]->pt_, (*selectedLeptonsMuScaleDown)[1]->eta_, (*selectedLeptonsMuScaleDown)[1]->phi_);
          recoW = recoNu + recoL2;
          recoTop = recoW + recoBjet;
          MVA_lep1Pt=(*selectedLeptonsMuScaleDown)[0]->pt_;
          MVA_lep2Pt=(*selectedLeptonsMuScaleDown)[1]->pt_;
          MVA_llM=((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).M();
          MVA_llPt=((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).Pt();
          MVA_llDr=deltaR((*selectedLeptonsMuScaleDown)[0]->eta_,(*selectedLeptonsMuScaleDown)[0]->phi_,(*selectedLeptonsMuScaleDown)[1]->eta_,(*selectedLeptonsMuScaleDown)[1]->phi_);
          MVA_llDphi=abs(deltaPhi((*selectedLeptonsMuScaleDown)[0]->phi_,(*selectedLeptonsMuScaleDown)[1]->phi_));
          MVA_topL1Dphi=abs(deltaPhi((*selectedLeptonsMuScaleDown)[0]->phi_,recoTop.Phi()));
          MVA_topL1Dr=deltaR((*selectedLeptonsMuScaleDown)[0]->eta_,(*selectedLeptonsMuScaleDown)[0]->phi_,recoTop.Eta(),recoTop.Phi());
          MVA_topL1DptOsumPt=(abs((*selectedLeptonsMuScaleDown)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuScaleDown)[0]->pt_ + recoTop.Pt());
          MVA_topPt=recoTop.Pt();
          MVAoutput = readerMVA->EvaluateMVA( "BDTG");
  
          HistsSysDown[ch][2][0][11]->Fill((*selectedLeptonsMuScaleDown)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][1][11]->Fill((*selectedLeptonsMuScaleDown)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][2][11]->Fill((*selectedLeptonsMuScaleDown)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][3][11]->Fill((*selectedLeptonsMuScaleDown)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][4][11]->Fill((*selectedLeptonsMuScaleDown)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][5][11]->Fill((*selectedLeptonsMuScaleDown)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][6][11]->Fill(((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][7][11]->Fill(((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][8][11]->Fill(deltaR((*selectedLeptonsMuScaleDown)[0]->eta_,(*selectedLeptonsMuScaleDown)[0]->phi_,(*selectedLeptonsMuScaleDown)[1]->eta_,(*selectedLeptonsMuScaleDown)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][9][11]->Fill(deltaPhi((*selectedLeptonsMuScaleDown)[0]->phi_,(*selectedLeptonsMuScaleDown)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][10][11]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][11][11]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][12][11]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][13][11]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][14][11]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][15][11]->Fill(MET_T1_pt,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][16][11]->Fill(MET_T1_phi,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][17][11]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][18][11]->Fill(((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][19][11]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][20][11]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][21][11]->Fill(abs(deltaPhi((*selectedLeptonsMuScaleDown)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][22][11]->Fill(deltaR((*selectedLeptonsMuScaleDown)[0]->eta_,(*selectedLeptonsMuScaleDown)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][23][11]->Fill((abs((*selectedLeptonsMuScaleDown)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuScaleDown)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][24][11]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        if(nbjet>1 && DyPass){
          HistsSysDown[ch][3][0][11]->Fill((*selectedLeptonsMuScaleDown)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][1][11]->Fill((*selectedLeptonsMuScaleDown)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][2][11]->Fill((*selectedLeptonsMuScaleDown)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][3][11]->Fill((*selectedLeptonsMuScaleDown)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][4][11]->Fill((*selectedLeptonsMuScaleDown)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][5][11]->Fill((*selectedLeptonsMuScaleDown)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][6][11]->Fill(((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][7][11]->Fill(((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][8][11]->Fill(deltaR((*selectedLeptonsMuScaleDown)[0]->eta_,(*selectedLeptonsMuScaleDown)[0]->phi_,(*selectedLeptonsMuScaleDown)[1]->eta_,(*selectedLeptonsMuScaleDown)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][9][11]->Fill(deltaPhi((*selectedLeptonsMuScaleDown)[0]->phi_,(*selectedLeptonsMuScaleDown)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][10][11]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][11][11]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][12][11]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][13][11]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][14][11]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][15][11]->Fill(MET_T1_pt,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][16][11]->Fill(MET_T1_phi,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][17][11]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][18][11]->Fill(((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][19][11]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][20][11]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][21][11]->Fill(abs(deltaPhi((*selectedLeptonsMuScaleDown)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][22][11]->Fill(deltaR((*selectedLeptonsMuScaleDown)[0]->eta_,(*selectedLeptonsMuScaleDown)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][23][11]->Fill((abs((*selectedLeptonsMuScaleDown)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuScaleDown)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][24][11]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        delete eft_fit;
      } 

      leptonPass=false;
      triggerPass=false;
      DyPass = true;
      ch=10;
      if(selectedLeptonsEleScaleUp->size() ==2){
        if (((*selectedLeptonsEleScaleUp)[0]->pt_ > 25) &&
        ((*selectedLeptonsEleScaleUp)[0]->charge_ * (*selectedLeptonsEleScaleUp)[1]->charge_ == -1) &&
        ((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).M()>20) leptonPass=true; 
        if ((*selectedLeptonsEleScaleUp)[0]->lep_ + (*selectedLeptonsEleScaleUp)[1]->lep_ == 2) ch = 0;
        if ((*selectedLeptonsEleScaleUp)[0]->lep_ + (*selectedLeptonsEleScaleUp)[1]->lep_ == 11) ch = 1;
        if ((*selectedLeptonsEleScaleUp)[0]->lep_ + (*selectedLeptonsEleScaleUp)[1]->lep_ == 20) ch = 2;
        if(ch ==0 && triggerPassEE) triggerPass=true;
        if(ch ==1 && triggerPassEMu) triggerPass=true;
        if(ch ==2 && triggerPassMuMu) triggerPass=true;
        if ((*selectedLeptonsEleScaleUp)[0]->lep_ + (*selectedLeptonsEleScaleUp)[1]->lep_ != 11 && ((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).M()<106 ) DyPass = false;
        if (data == "mc" && ch==0) sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptonsEleScaleUp)[0]->pt_, (*selectedLeptonsEleScaleUp)[1]->pt_,"central",false, true);
        if (data == "mc" && ch==1) sf_Trigger = scale_factor(&sf_triggeremu_H, (*selectedLeptonsEleScaleUp)[0]->pt_, (*selectedLeptonsEleScaleUp)[1]->pt_,"central",false, true);
        if (data == "mc" && ch==2) sf_Trigger = scale_factor(&sf_triggermumu_H, (*selectedLeptonsEleScaleUp)[0]->pt_, (*selectedLeptonsEleScaleUp)[1]->pt_,"central",false, true);
      }
      if(leptonPass && triggerPass){
        weight_lepB = weight_Lumi * signnum_typical(LHEWeight_originalXWGTUP) * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowheg * ttKFactor;
        weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowheg * ttKFactor;
        if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
        else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
        if(nbjet==1 && DyPass){
          recoL2 = (*selectedLeptonsEleScaleUp)[1]->p4_;
          for (int l=0;l<selectedJets->size();l++){
            if((*selectedJets)[l]->btag_){
              recoBjet = (*selectedJets)[l]->p4_;
              break;
            }
          }
          recoNu = Wneutrino(MET_T1_pt, MET_T1_phi, (*selectedLeptonsEleScaleUp)[1]->pt_, (*selectedLeptonsEleScaleUp)[1]->eta_, (*selectedLeptonsEleScaleUp)[1]->phi_);
          recoW = recoNu + recoL2;
          recoTop = recoW + recoBjet;
          MVA_lep1Pt=(*selectedLeptonsEleScaleUp)[0]->pt_;
          MVA_lep2Pt=(*selectedLeptonsEleScaleUp)[1]->pt_;
          MVA_llM=((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).M();
          MVA_llPt=((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).Pt();
          MVA_llDr=deltaR((*selectedLeptonsEleScaleUp)[0]->eta_,(*selectedLeptonsEleScaleUp)[0]->phi_,(*selectedLeptonsEleScaleUp)[1]->eta_,(*selectedLeptonsEleScaleUp)[1]->phi_);
          MVA_llDphi=abs(deltaPhi((*selectedLeptonsEleScaleUp)[0]->phi_,(*selectedLeptonsEleScaleUp)[1]->phi_));
          MVA_topL1Dphi=abs(deltaPhi((*selectedLeptonsEleScaleUp)[0]->phi_,recoTop.Phi()));
          MVA_topL1Dr=deltaR((*selectedLeptonsEleScaleUp)[0]->eta_,(*selectedLeptonsEleScaleUp)[0]->phi_,recoTop.Eta(),recoTop.Phi());
          MVA_topL1DptOsumPt=(abs((*selectedLeptonsEleScaleUp)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsEleScaleUp)[0]->pt_ + recoTop.Pt());
          MVA_topPt=recoTop.Pt();
          MVAoutput = readerMVA->EvaluateMVA( "BDTG");
  
          HistsSysUp[ch][2][0][12]->Fill((*selectedLeptonsEleScaleUp)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][1][12]->Fill((*selectedLeptonsEleScaleUp)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][2][12]->Fill((*selectedLeptonsEleScaleUp)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][3][12]->Fill((*selectedLeptonsEleScaleUp)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][4][12]->Fill((*selectedLeptonsEleScaleUp)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][5][12]->Fill((*selectedLeptonsEleScaleUp)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][6][12]->Fill(((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][7][12]->Fill(((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][8][12]->Fill(deltaR((*selectedLeptonsEleScaleUp)[0]->eta_,(*selectedLeptonsEleScaleUp)[0]->phi_,(*selectedLeptonsEleScaleUp)[1]->eta_,(*selectedLeptonsEleScaleUp)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][9][12]->Fill(deltaPhi((*selectedLeptonsEleScaleUp)[0]->phi_,(*selectedLeptonsEleScaleUp)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][10][12]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][11][12]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][12][12]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][13][12]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][14][12]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][15][12]->Fill(MET_T1_pt,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][16][12]->Fill(MET_T1_phi,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][17][12]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][18][12]->Fill(((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][19][12]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][20][12]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][21][12]->Fill(abs(deltaPhi((*selectedLeptonsEleScaleUp)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][22][12]->Fill(deltaR((*selectedLeptonsEleScaleUp)[0]->eta_,(*selectedLeptonsEleScaleUp)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][23][12]->Fill((abs((*selectedLeptonsEleScaleUp)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsEleScaleUp)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][24][12]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        if(nbjet>1 && DyPass){
          HistsSysUp[ch][3][0][12]->Fill((*selectedLeptonsEleScaleUp)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][1][12]->Fill((*selectedLeptonsEleScaleUp)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][2][12]->Fill((*selectedLeptonsEleScaleUp)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][3][12]->Fill((*selectedLeptonsEleScaleUp)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][4][12]->Fill((*selectedLeptonsEleScaleUp)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][5][12]->Fill((*selectedLeptonsEleScaleUp)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][6][12]->Fill(((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][7][12]->Fill(((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][8][12]->Fill(deltaR((*selectedLeptonsEleScaleUp)[0]->eta_,(*selectedLeptonsEleScaleUp)[0]->phi_,(*selectedLeptonsEleScaleUp)[1]->eta_,(*selectedLeptonsEleScaleUp)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][9][12]->Fill(deltaPhi((*selectedLeptonsEleScaleUp)[0]->phi_,(*selectedLeptonsEleScaleUp)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][10][12]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][11][12]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][12][12]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][13][12]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][14][12]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][15][12]->Fill(MET_T1_pt,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][16][12]->Fill(MET_T1_phi,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][17][12]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][18][12]->Fill(((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][19][12]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][20][12]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][21][12]->Fill(abs(deltaPhi((*selectedLeptonsEleScaleUp)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][22][12]->Fill(deltaR((*selectedLeptonsEleScaleUp)[0]->eta_,(*selectedLeptonsEleScaleUp)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][23][12]->Fill((abs((*selectedLeptonsEleScaleUp)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsEleScaleUp)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][24][12]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        delete eft_fit;
      } 

      leptonPass=false;
      triggerPass=false;
      DyPass = true;
      ch=10;
      if(selectedLeptonsEleScaleDown->size() ==2){
        if (((*selectedLeptonsEleScaleDown)[0]->pt_ > 25) &&
        ((*selectedLeptonsEleScaleDown)[0]->charge_ * (*selectedLeptonsEleScaleDown)[1]->charge_ == -1) &&
        ((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).M()>20) leptonPass=true; 
        if ((*selectedLeptonsEleScaleDown)[0]->lep_ + (*selectedLeptonsEleScaleDown)[1]->lep_ == 2) ch = 0;
        if ((*selectedLeptonsEleScaleDown)[0]->lep_ + (*selectedLeptonsEleScaleDown)[1]->lep_ == 11) ch = 1;
        if ((*selectedLeptonsEleScaleDown)[0]->lep_ + (*selectedLeptonsEleScaleDown)[1]->lep_ == 20) ch = 2;
        if(ch ==0 && triggerPassEE) triggerPass=true;
        if(ch ==1 && triggerPassEMu) triggerPass=true;
        if(ch ==2 && triggerPassMuMu) triggerPass=true;
        if ((*selectedLeptonsEleScaleDown)[0]->lep_ + (*selectedLeptonsEleScaleDown)[1]->lep_ != 11 && ((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).M()<106 ) DyPass = false;
        if (data == "mc" && ch==0) sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptonsEleScaleDown)[0]->pt_, (*selectedLeptonsEleScaleDown)[1]->pt_,"central",false, true);
        if (data == "mc" && ch==1) sf_Trigger = scale_factor(&sf_triggeremu_H, (*selectedLeptonsEleScaleDown)[0]->pt_, (*selectedLeptonsEleScaleDown)[1]->pt_,"central",false, true);
        if (data == "mc" && ch==2) sf_Trigger = scale_factor(&sf_triggermumu_H, (*selectedLeptonsEleScaleDown)[0]->pt_, (*selectedLeptonsEleScaleDown)[1]->pt_,"central",false, true);
      }
      if(leptonPass && triggerPass){
        weight_lepB = weight_Lumi * signnum_typical(LHEWeight_originalXWGTUP) * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowheg * ttKFactor * sf_JetPuId;
        weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowheg * ttKFactor * sf_JetPuId;
        if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
        else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
        if(nbjet==1 && DyPass){
          recoL2 = (*selectedLeptonsEleScaleDown)[1]->p4_;
          for (int l=0;l<selectedJets->size();l++){
            if((*selectedJets)[l]->btag_){
              recoBjet = (*selectedJets)[l]->p4_;
              break;
            }
          }
          recoNu = Wneutrino(MET_T1_pt, MET_T1_phi, (*selectedLeptonsEleScaleDown)[1]->pt_, (*selectedLeptonsEleScaleDown)[1]->eta_, (*selectedLeptonsEleScaleDown)[1]->phi_);
          recoW = recoNu + recoL2;
          recoTop = recoW + recoBjet;
          MVA_lep1Pt=(*selectedLeptonsEleScaleDown)[0]->pt_;
          MVA_lep2Pt=(*selectedLeptonsEleScaleDown)[1]->pt_;
          MVA_llM=((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).M();
          MVA_llPt=((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).Pt();
          MVA_llDr=deltaR((*selectedLeptonsEleScaleDown)[0]->eta_,(*selectedLeptonsEleScaleDown)[0]->phi_,(*selectedLeptonsEleScaleDown)[1]->eta_,(*selectedLeptonsEleScaleDown)[1]->phi_);
          MVA_llDphi=abs(deltaPhi((*selectedLeptonsEleScaleDown)[0]->phi_,(*selectedLeptonsEleScaleDown)[1]->phi_));
          MVA_topL1Dphi=abs(deltaPhi((*selectedLeptonsEleScaleDown)[0]->phi_,recoTop.Phi()));
          MVA_topL1Dr=deltaR((*selectedLeptonsEleScaleDown)[0]->eta_,(*selectedLeptonsEleScaleDown)[0]->phi_,recoTop.Eta(),recoTop.Phi());
          MVA_topL1DptOsumPt=(abs((*selectedLeptonsEleScaleDown)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsEleScaleDown)[0]->pt_ + recoTop.Pt());
          MVA_topPt=recoTop.Pt();
          MVAoutput = readerMVA->EvaluateMVA( "BDTG");
  
          HistsSysDown[ch][2][0][12]->Fill((*selectedLeptonsEleScaleDown)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][1][12]->Fill((*selectedLeptonsEleScaleDown)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][2][12]->Fill((*selectedLeptonsEleScaleDown)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][3][12]->Fill((*selectedLeptonsEleScaleDown)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][4][12]->Fill((*selectedLeptonsEleScaleDown)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][5][12]->Fill((*selectedLeptonsEleScaleDown)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][6][12]->Fill(((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][7][12]->Fill(((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][8][12]->Fill(deltaR((*selectedLeptonsEleScaleDown)[0]->eta_,(*selectedLeptonsEleScaleDown)[0]->phi_,(*selectedLeptonsEleScaleDown)[1]->eta_,(*selectedLeptonsEleScaleDown)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][9][12]->Fill(deltaPhi((*selectedLeptonsEleScaleDown)[0]->phi_,(*selectedLeptonsEleScaleDown)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][10][12]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][11][12]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][12][12]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][13][12]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][14][12]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][15][12]->Fill(MET_T1_pt,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][16][12]->Fill(MET_T1_phi,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][17][12]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][18][12]->Fill(((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][19][12]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][20][12]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][21][12]->Fill(abs(deltaPhi((*selectedLeptonsEleScaleDown)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][22][12]->Fill(deltaR((*selectedLeptonsEleScaleDown)[0]->eta_,(*selectedLeptonsEleScaleDown)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][23][12]->Fill((abs((*selectedLeptonsEleScaleDown)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsEleScaleDown)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][24][12]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        if(nbjet>1 && DyPass){
          HistsSysDown[ch][3][0][12]->Fill((*selectedLeptonsEleScaleDown)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][1][12]->Fill((*selectedLeptonsEleScaleDown)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][2][12]->Fill((*selectedLeptonsEleScaleDown)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][3][12]->Fill((*selectedLeptonsEleScaleDown)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][4][12]->Fill((*selectedLeptonsEleScaleDown)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][5][12]->Fill((*selectedLeptonsEleScaleDown)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][6][12]->Fill(((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][7][12]->Fill(((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][8][12]->Fill(deltaR((*selectedLeptonsEleScaleDown)[0]->eta_,(*selectedLeptonsEleScaleDown)[0]->phi_,(*selectedLeptonsEleScaleDown)[1]->eta_,(*selectedLeptonsEleScaleDown)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][9][12]->Fill(deltaPhi((*selectedLeptonsEleScaleDown)[0]->phi_,(*selectedLeptonsEleScaleDown)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][10][12]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][11][12]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][12][12]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][13][12]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][14][12]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][15][12]->Fill(MET_T1_pt,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][16][12]->Fill(MET_T1_phi,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][17][12]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][18][12]->Fill(((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][19][12]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][20][12]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][21][12]->Fill(abs(deltaPhi((*selectedLeptonsEleScaleDown)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][22][12]->Fill(deltaR((*selectedLeptonsEleScaleDown)[0]->eta_,(*selectedLeptonsEleScaleDown)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][23][12]->Fill((abs((*selectedLeptonsEleScaleDown)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsEleScaleDown)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][24][12]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        delete eft_fit;
      } 

      leptonPass=false;
      triggerPass=false;
      DyPass = true;
      ch=10;
      if(selectedLeptonsMuResUp->size() ==2){
        if (((*selectedLeptonsMuResUp)[0]->pt_ > 25) &&
        ((*selectedLeptonsMuResUp)[0]->charge_ * (*selectedLeptonsMuResUp)[1]->charge_ == -1) &&
        ((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).M()>20) leptonPass=true; 
        if ((*selectedLeptonsMuResUp)[0]->lep_ + (*selectedLeptonsMuResUp)[1]->lep_ == 2) ch = 0;
        if ((*selectedLeptonsMuResUp)[0]->lep_ + (*selectedLeptonsMuResUp)[1]->lep_ == 11) ch = 1;
        if ((*selectedLeptonsMuResUp)[0]->lep_ + (*selectedLeptonsMuResUp)[1]->lep_ == 20) ch = 2;
        if(ch ==0 && triggerPassEE) triggerPass=true;
        if(ch ==1 && triggerPassEMu) triggerPass=true;
        if(ch ==2 && triggerPassMuMu) triggerPass=true;
        if ((*selectedLeptonsMuResUp)[0]->lep_ + (*selectedLeptonsMuResUp)[1]->lep_ != 11 && ((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).M()<106 ) DyPass = false;
        if (data == "mc" && ch==0) sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptonsMuResUp)[0]->pt_, (*selectedLeptonsMuResUp)[1]->pt_,"central",false, true);
        if (data == "mc" && ch==1) sf_Trigger = scale_factor(&sf_triggeremu_H, (*selectedLeptonsMuResUp)[0]->pt_, (*selectedLeptonsMuResUp)[1]->pt_,"central",false, true);
        if (data == "mc" && ch==2) sf_Trigger = scale_factor(&sf_triggermumu_H, (*selectedLeptonsMuResUp)[0]->pt_, (*selectedLeptonsMuResUp)[1]->pt_,"central",false, true);
      }
      if(leptonPass && triggerPass){
        weight_lepB = weight_Lumi * signnum_typical(LHEWeight_originalXWGTUP) * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowheg * ttKFactor * sf_JetPuId;
        weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowheg * ttKFactor * sf_JetPuId;
        if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
        else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
        if(nbjet==1 && DyPass){
          recoL2 = (*selectedLeptonsMuResUp)[1]->p4_;
          for (int l=0;l<selectedJets->size();l++){
            if((*selectedJets)[l]->btag_){
              recoBjet = (*selectedJets)[l]->p4_;
              break;
            }
          }
          recoNu = Wneutrino(MET_T1_pt, MET_T1_phi, (*selectedLeptonsMuResUp)[1]->pt_, (*selectedLeptonsMuResUp)[1]->eta_, (*selectedLeptonsMuResUp)[1]->phi_);
          recoW = recoNu + recoL2;
          recoTop = recoW + recoBjet;
          MVA_lep1Pt=(*selectedLeptonsMuResUp)[0]->pt_;
          MVA_lep2Pt=(*selectedLeptonsMuResUp)[1]->pt_;
          MVA_llM=((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).M();
          MVA_llPt=((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).Pt();
          MVA_llDr=deltaR((*selectedLeptonsMuResUp)[0]->eta_,(*selectedLeptonsMuResUp)[0]->phi_,(*selectedLeptonsMuResUp)[1]->eta_,(*selectedLeptonsMuResUp)[1]->phi_);
          MVA_llDphi=abs(deltaPhi((*selectedLeptonsMuResUp)[0]->phi_,(*selectedLeptonsMuResUp)[1]->phi_));
          MVA_topL1Dphi=abs(deltaPhi((*selectedLeptonsMuResUp)[0]->phi_,recoTop.Phi()));
          MVA_topL1Dr=deltaR((*selectedLeptonsMuResUp)[0]->eta_,(*selectedLeptonsMuResUp)[0]->phi_,recoTop.Eta(),recoTop.Phi());
          MVA_topL1DptOsumPt=(abs((*selectedLeptonsMuResUp)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuResUp)[0]->pt_ + recoTop.Pt());
          MVA_topPt=recoTop.Pt();
          MVAoutput = readerMVA->EvaluateMVA( "BDTG");
  
          HistsSysUp[ch][2][0][13]->Fill((*selectedLeptonsMuResUp)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][1][13]->Fill((*selectedLeptonsMuResUp)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][2][13]->Fill((*selectedLeptonsMuResUp)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][3][13]->Fill((*selectedLeptonsMuResUp)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][4][13]->Fill((*selectedLeptonsMuResUp)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][5][13]->Fill((*selectedLeptonsMuResUp)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][6][13]->Fill(((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][7][13]->Fill(((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][8][13]->Fill(deltaR((*selectedLeptonsMuResUp)[0]->eta_,(*selectedLeptonsMuResUp)[0]->phi_,(*selectedLeptonsMuResUp)[1]->eta_,(*selectedLeptonsMuResUp)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][9][13]->Fill(deltaPhi((*selectedLeptonsMuResUp)[0]->phi_,(*selectedLeptonsMuResUp)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][10][13]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][11][13]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][12][13]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][13][13]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][14][13]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][15][13]->Fill(MET_T1_pt,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][16][13]->Fill(MET_T1_phi,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][17][13]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][18][13]->Fill(((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][19][13]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysUp[ch][2][20][13]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][21][13]->Fill(abs(deltaPhi((*selectedLeptonsMuResUp)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][22][13]->Fill(deltaR((*selectedLeptonsMuResUp)[0]->eta_,(*selectedLeptonsMuResUp)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][23][13]->Fill((abs((*selectedLeptonsMuResUp)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuResUp)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysUp[ch][2][24][13]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        if(nbjet>1 && DyPass){
          HistsSysUp[ch][3][0][13]->Fill((*selectedLeptonsMuResUp)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][1][13]->Fill((*selectedLeptonsMuResUp)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][2][13]->Fill((*selectedLeptonsMuResUp)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][3][13]->Fill((*selectedLeptonsMuResUp)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][4][13]->Fill((*selectedLeptonsMuResUp)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][5][13]->Fill((*selectedLeptonsMuResUp)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][6][13]->Fill(((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][7][13]->Fill(((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][8][13]->Fill(deltaR((*selectedLeptonsMuResUp)[0]->eta_,(*selectedLeptonsMuResUp)[0]->phi_,(*selectedLeptonsMuResUp)[1]->eta_,(*selectedLeptonsMuResUp)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][9][13]->Fill(deltaPhi((*selectedLeptonsMuResUp)[0]->phi_,(*selectedLeptonsMuResUp)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][10][13]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][11][13]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][12][13]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][13][13]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][14][13]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][15][13]->Fill(MET_T1_pt,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][16][13]->Fill(MET_T1_phi,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][17][13]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][18][13]->Fill(((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][19][13]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysUp[ch][3][20][13]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][21][13]->Fill(abs(deltaPhi((*selectedLeptonsMuResUp)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][22][13]->Fill(deltaR((*selectedLeptonsMuResUp)[0]->eta_,(*selectedLeptonsMuResUp)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][23][13]->Fill((abs((*selectedLeptonsMuResUp)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuResUp)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysUp[ch][3][24][13]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        delete eft_fit;
      } 

      leptonPass=false;
      triggerPass=false;
      DyPass = true;
      ch=10;
      if(selectedLeptonsMuResDown->size() ==2){
        if (((*selectedLeptonsMuResDown)[0]->pt_ > 25) &&
        ((*selectedLeptonsMuResDown)[0]->charge_ * (*selectedLeptonsMuResDown)[1]->charge_ == -1) &&
        ((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).M()>20) leptonPass=true; 
        if ((*selectedLeptonsMuResDown)[0]->lep_ + (*selectedLeptonsMuResDown)[1]->lep_ == 2) ch = 0;
        if ((*selectedLeptonsMuResDown)[0]->lep_ + (*selectedLeptonsMuResDown)[1]->lep_ == 11) ch = 1;
        if ((*selectedLeptonsMuResDown)[0]->lep_ + (*selectedLeptonsMuResDown)[1]->lep_ == 20) ch = 2;
        if(ch ==0 && triggerPassEE) triggerPass=true;
        if(ch ==1 && triggerPassEMu) triggerPass=true;
        if(ch ==2 && triggerPassMuMu) triggerPass=true;
        if ((*selectedLeptonsMuResDown)[0]->lep_ + (*selectedLeptonsMuResDown)[1]->lep_ != 11 && ((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).M()<106 ) DyPass = false;
        if (data == "mc" && ch==0) sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptonsMuResDown)[0]->pt_, (*selectedLeptonsMuResDown)[1]->pt_,"central",false, true);
        if (data == "mc" && ch==1) sf_Trigger = scale_factor(&sf_triggeremu_H, (*selectedLeptonsMuResDown)[0]->pt_, (*selectedLeptonsMuResDown)[1]->pt_,"central",false, true);
        if (data == "mc" && ch==2) sf_Trigger = scale_factor(&sf_triggermumu_H, (*selectedLeptonsMuResDown)[0]->pt_, (*selectedLeptonsMuResDown)[1]->pt_,"central",false, true);
      }
      if(leptonPass && triggerPass){
        weight_lepB = weight_Lumi * signnum_typical(LHEWeight_originalXWGTUP) * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowheg * ttKFactor * sf_JetPuId;
        weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowheg * ttKFactor * sf_JetPuId;
        if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
        else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
        if(nbjet==1 && DyPass){
          recoL2 = (*selectedLeptonsMuResDown)[1]->p4_;
          for (int l=0;l<selectedJets->size();l++){
            if((*selectedJets)[l]->btag_){
              recoBjet = (*selectedJets)[l]->p4_;
              break;
            }
          }
          recoNu = Wneutrino(MET_T1_pt, MET_T1_phi, (*selectedLeptonsMuResDown)[1]->pt_, (*selectedLeptonsMuResDown)[1]->eta_, (*selectedLeptonsMuResDown)[1]->phi_);
          recoW = recoNu + recoL2;
          recoTop = recoW + recoBjet;
          MVA_lep1Pt=(*selectedLeptonsMuResDown)[0]->pt_;
          MVA_lep2Pt=(*selectedLeptonsMuResDown)[1]->pt_;
          MVA_llM=((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).M();
          MVA_llPt=((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).Pt();
          MVA_llDr=deltaR((*selectedLeptonsMuResDown)[0]->eta_,(*selectedLeptonsMuResDown)[0]->phi_,(*selectedLeptonsMuResDown)[1]->eta_,(*selectedLeptonsMuResDown)[1]->phi_);
          MVA_llDphi=abs(deltaPhi((*selectedLeptonsMuResDown)[0]->phi_,(*selectedLeptonsMuResDown)[1]->phi_));
          MVA_topL1Dphi=abs(deltaPhi((*selectedLeptonsMuResDown)[0]->phi_,recoTop.Phi()));
          MVA_topL1Dr=deltaR((*selectedLeptonsMuResDown)[0]->eta_,(*selectedLeptonsMuResDown)[0]->phi_,recoTop.Eta(),recoTop.Phi());
          MVA_topL1DptOsumPt=(abs((*selectedLeptonsMuResDown)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuResDown)[0]->pt_ + recoTop.Pt());
          MVA_topPt=recoTop.Pt();
          MVAoutput = readerMVA->EvaluateMVA( "BDTG");
  
          HistsSysDown[ch][2][0][13]->Fill((*selectedLeptonsMuResDown)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][1][13]->Fill((*selectedLeptonsMuResDown)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][2][13]->Fill((*selectedLeptonsMuResDown)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][3][13]->Fill((*selectedLeptonsMuResDown)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][4][13]->Fill((*selectedLeptonsMuResDown)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][5][13]->Fill((*selectedLeptonsMuResDown)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][6][13]->Fill(((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][7][13]->Fill(((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][8][13]->Fill(deltaR((*selectedLeptonsMuResDown)[0]->eta_,(*selectedLeptonsMuResDown)[0]->phi_,(*selectedLeptonsMuResDown)[1]->eta_,(*selectedLeptonsMuResDown)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][9][13]->Fill(deltaPhi((*selectedLeptonsMuResDown)[0]->phi_,(*selectedLeptonsMuResDown)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][10][13]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][11][13]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][12][13]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][13][13]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][14][13]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][15][13]->Fill(MET_T1_pt,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][16][13]->Fill(MET_T1_phi,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][17][13]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][18][13]->Fill(((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][19][13]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysDown[ch][2][20][13]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][21][13]->Fill(abs(deltaPhi((*selectedLeptonsMuResDown)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][22][13]->Fill(deltaR((*selectedLeptonsMuResDown)[0]->eta_,(*selectedLeptonsMuResDown)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][23][13]->Fill((abs((*selectedLeptonsMuResDown)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuResDown)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysDown[ch][2][24][13]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        if(nbjet>1 && DyPass){
          HistsSysDown[ch][3][0][13]->Fill((*selectedLeptonsMuResDown)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][1][13]->Fill((*selectedLeptonsMuResDown)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][2][13]->Fill((*selectedLeptonsMuResDown)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][3][13]->Fill((*selectedLeptonsMuResDown)[1]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][4][13]->Fill((*selectedLeptonsMuResDown)[1]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][5][13]->Fill((*selectedLeptonsMuResDown)[1]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][6][13]->Fill(((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][7][13]->Fill(((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).Pt(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][8][13]->Fill(deltaR((*selectedLeptonsMuResDown)[0]->eta_,(*selectedLeptonsMuResDown)[0]->phi_,(*selectedLeptonsMuResDown)[1]->eta_,(*selectedLeptonsMuResDown)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][9][13]->Fill(deltaPhi((*selectedLeptonsMuResDown)[0]->phi_,(*selectedLeptonsMuResDown)[1]->phi_),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][10][13]->Fill((*selectedJets)[0]->pt_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][11][13]->Fill((*selectedJets)[0]->eta_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][12][13]->Fill((*selectedJets)[0]->phi_,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][13][13]->Fill(selectedJets->size(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][14][13]->Fill(nbjet,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][15][13]->Fill(MET_T1_pt,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][16][13]->Fill(MET_T1_phi,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][17][13]->Fill(PV_npvs,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][18][13]->Fill(((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][19][13]->Fill(MVAoutput,weight_lepB, *eft_fit);
          HistsSysDown[ch][3][20][13]->Fill(recoTop.M(),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][21][13]->Fill(abs(deltaPhi((*selectedLeptonsMuResDown)[0]->phi_,recoTop.Phi())),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][22][13]->Fill(deltaR((*selectedLeptonsMuResDown)[0]->eta_,(*selectedLeptonsMuResDown)[0]->phi_,recoTop.Eta(),recoTop.Phi()),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][23][13]->Fill((abs((*selectedLeptonsMuResDown)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuResDown)[0]->pt_ + recoTop.Pt()),weight_lepB, *eft_fit);
          HistsSysDown[ch][3][24][13]->Fill(recoTop.Pt(),weight_lepB, *eft_fit);
        }
        delete eft_fit;
      } 


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

    for (int l=0;l<JECsysUp->size();l++){
      for (int n=0;n<(*JECsysUp)[l].size();n++){
        delete (*JECsysUp)[l][n];
      }
    }
    for (int l=0;l<JECsysDown->size();l++){
      for (int n=0;n<(*JECsysDown)[l].size();n++){
        delete (*JECsysDown)[l][n];
      }
    }

    JECsysUp->clear();
    JECsysUp->shrink_to_fit();
    JECsysDown->clear();
    JECsysDown->shrink_to_fit();
    JECsysNbtagUp->clear();
    JECsysNbtagDown->clear();
    JECsysMETUp->clear();
    JECsysMETDown->clear();
    JECsysMETUp->shrink_to_fit();
    JECsysMETDown->shrink_to_fit();
    delete JECsysUp;
    delete JECsysDown;
    delete JECsysNbtagUp;
    delete JECsysNbtagDown;
    delete JECsysMETUp;
    delete JECsysMETDown;

   }
  cout<<"Loop is completed"<<endl;
  cout<<"from "<<ntr<<" events, "<<nAccept<<" events are accepted"<<endl;
  cout<<"sumPuWeight "<<sumPuWeight<<endl;
  cout<<"sumPreFireWeight "<<sumPreFireWeight<<endl;
  cout<<"sumWeightMuID "<<sumWeightMuID<<endl;
  cout<<"sumWeightMuIso "<<sumWeightMuIso<<endl;
  cout<<"sumWeighttTrigger "<<sumWeighttTrigger<<endl;
  cout<<"sumWeightEleReco "<<sumWeightEleReco<<endl;
  cout<<"sumWeightEleID "<<sumWeightEleID<<endl;  
  cout<<"sumWeightBtag "<<sumWeightBtag<<endl;

  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int l=0;l<vars.size();++l){
        Hists[i][k][l]  ->Write("",TObject::kOverwrite);
        if(data=="mc"){
          for (int n=0;n<sys.size();++n){
            if(k<2) continue;
            HistsSysUp[i][k][l][n]->Write("",TObject::kOverwrite);
            HistsSysDown[i][k][l][n]->Write("",TObject::kOverwrite);
          }
        }
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

  if(data=="mc"){
    file_out.mkdir("JECSys");
    file_out.cd("JECSys/");
    for (int i=0;i<channels.size();++i){
      for (int k=2;k<regions.size();++k){
        for (int l=0;l<1;++l){
          for (int n=0;n<sysJecNames.size();++n){
            HistsJECUp[i][k-2][l][n]->Write("",TObject::kOverwrite);
            HistsJECDown[i][k-2][l][n]->Write("",TObject::kOverwrite);
          }
        }
      }
    }
  }
  file_out.cd("");

  if ((fname.Contains("TTTo2L2Nu")&& !fname.Contains("sys")) || fname.Contains("BNV")){
    file_out.mkdir("reweightingSys");
    file_out.cd("reweightingSys/");
    for (int i=0;i<channels.size();++i){
      for (int k=2;k<regions.size();++k){
        for (int l=0;l<1;++l){
          for (int n=0;n<nScale;++n){
            HistsSysReweightsQscale[i][k-2][l][n]->Write("",TObject::kOverwrite);
          }
          for (int n=0;n<nPdf;++n){
            HistsSysReweightsPDF[i][k-2][l][n]->Write("",TObject::kOverwrite);
          }
          for (int n=0;n<nPS;++n){
            HistsSysReweightsPS[i][k-2][l][n]->Write("",TObject::kOverwrite);
          }
        }
      }
    }
  }

  file_out.cd("");
  tree_out.Write() ;
  file_out.Close() ;
cout<<"Cleaning the memory"<<endl;

Hists.clear();
cout<<"Hists cleaned"<<endl;
HistsSysUp.clear();
cout<<"HistsSysUp cleaned"<<endl;
HistsSysDown.clear();
cout<<"HistsSysDown cleaned"<<endl;
HistsJECUp.clear();
cout<<"HistsJECUp cleaned"<<endl;
HistsJECDown.clear();
cout<<"HistsJECDown cleaned"<<endl;
HistsSysReweightsPDF.clear();
cout<<"HistsSysReweightPDF cleaned"<<endl;
HistsSysReweightsQscale.clear();
cout<<"HistsSysReweightsQscale cleaned"<<endl;
HistsSysReweightsPS.clear();
cout<<"HistsSysReweightPS cleaned"<<endl;

  cout<<"Job is finished"<<endl;
}

