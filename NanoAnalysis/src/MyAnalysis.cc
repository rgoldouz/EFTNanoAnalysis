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
#include "XYMETCorrection_withUL17andUL18andUL16.h"
#include "lester_mt2_bisect.h"
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
#include "sys/types.h"
#include "sys/sysinfo.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "TMath.h"
#include <math.h>     
#endif

using namespace correction;
using namespace std;

void resetVec(std::vector<int> &K){
  for (int i=0;i<K.size();++i){  
    K[i]=-1;
  }
}

bool useRegions(std::vector<int> K){
  bool ifuse=false;
  for (int i=0;i<K.size();++i){
    if(K[i]>=0) ifuse=true;
  }
  return ifuse;
}

float DYMassPowheg(float M){
  return (1.0678 - 0.000120666*M + 0.0000000322*pow(M,2) - 0.0000000000039*pow(M,3));
}

float topPtPowhegData(float pt){
//Data top pt reweighting
  if (pt<500) return TMath::Exp(0.0615- (0.0005*pt));
  else return TMath::Exp(0.0615- (0.0005*500));
//MC top pt reweighting
//  return (0.973 - (0.000134 * pt) + (0.103 * exp(pt * (-0.0118))));
}

float topPtPowhegMC(float pt){
//MC top pt reweighting
  return (0.973 - (0.000134 * pt) + (0.103 * exp(pt * (-0.0118))));
}

float topPtMGLO(float x){
  return (0.688 -  0.0000174*x + 0.824*exp(-0.0000253*x)/(pow(x,0.2185)));
}

int vInd(std::map<TString, std::vector<float>> V, TString name, int i){
  return V.find(name)->second.at(i);
}

bool ifSysNeeded(std::vector<lepton_candidate*> *lep, float cut){
  bool pass=false;
  if (lep->size()==2){
    if(((*lep)[0]->p4_ + (*lep)[1]->p4_).M()> cut) pass=true;
  }
  return pass;
}     

int getVecPos(std::vector<TString> vec, string element){
    int i;
    for(i = 0; i < vec.size(); i++){
        if(vec[i] == element){
            break;
        }
    }
    if(i == vec.size()){
        std::cout<<"No such element as "<<element<<" found. Please enter again: ";
        std::cin>>element;
        i = getVecPos(vec, element);
    }

    return i;
}

void MyAnalysis::Loop(TString fname, TString data, TString dataset ,string year, TString RUN, float xs, float lumi, float Nevent, int iseft, int nRuns){
  bool ifSys=true;
  bool ifCR=false;
  double memoryInit=getValue();
  TH1EFT  *crossSection = new TH1EFT("crossSection","crossSection",1,0,1);
  TH1F  *eleEffNum = new TH1F("eleEffNum","eleEffNum",10,0,1000);
  TH1F  *eleEffDen = new TH1F("eleEffDen","eleEffDen",10,0,1000);
  TH1F  *truePU = new TH1F("truePU","truePU",100,0,100);

  TH2F  btagEff_b_H;
  TH2F  btagEff_c_H;
  TH2F  btagEff_udsg_H;
  TH2F  sf_triggeree_H;
  TH2F  sf_triggeremu_H;
  TH2F  sf_triggermumu_H;
  TH2F  jetVetoMaps_H;
  TH2F  highPtMuRecoSF_pVsAbsEta_H;
  TH2F  eff_triggermumu_mc_H;
  TH2F  eff_triggermumu_data_H;
  TH2F  eff_triggermumu_dataUp_H;
  TH2F  eff_triggermumu_dataDown_H;
  TGraphAsymmErrors* TTptQS_TG = new TGraphAsymmErrors();;
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

    TFile *f_HighPtMuRecoSF = new TFile(("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/HighPtMuRecoSF_" + year + ".root").c_str());
    highPtMuRecoSF_pVsAbsEta_H = *(TH2F*)f_HighPtMuRecoSF->Get("h2_HighPtMuRecoSF_pVsAbsEta");

    TFile *TTptQS = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/TTSYS.root");
    TTptQS_TG = (TGraphAsymmErrors*) TTptQS->Get("PTtop_QS")->Clone();

    f_btagEff_Map->Close();
    f_trigger->Close();
    f_HighPtMuRecoSF->Close();
    TTptQS->Close();
    delete f_btagEff_Map;
    delete f_trigger;
    delete f_HighPtMuRecoSF;
    delete TTptQS;
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
    delete Map2016preVFP;
    TFile *f_triggerEffMuMu = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root");
    eff_triggermumu_mc_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_DY_var");
    eff_triggermumu_data_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_Data_var");
    eff_triggermumu_dataUp_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_Data_errorUpper");
    eff_triggermumu_dataDown_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_Data_errorLower");
    f_triggerEffMuMu->Close();
    delete f_triggerEffMuMu;
  }
  if(year == "2016postVFP"){
    eleSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/EGM/2016postVFP_UL/electron.json.gz";
    muSF= "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/MUO/2016postVFP_UL/muon_Z.json.gz";
    bSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/BTV/2016postVFP_UL/btagging.json.gz";
    jetSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/JME/2016postVFP_UL/UL16postVFP_jmar.json.gz";
    TFile *Map2016postVFP = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/hotjets-UL16.root");
    jetVetoMaps_H = *(TH2F*)Map2016postVFP->Get("h2hot_ul16_plus_hbm2_hbp12_qie11");
    Map2016postVFP->Close();
    delete Map2016postVFP;
    TFile *f_triggerEffMuMu = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root");
    eff_triggermumu_mc_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_DY_var");
    eff_triggermumu_data_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_Data_var");
    eff_triggermumu_dataUp_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_Data_errorUpper");
    eff_triggermumu_dataDown_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2016_Data_errorLower");
    f_triggerEffMuMu->Close();
    delete f_triggerEffMuMu;
  }
  if(year == "2017"){
    eleSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/EGM/2017_UL/electron.json.gz";
    muSF= "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/MUO/2017_UL/muon_Z.json.gz";
    bSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/BTV/2017_UL/btagging.json.gz";
    jetSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/JME/2017_UL/UL17_jmar.json.gz";
    TFile *Map2017 = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/hotjets-UL17_v2.root");
    jetVetoMaps_H = *(TH2F*)Map2017->Get("h2hot_ul17_plus_hep17_plus_hbpw89");
    Map2017->Close();
    delete Map2017;
    TFile *f_triggerEffMuMu = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root");
    eff_triggermumu_mc_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2017_DY_var");
    eff_triggermumu_data_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2017_Data_var");
    eff_triggermumu_dataUp_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2017_Data_errorUpper");
    eff_triggermumu_dataDown_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2017_Data_errorLower");
    f_triggerEffMuMu->Close();
    delete f_triggerEffMuMu;
  }
  if(year == "2018"){
    eleSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/EGM/2018_UL/electron.json.gz";
    muSF= "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/MUO/2018_UL/muon_Z.json.gz";
    bSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/BTV/2018_UL/btagging.json.gz";
    jetSF="/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/data/POG/JME/2018_UL/UL18_jmar.json.gz";
    TFile *Map2018 = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/hotjets-UL18.root");
    jetVetoMaps_H = *(TH2F*)Map2018->Get("h2hot_ul18_plus_hem1516_and_hbp2m1");
    Map2018->Close();
    delete Map2018;
    TFile *f_triggerEffMuMu = new TFile("/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root");
    eff_triggermumu_mc_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2018_DY_var");
    eff_triggermumu_data_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2018_Data_var");
    eff_triggermumu_dataUp_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2018_Data_errorUpper");
    eff_triggermumu_dataDown_H =  *(TH2F*)f_triggerEffMuMu->Get("Eff_2018_Data_errorLower");
    f_triggerEffMuMu->Close();
    delete f_triggerEffMuMu;
  }

  map<string,float> HEEPSF_B{ {"2016preVFP",0.985}, { "2016postVFP",0.985},{"2017",0.979},{"2018",0.973} };
  map<string,float> HEEPSF_E{ {"2016preVFP",0.990}, { "2016postVFP",0.990},{"2017",0.987},{"2018",0.980} };

  auto csetFileEleSF = CorrectionSet::from_file(eleSF);
  auto csetEleIdReco = csetFileEleSF->at("UL-Electron-ID-SF");

  auto csetFileMuSF = CorrectionSet::from_file(muSF);
  auto csetMuTightId = csetFileMuSF->at("NUM_HighPtID_DEN_TrackerMuons");
  auto csetMuTightRelIso = csetFileMuSF->at("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut");

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

  std::vector<TString> regions{"ll","llOffZ","llB1", "llBg1"};
  if (ifCR && !ifSys ) {
    regions.push_back("llB1InZ");
    regions.push_back("llB1-BDTm1to0");
    regions.push_back("llB1-BDT0to0p8");
    regions.push_back("llB1-BDT0p8to1");
    regions.push_back("llptl1G300");
    regions.push_back("llB1ptl1G300");
    regions.push_back("ll-BB");
    regions.push_back("ll-BE");
    regions.push_back("ll-EE");
    regions.push_back("ll-leadMuP");
    regions.push_back("ll-leadMuN");
  }
  std::vector<TString> channels{"ee", "emu", "mumu","ll"};

  std::map<TString, std::vector<float>> vars =
  {
    {"lep1Pt",                         {0,      60,   0,  1500, -1}},
    {"lep1Eta",                        {1,      20,   -3, 3   , -1}},
    {"lep1Phi",                        {2,      25,   -4, 4   , -1}},
    {"lep2Pt",                         {3,      25,   0,  1000, -1}},
    {"lep2Eta",                        {4,      20,   -3, 3   , -1}},
    {"lep2Phi",                        {5,      25,   -4, 4   , -1}},
    {"llM",                            {6,      100,    0, 2000 , -1}},
    {"llPt",                           {7,      20,    0, 200 , -1}},
    {"llDr",                           {8,      25,    0, 7   , -1}},
    {"llDphi",                         {9,      15,    0, 4   , -1}},
    {"jet1Pt",                         {10,     20,    0, 300 , -1}},
    {"jet1Eta",                        {11,     20,    -3, 3  , -1}},
    {"jet1Phi",                        {12,     25,    -4, 4  , -1}},
    {"njet",                           {13,     10,    0, 10  , -1}},
    {"nbjet",                          {14,     6,     0, 6   , -1}},
    {"Met",                            {15,     30,    0, 210 , -1}},
    {"MetPhi",                         {16,     20,    -4, 4  , -1}},
    {"nVtx",                           {17,     70,    0, 70  , -1}},
    {"llMZw",                          {18,     80,    70, 110, -1}},
    {"BDT",                            {19,     100,    -1, 1  , -1}},
    {"topMass",                        {20,     25,    0, 500 , -1}},
    {"topL1Dphi",                      {21,     15,    0, 4   , -1}},
    {"topL1Dr",                        {22,     25,    0, 7   , -1}},
    {"topL1DptOsumPt",                 {23,     20,    0, 1   , -1}},
    {"topPt",                          {24,     25,    0, 500 , -1}},
  };
  if (ifCR && !ifSys ) {
    vars["muTkIso"]= std::vector<float> {25,     20,    0, 0.2, -1};
    vars["leadElePt"]= std::vector<float>   {26,     60,    0, 1500, -1};
    vars["leadEleEta"]= std::vector<float>  {27,     20,   -3, 3, -1};
    vars["leadElePhi"]= std::vector<float>  {28,     25,   -4, 4, -1};
    vars["leadMuPt"]= std::vector<float>    {29,     60,    0, 1500, -1};
    vars["leadMuEta"]= std::vector<float>   {30,     20,   -3, 3, -1};
    vars["leadMuPhi"]= std::vector<float>   {31,     25,   -4, 4, -1};
    vars["mu1bDr"]= std::vector<float>       {32,      25,    0, 7, -1};
    vars["mu2bDr"]= std::vector<float>       {33,      25,    0, 7, -1};
    vars["mu1jDrMin"]= std::vector<float>    {34,      25,    0, 1, -1};
    vars["mu2jDrMin"]= std::vector<float>    {35,      25,    0, 1, -1};
    vars["PtRelMu1jet"]= std::vector<float>  {36,     30,    0, 300, -1};
    vars["PtRelMu2jet"]= std::vector<float>  {37,     30,    0, 300 , -1};
    vars["ele1bDr"]= std::vector<float>      {38,      25,    0, 7, -1};
    vars["ele2bDr"]= std::vector<float>      {39,      25,    0, 7, -1};
    vars["ele1jDrMin"]= std::vector<float>   {40,      25,    0, 1, -1};
    vars["ele2jDrMin"]= std::vector<float>   {41,      25,    0, 1, -1};
    vars["PtRelEle1jet"]= std::vector<float> {42,     30,    0, 300, -1};
    vars["PtRelEle2jet"]= std::vector<float> {43,     30,    0, 300, -1};
  };

//  D3HistsContainer Hists;
  Hists.resize(channels.size());
  for (int i=0;i<channels.size();++i){
    Hists[i].resize(regions.size());
    for (int k=0;k<regions.size();++k){
      Hists[i][k].resize(vars.size());
    }
  }

  std::vector<TString> varsWithSys{"BDT"};
//  std::vector<TString> varsWithSys{"BDT","lep1Pt","lep1Eta","lep2Pt","lep2Eta","llM","llPt","llDr","llDphi","njet","nbjet","Met","topL1Dphi","topL1Dr","topL1DptOsumPt","topPt"};
  std::stringstream name;
  TH1EFT *h_test;

  int ii=0;
  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for( auto it = vars.cbegin() ; it != vars.cend() ; ++it ){
//        if(i==0 && k==0) cout<<"Histogram is made for variable "<<it->first<<" at "<<it->second.at(0)<<endl;
        if(i==0 && k==0 && (std::find(varsWithSys.begin(), varsWithSys.end(), it->first) != varsWithSys.end())) {
          vars[it->first][4]=ii;
          ii++;
        }
        name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first;
        h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
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

  std::vector<TString> sys{"eleRecoSf", "eleIDSf", "JetPuID", "muIdIsoSf", "bcTagSf", "LTagSf","pu", "prefiring", "trigSF","jes", "jer","muonScale","electronScale","muonRes", "unclusMET", "bcTagSfUnCorr", "LTagSfUnCorr","topPt"};
  std::vector<TString> sysNotWeight{"jes", "jer","muonScale","electronScale","muonRes", "unclusMET"};

  if(data == "mc" && !fname.Contains("sys")  && ifSys){
    HistsSysUp.resize(channels.size());
    for (int i=0;i<channels.size();++i){
      HistsSysUp[i].resize(regions.size());
      for (int k=0;k<regions.size();++k){
        HistsSysUp[i][k].resize(varsWithSys.size());
        for (int n=0;n<varsWithSys.size();++n){
          HistsSysUp[i][k][n].resize(sys.size());
        }
      }
    }
  
    HistsSysDown.resize(channels.size());
    for (int i=0;i<channels.size();++i){
      HistsSysDown[i].resize(regions.size());
      for (int k=0;k<regions.size();++k){
        HistsSysDown[i][k].resize(varsWithSys.size());
        for (int n=0;n<varsWithSys.size();++n){
          HistsSysDown[i][k][n].resize(sys.size());
        }
      }
    }

    for (int i=0;i<channels.size();++i){
      for (int k=0;k<regions.size();++k){
        for( auto it = vars.cbegin() ; it != vars.cend() ; ++it ){
          for (int n=0;n<sys.size();++n){
            if (!(std::find(varsWithSys.begin(), varsWithSys.end(), it->first) != varsWithSys.end())) continue;
            name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_"<<sys[n]<<"_Up";
            h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
            h_test->StatOverflows(kTRUE);
            h_test->Sumw2(kTRUE);
            HistsSysUp[i][k][it->second.at(4)][n] = h_test;
            name.str("");
            name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_"<<sys[n]<<"_Down";
            h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
            h_test->StatOverflows(kTRUE);
            h_test->Sumw2(kTRUE);
            HistsSysDown[i][k][it->second.at(4)][n] = h_test;
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

  if(data == "mc" && !fname.Contains("sys") && ifSys){
    HistsJecUp.resize(channels.size());
    for (int i=0;i<channels.size();++i){
      HistsJecUp[i].resize(regions.size()-2);
      for (int k=0;k<regions.size()-2;++k){
        HistsJecUp[i][k].resize(1);
        for (int n=0;n<1;++n){
          HistsJecUp[i][k][n].resize(sysJecNames.size());
        }
      }
    }
  
    HistsJecDown.resize(channels.size());
    for (int i=0;i<channels.size();++i){
      HistsJecDown[i].resize(regions.size()-2);
      for (int k=0;k<regions.size()-2;++k){
        HistsJecDown[i][k].resize(1);
        for (int n=0;n<1;++n){
          HistsJecDown[i][k][n].resize(sysJecNames.size());
        }
      }
    }
    for (int i=0;i<channels.size();++i){
      for (int k=0;k<regions.size()-2;++k){
        for( auto it = vars.cbegin() ; it != vars.cend() ; ++it ){
          if(it->first !="BDT") continue;
          for (int n=0;n<sysJecNames.size();++n){
            name<<channels[i]<<"_"<<regions[k+2]<<"_"<<it->first<<"_"<<sysJecNames[n]<<"_Up";
            h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
            h_test->StatOverflows(kTRUE);
            h_test->Sumw2(kTRUE);
            HistsJecUp[i][k][0][n] = h_test;
            name.str("");
            name<<channels[i]<<"_"<<regions[k+2]<<"_"<<it->first<<"_"<<sysJecNames[n]<<"_Down";
            h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
            h_test->StatOverflows(kTRUE);
            h_test->Sumw2(kTRUE);
            HistsJecDown[i][k][0][n] = h_test;
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

  if(data == "mc" && ifSys){
    if ((fname.Contains("TTTo2L2Nu") && !fname.Contains("sys")) || fname.Contains("BNV")){
    HistsSysReweightsQscale.resize(channels.size());
    for (int i=0;i<channels.size();++i){
      HistsSysReweightsQscale[i].resize(2);
      for (int k=0;k<2;++k){
        HistsSysReweightsQscale[i][k].resize(1);
        for (int n=0;n<1;++n){
          HistsSysReweightsQscale[i][k][n].resize(nScale);
        }
      }
    }
  
    HistsSysReweightsPDF.resize(channels.size());
    for (int i=0;i<channels.size();++i){
      HistsSysReweightsPDF[i].resize(2);
      for (int k=0;k<2;++k){
        HistsSysReweightsPDF[i][k].resize(1);
        for (int n=0;n<1;++n){
          HistsSysReweightsPDF[i][k][n].resize(nPdf);
        }
      }
    }
  
    HistsSysReweightsPS.resize(channels.size());
    for (int i=0;i<channels.size();++i){
      HistsSysReweightsPS[i].resize(2);
      for (int k=0;k<2;++k){
        HistsSysReweightsPS[i][k].resize(1);
        for (int n=0;n<1;++n){
          HistsSysReweightsPS[i][k][n].resize(nPS);
        }
      }
    }
      for (int i=0;i<channels.size();++i){
        for (int k=2;k<regions.size();++k){
        for( auto it = vars.cbegin() ; it != vars.cend() ; ++it ){
            if(it->first !="BDT") continue;
            for (int n=0;n<nScale;++n){
              name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_Qscale_"<<n;
              h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
              h_test->StatOverflows(kTRUE);
              h_test->Sumw2(kTRUE);
              HistsSysReweightsQscale[i][k-2][0][n] = h_test;
              name.str("");
            }
            for (int n=0;n<nPdf;++n){
              name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_PDF_"<<n;
              h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
              h_test->StatOverflows(kTRUE);
              h_test->Sumw2(kTRUE);
              HistsSysReweightsPDF[i][k-2][0][n] = h_test;
              name.str("");
            }
            for (int n=0;n<nPS;++n){
              name<<channels[i]<<"_"<<regions[k]<<"_"<<it->first<<"_PS_"<<n;
              h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),it->second.at(1), it->second.at(2), it->second.at(3));
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
//  TFile file_out ("ANoutput.root","RECREATE");
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
  float BDToutput_;
  int ch_;

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
  tree_out.Branch("BDToutput"      , &BDToutput_, "BDToutput/F" ) ;
  tree_out.Branch("ch"      , &ch_, "ch/I" ) ;

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
//  readerMVA->BookMVA( "BDTG", "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/TMVAClassification_BDTG.weights.xml");
  readerMVA->BookMVA( "BDTG", "/afs/crc.nd.edu/user/r/rgoldouz/BNV/NanoAnalysis/input/TMVAClassification_BDTG_NEW.weights.xml");

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
  std::pair<double,double> METXYCorr;

  TLorentzVector wp, wm, b, ab, top, atop, lep, alep;
  std::vector<float> nominalWeights;
  TLorentzVector recoTop, recoBjet, recoW, recoNu, recoL1, recoL2, highPtMu;
  nominalWeights.assign(sys.size(), 1);
  std::vector<float> sysUpWeights;
  sysUpWeights.assign(sys.size(), 1);
  std::vector<float> sysDownWeights;
  sysDownWeights.assign(sys.size(), 1);
  bool leptonPass;
  bool triggerPass;
  bool DyPass;
  bool DyInZPass;
  bool MetPass;
  bool triggerPassEE;
  bool triggerPassEMu;
  bool triggerPassMuMu;
  bool metFilterPass;
  bool ifTopPt=false;
  float ttKFactor;
  int ch;
  float MetCut=60;
  float sf_Ele_Reco;
  float sf_Ele_ID;
  float sf_Mu_ID;
  float sf_Mu_ISO;
  float sf_Mu_RECO;
  float sf_Trigger;
  float sf_JetPuId;
  float weight_PU;
  float weight_Lumi;
  float weight_lep;
  float weight_lepB;
  float weight_EFT;
  float weight_prefiring;
  float weight_topPtPowhegData;
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
  float mT2output;
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
  int lbjet;
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
  float myMETpt;
  float myMETphi;
  int nLHEl;
  float mllLHE;
  int weightSign;
  float effTmc1;   
  float effTmc2;   
  float effTdata1; 
  float effTdata2; 
  float effTdata1Up;
  float effTdata2Up;
  float effTdata1Down;
  float effTdata2Down;
  std::vector<int> reg(regions.size());
  std::vector<float> wgt(regions.size());
  std::vector<WCFit> wcfit(regions.size());
  int JECnbUp;
  int JECnbDown;
  float nMet=0;
  float nMetUp=0;
  float nMetDown=0;

  if (fname.Contains("TTTo2L2Nu") || fname.Contains("sys") || fname.Contains("TTBNV")) ifTopPt=true;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  Long64_t ntr = fChain->GetEntries ();

//Loop over events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
//  for (Long64_t jentry=0; jentry<10000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    displayProgress(jentry, ntr) ;
    nLHEl=0;
    if(data == "mc" && !fname.Contains("pythia")){
      for (int l=0;l<nLHEPart;l++){
        if(abs(LHEPart_pdgId[l]) ==11 || abs(LHEPart_pdgId[l]) ==13 ){
          if(LHEPart_pt[l]>20 && abs(LHEPart_eta[l])<2.4) nLHEl++;
        }
      }
      for (int l=0;l<nGenPart;l++){
        if(abs(GenPart_pdgId[l])==11 || GenPart_pdgId[l]==13){
          if(abs(GenPart_pdgId[GenPart_genPartIdxMother[l]])==24 && GenPart_pt[l]>20 &&  abs(GenPart_eta[l])<2.4) nLHEl++;
        }
      }
    }

    if(nLHEl >1) nAccept++;

    triggerPassEE = false;
    triggerPassEMu = false;
    triggerPassMuMu = false;
    metFilterPass = false;
    leptonPass = false;
    triggerPass = false;
    DyPass = false;
    DyInZPass  = false;
    MetPass = false;
    ch =10;
    ttKFactor=1;
    sf_Ele_Reco =1;
    sf_Ele_ID =1;
    sf_Mu_ID =1;
    sf_Mu_ISO =1;
    sf_Mu_RECO =1;
    sf_Trigger =1;
    sf_JetPuId =1;
    muPtSFRochester=1;
    weight_PU =1;
    weight_Lumi =1;
    weight_lep =1;
    weight_lepB =1;
    weight_EFT =1;
    weight_prefiring =1;
    weight_topPtPowhegData =1;
    weight_topPtMGLO =1;
    P_bjet_data =1;
    P_bjet_mc =1;
    MVAoutput=0;
    mT2output=0;
    nbjetGen=0;
    nbjet=0;
    lbjet=0;
    nbjetJerUp=0;
    nbjetJerDown=0;
    nbjetJesUp=0;
    nbjetJesDown=0;
    metUnclusMETUp=0;
    metUnclusMETDown=0;
    metUnclusMETPhiUp=0;
    metUnclusMETPhiDown=0;
    if(year == "2016preVFP") METXYCorr = METXYCorr_Met_MetPhi(MET_pt, MET_phi, run, "2016APV", data == "mc", PV_npvs, true, false);
    else if(year == "2016postVFP") METXYCorr = METXYCorr_Met_MetPhi(MET_pt, MET_phi, run, "2016nonAPV", data == "mc", PV_npvs, true, false);
    else METXYCorr = METXYCorr_Met_MetPhi(MET_pt, MET_phi, run, year, data == "mc", PV_npvs, true, false);
//    myMETpt= METXYCorr.first;
//    myMETphi= METXYCorr.second;
    myMETpt=MET_pt;
    myMETphi=MET_phi;

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
    mllLHE=0;
    weightSign=1;
    effTmc1=1;
    effTmc2=1;
    effTdata1=1;
    effTdata2=1;
   
    if (data == "mc" && fname.Contains("DY50")){
      for (int l=0;l<nLHEPart;l++){
        if(LHEPart_pdgId[l]==11 || LHEPart_pdgId[l]==13 || LHEPart_pdgId[l]==15) lep.SetPtEtaPhiM(LHEPart_pt[l], LHEPart_eta[l], LHEPart_phi[l], LHEPart_mass[l]) ;
        if(LHEPart_pdgId[l]==-11 || LHEPart_pdgId[l]==-13 || LHEPart_pdgId[l]==-15) alep.SetPtEtaPhiM(LHEPart_pt[l], LHEPart_eta[l], LHEPart_phi[l], LHEPart_mass[l]) ;
      }
      mllLHE=(lep+alep).M();
      if(year == "2016preVFP" && mllLHE>200) continue;
      if(year != "2016preVFP" && year != "2016postVFP" && mllLHE>100) continue;
    }
 
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
      if ( Flag_goodVertices  &&  Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter &&  Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_eeBadScFilter && Flag_ecalBadCalibFilter && Flag_BadPFMuonDzFilter) metFilterPass = true;
    }
    else{
      if ( Flag_goodVertices  &&  Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter &&  Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_eeBadScFilter && Flag_BadPFMuonDzFilter) metFilterPass = true;
    }

//trigger MC
    if(data == "mc" && year == "2016preVFP"){
      if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele27_WPTight_Gsf || HLT_Photon175 ) triggerPassEE =true;
      if(HLT_Mu50 || HLT_TkMu50) triggerPassEMu =true;
      if(HLT_Mu50 || HLT_TkMu50) triggerPassMuMu =true;
    }

    if(data == "mc" && year == "2016postVFP"){
      if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele27_WPTight_Gsf || HLT_Photon175 ) triggerPassEE =true;
      if(HLT_Mu50 || HLT_TkMu50) triggerPassEMu =true;
      if(HLT_Mu50 || HLT_TkMu50) triggerPassMuMu =true;
    }
    if(data == "mc" && year == "2017"){
      if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele35_WPTight_Gsf || HLT_Photon200) triggerPassEE =true;
      if(HLT_Mu50 || HLT_TkMu100 || HLT_OldMu100) triggerPassEMu =true;
      if(HLT_Mu50 || HLT_TkMu100 || HLT_OldMu100) triggerPassMuMu =true;
    }

    if(data == "mc" && year == "2018"){
      if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele32_WPTight_Gsf || HLT_Photon200) triggerPassEE =true;
      if(HLT_Mu50 || HLT_TkMu100 || HLT_OldMu100) triggerPassEMu =true;
      if(HLT_Mu50 || HLT_TkMu100 || HLT_OldMu100) triggerPassMuMu =true;
    }

//trigger DATA
    if(data == "data"){
      if(year == "2016preVFP"){
        if(dataset=="SingleElectron"){
          if(!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && (HLT_Ele27_WPTight_Gsf || HLT_Photon175)) triggerPassEE =true;
        }
        if(dataset=="SingleMuon"){
          if(HLT_Mu50 || HLT_TkMu50) triggerPassMuMu =true;
          if(HLT_Mu50 || HLT_TkMu50) triggerPassEMu =true;
        }
        if(dataset=="DoubleEG"){
          if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) triggerPassEE =true;
        }
      }
      if(year == "2016postVFP"){
        if(dataset=="SingleElectron"){
          if(!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && (HLT_Ele27_WPTight_Gsf||HLT_Photon175)) triggerPassEE =true;
        }
        if(dataset=="SingleMuon"){
          if(HLT_Mu50 || HLT_TkMu50) triggerPassMuMu =true;
          if(HLT_Mu50 || HLT_TkMu50) triggerPassEMu =true;
        }
        if(dataset=="DoubleEG"){
          if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) triggerPassEE =true;
        }
      }
      if(year == "2017"){
        if(dataset=="SingleElectron"){
          if(!HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL && (HLT_Ele35_WPTight_Gsf || HLT_Photon200)) triggerPassEE =true;
        }
        if(dataset=="SingleMuon"){
          if(HLT_Mu50 || HLT_TkMu100 || HLT_OldMu100) triggerPassMuMu =true;
          if(HLT_Mu50 || HLT_TkMu100 || HLT_OldMu100) triggerPassEMu =true;
        }
        if(dataset=="DoubleEG"){
          if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL) triggerPassEE =true;
        }
      }
      if(year == "2018"){
        if(dataset=="EGamma"){
          if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele32_WPTight_Gsf || HLT_Photon200) triggerPassEE =true;
        }
        if(dataset=="SingleMuon"){
          if(HLT_Mu50 || HLT_TkMu100 || HLT_OldMu100) triggerPassMuMu =true;
          if(HLT_Mu50 || HLT_TkMu100 || HLT_OldMu100) triggerPassEMu =true;
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
      if(Electron_pt[l] > 35) eleEffDen->Fill(Electron_pt[l]);
//      if(Electron_cutBased[l] < 4) continue;
      if(!Electron_cutBased_HEEP[l]) continue;
//      if(!Electron_convVeto[l]) continue;
      //remove HEM affected region
      if (year == "2018" && Electron_eta[l] < -1.3 && Electron_phi[l] < -0.87 && Electron_phi[l] > -1.57) continue;
      if (data == "mc"){
        if((Electron_pt[l] + 0.004*Electron_pt[l]) > 35 && abs(Electron_eta[l])< 1.4442) selectedLeptonsEleScaleUp->push_back(new lepton_candidate(Electron_pt[l]+ 0.004*Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
        if((Electron_pt[l] + 0.008*Electron_pt[l]) > 35 && abs(Electron_eta[l])> 1.4442) selectedLeptonsEleScaleUp->push_back(new lepton_candidate(Electron_pt[l]+ 0.008*Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
        if((Electron_pt[l] - 0.004*Electron_pt[l]) > 35 && abs(Electron_eta[l])< 1.4442) selectedLeptonsEleScaleDown->push_back(new lepton_candidate(Electron_pt[l]- 0.004*Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
        if((Electron_pt[l] - 0.008*Electron_pt[l]) > 35 && abs(Electron_eta[l])> 1.4442) selectedLeptonsEleScaleDown->push_back(new lepton_candidate(Electron_pt[l]- 0.008*Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
      }
      if(Electron_pt[l] <35) continue;
      eleEffNum->Fill(Electron_pt[l]);
      selectedLeptons->push_back(new lepton_candidate(Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
      if (data == "mc"){
        selectedLeptonsMuScaleUp->push_back(new lepton_candidate(Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
        selectedLeptonsMuScaleDown->push_back(new lepton_candidate(Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
        selectedLeptonsMuResUp->push_back(new lepton_candidate(Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
        selectedLeptonsMuResDown->push_back(new lepton_candidate(Electron_pt[l],Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
        sf_Ele_Reco = sf_Ele_Reco * csetEleIdReco->evaluate({year, "sf", "RecoAbove20", Electron_eta[l],Electron_pt[l]}); 
       
        nominalWeights[getVecPos(sys,"eleRecoSf")] = nominalWeights[getVecPos(sys,"eleRecoSf")] * csetEleIdReco->evaluate({year, "sf", "RecoAbove20", Electron_eta[l],Electron_pt[l]});
        sysUpWeights[getVecPos(sys,"eleRecoSf")] = sysUpWeights[getVecPos(sys,"eleRecoSf")] * csetEleIdReco->evaluate({year, "sfup", "RecoAbove20", Electron_eta[l],Electron_pt[l]});
        sysDownWeights[getVecPos(sys,"eleRecoSf")] = sysDownWeights[getVecPos(sys,"eleRecoSf")] * csetEleIdReco->evaluate({year, "sfdown", "RecoAbove20", Electron_eta[l],Electron_pt[l]});

        if (abs(Electron_eta[l])< 1.4442) {
          sf_Ele_ID = sf_Ele_ID *HEEPSF_B[year];
          nominalWeights[getVecPos(sys,"eleIDSf")] = nominalWeights[getVecPos(sys,"eleIDSf")] * HEEPSF_B[year];
          if(Electron_pt[l] <100){
            sysUpWeights[getVecPos(sys,"eleIDSf")] = sysUpWeights[getVecPos(sys,"eleIDSf")] * HEEPSF_B[year] * 1.01;
            sysDownWeights[getVecPos(sys,"eleIDSf")] = sysDownWeights[getVecPos(sys,"eleIDSf")] * HEEPSF_B[year] * 0.99;
          }
          else{
            sysUpWeights[getVecPos(sys,"eleIDSf")] = sysUpWeights[getVecPos(sys,"eleIDSf")] * HEEPSF_B[year] * (1 + (0.01 + (0.02/900)*(Electron_pt[l]-100)));
            sysDownWeights[getVecPos(sys,"eleIDSf")] = sysDownWeights[getVecPos(sys,"eleIDSf")] * HEEPSF_B[year] * (1 - (0.01 + (0.02/900)*(Electron_pt[l]-100)));
          }
        }
        else {
          sf_Ele_ID = sf_Ele_ID * HEEPSF_E[year];
          nominalWeights[getVecPos(sys,"eleIDSf")] = nominalWeights[getVecPos(sys,"eleIDSf")] * HEEPSF_E[year];
          if(Electron_pt[l] <100){
            sysUpWeights[getVecPos(sys,"eleIDSf")] = sysUpWeights[getVecPos(sys,"eleIDSf")] * HEEPSF_B[year] * 1.02;
            sysDownWeights[getVecPos(sys,"eleIDSf")] = sysDownWeights[getVecPos(sys,"eleIDSf")] * HEEPSF_B[year] * 0.98;
          }
          else if(Electron_pt[l] <300){
            sysUpWeights[getVecPos(sys,"eleIDSf")] = sysUpWeights[getVecPos(sys,"eleIDSf")] * HEEPSF_B[year] * (1 + (0.02 + (0.03/200)*(Electron_pt[l]-100)));
            sysDownWeights[getVecPos(sys,"eleIDSf")] = sysDownWeights[getVecPos(sys,"eleIDSf")] * HEEPSF_B[year] * (1 - (0.02 + (0.03/200)*(Electron_pt[l]-100)));
          }
          else{
            sysUpWeights[getVecPos(sys,"eleIDSf")] = sysUpWeights[getVecPos(sys,"eleIDSf")] * HEEPSF_B[year] * 1.05;
            sysDownWeights[getVecPos(sys,"eleIDSf")] = sysDownWeights[getVecPos(sys,"eleIDSf")] * HEEPSF_B[year] * 0.95;
          }
        }
      }
    }
// Muon selection
    for (int l=0;l<nMuon;l++){
      if(abs(Muon_eta[l]) > 2.4) continue;
      if(Muon_highPtId[l]<2) continue;
      if(Muon_tkIsoId[l]<1) continue;
      muPtSFRochester =1;
      if (data == "data" && Muon_tunepRelPt[l]*Muon_pt[l]>53) selectedLeptons->push_back(new lepton_candidate(Muon_tunepRelPt[l]*Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
        if(data == "mc"){
          if (Muon_tunepRelPt[l]*Muon_pt[l]<53) continue;
          selectedLeptons->push_back(new lepton_candidate(Muon_tunepRelPt[l]*Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          pt_res = 0.02;
          if(abs(Muon_eta[l]) < 1.2) pt_res = 0.01;
          SmearedMuonPt = Tr.Gaus(0, pt_res);
          correctMuonPt = GE->GEScaleCorrPt(1600, Muon_tunepRelPt[l]*Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]);
          correctMuonPtUp = GE->GEScaleCorrPt(1601, Muon_tunepRelPt[l]*Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]);
          correctMuonPtDown = GE->GEScaleCorrPt(1602, Muon_tunepRelPt[l]*Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l]);
          selectedLeptonsEleScaleUp->push_back(new lepton_candidate(Muon_tunepRelPt[l]* Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          selectedLeptonsEleScaleDown->push_back(new lepton_candidate(Muon_tunepRelPt[l]*Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          if((1+abs(SmearedMuonPt))*Muon_tunepRelPt[l]*Muon_pt[l]>53) selectedLeptonsMuResUp->push_back(new lepton_candidate((1+abs(SmearedMuonPt))*Muon_tunepRelPt[l]*Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          if((1-abs(SmearedMuonPt))*Muon_tunepRelPt[l]*Muon_pt[l]>53) selectedLeptonsMuResDown->push_back(new lepton_candidate((1-abs(SmearedMuonPt))*Muon_tunepRelPt[l]*Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          if(Muon_tunepRelPt[l]*Muon_pt[l]*(1+abs((correctMuonPtUp - correctMuonPt)/correctMuonPt))>53) selectedLeptonsMuScaleUp->push_back(new lepton_candidate(Muon_tunepRelPt[l]*Muon_pt[l]*(1+abs((correctMuonPtUp - correctMuonPt)/correctMuonPt)),Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          if(Muon_tunepRelPt[l]*Muon_pt[l]*(1-abs((correctMuonPtUp - correctMuonPt)/correctMuonPt))>53) selectedLeptonsMuScaleDown->push_back(new lepton_candidate(Muon_tunepRelPt[l]*Muon_pt[l]*(1-abs((correctMuonPtUp - correctMuonPt)/correctMuonPt)) ,Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          highPtMu.SetPtEtaPhiM(Muon_tunepRelPt[l]*Muon_pt[l], Muon_eta[l],Muon_phi[l], 0.10566) ;
          sf_Mu_RECO=scale_factor(&highPtMuRecoSF_pVsAbsEta_H, highPtMu.E(), abs(highPtMu.Eta()),"central",false, true);
          sf_Mu_ID = sf_Mu_ID * csetMuTightId->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_tunepRelPt[l]*Muon_pt[l], "sf"})*sf_Mu_RECO;
          sf_Mu_ISO = sf_Mu_ISO * csetMuTightRelIso->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_tunepRelPt[l]*Muon_pt[l], "sf"});
          nominalWeights[getVecPos(sys,"muIdIsoSf")] = nominalWeights[getVecPos(sys,"muIdIsoSf")] * csetMuTightRelIso->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_tunepRelPt[l]*Muon_pt[l], "sf"}) * csetMuTightId->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_tunepRelPt[l]*Muon_pt[l], "sf"})*sf_Mu_RECO;
          sysUpWeights[getVecPos(sys,"muIdIsoSf")] = sysUpWeights[getVecPos(sys,"muIdIsoSf")] * csetMuTightRelIso->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_tunepRelPt[l]*Muon_pt[l], "systup"}) * csetMuTightId->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_tunepRelPt[l]*Muon_pt[l], "systup"})*sf_Mu_RECO;
          sysDownWeights[getVecPos(sys,"muIdIsoSf")] = sysDownWeights[getVecPos(sys,"muIdIsoSf")] * csetMuTightRelIso->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_tunepRelPt[l]*Muon_pt[l], "systdown"}) * csetMuTightId->evaluate({year + "_UL", abs(Muon_eta[l]),  Muon_tunepRelPt[l]*Muon_pt[l], "systdown"})*sf_Mu_RECO;
        }
//correct MET if the muon is not a PF muon
      if(!Muon_isPFcand) myMETpt = sqrt(pow(myMETpt*cos(myMETphi) - Muon_tunepRelPt[l]*Muon_pt[l]*cos(Muon_phi[l]),2) + pow(myMETpt*sin(myMETphi) - Muon_tunepRelPt[l]*Muon_pt[l]*sin(Muon_phi[l]),2));
      if(Muon_isPFcand &&  Muon_tunepRelPt[l]!=1) myMETpt = sqrt(pow(myMETpt*cos(myMETphi) + ((1-Muon_tunepRelPt[l])*Muon_pt[l]*cos(Muon_phi[l])),2) + pow(myMETpt*sin(myMETphi) + ((1-Muon_tunepRelPt[l])*Muon_pt[l]*sin(Muon_phi[l])),2));
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
//        if(deltaR((*selectedLeptons)[i]->eta_,(*selectedLeptons)[i]->phi_,Jet_eta[l],Jet_phi[l]) < 0.4 && (*selectedLeptons)[i]->lep_ ==10) cout<<" A jet with pt "<< Jet_pt_nom[l] << " is "<<deltaR((*selectedLeptons)[i]->eta_,(*selectedLeptons)[i]->phi_,Jet_eta[l],Jet_phi[l]) <<" close to a muon with pt "<< (*selectedLeptons)[i]->pt_<<" and jet rel pt "<< Muon_jetPtRelv2[i]<<" and trkiso "<<Muon_tkRelIso[(*selectedLeptons)[i]->indice_]<<" isbtagged:"<<(Jet_btagDeepFlavB[l] > 0.2783)<<endl;
      }
//      if(jetlepfail) continue;
      if(data == "mc" && abs(Jet_eta[l]) < 2.4){
        if(Jet_pt_nom[l] >30) {
          selectedJets->push_back(new jet_candidate(Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l], year,Jet_partonFlavour[l]));
          if(Jet_pt_nom[l] <50){
            sf_JetPuId = sf_JetPuId * csetJetPuID->evaluate({Jet_eta[l],Jet_pt_nom[l],"nom","L"});
            nominalWeights[getVecPos(sys,"JetPuID")] = nominalWeights[getVecPos(sys,"JetPuID")] * csetJetPuID->evaluate({Jet_eta[l],Jet_pt_nom[l],"nom","L"});
            sysUpWeights[getVecPos(sys,"JetPuID")] = sysUpWeights[getVecPos(sys,"JetPuID")] * csetJetPuID->evaluate({Jet_eta[l],Jet_pt_nom[l],"up","L"});
            sysDownWeights[getVecPos(sys,"JetPuID")] = sysDownWeights[getVecPos(sys,"JetPuID")] * csetJetPuID->evaluate({Jet_eta[l],Jet_pt_nom[l],"down","L"});
          }
        }
        uncJes->setJetPt(Jet_pt_nom[l]);
        uncJes->setJetEta(Jet_eta[l]);
        sup = abs(uncJes->getUncertainty(true));
        if ((1+sup)*Jet_pt_nom[l]>30) {
          selectedJetsJesUp->push_back(new jet_candidate((1+sup)*Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l], year,Jet_partonFlavour[l]));
          MetXJetsJesUp = MetXJetsJesUp + (sup * Jet_pt_nom[l] * cos(Jet_phi[l]));
          MetYJetsJesUp = MetYJetsJesUp + (sup * Jet_pt_nom[l] * sin(Jet_phi[l]));
        }

        uncJes->setJetPt(Jet_pt_nom[l]);
        uncJes->setJetEta(Jet_eta[l]);
        sdw = abs(uncJes->getUncertainty(false));
        if ((1-sdw)*Jet_pt_nom[l]>30) {
          selectedJetsJesDown->push_back(new jet_candidate((1-sdw)*Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l], year,Jet_partonFlavour[l]));
          MetXJetsJesDown = MetXJetsJesDown - (sdw * Jet_pt_nom[l] * cos(Jet_phi[l]));
          MetYJetsJesDown = MetYJetsJesDown - (sdw * Jet_pt_nom[l] * sin(Jet_phi[l]));
        }

        uncJer->setJetPt(Jet_pt_nom[l]);
        uncJer->setJetEta(Jet_eta[l]);
        sup = abs(uncJer->getUncertainty(true));
        if ((1+sup)*Jet_pt_nom[l]>30) {
          selectedJetsJerUp->push_back(new jet_candidate((1+sup)*Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l], year,Jet_partonFlavour[l]));
          MetXJetsJerUp = MetXJetsJerUp + (sup * Jet_pt_nom[l] * cos(Jet_phi[l]));
          MetYJetsJerUp = MetYJetsJerUp + (sup * Jet_pt_nom[l] * sin(Jet_phi[l]));
        }

        uncJer->setJetPt(Jet_pt_nom[l]);
        uncJer->setJetEta(Jet_eta[l]);
        sdw = abs(uncJer->getUncertainty(false));
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

    MetJetsPtJesUp = sqrt(pow(myMETpt * cos(myMETphi) - MetXJetsJesUp,2)+pow(myMETpt * sin(myMETphi) - MetYJetsJesUp,2));
    MetJetsPhiJesUp = atan2 ((myMETpt * sin(myMETphi) - MetYJetsJesUp),(myMETpt * cos(myMETphi) - MetXJetsJesUp));
    MetJetsPtJesDown = sqrt(pow(myMETpt * cos(myMETphi) - MetXJetsJesDown,2)+pow(myMETpt * sin(myMETphi) - MetYJetsJesDown,2));
    MetJetsPhiJesDown = atan2 ((myMETpt * sin(myMETphi) - MetYJetsJesDown),(myMETpt * cos(myMETphi) - MetXJetsJesDown));
    MetJetsPtJerUp = sqrt(pow(myMETpt * cos(myMETphi) - MetXJetsJerUp,2)+pow(myMETpt * sin(myMETphi) - MetYJetsJerUp,2));
    MetJetsPhiJerUp = atan2 ((myMETpt * sin(myMETphi) - MetYJetsJerUp),(myMETpt * cos(myMETphi) - MetXJetsJerUp));
    MetJetsPtJerDown = sqrt(pow(myMETpt * cos(myMETphi) - MetXJetsJerDown,2)+pow(myMETpt * sin(myMETphi) - MetYJetsJerDown,2));
    MetJetsPhiJerDown = atan2 ((myMETpt * sin(myMETphi) - MetYJetsJerDown),(myMETpt * cos(myMETphi) - MetXJetsJerDown));

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
//cout<<Jet_jetId[l]<<":"<<(Jet_jetId[l]>>0&1)<<":"<<(Jet_jetId[l]>>1&1)<<":"<<(Jet_jetId[l]>>2&1)<<endl;
          if(Jet_jetId[l]<6) continue;
          if (jetVetoMaps_H.GetBinContent(jetVetoMaps_H.GetXaxis()->FindBin(Jet_eta[l]),jetVetoMaps_H.GetYaxis()->FindBin(Jet_phi[l]))>0) continue;
          jetlepfail = false;
          for (int i=0;i<selectedLeptons->size();i++){
            if(deltaR((*selectedLeptons)[i]->eta_,(*selectedLeptons)[i]->phi_,Jet_eta[l],Jet_phi[l]) < 0.4 ) jetlepfail=true;
          }
//          if(jetlepfail) continue;
//          JetCorrectionUncertainty *unc = vsrc[n];
          vsrc[n]->setJetPt(Jet_pt_nom[l]);
          vsrc[n]->setJetEta(Jet_eta[l]);
          sup = abs(vsrc[n]->getUncertainty(true));
          if ((1+sup)*Jet_pt_nom[l]>30) {
            JECJetsUp->push_back(new jet_candidate((1+sup)*Jet_pt_nom[l],Jet_eta[l],Jet_phi[l],Jet_mass[l],Jet_btagDeepFlavB[l], year,Jet_partonFlavour[l]));
            JECMETUpx = JECMETUpx + (sup * Jet_pt_nom[l] * cos(Jet_phi[l]));
            JECMETUpy = JECMETUpy + (sup * Jet_pt_nom[l] * sin(Jet_phi[l]));
          }
          vsrc[n]->setJetPt(Jet_pt_nom[l]);
          vsrc[n]->setJetEta(Jet_eta[l]);
          sdw = abs(vsrc[n]->getUncertainty(false));
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
        JECsysMETUp->push_back(sqrt(pow(myMETpt * cos(myMETphi) - JECMETUpx,2)+pow(myMETpt * sin(myMETphi) - JECMETUpy,2)));
        JECsysMETDown->push_back(sqrt(pow(myMETpt * cos(myMETphi) - JECMETDownx,2)+pow(myMETpt * sin(myMETphi) - JECMETDowny,2)));
        JECsysMETPhiUp->push_back(atan2 ((myMETpt * sin(myMETphi) - JECMETUpx),(myMETpt * cos(myMETphi) - JECMETUpy)));
        JECsysMETPhiDown->push_back(atan2 ((myMETpt * sin(myMETphi) - JECMETDownx),(myMETpt * cos(myMETphi) - JECMETDowny)));
      }
    
      for (int n=0;n<sysJecNames.size();++n){
        JECnbUp = 0;
        JECnbDown = 0;
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
      if((*selectedJets)[l]->btag_ && lbjet==0) lbjet=l;
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
          nominalWeights[getVecPos(sys,"bcTagSf")] = nominalWeights[getVecPos(sys,"bcTagSf")] * BJetEff * BJetSF;
          sysUpWeights[getVecPos(sys,"bcTagSf")] = sysUpWeights[getVecPos(sys,"bcTagSf")] * BJetEff * BJetSF_UpCorr;
          sysDownWeights[getVecPos(sys,"bcTagSf")] = sysDownWeights[getVecPos(sys,"bcTagSf")] * BJetEff * BJetSF_DownCorr;
          nominalWeights[getVecPos(sys,"LTagSf")] = nominalWeights[getVecPos(sys,"LTagSf")] * BJetEff * BJetSF;
          sysUpWeights[getVecPos(sys,"LTagSf")] = sysUpWeights[getVecPos(sys,"LTagSf")] * BJetEff * BJetSF;
          sysDownWeights[getVecPos(sys,"LTagSf")] = sysDownWeights[getVecPos(sys,"LTagSf")] * BJetEff * BJetSF;

          nominalWeights[getVecPos(sys,"bcTagSfUnCorr")] = nominalWeights[getVecPos(sys,"bcTagSfUnCorr")] * BJetEff * BJetSF;
          sysUpWeights[getVecPos(sys,"bcTagSfUnCorr")] = sysUpWeights[getVecPos(sys,"bcTagSfUnCorr")] * BJetEff * BJetSF_UpUnCorr;
          sysDownWeights[getVecPos(sys,"bcTagSfUnCorr")] = sysDownWeights[getVecPos(sys,"bcTagSfUnCorr")] * BJetEff * BJetSF_DownUnCorr;
          nominalWeights[getVecPos(sys,"LTagSfUnCorr")] = nominalWeights[getVecPos(sys,"LTagSfUnCorr")] * BJetEff * BJetSF;
          sysUpWeights[getVecPos(sys,"LTagSfUnCorr")] = sysUpWeights[getVecPos(sys,"LTagSfUnCorr")] * BJetEff * BJetSF;
          sysDownWeights[getVecPos(sys,"LTagSfUnCorr")] = sysDownWeights[getVecPos(sys,"LTagSfUnCorr")] * BJetEff * BJetSF;
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - BJetEff);
          P_bjet_data = P_bjet_data * (1- (BJetEff * BJetSF));
          nominalWeights[getVecPos(sys,"bcTagSf")] = nominalWeights[getVecPos(sys,"bcTagSf")]* (1- (BJetEff * BJetSF));
          sysUpWeights[getVecPos(sys,"bcTagSf")] = sysUpWeights[getVecPos(sys,"bcTagSf")]* (1- (BJetEff * BJetSF_UpCorr));
          sysDownWeights[getVecPos(sys,"bcTagSf")] = sysDownWeights[getVecPos(sys,"bcTagSf")]* (1- (BJetEff * BJetSF_DownCorr));
          nominalWeights[getVecPos(sys,"LTagSf")] = nominalWeights[getVecPos(sys,"LTagSf")]* (1- (BJetEff * BJetSF));
          sysUpWeights[getVecPos(sys,"LTagSf")] = sysUpWeights[getVecPos(sys,"LTagSf")]* (1- (BJetEff * BJetSF));
          sysDownWeights[getVecPos(sys,"LTagSf")] = sysDownWeights[getVecPos(sys,"LTagSf")]* (1- (BJetEff * BJetSF));

          nominalWeights[getVecPos(sys,"bcTagSfUnCorr")] = nominalWeights[getVecPos(sys,"bcTagSfUnCorr")]* (1- (BJetEff * BJetSF));
          sysUpWeights[getVecPos(sys,"bcTagSfUnCorr")] = sysUpWeights[getVecPos(sys,"bcTagSfUnCorr")]* (1- (BJetEff * BJetSF_UpUnCorr));
          sysDownWeights[getVecPos(sys,"bcTagSfUnCorr")] = sysDownWeights[getVecPos(sys,"bcTagSfUnCorr")]* (1- (BJetEff * BJetSF_DownUnCorr));
          nominalWeights[getVecPos(sys,"LTagSfUnCorr")] = nominalWeights[getVecPos(sys,"LTagSfUnCorr")]* (1- (BJetEff * BJetSF));
          sysUpWeights[getVecPos(sys,"LTagSfUnCorr")] = sysUpWeights[getVecPos(sys,"LTagSfUnCorr")]* (1- (BJetEff * BJetSF));
          sysDownWeights[getVecPos(sys,"LTagSfUnCorr")] = sysDownWeights[getVecPos(sys,"LTagSfUnCorr")]* (1- (BJetEff * BJetSF));
        }
      }
//c-quark
      if( abs((*selectedJets)[l]->flavor_) == 4){
        h2_BTaggingEff_Denom_c->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
        if( (*selectedJets)[l]->btag_) {
          h2_BTaggingEff_Num_c->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
          P_bjet_mc = P_bjet_mc * CJetEff;
          P_bjet_data = P_bjet_data * CJetEff * CJetSF;
          nominalWeights[getVecPos(sys,"bcTagSf")] = nominalWeights[getVecPos(sys,"bcTagSf")] * CJetEff * CJetSF;
          sysUpWeights[getVecPos(sys,"bcTagSf")] = sysUpWeights[getVecPos(sys,"bcTagSf")] * CJetEff * CJetSF_UpCorr;
          sysDownWeights[getVecPos(sys,"bcTagSf")] = sysDownWeights[getVecPos(sys,"bcTagSf")] * CJetEff * CJetSF_DownCorr;
          nominalWeights[getVecPos(sys,"LTagSf")] = nominalWeights[getVecPos(sys,"LTagSf")] * CJetEff * CJetSF;
          sysUpWeights[getVecPos(sys,"LTagSf")] = sysUpWeights[getVecPos(sys,"LTagSf")] * CJetEff * CJetSF;
          sysDownWeights[getVecPos(sys,"LTagSf")] = sysDownWeights[getVecPos(sys,"LTagSf")] * CJetEff * CJetSF;

          nominalWeights[getVecPos(sys,"bcTagSfUnCorr")] = nominalWeights[getVecPos(sys,"bcTagSfUnCorr")] * CJetEff * CJetSF;
          sysUpWeights[getVecPos(sys,"bcTagSfUnCorr")] = sysUpWeights[getVecPos(sys,"bcTagSfUnCorr")] * CJetEff * CJetSF_UpUnCorr;
          sysDownWeights[getVecPos(sys,"bcTagSfUnCorr")] = sysDownWeights[getVecPos(sys,"bcTagSfUnCorr")] * CJetEff * CJetSF_DownUnCorr;
          nominalWeights[getVecPos(sys,"LTagSfUnCorr")] = nominalWeights[getVecPos(sys,"LTagSfUnCorr")] * CJetEff * CJetSF;
          sysUpWeights[getVecPos(sys,"LTagSfUnCorr")] = sysUpWeights[getVecPos(sys,"LTagSfUnCorr")] * CJetEff * CJetSF;
          sysDownWeights[getVecPos(sys,"LTagSfUnCorr")] = sysDownWeights[getVecPos(sys,"LTagSfUnCorr")] * CJetEff * CJetSF;
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - CJetEff);
          P_bjet_data = P_bjet_data * (1- (CJetEff * CJetSF));
          nominalWeights[getVecPos(sys,"bcTagSf")] = nominalWeights[getVecPos(sys,"bcTagSf")]* (1- (CJetEff * CJetSF));
          sysUpWeights[getVecPos(sys,"bcTagSf")] = sysUpWeights[getVecPos(sys,"bcTagSf")]* (1- (CJetEff * CJetSF_UpCorr));
          sysDownWeights[getVecPos(sys,"bcTagSf")] = sysDownWeights[getVecPos(sys,"bcTagSf")]* (1- (CJetEff * CJetSF_DownCorr));
          nominalWeights[getVecPos(sys,"LTagSf")] = nominalWeights[getVecPos(sys,"LTagSf")]* (1- (CJetEff * CJetSF));
          sysUpWeights[getVecPos(sys,"LTagSf")] = sysUpWeights[getVecPos(sys,"LTagSf")]* (1- (CJetEff * CJetSF));
          sysDownWeights[getVecPos(sys,"LTagSf")] = sysDownWeights[getVecPos(sys,"LTagSf")]* (1- (CJetEff * CJetSF));

          nominalWeights[getVecPos(sys,"bcTagSfUnCorr")] = nominalWeights[getVecPos(sys,"bcTagSfUnCorr")]* (1- (CJetEff * CJetSF));
          sysUpWeights[getVecPos(sys,"bcTagSfUnCorr")] = sysUpWeights[getVecPos(sys,"bcTagSfUnCorr")]* (1- (CJetEff * CJetSF_UpUnCorr));
          sysDownWeights[getVecPos(sys,"bcTagSfUnCorr")] = sysDownWeights[getVecPos(sys,"bcTagSfUnCorr")]* (1- (CJetEff * CJetSF_DownUnCorr));
          nominalWeights[getVecPos(sys,"LTagSfUnCorr")] = nominalWeights[getVecPos(sys,"LTagSfUnCorr")]* (1- (CJetEff * CJetSF));
          sysUpWeights[getVecPos(sys,"LTagSfUnCorr")] = sysUpWeights[getVecPos(sys,"LTagSfUnCorr")]* (1- (CJetEff * CJetSF));
          sysDownWeights[getVecPos(sys,"LTagSfUnCorr")] = sysDownWeights[getVecPos(sys,"LTagSfUnCorr")]* (1- (CJetEff * CJetSF));
        }
      }
//light-quark
      if( abs((*selectedJets)[l]->flavor_) != 4 && abs((*selectedJets)[l]->flavor_) != 5){
        h2_BTaggingEff_Denom_udsg->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
        if( (*selectedJets)[l]->btag_) {
          h2_BTaggingEff_Num_udsg->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
          P_bjet_mc = P_bjet_mc * LJetEff;
          P_bjet_data = P_bjet_data * LJetEff * LJetSF;
          nominalWeights[getVecPos(sys,"bcTagSf")] = nominalWeights[getVecPos(sys,"bcTagSf")]* LJetEff * LJetSF;
          sysUpWeights[getVecPos(sys,"bcTagSf")] = sysUpWeights[getVecPos(sys,"bcTagSf")]* LJetEff * LJetSF;
          sysDownWeights[getVecPos(sys,"bcTagSf")] = sysDownWeights[getVecPos(sys,"bcTagSf")]* LJetEff * LJetSF;
          nominalWeights[getVecPos(sys,"LTagSf")] = nominalWeights[getVecPos(sys,"LTagSf")] * LJetEff * LJetSF;
          sysUpWeights[getVecPos(sys,"LTagSf")] = sysUpWeights[getVecPos(sys,"LTagSf")] * LJetEff * LJetSF_UpCorr;
          sysDownWeights[getVecPos(sys,"LTagSf")] = sysDownWeights[getVecPos(sys,"LTagSf")] * LJetEff * LJetSF_DownCorr;

          nominalWeights[getVecPos(sys,"bcTagSfUnCorr")] = nominalWeights[getVecPos(sys,"bcTagSfUnCorr")]* LJetEff * LJetSF;
          sysUpWeights[getVecPos(sys,"bcTagSfUnCorr")] = sysUpWeights[getVecPos(sys,"bcTagSfUnCorr")]* LJetEff * LJetSF;
          sysDownWeights[getVecPos(sys,"bcTagSfUnCorr")] = sysDownWeights[getVecPos(sys,"bcTagSfUnCorr")]* LJetEff * LJetSF;
          nominalWeights[getVecPos(sys,"LTagSfUnCorr")] = nominalWeights[getVecPos(sys,"LTagSfUnCorr")] * LJetEff * LJetSF;
          sysUpWeights[getVecPos(sys,"LTagSfUnCorr")] = sysUpWeights[getVecPos(sys,"LTagSfUnCorr")] * LJetEff * LJetSF_UpUnCorr;
          sysDownWeights[getVecPos(sys,"LTagSfUnCorr")] = sysDownWeights[getVecPos(sys,"LTagSfUnCorr")] * LJetEff * LJetSF_DownUnCorr;
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - LJetEff);
          P_bjet_data = P_bjet_data * (1- (LJetEff * LJetSF));
          nominalWeights[getVecPos(sys,"bcTagSf")] = nominalWeights[getVecPos(sys,"bcTagSf")]* (1- (LJetEff * LJetSF));
          sysUpWeights[getVecPos(sys,"bcTagSf")] = sysUpWeights[getVecPos(sys,"bcTagSf")]* (1- (LJetEff * LJetSF));
          sysDownWeights[getVecPos(sys,"bcTagSf")] = sysDownWeights[getVecPos(sys,"bcTagSf")]* (1- (LJetEff * LJetSF));
          nominalWeights[getVecPos(sys,"LTagSf")] = nominalWeights[getVecPos(sys,"LTagSf")]* (1- (LJetEff * LJetSF));
          sysUpWeights[getVecPos(sys,"LTagSf")] = sysUpWeights[getVecPos(sys,"LTagSf")]* (1- (LJetEff * LJetSF_UpCorr));
          sysDownWeights[getVecPos(sys,"LTagSf")] = sysDownWeights[getVecPos(sys,"LTagSf")]* (1- (LJetEff * LJetSF_DownCorr));

          nominalWeights[getVecPos(sys,"bcTagSfUnCorr")] = nominalWeights[getVecPos(sys,"bcTagSfUnCorr")]* (1- (LJetEff * LJetSF));
          sysUpWeights[getVecPos(sys,"bcTagSfUnCorr")] = sysUpWeights[getVecPos(sys,"bcTagSfUnCorr")]* (1- (LJetEff * LJetSF));
          sysDownWeights[getVecPos(sys,"bcTagSfUnCorr")] = sysDownWeights[getVecPos(sys,"bcTagSfUnCorr")]* (1- (LJetEff * LJetSF));
          nominalWeights[getVecPos(sys,"LTagSfUnCorr")] = nominalWeights[getVecPos(sys,"LTagSfUnCorr")]* (1- (LJetEff * LJetSF));
          sysUpWeights[getVecPos(sys,"LTagSfUnCorr")] = sysUpWeights[getVecPos(sys,"LTagSfUnCorr")]* (1- (LJetEff * LJetSF_UpUnCorr));
          sysDownWeights[getVecPos(sys,"LTagSfUnCorr")] = sysDownWeights[getVecPos(sys,"LTagSfUnCorr")]* (1- (LJetEff * LJetSF_DownUnCorr));
        }
      }
    }

//MET Unc 
    metUnclusMETUp = sqrt(pow(myMETpt * cos(myMETphi) + MET_MetUnclustEnUpDeltaX,2)+pow(myMETpt * sin(myMETphi) + MET_MetUnclustEnUpDeltaY,2));
    metUnclusMETDown = sqrt(pow(myMETpt * cos(myMETphi) - MET_MetUnclustEnUpDeltaX,2)+pow(myMETpt * sin(myMETphi) - MET_MetUnclustEnUpDeltaY,2));
    metUnclusMETPhiUp = atan2 ((myMETpt * sin(myMETphi) + MET_MetUnclustEnUpDeltaY),(myMETpt * cos(myMETphi) + MET_MetUnclustEnUpDeltaX));
    metUnclusMETPhiDown = atan2 ((myMETpt * sin(myMETphi) - MET_MetUnclustEnUpDeltaY),(myMETpt * cos(myMETphi) - MET_MetUnclustEnUpDeltaX));
//PU reweighting
    if (data == "mc" && year == "2016preVFP") {
      weight_PU = wPU.PU_2016preVFP(int(Pileup_nTrueInt),"nominal");
      nominalWeights[getVecPos(sys,"pu")] = wPU.PU_2016preVFP(int(Pileup_nTrueInt),"nominal");
      sysUpWeights[getVecPos(sys,"pu")] = wPU.PU_2016preVFP(int(Pileup_nTrueInt),"up");
      sysDownWeights[getVecPos(sys,"pu")] = wPU.PU_2016preVFP(int(Pileup_nTrueInt),"down");
    }
    if (data == "mc" && year == "2016postVFP") {
      weight_PU = wPU.PU_2016postVFP(int(Pileup_nTrueInt),"nominal");
      nominalWeights[getVecPos(sys,"pu")] = wPU.PU_2016postVFP(int(Pileup_nTrueInt),"nominal");
      sysUpWeights[getVecPos(sys,"pu")] = wPU.PU_2016postVFP(int(Pileup_nTrueInt),"up");
      sysDownWeights[getVecPos(sys,"pu")] = wPU.PU_2016postVFP(int(Pileup_nTrueInt),"down");
    }
    if (data == "mc" && year == "2017") {
      weight_PU = wPU.PU_2017(int(Pileup_nTrueInt),"nominal");
      nominalWeights[getVecPos(sys,"pu")] = wPU.PU_2017(int(Pileup_nTrueInt),"nominal");
      sysUpWeights[getVecPos(sys,"pu")] = wPU.PU_2017(int(Pileup_nTrueInt),"up");
      sysDownWeights[getVecPos(sys,"pu")] = wPU.PU_2017(int(Pileup_nTrueInt),"down");
    }
    if (data == "mc" && year == "2018") {
      weight_PU = wPU.PU_2018(int(Pileup_nTrueInt),"nominal");
      nominalWeights[getVecPos(sys,"pu")] = wPU.PU_2018(int(Pileup_nTrueInt),"nominal");
      sysUpWeights[getVecPos(sys,"pu")] = wPU.PU_2018(int(Pileup_nTrueInt),"up");
      sysDownWeights[getVecPos(sys,"pu")] = wPU.PU_2018(int(Pileup_nTrueInt),"down");
    }
    if (data == "mc") weight_Lumi = (1000*xs*lumi)/Nevent;
    if (data == "mc"){
        weight_prefiring = L1PreFiringWeight_Nom;
        nominalWeights[getVecPos(sys,"prefiring")] = L1PreFiringWeight_Nom;
        sysUpWeights[getVecPos(sys,"prefiring")] = L1PreFiringWeight_Up;
        sysDownWeights[getVecPos(sys,"prefiring")] = L1PreFiringWeight_Dn;
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
//      weight_topPtPowhegData = sqrt(topPtPowhegData((wp + b).Pt()) * topPtPowhegData((wm + ab).Pt()));
      weight_topPtPowhegData = sqrt(topPtPowhegMC((wp + b).Pt()) * topPtPowhegMC((wm + ab).Pt()));
      weight_topPtMGLO = sqrt(topPtMGLO((atop).Pt()) * topPtMGLO((top).Pt()));
//      weight_topPtPowhegData = weight_topPtMGLO;
      nominalWeights[getVecPos(sys,"topPt")] = sqrt(topPtPowhegMC((wp + b).Pt()) * topPtPowhegMC((wm + ab).Pt()));
      sysUpWeights[getVecPos(sys,"topPt")]   = sqrt(topPtPowhegMC((wp + b).Pt()) * topPtPowhegMC((wm + ab).Pt()));
      sysDownWeights[getVecPos(sys,"topPt")] = sqrt(topPtPowhegMC((wp + b).Pt()) * topPtPowhegMC((wm + ab).Pt()));
      if(sqrt((wp + b).Pt() * (wm + ab).Pt())>300){
      sysUpWeights[getVecPos(sys,"topPt")]   = sqrt(topPtPowhegMC((wp + b).Pt()) * topPtPowhegMC((wm + ab).Pt()))*(1+GetGraphYError(TTptQS_TG, sqrt((wp + b).Pt() * (wm + ab).Pt())));
      sysDownWeights[getVecPos(sys,"topPt")] = sqrt(topPtPowhegMC((wp + b).Pt()) * topPtPowhegMC((wm + ab).Pt()))*(1-GetGraphYError(TTptQS_TG, sqrt((wp + b).Pt() * (wm + ab).Pt())));
      }
    }

    if (fname.Contains("TTBNV")) {
      weight_topPtPowhegData = weight_topPtMGLO;
      ttKFactor=831.7/445.0;
    }

//Check if analysis cuts are passed and then categorize dilepton channels
    if(selectedLeptons->size() ==2){ 
      if (((*selectedLeptons)[0]->pt_ > 25) && ((*selectedLeptons)[0]->charge_ * (*selectedLeptons)[1]->charge_ == -1) && ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()>20 && ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()<1200) leptonPass=true;  
      if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 2 && triggerPassEE) ch = 0;
      if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 11 && triggerPassEMu) ch = 1;
      if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 20 && triggerPassMuMu) ch = 2;
      if (ch ==0 && triggerPassEE) triggerPass=true;
      if (ch ==1 && triggerPassEMu) triggerPass=true;
      if (ch ==2 && triggerPassMuMu) triggerPass=true;
      if (((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()>106) DyPass = true;
      if (((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()<106 && ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()>76) DyInZPass=true;
      if (myMETpt > MetCut) MetPass = true;
      if (data == "mc" && ch==0) {
        sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"central",false, true);
        nominalWeights[getVecPos(sys,"trigSF")] = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"central",false, true);
        sysUpWeights[getVecPos(sys,"trigSF")] = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"up",false, true);
        sysDownWeights[getVecPos(sys,"trigSF")] = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"down",false, true);
      }
      if (data == "mc" && ch==1) {
          if((*selectedLeptons)[0]->lep_==10) {
            effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptons)[0]->pt_, abs((*selectedLeptons)[0]->eta_),"central",false, true);
            effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptons)[0]->pt_, abs((*selectedLeptons)[0]->eta_),"central",false, true);
            effTdata1Up = effTdata1 + scale_factor(&eff_triggermumu_dataUp_H, (*selectedLeptons)[0]->pt_, abs((*selectedLeptons)[0]->eta_),"central",false, true);
            effTdata1Down = effTdata1 - scale_factor(&eff_triggermumu_dataDown_H, (*selectedLeptons)[0]->pt_, abs((*selectedLeptons)[0]->eta_),"central",false, true);
          }
          else{
            effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptons)[1]->pt_, abs((*selectedLeptons)[1]->eta_),"central",false, true);
            effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptons)[1]->pt_, abs((*selectedLeptons)[1]->eta_),"central",false, true);
            effTdata1Up = effTdata1 + scale_factor(&eff_triggermumu_dataUp_H, (*selectedLeptons)[1]->pt_, abs((*selectedLeptons)[1]->eta_),"central",false, true);
            effTdata1Down = effTdata1 - scale_factor(&eff_triggermumu_dataDown_H, (*selectedLeptons)[1]->pt_, abs((*selectedLeptons)[1]->eta_),"central",false, true);
          }
          sf_Trigger = effTdata1/effTmc1;
          nominalWeights[getVecPos(sys,"trigSF")] = effTdata1/effTmc1;
          sysUpWeights[getVecPos(sys,"trigSF")] = effTdata1Up/effTmc1;
          sysDownWeights[getVecPos(sys,"trigSF")] = effTdata1Down/effTmc1;
//if ((*selectedLeptons)[0]->pt_>300 && (*selectedLeptons)[0]->pt_<500) 
//cout<<(*selectedLeptons)[0]->pt_<<":"<<TMath::Exp(0.0615- (0.0005*500))<<":my"<<(0.973 - (0.000134 * (*selectedLeptons)[0]->pt_) + (0.103 * exp((*selectedLeptons)[0]->pt_ * (-0.0118))));
      }
      if (data == "mc" && ch==2) {
        effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptons)[0]->pt_, abs((*selectedLeptons)[0]->eta_),"central",false, true);
        effTmc2 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptons)[1]->pt_, abs((*selectedLeptons)[1]->eta_),"central",false, true);
        effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptons)[0]->pt_, abs((*selectedLeptons)[0]->eta_),"central",false, true);
        effTdata2 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptons)[1]->pt_, abs((*selectedLeptons)[1]->eta_),"central",false, true);
        effTdata1Up = effTdata1 + scale_factor(&eff_triggermumu_dataUp_H, (*selectedLeptons)[0]->pt_, abs((*selectedLeptons)[0]->eta_),"central",false, true);
        effTdata1Down = effTdata1 - scale_factor(&eff_triggermumu_dataDown_H, (*selectedLeptons)[0]->pt_, abs((*selectedLeptons)[0]->eta_),"central",false, true);
        effTdata2Up = effTdata2 + scale_factor(&eff_triggermumu_dataUp_H, (*selectedLeptons)[1]->pt_, abs((*selectedLeptons)[1]->eta_),"central",false, true);
        effTdata2Down = effTdata2 - scale_factor(&eff_triggermumu_dataDown_H, (*selectedLeptons)[1]->pt_, abs((*selectedLeptons)[1]->eta_),"central",false, true);
        sf_Trigger = ((1-((1-effTdata1)*(1-effTdata2)))/(1-((1-effTmc1)*(1-effTmc2))));
        nominalWeights[getVecPos(sys,"trigSF")] = (1-((1-effTdata1)*(1-effTdata2)))/(1-((1-effTmc1)*(1-effTmc2)));
        sysUpWeights[getVecPos(sys,"trigSF")] = (1-((1-effTdata1Up)*(1-effTdata2Up)))/(1-((1-effTmc1)*(1-effTmc2)));
        sysDownWeights[getVecPos(sys,"trigSF")] = (1-((1-effTdata1Down)*(1-effTdata2Down)))/(1-((1-effTmc1)*(1-effTmc2)));
      }
//cout<<ch<<":"<<(*selectedLeptons)[0]->pt_<<":"<<exp(0.0615-(0.0005*(*selectedLeptons)[0]->pt_))<<":my"<<(0.973 - (0.000134 * (*selectedLeptons)[0]->pt_) + (0.103 * exp((*selectedLeptons)[0]->pt_ * (-0.0118))))<<endl;
//Fill histograms
      if (data == "mc"){
        if(!fname.Contains("pythia")) weightSign = signnum_typical(LHEWeight_originalXWGTUP);
        weight_lep  = weight_Lumi * weightSign * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
        weight_lepB = weight_lep * (P_bjet_data/P_bjet_mc);
        weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      }
      if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
      else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
      if(nbjet==1){
        recoL2 = (*selectedLeptons)[1]->p4_;
        for (int l=0;l<selectedJets->size();l++){
          if((*selectedJets)[l]->btag_){
            recoBjet = (*selectedJets)[l]->p4_;
            break;
          }
        }
        recoNu = Wneutrino(myMETpt, myMETphi, (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
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
      }
      resetVec(reg);
      if(leptonPass && triggerPass){
        reg[0]=0;
        wgt[0]=weight_lep;
        wcfit[0]=*eft_fit;
        if (data == "mc") truePU->Fill(Pileup_nTrueInt);
      }
      if(leptonPass && triggerPass && MetPass && DyPass){
        reg[1]=1;
        wgt[1]=weight_lep;
        wcfit[1]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjet==1){
        reg[2]=2;
        wgt[2]=weight_lepB;
        wcfit[2]=*eft_fit;
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
        BDToutput_=MVAoutput;
        ch_=ch;
        tree_out.Fill();
//        if(MVAoutput>0.9) {
//          cout<<(*selectedLeptons)[0]->pt_<<","<<(*selectedLeptons)[1]->pt_<<":"<< (*selectedLeptons)[0]->eta_<<","<<(*selectedLeptons)[1]->eta_<<" -weight="<<weight_lepB<<endl;
//          cout<< weight_Lumi<<","<<weightSign<<","<<sf_Ele_Reco<<","<<sf_Ele_ID<<","<<sf_Mu_ID<<","<<sf_Mu_ISO<<","<<sf_Trigger<<","<<weight_PU<<","<<weight_prefiring<<","<<weight_topPtPowhegData<<","<<ttKFactor<<","<<sf_JetPuId<<","<<(P_bjet_data/P_bjet_mc)<<endl;
//        }
      }
  
      if(leptonPass && triggerPass && MetPass && DyPass && nbjet>1){
        reg[3]=3;
        wgt[3]=weight_lepB;
        wcfit[3]=*eft_fit;
      }
      if(leptonPass && triggerPass && MetPass && DyInZPass && nbjet==1 && ifCR && !ifSys){
        reg[getVecPos(regions,"llB1InZ")]=getVecPos(regions,"llB1InZ");
        wgt[getVecPos(regions,"llB1InZ")]=weight_lepB;
        wcfit[getVecPos(regions,"llB1InZ")]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjet==1 && MVAoutput<0 && ifCR && !ifSys){
        reg[getVecPos(regions,"llB1-BDTm1to0")]=getVecPos(regions,"llB1-BDTm1to0");
        wgt[getVecPos(regions,"llB1-BDTm1to0")]=weight_lepB;
        wcfit[getVecPos(regions,"llB1-BDTm1to0")]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjet==1 && MVAoutput>0 && MVAoutput<0.8 && ifCR && !ifSys){
        reg[getVecPos(regions,"llB1-BDT0to0p8")]=getVecPos(regions,"llB1-BDT0to0p8");
        wgt[getVecPos(regions,"llB1-BDT0to0p8")]=weight_lepB;
        wcfit[getVecPos(regions,"llB1-BDT0to0p8")]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjet==1 && MVAoutput>0.8 && ifCR && !ifSys){
        reg[getVecPos(regions,"llB1-BDT0p8to1")]=getVecPos(regions,"llB1-BDT0p8to1");
        wgt[getVecPos(regions,"llB1-BDT0p8to1")]=weight_lepB;
        wcfit[getVecPos(regions,"llB1-BDT0p8to1")]=*eft_fit;
      }
      if(leptonPass && triggerPass && nbjet==1 && (*selectedLeptons)[0]->pt_>300  && ifCR && !ifSys){
        reg[getVecPos(regions,"llB1ptl1G300")]=getVecPos(regions,"llB1ptl1G300");
        wgt[getVecPos(regions,"llB1ptl1G300")]=weight_lepB;
        wcfit[getVecPos(regions,"llB1ptl1G300")]=*eft_fit;
      }
      if(leptonPass && triggerPass && (*selectedLeptons)[0]->pt_>300  && ifCR && !ifSys){
        reg[getVecPos(regions,"llptl1G300")]=getVecPos(regions,"llptl1G300");
        wgt[getVecPos(regions,"llptl1G300")]=weight_lepB;
        wcfit[getVecPos(regions,"llptl1G300")]=*eft_fit;
      }
      if(leptonPass && triggerPass && abs((*selectedLeptons)[0]->eta_)<1.4 && abs((*selectedLeptons)[1]->eta_)<1.4 && ifCR && !ifSys){
        reg[getVecPos(regions,"ll-BB")]=getVecPos(regions,"ll-BB");
        wgt[getVecPos(regions,"ll-BB")]=weight_lepB;
        wcfit[getVecPos(regions,"ll-BB")]=*eft_fit;
      }
      if(leptonPass && triggerPass && ((abs((*selectedLeptons)[0]->eta_)<1.4 && abs((*selectedLeptons)[1]->eta_)>1.4) || (abs((*selectedLeptons)[0]->eta_)>1.4 && abs((*selectedLeptons)[1]->eta_)<1.4)) && ifCR && !ifSys){
        reg[getVecPos(regions,"ll-BE")]=getVecPos(regions,"ll-BE");
        wgt[getVecPos(regions,"ll-BE")]=weight_lepB;
        wcfit[getVecPos(regions,"ll-BE")]=*eft_fit;
      }
      if(leptonPass && triggerPass && abs((*selectedLeptons)[0]->eta_)>1.4 && abs((*selectedLeptons)[1]->eta_)>1.4 && ifCR && !ifSys){
        reg[getVecPos(regions,"ll-EE")]=getVecPos(regions,"ll-EE");
        wgt[getVecPos(regions,"ll-EE")]=weight_lepB;
        wcfit[getVecPos(regions,"ll-EE")]=*eft_fit;
      }
      if(leptonPass && triggerPass && ifCR && !ifSys){
        if((*selectedLeptons)[0]->lep_==10 && (*selectedLeptons)[0]->charge_>0){
          reg[getVecPos(regions,"ll-leadMuP")]=getVecPos(regions,"ll-leadMuP");
          wgt[getVecPos(regions,"ll-leadMuP")]=weight_lepB;
          wcfit[getVecPos(regions,"ll-leadMuP")]=*eft_fit;
        }
      }
      if(leptonPass && triggerPass && ifCR && !ifSys){
        if((*selectedLeptons)[0]->lep_==10 && (*selectedLeptons)[0]->charge_<0){
          reg[getVecPos(regions,"ll-leadMuN")]=getVecPos(regions,"ll-leadMuN");
          wgt[getVecPos(regions,"ll-leadMuN")]=weight_lepB;
          wcfit[getVecPos(regions,"ll-leadMuN")]=*eft_fit;
        }
      }

      FillD3Hists(Hists, ch, reg, vInd(vars,"lep1Pt",0), (*selectedLeptons)[0]->pt_ ,wgt, wcfit);
      FillD3Hists(Hists, ch, reg, vInd(vars,"lep1Eta",0), (*selectedLeptons)[0]->eta_ ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"lep1Phi",0), (*selectedLeptons)[0]->phi_ ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"lep2Pt",0), (*selectedLeptons)[1]->pt_ ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"lep2Eta",0), (*selectedLeptons)[1]->eta_ ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"lep2Phi",0), (*selectedLeptons)[1]->phi_ ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"llM",0), ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"llPt",0), ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt() ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"llDr",0), deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_) ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"llDphi",0), abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_)) ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"njet",0), selectedJets->size() ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"nbjet",0), nbjet ,wgt, wcfit); 
      FillD3Hists(Hists, ch, reg, vInd(vars,"Met",0), myMETpt ,wgt, wcfit);
      FillD3Hists(Hists, ch, reg, vInd(vars,"MetPhi",0), myMETphi ,wgt, wcfit);
      FillD3Hists(Hists, ch, reg, vInd(vars,"nVtx",0), PV_npvs ,wgt, wcfit);
      FillD3Hists(Hists, ch, reg, vInd(vars,"llMZw",0), ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
      FillD3Hists(Hists, ch, reg, vInd(vars,"BDT",0), MVAoutput ,wgt, wcfit);
      FillD3Hists(Hists, ch, reg, vInd(vars,"topMass",0), recoTop.M() ,wgt, wcfit);
      FillD3Hists(Hists, ch, reg, vInd(vars,"topL1Dphi",0), abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
      FillD3Hists(Hists, ch, reg, vInd(vars,"topL1Dr",0), deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
      FillD3Hists(Hists, ch, reg, vInd(vars,"topL1DptOsumPt",0), (abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
      FillD3Hists(Hists, ch, reg, vInd(vars,"topPt",0), recoTop.Pt() ,wgt, wcfit);
      if(selectedJets->size()>0){
        FillD3Hists(Hists, ch, reg, vInd(vars,"jet1Pt",0), (*selectedJets)[0]->pt_ ,wgt, wcfit);
        FillD3Hists(Hists, ch, reg, vInd(vars,"jet1Eta",0), (*selectedJets)[0]->eta_ ,wgt, wcfit);
        FillD3Hists(Hists, ch, reg, vInd(vars,"jet1Phi",0), (*selectedJets)[0]->phi_ ,wgt, wcfit);
      }
      if (ifCR && !ifSys ) {
        if((*selectedLeptons)[0]->lep_==10) {
          FillD3Hists(Hists, ch, reg, vInd(vars,"muTkIso",0), Muon_tkRelIso[(*selectedLeptons)[0]->indice_] ,wgt, wcfit);
          FillD3Hists(Hists, ch, reg, vInd(vars,"leadMuPt",0), (*selectedLeptons)[0]->pt_ ,wgt, wcfit);
          FillD3Hists(Hists, ch, reg, vInd(vars,"leadMuEta",0), (*selectedLeptons)[0]->eta_ ,wgt, wcfit);
          FillD3Hists(Hists, ch, reg, vInd(vars,"leadMuPhi",0), (*selectedLeptons)[0]->phi_ ,wgt, wcfit);
          if(Muon_jetIdx[(*selectedLeptons)[0]->indice_]>=0) FillD3Hists(Hists, ch, reg, vInd(vars,"mu1jDrMin",0), deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,Jet_eta[Muon_jetIdx[(*selectedLeptons)[0]->indice_]],Jet_phi[Muon_jetIdx[(*selectedLeptons)[0]->indice_]]) ,wgt, wcfit);
          FillD3Hists(Hists, ch, reg, vInd(vars,"PtRelMu1jet",0),Muon_jetPtRelv2[(*selectedLeptons)[0]->indice_],wgt, wcfit);
        }
        if((*selectedLeptons)[1]->lep_==10) {
          if(Muon_jetIdx[(*selectedLeptons)[1]->indice_]>=0) FillD3Hists(Hists, ch, reg, vInd(vars,"mu2jDrMin",0), deltaR((*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_,Jet_eta[Muon_jetIdx[(*selectedLeptons)[1]->indice_]],Jet_phi[Muon_jetIdx[(*selectedLeptons)[1]->indice_]]) ,wgt, wcfit);
          FillD3Hists(Hists, ch, reg, vInd(vars,"PtRelMu2jet",0),Muon_jetPtRelv2[(*selectedLeptons)[1]->indice_],wgt, wcfit);
        }
        if((*selectedLeptons)[0]->lep_==1) {
          FillD3Hists(Hists, ch, reg, vInd(vars,"leadElePt",0), (*selectedLeptons)[0]->pt_ ,wgt, wcfit);
          FillD3Hists(Hists, ch, reg, vInd(vars,"leadEleEta",0), (*selectedLeptons)[0]->eta_ ,wgt, wcfit);
          FillD3Hists(Hists, ch, reg, vInd(vars,"leadElePhi",0), (*selectedLeptons)[0]->phi_ ,wgt, wcfit);
          if(Electron_jetIdx[(*selectedLeptons)[0]->indice_]>=0) FillD3Hists(Hists, ch, reg, vInd(vars,"ele1jDrMin",0), deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,Jet_eta[Electron_jetIdx[(*selectedLeptons)[0]->indice_]],Jet_phi[Electron_jetIdx[(*selectedLeptons)[0]->indice_]]) ,wgt, wcfit);
          FillD3Hists(Hists, ch, reg, vInd(vars,"PtRelEle1jet",0),Electron_jetPtRelv2[(*selectedLeptons)[0]->indice_],wgt, wcfit);
        }
        if((*selectedLeptons)[1]->lep_==1) {
          if(Electron_jetIdx[(*selectedLeptons)[1]->indice_]>=0) FillD3Hists(Hists, ch, reg, vInd(vars,"ele2jDrMin",0), deltaR((*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_,Jet_eta[Electron_jetIdx[(*selectedLeptons)[1]->indice_]],Jet_phi[Electron_jetIdx[(*selectedLeptons)[1]->indice_]]) ,wgt, wcfit);
          FillD3Hists(Hists, ch, reg, vInd(vars,"PtRelEle2jet",0),Electron_jetPtRelv2[(*selectedLeptons)[1]->indice_],wgt, wcfit);
        }
        if(selectedJets->size()>0){
          if((*selectedLeptons)[0]->lep_==10) FillD3Hists(Hists, ch, reg, vInd(vars,"mu1bDr",0), deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedJets)[lbjet]->eta_,(*selectedJets)[lbjet]->phi_) ,wgt, wcfit);
          if((*selectedLeptons)[1]->lep_==10) FillD3Hists(Hists, ch, reg, vInd(vars,"mu2bDr",0), deltaR((*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_,(*selectedJets)[lbjet]->eta_,(*selectedJets)[lbjet]->phi_) ,wgt, wcfit);
          if((*selectedLeptons)[0]->lep_==1) FillD3Hists(Hists, ch, reg, vInd(vars,"ele1bDr",0), deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedJets)[lbjet]->eta_,(*selectedJets)[lbjet]->phi_) ,wgt, wcfit);
          if((*selectedLeptons)[1]->lep_==1) FillD3Hists(Hists, ch, reg, vInd(vars,"ele2bDr",0), deltaR((*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_,(*selectedJets)[lbjet]->eta_,(*selectedJets)[lbjet]->phi_) ,wgt, wcfit);
        }
      }  
      delete eft_fit;
    }
    if (data != "mc" || fname.Contains("sys") || !ifSys){
      for (int l=0;l<selectedLeptons->size();l++) delete (*selectedLeptons)[l];
      for (int l=0;l<selectedLeptonsMuScaleUp->size();l++) delete (*selectedLeptonsMuScaleUp)[l];
      for (int l=0;l<selectedLeptonsMuScaleDown->size();l++) delete (*selectedLeptonsMuScaleDown)[l];
      for (int l=0;l<selectedLeptonsEleScaleUp->size();l++) delete (*selectedLeptonsEleScaleUp)[l];
      for (int l=0;l<selectedLeptonsEleScaleDown->size();l++) delete (*selectedLeptonsEleScaleDown)[l];
      for (int l=0;l<selectedLeptonsMuResUp->size();l++) delete (*selectedLeptonsMuResUp)[l];
      for (int l=0;l<selectedLeptonsMuResDown->size();l++) delete (*selectedLeptonsMuResDown)[l];
      
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
  
      for (int l=0;l<selectedJets->size();l++) delete (*selectedJets)[l];
      for (int l=0;l<selectedJetsJerUp->size();l++) delete (*selectedJetsJerUp)[l];
      for (int l=0;l<selectedJetsJerDown->size();l++) delete (*selectedJetsJerDown)[l];
      for (int l=0;l<selectedJetsJesUp->size();l++) delete (*selectedJetsJesUp)[l];
      for (int l=0;l<selectedJetsJesDown->size();l++) delete (*selectedJetsJesDown)[l];
  
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
        (*JECsysUp)[l].clear();
        (*JECsysUp)[l].shrink_to_fit();
      }
      for (int l=0;l<JECsysDown->size();l++){
        for (int n=0;n<(*JECsysDown)[l].size();n++){
          delete (*JECsysDown)[l][n];
        }
        (*JECsysDown)[l].clear();
        (*JECsysDown)[l].shrink_to_fit();
      }
  
      JECsysUp->clear();
      JECsysDown->clear();
      JECsysNbtagUp->clear();
      JECsysNbtagDown->clear();
      JECsysMETUp->clear();
      JECsysMETDown->clear();
      JECsysMETPhiUp->clear();
      JECsysMETPhiDown->clear();
      JECsysUp->shrink_to_fit();
      JECsysDown->shrink_to_fit();
      JECsysNbtagUp->shrink_to_fit();
      JECsysNbtagDown->shrink_to_fit();
      JECsysMETUp->shrink_to_fit();
      JECsysMETDown->shrink_to_fit();
      JECsysMETPhiUp->shrink_to_fit();
      JECsysMETPhiDown->shrink_to_fit();
      delete JECsysUp;
      delete JECsysDown;
      delete JECsysNbtagUp;
      delete JECsysNbtagDown;
      delete JECsysMETUp;
      delete JECsysMETDown;
      delete JECsysMETPhiUp;
      delete JECsysMETPhiDown;
      continue;
    }
//include systematic histograms
    if(ch<10 && selectedLeptons->size() ==2){
      if(useRegions(reg)){
        for (int n=0;n<sys.size();++n){
          if (std::find(sysNotWeight.begin(), sysNotWeight.end(), sys[n]) != sysNotWeight.end()) continue;
          if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* (sysUpWeights[n]/nominalWeights[n]));
          else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
          for(int l=0;l<reg.size();++l){
            if(l<2)  wgt[l]=(weight_lep* (sysUpWeights[n]/nominalWeights[n]));
            else wgt[l]=(weight_lepB* (sysUpWeights[n]/nominalWeights[n]));
            wcfit[l]=*eft_fit;
          }
  
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Pt",4), n, (*selectedLeptons)[0]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Eta",4), n, (*selectedLeptons)[0]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Phi",4), n, (*selectedLeptons)[0]->phi_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Pt",4), n, (*selectedLeptons)[1]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Eta",4), n, (*selectedLeptons)[1]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Phi",4), n, (*selectedLeptons)[1]->phi_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llM",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llPt",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt() ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llDr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_) ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llDphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_)) ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"njet",4), n, selectedJets->size() ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"nbjet",4), n, nbjet ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"Met",4), n, myMETpt ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"MetPhi",4), n, myMETphi ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"nVtx",4), n, PV_npvs ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llMZw",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"BDT",4), n, MVAoutput ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topMass",4), n, recoTop.M() ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1Dphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1Dr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1DptOsumPt",4), n, (abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topPt",4), n, recoTop.Pt() ,wgt, wcfit);
          if(selectedJets->size()>0){
            FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Pt",4), n, (*selectedJets)[0]->pt_ ,wgt, wcfit);
            FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Eta",4), n, (*selectedJets)[0]->eta_ ,wgt, wcfit);
            FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Phi",4), n, (*selectedJets)[0]->phi_ ,wgt, wcfit);
          }
          delete eft_fit;
          if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* (sysDownWeights[n]/nominalWeights[n]));
          else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
          for(int l=0;l<reg.size();++l){
            if(l<2) wgt[l]=(weight_lep* (sysDownWeights[n]/nominalWeights[n]));
            else wgt[l]=(weight_lepB* (sysDownWeights[n]/nominalWeights[n]));
            wcfit[l]=*eft_fit;
          }
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Pt",4), n, (*selectedLeptons)[0]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Eta",4), n, (*selectedLeptons)[0]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Phi",4), n, (*selectedLeptons)[0]->phi_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Pt",4), n, (*selectedLeptons)[1]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Eta",4), n, (*selectedLeptons)[1]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Phi",4), n, (*selectedLeptons)[1]->phi_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llM",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llPt",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt() ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llDr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_) ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llDphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_)) ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"njet",4), n, selectedJets->size() ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"nbjet",4), n, nbjet ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"Met",4), n, myMETpt ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"MetPhi",4), n, myMETphi ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"nVtx",4), n, PV_npvs ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llMZw",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"BDT",4), n, MVAoutput ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topMass",4), n, recoTop.M() ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1Dphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1Dr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1DptOsumPt",4), n, (abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topPt",4), n, recoTop.Pt() ,wgt, wcfit);
          if(selectedJets->size()>0){
            FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Pt",4), n, (*selectedJets)[0]->pt_ ,wgt, wcfit);
            FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Eta",4), n, (*selectedJets)[0]->eta_ ,wgt, wcfit);
            FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Phi",4), n, (*selectedJets)[0]->phi_ ,wgt, wcfit);
          }
          delete eft_fit;
        }
  //PDF, QS, PS
        if ((fname.Contains("TTTo2L2Nu")&& !fname.Contains("sys")) || fname.Contains("BNV")){
           resetVec(reg);
          if(leptonPass && triggerPass && DyPass && MetPass && nbjet==1) reg[0]=0;
          if(leptonPass && triggerPass && MetPass && DyPass && nbjet>1) reg[1]=1;
          for (int n=0;n<nScale;++n){
            if(isnan(LHEScaleWeight[n]) || isinf(LHEScaleWeight[n])) continue;
            if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* LHEScaleWeight[n]);
            else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
            for(int l=0;l<reg.size();++l){
              wgt[l]=(weight_lepB * LHEScaleWeight[n]);
              wcfit[l]=*eft_fit;
            }
            FillD4Hists(HistsSysReweightsQscale, ch, reg, 0, n, MVAoutput ,wgt, wcfit);
            delete eft_fit;
          }
          for (int n=0;n<nPdf;++n){
            if(isnan(LHEPdfWeight[n]) || isinf(LHEPdfWeight[n])) continue;
            if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* LHEPdfWeight[n]);
            else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
            for(int l=0;l<reg.size();++l){
              wgt[l]=(weight_lepB * LHEPdfWeight[n]);
              wcfit[l]=*eft_fit;
            }
            FillD4Hists(HistsSysReweightsPDF, ch, reg, 0, n, MVAoutput ,wgt, wcfit);
            delete eft_fit;
          }
          for (int n=0;n<nPS;++n){
            if(isnan(PSWeight[n]) || isinf(PSWeight[n])) continue;
            if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT* PSWeight[n]);
            else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
            for(int l=0;l<reg.size();++l){
              wgt[l]=(weight_lepB * PSWeight[n]);
              wcfit[l]=*eft_fit;
            }
            FillD4Hists(HistsSysReweightsPS, ch, reg, 0, n, MVAoutput ,wgt, wcfit);
            delete eft_fit;
          }
        }
      }
//Uncluster MET uncertainty      
      if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
      else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
      //if (leptonPass && triggerPass && myMETpt > MetCut && ch==1 && DyPass && nbjet==1) nMet = nMet +weight_lepB;
      //if (leptonPass && triggerPass && metUnclusMETUp > MetCut && ch==1  && DyPass && nbjet==1) nMetUp=nMetUp+weight_lepB;
      //if (leptonPass && triggerPass && metUnclusMETDown > MetCut && ch==1  && DyPass && nbjet==1) nMetDown=nMetDown+weight_lepB;
      MetPass = (metUnclusMETUp > MetCut) ? true : false;
      resetVec(reg);
      if(leptonPass && triggerPass){
        reg[0]=0;
        wgt[0]=weight_lep;
        wcfit[0]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass){
        reg[1]=1;
        wgt[1]=weight_lep;
        wcfit[1]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjet==1){
        reg[2]=2;
        wgt[2]=weight_lepB;
        wcfit[2]=*eft_fit;
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
      }
      if(leptonPass && triggerPass && MetPass && DyPass && nbjet>1){
        reg[3]=3;
        wgt[3]=weight_lepB;
        wcfit[3]=*eft_fit;
      }
      for (int n=0;n<sys.size();++n){
        if (sys[n]!="unclusMET") continue;
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Pt",4), n, (*selectedLeptons)[0]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Eta",4), n, (*selectedLeptons)[0]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Phi",4), n, (*selectedLeptons)[0]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Pt",4), n, (*selectedLeptons)[1]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Eta",4), n, (*selectedLeptons)[1]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Phi",4), n, (*selectedLeptons)[1]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llM",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llPt",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llDr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llDphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_)) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"njet",4), n, selectedJets->size() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"nbjet",4), n, nbjet ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"Met",4), n, metUnclusMETUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"MetPhi",4), n, metUnclusMETPhiUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"nVtx",4), n, PV_npvs ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llMZw",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"BDT",4), n, MVAoutput ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topMass",4), n, recoTop.M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1Dphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1Dr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1DptOsumPt",4), n, (abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topPt",4), n, recoTop.Pt() ,wgt, wcfit);
        if(selectedJets->size()>0){
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Pt",4), n, (*selectedJets)[0]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Eta",4), n, (*selectedJets)[0]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Phi",4), n, (*selectedJets)[0]->phi_ ,wgt, wcfit);
        }
      }

      MetPass = (metUnclusMETDown > MetCut) ? true : false;
      resetVec(reg);
      if(leptonPass && triggerPass){
        reg[0]=0;
        wgt[0]=weight_lep;
        wcfit[0]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass){
        reg[1]=1;
        wgt[1]=weight_lep;
        wcfit[1]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjet==1){
        reg[2]=2;
        wgt[2]=weight_lepB;
        wcfit[2]=*eft_fit;
        for (int l=0;l<selectedJets->size();l++){
          if((*selectedJets)[l]->btag_){
            recoBjet = (*selectedJets)[l]->p4_;
            break;
          }
        }
        recoNu = Wneutrino(metUnclusMETDown, metUnclusMETPhiDown, (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
        recoW = recoNu + recoL2;
        recoTop = recoW + recoBjet;
        MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
        MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
        MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
        MVA_topPt=recoTop.Pt();
        MVAoutput = readerMVA->EvaluateMVA( "BDTG");
      }
      if(leptonPass && triggerPass && MetPass && DyPass && nbjet>1){
        reg[3]=3;
        wgt[3]=weight_lepB;
        wcfit[3]=*eft_fit;
      }
      for (int n=0;n<sys.size();++n){
        if (sys[n]!="unclusMET") continue;
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Pt",4), n, (*selectedLeptons)[0]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Eta",4), n, (*selectedLeptons)[0]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Phi",4), n, (*selectedLeptons)[0]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Pt",4), n, (*selectedLeptons)[1]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Eta",4), n, (*selectedLeptons)[1]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Phi",4), n, (*selectedLeptons)[1]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llM",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llPt",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llDr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llDphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_)) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"njet",4), n, selectedJets->size() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"nbjet",4), n, nbjet ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"Met",4), n, metUnclusMETDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"MetPhi",4), n, metUnclusMETPhiDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"nVtx",4), n, PV_npvs ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llMZw",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"BDT",4), n, MVAoutput ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topMass",4), n, recoTop.M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1Dphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1Dr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1DptOsumPt",4), n, (abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topPt",4), n, recoTop.Pt() ,wgt, wcfit);
        if(selectedJets->size()>0){
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Pt",4), n, (*selectedJets)[0]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Eta",4), n, (*selectedJets)[0]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Phi",4), n, (*selectedJets)[0]->phi_ ,wgt, wcfit);
        }
      }
//JES uncertainty
      MetPass = (MetJetsPtJesUp > MetCut) ? true : false;
      resetVec(reg);
      if(leptonPass && triggerPass){
        reg[0]=0;
        wgt[0]=weight_lep;
        wcfit[0]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass){
        reg[1]=1;
        wgt[1]=weight_lep;
        wcfit[1]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjetJesUp==1){
        reg[2]=2;
        wgt[2]=weight_lepB;
        wcfit[2]=*eft_fit;
        for (int l=0;l<selectedJetsJesUp->size();l++){
          if((*selectedJetsJesUp)[l]->btag_){
            recoBjet = (*selectedJetsJesUp)[l]->p4_;
            break;
          }
        }
        recoNu = Wneutrino(MetJetsPtJesUp, MetJetsPhiJesUp, (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
        recoW = recoNu + recoL2;
        recoTop = recoW + recoBjet;
        MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
        MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
        MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
        MVA_topPt=recoTop.Pt();
        MVAoutput = readerMVA->EvaluateMVA( "BDTG");
      }
      if(leptonPass && triggerPass && MetPass && DyPass && nbjetJesUp>1){
        reg[3]=3;
        wgt[3]=weight_lepB;
        wcfit[3]=*eft_fit;
      }
      for (int n=0;n<sys.size();++n){
        if (sys[n]!="jes") continue;
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Pt",4), n, (*selectedLeptons)[0]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Eta",4), n, (*selectedLeptons)[0]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Phi",4), n, (*selectedLeptons)[0]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Pt",4), n, (*selectedLeptons)[1]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Eta",4), n, (*selectedLeptons)[1]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Phi",4), n, (*selectedLeptons)[1]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llM",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llPt",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llDr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llDphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_)) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"njet",4), n, selectedJetsJesUp->size() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"nbjet",4), n, nbjetJesUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"Met",4), n, MetJetsPtJesUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"MetPhi",4), n, MetJetsPhiJesUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"nVtx",4), n, PV_npvs ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llMZw",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"BDT",4), n, MVAoutput ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topMass",4), n, recoTop.M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1Dphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1Dr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1DptOsumPt",4), n, (abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topPt",4), n, recoTop.Pt() ,wgt, wcfit);
        if(selectedJetsJesUp->size()>0){
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Pt",4), n, (*selectedJetsJesUp)[0]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Eta",4), n, (*selectedJetsJesUp)[0]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Phi",4), n, (*selectedJetsJesUp)[0]->phi_ ,wgt, wcfit);
        }
      }
      MetPass = (MetJetsPtJesDown > MetCut) ? true : false;
      resetVec(reg);
      if(leptonPass && triggerPass){
        reg[0]=0;
        wgt[0]=weight_lep;
        wcfit[0]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass){
        reg[1]=1;
        wgt[1]=weight_lep;
        wcfit[1]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjetJesDown==1){
        reg[2]=2;
        wgt[2]=weight_lepB;
        wcfit[2]=*eft_fit;
        for (int l=0;l<selectedJetsJesDown->size();l++){
          if((*selectedJetsJesDown)[l]->btag_){
            recoBjet = (*selectedJetsJesDown)[l]->p4_;
            break;
          }
        }
        recoNu = Wneutrino(MetJetsPtJesDown, MetJetsPhiJesDown, (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
        recoW = recoNu + recoL2;
        recoTop = recoW + recoBjet;
        MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
        MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
        MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
        MVA_topPt=recoTop.Pt();
        MVAoutput = readerMVA->EvaluateMVA( "BDTG");
      }
      if(leptonPass && triggerPass && MetPass && DyPass && nbjetJesDown>1){
        reg[3]=3;
        wgt[3]=weight_lepB;
        wcfit[3]=*eft_fit;
      }
      for (int n=0;n<sys.size();++n){
        if (sys[n]!="jes") continue;
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Pt",4), n, (*selectedLeptons)[0]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Eta",4), n, (*selectedLeptons)[0]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Phi",4), n, (*selectedLeptons)[0]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Pt",4), n, (*selectedLeptons)[1]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Eta",4), n, (*selectedLeptons)[1]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Phi",4), n, (*selectedLeptons)[1]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llM",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llPt",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llDr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llDphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_)) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"njet",4), n, selectedJetsJesDown->size() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"nbjet",4), n, nbjetJesDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"Met",4), n, MetJetsPtJesDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"MetPhi",4), n, MetJetsPhiJesDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"nVtx",4), n, PV_npvs ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llMZw",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"BDT",4), n, MVAoutput ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topMass",4), n, recoTop.M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1Dphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1Dr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1DptOsumPt",4), n, (abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topPt",4), n, recoTop.Pt() ,wgt, wcfit);
        if(selectedJetsJesDown->size()>0){
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Pt",4), n, (*selectedJetsJesDown)[0]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Eta",4), n, (*selectedJetsJesDown)[0]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Phi",4), n, (*selectedJetsJesDown)[0]->phi_ ,wgt, wcfit);
        }
      }
//JER uncertainty
      MetPass = (MetJetsPtJerUp > MetCut) ? true : false;
      resetVec(reg);
      if(leptonPass && triggerPass){
        reg[0]=0;
        wgt[0]=weight_lep;
        wcfit[0]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass){
        reg[1]=1;
        wgt[1]=weight_lep;
        wcfit[1]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjetJerUp==1){
        reg[2]=2;
        wgt[2]=weight_lepB;
        wcfit[2]=*eft_fit;
        for (int l=0;l<selectedJetsJerUp->size();l++){
          if((*selectedJetsJerUp)[l]->btag_){
            recoBjet = (*selectedJetsJerUp)[l]->p4_;
            break;
          }
        }
        recoNu = Wneutrino(MetJetsPtJerUp, MetJetsPhiJerUp, (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
        recoW = recoNu + recoL2;
        recoTop = recoW + recoBjet;
        MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
        MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
        MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
        MVA_topPt=recoTop.Pt();
        MVAoutput = readerMVA->EvaluateMVA( "BDTG");
      }
      if(leptonPass && triggerPass && MetPass && DyPass && nbjetJerUp>1){
        reg[3]=3;
        wgt[3]=weight_lepB;
        wcfit[3]=*eft_fit;
      }
      for (int n=0;n<sys.size();++n){
        if (sys[n]!="jer") continue;
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Pt",4), n, (*selectedLeptons)[0]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Eta",4), n, (*selectedLeptons)[0]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Phi",4), n, (*selectedLeptons)[0]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Pt",4), n, (*selectedLeptons)[1]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Eta",4), n, (*selectedLeptons)[1]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Phi",4), n, (*selectedLeptons)[1]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llM",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llPt",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llDr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llDphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_)) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"njet",4), n, selectedJetsJerUp->size() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"nbjet",4), n, nbjetJerUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"Met",4), n, MetJetsPtJerUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"MetPhi",4), n, MetJetsPhiJerUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"nVtx",4), n, PV_npvs ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llMZw",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"BDT",4), n, MVAoutput ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topMass",4), n, recoTop.M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1Dphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1Dr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1DptOsumPt",4), n, (abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topPt",4), n, recoTop.Pt() ,wgt, wcfit);
        if(selectedJetsJerUp->size()>0){
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Pt",4), n, (*selectedJetsJerUp)[0]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Eta",4), n, (*selectedJetsJerUp)[0]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Phi",4), n, (*selectedJetsJerUp)[0]->phi_ ,wgt, wcfit);
        }
      }

      MetPass = (MetJetsPtJerDown > MetCut) ? true : false;
      resetVec(reg);
      if(leptonPass && triggerPass){
        reg[0]=0;
        wgt[0]=weight_lep;
        wcfit[0]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass){
        reg[1]=1;
        wgt[1]=weight_lep;
        wcfit[1]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjetJerDown==1){
        reg[2]=2;
        wgt[2]=weight_lepB;
        wcfit[2]=*eft_fit;
        for (int l=0;l<selectedJetsJerDown->size();l++){
          if((*selectedJetsJerDown)[l]->btag_){
            recoBjet = (*selectedJetsJerDown)[l]->p4_;
            break;
          }
        }
        recoNu = Wneutrino(MetJetsPtJerDown, MetJetsPhiJerDown, (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
        recoW = recoNu + recoL2;
        recoTop = recoW + recoBjet;
        MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
        MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
        MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
        MVA_topPt=recoTop.Pt();
        MVAoutput = readerMVA->EvaluateMVA( "BDTG");
      }
      if(leptonPass && triggerPass && MetPass && DyPass && nbjetJerDown>1){
        reg[3]=3;
        wgt[3]=weight_lepB;
        wcfit[3]=*eft_fit;
      }
      for (int n=0;n<sys.size();++n){
        if (sys[n]!="jer") continue;
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Pt",4), n, (*selectedLeptons)[0]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Eta",4), n, (*selectedLeptons)[0]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Phi",4), n, (*selectedLeptons)[0]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Pt",4), n, (*selectedLeptons)[1]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Eta",4), n, (*selectedLeptons)[1]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Phi",4), n, (*selectedLeptons)[1]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llM",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llPt",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llDr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llDphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_)) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"njet",4), n, selectedJetsJerDown->size() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"nbjet",4), n, nbjetJerDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"Met",4), n, MetJetsPtJerDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"MetPhi",4), n, MetJetsPhiJerDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"nVtx",4), n, PV_npvs ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llMZw",4), n, ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"BDT",4), n, MVAoutput ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topMass",4), n, recoTop.M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1Dphi",4), n, abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1Dr",4), n, deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1DptOsumPt",4), n, (abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topPt",4), n, recoTop.Pt() ,wgt, wcfit);
        if(selectedJetsJerDown->size()>0){
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Pt",4), n, (*selectedJetsJerDown)[0]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Eta",4), n, (*selectedJetsJerDown)[0]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Phi",4), n, (*selectedJetsJerDown)[0]->phi_ ,wgt, wcfit);
        }
      }
//JES uncertainty sub sources
      if(leptonPass && triggerPass && DyPass){
        for (int n=0;n<sysJecNames.size();++n){
          if((*JECsysNbtagUp)[n]==0) continue;
          MetPass = ((*JECsysMETUp)[n] > MetCut) ? true : false;
          resetVec(reg);
          if(leptonPass && triggerPass && DyPass && MetPass && (*JECsysNbtagUp)[n]==1){
            reg[0]=0;
            wgt[0]=weight_lepB;
            wcfit[0]=*eft_fit;
            for (int l=0;l<(*JECsysUp)[n].size();l++){
              if((*JECsysUp)[n][l]->btag_){
                recoBjet = (*JECsysUp)[n][l]->p4_;
                break;
              }
            }
            recoNu = Wneutrino((*JECsysMETUp)[n], (*JECsysMETPhiUp)[n], (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
            recoW = recoNu + recoL2;
            recoTop = recoW + recoBjet;
            MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
            MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
            MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
            MVA_topPt=recoTop.Pt();
            MVAoutput = readerMVA->EvaluateMVA( "BDTG");
          }
          if(leptonPass && triggerPass && MetPass && DyPass && (*JECsysNbtagUp)[n]>1){
            reg[1]=1;
            wgt[1]=weight_lepB;
            wcfit[1]=*eft_fit;
          }
          FillD4Hists(HistsJecUp, ch, reg, 0, n, MVAoutput ,wgt, wcfit);
        }
        for (int n=0;n<sysJecNames.size();++n){
          if((*JECsysNbtagDown)[n]==0) continue;
          MetPass = ((*JECsysMETDown)[n] > MetCut) ? true : false;
          resetVec(reg);
          if(leptonPass && triggerPass && DyPass && MetPass && (*JECsysNbtagDown)[n]==1){
            reg[0]=0;
            wgt[0]=weight_lepB;
            wcfit[0]=*eft_fit;
            for (int l=0;l<(*JECsysDown)[n].size();l++){
              if((*JECsysDown)[n][l]->btag_){
                recoBjet = (*JECsysDown)[n][l]->p4_;
                break;
              }
            }
            recoNu = Wneutrino((*JECsysMETDown)[n], (*JECsysMETPhiDown)[n], (*selectedLeptons)[1]->pt_, (*selectedLeptons)[1]->eta_, (*selectedLeptons)[1]->phi_);
            recoW = recoNu + recoL2;
            recoTop = recoW + recoBjet;
            MVA_topL1Dphi=abs(deltaPhi((*selectedLeptons)[0]->phi_,recoTop.Phi()));
            MVA_topL1Dr=deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,recoTop.Eta(),recoTop.Phi());
            MVA_topL1DptOsumPt=(abs((*selectedLeptons)[0]->pt_ - recoTop.Pt()))/((*selectedLeptons)[0]->pt_ + recoTop.Pt());
            MVA_topPt=recoTop.Pt();
            MVAoutput = readerMVA->EvaluateMVA( "BDTG");
          }
          if(leptonPass && triggerPass && MetPass && DyPass && (*JECsysNbtagDown)[n]>1){
            reg[1]=1;
            wgt[1]=weight_lepB;
            wcfit[1]=*eft_fit;
          }
          FillD4Hists(HistsJecDown, ch, reg, 0, n, MVAoutput ,wgt, wcfit);
        }
      }
      delete eft_fit;
    }

// Lepton scale resolution uncertainties
//Muon scale
    if(selectedLeptonsMuScaleUp->size() ==2){
      leptonPass=false;
      triggerPass=false;
      DyPass = false;
      MetPass = false;
      ch=100;
      if (((*selectedLeptonsMuScaleUp)[0]->pt_ > 25) && ((*selectedLeptonsMuScaleUp)[0]->charge_ * (*selectedLeptonsMuScaleUp)[1]->charge_ == -1) && ((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).M()>20 && ((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).M()<1200) leptonPass=true;
      if ((*selectedLeptonsMuScaleUp)[0]->lep_ + (*selectedLeptonsMuScaleUp)[1]->lep_ == 2) ch = 0;
      if ((*selectedLeptonsMuScaleUp)[0]->lep_ + (*selectedLeptonsMuScaleUp)[1]->lep_ == 11) ch = 1;
      if ((*selectedLeptonsMuScaleUp)[0]->lep_ + (*selectedLeptonsMuScaleUp)[1]->lep_ == 20) ch = 2;
      if (ch ==0 && triggerPassEE) triggerPass=true;
      if (ch ==1 && triggerPassEMu) triggerPass=true;
      if (ch ==2 && triggerPassMuMu) triggerPass=true;
      if (((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).M()>106 ) DyPass = true;
      if (myMETpt > MetCut) MetPass = true;
      if (ch==0)  sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptonsMuScaleUp)[0]->pt_, (*selectedLeptonsMuScaleUp)[1]->pt_,"central",false, true);
      if (ch==1) {
        if((*selectedLeptonsMuScaleUp)[0]->lep_==10) {
          effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuScaleUp)[0]->pt_, abs((*selectedLeptonsMuScaleUp)[0]->eta_),"central",false, true);
          effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuScaleUp)[0]->pt_, abs((*selectedLeptonsMuScaleUp)[0]->eta_),"central",false, true);
        }
        else{
          effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuScaleUp)[1]->pt_, abs((*selectedLeptonsMuScaleUp)[1]->eta_),"central",false, true);
          effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuScaleUp)[1]->pt_, abs((*selectedLeptonsMuScaleUp)[1]->eta_),"central",false, true);
        }
        sf_Trigger = effTdata1/effTmc1;
      }
      if (data == "mc" && ch==2) {
        effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuScaleUp)[0]->pt_, abs((*selectedLeptonsMuScaleUp)[0]->eta_),"central",false, true);
        effTmc2 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuScaleUp)[1]->pt_, abs((*selectedLeptonsMuScaleUp)[1]->eta_),"central",false, true);
        effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuScaleUp)[0]->pt_, abs((*selectedLeptonsMuScaleUp)[0]->eta_),"central",false, true);
        effTdata2 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuScaleUp)[1]->pt_, abs((*selectedLeptonsMuScaleUp)[1]->eta_),"central",false, true);
        sf_Trigger = ((1-((1-effTdata1)*(1-effTdata2)))/(1-((1-effTmc1)*(1-effTmc2))));
      }
      if(!fname.Contains("pythia")) weightSign = signnum_typical(LHEWeight_originalXWGTUP);
      weight_lep = weight_Lumi * weightSign * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      weight_lepB = weight_Lumi * weightSign * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
      else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
      resetVec(reg);
      if(leptonPass && triggerPass){
        reg[0]=0;
        wgt[0]=weight_lep;
        wcfit[0]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass){
        reg[1]=1;
        wgt[1]=weight_lep;
        wcfit[1]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjet==1){
        reg[2]=2;
        wgt[2]=weight_lepB;
        wcfit[2]=*eft_fit;
        for (int l=0;l<selectedJets->size();l++){
          if((*selectedJets)[l]->btag_){
            recoBjet = (*selectedJets)[l]->p4_;
            break;
          }
        }
        recoL2 = (*selectedLeptonsMuScaleUp)[1]->p4_;
        recoNu = Wneutrino(myMETpt, myMETphi, (*selectedLeptonsMuScaleUp)[1]->pt_, (*selectedLeptonsMuScaleUp)[1]->eta_, (*selectedLeptonsMuScaleUp)[1]->phi_);
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
      }
      if(leptonPass && triggerPass && MetPass && DyPass && nbjet>1){
        reg[3]=3;
        wgt[3]=weight_lepB;
        wcfit[3]=*eft_fit;
      }
      for (int n=0;n<sys.size();++n){
        if (sys[n]!="muonScale") continue;
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Pt",4), n, (*selectedLeptonsMuScaleUp)[0]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Eta",4), n, (*selectedLeptonsMuScaleUp)[0]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Phi",4), n, (*selectedLeptonsMuScaleUp)[0]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Pt",4), n, (*selectedLeptonsMuScaleUp)[1]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Eta",4), n, (*selectedLeptonsMuScaleUp)[1]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Phi",4), n, (*selectedLeptonsMuScaleUp)[1]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llM",4), n, ((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llPt",4), n, ((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).Pt() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llDr",4), n, deltaR((*selectedLeptonsMuScaleUp)[0]->eta_,(*selectedLeptonsMuScaleUp)[0]->phi_,(*selectedLeptonsMuScaleUp)[1]->eta_,(*selectedLeptonsMuScaleUp)[1]->phi_) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llDphi",4), n, abs(deltaPhi((*selectedLeptonsMuScaleUp)[0]->phi_,(*selectedLeptonsMuScaleUp)[1]->phi_)) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"njet",4), n, selectedJets->size() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"nbjet",4), n, nbjetJerUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"Met",4), n, MetJetsPtJerUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"MetPhi",4), n, MetJetsPhiJerUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"nVtx",4), n, PV_npvs ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llMZw",4), n, ((*selectedLeptonsMuScaleUp)[0]->p4_ + (*selectedLeptonsMuScaleUp)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"BDT",4), n, MVAoutput ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topMass",4), n, recoTop.M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1Dphi",4), n, abs(deltaPhi((*selectedLeptonsMuScaleUp)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1Dr",4), n, deltaR((*selectedLeptonsMuScaleUp)[0]->eta_,(*selectedLeptonsMuScaleUp)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1DptOsumPt",4), n, (abs((*selectedLeptonsMuScaleUp)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuScaleUp)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topPt",4), n, recoTop.Pt() ,wgt, wcfit);
        if(selectedJets->size()>0){
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Pt",4), n, (*selectedJets)[0]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Eta",4), n, (*selectedJets)[0]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Phi",4), n, (*selectedJets)[0]->phi_ ,wgt, wcfit);
        }
      }
    }

    if(selectedLeptonsMuScaleDown->size() ==2){
      leptonPass=false;
      triggerPass=false;
      DyPass = false;
      MetPass = false;
      ch=100;
      if (((*selectedLeptonsMuScaleDown)[0]->pt_ > 25) && ((*selectedLeptonsMuScaleDown)[0]->charge_ * (*selectedLeptonsMuScaleDown)[1]->charge_ == -1) && ((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).M()>20 && ((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).M()<1200) leptonPass=true;
      if ((*selectedLeptonsMuScaleDown)[0]->lep_ + (*selectedLeptonsMuScaleDown)[1]->lep_ == 2) ch = 0;
      if ((*selectedLeptonsMuScaleDown)[0]->lep_ + (*selectedLeptonsMuScaleDown)[1]->lep_ == 11) ch = 1;
      if ((*selectedLeptonsMuScaleDown)[0]->lep_ + (*selectedLeptonsMuScaleDown)[1]->lep_ == 20) ch = 2;
      if (ch ==0 && triggerPassEE) triggerPass=true;
      if (ch ==1 && triggerPassEMu) triggerPass=true;
      if (ch ==2 && triggerPassMuMu) triggerPass=true;
      if (((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).M()>106 ) DyPass = true;
      if (myMETpt > MetCut) MetPass = true;
      if (ch==0)  sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptonsMuScaleDown)[0]->pt_, (*selectedLeptonsMuScaleDown)[1]->pt_,"central",false, true);
      if (ch==1) {
        if((*selectedLeptonsMuScaleDown)[0]->lep_==10) {
          effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuScaleDown)[0]->pt_, abs((*selectedLeptonsMuScaleDown)[0]->eta_),"central",false, true);
          effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuScaleDown)[0]->pt_, abs((*selectedLeptonsMuScaleDown)[0]->eta_),"central",false, true);
        }
        else{
          effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuScaleDown)[1]->pt_, abs((*selectedLeptonsMuScaleDown)[1]->eta_),"central",false, true);
          effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuScaleDown)[1]->pt_, abs((*selectedLeptonsMuScaleDown)[1]->eta_),"central",false, true);
        }
        sf_Trigger = effTdata1/effTmc1;
      }
      if (data == "mc" && ch==2) {
        effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuScaleDown)[0]->pt_, abs((*selectedLeptonsMuScaleDown)[0]->eta_),"central",false, true);
        effTmc2 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuScaleDown)[1]->pt_, abs((*selectedLeptonsMuScaleDown)[1]->eta_),"central",false, true);
        effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuScaleDown)[0]->pt_, abs((*selectedLeptonsMuScaleDown)[0]->eta_),"central",false, true);
        effTdata2 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuScaleDown)[1]->pt_, abs((*selectedLeptonsMuScaleDown)[1]->eta_),"central",false, true);
        sf_Trigger = ((1-((1-effTdata1)*(1-effTdata2)))/(1-((1-effTmc1)*(1-effTmc2))));
      }
      if(!fname.Contains("pythia")) weightSign = signnum_typical(LHEWeight_originalXWGTUP);
      weight_lep = weight_Lumi * weightSign * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      weight_lepB = weight_Lumi * weightSign * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
      else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
      resetVec(reg);
      if(leptonPass && triggerPass){
        reg[0]=0;
        wgt[0]=weight_lep;
        wcfit[0]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass){
        reg[1]=1;
        wgt[1]=weight_lep;
        wcfit[1]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjet==1){
        reg[2]=2;
        wgt[2]=weight_lepB;
        wcfit[2]=*eft_fit;
        for (int l=0;l<selectedJets->size();l++){
          if((*selectedJets)[l]->btag_){
            recoBjet = (*selectedJets)[l]->p4_;
            break;
          }
        }
        recoL2 = (*selectedLeptonsMuScaleDown)[1]->p4_;
        recoNu = Wneutrino(myMETpt, myMETphi, (*selectedLeptonsMuScaleDown)[1]->pt_, (*selectedLeptonsMuScaleDown)[1]->eta_, (*selectedLeptonsMuScaleDown)[1]->phi_);
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
      }
      if(leptonPass && triggerPass && MetPass && DyPass && nbjet>1){
        reg[3]=3;
        wgt[3]=weight_lepB;
        wcfit[3]=*eft_fit;
      }
      for (int n=0;n<sys.size();++n){
        if (sys[n]!="muonScale") continue;
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Pt",4), n, (*selectedLeptonsMuScaleDown)[0]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Eta",4), n, (*selectedLeptonsMuScaleDown)[0]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Phi",4), n, (*selectedLeptonsMuScaleDown)[0]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Pt",4), n, (*selectedLeptonsMuScaleDown)[1]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Eta",4), n, (*selectedLeptonsMuScaleDown)[1]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Phi",4), n, (*selectedLeptonsMuScaleDown)[1]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llM",4), n, ((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llPt",4), n, ((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).Pt() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llDr",4), n, deltaR((*selectedLeptonsMuScaleDown)[0]->eta_,(*selectedLeptonsMuScaleDown)[0]->phi_,(*selectedLeptonsMuScaleDown)[1]->eta_,(*selectedLeptonsMuScaleDown)[1]->phi_) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llDphi",4), n, abs(deltaPhi((*selectedLeptonsMuScaleDown)[0]->phi_,(*selectedLeptonsMuScaleDown)[1]->phi_)) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"njet",4), n, selectedJets->size() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"nbjet",4), n, nbjetJerDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"Met",4), n, MetJetsPtJerDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"MetPhi",4), n, MetJetsPhiJerDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"nVtx",4), n, PV_npvs ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llMZw",4), n, ((*selectedLeptonsMuScaleDown)[0]->p4_ + (*selectedLeptonsMuScaleDown)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"BDT",4), n, MVAoutput ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topMass",4), n, recoTop.M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1Dphi",4), n, abs(deltaPhi((*selectedLeptonsMuScaleDown)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1Dr",4), n, deltaR((*selectedLeptonsMuScaleDown)[0]->eta_,(*selectedLeptonsMuScaleDown)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1DptOsumPt",4), n, (abs((*selectedLeptonsMuScaleDown)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuScaleDown)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topPt",4), n, recoTop.Pt() ,wgt, wcfit);
        if(selectedJets->size()>0){
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Pt",4), n, (*selectedJets)[0]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Eta",4), n, (*selectedJets)[0]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Phi",4), n, (*selectedJets)[0]->phi_ ,wgt, wcfit);
        }
      }
      delete eft_fit;
    }

    if(selectedLeptonsEleScaleUp->size() ==2){
      leptonPass=false;
      triggerPass=false;
      DyPass = false;
      MetPass = false;
      ch=100;
      if (((*selectedLeptonsEleScaleUp)[0]->pt_ > 25) && ((*selectedLeptonsEleScaleUp)[0]->charge_ * (*selectedLeptonsEleScaleUp)[1]->charge_ == -1) && ((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).M()>20 && ((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).M()<1200) leptonPass=true;
      if ((*selectedLeptonsEleScaleUp)[0]->lep_ + (*selectedLeptonsEleScaleUp)[1]->lep_ == 2) ch = 0;
      if ((*selectedLeptonsEleScaleUp)[0]->lep_ + (*selectedLeptonsEleScaleUp)[1]->lep_ == 11) ch = 1;
      if ((*selectedLeptonsEleScaleUp)[0]->lep_ + (*selectedLeptonsEleScaleUp)[1]->lep_ == 20) ch = 2;
      if (ch ==0 && triggerPassEE) triggerPass=true;
      if (ch ==1 && triggerPassEMu) triggerPass=true;
      if (ch ==2 && triggerPassMuMu) triggerPass=true;
      if (((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).M()>106 ) DyPass = true;
      if (myMETpt > MetCut) MetPass = true;
      if (ch==0)  sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptonsEleScaleUp)[0]->pt_, (*selectedLeptonsEleScaleUp)[1]->pt_,"central",false, true);
      if (ch==1) {
        if((*selectedLeptonsEleScaleUp)[0]->lep_==10) {
          effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsEleScaleUp)[0]->pt_, abs((*selectedLeptonsEleScaleUp)[0]->eta_),"central",false, true);
          effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsEleScaleUp)[0]->pt_, abs((*selectedLeptonsEleScaleUp)[0]->eta_),"central",false, true);
        }
        else{
          effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsEleScaleUp)[1]->pt_, abs((*selectedLeptonsEleScaleUp)[1]->eta_),"central",false, true);
          effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsEleScaleUp)[1]->pt_, abs((*selectedLeptonsEleScaleUp)[1]->eta_),"central",false, true);
        }
        sf_Trigger = effTdata1/effTmc1;
      }
      if (data == "mc" && ch==2) {
        effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsEleScaleUp)[0]->pt_, abs((*selectedLeptonsEleScaleUp)[0]->eta_),"central",false, true);
        effTmc2 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsEleScaleUp)[1]->pt_, abs((*selectedLeptonsEleScaleUp)[1]->eta_),"central",false, true);
        effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsEleScaleUp)[0]->pt_, abs((*selectedLeptonsEleScaleUp)[0]->eta_),"central",false, true);
        effTdata2 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsEleScaleUp)[1]->pt_, abs((*selectedLeptonsEleScaleUp)[1]->eta_),"central",false, true);
        sf_Trigger = ((1-((1-effTdata1)*(1-effTdata2)))/(1-((1-effTmc1)*(1-effTmc2))));
      }
      if(!fname.Contains("pythia")) weightSign = signnum_typical(LHEWeight_originalXWGTUP);
      weight_lep = weight_Lumi * weightSign * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      weight_lepB = weight_Lumi * weightSign * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
      else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
      resetVec(reg);
      if(leptonPass && triggerPass){
        reg[0]=0;
        wgt[0]=weight_lep;
        wcfit[0]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass){
        reg[1]=1;
        wgt[1]=weight_lep;
        wcfit[1]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjet==1){
        reg[2]=2;
        wgt[2]=weight_lepB;
        wcfit[2]=*eft_fit;
        for (int l=0;l<selectedJets->size();l++){
          if((*selectedJets)[l]->btag_){
            recoBjet = (*selectedJets)[l]->p4_;
            break;
          }
        }
        recoL2 = (*selectedLeptonsEleScaleUp)[1]->p4_;
        recoNu = Wneutrino(myMETpt, myMETphi, (*selectedLeptonsEleScaleUp)[1]->pt_, (*selectedLeptonsEleScaleUp)[1]->eta_, (*selectedLeptonsEleScaleUp)[1]->phi_);
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
      }
      if(leptonPass && triggerPass && MetPass && DyPass && nbjet>1){
        reg[3]=3;
        wgt[3]=weight_lepB;
        wcfit[3]=*eft_fit;
      }
      for (int n=0;n<sys.size();++n){
        if (sys[n]!="electronScale") continue;
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Pt",4), n, (*selectedLeptonsEleScaleUp)[0]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Eta",4), n, (*selectedLeptonsEleScaleUp)[0]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Phi",4), n, (*selectedLeptonsEleScaleUp)[0]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Pt",4), n, (*selectedLeptonsEleScaleUp)[1]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Eta",4), n, (*selectedLeptonsEleScaleUp)[1]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Phi",4), n, (*selectedLeptonsEleScaleUp)[1]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llM",4), n, ((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llPt",4), n, ((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).Pt() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llDr",4), n, deltaR((*selectedLeptonsEleScaleUp)[0]->eta_,(*selectedLeptonsEleScaleUp)[0]->phi_,(*selectedLeptonsEleScaleUp)[1]->eta_,(*selectedLeptonsEleScaleUp)[1]->phi_) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llDphi",4), n, abs(deltaPhi((*selectedLeptonsEleScaleUp)[0]->phi_,(*selectedLeptonsEleScaleUp)[1]->phi_)) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"njet",4), n, selectedJets->size() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"nbjet",4), n, nbjetJerUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"Met",4), n, MetJetsPtJerUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"MetPhi",4), n, MetJetsPhiJerUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"nVtx",4), n, PV_npvs ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llMZw",4), n, ((*selectedLeptonsEleScaleUp)[0]->p4_ + (*selectedLeptonsEleScaleUp)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"BDT",4), n, MVAoutput ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topMass",4), n, recoTop.M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1Dphi",4), n, abs(deltaPhi((*selectedLeptonsEleScaleUp)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1Dr",4), n, deltaR((*selectedLeptonsEleScaleUp)[0]->eta_,(*selectedLeptonsEleScaleUp)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1DptOsumPt",4), n, (abs((*selectedLeptonsEleScaleUp)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsEleScaleUp)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topPt",4), n, recoTop.Pt() ,wgt, wcfit);
        if(selectedJets->size()>0){
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Pt",4), n, (*selectedJets)[0]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Eta",4), n, (*selectedJets)[0]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Phi",4), n, (*selectedJets)[0]->phi_ ,wgt, wcfit);
        }
      }
    }

    if(selectedLeptonsEleScaleDown->size() ==2){
      leptonPass=false;
      triggerPass=false;
      DyPass = false;
      MetPass = false;
      ch=100;
      if (((*selectedLeptonsEleScaleDown)[0]->pt_ > 25) && ((*selectedLeptonsEleScaleDown)[0]->charge_ * (*selectedLeptonsEleScaleDown)[1]->charge_ == -1) && ((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).M()>20 && ((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).M()<1200) leptonPass=true;
      if ((*selectedLeptonsEleScaleDown)[0]->lep_ + (*selectedLeptonsEleScaleDown)[1]->lep_ == 2) ch = 0;
      if ((*selectedLeptonsEleScaleDown)[0]->lep_ + (*selectedLeptonsEleScaleDown)[1]->lep_ == 11) ch = 1;
      if ((*selectedLeptonsEleScaleDown)[0]->lep_ + (*selectedLeptonsEleScaleDown)[1]->lep_ == 20) ch = 2;
      if (ch ==0 && triggerPassEE) triggerPass=true;
      if (ch ==1 && triggerPassEMu) triggerPass=true;
      if (ch ==2 && triggerPassMuMu) triggerPass=true;
      if (((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).M()>106 ) DyPass = true;
      if (myMETpt > MetCut) MetPass = true;
      if (ch==0)  sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptonsEleScaleDown)[0]->pt_, (*selectedLeptonsEleScaleDown)[1]->pt_,"central",false, true);
      if (ch==1) {
        if((*selectedLeptonsEleScaleDown)[0]->lep_==10) {
          effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsEleScaleDown)[0]->pt_, abs((*selectedLeptonsEleScaleDown)[0]->eta_),"central",false, true);
          effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsEleScaleDown)[0]->pt_, abs((*selectedLeptonsEleScaleDown)[0]->eta_),"central",false, true);
        }
        else{
          effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsEleScaleDown)[1]->pt_, abs((*selectedLeptonsEleScaleDown)[1]->eta_),"central",false, true);
          effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsEleScaleDown)[1]->pt_, abs((*selectedLeptonsEleScaleDown)[1]->eta_),"central",false, true);
        }
        sf_Trigger = effTdata1/effTmc1;
      }
      if (data == "mc" && ch==2) {
        effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsEleScaleDown)[0]->pt_, abs((*selectedLeptonsEleScaleDown)[0]->eta_),"central",false, true);
        effTmc2 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsEleScaleDown)[1]->pt_, abs((*selectedLeptonsEleScaleDown)[1]->eta_),"central",false, true);
        effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsEleScaleDown)[0]->pt_, abs((*selectedLeptonsEleScaleDown)[0]->eta_),"central",false, true);
        effTdata2 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsEleScaleDown)[1]->pt_, abs((*selectedLeptonsEleScaleDown)[1]->eta_),"central",false, true);
        sf_Trigger = ((1-((1-effTdata1)*(1-effTdata2)))/(1-((1-effTmc1)*(1-effTmc2))));
      }
      if(!fname.Contains("pythia")) weightSign = signnum_typical(LHEWeight_originalXWGTUP);
      weight_lep = weight_Lumi * weightSign * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      weight_lepB = weight_Lumi * weightSign * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
      else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
      resetVec(reg);
      if(leptonPass && triggerPass){
        reg[0]=0;
        wgt[0]=weight_lep;
        wcfit[0]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass){
        reg[1]=1;
        wgt[1]=weight_lep;
        wcfit[1]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjet==1){
        reg[2]=2;
        wgt[2]=weight_lepB;
        wcfit[2]=*eft_fit;
        for (int l=0;l<selectedJets->size();l++){
          if((*selectedJets)[l]->btag_){
            recoBjet = (*selectedJets)[l]->p4_;
            break;
          }
        }
        recoL2 = (*selectedLeptonsEleScaleDown)[1]->p4_;
        recoNu = Wneutrino(myMETpt, myMETphi, (*selectedLeptonsEleScaleDown)[1]->pt_, (*selectedLeptonsEleScaleDown)[1]->eta_, (*selectedLeptonsEleScaleDown)[1]->phi_);
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
      }
      if(leptonPass && triggerPass && MetPass && DyPass && nbjet>1){
        reg[3]=3;
        wgt[3]=weight_lepB;
        wcfit[3]=*eft_fit;
      }
      for (int n=0;n<sys.size();++n){
        if (sys[n]!="electronScale") continue;
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Pt",4), n, (*selectedLeptonsEleScaleDown)[0]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Eta",4), n, (*selectedLeptonsEleScaleDown)[0]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Phi",4), n, (*selectedLeptonsEleScaleDown)[0]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Pt",4), n, (*selectedLeptonsEleScaleDown)[1]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Eta",4), n, (*selectedLeptonsEleScaleDown)[1]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Phi",4), n, (*selectedLeptonsEleScaleDown)[1]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llM",4), n, ((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llPt",4), n, ((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).Pt() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llDr",4), n, deltaR((*selectedLeptonsEleScaleDown)[0]->eta_,(*selectedLeptonsEleScaleDown)[0]->phi_,(*selectedLeptonsEleScaleDown)[1]->eta_,(*selectedLeptonsEleScaleDown)[1]->phi_) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llDphi",4), n, abs(deltaPhi((*selectedLeptonsEleScaleDown)[0]->phi_,(*selectedLeptonsEleScaleDown)[1]->phi_)) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"njet",4), n, selectedJets->size() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"nbjet",4), n, nbjetJerDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"Met",4), n, MetJetsPtJerDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"MetPhi",4), n, MetJetsPhiJerDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"nVtx",4), n, PV_npvs ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llMZw",4), n, ((*selectedLeptonsEleScaleDown)[0]->p4_ + (*selectedLeptonsEleScaleDown)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"BDT",4), n, MVAoutput ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topMass",4), n, recoTop.M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1Dphi",4), n, abs(deltaPhi((*selectedLeptonsEleScaleDown)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1Dr",4), n, deltaR((*selectedLeptonsEleScaleDown)[0]->eta_,(*selectedLeptonsEleScaleDown)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1DptOsumPt",4), n, (abs((*selectedLeptonsEleScaleDown)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsEleScaleDown)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topPt",4), n, recoTop.Pt() ,wgt, wcfit);
        if(selectedJets->size()>0){
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Pt",4), n, (*selectedJets)[0]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Eta",4), n, (*selectedJets)[0]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Phi",4), n, (*selectedJets)[0]->phi_ ,wgt, wcfit);
        }
      }
      delete eft_fit;
    }


    if(selectedLeptonsMuResUp->size() ==2){
      leptonPass=false;
      triggerPass=false;
      DyPass = false;
      MetPass = false;
      ch=100;
      if (((*selectedLeptonsMuResUp)[0]->pt_ > 25) && ((*selectedLeptonsMuResUp)[0]->charge_ * (*selectedLeptonsMuResUp)[1]->charge_ == -1) && ((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).M()>20 && ((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).M()<1200) leptonPass=true;
      if ((*selectedLeptonsMuResUp)[0]->lep_ + (*selectedLeptonsMuResUp)[1]->lep_ == 2) ch = 0;
      if ((*selectedLeptonsMuResUp)[0]->lep_ + (*selectedLeptonsMuResUp)[1]->lep_ == 11) ch = 1;
      if ((*selectedLeptonsMuResUp)[0]->lep_ + (*selectedLeptonsMuResUp)[1]->lep_ == 20) ch = 2;
      if (ch ==0 && triggerPassEE) triggerPass=true;
      if (ch ==1 && triggerPassEMu) triggerPass=true;
      if (ch ==2 && triggerPassMuMu) triggerPass=true;
      if (((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).M()>106 ) DyPass = true;
      if (myMETpt > MetCut) MetPass = true;
      if (ch==0)  sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptonsMuResUp)[0]->pt_, (*selectedLeptonsMuResUp)[1]->pt_,"central",false, true);
      if (ch==1) {
        if((*selectedLeptonsMuResUp)[0]->lep_==10) {
          effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuResUp)[0]->pt_, abs((*selectedLeptonsMuResUp)[0]->eta_),"central",false, true);
          effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuResUp)[0]->pt_, abs((*selectedLeptonsMuResUp)[0]->eta_),"central",false, true);
        }
        else{
          effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuResUp)[1]->pt_, abs((*selectedLeptonsMuResUp)[1]->eta_),"central",false, true);
          effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuResUp)[1]->pt_, abs((*selectedLeptonsMuResUp)[1]->eta_),"central",false, true);
        }
        sf_Trigger = effTdata1/effTmc1;
      }
      if (data == "mc" && ch==2) {
        effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuResUp)[0]->pt_, abs((*selectedLeptonsMuResUp)[0]->eta_),"central",false, true);
        effTmc2 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuResUp)[1]->pt_, abs((*selectedLeptonsMuResUp)[1]->eta_),"central",false, true);
        effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuResUp)[0]->pt_, abs((*selectedLeptonsMuResUp)[0]->eta_),"central",false, true);
        effTdata2 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuResUp)[1]->pt_, abs((*selectedLeptonsMuResUp)[1]->eta_),"central",false, true);
        sf_Trigger = ((1-((1-effTdata1)*(1-effTdata2)))/(1-((1-effTmc1)*(1-effTmc2))));
      }
      if(!fname.Contains("pythia")) weightSign = signnum_typical(LHEWeight_originalXWGTUP);
      weight_lep = weight_Lumi * weightSign * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      weight_lepB = weight_Lumi * weightSign * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
      else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
      resetVec(reg);
      if(leptonPass && triggerPass){
        reg[0]=0;
        wgt[0]=weight_lep;
        wcfit[0]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass){
        reg[1]=1;
        wgt[1]=weight_lep;
        wcfit[1]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjet==1){
        reg[2]=2;
        wgt[2]=weight_lepB;
        wcfit[2]=*eft_fit;
        for (int l=0;l<selectedJets->size();l++){
          if((*selectedJets)[l]->btag_){
            recoBjet = (*selectedJets)[l]->p4_;
            break;
          }
        }
        recoL2 = (*selectedLeptonsMuResUp)[1]->p4_;
        recoNu = Wneutrino(myMETpt, myMETphi, (*selectedLeptonsMuResUp)[1]->pt_, (*selectedLeptonsMuResUp)[1]->eta_, (*selectedLeptonsMuResUp)[1]->phi_);
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
      }
      if(leptonPass && triggerPass && MetPass && DyPass && nbjet>1){
        reg[3]=3;
        wgt[3]=weight_lepB;
        wcfit[3]=*eft_fit;
      }
      for (int n=0;n<sys.size();++n){
        if (sys[n]!="muonRes") continue;
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Pt",4), n, (*selectedLeptonsMuResUp)[0]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Eta",4), n, (*selectedLeptonsMuResUp)[0]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep1Phi",4), n, (*selectedLeptonsMuResUp)[0]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Pt",4), n, (*selectedLeptonsMuResUp)[1]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Eta",4), n, (*selectedLeptonsMuResUp)[1]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"lep2Phi",4), n, (*selectedLeptonsMuResUp)[1]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llM",4), n, ((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llPt",4), n, ((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).Pt() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llDr",4), n, deltaR((*selectedLeptonsMuResUp)[0]->eta_,(*selectedLeptonsMuResUp)[0]->phi_,(*selectedLeptonsMuResUp)[1]->eta_,(*selectedLeptonsMuResUp)[1]->phi_) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llDphi",4), n, abs(deltaPhi((*selectedLeptonsMuResUp)[0]->phi_,(*selectedLeptonsMuResUp)[1]->phi_)) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"njet",4), n, selectedJets->size() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"nbjet",4), n, nbjetJerUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"Met",4), n, MetJetsPtJerUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"MetPhi",4), n, MetJetsPhiJerUp ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"nVtx",4), n, PV_npvs ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"llMZw",4), n, ((*selectedLeptonsMuResUp)[0]->p4_ + (*selectedLeptonsMuResUp)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"BDT",4), n, MVAoutput ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topMass",4), n, recoTop.M() ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1Dphi",4), n, abs(deltaPhi((*selectedLeptonsMuResUp)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1Dr",4), n, deltaR((*selectedLeptonsMuResUp)[0]->eta_,(*selectedLeptonsMuResUp)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topL1DptOsumPt",4), n, (abs((*selectedLeptonsMuResUp)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuResUp)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
        FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"topPt",4), n, recoTop.Pt() ,wgt, wcfit);
        if(selectedJets->size()>0){
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Pt",4), n, (*selectedJets)[0]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Eta",4), n, (*selectedJets)[0]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysUp, ch, reg, vInd(vars,"jet1Phi",4), n, (*selectedJets)[0]->phi_ ,wgt, wcfit);
        }
      }
    }

    if(selectedLeptonsMuResDown->size() ==2){
      leptonPass=false;
      triggerPass=false;
      DyPass = false;
      MetPass = false;
      ch=100;
      if (((*selectedLeptonsMuResDown)[0]->pt_ > 25) && ((*selectedLeptonsMuResDown)[0]->charge_ * (*selectedLeptonsMuResDown)[1]->charge_ == -1) && ((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).M()>20 && ((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).M()<1200) leptonPass=true;
      if ((*selectedLeptonsMuResDown)[0]->lep_ + (*selectedLeptonsMuResDown)[1]->lep_ == 2) ch = 0;
      if ((*selectedLeptonsMuResDown)[0]->lep_ + (*selectedLeptonsMuResDown)[1]->lep_ == 11) ch = 1;
      if ((*selectedLeptonsMuResDown)[0]->lep_ + (*selectedLeptonsMuResDown)[1]->lep_ == 20) ch = 2;
      if (ch ==0 && triggerPassEE) triggerPass=true;
      if (ch ==1 && triggerPassEMu) triggerPass=true;
      if (ch ==2 && triggerPassMuMu) triggerPass=true;
      if (((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).M()>106 ) DyPass = true;
      if (myMETpt > MetCut) MetPass = true;
      if (ch==0)  sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptonsMuResDown)[0]->pt_, (*selectedLeptonsMuResDown)[1]->pt_,"central",false, true);
      if (ch==1) {
        if((*selectedLeptonsMuResDown)[0]->lep_==10) {
          effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuResDown)[0]->pt_, abs((*selectedLeptonsMuResDown)[0]->eta_),"central",false, true);
          effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuResDown)[0]->pt_, abs((*selectedLeptonsMuResDown)[0]->eta_),"central",false, true);
        }
        else{
          effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuResDown)[1]->pt_, abs((*selectedLeptonsMuResDown)[1]->eta_),"central",false, true);
          effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuResDown)[1]->pt_, abs((*selectedLeptonsMuResDown)[1]->eta_),"central",false, true);
        }
        sf_Trigger = effTdata1/effTmc1;
      }
      if (data == "mc" && ch==2) {
        effTmc1 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuResDown)[0]->pt_, abs((*selectedLeptonsMuResDown)[0]->eta_),"central",false, true);
        effTmc2 = scale_factor(&eff_triggermumu_mc_H, (*selectedLeptonsMuResDown)[1]->pt_, abs((*selectedLeptonsMuResDown)[1]->eta_),"central",false, true);
        effTdata1 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuResDown)[0]->pt_, abs((*selectedLeptonsMuResDown)[0]->eta_),"central",false, true);
        effTdata2 = scale_factor(&eff_triggermumu_data_H, (*selectedLeptonsMuResDown)[1]->pt_, abs((*selectedLeptonsMuResDown)[1]->eta_),"central",false, true);
        sf_Trigger = ((1-((1-effTdata1)*(1-effTdata2)))/(1-((1-effTmc1)*(1-effTmc2))));
      }
      if(!fname.Contains("pythia")) weightSign = signnum_typical(LHEWeight_originalXWGTUP);
      weight_lep = weight_Lumi * weightSign * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      weight_lepB = weight_Lumi * weightSign * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      weight_EFT = lumi * (1000.0/nRuns)  * sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO* (P_bjet_data/P_bjet_mc) * sf_Trigger * weight_PU * weight_prefiring * weight_topPtPowhegData * ttKFactor * sf_JetPuId;
      if (iseft) eft_fit = new WCFit(nWCnames, wc_names_lst, nEFTfitCoefficients, EFTfitCoefficients, weight_EFT);
      else eft_fit = new WCFit(0,wc_names_lst,1, &genWeight, 1.0);
      resetVec(reg);
      if(leptonPass && triggerPass){
        reg[0]=0;
        wgt[0]=weight_lep;
        wcfit[0]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass){
        reg[1]=1;
        wgt[1]=weight_lep;
        wcfit[1]=*eft_fit;
      }
      if(leptonPass && triggerPass && DyPass && MetPass && nbjet==1){
        reg[2]=2;
        wgt[2]=weight_lepB;
        wcfit[2]=*eft_fit;
        for (int l=0;l<selectedJets->size();l++){
          if((*selectedJets)[l]->btag_){
            recoBjet = (*selectedJets)[l]->p4_;
            break;
          }
        }
        recoL2 = (*selectedLeptonsMuResDown)[1]->p4_;
        recoNu = Wneutrino(myMETpt, myMETphi, (*selectedLeptonsMuResDown)[1]->pt_, (*selectedLeptonsMuResDown)[1]->eta_, (*selectedLeptonsMuResDown)[1]->phi_);
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
      }
      if(leptonPass && triggerPass && MetPass && DyPass && nbjet>1){
        reg[3]=3;
        wgt[3]=weight_lepB;
        wcfit[3]=*eft_fit;
      }
      for (int n=0;n<sys.size();++n){
        if (sys[n]!="muonRes") continue;
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Pt",4), n, (*selectedLeptonsMuResDown)[0]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Eta",4), n, (*selectedLeptonsMuResDown)[0]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep1Phi",4), n, (*selectedLeptonsMuResDown)[0]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Pt",4), n, (*selectedLeptonsMuResDown)[1]->pt_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Eta",4), n, (*selectedLeptonsMuResDown)[1]->eta_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"lep2Phi",4), n, (*selectedLeptonsMuResDown)[1]->phi_ ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llM",4), n, ((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llPt",4), n, ((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).Pt() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llDr",4), n, deltaR((*selectedLeptonsMuResDown)[0]->eta_,(*selectedLeptonsMuResDown)[0]->phi_,(*selectedLeptonsMuResDown)[1]->eta_,(*selectedLeptonsMuResDown)[1]->phi_) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llDphi",4), n, abs(deltaPhi((*selectedLeptonsMuResDown)[0]->phi_,(*selectedLeptonsMuResDown)[1]->phi_)) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"njet",4), n, selectedJets->size() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"nbjet",4), n, nbjetJerDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"Met",4), n, MetJetsPtJerDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"MetPhi",4), n, MetJetsPhiJerDown ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"nVtx",4), n, PV_npvs ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"llMZw",4), n, ((*selectedLeptonsMuResDown)[0]->p4_ + (*selectedLeptonsMuResDown)[1]->p4_).M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"BDT",4), n, MVAoutput ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topMass",4), n, recoTop.M() ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1Dphi",4), n, abs(deltaPhi((*selectedLeptonsMuResDown)[0]->phi_,recoTop.Phi())) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1Dr",4), n, deltaR((*selectedLeptonsMuResDown)[0]->eta_,(*selectedLeptonsMuResDown)[0]->phi_,recoTop.Eta(),recoTop.Phi()) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topL1DptOsumPt",4), n, (abs((*selectedLeptonsMuResDown)[0]->pt_ - recoTop.Pt()))/((*selectedLeptonsMuResDown)[0]->pt_ + recoTop.Pt()) ,wgt, wcfit);
        FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"topPt",4), n, recoTop.Pt() ,wgt, wcfit);
        if(selectedJets->size()>0){
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Pt",4), n, (*selectedJets)[0]->pt_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Eta",4), n, (*selectedJets)[0]->eta_ ,wgt, wcfit);
          FillD4Hists(HistsSysDown, ch, reg, vInd(vars,"jet1Phi",4), n, (*selectedJets)[0]->phi_ ,wgt, wcfit);
        }
      }
      delete eft_fit;

    }
/*
      sumPuWeight+=weight_PU;
      sumPreFireWeight+=weight_prefiring;
      sumWeightMuID+=sf_Mu_ID;
      sumWeightMuIso+=sf_Mu_ISO;
      sumWeighttTrigger+=sf_Trigger;
      sumWeightEleID+=sf_Ele_ID;
      sumWeightEleReco+=sf_Ele_Reco;
      sumWeightBtag+=(P_bjet_data/P_bjet_mc);
*/

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
    JECsysDown->clear();
    JECsysNbtagUp->clear();
    JECsysNbtagDown->clear();
    JECsysMETUp->clear();
    JECsysMETDown->clear();
    JECsysMETPhiUp->clear();
    JECsysMETPhiDown->clear();
    JECsysUp->shrink_to_fit();
    JECsysDown->shrink_to_fit();
    JECsysNbtagUp->shrink_to_fit();
    JECsysNbtagDown->shrink_to_fit();
    JECsysMETUp->shrink_to_fit();
    JECsysMETDown->shrink_to_fit();
    JECsysMETPhiUp->shrink_to_fit();
    JECsysMETPhiDown->shrink_to_fit();
    delete JECsysUp;
    delete JECsysDown;
    delete JECsysNbtagUp;
    delete JECsysNbtagDown;
    delete JECsysMETUp;
    delete JECsysMETDown;
    delete JECsysMETPhiUp;
    delete JECsysMETPhiDown;
  }
  cout<<"Loop is completed"<<endl;
  cout<<"from "<<ntr<<" events, "<<nAccept<<" events are accepted"<<endl;
  cout<<"Total Virtual Memory just after the loop: "<<(getValue()-memoryInit)/1000.0<<" MB"<<endl;
//  cout<<"sumPuWeight "<<sumPuWeight<<endl;
//  cout<<"sumPreFireWeight "<<sumPreFireWeight<<endl;
//  cout<<"sumWeightMuID "<<sumWeightMuID<<endl;
//  cout<<"sumWeightMuIso "<<sumWeightMuIso<<endl;
//  cout<<"sumWeighttTrigger "<<sumWeighttTrigger<<endl;
//  cout<<"sumWeightEleReco "<<sumWeightEleReco<<endl;
//  cout<<"sumWeightEleID "<<sumWeightEleID<<endl;  
//  cout<<"sumWeightBtag "<<sumWeightBtag<<endl;
//cout<<"Integral of the region "<<regions[getVecPos(regions,"llB1-BmuExc")]<<":"<<Hists[3][getVecPos(regions,"llB1-BmuExc")][vInd(vars,"lep1Pt")]->Integral()<<endl;
//cout<<"Integral of the region "<<regions[getVecPos(regions,"llB1")]<<":"<<Hists[3][getVecPos(regions,"llB1")][vInd(vars,"lep1Pt")]->Integral()<<endl;

  TFile file_out ("ANoutput.root","RECREATE");
  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int l=0;l<vars.size();++l){
        Hists[i][k][l]  ->Write("",TObject::kOverwrite);
      }
    }
  }
  cout<<"Cleaning the main hists"<<endl;
  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int l=0;l<vars.size();++l){
        delete Hists[i][k][l];
      }
    }
  }

  if(data=="mc" && !fname.Contains("sys")  && ifSys){
    for (int i=0;i<channels.size();++i){
      file_out.mkdir("sys"+channels[i]);
      file_out.cd("sys"+channels[i]+"/");
      for (int k=0;k<regions.size();++k){
        for (int l=0;l<varsWithSys.size();++l){
          for (int n=0;n<sys.size();++n){
            HistsSysUp[i][k][l][n]->Write("",TObject::kOverwrite);
            HistsSysDown[i][k][l][n]->Write("",TObject::kOverwrite);
          }
        }
      }
      file_out.cd("");
    }
  }
  file_out.cd("");
  h2_BTaggingEff_Denom_b   ->Write("",TObject::kOverwrite);
  h2_BTaggingEff_Denom_c   ->Write("",TObject::kOverwrite);
  h2_BTaggingEff_Denom_udsg->Write("",TObject::kOverwrite);
  h2_BTaggingEff_Num_b     ->Write("",TObject::kOverwrite);
  h2_BTaggingEff_Num_c     ->Write("",TObject::kOverwrite);
  h2_BTaggingEff_Num_udsg  ->Write("",TObject::kOverwrite);
  crossSection             ->Write("",TObject::kOverwrite);
  eleEffNum                ->Write("",TObject::kOverwrite);
  eleEffDen                ->Write("",TObject::kOverwrite);
  truePU                   ->Write("",TObject::kOverwrite);

  if(data=="mc" && !fname.Contains("sys")  && ifSys){
    file_out.mkdir("JECSys");
    file_out.cd("JECSys/");
    for (int i=0;i<channels.size();++i){
      for (int k=2;k<regions.size();++k){
        for (int n=0;n<sysJecNames.size();++n){
          HistsJecUp[i][k-2][0][n]->Write("",TObject::kOverwrite);
          HistsJecDown[i][k-2][0][n]->Write("",TObject::kOverwrite);
        }
      }
    }
  }
  file_out.cd("");
  if(ifSys){
    if ((fname.Contains("TTTo2L2Nu")&& !fname.Contains("sys")) || fname.Contains("BNV")){
      file_out.mkdir("reweightingSys");
      file_out.cd("reweightingSys/");
      for (int i=0;i<channels.size();++i){
        for (int k=2;k<regions.size();++k){
          for (int n=0;n<nScale;++n){
            HistsSysReweightsQscale[i][k-2][0][n]->Write("",TObject::kOverwrite);
          }
          for (int n=0;n<nPdf;++n){
            HistsSysReweightsPDF[i][k-2][0][n]->Write("",TObject::kOverwrite);
          }
          for (int n=0;n<nPS;++n){
            HistsSysReweightsPS[i][k-2][0][n]->Write("",TObject::kOverwrite);
          }
        }
      }
    }
    for (int i=0;i<channels.size();++i){
      for (int k=0;k<regions.size();++k){
        for (int l=0;l<varsWithSys.size();++l){
          for (int n=0;n<sys.size();++n){
            delete HistsSysUp[i][k][l][n];
            delete HistsSysDown[i][k][l][n];
          }
        }
      }
    }

    if(data=="mc" && !fname.Contains("sys")){
      for (int i=0;i<channels.size();++i){
        for (int k=2;k<regions.size();++k){
          for (int n=0;n<sysJecNames.size();++n){
            delete HistsJecUp[i][k-2][0][n];
            delete HistsJecDown[i][k-2][0][n];
          }
        }
      }
    }
    if ((fname.Contains("TTTo2L2Nu")&& !fname.Contains("sys")) || fname.Contains("BNV")){
      for (int i=0;i<channels.size();++i){
        for (int k=2;k<regions.size();++k){
          for (int n=0;n<nScale;++n){
            delete HistsSysReweightsQscale[i][k-2][0][n];
          }
          for (int n=0;n<nPdf;++n){
            delete HistsSysReweightsPDF[i][k-2][0][n];
          }
          for (int n=0;n<nPS;++n){
            delete HistsSysReweightsPS[i][k-2][0][n];
          }
        }
      }
    }
  }
  file_out.cd("");
  tree_out.Write() ;
  file_out.Close() ;
  Hists.clear();
  Hists.shrink_to_fit();
//  TTptQS_TG.~TGraphAsymmErrors();
  cout<<"Hists cleaned"<<endl;
  if(ifSys){
    HistsSysUp.clear();
    HistsSysUp.shrink_to_fit();
    cout<<"HistsSysUp cleaned"<<endl;
    HistsSysDown.clear();
    HistsSysDown.shrink_to_fit();
    cout<<"HistsSysDown cleaned"<<endl;
    HistsJecUp.clear();
    HistsJecUp.shrink_to_fit();
    cout<<"HistsJecUp cleaned"<<endl;
    HistsJecDown.clear();
    HistsJecDown.shrink_to_fit();
    cout<<"HistsJecDown cleaned"<<endl;
    HistsSysReweightsPDF.clear();
    HistsSysReweightsPDF.shrink_to_fit();
    cout<<"HistsSysReweightPDF cleaned"<<endl;
    HistsSysReweightsQscale.clear();
    HistsSysReweightsQscale.shrink_to_fit();
    cout<<"HistsSysReweightsQscale cleaned"<<endl;
    HistsSysReweightsPS.clear();
    HistsSysReweightsPS.shrink_to_fit();
    cout<<"HistsSysReweightPS cleaned"<<endl;
    cout<<"Job is finished"<<endl;
  }
cout<<"Total Virtual Memory: "<<(getValue()-memoryInit)/1000.0<<" MB"<<endl;

}


void MyAnalysis::FillD3Hists(D3HistsContainer H3, int v1, std::vector<int> v2, int v3, float value, std::vector<float> weight, std::vector<WCFit> wcfit){
  for (int i = 0; i < v2.size(); ++i) {
    if(v2[i]>=0) H3[v1][v2[i]][v3]->Fill(value, weight[i], wcfit[i]);
    if(v2[i]>=0) H3[3][v2[i]][v3]->Fill(value, weight[i], wcfit[i]);
  }
}

void MyAnalysis::FillD4Hists(D4HistsContainer H4, int v1, std::vector<int> v2, int v3, int v4, float value, std::vector<float> weight, std::vector<WCFit> wcfit){
  for (int i = 0; i < v2.size(); ++i) {
    if(v2[i]>=0 && v3>=0) H4[v1][v2[i]][v3][v4]->Fill(value, weight[i], wcfit[i]);
    if(v2[i]>=0 && v3>=0) H4[3][v2[i]][v3][v4]->Fill(value, weight[i], wcfit[i]);
  }
}
