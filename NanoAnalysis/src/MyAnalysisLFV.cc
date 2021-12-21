#define MyAnalysis_cxx
#include "MyAnalysis.h"
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
#include "GEScaleSyst.h"
#endif

void displayProgress(long current, long max){
  using std::cerr;
  if (max<2500) return;
  if (current%(max/2500)!=0 && current<max-1) return;

  int width = 52; // Hope the terminal is at least that wide.
  int barWidth = width - 2;
  cerr << "\x1B[2K"; // Clear line
  cerr << "\x1B[2000D"; // Cursor left
  cerr << '[';
  for(int i=0 ; i<barWidth ; ++i){ if(i<barWidth*current/max){ cerr << '=' ; }else{ cerr << ' ' ; } }
  cerr << ']';
  cerr << " " << Form("%8d/%8d (%5.2f%%)", (int)current, (int)max, 100.0*current/max) ;
  cerr.flush();
}

Double_t deltaPhi(Double_t phi1, Double_t phi2) {
  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
  return dPhi;
}


Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {
  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);
  return sqrt(dEta*dEta+dPhi*dPhi);
}

bool ComparePtLep(lepton_candidate *a, lepton_candidate *b) { return a->pt_ > b->pt_; }
bool ComparePtJet(jet_candidate *a, jet_candidate *b) { return a->pt_ > b->pt_; }


float scale_factor( TH2F* h, float X, float Y , TString uncert, bool eff=false, bool out=false){
  int NbinsX=h->GetXaxis()->GetNbins();
  int NbinsY=h->GetYaxis()->GetNbins();
  float x_min=h->GetXaxis()->GetBinLowEdge(1);
  float x_max=h->GetXaxis()->GetBinLowEdge(NbinsX)+h->GetXaxis()->GetBinWidth(NbinsX);
  float y_min=h->GetYaxis()->GetBinLowEdge(1);
  float y_max=h->GetYaxis()->GetBinLowEdge(NbinsY)+h->GetYaxis()->GetBinWidth(NbinsY);
  TAxis *Xaxis = h->GetXaxis();
  TAxis *Yaxis = h->GetYaxis();
  Int_t binx=1;
  Int_t biny=1;
  if(x_min < X && X < x_max) binx = Xaxis->FindBin(X);
  else binx= (X<=x_min) ? 1 : NbinsX ;
  if(y_min < Y && Y < y_max) biny = Yaxis->FindBin(Y);
  else biny= (Y<=y_min) ? 1 : NbinsY ;
// for efficiency return central value
//  if(eff) return  h->GetBinContent(binx, biny);
// if out of range return 1 +- last bin error
//  if(out){
//    if ((X<x_min || X > x_max || Y < y_min || Y > y_max) && uncert=="central") return 1;
//    if ((X<x_min || X > x_max || Y < y_min || Y > y_max) && uncert=="up") return (1+h->GetBinError(binx, biny));
//    if ((X<x_min || X > x_max || Y < y_min || Y > y_max) && uncert=="down") return (1-h->GetBinError(binx, biny));
//  }
// use SF from the histogram and send out SF +_ error
  if(uncert=="up") return (h->GetBinContent(binx, biny)+h->GetBinError(binx, biny));
  if(uncert=="down") return (h->GetBinContent(binx, biny)-h->GetBinError(binx, biny));
  if(uncert=="central") return  h->GetBinContent(binx, biny);
}

float topPtPowheg(float pt){
  return (0.973 - (0.000134 * pt) + (0.103 * exp(pt * (-0.0118))));  
}

//float topPtPowheg(float pt){
//  return (exp(0.0615 - pt*0.0005));
//}

float topPtMGLO(float x){
  return (0.688 -  0.0000174*x + 0.824*exp(-0.0000253*x)/(pow(x,0.2185)));
}

void MyAnalysis::Loop(TString fname, TString sname, TString data, TString dataset ,TString year, TString run, float xs, float lumi, float Nevent)
{
  TRandom3 Tr;

//JEC sources 
  std::string JECFile;
  if(year == "2016")    JECFile = "/user/rgoldouz/NewAnalysis2020/Analysis/input/Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt";
  if(year == "2017")    JECFile = "/user/rgoldouz/NewAnalysis2020/Analysis/input/Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt";
  if(year == "2018")    JECFile = "/user/rgoldouz/NewAnalysis2020/Analysis/input/Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt";

    std::vector<TString> sysJecNames{"AbsoluteMPFBias","AbsoluteScale","AbsoluteStat","FlavorQCD","Fragmentation","PileUpDataMC","PileUpPtBB","PileUpPtEC1","PileUpPtEC2","PileUpPtHF","PileUpPtRef","RelativeFSR","RelativeJEREC1","RelativeJEREC2","RelativeJERHF","RelativePtBB","RelativePtEC1","RelativePtEC2","RelativePtHF","RelativeBal","RelativeSample","RelativeStatEC","RelativeStatFSR","RelativeStatHF","SinglePionECAL","SinglePionHCAL","TimePtEta"};
  const int nsrc = 27;
  const char* srcnames[nsrc] = {"AbsoluteMPFBias","AbsoluteScale","AbsoluteStat","FlavorQCD","Fragmentation","PileUpDataMC","PileUpPtBB","PileUpPtEC1","PileUpPtEC2","PileUpPtHF","PileUpPtRef","RelativeFSR","RelativeJEREC1","RelativeJEREC2","RelativeJERHF","RelativePtBB","RelativePtEC1","RelativePtEC2","RelativePtHF","RelativeBal","RelativeSample","RelativeStatEC","RelativeStatFSR","RelativeStatHF","SinglePionECAL","SinglePionHCAL","TimePtEta"};
  std::vector<JetCorrectionUncertainty*> vsrc(nsrc);
  for (int isrc = 0; isrc < nsrc; isrc++) {
    JetCorrectorParameters *p = new JetCorrectorParameters(JECFile, srcnames[isrc]);
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
    vsrc[isrc] = unc;
  }

//MVA setting
   TMVA::Tools::Instance();
   TMVA::Reader *readerMVA = new TMVA::Reader( "!Color:!Silent" );
   Float_t leading_pt, jet_leading_pt, deltaR_ll, MET, n_jet;
 
   readerMVA->AddVariable ("leading_pt", &leading_pt);
   readerMVA->AddVariable ("jet_leading_pt", &jet_leading_pt);
   readerMVA->AddVariable ("deltaR_ll", &deltaR_ll);
   readerMVA->AddVariable ("MET", &MET);
   readerMVA->AddVariable ("n_jet", &n_jet);
   readerMVA->BookMVA( "BDT_1b_tc_BDT", "/user/rgoldouz/NewAnalysis2020/Analysis/input/TMVA_BDT_1b_BDT.weights.xml");

   Double_t ptBins[11] = {30., 40., 60., 80., 100., 150., 200., 300., 400., 500., 1000.};
   Double_t etaBins [4]= {0., 0.6, 1.2, 2.4};
   TH2D *h2_BTaggingEff_Denom_b    = new TH2D("h2_BTaggingEff_Denom_b"   , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Denom_c    = new TH2D("h2_BTaggingEff_Denom_c"   , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Denom_udsg = new TH2D("h2_BTaggingEff_Denom_udsg", ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Num_b      = new TH2D("h2_BTaggingEff_Num_b"     , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Num_c      = new TH2D("h2_BTaggingEff_Num_c"     , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Num_udsg   = new TH2D("h2_BTaggingEff_Num_udsg"  , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins); 


  typedef vector<std::shared_ptr<TH1F>> Dim1;
  typedef vector<Dim1> Dim2;
  typedef vector<Dim2> Dim3;
  typedef vector<Dim3> Dim4;

  std::vector<TString> regions{"ll","llOffZ","llB1", "llBg1", };
  std::vector<TString> channels{"ee", "emu", "mumu"};
  std::vector<TString> vars   {"lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","llM","llPt","llDr","llDphi","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx", "llMZw","BDT","muPt","elePt","muEta","eleEta"};
  std::vector<int>    nbins   {60      ,20       ,25       ,25      ,20       ,25       ,30   ,20    ,25    ,15      ,20      ,20       ,25       ,10    ,6      ,30   ,20      ,70    ,80      ,100, 70,70,20,20};   
  std::vector<float> lowEdge  {0       ,-3       ,-4       ,0       ,-3       ,-4       ,0    ,0     ,0     ,0       ,0       ,-3       ,-4       ,0     ,0      ,0    ,-4      ,0     ,70      ,-0.4, 0,0,-3,-3};
  std::vector<float> highEdge {1500     ,3        ,4        ,1000     ,3        ,4        ,500  ,200   ,7     ,4       ,300     ,3        ,4        ,10    ,6      ,210  ,4       ,70    ,110,     0.6 ,1500,1500,3,3};       


  Dim3 Hists(channels.size(),Dim2(regions.size(),Dim1(vars.size())));  
  std::stringstream name;
//  TH1F *h_test;
  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int l=0;l<vars.size();++l){
        name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l];
        std::shared_ptr<TH1F> h_test(new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]));        
        h_test->StatOverflows(kTRUE);
        h_test->Sumw2(kTRUE);
        Hists[i][k][l] = h_test;
        name.str("");
      }
    }
  }

  std::vector<TString> sys{"eleRecoSf", "eleIDSf", "muIdSf", "muIsoSf", "bcTagSF", "udsgTagSF","pu", "prefiring", "trigSF","jes", "jer","unclusMET","muonScale","electronScale","muonRes" };
  Dim4 HistsSysUp(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(sys.size()))));
  Dim4 HistsSysDown(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(sys.size()))));

  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int l=0;l<vars.size();++l){
        for (int n=0;n<sys.size();++n){
          name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l]<<"_"<<sys[n]<<"_Up";
          std::shared_ptr<TH1F> h_test(new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]));
          h_test->StatOverflows(kTRUE);
          h_test->Sumw2(kTRUE);
          HistsSysUp[i][k][l][n] = h_test;
          name.str("");
          name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l]<<"_"<<sys[n]<<"_Down";
          std::shared_ptr<TH1F> h_test2(new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]));
          h_test2->StatOverflows(kTRUE);
          h_test2->Sumw2(kTRUE);
          HistsSysDown[i][k][l][n] = h_test2;
          name.str("");
        }
      }
    }
  }

  Dim4 HistsJECUp(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(sysJecNames.size()))));
  Dim4 HistsJECDown(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(sysJecNames.size()))));

  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int l=0;l<vars.size();++l){
        for (int n=0;n<sysJecNames.size();++n){
          name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l]<<"_"<<sysJecNames[n]<<"_Up";
          std::shared_ptr<TH1F> h_test(new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]));
          h_test->StatOverflows(kTRUE);
          h_test->Sumw2(kTRUE);
          HistsJECUp[i][k][l][n] = h_test;
          name.str("");
          name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l]<<"_"<<sysJecNames[n]<<"_Down";
          std::shared_ptr<TH1F> h_test2(new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]));
          h_test2->StatOverflows(kTRUE);
          h_test2->Sumw2(kTRUE);
          HistsJECDown[i][k][l][n] = h_test2;
          name.str("");
        }
      }
    }
  }


  int reweightSizeQscalePDF = 150;
  int reweightSizePS = 14;
  Dim4 HistsSysReweightsQscalePDF(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(reweightSizeQscalePDF))));
  Dim4 HistsSysReweightsPS(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(reweightSizePS))));
  sumOfWeights SW;
  sumOfWeightsSignal SWS;
  std::vector<float> SLW;
  std::vector<float> SGW;

  if (fname.Contains("TTTo2L2Nu")){
     SLW = SW.LHEWeight(year);
     SGW = SW.GenWeight(year);
  }
  if (fname.Contains("LFV")){
     SLW = SWS.LHEWeightSignal(sname);
     SGW = SWS.GenWeightSignal(sname);
  }

  if (fname.Contains("TTTo2L2Nu") || fname.Contains("LFV")){
    for (int i=0;i<channels.size();++i){
      for (int k=0;k<regions.size();++k){
        for (int l=0;l<vars.size();++l){
          for (int n=0;n<reweightSizeQscalePDF;++n){
            name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l]<<"_QscalePDF_"<<n;
            std::shared_ptr<TH1F> h_test(new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]));
            h_test->StatOverflows(kTRUE);
            h_test->Sumw2(kTRUE);
            HistsSysReweightsQscalePDF[i][k][l][n] = h_test;
            name.str("");
          }
          for (int n=0;n<reweightSizePS;++n){
            name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l]<<"_PS_"<<n;
            std::shared_ptr<TH1F> h_test(new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]));
            h_test->StatOverflows(kTRUE);
            h_test->Sumw2(kTRUE);
            HistsSysReweightsPS[i][k][l][n] = h_test;
            name.str("");
          }
        }
      }
    }
  }


////Get scale factor and weight histograms
  TH2F  sf_Ele_Reco_H;
  TH2F  sf_Ele_ID_H;
  TH2F  sf_Mu_ID_H;
  TH2F  sf_Mu_ISO_H;
  TH2F  sf_triggeree_H;
  TH2F  sf_triggeremu_H;
  TH2F  sf_triggermumu_H;
  TH2F  btagEff_b_H;
  TH2F  btagEff_c_H;
  TH2F  btagEff_udsg_H;
  PU wPU;
  std::string rochesterFile;
  std::string btagFile;
  BTagCalibrationReader reader(BTagEntry::OP_MEDIUM, "central", {"up", "down"});
  GEScaleSyst *GE = new GEScaleSyst();

  RoccoR  rc;
  if(year == "2016")    rochesterFile = "/user/rgoldouz/NewAnalysis2020/Analysis/input/RoccoR2016.txt";
  if(year == "2017")    rochesterFile = "/user/rgoldouz/NewAnalysis2020/Analysis/input/RoccoR2017.txt";
  if(year == "2018")    rochesterFile = "/user/rgoldouz/NewAnalysis2020/Analysis/input/RoccoR2018.txt";
  rc.init(rochesterFile);

  if(data == "mc"){
    TFile *f_btagEff_Map = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/btagEff.root");
    if(year == "2016"){
      btagEff_b_H = *(TH2F*)f_btagEff_Map->Get("2016_h2_BTaggingEff_b");
      btagEff_c_H = *(TH2F*)f_btagEff_Map->Get("2016_h2_BTaggingEff_c");
      btagEff_udsg_H = *(TH2F*)f_btagEff_Map->Get("2016_h2_BTaggingEff_udsg");
    }
    if(year == "2017"){
      btagEff_b_H = *(TH2F*)f_btagEff_Map->Get("2017_h2_BTaggingEff_b");
      btagEff_c_H = *(TH2F*)f_btagEff_Map->Get("2017_h2_BTaggingEff_c");
      btagEff_udsg_H = *(TH2F*)f_btagEff_Map->Get("2017_h2_BTaggingEff_udsg");
    }
    if(year == "2018"){
      btagEff_b_H = *(TH2F*)f_btagEff_Map->Get("2018_h2_BTaggingEff_b");
      btagEff_c_H = *(TH2F*)f_btagEff_Map->Get("2018_h2_BTaggingEff_c");
      btagEff_udsg_H = *(TH2F*)f_btagEff_Map->Get("2018_h2_BTaggingEff_udsg");
    }

    if(year == "2016")    btagFile = "/user/rgoldouz/NewAnalysis2020/Analysis/input/DeepCSV_2016LegacySF_WP_V1.csv";
    if(year == "2017")    btagFile = "/user/rgoldouz/NewAnalysis2020/Analysis/input/DeepCSV_94XSF_WP_V4_B_F.csv";
    if(year == "2018")    btagFile = "/user/rgoldouz/NewAnalysis2020/Analysis/input/DeepCSV_102XSF_WP_V1.csv";

    BTagCalibration calib("DeepCSV",btagFile);
    reader.load(calib,BTagEntry::FLAV_B,"mujets"); 
    reader.load(calib,BTagEntry::FLAV_C,"mujets");
    reader.load(calib,BTagEntry::FLAV_UDSG,"incl");

    if(year == "2016"){
      TFile *f_Ele_Reco_Map = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root");
      sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

      TFile *f_Ele_ID_Map = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/2016LegacyReReco_ElectronTight_Fall17V2.root");
      sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

      TFile *f_Mu_ID_Map_1 = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/2016_RunBCDEF_SF_ID.root");
      TH2F *sf_Mu_ID_H_1 = (TH2F*)f_Mu_ID_Map_1->Get("NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt");
      TFile *f_Mu_ID_Map_2 = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/2016_RunGH_SF_ID.root");
      TH2F *sf_Mu_ID_H_2 = (TH2F*)f_Mu_ID_Map_2->Get("NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt");

      sf_Mu_ID_H_1->Scale(0.55);
      sf_Mu_ID_H_2->Scale(0.45);
      sf_Mu_ID_H_1->Add(sf_Mu_ID_H_2);
      sf_Mu_ID_H = *sf_Mu_ID_H_1;

      TFile *f_Mu_ISO_Map_1 = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/2016_RunBCDEF_SF_ISO.root");
      TH2F *sf_Mu_ISO_H_1 = (TH2F*)f_Mu_ISO_Map_1->Get("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_eta_pair_newTuneP_probe_pt");
      TFile *f_Mu_ISO_Map_2 = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/2016_RunGH_SF_ISO.root");
      TH2F *sf_Mu_ISO_H_2 = (TH2F*)f_Mu_ISO_Map_2->Get("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_eta_pair_newTuneP_probe_pt");
      sf_Mu_ISO_H_1->Scale(0.55);
      sf_Mu_ISO_H_2->Scale(0.45);
      sf_Mu_ISO_H_1->Add(sf_Mu_ISO_H_2);
      sf_Mu_ISO_H = *sf_Mu_ISO_H_1;

      TFile *f_triggeree = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/TriggerSF_ee2016_pt.root");
      sf_triggeree_H = *(TH2F*)f_triggeree->Get("h_lep1Pt_lep2Pt_Step6");
      TFile *f_triggeremu = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/TriggerSF_emu2016_Pt.root");
      sf_triggeremu_H = *(TH2F*)f_triggeremu->Get("SF_lep1Pt_lep2Pt");
      TFile *f_triggermumu = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/TriggerSF_mumu2016_pt.root");
      sf_triggermumu_H = *(TH2F*)f_triggermumu->Get("h_lep1Pt_lep2Pt_Step9");

      f_Ele_Reco_Map->Close();
      f_Ele_ID_Map->Close();
      f_Mu_ID_Map_1->Close();
      f_Mu_ID_Map_2->Close();
      f_Mu_ISO_Map_1->Close();
      f_Mu_ISO_Map_2->Close();
      f_triggeree->Close();
      f_triggeremu->Close();
      f_triggermumu->Close();
    }
    if(year == "2017"){
      TFile *f_Ele_Reco_Map = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root");
      sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

      TFile *f_Ele_ID_Map = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/2017_ElectronTight.root");
      sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

      TFile *f_Mu_ID_Map = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/2017_RunBCDEF_SF_ID_syst.root");
      sf_Mu_ID_H = *(TH2F*)f_Mu_ID_Map->Get("NUM_HighPtID_DEN_genTracks_pair_newTuneP_probe_pt_abseta");

      TFile *f_Mu_ISO_Map = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/2017_RunBCDEF_SF_ISO_syst.root");
      sf_Mu_ISO_H = *(TH2F*)f_Mu_ISO_Map->Get("NUM_LooseRelTkIso_DEN_TrkHighPtID_pair_newTuneP_probe_pt_abseta");

      TFile *f_triggeree = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/TriggerSF_ee2017_pt.root");
      sf_triggeree_H = *(TH2F*)f_triggeree->Get("h_lep1Pt_lep2Pt_Step6");
      TFile *f_triggeremu = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/TriggerSF_emu2017_Pt.root");
      sf_triggeremu_H = *(TH2F*)f_triggeremu->Get("SF_lep1Pt_lep2Pt");
      TFile *f_triggermumu = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/TriggerSF_mumu2017_pt.root");
      sf_triggermumu_H = *(TH2F*)f_triggermumu->Get("h_lep1Pt_lep2Pt_Step9");

      f_Ele_Reco_Map->Close();
      f_Ele_ID_Map->Close();
      f_Mu_ID_Map->Close();
      f_Mu_ISO_Map->Close();
      f_triggeree->Close();
      f_triggeremu->Close();
      f_triggermumu->Close();
    }
    if(year == "2018"){
      TFile *f_Ele_Reco_Map = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/egammaEffi.txt_EGM2D_updatedAll.root");
      sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

      TFile *f_Ele_ID_Map = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/2018_ElectronTight.root");
      sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

      TFile *f_Mu_ID_Map = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/2018_RunABCD_SF_ID.root");
      sf_Mu_ID_H = *(TH2F*)f_Mu_ID_Map->Get("NUM_HighPtID_DEN_TrackerMuons_pair_newTuneP_probe_pt_abseta");

      TFile *f_Mu_ISO_Map = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/2018_RunABCD_SF_ISO.root");
      sf_Mu_ISO_H = *(TH2F*)f_Mu_ISO_Map->Get("NUM_LooseRelTkIso_DEN_TrkHighPtID_pair_newTuneP_probe_pt_abseta");

      TFile *f_triggeree = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/TriggerSF_ee2018_pt.root");
      sf_triggeree_H = *(TH2F*)f_triggeree->Get("h_lep1Pt_lep2Pt_Step6");
      TFile *f_triggeremu = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/TriggerSF_emu2018_Pt.root");
      sf_triggeremu_H = *(TH2F*)f_triggeremu->Get("SF_lep1Pt_lep2Pt");
      TFile *f_triggermumu = new TFile("/user/rgoldouz/NewAnalysis2020/Analysis/input/TriggerSF_mumu2018_pt.root");
      sf_triggermumu_H = *(TH2F*)f_triggermumu->Get("h_lep1Pt_lep2Pt_Step9");

      f_Ele_Reco_Map->Close();
      f_Ele_ID_Map->Close();
      f_Mu_ID_Map->Close();
      f_Mu_ISO_Map->Close();
      f_triggeree->Close();
      f_triggeremu->Close();
      f_triggermumu->Close();
    }
  }

  TFile file_out (fname,"RECREATE");
  TTree tree_out("analysis","main analysis") ;

//    cout<<"ev_event"<<"   "<<"sf_Ele_Reco"<<"   "<<"sf_Ele_ID"<<"      "<<"sf_Mu_ID"<<"   "<<"sf_Mu_ISO"<<"   "<<"sf_trigger"<<"   "<<"PU weight"<<endl;
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

  TLorentzVector wp, wm, b, ab, top, atop;
  std::vector<float> nominalWeights;
  nominalWeights.assign(sys.size(), 1);
  std::vector<float> sysUpWeights;
  sysUpWeights.assign(sys.size(), 1);
  std::vector<float> sysDownWeights;
  sysDownWeights.assign(sys.size(), 1);
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
  float elePt;
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
//  for (Long64_t jentry=0; jentry<10;jentry++) {
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
    for (int n=0;n<sys.size();++n){
      nominalWeights[n] =1;
      sysUpWeights[n] =1;
      sysDownWeights[n] =1;
    }
//MET filters

    if(data == "mc"){
      if(year == "2016" || year == "2018" ){
        if ( trig_Flag_goodVertices_accept==1  &&  trig_Flag_globalSuperTightHalo2016Filter_accept==1 && trig_Flag_HBHENoiseFilter_accept==1 &&  trig_Flag_HBHENoiseIsoFilter_accept==1 && trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept==1 && trig_Flag_BadPFMuonFilter_accept==1)
        metFilterPass = true;
        }
        if(year == "2017"){
        if ( trig_Flag_goodVertices_accept==1  &&  trig_Flag_globalSuperTightHalo2016Filter_accept==1 && trig_Flag_HBHENoiseFilter_accept==1 &&  trig_Flag_HBHENoiseIsoFilter_accept==1 && trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept==1 && trig_Flag_BadPFMuonFilter_accept==1 && trig_Flag_ecalBadCalibReduced ==1)
        metFilterPass = true;
        }
    }

    if(data == "data"){
      if(year == "2016" || year == "2018"){
        if ( trig_Flag_goodVertices_accept==1  &&  trig_Flag_globalSuperTightHalo2016Filter_accept==1 && trig_Flag_HBHENoiseFilter_accept==1 &&  trig_Flag_HBHENoiseIsoFilter_accept==1 && trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept==1 && trig_Flag_BadPFMuonFilter_accept==1 &&  trig_Flag_eeBadScFilter_accept==1)
        metFilterPass = true;
        }
        if(year == "2017"){
        if ( trig_Flag_goodVertices_accept==1  &&  trig_Flag_globalSuperTightHalo2016Filter_accept==1 && trig_Flag_HBHENoiseFilter_accept==1 &&  trig_Flag_HBHENoiseIsoFilter_accept==1 && trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept==1 && trig_Flag_BadPFMuonFilter_accept==1 && trig_Flag_ecalBadCalibReduced ==1)
        metFilterPass = true;
        }
    }

//trigger
////MC
      if(data == "mc" && year == "2016"){
        if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Ele27_WPTight_Gsf_accept ) triggerPassEE =true;
        if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele27_WPTight_Gsf_accept || trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept) triggerPassEMu =true;
        if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept || trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept) triggerPassMuMu =true;
      }

      if(data == "mc" && year == "2017"){
        if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele35_WPTight_Gsf_accept) triggerPassEE =true;
        if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele35_WPTight_Gsf_accept || trig_HLT_IsoMu27_accept) triggerPassEMu =true;
        if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept || trig_HLT_IsoMu27_accept) triggerPassMuMu =true;
      }

      if(data == "mc" && year == "2018"){
        if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept) triggerPassEE =true;
        if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept || trig_HLT_IsoMu24_accept) triggerPassEMu =true;
        if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept || trig_HLT_IsoMu24_accept) triggerPassMuMu =true;
      } 

////DATA
    if(data == "data"){
      if(year == "2016"){
        if(run == "H"){
          if(dataset=="MuonEG"){
            if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) triggerPassEMu =true;
          }
          if(dataset=="SingleElectron"){
            if(!(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEE =true;
            if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEMu =true;
          }
          if(dataset=="SingleMuon"){
            if(!(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassMuMu =true;
            if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Ele27_WPTight_Gsf_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassEMu =true; 
          }
          if(dataset=="DoubleEG"){
            if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) triggerPassEE =true;
          }
          if(dataset=="DoubleMu"){
            if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept) triggerPassMuMu =true;
          }
        }
        if(run != "H"){
          if(dataset=="MuonEG"){
            if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) triggerPassEMu =true;
          }
          if(dataset=="SingleElectron"){
            if(!(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEE =true;
            if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEMu =true;
          }
          if(dataset=="SingleMuon"){
            if(!(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassMuMu =true;
            if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele27_WPTight_Gsf_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassEMu =true;
          }
          if(dataset=="DoubleEG"){
            if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) triggerPassEE =true;
          }
          if(dataset=="DoubleMu"){
            if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept) triggerPassMuMu =true;
          }
        }
      }
      if(year == "2017"){
        if(dataset=="MuonEG"){
          if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept ) triggerPassEMu =true;
        }
        if(dataset=="SingleElectron"){
          if(!trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept && trig_HLT_Ele35_WPTight_Gsf_accept) triggerPassEE =true;
          if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) && trig_HLT_Ele35_WPTight_Gsf_accept) triggerPassEMu =true;
        }
        if(dataset=="SingleMuon"){
           if(!trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept && trig_HLT_IsoMu27_accept) triggerPassMuMu =true;
          if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele35_WPTight_Gsf_accept) && trig_HLT_IsoMu27_accept) triggerPassEMu =true;
        }
        if(dataset=="DoubleEG"){
          if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept) triggerPassEE =true;
        }
        if(dataset=="DoubleMu"){
          if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept)triggerPassMuMu =true;
        }
      }
      if(year == "2018"){
        if(dataset=="MuonEG"){
          if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) triggerPassEMu =true;
        }
        if(dataset=="EGamma"){
          if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept) triggerPassEE =true;
          if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) && trig_HLT_Ele32_WPTight_Gsf_accept) triggerPassEMu =true;
      }
      if(dataset=="SingleMuon"){
        if(!trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept && trig_HLT_IsoMu24_accept) triggerPassMuMu =true;
        if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept) && trig_HLT_IsoMu24_accept) triggerPassEMu =true;
      }
      if(dataset=="DoubleMu"){
        if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept) triggerPassMuMu =true;
      }
    }
  }
 


//cout<<ev_event<<"  "<<trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept <<"  "<< trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept <<"  "<< trig_HLT_Ele27_WPTight_Gsf_accept  <<"  "<<trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept <<"  "<< trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept <<"  "<< trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept<<"  "<<trig_HLT_IsoMu24_accept<<"  "<<trig_HLT_IsoTkMu24_accept<<endl;
//cout<<"trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept "<< "trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept "<< "trig_HLT_Ele27_WPTight_Gsf_accept " <<"trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept "<< "trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept "<< "trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept "  <<"trig_HLT_IsoMu24_accept "<<"trig_HLT_IsoTkMu24_accept"<<endl;
    if(!(triggerPassEE || triggerPassEMu || triggerPassMuMu)) continue;
    if(!metFilterPass) continue;

// lepton selection
  selectedLeptons = new std::vector<lepton_candidate*>();
  selectedLeptonsMuScaleUp = new std::vector<lepton_candidate*>();
  selectedLeptonsMuScaleDown = new std::vector<lepton_candidate*>();
  selectedLeptonsEleScaleUp = new std::vector<lepton_candidate*>();
  selectedLeptonsEleScaleDown = new std::vector<lepton_candidate*>();
  selectedLeptonsMuResUp = new std::vector<lepton_candidate*>();
  selectedLeptonsMuResDown = new std::vector<lepton_candidate*>();

// electron
    for (int l=0;l<gsf_pt->size();l++){
      elePt = (*gsf_ecalTrkEnergyPostCorr)[l]*sin(2.*atan(exp(-1.*(*gsf_eta)[l]))) ;
      if(abs((*gsf_eta)[l]) > 2.4 || (abs((*gsf_sc_eta)[l])> 1.4442 && (abs((*gsf_sc_eta)[l])< 1.566))) continue;
      if(!(*gsf_VID_cutBasedElectronID_Fall17_94X_V2_tight)[l]) continue;
      if (data == "mc"){
        if((elePt + 0.004*elePt) > 20 && abs((*gsf_sc_eta)[l])< 1.4442) selectedLeptonsEleScaleUp->push_back(new lepton_candidate(elePt+ 0.004*elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
        if((elePt + 0.008*elePt) > 20 && abs((*gsf_sc_eta)[l])> 1.4442) selectedLeptonsEleScaleUp->push_back(new lepton_candidate(elePt+ 0.008*elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
        if((elePt - 0.004*elePt) > 20 && abs((*gsf_sc_eta)[l])< 1.4442) selectedLeptonsEleScaleDown->push_back(new lepton_candidate(elePt- 0.004*elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
        if((elePt - 0.008*elePt) > 20 && abs((*gsf_sc_eta)[l])> 1.4442) selectedLeptonsEleScaleDown->push_back(new lepton_candidate(elePt- 0.008*elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
      }
      if(elePt <20) continue;
      selectedLeptons->push_back(new lepton_candidate(elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
      if (data == "mc"){
        selectedLeptonsMuScaleUp->push_back(new lepton_candidate(elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
        selectedLeptonsMuScaleDown->push_back(new lepton_candidate(elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
        selectedLeptonsMuResUp->push_back(new lepton_candidate(elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
        selectedLeptonsMuResDown->push_back(new lepton_candidate(elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
        sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"central",false, true);
        nominalWeights[0] = nominalWeights[0] * scale_factor(&sf_Ele_Reco_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"central",false, true);
        sysUpWeights[0] = sysUpWeights[0] * scale_factor(&sf_Ele_Reco_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"up",false, true);
        sysDownWeights[0] = sysDownWeights[0] * scale_factor(&sf_Ele_Reco_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"down",false, true);

        sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"central",false, true);
        nominalWeights[1] = nominalWeights[1] * scale_factor(&sf_Ele_ID_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"central",false, true);
        sysUpWeights[1] = sysUpWeights[1] * scale_factor(&sf_Ele_ID_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"up",false, true);
        sysDownWeights[1] = sysDownWeights[1] * scale_factor(&sf_Ele_ID_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"down",false, true);
      }
    }

// Muon selection
    for (int l=0;l<mu_ibt_pt->size();l++){
      if(abs((*mu_ibt_eta)[l]) > 2.4) continue;
      if(!(*mu_isHighPtMuon)[l]) continue;
      if(!(*mu_TkIsoLoose)[l]) continue;
//is normal muon
      if((*mu_ibt_pt)[l] == (*mu_it_pt)[l]){
        if(data == "data") {
          muPtSFRochester = rc.kScaleDT((*mu_ibt_charge)[l], (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l], 0, 0);
          if(muPtSFRochester * (*mu_ibt_pt)[l] > 20) selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));
        }
        if (data == "mc"){
          if ((*mu_mc_index)[l]!=-1 && abs((*mc_pdgId)[(*mu_mc_index)[l]]) == 13) muPtSFRochester = rc.kSpreadMC((*mu_ibt_charge)[l], (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l], (*mc_pt)[(*mu_mc_index)[l]],0, 0);
          if ((*mu_mc_index)[l]<0) muPtSFRochester = rc.kSmearMC((*mu_ibt_charge)[l], (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l], (*mu_trackerLayersWithMeasurement)[l] , gRandom->Rndm(),0, 0);
          if(muPtSFRochester * (*mu_ibt_pt)[l] > 20) selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));
          if(muPtSFRochester * (*mu_ibt_pt)[l] > 20) selectedLeptonsEleScaleUp->push_back(new lepton_candidate(muPtSFRochester * (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));
          if(muPtSFRochester * (*mu_ibt_pt)[l] > 20) selectedLeptonsEleScaleDown->push_back(new lepton_candidate(muPtSFRochester * (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));
          if(muPtSFRochester * (*mu_ibt_pt)[l] > 20) selectedLeptonsMuResUp->push_back(new lepton_candidate(muPtSFRochester * (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10)); 
          if(muPtSFRochester * (*mu_ibt_pt)[l] > 20) selectedLeptonsMuResDown->push_back(new lepton_candidate(muPtSFRochester * (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));
          if(muPtSFRochester * (*mu_ibt_pt)[l] > 20) selectedLeptonsMuScaleUp->push_back(new lepton_candidate(muPtSFRochester * (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));
          if(muPtSFRochester * (*mu_ibt_pt)[l] > 20) selectedLeptonsMuScaleDown->push_back(new lepton_candidate(muPtSFRochester * (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));
        }
      }
      if((*mu_ibt_pt)[l] != (*mu_it_pt)[l]){
        pt_res = 0.02;
        if(abs((*mu_ibt_eta)[l]) < 1.2) pt_res = 0.01;
        if (data == "data" && (*mu_ibt_pt)[l]>20) selectedLeptons->push_back(new lepton_candidate((*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));
        if(data == "mc"){
          if (GE->GEScaleCorrPt(1600, (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l]) > 20) selectedLeptons->push_back(new lepton_candidate(GE->GEScaleCorrPt(1600, (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l]),(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));
          if (GE->GEScaleCorrPt(1600, (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l]) > 20) selectedLeptonsEleScaleUp->push_back(new lepton_candidate(GE->GEScaleCorrPt(1600, (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l]),(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));
          if (GE->GEScaleCorrPt(1600, (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l]) > 20) selectedLeptonsEleScaleDown->push_back(new lepton_candidate(GE->GEScaleCorrPt(1600, (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l]),(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));

          if (GE->GEScaleCorrPt(1600, (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l])*(Tr.Gaus(1, pt_res)) > 20) selectedLeptonsMuResUp->push_back(new lepton_candidate((Tr.Gaus(1, pt_res))*GE->GEScaleCorrPt(1600, (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l]),(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));
          if (GE->GEScaleCorrPt(1600, (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l])*(Tr.Gaus(1, pt_res)) > 20) selectedLeptonsMuResDown->push_back(new lepton_candidate((Tr.Gaus(1, pt_res))*GE->GEScaleCorrPt(1600, (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l]),(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));

          if (GE->GEScaleCorrPt(1601, (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l]) > 20) selectedLeptonsMuScaleUp->push_back(new lepton_candidate(GE->GEScaleCorrPt(1601, (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l]),(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));
          if (GE->GEScaleCorrPt(1602, (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l]) > 20) selectedLeptonsMuScaleDown->push_back(new lepton_candidate(GE->GEScaleCorrPt(1602, (*mu_ibt_pt)[l],(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l]),(*mu_ibt_eta)[l],(*mu_ibt_phi)[l],(*mu_ibt_charge)[l],l,10));
        }
      }
      if((*mu_ibt_pt)[l] <20) continue;
      if (data == "mc" && year == "2016") {
        sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_ibt_eta)[l], (*mu_ibt_pt)[l],"central",false, true);
        nominalWeights[2] = nominalWeights[2] * scale_factor(&sf_Mu_ID_H, (*mu_ibt_eta)[l], (*mu_ibt_pt)[l],"central",false, true);
        sysUpWeights[2] = sysUpWeights[2] * scale_factor(&sf_Mu_ID_H, (*mu_ibt_eta)[l], (*mu_ibt_pt)[l],"up",false, true);
        sysDownWeights[2] = sysDownWeights[2] * scale_factor(&sf_Mu_ID_H, (*mu_ibt_eta)[l], (*mu_ibt_pt)[l],"down",false, true);

        sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_ibt_eta)[l], (*mu_ibt_pt)[l],"central",false, true);
        nominalWeights[3] = nominalWeights[3] * scale_factor(&sf_Mu_ISO_H, (*mu_ibt_eta)[l], (*mu_ibt_pt)[l],"central",false, true);
        sysUpWeights[3] = sysUpWeights[3] * scale_factor(&sf_Mu_ISO_H, (*mu_ibt_eta)[l], (*mu_ibt_pt)[l],"up",false, true);
        sysDownWeights[3] = sysDownWeights[3] * scale_factor(&sf_Mu_ISO_H, (*mu_ibt_eta)[l], (*mu_ibt_pt)[l],"down",false, true);
      }
      if (data == "mc" && year != "2016") {
        sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_ibt_pt)[l], abs((*mu_ibt_eta)[l]),"central",false, true);
        nominalWeights[2] = nominalWeights[2] * scale_factor(&sf_Mu_ID_H, (*mu_ibt_pt)[l], abs((*mu_ibt_eta)[l]),"central",false, true);
        sysUpWeights[2] = sysUpWeights[2] * scale_factor(&sf_Mu_ID_H, (*mu_ibt_pt)[l], abs((*mu_ibt_eta)[l]),"up",false, true);
        sysDownWeights[2] = sysDownWeights[2] * scale_factor(&sf_Mu_ID_H, (*mu_ibt_pt)[l], abs((*mu_ibt_eta)[l]),"down",false, true);

        sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_ibt_pt)[l], abs((*mu_ibt_eta)[l]),"central",false, true);
        nominalWeights[3] = nominalWeights[3] * scale_factor(&sf_Mu_ISO_H, (*mu_ibt_pt)[l], abs((*mu_ibt_eta)[l]),"central",false, true);
        sysUpWeights[3] = sysUpWeights[3] * scale_factor(&sf_Mu_ISO_H, (*mu_ibt_pt)[l], abs((*mu_ibt_eta)[l]),"up",false, true);
        sysDownWeights[3] = sysDownWeights[3] * scale_factor(&sf_Mu_ISO_H, (*mu_ibt_pt)[l], abs((*mu_ibt_eta)[l]),"down",false, true);
      }
    }
    sort(selectedLeptons->begin(), selectedLeptons->end(), ComparePtLep);
    sort(selectedLeptonsMuScaleUp->begin(), selectedLeptonsMuScaleUp->end(), ComparePtLep);
    sort(selectedLeptonsMuScaleDown->begin(), selectedLeptonsMuScaleDown->end(), ComparePtLep);
    sort(selectedLeptonsMuResDown->begin(), selectedLeptonsMuResDown->end(), ComparePtLep);
    sort(selectedLeptonsMuResUp->begin(), selectedLeptonsMuResUp->end(), ComparePtLep);
    sort(selectedLeptonsEleScaleUp->begin(), selectedLeptonsEleScaleUp->end(), ComparePtLep);
    sort(selectedLeptonsEleScaleDown->begin(), selectedLeptonsEleScaleDown->end(), ComparePtLep);
// dilepton selection
//cout<<ev_event<<"  "<<triggerPass<<"  "<<metFilterPass<<"  "<<selectedLeptons->size()<<endl;
    if(selectedLeptons->size()!=2 ||
      ((*selectedLeptons)[0]->pt_ <25) ||
      ((*selectedLeptons)[0]->charge_ * (*selectedLeptons)[1]->charge_ == 1) ||
//      ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ != 11 && ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()<106 && ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()>76) ||
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
    for (int l=0;l<jet_pt->size();l++){
      if(year == "2016" && !(*jet_isJetIDTightLepVeto_2016)[l]) continue;
      if(year == "2017" && !(*jet_isJetIDLepVeto_2017)[l]) continue;
      if(year == "2018" && !(*jet_isJetIDLepVeto_2018)[l]) continue;
      jetlepfail = false;
      for (int i=0;i<selectedLeptons->size();i++){
        if(deltaR((*selectedLeptons)[i]->eta_,(*selectedLeptons)[i]->phi_,(*jet_eta)[l],(*jet_phi)[l]) < 0.4 ) jetlepfail=true;
      }
      if(jetlepfail) continue; 
      if(data == "mc" && abs((*jet_eta)[l]) < 2.4){
        if((*jet_Smeared_pt)[l] >30) selectedJets->push_back(new jet_candidate((*jet_Smeared_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l], year,(*jet_partonFlavour)[l]));
        if((*jet_SmearedJetResUp_pt)[l] >30) selectedJetsJerUp->push_back(new jet_candidate((*jet_SmearedJetResUp_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l], year,(*jet_partonFlavour)[l]));
        if((*jet_SmearedJetResDown_pt)[l] >30) selectedJetsJerDown->push_back(new jet_candidate((*jet_SmearedJetResDown_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l], year,(*jet_partonFlavour)[l]));
        if((*jet_SmearedJetEnUp_pt)[l] >30) selectedJetsJesUp->push_back(new jet_candidate((*jet_SmearedJetEnUp_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l], year,(*jet_partonFlavour)[l]));
        if((*jet_SmearedJetEnDown_pt)[l] >30) selectedJetsJesDown->push_back(new jet_candidate((*jet_SmearedJetEnDown_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l], year,(*jet_partonFlavour)[l]));
      }
      if(data == "data" && (*jet_pt)[l] >30 && abs((*jet_eta)[l]) < 2.4){
        selectedJets->push_back(new jet_candidate((*jet_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l],year,0));
        selectedJetsJerUp->push_back(new jet_candidate((*jet_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l],year,0));
        selectedJetsJerDown->push_back(new jet_candidate((*jet_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l],year,0));
        selectedJetsJesUp->push_back(new jet_candidate((*jet_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l],year,0));
        selectedJetsJesDown->push_back(new jet_candidate((*jet_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l],year,0));
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

    JECsysUp = new std::vector<std::vector<jet_candidate*>>();
    JECsysDown = new std::vector<std::vector<jet_candidate*>>();
    JECsysMETUp = new std::vector<float>();
    JECsysMETDown = new std::vector<float>();
    JECsysMVAUp = new std::vector<float>();
    JECsysMVADown = new std::vector<float>();

    for (int n=0;n<sysJecNames.size();++n){
      JECJetsUp= new std::vector<jet_candidate*>();
      JECJetsDown= new std::vector<jet_candidate*>();
      double sup = 0;
      double sdw = 0;
      JECMETUpx =   0;
      JECMETUpy =   0;
      JECMETDownx = 0;
      JECMETDowny = 0;
      for (int l=0;l<jet_pt->size();l++){
        if(data == "data" || abs((*jet_eta)[l]) > 2.4) continue;
        if(year == "2016" && !(*jet_isJetIDTightLepVeto_2016)[l]) continue;
        if(year == "2017" && !(*jet_isJetIDLepVeto_2017)[l]) continue;
        if(year == "2018" && !(*jet_isJetIDLepVeto_2018)[l]) continue;
        jetlepfail = false;
        for (int i=0;i<selectedLeptons->size();i++){
          if(deltaR((*selectedLeptons)[i]->eta_,(*selectedLeptons)[i]->phi_,(*jet_eta)[l],(*jet_phi)[l]) < 0.4 ) jetlepfail=true;
        }
        if(jetlepfail) continue;
        JetCorrectionUncertainty *unc = vsrc[n];
        unc->setJetPt((*jet_Smeared_pt)[l]);
        unc->setJetEta((*jet_eta)[l]);
        sup = unc->getUncertainty(true); // up variation
        if ((1+sup)*(*jet_Smeared_pt)[l]>30) {
          JECJetsUp->push_back(new jet_candidate((1+sup)*(*jet_Smeared_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l], year,(*jet_partonFlavour)[l]));
          JECMETUpx = JECMETUpx + (sup * (*jet_Smeared_pt)[l] * cos((*jet_phi)[l]));
          JECMETUpy = JECMETUpy + (sup * (*jet_Smeared_pt)[l] * sin((*jet_phi)[l])); 
        }
        unc->setJetPt((*jet_Smeared_pt)[l]);
        unc->setJetEta((*jet_eta)[l]);
        sdw = unc->getUncertainty(false); // down variation
        if ((1-sdw)*(*jet_Smeared_pt)[l]>30){
          JECJetsDown->push_back(new jet_candidate((1-sdw)*(*jet_Smeared_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l], year,(*jet_partonFlavour)[l]));
          JECMETDownx = JECMETDownx - (sdw * (*jet_Smeared_pt)[l] * cos((*jet_phi)[l]));
          JECMETDowny = JECMETDowny - (sdw * (*jet_Smeared_pt)[l] * sin((*jet_phi)[l]));
        }
      }
      sort(JECJetsUp->begin(), JECJetsUp->end(), ComparePtJet);
      sort(JECJetsDown->begin(), JECJetsDown->end(), ComparePtJet);
      JECsysUp->push_back(*JECJetsUp);
      JECsysDown->push_back(*JECJetsDown);
      JECsysMETUp->push_back(sqrt(pow(MET_FinalCollection_Pt * cos(MET_FinalCollection_phi) - JECMETUpx,2)+pow(MET_FinalCollection_Pt * sin(MET_FinalCollection_phi) - JECMETUpy,2)));
      JECsysMETDown->push_back(sqrt(pow(MET_FinalCollection_Pt * cos(MET_FinalCollection_phi) - JECMETDownx,2)+pow(MET_FinalCollection_Pt * sin(MET_FinalCollection_phi) - JECMETDowny,2)));
//      cout<<n<<"-nominal="<<MET_FinalCollection_Pt<<"["<<(*JECsysMETDown)[n]<<","<<(*JECsysMETUp)[n]<<"]"<<endl;
   }

// MVA output
// Nominal output
    leading_pt = (*selectedLeptons)[0]->pt_;
    jet_leading_pt = 0;
    if (selectedJets->size()>0) jet_leading_pt = (*selectedJets)[0]->pt_;
    deltaR_ll = deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_);
    MET = MET_FinalCollection_Pt;
    n_jet = selectedJets->size();
    MVAoutput = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");

//MVA systematic output
    if (data == "mc"){
//jes Up, Down
      jet_leading_pt = 0;
      if (selectedJetsJesUp->size()>0) jet_leading_pt = (*selectedJetsJesUp)[0]->pt_;
      MET = MET_T1SmearJetEnUp_Pt;
      n_jet = selectedJetsJesUp->size();
      MVAoutputJesUp = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");
      jet_leading_pt = 0;
      if (selectedJetsJesDown->size()>0) jet_leading_pt = (*selectedJetsJesDown)[0]->pt_;
      MET = MET_T1SmearJetEnDown_Pt;
      n_jet = selectedJetsJesDown->size();
      MVAoutputJesDown = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");
//Jer Up, Down
      jet_leading_pt = 0;
      if (selectedJetsJerUp->size()>0) jet_leading_pt = (*selectedJetsJerUp)[0]->pt_;
      MET = MET_T1SmearJetResUp_Pt;
      n_jet = selectedJetsJerUp->size();
      MVAoutputJerUp = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");
      jet_leading_pt = 0;
      if (selectedJetsJerDown->size()>0) jet_leading_pt = (*selectedJetsJerDown)[0]->pt_;
      MET = MET_T1SmearJetResDown_Pt;
      n_jet = selectedJetsJerDown->size();
      MVAoutputJerDown = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");
//lepton scale/resolution Up, Down
      if (selectedJets->size()>0) jet_leading_pt = (*selectedJets)[0]->pt_;
      MET = MET_FinalCollection_Pt;
      n_jet = selectedJets->size();
      if (selectedLeptonsMuScaleUp->size()>0) leading_pt = (*selectedLeptonsMuScaleUp)[0]->pt_;
      MVAoutputMuScaleUp = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");
      if (selectedLeptonsMuScaleDown->size()>0) leading_pt = (*selectedLeptonsMuScaleDown)[0]->pt_;
      MVAoutputMuScaleDown = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");
      if (selectedLeptonsEleScaleUp->size()>0) leading_pt = (*selectedLeptonsEleScaleUp)[0]->pt_;
      MVAoutputEleScaleUp = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");
      if (selectedLeptonsEleScaleDown->size()>0) leading_pt = (*selectedLeptonsEleScaleDown)[0]->pt_;
      MVAoutputEleScaleDown = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");
      if (selectedLeptonsMuResUp->size()>0) leading_pt = (*selectedLeptonsMuResUp)[0]->pt_;
      MVAoutputMuResUp = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");
      if (selectedLeptonsMuResDown->size()>0) leading_pt = (*selectedLeptonsMuResDown)[0]->pt_;
      MVAoutputMuResDown = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");
//Uncluster MET Up, Down  
      leading_pt = (*selectedLeptons)[0]->pt_;
      MET = (*MET_Type1Unc_Pt)[10];
      MVAoutputUnclusMETUp = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");
      MET = (*MET_Type1Unc_Pt)[11];
      MVAoutputUnclusMETDown = readerMVA->EvaluateMVA( "BDT_1b_tc_BDT");
      metUnclusMETUp = (*MET_Type1Unc_Pt)[10];
      metUnclusMETDown = (*MET_Type1Unc_Pt)[11];
  }

    for (int n=0;n<sysJecNames.size();++n){
      jet_leading_pt = 0;
      if ((*JECsysUp)[n].size()>0) jet_leading_pt = (*JECsysUp)[n][0]->pt_;
      MET = (*JECsysMETUp)[n];
      n_jet = (*JECsysUp)[n].size();
      JECsysMVAUp->push_back(readerMVA->EvaluateMVA( "BDT_1b_tc_BDT"));
      jet_leading_pt = 0;
      if ((*JECsysDown)[n].size()>0) jet_leading_pt = (*JECsysDown)[n][0]->pt_;
      MET = (*JECsysMETDown)[n];
      n_jet = (*JECsysDown)[n].size();
      JECsysMVADown->push_back(readerMVA->EvaluateMVA( "BDT_1b_tc_BDT"));
  }


    JECsysNbtagUp = new std::vector<int>();
    JECsysNbtagDown = new std::vector<int>();
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

    for (int l=0;l<selectedJets->size();l++){
      if((*selectedJets)[l]->btag_) nbjet++;
      if(data == "data") continue;
      if( abs((*selectedJets)[l]->flavor_) == 5){
        nbjetGen++;
        h2_BTaggingEff_Denom_b->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
        if( (*selectedJets)[l]->btag_ ) {
          h2_BTaggingEff_Num_b->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false);
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          nominalWeights[4] = nominalWeights[4] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysUpWeights[4] = sysUpWeights[4] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("up", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysDownWeights[4] = sysDownWeights[4] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("down", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
 
          nominalWeights[5] = nominalWeights[5] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysUpWeights[5] = sysUpWeights[5] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysDownWeights[5] = sysDownWeights[5] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          nominalWeights[4] = nominalWeights[4]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysUpWeights[4] = sysUpWeights[4]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("up", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysDownWeights[4] = sysDownWeights[4]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("down", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));

          nominalWeights[5] = nominalWeights[5]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysUpWeights[5] = sysUpWeights[5]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysDownWeights[5] = sysDownWeights[5]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }  
      }
      if( abs((*selectedJets)[l]->flavor_) == 4){
        h2_BTaggingEff_Denom_c->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
        if( (*selectedJets)[l]->btag_) {
          h2_BTaggingEff_Num_c->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false);
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          nominalWeights[4] = nominalWeights[4] * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysUpWeights[4] = sysUpWeights[4] * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("up", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysDownWeights[4] = sysDownWeights[4] * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("down", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);

          nominalWeights[5] = nominalWeights[5] * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysUpWeights[5] = sysUpWeights[5] * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysDownWeights[5] = sysDownWeights[5] * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          nominalWeights[4] = nominalWeights[4]* (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysUpWeights[4] = sysUpWeights[4]* (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("up", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysDownWeights[4] = sysDownWeights[4]* (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("down", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));

          nominalWeights[5] = nominalWeights[5]* (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysUpWeights[5] = sysUpWeights[5]* (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysDownWeights[5] = sysDownWeights[5]* (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }
      }
      if( abs((*selectedJets)[l]->flavor_) != 4 && abs((*selectedJets)[l]->flavor_) != 5){
        h2_BTaggingEff_Denom_udsg->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
        if( (*selectedJets)[l]->btag_) {
          h2_BTaggingEff_Num_udsg->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false);
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          nominalWeights[4] = nominalWeights[4]* scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysUpWeights[4] = sysUpWeights[4]* scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysDownWeights[4] = sysDownWeights[4]* scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);

          nominalWeights[5] = nominalWeights[5] * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysUpWeights[5] = sysUpWeights[5] * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("up", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysDownWeights[5] = sysDownWeights[5] * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("down", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          nominalWeights[4] = nominalWeights[4]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysUpWeights[4] = sysUpWeights[4]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysDownWeights[4] = sysDownWeights[4]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));

          nominalWeights[5] = nominalWeights[5]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysUpWeights[5] = sysUpWeights[5]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("up", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysDownWeights[5] = sysDownWeights[5]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"central", true, false) * reader.eval_auto_bounds("down", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }
      }
    }

//PU reweighting
    if (data == "mc" && year == "2016") {
      weight_PU = wPU.PU_2016(mc_trueNumInteractions,"nominal");
      nominalWeights[6] = wPU.PU_2016(mc_trueNumInteractions,"nominal");
      sysUpWeights[6] = wPU.PU_2016(mc_trueNumInteractions,"up");
      sysDownWeights[6] = wPU.PU_2016(mc_trueNumInteractions,"down");
    }
    if (data == "mc" && year == "2017") {
      weight_PU = wPU.PU_2017(mc_trueNumInteractions,"nominal");
      nominalWeights[6] = wPU.PU_2017(mc_trueNumInteractions,"nominal");
      sysUpWeights[6] = wPU.PU_2017(mc_trueNumInteractions,"up");
      sysDownWeights[6] = wPU.PU_2017(mc_trueNumInteractions,"down");
    }
    if (data == "mc" && year == "2018") {
      weight_PU = wPU.PU_2018(mc_trueNumInteractions,"nominal");
      nominalWeights[6] = wPU.PU_2018(mc_trueNumInteractions,"nominal");
      sysUpWeights[6] = wPU.PU_2018(mc_trueNumInteractions,"up");
      sysDownWeights[6] = wPU.PU_2018(mc_trueNumInteractions,"down");
    }

//lumi xs weights
    if (data == "mc") weight_Lumi = (1000*xs*lumi)/Nevent;

    if (data == "mc" && (year == "2016" || year == "2017")) {
      weight_prefiring = ev_prefiringweight;
      nominalWeights[7] = ev_prefiringweight;
      sysUpWeights[7] = ev_prefiringweightup;
      sysDownWeights[7] = ev_prefiringweightdown;
    }

//trigger SF
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


    if (data == "mc" && ifTopPt) {
      for (int l=0;l<mc_status->size();l++){
        if((*mc_status)[l]<30 && (*mc_status)[l]>20 && (*mc_pdgId)[l]==24) wp.SetPtEtaPhiE((*mc_pt)[l], (*mc_eta)[l], (*mc_phi)[l], (*mc_energy)[l]) ;
        if((*mc_status)[l]<30 && (*mc_status)[l]>20 && (*mc_pdgId)[l]==-24) wm.SetPtEtaPhiE((*mc_pt)[l], (*mc_eta)[l], (*mc_phi)[l], (*mc_energy)[l]) ;
        if((*mc_status)[l]<30 && (*mc_status)[l]>20 && (*mc_pdgId)[l]==5) b.SetPtEtaPhiE((*mc_pt)[l], (*mc_eta)[l], (*mc_phi)[l], (*mc_energy)[l]) ;
        if((*mc_status)[l]<30 && (*mc_status)[l]>20 && (*mc_pdgId)[l]==-5) ab.SetPtEtaPhiE((*mc_pt)[l], (*mc_eta)[l], (*mc_phi)[l], (*mc_energy)[l]) ;
        if((*mc_status)[l]<30 && (*mc_status)[l]>20 && (*mc_pdgId)[l]==6) top.SetPtEtaPhiE((*mc_pt)[l], (*mc_eta)[l], (*mc_phi)[l], (*mc_energy)[l]) ;
        if((*mc_status)[l]<30 && (*mc_status)[l]>20 && (*mc_pdgId)[l]==-6) atop.SetPtEtaPhiE((*mc_pt)[l], (*mc_eta)[l], (*mc_phi)[l], (*mc_energy)[l]) ;
      }
    weight_topPtPowheg = sqrt(topPtPowheg((wp + b).Pt()) * topPtPowheg((wm + ab).Pt()));
    weight_topPtMGLO = sqrt(topPtMGLO((atop).Pt()) * topPtMGLO((top).Pt()));
    }

    if (fname.Contains("LFVTt")) weight_topPtPowheg = weight_topPtMGLO;
    if (data == "mc") weight_lep = sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO * sf_Trigger * weight_PU * weight_Lumi  * mc_w_sign * weight_prefiring * weight_topPtPowheg;
    if (data == "mc") weight_lepB = sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO * sf_Trigger * weight_PU * weight_Lumi * mc_w_sign *  weight_prefiring * weight_topPtPowheg * (P_bjet_data/P_bjet_mc);

if(isnan(weight_lepB) || isinf(weight_lepB)) cout<<weight_topPtMGLO<<"  "<<weight_topPtPowheg<<endl;

if(ch==1 && (*selectedLeptons)[0]->pt_>100 && (*selectedLeptons)[1]->pt_>100 && nbjet>0) cout<<"REZA :"<<MVAoutput<<" - "<<ev_run<<","<<ev_luminosityBlock<<","<<ev_event<<endl;
//fill histograms
//dilepton step
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
    Hists[ch][0][15]->Fill(MET_FinalCollection_Pt,weight_lep);
    Hists[ch][0][16]->Fill(MET_FinalCollection_phi,weight_lep);
    Hists[ch][0][17]->Fill(pv_n,weight_lep);
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

    for (int n=0;n<9;++n){
      HistsSysUp[ch][0][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      if(selectedJets->size()>0) HistsSysUp[ch][0][10][n]->Fill((*selectedJets)[0]->pt_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      if(selectedJets->size()>0) HistsSysUp[ch][0][11][n]->Fill((*selectedJets)[0]->eta_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      if(selectedJets->size()>0) HistsSysUp[ch][0][12][n]->Fill((*selectedJets)[0]->phi_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][13][n]->Fill(selectedJets->size(),weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][14][n]->Fill(nbjet,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][15][n]->Fill(MET_FinalCollection_Pt,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][16][n]->Fill(MET_FinalCollection_phi,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][17][n]->Fill(pv_n,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][0][19][n]->Fill(MVAoutput,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][0][20]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep);
      if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][0][22]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep);
      if ((*selectedLeptons)[0]->lep_ == 1)  Hists[ch][0][21]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep);
      if ((*selectedLeptons)[0]->lep_ == 1)  Hists[ch][0][23]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][0][20]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][0][22]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 1)  Hists[ch][0][21]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 1)  Hists[ch][0][23]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep);

      HistsSysDown[ch][0][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      if(selectedJets->size()>0) HistsSysDown[ch][0][10][n]->Fill((*selectedJets)[0]->pt_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      if(selectedJets->size()>0) HistsSysDown[ch][0][11][n]->Fill((*selectedJets)[0]->eta_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      if(selectedJets->size()>0) HistsSysDown[ch][0][12][n]->Fill((*selectedJets)[0]->phi_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][13][n]->Fill(selectedJets->size(),weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][14][n]->Fill(nbjet,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][15][n]->Fill(MET_FinalCollection_Pt,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][16][n]->Fill(MET_FinalCollection_phi,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][17][n]->Fill(pv_n,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][0][19][n]->Fill(MVAoutput,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][0][20]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep);
      if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][0][22]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep);
      if ((*selectedLeptons)[0]->lep_ == 1)  Hists[ch][0][21]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep);
      if ((*selectedLeptons)[0]->lep_ == 1)  Hists[ch][0][23]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][0][20]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][0][22]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 1)  Hists[ch][0][21]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 1)  Hists[ch][0][23]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep);
    }
    if(selectedJetsJesUp->size()>0) HistsSysUp[ch][0][10][9]->Fill((*selectedJetsJesUp)[0]->pt_,weight_lep );
    if(selectedJetsJesDown->size()>0) HistsSysDown[ch][0][10][9]->Fill((*selectedJetsJesDown)[0]->pt_,weight_lep );
    if(selectedJetsJerUp->size()>0) HistsSysUp[ch][0][10][10]->Fill((*selectedJetsJerUp)[0]->pt_,weight_lep );
    if(selectedJetsJerDown->size()>0) HistsSysDown[ch][0][10][10]->Fill((*selectedJetsJerDown)[0]->pt_,weight_lep );
    HistsSysUp[ch][0][13][9]->Fill(selectedJetsJesUp->size(),weight_lep );
    HistsSysDown[ch][0][13][9]->Fill(selectedJetsJesDown->size(),weight_lep );
    HistsSysUp[ch][0][13][10]->Fill(selectedJetsJerUp->size(),weight_lep );
    HistsSysDown[ch][0][13][10]->Fill(selectedJetsJerDown->size(),weight_lep );

  if (fname.Contains("TTTo2L2Nu")|| fname.Contains("LFV")){
     for (int n=0;n<reweightSizeQscalePDF;++n){
       HistsSysReweightsQscalePDF[ch][0][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       if(selectedJets->size()>0) HistsSysReweightsQscalePDF[ch][0][10][n]->Fill((*selectedJets)[0]->pt_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       if(selectedJets->size()>0) HistsSysReweightsQscalePDF[ch][0][11][n]->Fill((*selectedJets)[0]->eta_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       if(selectedJets->size()>0) HistsSysReweightsQscalePDF[ch][0][12][n]->Fill((*selectedJets)[0]->phi_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][13][n]->Fill(selectedJets->size(),weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][14][n]->Fill(nbjet,weight_lepB*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][15][n]->Fill(MET_FinalCollection_Pt,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][16][n]->Fill(MET_FinalCollection_phi,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][17][n]->Fill(pv_n,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][0][19][n]->Fill(MVAoutput,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
      if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][0][20]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep);
      if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][0][22]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep);
      if ((*selectedLeptons)[0]->lep_ == 1)  Hists[ch][0][21]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep);
      if ((*selectedLeptons)[0]->lep_ == 1)  Hists[ch][0][23]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][0][20]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][0][22]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 1)  Hists[ch][0][21]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 1)  Hists[ch][0][23]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep);
     }
     for (int n=0;n<reweightSizePS;++n){
       HistsSysReweightsPS[ch][0][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lep * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       if(selectedJets->size()>0) HistsSysReweightsPS[ch][0][10][n]->Fill((*selectedJets)[0]->pt_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       if(selectedJets->size()>0) HistsSysReweightsPS[ch][0][11][n]->Fill((*selectedJets)[0]->eta_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       if(selectedJets->size()>0) HistsSysReweightsPS[ch][0][12][n]->Fill((*selectedJets)[0]->phi_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][13][n]->Fill(selectedJets->size(),weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][14][n]->Fill(nbjet,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][15][n]->Fill(MET_FinalCollection_Pt,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][16][n]->Fill(MET_FinalCollection_phi,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][17][n]->Fill(pv_n,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][0][19][n]->Fill(MVAoutput,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
      if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][0][20]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep);
      if ((*selectedLeptons)[0]->lep_ == 10) Hists[ch][0][22]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep);
      if ((*selectedLeptons)[0]->lep_ == 1)  Hists[ch][0][21]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep);
      if ((*selectedLeptons)[0]->lep_ == 1)  Hists[ch][0][23]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][0][20]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 10) Hists[ch][0][22]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 1)  Hists[ch][0][21]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep);
      if ((*selectedLeptons)[1]->lep_ == 1)  Hists[ch][0][23]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep);
    }
  }

//exclude events in zmass window
  if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ != 11 && ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()<106 && ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()>76) continue;
//ll off Z region
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
  Hists[ch][1][15]->Fill(MET_FinalCollection_Pt,weight_lep);
  Hists[ch][1][16]->Fill(MET_FinalCollection_phi,weight_lep);
  Hists[ch][1][17]->Fill(pv_n,weight_lep);
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

  for (int n=0;n<9;++n){
    HistsSysUp[ch][1][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysUp[ch][1][10][n]->Fill((*selectedJets)[0]->pt_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysUp[ch][1][11][n]->Fill((*selectedJets)[0]->eta_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysUp[ch][1][12][n]->Fill((*selectedJets)[0]->phi_,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][13][n]->Fill(selectedJets->size(),weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][14][n]->Fill(nbjet,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][15][n]->Fill(MET_FinalCollection_Pt,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][16][n]->Fill(MET_FinalCollection_phi,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][17][n]->Fill(pv_n,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][19][n]->Fill(MVAoutput,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysUp[ch][1][20][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysUp[ch][1][22][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysUp[ch][1][21][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysUp[ch][1][23][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysUp[ch][1][20][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysUp[ch][1][22][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysUp[ch][1][21][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep * (sysUpWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysUp[ch][1][23][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep * (sysUpWeights[n]/nominalWeights[n]));


    HistsSysDown[ch][1][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysDown[ch][1][10][n]->Fill((*selectedJets)[0]->pt_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysDown[ch][1][11][n]->Fill((*selectedJets)[0]->eta_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysDown[ch][1][12][n]->Fill((*selectedJets)[0]->phi_,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][13][n]->Fill(selectedJets->size(),weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][14][n]->Fill(nbjet,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][15][n]->Fill(MET_FinalCollection_Pt,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][16][n]->Fill(MET_FinalCollection_phi,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][17][n]->Fill(pv_n,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][19][n]->Fill(MVAoutput,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysDown[ch][1][20][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysDown[ch][1][22][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysDown[ch][1][21][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysDown[ch][1][23][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysDown[ch][1][20][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysDown[ch][1][22][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysDown[ch][1][21][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysDown[ch][1][23][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lep * (sysDownWeights[n]/nominalWeights[n]));
  }
  if(selectedJetsJesUp->size()>0) {
    HistsSysUp[ch][1][10][9]->Fill((*selectedJetsJesUp)[0]->pt_,weight_lep );
    HistsSysUp[ch][1][15][9]->Fill(MET_T1SmearJetEnUp_Pt,weight_lep );
    HistsSysUp[ch][1][19][9]->Fill(MVAoutput,weight_lep );
  }
  if(selectedJetsJesDown->size()>0){
    HistsSysDown[ch][1][10][9]->Fill((*selectedJetsJesDown)[0]->pt_,weight_lep );
    HistsSysDown[ch][1][15][9]->Fill(MET_T1SmearJetEnDown_Pt,weight_lep );
    HistsSysDown[ch][1][19][9]->Fill(MVAoutput,weight_lep );
  }
  if(selectedJetsJerUp->size()>0){
    HistsSysUp[ch][1][10][10]->Fill((*selectedJetsJerUp)[0]->pt_,weight_lep );
    HistsSysUp[ch][1][15][10]->Fill(MET_T1SmearJetResUp_Pt,weight_lep );
    HistsSysUp[ch][1][19][10]->Fill(MVAoutput,weight_lep );
  }
  if(selectedJetsJerDown->size()>0){
    HistsSysDown[ch][1][10][10]->Fill((*selectedJetsJerDown)[0]->pt_,weight_lep );
    HistsSysDown[ch][1][15][10]->Fill(MET_T1SmearJetResDown_Pt,weight_lep );
    HistsSysDown[ch][1][19][10]->Fill(MVAoutput,weight_lep );
  }
  HistsSysUp[ch][1][13][9]->Fill(selectedJetsJesUp->size(),weight_lep );
  HistsSysDown[ch][1][13][9]->Fill(selectedJetsJesDown->size(),weight_lep );
  HistsSysUp[ch][1][13][10]->Fill(selectedJetsJerUp->size(),weight_lep );
  HistsSysDown[ch][1][13][10]->Fill(selectedJetsJerDown->size(),weight_lep );
  HistsSysUp[ch][1][15][11]->Fill(metUnclusMETUp,weight_lep );
  HistsSysUp[ch][1][19][11]->Fill(MVAoutputUnclusMETUp,weight_lep );
  HistsSysDown[ch][1][15][11]->Fill(metUnclusMETDown,weight_lep );
  HistsSysDown[ch][1][19][11]->Fill(MVAoutputUnclusMETDown,weight_lep );
  if (selectedLeptonsMuScaleUp->size()>1) HistsSysUp[ch][1][0][12]->Fill((*selectedLeptonsMuScaleUp)[0]->pt_,weight_lep );
  if (selectedLeptonsMuScaleDown->size()>1) HistsSysDown[ch][1][0][12]->Fill((*selectedLeptonsMuScaleDown)[0]->pt_,weight_lep );
  if (selectedLeptonsMuScaleUp->size()>1) HistsSysUp[ch][1][3][12]->Fill((*selectedLeptonsMuScaleUp)[1]->pt_,weight_lep );
  if (selectedLeptonsMuScaleDown->size()>1) HistsSysDown[ch][1][3][12]->Fill((*selectedLeptonsMuScaleDown)[1]->pt_,weight_lep );

  if (selectedLeptonsEleScaleUp->size()>1) HistsSysUp[ch][1][0][13]->Fill((*selectedLeptonsEleScaleUp)[0]->pt_,weight_lep );
  if (selectedLeptonsEleScaleDown->size()>1) HistsSysDown[ch][1][0][13]->Fill((*selectedLeptonsEleScaleDown)[0]->pt_,weight_lep );
  if (selectedLeptonsEleScaleUp->size()>1) HistsSysUp[ch][1][3][13]->Fill((*selectedLeptonsEleScaleUp)[1]->pt_,weight_lep );
  if (selectedLeptonsEleScaleDown->size()>1) HistsSysDown[ch][1][3][13]->Fill((*selectedLeptonsEleScaleDown)[1]->pt_,weight_lep );

  if (selectedLeptonsMuResUp->size()>1) HistsSysUp[ch][1][0][14]->Fill((*selectedLeptonsMuResUp)[0]->pt_,weight_lep );
  if (selectedLeptonsMuResDown->size()>1) HistsSysDown[ch][1][0][14]->Fill((*selectedLeptonsMuResDown)[0]->pt_,weight_lep );
  if (selectedLeptonsMuResUp->size()>1) HistsSysUp[ch][1][3][14]->Fill((*selectedLeptonsMuResUp)[1]->pt_,weight_lep );
  if (selectedLeptonsMuResDown->size()>1) HistsSysDown[ch][1][3][14]->Fill((*selectedLeptonsMuResDown)[1]->pt_,weight_lep );

  HistsSysUp[ch][1][19][12]->Fill(MVAoutputMuScaleUp,weight_lep );
  HistsSysUp[ch][1][19][13]->Fill(MVAoutputEleScaleUp,weight_lep );
  HistsSysUp[ch][1][19][14]->Fill(MVAoutputMuResUp,weight_lep );

  HistsSysDown[ch][1][19][12]->Fill(MVAoutputMuScaleDown,weight_lep );
  HistsSysDown[ch][1][19][13]->Fill(MVAoutputEleScaleDown,weight_lep );
  HistsSysDown[ch][1][19][14]->Fill(MVAoutputMuResDown,weight_lep );

  if (fname.Contains("TTTo2L2Nu")|| fname.Contains("LFV")){
     for (int n=0;n<reweightSizeQscalePDF;++n){
       HistsSysReweightsQscalePDF[ch][1][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       if(selectedJets->size()>0) HistsSysReweightsQscalePDF[ch][1][10][n]->Fill((*selectedJets)[0]->pt_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       if(selectedJets->size()>0) HistsSysReweightsQscalePDF[ch][1][11][n]->Fill((*selectedJets)[0]->eta_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       if(selectedJets->size()>0) HistsSysReweightsQscalePDF[ch][1][12][n]->Fill((*selectedJets)[0]->phi_,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][13][n]->Fill(selectedJets->size(),weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][14][n]->Fill(nbjet,weight_lepB*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][15][n]->Fill(MET_FinalCollection_Pt,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][16][n]->Fill(MET_FinalCollection_phi,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][17][n]->Fill(pv_n,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
       HistsSysReweightsQscalePDF[ch][1][19][n]->Fill(MVAoutput,weight_lep*(SLW[0]/SLW[n])*((*LHE_weight_sys)[n]/(*LHE_weight_sys)[0]));
     }
     for (int n=0;n<reweightSizePS;++n){
       HistsSysReweightsPS[ch][1][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lep * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       if(selectedJets->size()>0) HistsSysReweightsPS[ch][1][10][n]->Fill((*selectedJets)[0]->pt_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       if(selectedJets->size()>0) HistsSysReweightsPS[ch][1][11][n]->Fill((*selectedJets)[0]->eta_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       if(selectedJets->size()>0) HistsSysReweightsPS[ch][1][12][n]->Fill((*selectedJets)[0]->phi_,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][13][n]->Fill(selectedJets->size(),weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][14][n]->Fill(nbjet,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][15][n]->Fill(MET_FinalCollection_Pt,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][16][n]->Fill(MET_FinalCollection_phi,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][17][n]->Fill(pv_n,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
       HistsSysReweightsPS[ch][1][19][n]->Fill(MVAoutput,weight_lep* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
    }
  }

//************* SIGNAL region
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
    Hists[ch][2][15]->Fill(MET_FinalCollection_Pt,weight_lepB);
    Hists[ch][2][16]->Fill(MET_FinalCollection_phi,weight_lepB);
    Hists[ch][2][17]->Fill(pv_n,weight_lepB);
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

    for (int n=0;n<9;++n){
      HistsSysUp[ch][2][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][10][n]->Fill((*selectedJets)[0]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][11][n]->Fill((*selectedJets)[0]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][12][n]->Fill((*selectedJets)[0]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][13][n]->Fill(selectedJets->size(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][14][n]->Fill(nbjet,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][15][n]->Fill(MET_FinalCollection_Pt,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][16][n]->Fill(MET_FinalCollection_phi,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][17][n]->Fill(pv_n,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][2][19][n]->Fill(MVAoutput,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 10) HistsSysUp[ch][2][20][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 10) HistsSysUp[ch][2][22][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysUp[ch][2][21][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysUp[ch][2][23][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 10) HistsSysUp[ch][2][20][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 10) HistsSysUp[ch][2][22][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysUp[ch][2][21][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysUp[ch][2][23][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));

      HistsSysDown[ch][2][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][10][n]->Fill((*selectedJets)[0]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][11][n]->Fill((*selectedJets)[0]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][12][n]->Fill((*selectedJets)[0]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][13][n]->Fill(selectedJets->size(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][14][n]->Fill(nbjet,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][15][n]->Fill(MET_FinalCollection_Pt,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][16][n]->Fill(MET_FinalCollection_phi,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][17][n]->Fill(pv_n,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][2][19][n]->Fill(MVAoutput,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 10) HistsSysDown[ch][2][20][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 10) HistsSysDown[ch][2][22][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysDown[ch][2][21][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysDown[ch][2][23][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 10) HistsSysDown[ch][2][20][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 10) HistsSysDown[ch][2][22][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysDown[ch][2][21][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysDown[ch][2][23][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    }
    HistsSysUp[ch][2][15][11]->Fill(metUnclusMETUp,weight_lepB );
    HistsSysUp[ch][2][19][11]->Fill(MVAoutputUnclusMETUp,weight_lepB );
    HistsSysDown[ch][2][15][11]->Fill(metUnclusMETDown,weight_lepB );
    HistsSysDown[ch][2][19][11]->Fill(MVAoutputUnclusMETDown,weight_lepB );
    if (selectedLeptonsMuScaleUp->size()>1) {
      HistsSysUp[ch][2][0][12]->Fill((*selectedLeptonsMuScaleUp)[0]->pt_,weight_lepB );
      HistsSysUp[ch][2][3][12]->Fill((*selectedLeptonsMuScaleUp)[1]->pt_,weight_lepB );
      if ((*selectedLeptonsMuScaleUp)[0]->lep_ == 10) HistsSysUp[ch][2][20][12]->Fill((*selectedLeptonsMuScaleUp)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuScaleUp)[0]->lep_ == 1)  HistsSysUp[ch][2][21][12]->Fill((*selectedLeptonsMuScaleUp)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuScaleUp)[1]->lep_ == 10) HistsSysUp[ch][2][20][12]->Fill((*selectedLeptonsMuScaleUp)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuScaleUp)[1]->lep_ == 1)  HistsSysUp[ch][2][21][12]->Fill((*selectedLeptonsMuScaleUp)[1]->pt_   ,weight_lepB);
    }

    if (selectedLeptonsMuScaleDown->size()>1) {
       HistsSysDown[ch][2][0][12]->Fill((*selectedLeptonsMuScaleDown)[0]->pt_,weight_lepB );
       HistsSysDown[ch][2][3][12]->Fill((*selectedLeptonsMuScaleDown)[1]->pt_,weight_lepB );
      if ((*selectedLeptonsMuScaleDown)[0]->lep_ == 10) HistsSysDown[ch][2][20][12]->Fill((*selectedLeptonsMuScaleDown)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuScaleDown)[0]->lep_ == 1)  HistsSysDown[ch][2][21][12]->Fill((*selectedLeptonsMuScaleDown)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuScaleDown)[1]->lep_ == 10) HistsSysDown[ch][2][20][12]->Fill((*selectedLeptonsMuScaleDown)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuScaleDown)[1]->lep_ == 1)  HistsSysDown[ch][2][21][12]->Fill((*selectedLeptonsMuScaleDown)[1]->pt_   ,weight_lepB);
    }

    if (selectedLeptonsEleScaleUp->size()>1) {
      HistsSysUp[ch][2][0][13]->Fill((*selectedLeptonsEleScaleUp)[0]->pt_,weight_lepB );
      HistsSysUp[ch][2][3][13]->Fill((*selectedLeptonsEleScaleUp)[1]->pt_,weight_lepB );
      if ((*selectedLeptonsEleScaleUp)[0]->lep_ == 10) HistsSysUp[ch][2][20][13]->Fill((*selectedLeptonsEleScaleUp)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsEleScaleUp)[0]->lep_ == 1)  HistsSysUp[ch][2][21][13]->Fill((*selectedLeptonsEleScaleUp)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsEleScaleUp)[1]->lep_ == 10) HistsSysUp[ch][2][20][13]->Fill((*selectedLeptonsEleScaleUp)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptonsEleScaleUp)[1]->lep_ == 1)  HistsSysUp[ch][2][21][13]->Fill((*selectedLeptonsEleScaleUp)[1]->pt_   ,weight_lepB);
    }

    if (selectedLeptonsEleScaleDown->size()>1) {
      HistsSysDown[ch][2][0][13]->Fill((*selectedLeptonsEleScaleDown)[0]->pt_,weight_lepB );
      HistsSysDown[ch][2][3][13]->Fill((*selectedLeptonsEleScaleDown)[1]->pt_,weight_lepB );
      if ((*selectedLeptonsEleScaleDown)[0]->lep_ == 10) HistsSysDown[ch][2][20][13]->Fill((*selectedLeptonsEleScaleDown)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsEleScaleDown)[0]->lep_ == 1)  HistsSysDown[ch][2][21][13]->Fill((*selectedLeptonsEleScaleDown)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsEleScaleDown)[1]->lep_ == 10) HistsSysDown[ch][2][20][13]->Fill((*selectedLeptonsEleScaleDown)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptonsEleScaleDown)[1]->lep_ == 1)  HistsSysDown[ch][2][21][13]->Fill((*selectedLeptonsEleScaleDown)[1]->pt_   ,weight_lepB);
    }

    if (selectedLeptonsMuResUp->size()>1) {
      HistsSysUp[ch][2][0][14]->Fill((*selectedLeptonsMuResUp)[0]->pt_,weight_lepB );
      HistsSysUp[ch][2][3][14]->Fill((*selectedLeptonsMuResUp)[1]->pt_,weight_lepB );
      if ((*selectedLeptonsMuResUp)[0]->lep_ == 10) HistsSysUp[ch][2][20][14]->Fill((*selectedLeptonsMuResUp)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuResUp)[0]->lep_ == 1)  HistsSysUp[ch][2][21][14]->Fill((*selectedLeptonsMuResUp)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuResUp)[1]->lep_ == 10) HistsSysUp[ch][2][20][14]->Fill((*selectedLeptonsMuResUp)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuResUp)[1]->lep_ == 1)  HistsSysUp[ch][2][21][14]->Fill((*selectedLeptonsMuResUp)[1]->pt_   ,weight_lepB);
    }

    if (selectedLeptonsMuResDown->size()>1) {
      HistsSysDown[ch][2][0][14]->Fill((*selectedLeptonsMuResDown)[0]->pt_,weight_lepB );
      HistsSysDown[ch][2][3][14]->Fill((*selectedLeptonsMuResDown)[1]->pt_,weight_lepB );
      if ((*selectedLeptonsMuResDown)[0]->lep_ == 10) HistsSysDown[ch][2][20][14]->Fill((*selectedLeptonsMuResDown)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuResDown)[0]->lep_ == 1)  HistsSysDown[ch][2][21][14]->Fill((*selectedLeptonsMuResDown)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuResDown)[1]->lep_ == 10) HistsSysDown[ch][2][20][14]->Fill((*selectedLeptonsMuResDown)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuResDown)[1]->lep_ == 1)  HistsSysDown[ch][2][21][14]->Fill((*selectedLeptonsMuResDown)[1]->pt_   ,weight_lepB);
    }
  
    HistsSysUp[ch][2][19][12]->Fill(MVAoutputMuScaleUp,weight_lepB );
    HistsSysUp[ch][2][19][13]->Fill(MVAoutputEleScaleUp,weight_lepB );
    HistsSysUp[ch][2][19][14]->Fill(MVAoutputMuResUp,weight_lepB );
  
    HistsSysDown[ch][2][19][12]->Fill(MVAoutputMuScaleDown,weight_lepB );
    HistsSysDown[ch][2][19][13]->Fill(MVAoutputEleScaleDown,weight_lepB );
    HistsSysDown[ch][2][19][14]->Fill(MVAoutputMuResDown,weight_lepB );
    if (fname.Contains("TTTo2L2Nu")|| fname.Contains("LFV")){
       for (int n=0;n<reweightSizeQscalePDF;++n){
         R = n;
//         if(fname.Contains("LFV") && n>49) R = 611 + n -50; 
         HistsSysReweightsQscalePDF[ch][2][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][10][n]->Fill((*selectedJets)[0]->pt_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][11][n]->Fill((*selectedJets)[0]->eta_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][12][n]->Fill((*selectedJets)[0]->phi_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][13][n]->Fill(selectedJets->size(),weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][14][n]->Fill(nbjet,weight_lepB*(SLW[0]/SLW[n])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[R]));
         HistsSysReweightsQscalePDF[ch][2][15][n]->Fill(MET_FinalCollection_Pt,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][16][n]->Fill(MET_FinalCollection_phi,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][17][n]->Fill(pv_n,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][2][19][n]->Fill(MVAoutput,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 10) HistsSysReweightsQscalePDF[ch][2][20][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 10) HistsSysReweightsQscalePDF[ch][2][22][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysReweightsQscalePDF[ch][2][21][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysReweightsQscalePDF[ch][2][23][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 10) HistsSysReweightsQscalePDF[ch][2][20][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 10) HistsSysReweightsQscalePDF[ch][2][22][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysReweightsQscalePDF[ch][2][21][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysReweightsQscalePDF[ch][2][23][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
       }
       for (int n=0;n<reweightSizePS;++n){
         if ((*gen_weight_sys)[n]/(*gen_weight_sys)[0]>20) continue;
         HistsSysReweightsPS[ch][2][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][10][n]->Fill((*selectedJets)[0]->pt_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][11][n]->Fill((*selectedJets)[0]->eta_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][12][n]->Fill((*selectedJets)[0]->phi_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][13][n]->Fill(selectedJets->size(),weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][14][n]->Fill(nbjet,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][15][n]->Fill(MET_FinalCollection_Pt,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][16][n]->Fill(MET_FinalCollection_phi,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][17][n]->Fill(pv_n,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][2][19][n]->Fill(MVAoutput,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 10) HistsSysReweightsPS[ch][2][20][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 10) HistsSysReweightsPS[ch][2][22][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysReweightsPS[ch][2][21][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysReweightsPS[ch][2][23][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 10) HistsSysReweightsPS[ch][2][20][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 10) HistsSysReweightsPS[ch][2][22][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysReweightsPS[ch][2][21][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysReweightsPS[ch][2][23][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
      }
    }
  }
  if(nbjetJesUp ==1) {
    HistsSysUp[ch][2][0][9]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    HistsSysUp[ch][2][1][9]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    HistsSysUp[ch][2][2][9]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    HistsSysUp[ch][2][3][9]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    HistsSysUp[ch][2][4][9]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    HistsSysUp[ch][2][5][9]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    HistsSysUp[ch][2][6][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysUp[ch][2][7][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB);
    HistsSysUp[ch][2][8][9]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysUp[ch][2][9][9]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysUp[ch][2][16][9]->Fill(MET_FinalCollection_phi,weight_lepB);
    HistsSysUp[ch][2][17][9]->Fill(pv_n,weight_lepB);
    HistsSysUp[ch][2][18][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysUp[ch][2][10][9]->Fill((*selectedJetsJesUp)[0]->pt_,weight_lepB );
    HistsSysUp[ch][2][11][9]->Fill((*selectedJetsJesUp)[0]->eta_,weight_lepB);
    HistsSysUp[ch][2][12][9]->Fill((*selectedJetsJesUp)[0]->phi_,weight_lepB);
    HistsSysUp[ch][2][13][9]->Fill(selectedJetsJesUp->size(),weight_lepB );
    HistsSysUp[ch][2][14][9]->Fill(nbjetJesUp,weight_lepB);
    HistsSysUp[ch][2][15][9]->Fill(MET_T1SmearJetEnUp_Pt,weight_lepB );
    HistsSysUp[ch][2][19][9]->Fill(MVAoutputJesUp,weight_lepB );
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysUp[ch][2][20][9]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysUp[ch][2][22][9]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysUp[ch][2][21][9]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysUp[ch][2][23][9]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysUp[ch][2][20][9]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysUp[ch][2][22][9]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysUp[ch][2][21][9]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysUp[ch][2][23][9]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
  }
  if(nbjetJesDown==1){
    HistsSysDown[ch][2][0][9]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    HistsSysDown[ch][2][1][9]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    HistsSysDown[ch][2][2][9]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    HistsSysDown[ch][2][3][9]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    HistsSysDown[ch][2][4][9]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    HistsSysDown[ch][2][5][9]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    HistsSysDown[ch][2][6][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysDown[ch][2][7][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB);
    HistsSysDown[ch][2][8][9]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysDown[ch][2][9][9]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysDown[ch][2][16][9]->Fill(MET_FinalCollection_phi,weight_lepB);
    HistsSysDown[ch][2][17][9]->Fill(pv_n,weight_lepB);
    HistsSysDown[ch][2][18][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysDown[ch][2][10][9]->Fill((*selectedJetsJesDown)[0]->pt_,weight_lepB );
    HistsSysDown[ch][2][11][9]->Fill((*selectedJetsJesDown)[0]->eta_,weight_lepB);
    HistsSysDown[ch][2][12][9]->Fill((*selectedJetsJesDown)[0]->phi_,weight_lepB);
    HistsSysDown[ch][2][13][9]->Fill(selectedJetsJesDown->size(),weight_lepB );
    HistsSysDown[ch][2][14][9]->Fill(nbjetJesDown,weight_lepB);
    HistsSysDown[ch][2][15][9]->Fill(MET_T1SmearJetEnDown_Pt,weight_lepB );
    HistsSysDown[ch][2][19][9]->Fill(MVAoutputJesDown,weight_lepB );
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysDown[ch][2][20][9]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysDown[ch][2][22][9]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysDown[ch][2][21][9]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysDown[ch][2][23][9]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysDown[ch][2][20][9]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysDown[ch][2][22][9]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysDown[ch][2][21][9]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysDown[ch][2][23][9]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
  }
  if(nbjetJerUp==1){
    HistsSysUp[ch][2][0][10]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    HistsSysUp[ch][2][1][10]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    HistsSysUp[ch][2][2][10]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    HistsSysUp[ch][2][3][10]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    HistsSysUp[ch][2][4][10]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    HistsSysUp[ch][2][5][10]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    HistsSysUp[ch][2][6][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysUp[ch][2][7][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB);
    HistsSysUp[ch][2][8][10]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysUp[ch][2][9][10]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysUp[ch][2][16][10]->Fill(MET_FinalCollection_phi,weight_lepB);
    HistsSysUp[ch][2][17][10]->Fill(pv_n,weight_lepB);
    HistsSysUp[ch][2][18][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysUp[ch][2][10][10]->Fill((*selectedJetsJerUp)[0]->pt_,weight_lepB );
    HistsSysUp[ch][2][11][10]->Fill((*selectedJetsJerUp)[0]->eta_,weight_lepB);
    HistsSysUp[ch][2][12][10]->Fill((*selectedJetsJerUp)[0]->phi_,weight_lepB);
    HistsSysUp[ch][2][13][10]->Fill(selectedJetsJerUp->size(),weight_lepB );
    HistsSysUp[ch][2][14][10]->Fill(nbjetJerUp,weight_lepB);
    HistsSysUp[ch][2][15][10]->Fill(MET_T1SmearJetResUp_Pt,weight_lepB );
    HistsSysUp[ch][2][19][10]->Fill(MVAoutputJerUp,weight_lepB );
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysUp[ch][2][20][10]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysUp[ch][2][22][10]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysUp[ch][2][21][10]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysUp[ch][2][23][10]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysUp[ch][2][20][10]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysUp[ch][2][22][10]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysUp[ch][2][21][10]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysUp[ch][2][23][10]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
  }
  if(nbjetJerDown==1){
    HistsSysDown[ch][2][0][10]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    HistsSysDown[ch][2][1][10]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    HistsSysDown[ch][2][2][10]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    HistsSysDown[ch][2][3][10]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    HistsSysDown[ch][2][4][10]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    HistsSysDown[ch][2][5][10]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    HistsSysDown[ch][2][6][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysDown[ch][2][7][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB);
    HistsSysDown[ch][2][8][10]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysDown[ch][2][9][10]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysDown[ch][2][16][10]->Fill(MET_FinalCollection_phi,weight_lepB);
    HistsSysDown[ch][2][17][10]->Fill(pv_n,weight_lepB);
    HistsSysDown[ch][2][18][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysDown[ch][2][10][10]->Fill((*selectedJetsJerDown)[0]->pt_,weight_lepB );
    HistsSysDown[ch][2][11][10]->Fill((*selectedJetsJerDown)[0]->eta_,weight_lepB);
    HistsSysDown[ch][2][12][10]->Fill((*selectedJetsJerDown)[0]->phi_,weight_lepB);
    HistsSysDown[ch][2][13][10]->Fill(selectedJetsJerDown->size(),weight_lepB );
    HistsSysDown[ch][2][14][10]->Fill(nbjetJerDown,weight_lepB);
    HistsSysDown[ch][2][15][10]->Fill(MET_T1SmearJetResDown_Pt,weight_lepB );
    HistsSysDown[ch][2][19][10]->Fill(MVAoutputJerDown,weight_lepB );
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysDown[ch][2][20][10]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysDown[ch][2][22][10]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysDown[ch][2][21][10]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysDown[ch][2][23][10]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysDown[ch][2][20][10]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysDown[ch][2][22][10]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysDown[ch][2][21][10]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysDown[ch][2][23][10]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
  }
  for (int n=0;n<sysJecNames.size();++n){
    if ((*JECsysNbtagUp)[n]==1) {
      HistsJECUp[ch][2][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
      HistsJECUp[ch][2][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
      HistsJECUp[ch][2][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
      HistsJECUp[ch][2][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      HistsJECUp[ch][2][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      HistsJECUp[ch][2][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
      HistsJECUp[ch][2][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
      HistsJECUp[ch][2][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB);
      HistsJECUp[ch][2][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB);
      HistsJECUp[ch][2][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB);
      HistsJECUp[ch][2][16][n]->Fill(MET_FinalCollection_phi,weight_lepB);
      HistsJECUp[ch][2][17][n]->Fill(pv_n,weight_lepB);
      HistsJECUp[ch][2][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
      HistsJECUp[ch][2][10][n]->Fill((*JECsysUp)[n][0]->pt_,weight_lepB );
      HistsJECUp[ch][2][11][n]->Fill((*JECsysUp)[n][0]->eta_,weight_lepB);
      HistsJECUp[ch][2][12][n]->Fill((*JECsysUp)[n][0]->phi_,weight_lepB);
      HistsJECUp[ch][2][13][n]->Fill((*JECsysUp)[n].size(),weight_lepB );
      HistsJECUp[ch][2][14][n]->Fill((*JECsysNbtagUp)[n],weight_lepB);
      HistsJECUp[ch][2][15][n]->Fill((*JECsysMETUp)[n],weight_lepB );
      HistsJECUp[ch][2][19][n]->Fill((*JECsysMVAUp)[n],weight_lepB );    
      if ((*selectedLeptons)[0]->lep_ == 10) HistsJECUp[ch][2][20][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[0]->lep_ == 10) HistsJECUp[ch][2][22][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsJECUp[ch][2][21][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsJECUp[ch][2][23][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 10) HistsJECUp[ch][2][20][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 10) HistsJECUp[ch][2][22][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsJECUp[ch][2][21][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsJECUp[ch][2][23][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB); 
      }
    if ((*JECsysNbtagDown)[n]==1) {
      HistsJECDown[ch][2][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
      HistsJECDown[ch][2][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
      HistsJECDown[ch][2][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
      HistsJECDown[ch][2][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      HistsJECDown[ch][2][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      HistsJECDown[ch][2][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
      HistsJECDown[ch][2][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
      HistsJECDown[ch][2][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB);
      HistsJECDown[ch][2][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB);
      HistsJECDown[ch][2][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB);
      HistsJECDown[ch][2][16][n]->Fill(MET_FinalCollection_phi,weight_lepB);
      HistsJECDown[ch][2][17][n]->Fill(pv_n,weight_lepB);
      HistsJECDown[ch][2][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
      HistsJECDown[ch][2][10][n]->Fill((*JECsysDown)[n][0]->pt_,weight_lepB );
      HistsJECDown[ch][2][11][n]->Fill((*JECsysDown)[n][0]->eta_,weight_lepB);
      HistsJECDown[ch][2][12][n]->Fill((*JECsysDown)[n][0]->phi_,weight_lepB);
      HistsJECDown[ch][2][13][n]->Fill((*JECsysDown)[n].size(),weight_lepB );
      HistsJECDown[ch][2][14][n]->Fill((*JECsysNbtagDown)[n],weight_lepB);
      HistsJECDown[ch][2][15][n]->Fill((*JECsysMETDown)[n],weight_lepB );
      HistsJECDown[ch][2][19][n]->Fill((*JECsysMVADown)[n],weight_lepB );
      if ((*selectedLeptons)[0]->lep_ == 10) HistsJECDown[ch][2][20][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[0]->lep_ == 10) HistsJECDown[ch][2][22][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsJECDown[ch][2][21][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsJECDown[ch][2][23][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 10) HistsJECDown[ch][2][20][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 10) HistsJECDown[ch][2][22][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsJECDown[ch][2][21][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsJECDown[ch][2][23][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
    }
  }


///////// TT control region
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
    Hists[ch][3][15]->Fill(MET_FinalCollection_Pt,weight_lepB);
    Hists[ch][3][16]->Fill(MET_FinalCollection_phi,weight_lepB);
    Hists[ch][3][17]->Fill(pv_n,weight_lepB);
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
    for (int n=0;n<9;++n){
      HistsSysUp[ch][3][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][10][n]->Fill((*selectedJets)[0]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][11][n]->Fill((*selectedJets)[0]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][12][n]->Fill((*selectedJets)[0]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][13][n]->Fill(selectedJets->size(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][14][n]->Fill(nbjet,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][15][n]->Fill(MET_FinalCollection_Pt,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][16][n]->Fill(MET_FinalCollection_phi,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][17][n]->Fill(pv_n,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      HistsSysUp[ch][3][19][n]->Fill(MVAoutput,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 10) HistsSysUp[ch][3][20][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 10) HistsSysUp[ch][3][22][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysUp[ch][3][21][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysUp[ch][3][23][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 10) HistsSysUp[ch][3][20][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 10) HistsSysUp[ch][3][22][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysUp[ch][3][21][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysUp[ch][3][23][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));

      HistsSysDown[ch][3][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][10][n]->Fill((*selectedJets)[0]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][11][n]->Fill((*selectedJets)[0]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][12][n]->Fill((*selectedJets)[0]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][13][n]->Fill(selectedJets->size(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][14][n]->Fill(nbjet,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][15][n]->Fill(MET_FinalCollection_Pt,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][16][n]->Fill(MET_FinalCollection_phi,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][17][n]->Fill(pv_n,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      HistsSysDown[ch][3][19][n]->Fill(MVAoutput,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 10) HistsSysDown[ch][3][20][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 10) HistsSysDown[ch][3][22][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysDown[ch][3][21][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysDown[ch][3][23][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 10) HistsSysDown[ch][3][20][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 10) HistsSysDown[ch][3][22][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysDown[ch][3][21][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysDown[ch][3][23][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    }

    HistsSysUp[ch][3][15][11]->Fill(metUnclusMETUp,weight_lepB );
    HistsSysUp[ch][3][19][11]->Fill(MVAoutputUnclusMETUp,weight_lepB );
    HistsSysDown[ch][3][15][11]->Fill(metUnclusMETDown,weight_lepB );
    HistsSysDown[ch][3][19][11]->Fill(MVAoutputUnclusMETDown,weight_lepB );

    if (selectedLeptonsMuScaleUp->size()>1) {
      HistsSysUp[ch][3][0][12]->Fill((*selectedLeptonsMuScaleUp)[0]->pt_,weight_lepB );
      HistsSysUp[ch][3][3][12]->Fill((*selectedLeptonsMuScaleUp)[1]->pt_,weight_lepB );
      if ((*selectedLeptonsMuScaleUp)[0]->lep_ == 10) HistsSysUp[ch][3][20][12]->Fill((*selectedLeptonsMuScaleUp)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuScaleUp)[0]->lep_ == 1)  HistsSysUp[ch][3][21][12]->Fill((*selectedLeptonsMuScaleUp)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuScaleUp)[1]->lep_ == 10) HistsSysUp[ch][3][20][12]->Fill((*selectedLeptonsMuScaleUp)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuScaleUp)[1]->lep_ == 1)  HistsSysUp[ch][3][21][12]->Fill((*selectedLeptonsMuScaleUp)[1]->pt_   ,weight_lepB);
    }

    if (selectedLeptonsMuScaleDown->size()>1) {
       HistsSysDown[ch][3][0][12]->Fill((*selectedLeptonsMuScaleDown)[0]->pt_,weight_lepB );
       HistsSysDown[ch][3][3][12]->Fill((*selectedLeptonsMuScaleDown)[1]->pt_,weight_lepB );
      if ((*selectedLeptonsMuScaleDown)[0]->lep_ == 10) HistsSysDown[ch][3][20][12]->Fill((*selectedLeptonsMuScaleDown)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuScaleDown)[0]->lep_ == 1)  HistsSysDown[ch][3][21][12]->Fill((*selectedLeptonsMuScaleDown)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuScaleDown)[1]->lep_ == 10) HistsSysDown[ch][3][20][12]->Fill((*selectedLeptonsMuScaleDown)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuScaleDown)[1]->lep_ == 1)  HistsSysDown[ch][3][21][12]->Fill((*selectedLeptonsMuScaleDown)[1]->pt_   ,weight_lepB);
    }

    if (selectedLeptonsEleScaleUp->size()>1) {
      HistsSysUp[ch][3][0][13]->Fill((*selectedLeptonsEleScaleUp)[0]->pt_,weight_lepB );
      HistsSysUp[ch][3][3][13]->Fill((*selectedLeptonsEleScaleUp)[1]->pt_,weight_lepB );
      if ((*selectedLeptonsEleScaleUp)[0]->lep_ == 10) HistsSysUp[ch][3][20][13]->Fill((*selectedLeptonsEleScaleUp)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsEleScaleUp)[0]->lep_ == 1)  HistsSysUp[ch][3][21][13]->Fill((*selectedLeptonsEleScaleUp)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsEleScaleUp)[1]->lep_ == 10) HistsSysUp[ch][3][20][13]->Fill((*selectedLeptonsEleScaleUp)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptonsEleScaleUp)[1]->lep_ == 1)  HistsSysUp[ch][3][21][13]->Fill((*selectedLeptonsEleScaleUp)[1]->pt_   ,weight_lepB);
    }

    if (selectedLeptonsEleScaleDown->size()>1) {
      HistsSysDown[ch][3][0][13]->Fill((*selectedLeptonsEleScaleDown)[0]->pt_,weight_lepB );
      HistsSysDown[ch][3][3][13]->Fill((*selectedLeptonsEleScaleDown)[1]->pt_,weight_lepB );
      if ((*selectedLeptonsEleScaleDown)[0]->lep_ == 10) HistsSysDown[ch][3][20][13]->Fill((*selectedLeptonsEleScaleDown)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsEleScaleDown)[0]->lep_ == 1)  HistsSysDown[ch][3][21][13]->Fill((*selectedLeptonsEleScaleDown)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsEleScaleDown)[1]->lep_ == 10) HistsSysDown[ch][3][20][13]->Fill((*selectedLeptonsEleScaleDown)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptonsEleScaleDown)[1]->lep_ == 1)  HistsSysDown[ch][3][21][13]->Fill((*selectedLeptonsEleScaleDown)[1]->pt_   ,weight_lepB);
    }

    if (selectedLeptonsMuResUp->size()>1) {
      HistsSysUp[ch][3][0][14]->Fill((*selectedLeptonsMuResUp)[0]->pt_,weight_lepB );
      HistsSysUp[ch][3][3][14]->Fill((*selectedLeptonsMuResUp)[1]->pt_,weight_lepB );
      if ((*selectedLeptonsMuResUp)[0]->lep_ == 10) HistsSysUp[ch][3][20][14]->Fill((*selectedLeptonsMuResUp)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuResUp)[0]->lep_ == 1)  HistsSysUp[ch][3][21][14]->Fill((*selectedLeptonsMuResUp)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuResUp)[1]->lep_ == 10) HistsSysUp[ch][3][20][14]->Fill((*selectedLeptonsMuResUp)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuResUp)[1]->lep_ == 1)  HistsSysUp[ch][3][21][14]->Fill((*selectedLeptonsMuResUp)[1]->pt_   ,weight_lepB);
    }

    if (selectedLeptonsMuResDown->size()>1) {
      HistsSysDown[ch][3][0][14]->Fill((*selectedLeptonsMuResDown)[0]->pt_,weight_lepB );
      HistsSysDown[ch][3][3][14]->Fill((*selectedLeptonsMuResDown)[1]->pt_,weight_lepB );
      if ((*selectedLeptonsMuResDown)[0]->lep_ == 10) HistsSysDown[ch][3][20][14]->Fill((*selectedLeptonsMuResDown)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuResDown)[0]->lep_ == 1)  HistsSysDown[ch][3][21][14]->Fill((*selectedLeptonsMuResDown)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuResDown)[1]->lep_ == 10) HistsSysDown[ch][3][20][14]->Fill((*selectedLeptonsMuResDown)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptonsMuResDown)[1]->lep_ == 1)  HistsSysDown[ch][3][21][14]->Fill((*selectedLeptonsMuResDown)[1]->pt_   ,weight_lepB);
    }

    HistsSysUp[ch][3][19][12]->Fill(MVAoutputMuScaleUp,weight_lepB );
    HistsSysUp[ch][3][19][13]->Fill(MVAoutputEleScaleUp,weight_lepB );
    HistsSysUp[ch][3][19][14]->Fill(MVAoutputMuResUp,weight_lepB );

    HistsSysDown[ch][3][19][12]->Fill(MVAoutputMuScaleDown,weight_lepB );
    HistsSysDown[ch][3][19][13]->Fill(MVAoutputEleScaleDown,weight_lepB );
    HistsSysDown[ch][3][19][14]->Fill(MVAoutputMuResDown,weight_lepB );

    if (fname.Contains("TTTo2L2Nu")|| fname.Contains("LFV")){
       for (int n=0;n<reweightSizeQscalePDF;++n){
         R = n;
//         if(fname.Contains("LFV") && n>49) R = 611 + n -50;
         HistsSysReweightsQscalePDF[ch][3][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][10][n]->Fill((*selectedJets)[0]->pt_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][11][n]->Fill((*selectedJets)[0]->eta_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][12][n]->Fill((*selectedJets)[0]->phi_,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][13][n]->Fill(selectedJets->size(),weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][14][n]->Fill(nbjet,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][15][n]->Fill(MET_FinalCollection_Pt,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][16][n]->Fill(MET_FinalCollection_phi,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][17][n]->Fill(pv_n,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         HistsSysReweightsQscalePDF[ch][3][19][n]->Fill(MVAoutput,weight_lepB*(SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 10) HistsSysReweightsQscalePDF[ch][3][20][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 10) HistsSysReweightsQscalePDF[ch][3][22][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysReweightsQscalePDF[ch][3][21][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysReweightsQscalePDF[ch][3][23][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 10) HistsSysReweightsQscalePDF[ch][3][20][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 10) HistsSysReweightsQscalePDF[ch][3][22][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysReweightsQscalePDF[ch][3][21][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysReweightsQscalePDF[ch][3][23][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (SLW[0]/SLW[R])*((*LHE_weight_sys)[R]/(*LHE_weight_sys)[0]));
       }
       for (int n=0;n<reweightSizePS;++n){
         if ((*gen_weight_sys)[n]/(*gen_weight_sys)[0]>20) continue;
         HistsSysReweightsPS[ch][3][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][10][n]->Fill((*selectedJets)[0]->pt_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][11][n]->Fill((*selectedJets)[0]->eta_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][12][n]->Fill((*selectedJets)[0]->phi_,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][13][n]->Fill(selectedJets->size(),weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][14][n]->Fill(nbjet,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][15][n]->Fill(MET_FinalCollection_Pt,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][16][n]->Fill(MET_FinalCollection_phi,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][17][n]->Fill(pv_n,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         HistsSysReweightsPS[ch][3][19][n]->Fill(MVAoutput,weight_lepB* (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 10) HistsSysReweightsPS[ch][3][20][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 10) HistsSysReweightsPS[ch][3][22][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysReweightsPS[ch][3][21][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysReweightsPS[ch][3][23][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 10) HistsSysReweightsPS[ch][3][20][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 10) HistsSysReweightsPS[ch][3][22][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysReweightsPS[ch][3][21][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
         if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysReweightsPS[ch][3][23][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB * (SGW[0]/SGW[n])*((*gen_weight_sys)[n]/(*gen_weight_sys)[0]));
      }
    }
  }
  if(nbjetJesUp >1) {
    HistsSysUp[ch][3][0][9]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    HistsSysUp[ch][3][1][9]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    HistsSysUp[ch][3][2][9]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    HistsSysUp[ch][3][3][9]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    HistsSysUp[ch][3][4][9]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    HistsSysUp[ch][3][5][9]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    HistsSysUp[ch][3][6][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysUp[ch][3][7][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB);
    HistsSysUp[ch][3][8][9]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysUp[ch][3][9][9]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysUp[ch][3][16][9]->Fill(MET_FinalCollection_phi,weight_lepB);
    HistsSysUp[ch][3][17][9]->Fill(pv_n,weight_lepB);
    HistsSysUp[ch][3][18][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysUp[ch][3][10][9]->Fill((*selectedJetsJesUp)[0]->pt_,weight_lepB );
    HistsSysUp[ch][3][11][9]->Fill((*selectedJetsJesUp)[0]->eta_,weight_lepB);
    HistsSysUp[ch][3][12][9]->Fill((*selectedJetsJesUp)[0]->phi_,weight_lepB);
    HistsSysUp[ch][3][13][9]->Fill(selectedJetsJesUp->size(),weight_lepB );
    HistsSysUp[ch][3][14][9]->Fill(nbjetJesUp,weight_lepB);
    HistsSysUp[ch][3][15][9]->Fill(MET_T1SmearJetEnUp_Pt,weight_lepB );
    HistsSysUp[ch][3][19][9]->Fill(MVAoutputJesUp,weight_lepB );
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysUp[ch][3][20][9]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysUp[ch][3][22][9]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysUp[ch][3][21][9]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysUp[ch][3][23][9]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysUp[ch][3][20][9]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysUp[ch][3][22][9]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysUp[ch][3][21][9]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysUp[ch][3][23][9]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
  }
  if(nbjetJesDown>1){
    HistsSysDown[ch][3][0][9]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    HistsSysDown[ch][3][1][9]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    HistsSysDown[ch][3][2][9]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    HistsSysDown[ch][3][3][9]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    HistsSysDown[ch][3][4][9]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    HistsSysDown[ch][3][5][9]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    HistsSysDown[ch][3][6][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysDown[ch][3][7][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB);
    HistsSysDown[ch][3][8][9]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysDown[ch][3][9][9]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysDown[ch][3][16][9]->Fill(MET_FinalCollection_phi,weight_lepB);
    HistsSysDown[ch][3][17][9]->Fill(pv_n,weight_lepB);
    HistsSysDown[ch][3][18][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysDown[ch][3][10][9]->Fill((*selectedJetsJesDown)[0]->pt_,weight_lepB );
    HistsSysDown[ch][3][11][9]->Fill((*selectedJetsJesDown)[0]->eta_,weight_lepB);
    HistsSysDown[ch][3][12][9]->Fill((*selectedJetsJesDown)[0]->phi_,weight_lepB);
    HistsSysDown[ch][3][13][9]->Fill(selectedJetsJesDown->size(),weight_lepB );
    HistsSysDown[ch][3][14][9]->Fill(nbjetJesDown,weight_lepB);
    HistsSysDown[ch][3][15][9]->Fill(MET_T1SmearJetEnDown_Pt,weight_lepB );
    HistsSysDown[ch][3][19][9]->Fill(MVAoutputJesDown,weight_lepB );
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysDown[ch][3][20][9]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysDown[ch][3][22][9]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysDown[ch][3][21][9]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysDown[ch][3][23][9]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysDown[ch][3][20][9]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysDown[ch][3][22][9]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysDown[ch][3][21][9]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysDown[ch][3][23][9]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
  }
  if(nbjetJerUp>1){
    HistsSysUp[ch][3][0][10]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    HistsSysUp[ch][3][1][10]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    HistsSysUp[ch][3][2][10]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    HistsSysUp[ch][3][3][10]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    HistsSysUp[ch][3][4][10]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    HistsSysUp[ch][3][5][10]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    HistsSysUp[ch][3][6][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysUp[ch][3][7][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB);
    HistsSysUp[ch][3][8][10]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysUp[ch][3][9][10]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysUp[ch][3][16][10]->Fill(MET_FinalCollection_phi,weight_lepB);
    HistsSysUp[ch][3][17][10]->Fill(pv_n,weight_lepB);
    HistsSysUp[ch][3][18][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysUp[ch][3][10][10]->Fill((*selectedJetsJerUp)[0]->pt_,weight_lepB );
    HistsSysUp[ch][3][11][10]->Fill((*selectedJetsJerUp)[0]->eta_,weight_lepB);
    HistsSysUp[ch][3][12][10]->Fill((*selectedJetsJerUp)[0]->phi_,weight_lepB);
    HistsSysUp[ch][3][13][10]->Fill(selectedJetsJerUp->size(),weight_lepB );
    HistsSysUp[ch][3][14][10]->Fill(nbjetJerUp,weight_lepB);
    HistsSysUp[ch][3][15][10]->Fill(MET_T1SmearJetResUp_Pt,weight_lepB );
    HistsSysUp[ch][3][19][10]->Fill(MVAoutputJerUp,weight_lepB );
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysUp[ch][3][20][10]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysUp[ch][3][22][10]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysUp[ch][3][21][10]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysUp[ch][3][23][10]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysUp[ch][3][20][10]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysUp[ch][3][22][10]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysUp[ch][3][21][10]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysUp[ch][3][23][10]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
  }

  if(nbjetJerDown>1){
    HistsSysDown[ch][3][0][10]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    HistsSysDown[ch][3][1][10]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    HistsSysDown[ch][3][2][10]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    HistsSysDown[ch][3][3][10]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    HistsSysDown[ch][3][4][10]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    HistsSysDown[ch][3][5][10]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    HistsSysDown[ch][3][6][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysDown[ch][3][7][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB);
    HistsSysDown[ch][3][8][10]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysDown[ch][3][9][10]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB);
    HistsSysDown[ch][3][16][10]->Fill(MET_FinalCollection_phi,weight_lepB);
    HistsSysDown[ch][3][17][10]->Fill(pv_n,weight_lepB);
    HistsSysDown[ch][3][18][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
    HistsSysDown[ch][3][10][10]->Fill((*selectedJetsJerDown)[0]->pt_,weight_lepB );
    HistsSysDown[ch][3][11][10]->Fill((*selectedJetsJerDown)[0]->eta_,weight_lepB);
    HistsSysDown[ch][3][12][10]->Fill((*selectedJetsJerDown)[0]->phi_,weight_lepB);
    HistsSysDown[ch][3][13][10]->Fill(selectedJetsJerDown->size(),weight_lepB );
    HistsSysDown[ch][3][14][10]->Fill(nbjetJerDown,weight_lepB);
    HistsSysDown[ch][3][15][10]->Fill(MET_T1SmearJetResDown_Pt,weight_lepB );
    HistsSysDown[ch][3][19][10]->Fill(MVAoutputJerDown,weight_lepB );
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysDown[ch][3][20][10]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 10) HistsSysDown[ch][3][22][10]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysDown[ch][3][21][10]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[0]->lep_ == 1)  HistsSysDown[ch][3][23][10]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysDown[ch][3][20][10]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 10) HistsSysDown[ch][3][22][10]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysDown[ch][3][21][10]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
    if ((*selectedLeptons)[1]->lep_ == 1)  HistsSysDown[ch][3][23][10]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
  }
  for (int n=0;n<sysJecNames.size();++n){
    if ((*JECsysNbtagUp)[n]>1) {
      HistsJECUp[ch][3][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
      HistsJECUp[ch][3][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
      HistsJECUp[ch][3][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
      HistsJECUp[ch][3][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      HistsJECUp[ch][3][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      HistsJECUp[ch][3][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
      HistsJECUp[ch][3][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
      HistsJECUp[ch][3][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB);
      HistsJECUp[ch][3][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB);
      HistsJECUp[ch][3][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB);
      HistsJECUp[ch][3][16][n]->Fill(MET_FinalCollection_phi,weight_lepB);
      HistsJECUp[ch][3][17][n]->Fill(pv_n,weight_lepB);
      HistsJECUp[ch][3][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
      HistsJECUp[ch][3][10][n]->Fill((*JECsysUp)[n][0]->pt_,weight_lepB );
      HistsJECUp[ch][3][11][n]->Fill((*JECsysUp)[n][0]->eta_,weight_lepB);
      HistsJECUp[ch][3][12][n]->Fill((*JECsysUp)[n][0]->phi_,weight_lepB);
      HistsJECUp[ch][3][13][n]->Fill((*JECsysUp)[n].size(),weight_lepB );
      HistsJECUp[ch][3][14][n]->Fill((*JECsysNbtagUp)[n],weight_lepB);
      HistsJECUp[ch][3][15][n]->Fill((*JECsysMETUp)[n],weight_lepB );
      HistsJECUp[ch][3][19][n]->Fill((*JECsysMVAUp)[n],weight_lepB );
      if ((*selectedLeptons)[0]->lep_ == 10) HistsJECUp[ch][3][20][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[0]->lep_ == 10) HistsJECUp[ch][3][22][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsJECUp[ch][3][21][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsJECUp[ch][3][23][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 10) HistsJECUp[ch][3][20][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 10) HistsJECUp[ch][3][22][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsJECUp[ch][3][21][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsJECUp[ch][3][23][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
      }
    if ((*JECsysNbtagDown)[n]>1) {
      HistsJECDown[ch][3][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
      HistsJECDown[ch][3][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
      HistsJECDown[ch][3][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
      HistsJECDown[ch][3][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      HistsJECDown[ch][3][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      HistsJECDown[ch][3][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
      HistsJECDown[ch][3][6][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
      HistsJECDown[ch][3][7][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight_lepB);
      HistsJECDown[ch][3][8][n]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight_lepB);
      HistsJECDown[ch][3][9][n]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight_lepB);
      HistsJECDown[ch][3][16][n]->Fill(MET_FinalCollection_phi,weight_lepB);
      HistsJECDown[ch][3][17][n]->Fill(pv_n,weight_lepB);
      HistsJECDown[ch][3][18][n]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight_lepB);
      HistsJECDown[ch][3][10][n]->Fill((*JECsysDown)[n][0]->pt_,weight_lepB );
      HistsJECDown[ch][3][11][n]->Fill((*JECsysDown)[n][0]->eta_,weight_lepB);
      HistsJECDown[ch][3][12][n]->Fill((*JECsysDown)[n][0]->phi_,weight_lepB);
      HistsJECDown[ch][3][13][n]->Fill((*JECsysDown)[n].size(),weight_lepB );
      HistsJECDown[ch][3][14][n]->Fill((*JECsysNbtagDown)[n],weight_lepB);
      HistsJECDown[ch][3][15][n]->Fill((*JECsysMETDown)[n],weight_lepB );
      HistsJECDown[ch][3][19][n]->Fill((*JECsysMVADown)[n],weight_lepB );
      if ((*selectedLeptons)[0]->lep_ == 10) HistsJECDown[ch][3][20][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[0]->lep_ == 10) HistsJECDown[ch][3][22][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsJECDown[ch][3][21][n]->Fill((*selectedLeptons)[0]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[0]->lep_ == 1)  HistsJECDown[ch][3][23][n]->Fill((*selectedLeptons)[0]->eta_  ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 10) HistsJECDown[ch][3][20][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 10) HistsJECDown[ch][3][22][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsJECDown[ch][3][21][n]->Fill((*selectedLeptons)[1]->pt_   ,weight_lepB);
      if ((*selectedLeptons)[1]->lep_ == 1)  HistsJECDown[ch][3][23][n]->Fill((*selectedLeptons)[1]->eta_  ,weight_lepB);
    }
  }

//clean the memory
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
    JECsysMVAUp->clear();
    JECsysMVADown->clear();
    JECsysMETUp->shrink_to_fit();
    JECsysMETDown->shrink_to_fit();
    JECsysMVAUp->shrink_to_fit();
    JECsysMVADown->shrink_to_fit();
    delete JECsysUp;
    delete JECsysDown;
    delete JECsysNbtagUp;
    delete JECsysNbtagDown;
    delete JECsysMETUp;
    delete JECsysMETDown;
    delete JECsysMVAUp;
    delete JECsysMVADown;
    nAccept++;
  } //end of event loop
  cout<<"from "<<ntr<<" events, "<<nAccept<<" events are accepted"<<endl;
  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int l=0;l<vars.size();++l){
        Hists[i][k][l]  ->Write("",TObject::kOverwrite);
        for (int n=0;n<sys.size();++n){
          HistsSysUp[i][k][l][n]->Write("",TObject::kOverwrite);
          HistsSysDown[i][k][l][n]->Write("",TObject::kOverwrite);
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

  file_out.mkdir("JECSys");
  file_out.cd("JECSys/");
  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int l=0;l<vars.size();++l){
        for (int n=0;n<sysJecNames.size();++n){
          if( channels[i] != "emu" || k<2) continue;
          HistsJECUp[i][k][l][n]->Write("",TObject::kOverwrite);
          HistsJECDown[i][k][l][n]->Write("",TObject::kOverwrite);
        }
      }
    }
  }

  file_out.cd("");  

  if (fname.Contains("TTTo2L2Nu") || fname.Contains("LFV")){
    file_out.mkdir("reweightingSys");
    file_out.cd("reweightingSys/");
    for (int i=0;i<channels.size();++i){
      for (int k=0;k<regions.size();++k){
        for (int l=0;l<vars.size();++l){
          if( channels[i] != "emu" || k<2) continue;
          for (int n=0;n<reweightSizeQscalePDF;++n){
            HistsSysReweightsQscalePDF[i][k][l][n]->Write("",TObject::kOverwrite);
          }
          for (int n=0;n<reweightSizePS;++n){
            HistsSysReweightsPS[i][k][l][n]->Write("",TObject::kOverwrite);
          }
        }
      }
    }
  }

  file_out.cd("");
  file_out.Close() ;
cout<<"Cleaning the memory"<<endl;

Hists.clear();
cout<<"1"<<endl;
HistsSysUp.clear();
cout<<"1"<<endl;
HistsSysDown.clear();
cout<<"1"<<endl;
HistsJECUp.clear();
cout<<"1"<<endl;
HistsJECDown.clear();
cout<<"1"<<endl;
HistsSysReweightsQscalePDF.clear();
cout<<"1"<<endl;
HistsSysReweightsPS.clear();

  cout<<"Job is finished"<<endl;
}
