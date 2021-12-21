#define MyAnalysis_cxx
#include "MyAnalysis.h"
#include "PU_reWeighting.h"
#include "lepton_candidate.h"
#include "jet_candidate.h"
#include "TRandom.h"
#include "TRandom3.h"
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
#include "WCPoint.h"
#include "WCFit.h"
#include "TH1EFT.h"
#include "utilities.h"
#include <memory>
#include <TLorentzVector.h>


void displayProgress(long current, long max){
  using std::cerr;
  if (max<500) return;
  if (current%(max/500)!=0 && current<max-1) return;

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

bool ComparePtLep(TLorentzVector a, TLorentzVector b) { return a.Pt() > b.Pt(); }
bool ComparePtJet(jet_candidate *a, jet_candidate *b) { return a->pt_ > b->pt_; }

float scale_factor( TH2F* h, float X, float Y , TString uncert){
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
  if(uncert=="up") return (h->GetBinContent(binx, biny)+h->GetBinError(binx, biny));
  else if(uncert=="down") return (h->GetBinContent(binx, biny)-h->GetBinError(binx, biny));
  else return  h->GetBinContent(binx, biny);
}

float topPt(float pt){
  return (0.973 - (0.000134 * pt) + (0.103 * exp(pt * (-0.0118))));  
}

void MyAnalysis::Loop(TString fname, TString data, TString dataset ,TString year, TString run, float xs, float lumi, float Nevent, int iseft, int nRuns)
{

  TFile *f = new TFile("ANoutput.root","RECREATE");
  f->cd();

  typedef vector<TH1EFT*> Dim1;
  typedef vector<Dim1> Dim2;
  typedef vector<Dim2> Dim3;
  typedef vector<Dim3> Dim4;

  std::vector<TString> vars   {"All","l1Pt","l2Pt","l1Eta","l2Eta", "Mll","Drll","Dphill","Njet","jet1pt","jet1eta", "MET"};
  std::vector<int>    nbins   {1    ,25    ,25    ,20     ,20     , 30   ,20    ,20      ,10    ,25      ,20      ,20};
  std::vector<float> lowEdge  {0    ,0     ,0     ,-3     ,-3     ,0     ,0     ,0       ,0     ,0       ,-3      ,0};
  std::vector<float> highEdge {1    ,500   ,500   ,3      ,3      ,500   ,7     ,7       ,10    ,500     ,3       ,400};
  int nn=0;
  float mtW=0;
  float sOw=0;
 
  std::stringstream name;
  TH1EFT *h_test;
  Dim1 HistsSysReweightsEFT(vars.size());
  for (int l=0;l<vars.size();++l){
    name<<vars[l];
    h_test = new TH1EFT((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]);
    HistsSysReweightsEFT[l] = h_test;
    name.str("");
  }

  TLorentzVector nu, ell, TLV;
  std::vector<TLorentzVector> *selectedlep, *selectedjet;
  float weight_Lumi;
  float weight;

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t ntr = fChain->GetEntries ();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<ntr;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    displayProgress(jentry, ntr) ;
    if (data == "mc") weight_Lumi = (1000*xs*lumi)/Nevent;
    if (data == "mc" && iseft) weight_Lumi = (1000*lumi)/nRuns;

// Add EFT weights
//    sOw +=LHE_weight_nominal;
    std::vector<WCPoint> wc_pts;

    if (iseft) {
        for (int s=0;s<mc_LHEweightsId->size();++s) {
            if((*mc_LHEweightsId)[s]=="1001") sOw +=(*LHE_weight_sys)[s]; //cout<<(*mc_LHEweightsId)[s]<<":"<<(*LHE_weight_sys)[s]<<endl;
            auto LHEwgtstr = std::string((*mc_LHEweightsId)[s]);
            std::size_t foundstr = LHEwgtstr.find("EFTrwgt"); // only save our EFT weights
            if (foundstr!=std::string::npos) {
                weight = weight_Lumi * (*LHE_weight_sys)[s];
                WCPoint wc_pt((*mc_LHEweightsId)[s],weight);
                wc_pts.push_back(wc_pt);
            }
        }
      weight=1;
    } else {
        WCPoint wc_pt("smpt",1);
        wc_pts.push_back(wc_pt);
        weight=weight_Lumi;
    }

    WCFit eft_fit(wc_pts,"");
    HistsSysReweightsEFT[0]->Fill(0.5,weight,eft_fit);
//lepton candidates
    selectedlep = new std::vector<TLorentzVector>();
    for (int l=0;l<pl_lep_pt->size();l++){
//      if(abs((*pl_lep_pt)[l]) <25 || abs((*pl_lep_eta)[l]) > 2.4 ) continue;
      TLV.SetPtEtaPhiM((*pl_lep_pt)[l], (*pl_lep_eta)[l], (*pl_lep_phi)[l], 0) ;
      selectedlep->push_back(TLV);
    }
    sort(selectedlep->begin(), selectedlep->end(), ComparePtLep);
//dilepton selection
    if(selectedlep->size() < 2) continue;
    nn+=1;
//Jet selection
    selectedjet = new std::vector<TLorentzVector>();
    for (int l=0;l<pl_jet_pt->size();l++){
      if(abs((*pl_jet_pt)[l]) <30 || abs((*pl_jet_eta)[l]) > 2.4 ) continue;
      TLV.SetPtEtaPhiM((*pl_lep_pt)[l], (*pl_lep_eta)[l], (*pl_lep_phi)[l], 0) ;
      selectedjet->push_back(TLV);
    }
    sort(selectedjet->begin(), selectedjet->end(), ComparePtLep);
    nu.SetPtEtaPhiM((*pl_MET_pt)[0], 0,(*pl_MET_phi)[0], 0) ;
    ell.SetPtEtaPhiM(selectedlep->at(0).Pt(), 0,selectedlep->at(0).Phi(), 0) ;
    HistsSysReweightsEFT[1]->Fill(selectedlep->at(0).Pt(),weight,eft_fit);
    HistsSysReweightsEFT[2]->Fill(selectedlep->at(1).Pt(),weight,eft_fit);
    HistsSysReweightsEFT[3]->Fill(selectedlep->at(0).Eta(),weight,eft_fit);
    HistsSysReweightsEFT[4]->Fill(selectedlep->at(1).Eta(),weight,eft_fit);
    HistsSysReweightsEFT[5]->Fill((selectedlep->at(0)+selectedlep->at(1)).M(),weight,eft_fit);
    HistsSysReweightsEFT[6]->Fill(abs(deltaR(selectedlep->at(0).Eta(),selectedlep->at(0).Phi(),selectedlep->at(1).Eta(),selectedlep->at(1).Phi())),weight,eft_fit);
    HistsSysReweightsEFT[7]->Fill(abs(deltaPhi(selectedlep->at(0).Phi(),selectedlep->at(1).Phi())),weight,eft_fit);
    HistsSysReweightsEFT[8]->Fill(selectedjet->size(),weight,eft_fit);
    if(selectedjet->size()>0){
      HistsSysReweightsEFT[9]->Fill(selectedjet->at(0).Pt(),weight,eft_fit);
      HistsSysReweightsEFT[10]->Fill(selectedjet->at(0).Eta(),weight,eft_fit);
    }
    HistsSysReweightsEFT[11]->Fill((*pl_MET_pt)[0],weight,eft_fit);

  delete selectedlep;
  delete selectedjet;
  }
cout<<nn<<" events************** "<<sOw<<endl;
  for (int l=0;l<vars.size();++l){
    HistsSysReweightsEFT[l]->Write("",TObject::kOverwrite);
  }

  f->Close();

}


