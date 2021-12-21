//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jun 19 03:34:14 2021 by ROOT version 6.12/07
// from TTree IIHEAnalysis/IIHEAnalysis
// found on file: /hadoop/store/user/rgoldouz/FullProduction/tt_EFT/ntuple_2l2q_tt_lnubDecay/outfile_573.root
//////////////////////////////////////////////////////////

#ifndef MyAnalysis_h
#define MyAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
using namespace std;

class MyAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *mc_sumofLHEWeights;
   vector<string>  *mc_LHEweightsId;
   Float_t         mc_nEventsWeighted;
   vector<float>   *mc_sumofgenWeights;
   Float_t         nEventsRaw;
   Float_t         nEventsStored;
   vector<float>   *nRuns;


   vector<float>   *LHE_Pt;
   vector<float>   *LHE_Eta;
   vector<float>   *LHE_Phi;
   vector<float>   *LHE_E;
   vector<int>     *LHE_pdgid;
   vector<int>     *LHE_status;
   Float_t         LHE_weight_nominal;
   vector<float>   *LHE_weight_sys;
   UInt_t          mc_n;
   Int_t           mc_nMEPartons;
   Int_t           mc_nMEPartonsFiltered;
   vector<float>   *mc_DJRValues;
   Float_t         mc_weight;
   Float_t         mc_w_sign;
   Int_t           mc_id_first;
   Int_t           mc_id_second;
   Float_t         mc_x_first;
   Float_t         mc_x_second;
   Float_t         mc_xPDF_first;
   Float_t         mc_xPDF_second;
   Float_t         mc_scalePDF;
   vector<int>     *mc_index;
   vector<int>     *mc_pdgId;
   vector<int>     *mc_charge;
   vector<int>     *mc_status;
   vector<int>     *mc_status_flags;
   vector<int>     *mc_status_tau_flags;
   vector<int>     *mc_tau_charge;
   vector<int>     *mc_tau_pdgId;
   vector<int>     *mc_tau_decay;
   vector<int>     *mc_tau_had_status;
   vector<int>     *mc_tau_had_charge;
   vector<int>     *mc_tau_had_pdgId;
   vector<float>   *mc_mass;
   vector<float>   *mc_px;
   vector<float>   *mc_py;
   vector<float>   *mc_pz;
   vector<float>   *mc_pt;
   vector<float>   *mc_eta;
   vector<float>   *mc_phi;
   vector<float>   *mc_energy;
   vector<float>   *mc_tau_pt;
   vector<float>   *mc_tau_eta;
   vector<float>   *mc_tau_phi;
   vector<float>   *mc_tau_energy;
   vector<float>   *mc_tau_had_pt;
   vector<float>   *mc_tau_had_eta;
   vector<float>   *mc_tau_had_phi;
   vector<float>   *mc_tau_had_energy;
   vector<unsigned int> *mc_numberOfDaughters;
   vector<unsigned int> *mc_numberOfMothers;
   vector<vector<int> > *mc_mother_index;
   vector<vector<int> > *mc_mother_pdgId;
   vector<vector<float> > *mc_mother_px;
   vector<vector<float> > *mc_mother_py;
   vector<vector<float> > *mc_mother_pz;
   vector<vector<float> > *mc_mother_pt;
   vector<vector<float> > *mc_mother_eta;
   vector<vector<float> > *mc_mother_phi;
   vector<vector<float> > *mc_mother_energy;
   vector<vector<float> > *mc_mother_mass;
   Int_t           mc_trueNumInteractions;
   Int_t           mc_PU_NumInteractions;
   vector<float>   *genjet_pt;
   vector<float>   *genjet_eta;
   vector<float>   *genjet_phi;
   vector<float>   *genjet_energy;
   vector<float>   *genjetAK8_pt;
   vector<float>   *genjetAK8_eta;
   vector<float>   *genjetAK8_phi;
   vector<float>   *genjetAK8_energy;
   vector<float>   *gen_weight_sys;
   vector<float>   *pl_jet_pt;
   vector<float>   *pl_jet_eta;
   vector<float>   *pl_jet_phi;
   vector<float>   *pl_lep_pt;
   vector<float>   *pl_lep_eta;
   vector<float>   *pl_lep_phi;
   vector<float>   *pl_ph_pt;
   vector<float>   *pl_ph_eta;
   vector<float>   *pl_ph_phi;
   vector<float>   *pl_MET_pt;
   vector<float>   *pl_MET_phi;
   vector<int>     *pl_lep_pdgid;
   vector<int>     *pl_lep_charge;
   vector<int>     *pl_jet_pdgid;

   // List of branches
   TBranch        *b_mc_sumofLHEWeights;   //!
   TBranch        *b_mc_LHEweightsId;   //!
   TBranch        *b_mc_nEventsWeighted;   //!
   TBranch        *b_mc_sumofgenWeights;   //!
   TBranch        *b_nEventsRaw;   //!
   TBranch        *b_nEventsStored;   //!
   TBranch        *b_nRuns;   //!

   TBranch        *b_LHE_Pt;   //!
   TBranch        *b_LHE_Eta;   //!
   TBranch        *b_LHE_Phi;   //!
   TBranch        *b_LHE_E;   //!
   TBranch        *b_LHE_pdgid;   //!
   TBranch        *b_LHE_status;   //!
   TBranch        *b_LHE_weight_nominal;   //!
   TBranch        *b_LHE_weight_sys;   //!
   TBranch        *b_mc_n;   //!
   TBranch        *b_mc_nMEPartons;   //!
   TBranch        *b_mc_nMEPartonsFiltered;   //!
   TBranch        *b_mc_DJRValues;   //!
   TBranch        *b_mc_weight;   //!
   TBranch        *b_mc_w_sign;   //!
   TBranch        *b_mc_id_first;   //!
   TBranch        *b_mc_id_second;   //!
   TBranch        *b_mc_x_first;   //!
   TBranch        *b_mc_x_second;   //!
   TBranch        *b_mc_xPDF_first;   //!
   TBranch        *b_mc_xPDF_second;   //!
   TBranch        *b_mc_scalePDF;   //!
   TBranch        *b_mc_index;   //!
   TBranch        *b_mc_pdgId;   //!
   TBranch        *b_mc_charge;   //!
   TBranch        *b_mc_status;   //!
   TBranch        *b_mc_status_flags;   //!
   TBranch        *b_mc_status_tau_flags;   //!
   TBranch        *b_mc_tau_charge;   //!
   TBranch        *b_mc_tau_pdgId;   //!
   TBranch        *b_mc_tau_decay;   //!
   TBranch        *b_mc_tau_had_status;   //!
   TBranch        *b_mc_tau_had_charge;   //!
   TBranch        *b_mc_tau_had_pdgId;   //!
   TBranch        *b_mc_mass;   //!
   TBranch        *b_mc_px;   //!
   TBranch        *b_mc_py;   //!
   TBranch        *b_mc_pz;   //!
   TBranch        *b_mc_pt;   //!
   TBranch        *b_mc_eta;   //!
   TBranch        *b_mc_phi;   //!
   TBranch        *b_mc_energy;   //!
   TBranch        *b_mc_tau_pt;   //!
   TBranch        *b_mc_tau_eta;   //!
   TBranch        *b_mc_tau_phi;   //!
   TBranch        *b_mc_tau_energy;   //!
   TBranch        *b_mc_tau_had_pt;   //!
   TBranch        *b_mc_tau_had_eta;   //!
   TBranch        *b_mc_tau_had_phi;   //!
   TBranch        *b_mc_tau_had_energy;   //!
   TBranch        *b_mc_numberOfDaughters;   //!
   TBranch        *b_mc_numberOfMothers;   //!
   TBranch        *b_mc_mother_index;   //!
   TBranch        *b_mc_mother_pdgId;   //!
   TBranch        *b_mc_mother_px;   //!
   TBranch        *b_mc_mother_py;   //!
   TBranch        *b_mc_mother_pz;   //!
   TBranch        *b_mc_mother_pt;   //!
   TBranch        *b_mc_mother_eta;   //!
   TBranch        *b_mc_mother_phi;   //!
   TBranch        *b_mc_mother_energy;   //!
   TBranch        *b_mc_mother_mass;   //!
   TBranch        *b_mc_trueNumInteractions;   //!
   TBranch        *b_mc_PU_NumInteractions;   //!
   TBranch        *b_genjet_pt;   //!
   TBranch        *b_genjet_eta;   //!
   TBranch        *b_genjet_phi;   //!
   TBranch        *b_genjet_energy;   //!
   TBranch        *b_genjetAK8_pt;   //!
   TBranch        *b_genjetAK8_eta;   //!
   TBranch        *b_genjetAK8_phi;   //!
   TBranch        *b_genjetAK8_energy;   //!
   TBranch        *b_gen_weight_sys;   //!
   TBranch        *b_pl_jet_pt;   //!
   TBranch        *b_pl_jet_eta;   //!
   TBranch        *b_pl_jet_phi;   //!
   TBranch        *b_pl_lep_pt;   //!
   TBranch        *b_pl_lep_eta;   //!
   TBranch        *b_pl_lep_phi;   //!
   TBranch        *b_pl_ph_pt;   //!
   TBranch        *b_pl_ph_eta;   //!
   TBranch        *b_pl_ph_phi;   //!
   TBranch        *b_pl_MET_pt;   //!
   TBranch        *b_pl_MET_phi;   //!
   TBranch        *b_pl_lep_pdgid;   //!
   TBranch        *b_pl_lep_charge;   //!
   TBranch        *b_pl_jet_pdgid;   //!

   MyAnalysis(TTree *tree=0);
   virtual ~MyAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString, TString, TString, TString, TString, float,float,float,int,int);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyAnalysis_cxx
MyAnalysis::MyAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/hadoop/store/user/rgoldouz/FullProduction/tt_EFT/ntuple_2l2q_tt_lnubDecay/outfile_573.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/hadoop/store/user/rgoldouz/FullProduction/tt_EFT/ntuple_2l2q_tt_lnubDecay/outfile_573.root");
      }
      f->GetObject("IIHEAnalysis",tree);

   }
   Init(tree);
}

MyAnalysis::~MyAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   mc_sumofLHEWeights = 0;
   mc_LHEweightsId = 0;
   mc_sumofgenWeights = 0;
   nRuns = 0;
   LHE_Pt = 0;
   LHE_Eta = 0;
   LHE_Phi = 0;
   LHE_E = 0;
   LHE_pdgid = 0;
   LHE_status = 0;
   LHE_weight_sys = 0;
   mc_DJRValues = 0;
   mc_index = 0;
   mc_pdgId = 0;
   mc_charge = 0;
   mc_status = 0;
   mc_status_flags = 0;
   mc_status_tau_flags = 0;
   mc_tau_charge = 0;
   mc_tau_pdgId = 0;
   mc_tau_decay = 0;
   mc_tau_had_status = 0;
   mc_tau_had_charge = 0;
   mc_tau_had_pdgId = 0;
   mc_mass = 0;
   mc_px = 0;
   mc_py = 0;
   mc_pz = 0;
   mc_pt = 0;
   mc_eta = 0;
   mc_phi = 0;
   mc_energy = 0;
   mc_tau_pt = 0;
   mc_tau_eta = 0;
   mc_tau_phi = 0;
   mc_tau_energy = 0;
   mc_tau_had_pt = 0;
   mc_tau_had_eta = 0;
   mc_tau_had_phi = 0;
   mc_tau_had_energy = 0;
   mc_numberOfDaughters = 0;
   mc_numberOfMothers = 0;
   mc_mother_index = 0;
   mc_mother_pdgId = 0;
   mc_mother_px = 0;
   mc_mother_py = 0;
   mc_mother_pz = 0;
   mc_mother_pt = 0;
   mc_mother_eta = 0;
   mc_mother_phi = 0;
   mc_mother_energy = 0;
   mc_mother_mass = 0;
   genjet_pt = 0;
   genjet_eta = 0;
   genjet_phi = 0;
   genjet_energy = 0;
   genjetAK8_pt = 0;
   genjetAK8_eta = 0;
   genjetAK8_phi = 0;
   genjetAK8_energy = 0;
   gen_weight_sys = 0;
   pl_jet_pt = 0;
   pl_jet_eta = 0;
   pl_jet_phi = 0;
   pl_lep_pt = 0;
   pl_lep_eta = 0;
   pl_lep_phi = 0;
   pl_ph_pt = 0;
   pl_ph_eta = 0;
   pl_ph_phi = 0;
   pl_MET_pt = 0;
   pl_MET_phi = 0;
   pl_lep_pdgid = 0;
   pl_lep_charge = 0;
   pl_jet_pdgid = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mc_sumofLHEWeights", &mc_sumofLHEWeights, &b_mc_sumofLHEWeights);
   fChain->SetBranchAddress("mc_LHEweightsId", &mc_LHEweightsId, &b_mc_LHEweightsId);
   fChain->SetBranchAddress("mc_nEventsWeighted", &mc_nEventsWeighted, &b_mc_nEventsWeighted);
   fChain->SetBranchAddress("mc_sumofgenWeights", &mc_sumofgenWeights, &b_mc_sumofgenWeights);
   fChain->SetBranchAddress("nEventsRaw", &nEventsRaw, &b_nEventsRaw);
   fChain->SetBranchAddress("nEventsStored", &nEventsStored, &b_nEventsStored);
   fChain->SetBranchAddress("nRuns", &nRuns, &b_nRuns);

   fChain->SetBranchAddress("LHE_Pt", &LHE_Pt, &b_LHE_Pt);
   fChain->SetBranchAddress("LHE_Eta", &LHE_Eta, &b_LHE_Eta);
   fChain->SetBranchAddress("LHE_Phi", &LHE_Phi, &b_LHE_Phi);
   fChain->SetBranchAddress("LHE_E", &LHE_E, &b_LHE_E);
   fChain->SetBranchAddress("LHE_pdgid", &LHE_pdgid, &b_LHE_pdgid);
   fChain->SetBranchAddress("LHE_status", &LHE_status, &b_LHE_status);
   fChain->SetBranchAddress("LHE_weight_nominal", &LHE_weight_nominal, &b_LHE_weight_nominal);
   fChain->SetBranchAddress("LHE_weight_sys", &LHE_weight_sys, &b_LHE_weight_sys);
   fChain->SetBranchAddress("mc_n", &mc_n, &b_mc_n);
   fChain->SetBranchAddress("mc_nMEPartons", &mc_nMEPartons, &b_mc_nMEPartons);
   fChain->SetBranchAddress("mc_nMEPartonsFiltered", &mc_nMEPartonsFiltered, &b_mc_nMEPartonsFiltered);
   fChain->SetBranchAddress("mc_DJRValues", &mc_DJRValues, &b_mc_DJRValues);
   fChain->SetBranchAddress("mc_weight", &mc_weight, &b_mc_weight);
   fChain->SetBranchAddress("mc_w_sign", &mc_w_sign, &b_mc_w_sign);
   fChain->SetBranchAddress("mc_id_first", &mc_id_first, &b_mc_id_first);
   fChain->SetBranchAddress("mc_id_second", &mc_id_second, &b_mc_id_second);
   fChain->SetBranchAddress("mc_x_first", &mc_x_first, &b_mc_x_first);
   fChain->SetBranchAddress("mc_x_second", &mc_x_second, &b_mc_x_second);
   fChain->SetBranchAddress("mc_xPDF_first", &mc_xPDF_first, &b_mc_xPDF_first);
   fChain->SetBranchAddress("mc_xPDF_second", &mc_xPDF_second, &b_mc_xPDF_second);
   fChain->SetBranchAddress("mc_scalePDF", &mc_scalePDF, &b_mc_scalePDF);
   fChain->SetBranchAddress("mc_index", &mc_index, &b_mc_index);
   fChain->SetBranchAddress("mc_pdgId", &mc_pdgId, &b_mc_pdgId);
   fChain->SetBranchAddress("mc_charge", &mc_charge, &b_mc_charge);
   fChain->SetBranchAddress("mc_status", &mc_status, &b_mc_status);
   fChain->SetBranchAddress("mc_status_flags", &mc_status_flags, &b_mc_status_flags);
   fChain->SetBranchAddress("mc_status_tau_flags", &mc_status_tau_flags, &b_mc_status_tau_flags);
   fChain->SetBranchAddress("mc_tau_charge", &mc_tau_charge, &b_mc_tau_charge);
   fChain->SetBranchAddress("mc_tau_pdgId", &mc_tau_pdgId, &b_mc_tau_pdgId);
   fChain->SetBranchAddress("mc_tau_decay", &mc_tau_decay, &b_mc_tau_decay);
   fChain->SetBranchAddress("mc_tau_had_status", &mc_tau_had_status, &b_mc_tau_had_status);
   fChain->SetBranchAddress("mc_tau_had_charge", &mc_tau_had_charge, &b_mc_tau_had_charge);
   fChain->SetBranchAddress("mc_tau_had_pdgId", &mc_tau_had_pdgId, &b_mc_tau_had_pdgId);
   fChain->SetBranchAddress("mc_mass", &mc_mass, &b_mc_mass);
   fChain->SetBranchAddress("mc_px", &mc_px, &b_mc_px);
   fChain->SetBranchAddress("mc_py", &mc_py, &b_mc_py);
   fChain->SetBranchAddress("mc_pz", &mc_pz, &b_mc_pz);
   fChain->SetBranchAddress("mc_pt", &mc_pt, &b_mc_pt);
   fChain->SetBranchAddress("mc_eta", &mc_eta, &b_mc_eta);
   fChain->SetBranchAddress("mc_phi", &mc_phi, &b_mc_phi);
   fChain->SetBranchAddress("mc_energy", &mc_energy, &b_mc_energy);
   fChain->SetBranchAddress("mc_tau_pt", &mc_tau_pt, &b_mc_tau_pt);
   fChain->SetBranchAddress("mc_tau_eta", &mc_tau_eta, &b_mc_tau_eta);
   fChain->SetBranchAddress("mc_tau_phi", &mc_tau_phi, &b_mc_tau_phi);
   fChain->SetBranchAddress("mc_tau_energy", &mc_tau_energy, &b_mc_tau_energy);
   fChain->SetBranchAddress("mc_tau_had_pt", &mc_tau_had_pt, &b_mc_tau_had_pt);
   fChain->SetBranchAddress("mc_tau_had_eta", &mc_tau_had_eta, &b_mc_tau_had_eta);
   fChain->SetBranchAddress("mc_tau_had_phi", &mc_tau_had_phi, &b_mc_tau_had_phi);
   fChain->SetBranchAddress("mc_tau_had_energy", &mc_tau_had_energy, &b_mc_tau_had_energy);
   fChain->SetBranchAddress("mc_numberOfDaughters", &mc_numberOfDaughters, &b_mc_numberOfDaughters);
   fChain->SetBranchAddress("mc_numberOfMothers", &mc_numberOfMothers, &b_mc_numberOfMothers);
   fChain->SetBranchAddress("mc_mother_index", &mc_mother_index, &b_mc_mother_index);
   fChain->SetBranchAddress("mc_mother_pdgId", &mc_mother_pdgId, &b_mc_mother_pdgId);
   fChain->SetBranchAddress("mc_mother_px", &mc_mother_px, &b_mc_mother_px);
   fChain->SetBranchAddress("mc_mother_py", &mc_mother_py, &b_mc_mother_py);
   fChain->SetBranchAddress("mc_mother_pz", &mc_mother_pz, &b_mc_mother_pz);
   fChain->SetBranchAddress("mc_mother_pt", &mc_mother_pt, &b_mc_mother_pt);
   fChain->SetBranchAddress("mc_mother_eta", &mc_mother_eta, &b_mc_mother_eta);
   fChain->SetBranchAddress("mc_mother_phi", &mc_mother_phi, &b_mc_mother_phi);
   fChain->SetBranchAddress("mc_mother_energy", &mc_mother_energy, &b_mc_mother_energy);
   fChain->SetBranchAddress("mc_mother_mass", &mc_mother_mass, &b_mc_mother_mass);
   fChain->SetBranchAddress("mc_trueNumInteractions", &mc_trueNumInteractions, &b_mc_trueNumInteractions);
   fChain->SetBranchAddress("mc_PU_NumInteractions", &mc_PU_NumInteractions, &b_mc_PU_NumInteractions);
   fChain->SetBranchAddress("genjet_pt", &genjet_pt, &b_genjet_pt);
   fChain->SetBranchAddress("genjet_eta", &genjet_eta, &b_genjet_eta);
   fChain->SetBranchAddress("genjet_phi", &genjet_phi, &b_genjet_phi);
   fChain->SetBranchAddress("genjet_energy", &genjet_energy, &b_genjet_energy);
   fChain->SetBranchAddress("genjetAK8_pt", &genjetAK8_pt, &b_genjetAK8_pt);
   fChain->SetBranchAddress("genjetAK8_eta", &genjetAK8_eta, &b_genjetAK8_eta);
   fChain->SetBranchAddress("genjetAK8_phi", &genjetAK8_phi, &b_genjetAK8_phi);
   fChain->SetBranchAddress("genjetAK8_energy", &genjetAK8_energy, &b_genjetAK8_energy);
   fChain->SetBranchAddress("gen_weight_sys", &gen_weight_sys, &b_gen_weight_sys);
   fChain->SetBranchAddress("pl_jet_pt", &pl_jet_pt, &b_pl_jet_pt);
   fChain->SetBranchAddress("pl_jet_eta", &pl_jet_eta, &b_pl_jet_eta);
   fChain->SetBranchAddress("pl_jet_phi", &pl_jet_phi, &b_pl_jet_phi);
   fChain->SetBranchAddress("pl_lep_pt", &pl_lep_pt, &b_pl_lep_pt);
   fChain->SetBranchAddress("pl_lep_eta", &pl_lep_eta, &b_pl_lep_eta);
   fChain->SetBranchAddress("pl_lep_phi", &pl_lep_phi, &b_pl_lep_phi);
   fChain->SetBranchAddress("pl_ph_pt", &pl_ph_pt, &b_pl_ph_pt);
   fChain->SetBranchAddress("pl_ph_eta", &pl_ph_eta, &b_pl_ph_eta);
   fChain->SetBranchAddress("pl_ph_phi", &pl_ph_phi, &b_pl_ph_phi);
   fChain->SetBranchAddress("pl_MET_pt", &pl_MET_pt, &b_pl_MET_pt);
   fChain->SetBranchAddress("pl_MET_phi", &pl_MET_phi, &b_pl_MET_phi);
   fChain->SetBranchAddress("pl_lep_pdgid", &pl_lep_pdgid, &b_pl_lep_pdgid);
   fChain->SetBranchAddress("pl_lep_charge", &pl_lep_charge, &b_pl_lep_charge);
   fChain->SetBranchAddress("pl_jet_pdgid", &pl_jet_pdgid, &b_pl_jet_pdgid);
   Notify();
}

Bool_t MyAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyAnalysis_cxx
