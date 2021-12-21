//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov 26 11:13:46 2021 by ROOT version 6.12/07
// from TTree Events/Events
// found on file: /hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_TTJets/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_UL17_TTJets/211106_221617/0000/tree_77.root
//////////////////////////////////////////////////////////

#ifndef MyAnalysis_h
#define MyAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
using namespace std;
class MyAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   Float_t         HTXS_Higgs_pt;
   Float_t         HTXS_Higgs_y;
   Int_t           HTXS_stage1_1_cat_pTjet25GeV;
   Int_t           HTXS_stage1_1_cat_pTjet30GeV;
   Int_t           HTXS_stage1_1_fine_cat_pTjet25GeV;
   Int_t           HTXS_stage1_1_fine_cat_pTjet30GeV;
   Int_t           HTXS_stage1_2_cat_pTjet25GeV;
   Int_t           HTXS_stage1_2_cat_pTjet30GeV;
   Int_t           HTXS_stage1_2_fine_cat_pTjet25GeV;
   Int_t           HTXS_stage1_2_fine_cat_pTjet30GeV;
   Int_t           HTXS_stage_0;
   Int_t           HTXS_stage_1_pTjet25;
   Int_t           HTXS_stage_1_pTjet30;
   UChar_t         HTXS_njets25;
   UChar_t         HTXS_njets30;
   UInt_t          nboostedTau;
   Float_t         boostedTau_chargedIso[5];   //[nboostedTau]
   Float_t         boostedTau_eta[5];   //[nboostedTau]
   Float_t         boostedTau_leadTkDeltaEta[5];   //[nboostedTau]
   Float_t         boostedTau_leadTkDeltaPhi[5];   //[nboostedTau]
   Float_t         boostedTau_leadTkPtOverTauPt[5];   //[nboostedTau]
   Float_t         boostedTau_mass[5];   //[nboostedTau]
   Float_t         boostedTau_neutralIso[5];   //[nboostedTau]
   Float_t         boostedTau_phi[5];   //[nboostedTau]
   Float_t         boostedTau_photonsOutsideSignalCone[5];   //[nboostedTau]
   Float_t         boostedTau_pt[5];   //[nboostedTau]
   Float_t         boostedTau_puCorr[5];   //[nboostedTau]
   Float_t         boostedTau_rawAntiEle2018[5];   //[nboostedTau]
   Float_t         boostedTau_rawIso[5];   //[nboostedTau]
   Float_t         boostedTau_rawIsodR03[5];   //[nboostedTau]
   Float_t         boostedTau_rawMVAnewDM2017v2[5];   //[nboostedTau]
   Float_t         boostedTau_rawMVAoldDM2017v2[5];   //[nboostedTau]
   Float_t         boostedTau_rawMVAoldDMdR032017v2[5];   //[nboostedTau]
   Int_t           boostedTau_charge[5];   //[nboostedTau]
   Int_t           boostedTau_decayMode[5];   //[nboostedTau]
   Int_t           boostedTau_jetIdx[5];   //[nboostedTau]
   Int_t           boostedTau_rawAntiEleCat2018[5];   //[nboostedTau]
   UChar_t         boostedTau_idAntiEle2018[5];   //[nboostedTau]
   UChar_t         boostedTau_idAntiMu[5];   //[nboostedTau]
   UChar_t         boostedTau_idMVAnewDM2017v2[5];   //[nboostedTau]
   UChar_t         boostedTau_idMVAoldDM2017v2[5];   //[nboostedTau]
   UChar_t         boostedTau_idMVAoldDMdR032017v2[5];   //[nboostedTau]
   Float_t         btagWeight_CSVV2;
   Float_t         btagWeight_DeepCSVB;
   Float_t         CaloMET_phi;
   Float_t         CaloMET_pt;
   Float_t         CaloMET_sumEt;
   Float_t         ChsMET_phi;
   Float_t         ChsMET_pt;
   Float_t         ChsMET_sumEt;
   UInt_t          nCorrT1METJet;
   Float_t         CorrT1METJet_area[27];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_eta[27];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_muonSubtrFactor[27];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_phi[27];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_rawPt[27];   //[nCorrT1METJet]
   Float_t         DeepMETResolutionTune_phi;
   Float_t         DeepMETResolutionTune_pt;
   Float_t         DeepMETResponseTune_phi;
   Float_t         DeepMETResponseTune_pt;
   UInt_t          nElectron;
   Float_t         Electron_dEscaleDown[9];   //[nElectron]
   Float_t         Electron_dEscaleUp[9];   //[nElectron]
   Float_t         Electron_dEsigmaDown[9];   //[nElectron]
   Float_t         Electron_dEsigmaUp[9];   //[nElectron]
   Float_t         Electron_deltaEtaSC[9];   //[nElectron]
   Float_t         Electron_dr03EcalRecHitSumEt[9];   //[nElectron]
   Float_t         Electron_dr03HcalDepth1TowerSumEt[9];   //[nElectron]
   Float_t         Electron_dr03TkSumPt[9];   //[nElectron]
   Float_t         Electron_dr03TkSumPtHEEP[9];   //[nElectron]
   Float_t         Electron_dxy[9];   //[nElectron]
   Float_t         Electron_dxyErr[9];   //[nElectron]
   Float_t         Electron_dz[9];   //[nElectron]
   Float_t         Electron_dzErr[9];   //[nElectron]
   Float_t         Electron_eCorr[9];   //[nElectron]
   Float_t         Electron_eInvMinusPInv[9];   //[nElectron]
   Float_t         Electron_energyErr[9];   //[nElectron]
   Float_t         Electron_eta[9];   //[nElectron]
   Float_t         Electron_hoe[9];   //[nElectron]
   Float_t         Electron_ip3d[9];   //[nElectron]
   Float_t         Electron_jetPtRelv2[9];   //[nElectron]
   Float_t         Electron_jetRelIso[9];   //[nElectron]
   Float_t         Electron_mass[9];   //[nElectron]
   Float_t         Electron_miniPFRelIso_all[9];   //[nElectron]
   Float_t         Electron_miniPFRelIso_chg[9];   //[nElectron]
   Float_t         Electron_mvaFall17V2Iso[9];   //[nElectron]
   Float_t         Electron_mvaFall17V2noIso[9];   //[nElectron]
   Float_t         Electron_pfRelIso03_all[9];   //[nElectron]
   Float_t         Electron_pfRelIso03_chg[9];   //[nElectron]
   Float_t         Electron_phi[9];   //[nElectron]
   Float_t         Electron_pt[9];   //[nElectron]
   Float_t         Electron_r9[9];   //[nElectron]
   Float_t         Electron_scEtOverPt[9];   //[nElectron]
   Float_t         Electron_sieie[9];   //[nElectron]
   Float_t         Electron_sip3d[9];   //[nElectron]
   Float_t         Electron_mvaTTH[9];   //[nElectron]
   Int_t           Electron_charge[9];   //[nElectron]
   Int_t           Electron_cutBased[9];   //[nElectron]
   Int_t           Electron_jetIdx[9];   //[nElectron]
   Int_t           Electron_pdgId[9];   //[nElectron]
   Int_t           Electron_photonIdx[9];   //[nElectron]
   Int_t           Electron_tightCharge[9];   //[nElectron]
   Int_t           Electron_vidNestedWPBitmap[9];   //[nElectron]
   Int_t           Electron_vidNestedWPBitmapHEEP[9];   //[nElectron]
   Bool_t          Electron_convVeto[9];   //[nElectron]
   Bool_t          Electron_cutBased_HEEP[9];   //[nElectron]
   Bool_t          Electron_isPFcand[9];   //[nElectron]
   UChar_t         Electron_jetNDauCharged[9];   //[nElectron]
   UChar_t         Electron_lostHits[9];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WP80[9];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WP90[9];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WPL[9];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WP80[9];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WP90[9];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WPL[9];   //[nElectron]
   UChar_t         Electron_seedGain[9];   //[nElectron]
   UInt_t          nFatJet;
   Float_t         FatJet_area[7];   //[nFatJet]
   Float_t         FatJet_btagCSVV2[7];   //[nFatJet]
   Float_t         FatJet_btagDDBvLV2[7];   //[nFatJet]
   Float_t         FatJet_btagDDCvBV2[7];   //[nFatJet]
   Float_t         FatJet_btagDDCvLV2[7];   //[nFatJet]
   Float_t         FatJet_btagDeepB[7];   //[nFatJet]
   Float_t         FatJet_btagHbb[7];   //[nFatJet]
   Float_t         FatJet_deepTagMD_H4qvsQCD[7];   //[nFatJet]
   Float_t         FatJet_deepTagMD_HbbvsQCD[7];   //[nFatJet]
   Float_t         FatJet_deepTagMD_TvsQCD[7];   //[nFatJet]
   Float_t         FatJet_deepTagMD_WvsQCD[7];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZHbbvsQCD[7];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZHccvsQCD[7];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZbbvsQCD[7];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZvsQCD[7];   //[nFatJet]
   Float_t         FatJet_deepTagMD_bbvsLight[7];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ccvsLight[7];   //[nFatJet]
   Float_t         FatJet_deepTag_H[7];   //[nFatJet]
   Float_t         FatJet_deepTag_QCD[7];   //[nFatJet]
   Float_t         FatJet_deepTag_QCDothers[7];   //[nFatJet]
   Float_t         FatJet_deepTag_TvsQCD[7];   //[nFatJet]
   Float_t         FatJet_deepTag_WvsQCD[7];   //[nFatJet]
   Float_t         FatJet_deepTag_ZvsQCD[7];   //[nFatJet]
   Float_t         FatJet_eta[7];   //[nFatJet]
   Float_t         FatJet_mass[7];   //[nFatJet]
   Float_t         FatJet_msoftdrop[7];   //[nFatJet]
   Float_t         FatJet_n2b1[7];   //[nFatJet]
   Float_t         FatJet_n3b1[7];   //[nFatJet]
   Float_t         FatJet_particleNetMD_QCD[7];   //[nFatJet]
   Float_t         FatJet_particleNetMD_Xbb[7];   //[nFatJet]
   Float_t         FatJet_particleNetMD_Xcc[7];   //[nFatJet]
   Float_t         FatJet_particleNetMD_Xqq[7];   //[nFatJet]
   Float_t         FatJet_particleNet_H4qvsQCD[7];   //[nFatJet]
   Float_t         FatJet_particleNet_HbbvsQCD[7];   //[nFatJet]
   Float_t         FatJet_particleNet_HccvsQCD[7];   //[nFatJet]
   Float_t         FatJet_particleNet_QCD[7];   //[nFatJet]
   Float_t         FatJet_particleNet_TvsQCD[7];   //[nFatJet]
   Float_t         FatJet_particleNet_WvsQCD[7];   //[nFatJet]
   Float_t         FatJet_particleNet_ZvsQCD[7];   //[nFatJet]
   Float_t         FatJet_particleNet_mass[7];   //[nFatJet]
   Float_t         FatJet_phi[7];   //[nFatJet]
   Float_t         FatJet_pt[7];   //[nFatJet]
   Float_t         FatJet_rawFactor[7];   //[nFatJet]
   Float_t         FatJet_tau1[7];   //[nFatJet]
   Float_t         FatJet_tau2[7];   //[nFatJet]
   Float_t         FatJet_tau3[7];   //[nFatJet]
   Float_t         FatJet_tau4[7];   //[nFatJet]
   Float_t         FatJet_lsf3[7];   //[nFatJet]
   Int_t           FatJet_jetId[7];   //[nFatJet]
   Int_t           FatJet_subJetIdx1[7];   //[nFatJet]
   Int_t           FatJet_subJetIdx2[7];   //[nFatJet]
   Int_t           FatJet_electronIdx3SJ[7];   //[nFatJet]
   Int_t           FatJet_muonIdx3SJ[7];   //[nFatJet]
   UChar_t         FatJet_nConstituents[7];   //[nFatJet]
   UInt_t          nFsrPhoton;
   Float_t         FsrPhoton_dROverEt2[4];   //[nFsrPhoton]
   Float_t         FsrPhoton_eta[4];   //[nFsrPhoton]
   Float_t         FsrPhoton_phi[4];   //[nFsrPhoton]
   Float_t         FsrPhoton_pt[4];   //[nFsrPhoton]
   Float_t         FsrPhoton_relIso03[4];   //[nFsrPhoton]
   Int_t           FsrPhoton_muonIdx[4];   //[nFsrPhoton]
   UInt_t          nGenJetAK8;
   Float_t         GenJetAK8_eta[8];   //[nGenJetAK8]
   Float_t         GenJetAK8_mass[8];   //[nGenJetAK8]
   Float_t         GenJetAK8_phi[8];   //[nGenJetAK8]
   Float_t         GenJetAK8_pt[8];   //[nGenJetAK8]
   UInt_t          nGenJet;
   Float_t         GenJet_eta[25];   //[nGenJet]
   Float_t         GenJet_mass[25];   //[nGenJet]
   Float_t         GenJet_phi[25];   //[nGenJet]
   Float_t         GenJet_pt[25];   //[nGenJet]
   UInt_t          nGenPart;
   Float_t         GenPart_eta[160];   //[nGenPart]
   Float_t         GenPart_mass[160];   //[nGenPart]
   Float_t         GenPart_phi[160];   //[nGenPart]
   Float_t         GenPart_pt[160];   //[nGenPart]
   Int_t           GenPart_genPartIdxMother[160];   //[nGenPart]
   Int_t           GenPart_pdgId[160];   //[nGenPart]
   Int_t           GenPart_status[160];   //[nGenPart]
   Int_t           GenPart_statusFlags[160];   //[nGenPart]
   UInt_t          nSubGenJetAK8;
   Float_t         SubGenJetAK8_eta[16];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_mass[16];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_phi[16];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_pt[16];   //[nSubGenJetAK8]
   Float_t         Generator_binvar;
   Float_t         Generator_scalePDF;
   Float_t         Generator_weight;
   Float_t         Generator_x1;
   Float_t         Generator_x2;
   Float_t         Generator_xpdf1;
   Float_t         Generator_xpdf2;
   Int_t           Generator_id1;
   Int_t           Generator_id2;
   Float_t         GenVtx_x;
   Float_t         GenVtx_y;
   Float_t         GenVtx_z;
   UInt_t          nGenVisTau;
   Float_t         GenVisTau_eta[4];   //[nGenVisTau]
   Float_t         GenVisTau_mass[4];   //[nGenVisTau]
   Float_t         GenVisTau_phi[4];   //[nGenVisTau]
   Float_t         GenVisTau_pt[4];   //[nGenVisTau]
   Int_t           GenVisTau_charge[4];   //[nGenVisTau]
   Int_t           GenVisTau_genPartIdxMother[4];   //[nGenVisTau]
   Int_t           GenVisTau_status[4];   //[nGenVisTau]
   Float_t         genWeight;
   UInt_t          nEFTfitCoefficients;
   Float_t         EFTfitCoefficients[6];   //[nEFTfitCoefficients]
   Float_t         LHEWeight_originalXWGTUP;
   UInt_t          nLHEPdfWeight;
   Float_t         LHEPdfWeight[103];   //[nLHEPdfWeight]
   UInt_t          nLHEReweightingWeight;
   Float_t         LHEReweightingWeight[1];   //[nLHEReweightingWeight]
   UInt_t          nLHEScaleWeight;
   Float_t         LHEScaleWeight[9];   //[nLHEScaleWeight]
   UInt_t          nPSWeight;
   Float_t         PSWeight[4];   //[nPSWeight]
   UInt_t          nWCnames;
   Int_t           WCnames[2];   //[nWCnames]
   UInt_t          nIsoTrack;
   Float_t         IsoTrack_dxy[20];   //[nIsoTrack]
   Float_t         IsoTrack_dz[20];   //[nIsoTrack]
   Float_t         IsoTrack_eta[20];   //[nIsoTrack]
   Float_t         IsoTrack_pfRelIso03_all[20];   //[nIsoTrack]
   Float_t         IsoTrack_pfRelIso03_chg[20];   //[nIsoTrack]
   Float_t         IsoTrack_phi[20];   //[nIsoTrack]
   Float_t         IsoTrack_pt[20];   //[nIsoTrack]
   Float_t         IsoTrack_miniPFRelIso_all[20];   //[nIsoTrack]
   Float_t         IsoTrack_miniPFRelIso_chg[20];   //[nIsoTrack]
   Int_t           IsoTrack_charge[20];   //[nIsoTrack]
   Int_t           IsoTrack_fromPV[20];   //[nIsoTrack]
   Int_t           IsoTrack_pdgId[20];   //[nIsoTrack]
   Bool_t          IsoTrack_isHighPurityTrack[20];   //[nIsoTrack]
   Bool_t          IsoTrack_isPFcand[20];   //[nIsoTrack]
   Bool_t          IsoTrack_isFromLostTrack[20];   //[nIsoTrack]
   UInt_t          nJet;
   Float_t         Jet_area[30];   //[nJet]
   Float_t         Jet_btagCSVV2[30];   //[nJet]
   Float_t         Jet_btagDeepB[30];   //[nJet]
   Float_t         Jet_btagDeepCvB[30];   //[nJet]
   Float_t         Jet_btagDeepCvL[30];   //[nJet]
   Float_t         Jet_btagDeepFlavB[30];   //[nJet]
   Float_t         Jet_btagDeepFlavCvB[30];   //[nJet]
   Float_t         Jet_btagDeepFlavCvL[30];   //[nJet]
   Float_t         Jet_btagDeepFlavQG[30];   //[nJet]
   Float_t         Jet_chEmEF[30];   //[nJet]
   Float_t         Jet_chFPV0EF[30];   //[nJet]
   Float_t         Jet_chHEF[30];   //[nJet]
   Float_t         Jet_eta[30];   //[nJet]
   Float_t         Jet_hfsigmaEtaEta[30];   //[nJet]
   Float_t         Jet_hfsigmaPhiPhi[30];   //[nJet]
   Float_t         Jet_mass[30];   //[nJet]
   Float_t         Jet_muEF[30];   //[nJet]
   Float_t         Jet_muonSubtrFactor[30];   //[nJet]
   Float_t         Jet_neEmEF[30];   //[nJet]
   Float_t         Jet_neHEF[30];   //[nJet]
   Float_t         Jet_phi[30];   //[nJet]
   Float_t         Jet_pt[30];   //[nJet]
   Float_t         Jet_puIdDisc[30];   //[nJet]
   Float_t         Jet_qgl[30];   //[nJet]
   Float_t         Jet_rawFactor[30];   //[nJet]
   Float_t         Jet_bRegCorr[30];   //[nJet]
   Float_t         Jet_bRegRes[30];   //[nJet]
   Float_t         Jet_cRegCorr[30];   //[nJet]
   Float_t         Jet_cRegRes[30];   //[nJet]
   Int_t           Jet_electronIdx1[30];   //[nJet]
   Int_t           Jet_electronIdx2[30];   //[nJet]
   Int_t           Jet_hfadjacentEtaStripsSize[30];   //[nJet]
   Int_t           Jet_hfcentralEtaStripSize[30];   //[nJet]
   Int_t           Jet_jetId[30];   //[nJet]
   Int_t           Jet_muonIdx1[30];   //[nJet]
   Int_t           Jet_muonIdx2[30];   //[nJet]
   Int_t           Jet_nElectrons[30];   //[nJet]
   Int_t           Jet_nMuons[30];   //[nJet]
   Int_t           Jet_puId[30];   //[nJet]
   UChar_t         Jet_nConstituents[30];   //[nJet]
   Float_t         L1PreFiringWeight_Dn;
   Float_t         L1PreFiringWeight_ECAL_Dn;
   Float_t         L1PreFiringWeight_ECAL_Nom;
   Float_t         L1PreFiringWeight_ECAL_Up;
   Float_t         L1PreFiringWeight_Muon_Nom;
   Float_t         L1PreFiringWeight_Muon_StatDn;
   Float_t         L1PreFiringWeight_Muon_StatUp;
   Float_t         L1PreFiringWeight_Muon_SystDn;
   Float_t         L1PreFiringWeight_Muon_SystUp;
   Float_t         L1PreFiringWeight_Nom;
   Float_t         L1PreFiringWeight_Up;
   Float_t         LHE_HT;
   Float_t         LHE_HTIncoming;
   Float_t         LHE_Vpt;
   Float_t         LHE_AlphaS;
   UChar_t         LHE_Njets;
   UChar_t         LHE_Nb;
   UChar_t         LHE_Nc;
   UChar_t         LHE_Nuds;
   UChar_t         LHE_Nglu;
   UChar_t         LHE_NpNLO;
   UChar_t         LHE_NpLO;
   UInt_t          nLHEPart;
   Float_t         LHEPart_pt[11];   //[nLHEPart]
   Float_t         LHEPart_eta[11];   //[nLHEPart]
   Float_t         LHEPart_phi[11];   //[nLHEPart]
   Float_t         LHEPart_mass[11];   //[nLHEPart]
   Float_t         LHEPart_incomingpz[11];   //[nLHEPart]
   Int_t           LHEPart_pdgId[11];   //[nLHEPart]
   Int_t           LHEPart_status[11];   //[nLHEPart]
   Int_t           LHEPart_spin[11];   //[nLHEPart]
   UInt_t          nLowPtElectron;
   Float_t         LowPtElectron_ID[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_convVtxRadius[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_deltaEtaSC[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_dxy[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_dxyErr[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_dz[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_dzErr[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_eInvMinusPInv[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_embeddedID[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_energyErr[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_eta[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_hoe[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_mass[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_miniPFRelIso_all[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_miniPFRelIso_chg[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_phi[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_pt[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_ptbiased[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_r9[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_scEtOverPt[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_sieie[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_unbiased[9];   //[nLowPtElectron]
   Int_t           LowPtElectron_charge[9];   //[nLowPtElectron]
   Int_t           LowPtElectron_convWP[9];   //[nLowPtElectron]
   Int_t           LowPtElectron_pdgId[9];   //[nLowPtElectron]
   Bool_t          LowPtElectron_convVeto[9];   //[nLowPtElectron]
   UChar_t         LowPtElectron_lostHits[9];   //[nLowPtElectron]
   Float_t         GenMET_phi;
   Float_t         GenMET_pt;
   Float_t         MET_MetUnclustEnUpDeltaX;
   Float_t         MET_MetUnclustEnUpDeltaY;
   Float_t         MET_covXX;
   Float_t         MET_covXY;
   Float_t         MET_covYY;
   Float_t         MET_phi;
   Float_t         MET_pt;
   Float_t         MET_significance;
   Float_t         MET_sumEt;
   Float_t         MET_sumPtUnclustered;
   UInt_t          nMuon;
   Float_t         Muon_dxy[14];   //[nMuon]
   Float_t         Muon_dxyErr[14];   //[nMuon]
   Float_t         Muon_dxybs[14];   //[nMuon]
   Float_t         Muon_dz[14];   //[nMuon]
   Float_t         Muon_dzErr[14];   //[nMuon]
   Float_t         Muon_eta[14];   //[nMuon]
   Float_t         Muon_ip3d[14];   //[nMuon]
   Float_t         Muon_jetPtRelv2[14];   //[nMuon]
   Float_t         Muon_jetRelIso[14];   //[nMuon]
   Float_t         Muon_mass[14];   //[nMuon]
   Float_t         Muon_miniPFRelIso_all[14];   //[nMuon]
   Float_t         Muon_miniPFRelIso_chg[14];   //[nMuon]
   Float_t         Muon_pfRelIso03_all[14];   //[nMuon]
   Float_t         Muon_pfRelIso03_chg[14];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[14];   //[nMuon]
   Float_t         Muon_phi[14];   //[nMuon]
   Float_t         Muon_pt[14];   //[nMuon]
   Float_t         Muon_ptErr[14];   //[nMuon]
   Float_t         Muon_segmentComp[14];   //[nMuon]
   Float_t         Muon_sip3d[14];   //[nMuon]
   Float_t         Muon_softMva[14];   //[nMuon]
   Float_t         Muon_tkRelIso[14];   //[nMuon]
   Float_t         Muon_tunepRelPt[14];   //[nMuon]
   Float_t         Muon_mvaLowPt[14];   //[nMuon]
   Float_t         Muon_mvaTTH[14];   //[nMuon]
   Int_t           Muon_charge[14];   //[nMuon]
   Int_t           Muon_jetIdx[14];   //[nMuon]
   Int_t           Muon_nStations[14];   //[nMuon]
   Int_t           Muon_nTrackerLayers[14];   //[nMuon]
   Int_t           Muon_pdgId[14];   //[nMuon]
   Int_t           Muon_tightCharge[14];   //[nMuon]
   Int_t           Muon_fsrPhotonIdx[14];   //[nMuon]
   UChar_t         Muon_highPtId[14];   //[nMuon]
   Bool_t          Muon_highPurity[14];   //[nMuon]
   Bool_t          Muon_inTimeMuon[14];   //[nMuon]
   Bool_t          Muon_isGlobal[14];   //[nMuon]
   Bool_t          Muon_isPFcand[14];   //[nMuon]
   Bool_t          Muon_isStandalone[14];   //[nMuon]
   Bool_t          Muon_isTracker[14];   //[nMuon]
   UChar_t         Muon_jetNDauCharged[14];   //[nMuon]
   Bool_t          Muon_looseId[14];   //[nMuon]
   Bool_t          Muon_mediumId[14];   //[nMuon]
   Bool_t          Muon_mediumPromptId[14];   //[nMuon]
   UChar_t         Muon_miniIsoId[14];   //[nMuon]
   UChar_t         Muon_multiIsoId[14];   //[nMuon]
   UChar_t         Muon_mvaId[14];   //[nMuon]
   UChar_t         Muon_mvaLowPtId[14];   //[nMuon]
   UChar_t         Muon_pfIsoId[14];   //[nMuon]
   UChar_t         Muon_puppiIsoId[14];   //[nMuon]
   Bool_t          Muon_softId[14];   //[nMuon]
   Bool_t          Muon_softMvaId[14];   //[nMuon]
   Bool_t          Muon_tightId[14];   //[nMuon]
   UChar_t         Muon_tkIsoId[14];   //[nMuon]
   Bool_t          Muon_triggerIdLoose[14];   //[nMuon]
   UInt_t          nPhoton;
   Float_t         Photon_dEscaleDown[10];   //[nPhoton]
   Float_t         Photon_dEscaleUp[10];   //[nPhoton]
   Float_t         Photon_dEsigmaDown[10];   //[nPhoton]
   Float_t         Photon_dEsigmaUp[10];   //[nPhoton]
   Float_t         Photon_eCorr[10];   //[nPhoton]
   Float_t         Photon_energyErr[10];   //[nPhoton]
   Float_t         Photon_eta[10];   //[nPhoton]
   Float_t         Photon_hoe[10];   //[nPhoton]
   Float_t         Photon_mass[10];   //[nPhoton]
   Float_t         Photon_mvaID[10];   //[nPhoton]
   Float_t         Photon_mvaID_Fall17V1p1[10];   //[nPhoton]
   Float_t         Photon_pfRelIso03_all[10];   //[nPhoton]
   Float_t         Photon_pfRelIso03_chg[10];   //[nPhoton]
   Float_t         Photon_phi[10];   //[nPhoton]
   Float_t         Photon_pt[10];   //[nPhoton]
   Float_t         Photon_r9[10];   //[nPhoton]
   Float_t         Photon_sieie[10];   //[nPhoton]
   Int_t           Photon_charge[10];   //[nPhoton]
   Int_t           Photon_cutBased[10];   //[nPhoton]
   Int_t           Photon_cutBased_Fall17V1Bitmap[10];   //[nPhoton]
   Int_t           Photon_electronIdx[10];   //[nPhoton]
   Int_t           Photon_jetIdx[10];   //[nPhoton]
   Int_t           Photon_pdgId[10];   //[nPhoton]
   Int_t           Photon_vidNestedWPBitmap[10];   //[nPhoton]
   Bool_t          Photon_electronVeto[10];   //[nPhoton]
   Bool_t          Photon_isScEtaEB[10];   //[nPhoton]
   Bool_t          Photon_isScEtaEE[10];   //[nPhoton]
   Bool_t          Photon_mvaID_WP80[10];   //[nPhoton]
   Bool_t          Photon_mvaID_WP90[10];   //[nPhoton]
   Bool_t          Photon_pixelSeed[10];   //[nPhoton]
   UChar_t         Photon_seedGain[10];   //[nPhoton]
   Float_t         Pileup_nTrueInt;
   Float_t         Pileup_pudensity;
   Float_t         Pileup_gpudensity;
   Int_t           Pileup_nPU;
   Int_t           Pileup_sumEOOT;
   Int_t           Pileup_sumLOOT;
   Float_t         PuppiMET_phi;
   Float_t         PuppiMET_phiJERDown;
   Float_t         PuppiMET_phiJERUp;
   Float_t         PuppiMET_phiJESDown;
   Float_t         PuppiMET_phiJESUp;
   Float_t         PuppiMET_phiUnclusteredDown;
   Float_t         PuppiMET_phiUnclusteredUp;
   Float_t         PuppiMET_pt;
   Float_t         PuppiMET_ptJERDown;
   Float_t         PuppiMET_ptJERUp;
   Float_t         PuppiMET_ptJESDown;
   Float_t         PuppiMET_ptJESUp;
   Float_t         PuppiMET_ptUnclusteredDown;
   Float_t         PuppiMET_ptUnclusteredUp;
   Float_t         PuppiMET_sumEt;
   Float_t         RawMET_phi;
   Float_t         RawMET_pt;
   Float_t         RawMET_sumEt;
   Float_t         RawPuppiMET_phi;
   Float_t         RawPuppiMET_pt;
   Float_t         RawPuppiMET_sumEt;
   Float_t         fixedGridRhoFastjetAll;
   Float_t         fixedGridRhoFastjetCentral;
   Float_t         fixedGridRhoFastjetCentralCalo;
   Float_t         fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         fixedGridRhoFastjetCentralNeutral;
   UInt_t          nGenDressedLepton;
   Float_t         GenDressedLepton_eta[4];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_mass[4];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_phi[4];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_pt[4];   //[nGenDressedLepton]
   Int_t           GenDressedLepton_pdgId[4];   //[nGenDressedLepton]
   Bool_t          GenDressedLepton_hasTauAnc[4];   //[nGenDressedLepton]
   UInt_t          nGenIsolatedPhoton;
   Float_t         GenIsolatedPhoton_eta[2];   //[nGenIsolatedPhoton]
   Float_t         GenIsolatedPhoton_mass[2];   //[nGenIsolatedPhoton]
   Float_t         GenIsolatedPhoton_phi[2];   //[nGenIsolatedPhoton]
   Float_t         GenIsolatedPhoton_pt[2];   //[nGenIsolatedPhoton]
   UInt_t          nSoftActivityJet;
   Float_t         SoftActivityJet_eta[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_phi[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_pt[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJetHT;
   Float_t         SoftActivityJetHT10;
   Float_t         SoftActivityJetHT2;
   Float_t         SoftActivityJetHT5;
   Int_t           SoftActivityJetNjets10;
   Int_t           SoftActivityJetNjets2;
   Int_t           SoftActivityJetNjets5;
   UInt_t          nSubJet;
   Float_t         SubJet_btagCSVV2[12];   //[nSubJet]
   Float_t         SubJet_btagDeepB[12];   //[nSubJet]
   Float_t         SubJet_eta[12];   //[nSubJet]
   Float_t         SubJet_mass[12];   //[nSubJet]
   Float_t         SubJet_n2b1[12];   //[nSubJet]
   Float_t         SubJet_n3b1[12];   //[nSubJet]
   Float_t         SubJet_phi[12];   //[nSubJet]
   Float_t         SubJet_pt[12];   //[nSubJet]
   Float_t         SubJet_rawFactor[12];   //[nSubJet]
   Float_t         SubJet_tau1[12];   //[nSubJet]
   Float_t         SubJet_tau2[12];   //[nSubJet]
   Float_t         SubJet_tau3[12];   //[nSubJet]
   Float_t         SubJet_tau4[12];   //[nSubJet]
   UInt_t          nTau;
   Float_t         Tau_chargedIso[7];   //[nTau]
   Float_t         Tau_dxy[7];   //[nTau]
   Float_t         Tau_dz[7];   //[nTau]
   Float_t         Tau_eta[7];   //[nTau]
   Float_t         Tau_leadTkDeltaEta[7];   //[nTau]
   Float_t         Tau_leadTkDeltaPhi[7];   //[nTau]
   Float_t         Tau_leadTkPtOverTauPt[7];   //[nTau]
   Float_t         Tau_mass[7];   //[nTau]
   Float_t         Tau_neutralIso[7];   //[nTau]
   Float_t         Tau_phi[7];   //[nTau]
   Float_t         Tau_photonsOutsideSignalCone[7];   //[nTau]
   Float_t         Tau_pt[7];   //[nTau]
   Float_t         Tau_puCorr[7];   //[nTau]
   Float_t         Tau_rawDeepTau2017v2p1VSe[7];   //[nTau]
   Float_t         Tau_rawDeepTau2017v2p1VSjet[7];   //[nTau]
   Float_t         Tau_rawDeepTau2017v2p1VSmu[7];   //[nTau]
   Float_t         Tau_rawIso[7];   //[nTau]
   Float_t         Tau_rawIsodR03[7];   //[nTau]
   Int_t           Tau_charge[7];   //[nTau]
   Int_t           Tau_decayMode[7];   //[nTau]
   Int_t           Tau_jetIdx[7];   //[nTau]
   Bool_t          Tau_idAntiEleDeadECal[7];   //[nTau]
   UChar_t         Tau_idAntiMu[7];   //[nTau]
   Bool_t          Tau_idDecayModeOldDMs[7];   //[nTau]
   UChar_t         Tau_idDeepTau2017v2p1VSe[7];   //[nTau]
   UChar_t         Tau_idDeepTau2017v2p1VSjet[7];   //[nTau]
   UChar_t         Tau_idDeepTau2017v2p1VSmu[7];   //[nTau]
   Float_t         TkMET_phi;
   Float_t         TkMET_pt;
   Float_t         TkMET_sumEt;
   UInt_t          nTrigObj;
   Float_t         TrigObj_pt[51];   //[nTrigObj]
   Float_t         TrigObj_eta[51];   //[nTrigObj]
   Float_t         TrigObj_phi[51];   //[nTrigObj]
   Float_t         TrigObj_l1pt[51];   //[nTrigObj]
   Float_t         TrigObj_l1pt_2[51];   //[nTrigObj]
   Float_t         TrigObj_l2pt[51];   //[nTrigObj]
   Int_t           TrigObj_id[51];   //[nTrigObj]
   Int_t           TrigObj_l1iso[51];   //[nTrigObj]
   Int_t           TrigObj_l1charge[51];   //[nTrigObj]
   Int_t           TrigObj_filterBits[51];   //[nTrigObj]
   Int_t           genTtbarId;
   UInt_t          nOtherPV;
   Float_t         OtherPV_z[3];   //[nOtherPV]
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   UInt_t          nSV;
   Float_t         SV_dlen[19];   //[nSV]
   Float_t         SV_dlenSig[19];   //[nSV]
   Float_t         SV_dxy[19];   //[nSV]
   Float_t         SV_dxySig[19];   //[nSV]
   Float_t         SV_pAngle[19];   //[nSV]
   Int_t           SV_charge[19];   //[nSV]
   Int_t           boostedTau_genPartIdx[5];   //[nboostedTau]
   UChar_t         boostedTau_genPartFlav[5];   //[nboostedTau]
   Int_t           Electron_genPartIdx[9];   //[nElectron]
   UChar_t         Electron_genPartFlav[9];   //[nElectron]
   Int_t           FatJet_genJetAK8Idx[7];   //[nFatJet]
   Int_t           FatJet_hadronFlavour[7];   //[nFatJet]
   UChar_t         FatJet_nBHadrons[7];   //[nFatJet]
   UChar_t         FatJet_nCHadrons[7];   //[nFatJet]
   Int_t           GenJetAK8_partonFlavour[8];   //[nGenJetAK8]
   UChar_t         GenJetAK8_hadronFlavour[8];   //[nGenJetAK8]
   Int_t           GenJet_partonFlavour[25];   //[nGenJet]
   UChar_t         GenJet_hadronFlavour[25];   //[nGenJet]
   Float_t         GenVtx_t0;
   Int_t           Jet_genJetIdx[30];   //[nJet]
   Int_t           Jet_hadronFlavour[30];   //[nJet]
   Int_t           Jet_partonFlavour[30];   //[nJet]
   Int_t           LowPtElectron_genPartIdx[9];   //[nLowPtElectron]
   UChar_t         LowPtElectron_genPartFlav[9];   //[nLowPtElectron]
   Int_t           Muon_genPartIdx[14];   //[nMuon]
   UChar_t         Muon_genPartFlav[14];   //[nMuon]
   Int_t           Photon_genPartIdx[10];   //[nPhoton]
   UChar_t         Photon_genPartFlav[10];   //[nPhoton]
   Float_t         MET_fiducialGenPhi;
   Float_t         MET_fiducialGenPt;
   UChar_t         Electron_cleanmask[9];   //[nElectron]
   UChar_t         Jet_cleanmask[30];   //[nJet]
   UChar_t         Muon_cleanmask[14];   //[nMuon]
   UChar_t         Photon_cleanmask[10];   //[nPhoton]
   UChar_t         Tau_cleanmask[7];   //[nTau]
   Int_t           SubJet_hadronFlavour[12];   //[nSubJet]
   UChar_t         SubJet_nBHadrons[12];   //[nSubJet]
   UChar_t         SubJet_nCHadrons[12];   //[nSubJet]
   Float_t         SV_chi2[19];   //[nSV]
   Float_t         SV_eta[19];   //[nSV]
   Float_t         SV_mass[19];   //[nSV]
   Float_t         SV_ndof[19];   //[nSV]
   Float_t         SV_phi[19];   //[nSV]
   Float_t         SV_pt[19];   //[nSV]
   Float_t         SV_x[19];   //[nSV]
   Float_t         SV_y[19];   //[nSV]
   Float_t         SV_z[19];   //[nSV]
   UChar_t         SV_ntracks[19];   //[nSV]
   Int_t           Tau_genPartIdx[7];   //[nTau]
   UChar_t         Tau_genPartFlav[7];   //[nTau]
   Bool_t          L1_AlwaysTrue;
   Bool_t          L1_BPTX_AND_Ref1_VME;
   Bool_t          L1_BPTX_AND_Ref3_VME;
   Bool_t          L1_BPTX_AND_Ref4_VME;
   Bool_t          L1_BPTX_BeamGas_B1_VME;
   Bool_t          L1_BPTX_BeamGas_B2_VME;
   Bool_t          L1_BPTX_BeamGas_Ref1_VME;
   Bool_t          L1_BPTX_BeamGas_Ref2_VME;
   Bool_t          L1_BPTX_NotOR_VME;
   Bool_t          L1_BPTX_OR_Ref3_VME;
   Bool_t          L1_BPTX_OR_Ref4_VME;
   Bool_t          L1_BPTX_RefAND_VME;
   Bool_t          L1_BptxMinus;
   Bool_t          L1_BptxOR;
   Bool_t          L1_BptxPlus;
   Bool_t          L1_BptxXOR;
   Bool_t          L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;
   Bool_t          L1_DoubleEG6_HTT240er;
   Bool_t          L1_DoubleEG6_HTT250er;
   Bool_t          L1_DoubleEG6_HTT255er;
   Bool_t          L1_DoubleEG6_HTT270er;
   Bool_t          L1_DoubleEG6_HTT300er;
   Bool_t          L1_DoubleEG8er2p6_HTT255er;
   Bool_t          L1_DoubleEG8er2p6_HTT270er;
   Bool_t          L1_DoubleEG8er2p6_HTT300er;
   Bool_t          L1_DoubleEG_15_10;
   Bool_t          L1_DoubleEG_18_17;
   Bool_t          L1_DoubleEG_20_18;
   Bool_t          L1_DoubleEG_22_10;
   Bool_t          L1_DoubleEG_22_12;
   Bool_t          L1_DoubleEG_22_15;
   Bool_t          L1_DoubleEG_23_10;
   Bool_t          L1_DoubleEG_24_17;
   Bool_t          L1_DoubleEG_25_12;
   Bool_t          L1_DoubleEG_25_13;
   Bool_t          L1_DoubleEG_25_14;
   Bool_t          L1_DoubleEG_LooseIso23_10;
   Bool_t          L1_DoubleEG_LooseIso24_10;
   Bool_t          L1_DoubleIsoTau28er2p1;
   Bool_t          L1_DoubleIsoTau30er2p1;
   Bool_t          L1_DoubleIsoTau32er2p1;
   Bool_t          L1_DoubleIsoTau33er2p1;
   Bool_t          L1_DoubleIsoTau34er2p1;
   Bool_t          L1_DoubleIsoTau35er2p1;
   Bool_t          L1_DoubleIsoTau36er2p1;
   Bool_t          L1_DoubleIsoTau38er2p1;
   Bool_t          L1_DoubleJet100er2p3_dEta_Max1p6;
   Bool_t          L1_DoubleJet100er2p7;
   Bool_t          L1_DoubleJet112er2p3_dEta_Max1p6;
   Bool_t          L1_DoubleJet112er2p7;
   Bool_t          L1_DoubleJet120er2p7;
   Bool_t          L1_DoubleJet150er2p7;
   Bool_t          L1_DoubleJet30_Mass_Min300_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min320_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min340_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min360_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min380_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min400_Mu10;
   Bool_t          L1_DoubleJet30_Mass_Min400_Mu6;
   Bool_t          L1_DoubleJet30_Mass_Min400_dEta_Max1p5;
   Bool_t          L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450;
   Bool_t          L1_DoubleJet40er2p7;
   Bool_t          L1_DoubleJet50er2p7;
   Bool_t          L1_DoubleJet60er2p7;
   Bool_t          L1_DoubleJet60er2p7_ETM100;
   Bool_t          L1_DoubleJet60er2p7_ETM60;
   Bool_t          L1_DoubleJet60er2p7_ETM70;
   Bool_t          L1_DoubleJet60er2p7_ETM80;
   Bool_t          L1_DoubleJet60er2p7_ETM90;
   Bool_t          L1_DoubleJet80er2p7;
   Bool_t          L1_DoubleJet_100_30_DoubleJet30_Mass_Min620;
   Bool_t          L1_DoubleJet_100_35_DoubleJet35_Mass_Min620;
   Bool_t          L1_DoubleJet_110_35_DoubleJet35_Mass_Min620;
   Bool_t          L1_DoubleJet_110_40_DoubleJet40_Mass_Min620;
   Bool_t          L1_DoubleJet_115_35_DoubleJet35_Mass_Min620;
   Bool_t          L1_DoubleJet_115_40_DoubleJet40_Mass_Min620;
   Bool_t          L1_DoubleJet_90_30_DoubleJet30_Mass_Min620;
   Bool_t          L1_DoubleLooseIsoEG22er2p1;
   Bool_t          L1_DoubleLooseIsoEG24er2p1;
   Bool_t          L1_DoubleMu0;
   Bool_t          L1_DoubleMu0_ETM40;
   Bool_t          L1_DoubleMu0_ETM55;
   Bool_t          L1_DoubleMu0_ETM60;
   Bool_t          L1_DoubleMu0_ETM65;
   Bool_t          L1_DoubleMu0_ETM70;
   Bool_t          L1_DoubleMu0_SQ;
   Bool_t          L1_DoubleMu0_SQ_OS;
   Bool_t          L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;
   Bool_t          L1_DoubleMu0er1p4_dEta_Max1p8_OS;
   Bool_t          L1_DoubleMu0er1p5_SQ_OS;
   Bool_t          L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;
   Bool_t          L1_DoubleMu0er1p5_SQ_dR_Max1p4;
   Bool_t          L1_DoubleMu0er2_SQ_dR_Max1p4;
   Bool_t          L1_DoubleMu18er2p1;
   Bool_t          L1_DoubleMu22er2p1;
   Bool_t          L1_DoubleMu3_OS_DoubleEG7p5Upsilon;
   Bool_t          L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_HTT100er;
   Bool_t          L1_DoubleMu3_SQ_HTT200er;
   Bool_t          L1_DoubleMu3_SQ_HTT220er;
   Bool_t          L1_DoubleMu3_SQ_HTT240er;
   Bool_t          L1_DoubleMu4_OS_EG12;
   Bool_t          L1_DoubleMu4_SQ_OS;
   Bool_t          L1_DoubleMu4_SQ_OS_dR_Max1p2;
   Bool_t          L1_DoubleMu4p5_SQ;
   Bool_t          L1_DoubleMu4p5_SQ_OS;
   Bool_t          L1_DoubleMu4p5_SQ_OS_dR_Max1p2;
   Bool_t          L1_DoubleMu4p5er2p0_SQ_OS;
   Bool_t          L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18;
   Bool_t          L1_DoubleMu5Upsilon_OS_DoubleEG3;
   Bool_t          L1_DoubleMu5_OS_EG12;
   Bool_t          L1_DoubleMu5_SQ_OS;
   Bool_t          L1_DoubleMu5_SQ_OS_Mass7to18;
   Bool_t          L1_DoubleMu6_SQ_OS;
   Bool_t          L1_DoubleMu7_EG7;
   Bool_t          L1_DoubleMu7_SQ_EG7;
   Bool_t          L1_DoubleMu8_SQ;
   Bool_t          L1_DoubleMu_10_0_dEta_Max1p8;
   Bool_t          L1_DoubleMu_11_4;
   Bool_t          L1_DoubleMu_12_5;
   Bool_t          L1_DoubleMu_12_8;
   Bool_t          L1_DoubleMu_13_6;
   Bool_t          L1_DoubleMu_15_5;
   Bool_t          L1_DoubleMu_15_5_SQ;
   Bool_t          L1_DoubleMu_15_7;
   Bool_t          L1_DoubleMu_15_7_SQ;
   Bool_t          L1_DoubleMu_15_7_SQ_Mass_Min4;
   Bool_t          L1_DoubleMu_20_2_SQ_Mass_Max20;
   Bool_t          L1_DoubleTau50er2p1;
   Bool_t          L1_DoubleTau70er2p1;
   Bool_t          L1_EG25er2p1_HTT125er;
   Bool_t          L1_EG27er2p1_HTT200er;
   Bool_t          L1_ETM100;
   Bool_t          L1_ETM100_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM105;
   Bool_t          L1_ETM110;
   Bool_t          L1_ETM110_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM115;
   Bool_t          L1_ETM120;
   Bool_t          L1_ETM150;
   Bool_t          L1_ETM30;
   Bool_t          L1_ETM40;
   Bool_t          L1_ETM50;
   Bool_t          L1_ETM60;
   Bool_t          L1_ETM70;
   Bool_t          L1_ETM75;
   Bool_t          L1_ETM75_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM80;
   Bool_t          L1_ETM80_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM85;
   Bool_t          L1_ETM90;
   Bool_t          L1_ETM90_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM95;
   Bool_t          L1_ETMHF100;
   Bool_t          L1_ETMHF100_HTT60er;
   Bool_t          L1_ETMHF100_Jet60_OR_DiJet30woTT28;
   Bool_t          L1_ETMHF100_Jet60_OR_DoubleJet30;
   Bool_t          L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETMHF110;
   Bool_t          L1_ETMHF110_HTT60er;
   Bool_t          L1_ETMHF110_Jet60_OR_DiJet30woTT28;
   Bool_t          L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETMHF120;
   Bool_t          L1_ETMHF120_HTT60er;
   Bool_t          L1_ETMHF120_Jet60_OR_DiJet30woTT28;
   Bool_t          L1_ETMHF150;
   Bool_t          L1_ETMHF70;
   Bool_t          L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETMHF80;
   Bool_t          L1_ETMHF80_HTT60er;
   Bool_t          L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETMHF90;
   Bool_t          L1_ETMHF90_HTT60er;
   Bool_t          L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETT100_BptxAND;
   Bool_t          L1_ETT110_BptxAND;
   Bool_t          L1_ETT40_BptxAND;
   Bool_t          L1_ETT50_BptxAND;
   Bool_t          L1_ETT60_BptxAND;
   Bool_t          L1_ETT70_BptxAND;
   Bool_t          L1_ETT75_BptxAND;
   Bool_t          L1_ETT80_BptxAND;
   Bool_t          L1_ETT85_BptxAND;
   Bool_t          L1_ETT90_BptxAND;
   Bool_t          L1_ETT95_BptxAND;
   Bool_t          L1_FirstBunchAfterTrain;
   Bool_t          L1_FirstBunchInTrain;
   Bool_t          L1_FirstCollisionInOrbit;
   Bool_t          L1_FirstCollisionInTrain;
   Bool_t          L1_HTT120er;
   Bool_t          L1_HTT160er;
   Bool_t          L1_HTT200er;
   Bool_t          L1_HTT220er;
   Bool_t          L1_HTT240er;
   Bool_t          L1_HTT250er_QuadJet_70_55_40_35_er2p5;
   Bool_t          L1_HTT255er;
   Bool_t          L1_HTT270er;
   Bool_t          L1_HTT280er;
   Bool_t          L1_HTT280er_QuadJet_70_55_40_35_er2p5;
   Bool_t          L1_HTT300er;
   Bool_t          L1_HTT300er_QuadJet_70_55_40_35_er2p5;
   Bool_t          L1_HTT320er;
   Bool_t          L1_HTT320er_QuadJet_70_55_40_40_er2p4;
   Bool_t          L1_HTT320er_QuadJet_70_55_40_40_er2p5;
   Bool_t          L1_HTT320er_QuadJet_70_55_45_45_er2p5;
   Bool_t          L1_HTT340er;
   Bool_t          L1_HTT340er_QuadJet_70_55_40_40_er2p5;
   Bool_t          L1_HTT340er_QuadJet_70_55_45_45_er2p5;
   Bool_t          L1_HTT380er;
   Bool_t          L1_HTT400er;
   Bool_t          L1_HTT450er;
   Bool_t          L1_HTT500er;
   Bool_t          L1_IsoEG33_Mt40;
   Bool_t          L1_IsoEG33_Mt44;
   Bool_t          L1_IsoEG33_Mt48;
   Bool_t          L1_IsoTau40er_ETM100;
   Bool_t          L1_IsoTau40er_ETM105;
   Bool_t          L1_IsoTau40er_ETM110;
   Bool_t          L1_IsoTau40er_ETM115;
   Bool_t          L1_IsoTau40er_ETM120;
   Bool_t          L1_IsoTau40er_ETM80;
   Bool_t          L1_IsoTau40er_ETM85;
   Bool_t          L1_IsoTau40er_ETM90;
   Bool_t          L1_IsoTau40er_ETM95;
   Bool_t          L1_IsoTau40er_ETMHF100;
   Bool_t          L1_IsoTau40er_ETMHF110;
   Bool_t          L1_IsoTau40er_ETMHF120;
   Bool_t          L1_IsoTau40er_ETMHF80;
   Bool_t          L1_IsoTau40er_ETMHF90;
   Bool_t          L1_IsolatedBunch;
   Bool_t          L1_LastCollisionInTrain;
   Bool_t          L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3;
   Bool_t          L1_LooseIsoEG24er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3;
   Bool_t          L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3;
   Bool_t          L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7;
   Bool_t          L1_LooseIsoEG26er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3;
   Bool_t          L1_LooseIsoEG28er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3;
   Bool_t          L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3;
   Bool_t          L1_MU20_EG15;
   Bool_t          L1_MinimumBiasHF0_AND_BptxAND;
   Bool_t          L1_MinimumBiasHF0_OR_BptxAND;
   Bool_t          L1_Mu10er2p1_ETM30;
   Bool_t          L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6;
   Bool_t          L1_Mu12_EG10;
   Bool_t          L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6;
   Bool_t          L1_Mu14er2p1_ETM30;
   Bool_t          L1_Mu15_HTT100er;
   Bool_t          L1_Mu18_HTT100er;
   Bool_t          L1_Mu18_Jet24er2p7;
   Bool_t          L1_Mu18er2p1_IsoTau26er2p1;
   Bool_t          L1_Mu18er2p1_Tau24er2p1;
   Bool_t          L1_Mu20_EG10;
   Bool_t          L1_Mu20_EG17;
   Bool_t          L1_Mu20_LooseIsoEG6;
   Bool_t          L1_Mu20er2p1_IsoTau26er2p1;
   Bool_t          L1_Mu20er2p1_IsoTau27er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau28er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau30er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau32er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau33er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau34er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau35er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau36er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau38er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau40er2p1;
   Bool_t          L1_Mu22er2p1_Tau50er2p1;
   Bool_t          L1_Mu22er2p1_Tau70er2p1;
   Bool_t          L1_Mu23_EG10;
   Bool_t          L1_Mu23_LooseIsoEG10;
   Bool_t          L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4;
   Bool_t          L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4;
   Bool_t          L1_Mu3_Jet30er2p5;
   Bool_t          L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4;
   Bool_t          L1_Mu5_EG15;
   Bool_t          L1_Mu5_EG20;
   Bool_t          L1_Mu5_EG23;
   Bool_t          L1_Mu5_LooseIsoEG18;
   Bool_t          L1_Mu5_LooseIsoEG20;
   Bool_t          L1_Mu6_DoubleEG10;
   Bool_t          L1_Mu6_DoubleEG17;
   Bool_t          L1_Mu6_HTT200er;
   Bool_t          L1_Mu6_HTT240er;
   Bool_t          L1_Mu6_HTT250er;
   Bool_t          L1_Mu7_EG23;
   Bool_t          L1_Mu7_LooseIsoEG20;
   Bool_t          L1_Mu7_LooseIsoEG23;
   Bool_t          L1_Mu8_HTT150er;
   Bool_t          L1_NotBptxOR;
   Bool_t          L1_QuadJet36er2p7_IsoTau52er2p1;
   Bool_t          L1_QuadJet36er2p7_Tau52;
   Bool_t          L1_QuadJet40er2p7;
   Bool_t          L1_QuadJet50er2p7;
   Bool_t          L1_QuadJet60er2p7;
   Bool_t          L1_QuadMu0;
   Bool_t          L1_SingleEG10;
   Bool_t          L1_SingleEG15;
   Bool_t          L1_SingleEG18;
   Bool_t          L1_SingleEG24;
   Bool_t          L1_SingleEG26;
   Bool_t          L1_SingleEG28;
   Bool_t          L1_SingleEG2_BptxAND;
   Bool_t          L1_SingleEG30;
   Bool_t          L1_SingleEG32;
   Bool_t          L1_SingleEG34;
   Bool_t          L1_SingleEG34er2p1;
   Bool_t          L1_SingleEG36;
   Bool_t          L1_SingleEG36er2p1;
   Bool_t          L1_SingleEG38;
   Bool_t          L1_SingleEG38er2p1;
   Bool_t          L1_SingleEG40;
   Bool_t          L1_SingleEG42;
   Bool_t          L1_SingleEG45;
   Bool_t          L1_SingleEG5;
   Bool_t          L1_SingleEG50;
   Bool_t          L1_SingleIsoEG18;
   Bool_t          L1_SingleIsoEG18er2p1;
   Bool_t          L1_SingleIsoEG20;
   Bool_t          L1_SingleIsoEG20er2p1;
   Bool_t          L1_SingleIsoEG22;
   Bool_t          L1_SingleIsoEG22er2p1;
   Bool_t          L1_SingleIsoEG24;
   Bool_t          L1_SingleIsoEG24er2p1;
   Bool_t          L1_SingleIsoEG26;
   Bool_t          L1_SingleIsoEG26er2p1;
   Bool_t          L1_SingleIsoEG28;
   Bool_t          L1_SingleIsoEG28er2p1;
   Bool_t          L1_SingleIsoEG30;
   Bool_t          L1_SingleIsoEG30er2p1;
   Bool_t          L1_SingleIsoEG32;
   Bool_t          L1_SingleIsoEG32er2p1;
   Bool_t          L1_SingleIsoEG33er2p1;
   Bool_t          L1_SingleIsoEG34;
   Bool_t          L1_SingleIsoEG34er2p1;
   Bool_t          L1_SingleIsoEG35;
   Bool_t          L1_SingleIsoEG35er2p1;
   Bool_t          L1_SingleIsoEG36;
   Bool_t          L1_SingleIsoEG36er2p1;
   Bool_t          L1_SingleIsoEG37;
   Bool_t          L1_SingleIsoEG38;
   Bool_t          L1_SingleIsoEG38er2p1;
   Bool_t          L1_SingleIsoEG40;
   Bool_t          L1_SingleIsoEG40er2p1;
   Bool_t          L1_SingleJet120;
   Bool_t          L1_SingleJet120_FWD;
   Bool_t          L1_SingleJet12_BptxAND;
   Bool_t          L1_SingleJet140;
   Bool_t          L1_SingleJet150;
   Bool_t          L1_SingleJet16;
   Bool_t          L1_SingleJet160;
   Bool_t          L1_SingleJet170;
   Bool_t          L1_SingleJet180;
   Bool_t          L1_SingleJet20;
   Bool_t          L1_SingleJet200;
   Bool_t          L1_SingleJet20er2p7_NotBptxOR;
   Bool_t          L1_SingleJet20er2p7_NotBptxOR_3BX;
   Bool_t          L1_SingleJet35;
   Bool_t          L1_SingleJet35_FWD;
   Bool_t          L1_SingleJet35_HFm;
   Bool_t          L1_SingleJet35_HFp;
   Bool_t          L1_SingleJet43er2p7_NotBptxOR_3BX;
   Bool_t          L1_SingleJet46er2p7_NotBptxOR_3BX;
   Bool_t          L1_SingleJet60;
   Bool_t          L1_SingleJet60_FWD;
   Bool_t          L1_SingleJet60_HFm;
   Bool_t          L1_SingleJet60_HFp;
   Bool_t          L1_SingleJet90;
   Bool_t          L1_SingleJet90_FWD;
   Bool_t          L1_SingleMu0_BMTF;
   Bool_t          L1_SingleMu0_EMTF;
   Bool_t          L1_SingleMu0_OMTF;
   Bool_t          L1_SingleMu10_LowQ;
   Bool_t          L1_SingleMu11_LowQ;
   Bool_t          L1_SingleMu12_LowQ_BMTF;
   Bool_t          L1_SingleMu12_LowQ_EMTF;
   Bool_t          L1_SingleMu12_LowQ_OMTF;
   Bool_t          L1_SingleMu14er2p1;
   Bool_t          L1_SingleMu16;
   Bool_t          L1_SingleMu16er2p1;
   Bool_t          L1_SingleMu18;
   Bool_t          L1_SingleMu18er2p1;
   Bool_t          L1_SingleMu20;
   Bool_t          L1_SingleMu20er2p1;
   Bool_t          L1_SingleMu22;
   Bool_t          L1_SingleMu22_BMTF;
   Bool_t          L1_SingleMu22_EMTF;
   Bool_t          L1_SingleMu22_OMTF;
   Bool_t          L1_SingleMu22er2p1;
   Bool_t          L1_SingleMu25;
   Bool_t          L1_SingleMu3;
   Bool_t          L1_SingleMu30;
   Bool_t          L1_SingleMu5;
   Bool_t          L1_SingleMu7;
   Bool_t          L1_SingleMuCosmics;
   Bool_t          L1_SingleMuCosmics_BMTF;
   Bool_t          L1_SingleMuCosmics_EMTF;
   Bool_t          L1_SingleMuCosmics_OMTF;
   Bool_t          L1_SingleMuOpen;
   Bool_t          L1_SingleMuOpen_NotBptxOR;
   Bool_t          L1_SingleMuOpen_NotBptxOR_3BX;
   Bool_t          L1_SingleTau100er2p1;
   Bool_t          L1_SingleTau120er2p1;
   Bool_t          L1_SingleTau130er2p1;
   Bool_t          L1_SingleTau140er2p1;
   Bool_t          L1_SingleTau20;
   Bool_t          L1_SingleTau80er2p1;
   Bool_t          L1_TripleEG_14_10_8;
   Bool_t          L1_TripleEG_18_17_8;
   Bool_t          L1_TripleEG_LooseIso20_10_5;
   Bool_t          L1_TripleJet_100_85_72_VBF;
   Bool_t          L1_TripleJet_105_85_76_VBF;
   Bool_t          L1_TripleJet_84_68_48_VBF;
   Bool_t          L1_TripleJet_88_72_56_VBF;
   Bool_t          L1_TripleJet_92_76_64_VBF;
   Bool_t          L1_TripleJet_98_83_71_VBF;
   Bool_t          L1_TripleMu0;
   Bool_t          L1_TripleMu0_OQ;
   Bool_t          L1_TripleMu3;
   Bool_t          L1_TripleMu3_SQ;
   Bool_t          L1_TripleMu_4_4_4;
   Bool_t          L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17;
   Bool_t          L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14;
   Bool_t          L1_TripleMu_5SQ_3SQ_0OQ;
   Bool_t          L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9;
   Bool_t          L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9;
   Bool_t          L1_TripleMu_5_0_0;
   Bool_t          L1_TripleMu_5_3_3;
   Bool_t          L1_TripleMu_5_3p5_2p5;
   Bool_t          L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17;
   Bool_t          L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17;
   Bool_t          L1_TripleMu_5_5_3;
   Bool_t          L1_UnpairedBunchBptxMinus;
   Bool_t          L1_UnpairedBunchBptxPlus;
   Bool_t          L1_ZeroBias;
   Bool_t          L1_ZeroBias_copy;
   Bool_t          L1_UnprefireableEvent;
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHENoiseIsoFilter;
   Bool_t          Flag_CSCTightHaloFilter;
   Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter;
   Bool_t          Flag_CSCTightHalo2015Filter;
   Bool_t          Flag_globalTightHalo2016Filter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_HcalStripHaloFilter;
   Bool_t          Flag_hcalLaserEventFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_ecalLaserCorrFilter;
   Bool_t          Flag_trkPOGFilters;
   Bool_t          Flag_chargedHadronTrackResolutionFilter;
   Bool_t          Flag_muonBadTrackFilter;
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_BadPFMuonDzFilter;
   Bool_t          Flag_hfNoisyHitsFilter;
   Bool_t          Flag_BadChargedCandidateSummer16Filter;
   Bool_t          Flag_BadPFMuonSummer16Filter;
   Bool_t          Flag_trkPOG_manystripclus53X;
   Bool_t          Flag_trkPOG_toomanystripclus53X;
   Bool_t          Flag_trkPOG_logErrorTooManyClusters;
   Bool_t          Flag_METFilters;
   Bool_t          L1Reco_step;
   Bool_t          HLTriggerFirstPath;
   Bool_t          HLT_AK8PFJet360_TrimMass30;
   Bool_t          HLT_AK8PFJet380_TrimMass30;
   Bool_t          HLT_AK8PFJet400_TrimMass30;
   Bool_t          HLT_AK8PFJet420_TrimMass30;
   Bool_t          HLT_AK8PFHT750_TrimMass50;
   Bool_t          HLT_AK8PFHT800_TrimMass50;
   Bool_t          HLT_AK8PFHT850_TrimMass50;
   Bool_t          HLT_AK8PFHT900_TrimMass50;
   Bool_t          HLT_CaloJet500_NoJetID;
   Bool_t          HLT_CaloJet550_NoJetID;
   Bool_t          HLT_Trimuon5_3p5_2_Upsilon_Muon;
   Bool_t          HLT_DoubleEle25_CaloIdL_MW;
   Bool_t          HLT_DoubleEle27_CaloIdL_MW;
   Bool_t          HLT_DoubleEle33_CaloIdL_MW;
   Bool_t          HLT_DoubleEle24_eta2p1_WPTight_Gsf;
   Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;
   Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350;
   Bool_t          HLT_Ele27_Ele37_CaloIdL_MW;
   Bool_t          HLT_Mu27_Ele37_CaloIdL_MW;
   Bool_t          HLT_Mu37_Ele27_CaloIdL_MW;
   Bool_t          HLT_Mu37_TkMu27;
   Bool_t          HLT_DoubleMu4_3_Bs;
   Bool_t          HLT_DoubleMu4_3_Jpsi_Displaced;
   Bool_t          HLT_DoubleMu4_JpsiTrk_Displaced;
   Bool_t          HLT_DoubleMu4_LowMassNonResonantTrk_Displaced;
   Bool_t          HLT_DoubleMu3_Trk_Tau3mu;
   Bool_t          HLT_DoubleMu4_PsiPrimeTrk_Displaced;
   Bool_t          HLT_DoubleMu4_Mass8_DZ_PFHT350;
   Bool_t          HLT_DoubleMu8_Mass8_PFHT350;
   Bool_t          HLT_Mu3_PFJet40;
   Bool_t          HLT_Mu7p5_L2Mu2_Jpsi;
   Bool_t          HLT_Mu7p5_L2Mu2_Upsilon;
   Bool_t          HLT_Mu7p5_Track2_Jpsi;
   Bool_t          HLT_Mu7p5_Track3p5_Jpsi;
   Bool_t          HLT_Mu7p5_Track7_Jpsi;
   Bool_t          HLT_Mu7p5_Track2_Upsilon;
   Bool_t          HLT_Mu7p5_Track3p5_Upsilon;
   Bool_t          HLT_Mu7p5_Track7_Upsilon;
   Bool_t          HLT_DoublePhoton33_CaloIdL;
   Bool_t          HLT_DoublePhoton70;
   Bool_t          HLT_DoublePhoton85;
   Bool_t          HLT_Ele20_WPTight_Gsf;
   Bool_t          HLT_Ele20_WPLoose_Gsf;
   Bool_t          HLT_Ele20_eta2p1_WPLoose_Gsf;
   Bool_t          HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;
   Bool_t          HLT_Ele27_WPTight_Gsf;
   Bool_t          HLT_Ele32_WPTight_Gsf;
   Bool_t          HLT_Ele35_WPTight_Gsf;
   Bool_t          HLT_Ele35_WPTight_Gsf_L1EGMT;
   Bool_t          HLT_Ele38_WPTight_Gsf;
   Bool_t          HLT_Ele40_WPTight_Gsf;
   Bool_t          HLT_Ele32_WPTight_Gsf_L1DoubleEG;
   Bool_t          HLT_HT450_Beamspot;
   Bool_t          HLT_HT300_Beamspot;
   Bool_t          HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1;
   Bool_t          HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1;
   Bool_t          HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1;
   Bool_t          HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1;
   Bool_t          HLT_IsoMu20;
   Bool_t          HLT_IsoMu24;
   Bool_t          HLT_IsoTkMu24;
   Bool_t          HLT_IsoMu24_eta2p1;
   Bool_t          HLT_IsoMu27;
   Bool_t          HLT_IsoMu30;
   Bool_t          HLT_UncorrectedJetE30_NoBPTX;
   Bool_t          HLT_UncorrectedJetE30_NoBPTX3BX;
   Bool_t          HLT_UncorrectedJetE60_NoBPTX3BX;
   Bool_t          HLT_UncorrectedJetE70_NoBPTX3BX;
   Bool_t          HLT_L1SingleMu18;
   Bool_t          HLT_L1SingleMu25;
   Bool_t          HLT_L2Mu10;
   Bool_t          HLT_L2Mu10_NoVertex_NoBPTX3BX;
   Bool_t          HLT_L2Mu10_NoVertex_NoBPTX;
   Bool_t          HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;
   Bool_t          HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;
   Bool_t          HLT_L2Mu50;
   Bool_t          HLT_DoubleL2Mu50;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
   Bool_t          HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Mu25_TkMu0_Onia;
   Bool_t          HLT_Mu30_TkMu0_Onia;
   Bool_t          HLT_Mu20_TkMu0_Phi;
   Bool_t          HLT_Mu25_TkMu0_Phi;
   Bool_t          HLT_Mu20;
   Bool_t          HLT_Mu27;
   Bool_t          HLT_Mu50;
   Bool_t          HLT_Mu55;
   Bool_t          HLT_OldMu100;
   Bool_t          HLT_TkMu100;
   Bool_t          HLT_DiPFJet15_NoCaloMatched;
   Bool_t          HLT_DiPFJet25_NoCaloMatched;
   Bool_t          HLT_DiPFJet15_FBEta3_NoCaloMatched;
   Bool_t          HLT_DiPFJet25_FBEta3_NoCaloMatched;
   Bool_t          HLT_DiPFJetAve40;
   Bool_t          HLT_DiPFJetAve60;
   Bool_t          HLT_DiPFJetAve80;
   Bool_t          HLT_DiPFJetAve140;
   Bool_t          HLT_DiPFJetAve200;
   Bool_t          HLT_DiPFJetAve260;
   Bool_t          HLT_DiPFJetAve320;
   Bool_t          HLT_DiPFJetAve400;
   Bool_t          HLT_DiPFJetAve500;
   Bool_t          HLT_DiPFJetAve15_HFJEC;
   Bool_t          HLT_DiPFJetAve25_HFJEC;
   Bool_t          HLT_DiPFJetAve35_HFJEC;
   Bool_t          HLT_DiPFJetAve60_HFJEC;
   Bool_t          HLT_DiPFJetAve80_HFJEC;
   Bool_t          HLT_DiPFJetAve100_HFJEC;
   Bool_t          HLT_DiPFJetAve160_HFJEC;
   Bool_t          HLT_DiPFJetAve220_HFJEC;
   Bool_t          HLT_DiPFJetAve300_HFJEC;
   Bool_t          HLT_AK8PFJet40;
   Bool_t          HLT_AK8PFJet60;
   Bool_t          HLT_AK8PFJet80;
   Bool_t          HLT_AK8PFJet140;
   Bool_t          HLT_AK8PFJet200;
   Bool_t          HLT_AK8PFJet260;
   Bool_t          HLT_AK8PFJet320;
   Bool_t          HLT_AK8PFJet400;
   Bool_t          HLT_AK8PFJet450;
   Bool_t          HLT_AK8PFJet500;
   Bool_t          HLT_AK8PFJet550;
   Bool_t          HLT_PFJet40;
   Bool_t          HLT_PFJet60;
   Bool_t          HLT_PFJet80;
   Bool_t          HLT_PFJet140;
   Bool_t          HLT_PFJet200;
   Bool_t          HLT_PFJet260;
   Bool_t          HLT_PFJet320;
   Bool_t          HLT_PFJet400;
   Bool_t          HLT_PFJet450;
   Bool_t          HLT_PFJet500;
   Bool_t          HLT_PFJet550;
   Bool_t          HLT_PFJetFwd40;
   Bool_t          HLT_PFJetFwd60;
   Bool_t          HLT_PFJetFwd80;
   Bool_t          HLT_PFJetFwd140;
   Bool_t          HLT_PFJetFwd200;
   Bool_t          HLT_PFJetFwd260;
   Bool_t          HLT_PFJetFwd320;
   Bool_t          HLT_PFJetFwd400;
   Bool_t          HLT_PFJetFwd450;
   Bool_t          HLT_PFJetFwd500;
   Bool_t          HLT_AK8PFJetFwd40;
   Bool_t          HLT_AK8PFJetFwd60;
   Bool_t          HLT_AK8PFJetFwd80;
   Bool_t          HLT_AK8PFJetFwd140;
   Bool_t          HLT_AK8PFJetFwd200;
   Bool_t          HLT_AK8PFJetFwd260;
   Bool_t          HLT_AK8PFJetFwd320;
   Bool_t          HLT_AK8PFJetFwd400;
   Bool_t          HLT_AK8PFJetFwd450;
   Bool_t          HLT_AK8PFJetFwd500;
   Bool_t          HLT_PFHT180;
   Bool_t          HLT_PFHT250;
   Bool_t          HLT_PFHT370;
   Bool_t          HLT_PFHT430;
   Bool_t          HLT_PFHT510;
   Bool_t          HLT_PFHT590;
   Bool_t          HLT_PFHT680;
   Bool_t          HLT_PFHT780;
   Bool_t          HLT_PFHT890;
   Bool_t          HLT_PFHT1050;
   Bool_t          HLT_PFHT500_PFMET100_PFMHT100_IDTight;
   Bool_t          HLT_PFHT500_PFMET110_PFMHT110_IDTight;
   Bool_t          HLT_PFHT700_PFMET85_PFMHT85_IDTight;
   Bool_t          HLT_PFHT700_PFMET95_PFMHT95_IDTight;
   Bool_t          HLT_PFHT800_PFMET75_PFMHT75_IDTight;
   Bool_t          HLT_PFHT800_PFMET85_PFMHT85_IDTight;
   Bool_t          HLT_PFMET110_PFMHT110_IDTight;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight;
   Bool_t          HLT_PFMET130_PFMHT130_IDTight;
   Bool_t          HLT_PFMET140_PFMHT140_IDTight;
   Bool_t          HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1;
   Bool_t          HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1;
   Bool_t          HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1;
   Bool_t          HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight_PFHT60;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
   Bool_t          HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60;
   Bool_t          HLT_PFMETTypeOne110_PFMHT110_IDTight;
   Bool_t          HLT_PFMETTypeOne120_PFMHT120_IDTight;
   Bool_t          HLT_PFMETTypeOne130_PFMHT130_IDTight;
   Bool_t          HLT_PFMETTypeOne140_PFMHT140_IDTight;
   Bool_t          HLT_PFMETNoMu110_PFMHTNoMu110_IDTight;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
   Bool_t          HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;
   Bool_t          HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;
   Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight;
   Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight;
   Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight;
   Bool_t          HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight;
   Bool_t          HLT_L1ETMHadSeeds;
   Bool_t          HLT_CaloMHT90;
   Bool_t          HLT_CaloMET80_NotCleaned;
   Bool_t          HLT_CaloMET90_NotCleaned;
   Bool_t          HLT_CaloMET100_NotCleaned;
   Bool_t          HLT_CaloMET110_NotCleaned;
   Bool_t          HLT_CaloMET250_NotCleaned;
   Bool_t          HLT_CaloMET70_HBHECleaned;
   Bool_t          HLT_CaloMET80_HBHECleaned;
   Bool_t          HLT_CaloMET90_HBHECleaned;
   Bool_t          HLT_CaloMET100_HBHECleaned;
   Bool_t          HLT_CaloMET250_HBHECleaned;
   Bool_t          HLT_CaloMET300_HBHECleaned;
   Bool_t          HLT_CaloMET350_HBHECleaned;
   Bool_t          HLT_PFMET200_NotCleaned;
   Bool_t          HLT_PFMET200_HBHECleaned;
   Bool_t          HLT_PFMET250_HBHECleaned;
   Bool_t          HLT_PFMET300_HBHECleaned;
   Bool_t          HLT_PFMET200_HBHE_BeamHaloCleaned;
   Bool_t          HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned;
   Bool_t          HLT_MET105_IsoTrk50;
   Bool_t          HLT_MET120_IsoTrk50;
   Bool_t          HLT_SingleJet30_Mu12_SinglePFJet40;
   Bool_t          HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33;
   Bool_t          HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33;
   Bool_t          HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33;
   Bool_t          HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33;
   Bool_t          HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33;
   Bool_t          HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33;
   Bool_t          HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33;
   Bool_t          HLT_DoublePFJets40_CaloBTagCSV_p33;
   Bool_t          HLT_DoublePFJets100_CaloBTagCSV_p33;
   Bool_t          HLT_DoublePFJets200_CaloBTagCSV_p33;
   Bool_t          HLT_DoublePFJets350_CaloBTagCSV_p33;
   Bool_t          HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33;
   Bool_t          HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33;
   Bool_t          HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33;
   Bool_t          HLT_Photon300_NoHE;
   Bool_t          HLT_Mu8_TrkIsoVVL;
   Bool_t          HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;
   Bool_t          HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
   Bool_t          HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;
   Bool_t          HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu17_TrkIsoVVL;
   Bool_t          HLT_Mu19_TrkIsoVVL;
   Bool_t          HLT_BTagMu_AK4DiJet20_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet40_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet70_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet110_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet170_Mu5;
   Bool_t          HLT_BTagMu_AK4Jet300_Mu5;
   Bool_t          HLT_BTagMu_AK8DiJet170_Mu5;
   Bool_t          HLT_BTagMu_AK8Jet300_Mu5;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu12_DoublePhoton20;
   Bool_t          HLT_TriplePhoton_20_20_20_CaloIdLV2;
   Bool_t          HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL;
   Bool_t          HLT_TriplePhoton_30_30_10_CaloIdLV2;
   Bool_t          HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL;
   Bool_t          HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL;
   Bool_t          HLT_Photon25;
   Bool_t          HLT_Photon33;
   Bool_t          HLT_Photon50;
   Bool_t          HLT_Photon75;
   Bool_t          HLT_Photon90;
   Bool_t          HLT_Photon120;
   Bool_t          HLT_Photon150;
   Bool_t          HLT_Photon175;
   Bool_t          HLT_Photon200;
   Bool_t          HLT_Photon50_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon90_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon120_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon165_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon90_CaloIdL_PFHT700;
   Bool_t          HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;
   Bool_t          HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95;
   Bool_t          HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;
   Bool_t          HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;
   Bool_t          HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;
   Bool_t          HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;
   Bool_t          HLT_Dimuon0_Jpsi_L1_NoOS;
   Bool_t          HLT_Dimuon0_Jpsi_NoVertexing_NoOS;
   Bool_t          HLT_Dimuon0_Jpsi;
   Bool_t          HLT_Dimuon0_Jpsi_NoVertexing;
   Bool_t          HLT_Dimuon0_Jpsi_L1_4R_0er1p5R;
   Bool_t          HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R;
   Bool_t          HLT_Dimuon0_Jpsi3p5_Muon2;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5;
   Bool_t          HLT_Dimuon0_Upsilon_L1_5;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5NoOS;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5er2p0;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5er2p0M;
   Bool_t          HLT_Dimuon0_Upsilon_NoVertexing;
   Bool_t          HLT_Dimuon0_Upsilon_L1_5M;
   Bool_t          HLT_Dimuon0_LowMass_L1_0er1p5R;
   Bool_t          HLT_Dimuon0_LowMass_L1_0er1p5;
   Bool_t          HLT_Dimuon0_LowMass;
   Bool_t          HLT_Dimuon0_LowMass_L1_4;
   Bool_t          HLT_Dimuon0_LowMass_L1_4R;
   Bool_t          HLT_Dimuon0_LowMass_L1_TM530;
   Bool_t          HLT_Dimuon0_Upsilon_Muon_L1_TM0;
   Bool_t          HLT_Dimuon0_Upsilon_Muon_NoL1Mass;
   Bool_t          HLT_TripleMu_5_3_3_Mass3p8to60_DZ;
   Bool_t          HLT_TripleMu_10_5_5_DZ;
   Bool_t          HLT_TripleMu_12_10_5;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1;
   Bool_t          HLT_DoubleMu3_DZ_PFMET50_PFMHT60;
   Bool_t          HLT_DoubleMu3_DZ_PFMET70_PFMHT70;
   Bool_t          HLT_DoubleMu3_DZ_PFMET90_PFMHT90;
   Bool_t          HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass;
   Bool_t          HLT_DoubleMu4_Jpsi_Displaced;
   Bool_t          HLT_DoubleMu4_Jpsi_NoVertexing;
   Bool_t          HLT_DoubleMu4_JpsiTrkTrk_Displaced;
   Bool_t          HLT_DoubleMu43NoFiltersNoVtx;
   Bool_t          HLT_DoubleMu48NoFiltersNoVtx;
   Bool_t          HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL;
   Bool_t          HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL;
   Bool_t          HLT_DoubleMu20_7_Mass0to30_L1_DM4;
   Bool_t          HLT_DoubleMu20_7_Mass0to30_L1_DM4EG;
   Bool_t          HLT_HT425;
   Bool_t          HLT_HT430_DisplacedDijet40_DisplacedTrack;
   Bool_t          HLT_HT430_DisplacedDijet60_DisplacedTrack;
   Bool_t          HLT_HT430_DisplacedDijet80_DisplacedTrack;
   Bool_t          HLT_HT400_DisplacedDijet40_DisplacedTrack;
   Bool_t          HLT_HT650_DisplacedDijet60_Inclusive;
   Bool_t          HLT_HT550_DisplacedDijet80_Inclusive;
   Bool_t          HLT_HT550_DisplacedDijet60_Inclusive;
   Bool_t          HLT_HT650_DisplacedDijet80_Inclusive;
   Bool_t          HLT_HT750_DisplacedDijet80_Inclusive;
   Bool_t          HLT_DiJet110_35_Mjj650_PFMET110;
   Bool_t          HLT_DiJet110_35_Mjj650_PFMET120;
   Bool_t          HLT_DiJet110_35_Mjj650_PFMET130;
   Bool_t          HLT_TripleJet110_35_35_Mjj650_PFMET110;
   Bool_t          HLT_TripleJet110_35_35_Mjj650_PFMET120;
   Bool_t          HLT_TripleJet110_35_35_Mjj650_PFMET130;
   Bool_t          HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg;
   Bool_t          HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg;
   Bool_t          HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg;
   Bool_t          HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned;
   Bool_t          HLT_Ele28_eta2p1_WPTight_Gsf_HT150;
   Bool_t          HLT_Ele28_HighEta_SC20_Mass55;
   Bool_t          HLT_DoubleMu20_7_Mass0to30_Photon23;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT450_PFMET50;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT450;
   Bool_t          HLT_Ele50_IsoVVVL_PFHT450;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT600;
   Bool_t          HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;
   Bool_t          HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT450_PFMET50;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT450;
   Bool_t          HLT_Mu50_IsoVVVL_PFHT450;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT600;
   Bool_t          HLT_Dimuon10_PsiPrime_Barrel_Seagulls;
   Bool_t          HLT_Dimuon20_Jpsi_Barrel_Seagulls;
   Bool_t          HLT_Dimuon10_Upsilon_Barrel_Seagulls;
   Bool_t          HLT_Dimuon12_Upsilon_eta1p5;
   Bool_t          HLT_Dimuon14_Phi_Barrel_Seagulls;
   Bool_t          HLT_Dimuon18_PsiPrime;
   Bool_t          HLT_Dimuon25_Jpsi;
   Bool_t          HLT_Dimuon18_PsiPrime_noCorrL1;
   Bool_t          HLT_Dimuon24_Upsilon_noCorrL1;
   Bool_t          HLT_Dimuon24_Phi_noCorrL1;
   Bool_t          HLT_Dimuon25_Jpsi_noCorrL1;
   Bool_t          HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
   Bool_t          HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
   Bool_t          HLT_DoubleIsoMu20_eta2p1;
   Bool_t          HLT_DoubleIsoMu24_eta2p1;
   Bool_t          HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx;
   Bool_t          HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx;
   Bool_t          HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx;
   Bool_t          HLT_Mu8;
   Bool_t          HLT_Mu17;
   Bool_t          HLT_Mu19;
   Bool_t          HLT_Mu17_Photon30_IsoCaloId;
   Bool_t          HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele8_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Ele17_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Ele23_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
   Bool_t          HLT_Ele115_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele135_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele145_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele200_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele250_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele300_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_PFHT300PT30_QuadPFJet_75_60_45_40;
   Bool_t          HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0;
   Bool_t          HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2;
   Bool_t          HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2;
   Bool_t          HLT_PFHT380_SixPFJet32;
   Bool_t          HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5;
   Bool_t          HLT_PFHT430_SixPFJet40;
   Bool_t          HLT_PFHT350;
   Bool_t          HLT_PFHT350MinPFJet15;
   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL;
   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL;
   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15;
   Bool_t          HLT_FullTrack_Multiplicity85;
   Bool_t          HLT_FullTrack_Multiplicity100;
   Bool_t          HLT_FullTrack_Multiplicity130;
   Bool_t          HLT_FullTrack_Multiplicity155;
   Bool_t          HLT_ECALHT800;
   Bool_t          HLT_DiSC30_18_EIso_AND_HE_Mass70;
   Bool_t          HLT_Physics;
   Bool_t          HLT_Physics_part0;
   Bool_t          HLT_Physics_part1;
   Bool_t          HLT_Physics_part2;
   Bool_t          HLT_Physics_part3;
   Bool_t          HLT_Physics_part4;
   Bool_t          HLT_Physics_part5;
   Bool_t          HLT_Physics_part6;
   Bool_t          HLT_Physics_part7;
   Bool_t          HLT_Random;
   Bool_t          HLT_ZeroBias;
   Bool_t          HLT_ZeroBias_part0;
   Bool_t          HLT_ZeroBias_part1;
   Bool_t          HLT_ZeroBias_part2;
   Bool_t          HLT_ZeroBias_part3;
   Bool_t          HLT_ZeroBias_part4;
   Bool_t          HLT_ZeroBias_part5;
   Bool_t          HLT_ZeroBias_part6;
   Bool_t          HLT_ZeroBias_part7;
   Bool_t          HLT_AK4CaloJet30;
   Bool_t          HLT_AK4CaloJet40;
   Bool_t          HLT_AK4CaloJet50;
   Bool_t          HLT_AK4CaloJet80;
   Bool_t          HLT_AK4CaloJet100;
   Bool_t          HLT_AK4CaloJet120;
   Bool_t          HLT_AK4PFJet30;
   Bool_t          HLT_AK4PFJet50;
   Bool_t          HLT_AK4PFJet80;
   Bool_t          HLT_AK4PFJet100;
   Bool_t          HLT_AK4PFJet120;
   Bool_t          HLT_HISinglePhoton10_Eta3p1ForPPRef;
   Bool_t          HLT_HISinglePhoton20_Eta3p1ForPPRef;
   Bool_t          HLT_HISinglePhoton30_Eta3p1ForPPRef;
   Bool_t          HLT_HISinglePhoton40_Eta3p1ForPPRef;
   Bool_t          HLT_HISinglePhoton50_Eta3p1ForPPRef;
   Bool_t          HLT_HISinglePhoton60_Eta3p1ForPPRef;
   Bool_t          HLT_Photon20_HoverELoose;
   Bool_t          HLT_Photon30_HoverELoose;
   Bool_t          HLT_Photon40_HoverELoose;
   Bool_t          HLT_Photon50_HoverELoose;
   Bool_t          HLT_Photon60_HoverELoose;
   Bool_t          HLT_EcalCalibration;
   Bool_t          HLT_HcalCalibration;
   Bool_t          HLT_L1UnpairedBunchBptxMinus;
   Bool_t          HLT_L1UnpairedBunchBptxPlus;
   Bool_t          HLT_L1NotBptxOR;
   Bool_t          HLT_L1MinimumBiasHF_OR;
   Bool_t          HLT_L1MinimumBiasHF0OR;
   Bool_t          HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;
   Bool_t          HLT_HcalNZS;
   Bool_t          HLT_HcalPhiSym;
   Bool_t          HLT_HcalIsolatedbunch;
   Bool_t          HLT_IsoTrackHB;
   Bool_t          HLT_IsoTrackHE;
   Bool_t          HLT_ZeroBias_FirstCollisionAfterAbortGap;
   Bool_t          HLT_ZeroBias_IsolatedBunches;
   Bool_t          HLT_ZeroBias_FirstCollisionInTrain;
   Bool_t          HLT_ZeroBias_LastCollisionInTrain;
   Bool_t          HLT_ZeroBias_FirstBXAfterTrain;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1;
   Bool_t          HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg;
   Bool_t          HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130;
   Bool_t          HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;
   Bool_t          HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr;
   Bool_t          HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1;
   Bool_t          HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1;
   Bool_t          HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
   Bool_t          HLT_Rsq0p35;
   Bool_t          HLT_Rsq0p40;
   Bool_t          HLT_RsqMR300_Rsq0p09_MR200;
   Bool_t          HLT_RsqMR320_Rsq0p09_MR200;
   Bool_t          HLT_RsqMR300_Rsq0p09_MR200_4jet;
   Bool_t          HLT_RsqMR320_Rsq0p09_MR200_4jet;
   Bool_t          HLT_L1_DoubleJet30_Mass_Min400_Mu10;
   Bool_t          HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1;
   Bool_t          HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1;
   Bool_t          HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1;
   Bool_t          HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3;
   Bool_t          HLT_PFMET100_PFMHT100_IDTight_PFHT60;
   Bool_t          HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;
   Bool_t          HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60;
   Bool_t          HLT_Mu18_Mu9_SameSign;
   Bool_t          HLT_Mu18_Mu9_SameSign_DZ;
   Bool_t          HLT_Mu18_Mu9;
   Bool_t          HLT_Mu18_Mu9_DZ;
   Bool_t          HLT_Mu20_Mu10_SameSign;
   Bool_t          HLT_Mu20_Mu10_SameSign_DZ;
   Bool_t          HLT_Mu20_Mu10;
   Bool_t          HLT_Mu20_Mu10_DZ;
   Bool_t          HLT_Mu23_Mu12_SameSign;
   Bool_t          HLT_Mu23_Mu12_SameSign_DZ;
   Bool_t          HLT_Mu23_Mu12;
   Bool_t          HLT_Mu23_Mu12_DZ;
   Bool_t          HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi;
   Bool_t          HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi;
   Bool_t          HLT_DoubleMu3_DCA_PFMET50_PFMHT60;
   Bool_t          HLT_TripleMu_5_3_3_Mass3p8to60_DCA;
   Bool_t          HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1;
   Bool_t          HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1;
   Bool_t          HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1;
   Bool_t          HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1;
   Bool_t          HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2;
   Bool_t          HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2;
   Bool_t          HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2;
   Bool_t          HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2;
   Bool_t          HLT_QuadPFJet98_83_71_15;
   Bool_t          HLT_QuadPFJet103_88_75_15;
   Bool_t          HLT_QuadPFJet105_88_76_15;
   Bool_t          HLT_QuadPFJet111_90_80_15;
   Bool_t          HLT_AK8PFJet330_PFAK8BTagCSV_p17;
   Bool_t          HLT_AK8PFJet330_PFAK8BTagCSV_p1;
   Bool_t          HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;
   Bool_t          HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;
   Bool_t          HLTriggerFinalPath;
   Bool_t          L1simulation_step;
   Float_t         Jet_pt_raw[30];   //[nJet]
   Float_t         Jet_pt_nom[30];   //[nJet]
   Float_t         Jet_mass_raw[30];   //[nJet]
   Float_t         Jet_mass_nom[30];   //[nJet]
   Float_t         Jet_corr_JEC[30];   //[nJet]
   Float_t         Jet_corr_JER[30];   //[nJet]
   Float_t         MET_T1_pt;
   Float_t         MET_T1_phi;
   Float_t         MET_T1Smear_pt;
   Float_t         MET_T1Smear_phi;
   Float_t         Jet_pt_jer0Up[30];   //[nJet]
   Float_t         Jet_mass_jer0Up[30];   //[nJet]
   Float_t         MET_T1_pt_jer0Up;
   Float_t         MET_T1_phi_jer0Up;
   Float_t         MET_T1Smear_pt_jer0Up;
   Float_t         MET_T1Smear_phi_jer0Up;
   Float_t         Jet_pt_jer1Up[30];   //[nJet]
   Float_t         Jet_mass_jer1Up[30];   //[nJet]
   Float_t         MET_T1_pt_jer1Up;
   Float_t         MET_T1_phi_jer1Up;
   Float_t         MET_T1Smear_pt_jer1Up;
   Float_t         MET_T1Smear_phi_jer1Up;
   Float_t         Jet_pt_jer2Up[30];   //[nJet]
   Float_t         Jet_mass_jer2Up[30];   //[nJet]
   Float_t         MET_T1_pt_jer2Up;
   Float_t         MET_T1_phi_jer2Up;
   Float_t         MET_T1Smear_pt_jer2Up;
   Float_t         MET_T1Smear_phi_jer2Up;
   Float_t         Jet_pt_jer3Up[30];   //[nJet]
   Float_t         Jet_mass_jer3Up[30];   //[nJet]
   Float_t         MET_T1_pt_jer3Up;
   Float_t         MET_T1_phi_jer3Up;
   Float_t         MET_T1Smear_pt_jer3Up;
   Float_t         MET_T1Smear_phi_jer3Up;
   Float_t         Jet_pt_jer4Up[30];   //[nJet]
   Float_t         Jet_mass_jer4Up[30];   //[nJet]
   Float_t         MET_T1_pt_jer4Up;
   Float_t         MET_T1_phi_jer4Up;
   Float_t         MET_T1Smear_pt_jer4Up;
   Float_t         MET_T1Smear_phi_jer4Up;
   Float_t         Jet_pt_jer5Up[30];   //[nJet]
   Float_t         Jet_mass_jer5Up[30];   //[nJet]
   Float_t         MET_T1_pt_jer5Up;
   Float_t         MET_T1_phi_jer5Up;
   Float_t         MET_T1Smear_pt_jer5Up;
   Float_t         MET_T1Smear_phi_jer5Up;
   Float_t         Jet_pt_jesAbsoluteStatUp[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteStatUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesAbsoluteStatUp;
   Float_t         MET_T1_phi_jesAbsoluteStatUp;
   Float_t         MET_T1Smear_pt_jesAbsoluteStatUp;
   Float_t         MET_T1Smear_phi_jesAbsoluteStatUp;
   Float_t         Jet_pt_jesAbsoluteScaleUp[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteScaleUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesAbsoluteScaleUp;
   Float_t         MET_T1_phi_jesAbsoluteScaleUp;
   Float_t         MET_T1Smear_pt_jesAbsoluteScaleUp;
   Float_t         MET_T1Smear_phi_jesAbsoluteScaleUp;
   Float_t         Jet_pt_jesAbsoluteFlavMapUp[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteFlavMapUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesAbsoluteFlavMapUp;
   Float_t         MET_T1_phi_jesAbsoluteFlavMapUp;
   Float_t         MET_T1Smear_pt_jesAbsoluteFlavMapUp;
   Float_t         MET_T1Smear_phi_jesAbsoluteFlavMapUp;
   Float_t         Jet_pt_jesAbsoluteMPFBiasUp[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteMPFBiasUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesAbsoluteMPFBiasUp;
   Float_t         MET_T1_phi_jesAbsoluteMPFBiasUp;
   Float_t         MET_T1Smear_pt_jesAbsoluteMPFBiasUp;
   Float_t         MET_T1Smear_phi_jesAbsoluteMPFBiasUp;
   Float_t         Jet_pt_jesFragmentationUp[30];   //[nJet]
   Float_t         Jet_mass_jesFragmentationUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesFragmentationUp;
   Float_t         MET_T1_phi_jesFragmentationUp;
   Float_t         MET_T1Smear_pt_jesFragmentationUp;
   Float_t         MET_T1Smear_phi_jesFragmentationUp;
   Float_t         Jet_pt_jesSinglePionECALUp[30];   //[nJet]
   Float_t         Jet_mass_jesSinglePionECALUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesSinglePionECALUp;
   Float_t         MET_T1_phi_jesSinglePionECALUp;
   Float_t         MET_T1Smear_pt_jesSinglePionECALUp;
   Float_t         MET_T1Smear_phi_jesSinglePionECALUp;
   Float_t         Jet_pt_jesSinglePionHCALUp[30];   //[nJet]
   Float_t         Jet_mass_jesSinglePionHCALUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesSinglePionHCALUp;
   Float_t         MET_T1_phi_jesSinglePionHCALUp;
   Float_t         MET_T1Smear_pt_jesSinglePionHCALUp;
   Float_t         MET_T1Smear_phi_jesSinglePionHCALUp;
   Float_t         Jet_pt_jesFlavorQCDUp[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorQCDUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesFlavorQCDUp;
   Float_t         MET_T1_phi_jesFlavorQCDUp;
   Float_t         MET_T1Smear_pt_jesFlavorQCDUp;
   Float_t         MET_T1Smear_phi_jesFlavorQCDUp;
   Float_t         Jet_pt_jesTimePtEtaUp[30];   //[nJet]
   Float_t         Jet_mass_jesTimePtEtaUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesTimePtEtaUp;
   Float_t         MET_T1_phi_jesTimePtEtaUp;
   Float_t         MET_T1Smear_pt_jesTimePtEtaUp;
   Float_t         MET_T1Smear_phi_jesTimePtEtaUp;
   Float_t         Jet_pt_jesRelativeJEREC1Up[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeJEREC1Up[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeJEREC1Up;
   Float_t         MET_T1_phi_jesRelativeJEREC1Up;
   Float_t         MET_T1Smear_pt_jesRelativeJEREC1Up;
   Float_t         MET_T1Smear_phi_jesRelativeJEREC1Up;
   Float_t         Jet_pt_jesRelativeJEREC2Up[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeJEREC2Up[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeJEREC2Up;
   Float_t         MET_T1_phi_jesRelativeJEREC2Up;
   Float_t         MET_T1Smear_pt_jesRelativeJEREC2Up;
   Float_t         MET_T1Smear_phi_jesRelativeJEREC2Up;
   Float_t         Jet_pt_jesRelativeJERHFUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeJERHFUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeJERHFUp;
   Float_t         MET_T1_phi_jesRelativeJERHFUp;
   Float_t         MET_T1Smear_pt_jesRelativeJERHFUp;
   Float_t         MET_T1Smear_phi_jesRelativeJERHFUp;
   Float_t         Jet_pt_jesRelativePtBBUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtBBUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativePtBBUp;
   Float_t         MET_T1_phi_jesRelativePtBBUp;
   Float_t         MET_T1Smear_pt_jesRelativePtBBUp;
   Float_t         MET_T1Smear_phi_jesRelativePtBBUp;
   Float_t         Jet_pt_jesRelativePtEC1Up[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtEC1Up[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativePtEC1Up;
   Float_t         MET_T1_phi_jesRelativePtEC1Up;
   Float_t         MET_T1Smear_pt_jesRelativePtEC1Up;
   Float_t         MET_T1Smear_phi_jesRelativePtEC1Up;
   Float_t         Jet_pt_jesRelativePtEC2Up[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtEC2Up[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativePtEC2Up;
   Float_t         MET_T1_phi_jesRelativePtEC2Up;
   Float_t         MET_T1Smear_pt_jesRelativePtEC2Up;
   Float_t         MET_T1Smear_phi_jesRelativePtEC2Up;
   Float_t         Jet_pt_jesRelativePtHFUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtHFUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativePtHFUp;
   Float_t         MET_T1_phi_jesRelativePtHFUp;
   Float_t         MET_T1Smear_pt_jesRelativePtHFUp;
   Float_t         MET_T1Smear_phi_jesRelativePtHFUp;
   Float_t         Jet_pt_jesRelativeBalUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeBalUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeBalUp;
   Float_t         MET_T1_phi_jesRelativeBalUp;
   Float_t         MET_T1Smear_pt_jesRelativeBalUp;
   Float_t         MET_T1Smear_phi_jesRelativeBalUp;
   Float_t         Jet_pt_jesRelativeSampleUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeSampleUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeSampleUp;
   Float_t         MET_T1_phi_jesRelativeSampleUp;
   Float_t         MET_T1Smear_pt_jesRelativeSampleUp;
   Float_t         MET_T1Smear_phi_jesRelativeSampleUp;
   Float_t         Jet_pt_jesRelativeFSRUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeFSRUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeFSRUp;
   Float_t         MET_T1_phi_jesRelativeFSRUp;
   Float_t         MET_T1Smear_pt_jesRelativeFSRUp;
   Float_t         MET_T1Smear_phi_jesRelativeFSRUp;
   Float_t         Jet_pt_jesRelativeStatFSRUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatFSRUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeStatFSRUp;
   Float_t         MET_T1_phi_jesRelativeStatFSRUp;
   Float_t         MET_T1Smear_pt_jesRelativeStatFSRUp;
   Float_t         MET_T1Smear_phi_jesRelativeStatFSRUp;
   Float_t         Jet_pt_jesRelativeStatECUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatECUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeStatECUp;
   Float_t         MET_T1_phi_jesRelativeStatECUp;
   Float_t         MET_T1Smear_pt_jesRelativeStatECUp;
   Float_t         MET_T1Smear_phi_jesRelativeStatECUp;
   Float_t         Jet_pt_jesRelativeStatHFUp[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatHFUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeStatHFUp;
   Float_t         MET_T1_phi_jesRelativeStatHFUp;
   Float_t         MET_T1Smear_pt_jesRelativeStatHFUp;
   Float_t         MET_T1Smear_phi_jesRelativeStatHFUp;
   Float_t         Jet_pt_jesPileUpDataMCUp[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpDataMCUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpDataMCUp;
   Float_t         MET_T1_phi_jesPileUpDataMCUp;
   Float_t         MET_T1Smear_pt_jesPileUpDataMCUp;
   Float_t         MET_T1Smear_phi_jesPileUpDataMCUp;
   Float_t         Jet_pt_jesPileUpPtRefUp[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtRefUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpPtRefUp;
   Float_t         MET_T1_phi_jesPileUpPtRefUp;
   Float_t         MET_T1Smear_pt_jesPileUpPtRefUp;
   Float_t         MET_T1Smear_phi_jesPileUpPtRefUp;
   Float_t         Jet_pt_jesPileUpPtBBUp[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtBBUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpPtBBUp;
   Float_t         MET_T1_phi_jesPileUpPtBBUp;
   Float_t         MET_T1Smear_pt_jesPileUpPtBBUp;
   Float_t         MET_T1Smear_phi_jesPileUpPtBBUp;
   Float_t         Jet_pt_jesPileUpPtEC1Up[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtEC1Up[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpPtEC1Up;
   Float_t         MET_T1_phi_jesPileUpPtEC1Up;
   Float_t         MET_T1Smear_pt_jesPileUpPtEC1Up;
   Float_t         MET_T1Smear_phi_jesPileUpPtEC1Up;
   Float_t         Jet_pt_jesPileUpPtEC2Up[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtEC2Up[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpPtEC2Up;
   Float_t         MET_T1_phi_jesPileUpPtEC2Up;
   Float_t         MET_T1Smear_pt_jesPileUpPtEC2Up;
   Float_t         MET_T1Smear_phi_jesPileUpPtEC2Up;
   Float_t         Jet_pt_jesPileUpPtHFUp[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtHFUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpPtHFUp;
   Float_t         MET_T1_phi_jesPileUpPtHFUp;
   Float_t         MET_T1Smear_pt_jesPileUpPtHFUp;
   Float_t         MET_T1Smear_phi_jesPileUpPtHFUp;
   Float_t         Jet_pt_jesPileUpMuZeroUp[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpMuZeroUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpMuZeroUp;
   Float_t         MET_T1_phi_jesPileUpMuZeroUp;
   Float_t         MET_T1Smear_pt_jesPileUpMuZeroUp;
   Float_t         MET_T1Smear_phi_jesPileUpMuZeroUp;
   Float_t         Jet_pt_jesPileUpEnvelopeUp[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpEnvelopeUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpEnvelopeUp;
   Float_t         MET_T1_phi_jesPileUpEnvelopeUp;
   Float_t         MET_T1Smear_pt_jesPileUpEnvelopeUp;
   Float_t         MET_T1Smear_phi_jesPileUpEnvelopeUp;
   Float_t         Jet_pt_jesSubTotalPileUpUp[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalPileUpUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesSubTotalPileUpUp;
   Float_t         MET_T1_phi_jesSubTotalPileUpUp;
   Float_t         MET_T1Smear_pt_jesSubTotalPileUpUp;
   Float_t         MET_T1Smear_phi_jesSubTotalPileUpUp;
   Float_t         Jet_pt_jesSubTotalRelativeUp[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalRelativeUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesSubTotalRelativeUp;
   Float_t         MET_T1_phi_jesSubTotalRelativeUp;
   Float_t         MET_T1Smear_pt_jesSubTotalRelativeUp;
   Float_t         MET_T1Smear_phi_jesSubTotalRelativeUp;
   Float_t         Jet_pt_jesSubTotalPtUp[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalPtUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesSubTotalPtUp;
   Float_t         MET_T1_phi_jesSubTotalPtUp;
   Float_t         MET_T1Smear_pt_jesSubTotalPtUp;
   Float_t         MET_T1Smear_phi_jesSubTotalPtUp;
   Float_t         Jet_pt_jesSubTotalScaleUp[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalScaleUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesSubTotalScaleUp;
   Float_t         MET_T1_phi_jesSubTotalScaleUp;
   Float_t         MET_T1Smear_pt_jesSubTotalScaleUp;
   Float_t         MET_T1Smear_phi_jesSubTotalScaleUp;
   Float_t         Jet_pt_jesSubTotalAbsoluteUp[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalAbsoluteUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesSubTotalAbsoluteUp;
   Float_t         MET_T1_phi_jesSubTotalAbsoluteUp;
   Float_t         MET_T1Smear_pt_jesSubTotalAbsoluteUp;
   Float_t         MET_T1Smear_phi_jesSubTotalAbsoluteUp;
   Float_t         Jet_pt_jesSubTotalMCUp[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalMCUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesSubTotalMCUp;
   Float_t         MET_T1_phi_jesSubTotalMCUp;
   Float_t         MET_T1Smear_pt_jesSubTotalMCUp;
   Float_t         MET_T1Smear_phi_jesSubTotalMCUp;
   Float_t         Jet_pt_jesTotalUp[30];   //[nJet]
   Float_t         Jet_mass_jesTotalUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesTotalUp;
   Float_t         MET_T1_phi_jesTotalUp;
   Float_t         MET_T1Smear_pt_jesTotalUp;
   Float_t         MET_T1Smear_phi_jesTotalUp;
   Float_t         Jet_pt_jesTotalNoFlavorUp[30];   //[nJet]
   Float_t         Jet_mass_jesTotalNoFlavorUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesTotalNoFlavorUp;
   Float_t         MET_T1_phi_jesTotalNoFlavorUp;
   Float_t         MET_T1Smear_pt_jesTotalNoFlavorUp;
   Float_t         MET_T1Smear_phi_jesTotalNoFlavorUp;
   Float_t         Jet_pt_jesTotalNoTimeUp[30];   //[nJet]
   Float_t         Jet_mass_jesTotalNoTimeUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesTotalNoTimeUp;
   Float_t         MET_T1_phi_jesTotalNoTimeUp;
   Float_t         MET_T1Smear_pt_jesTotalNoTimeUp;
   Float_t         MET_T1Smear_phi_jesTotalNoTimeUp;
   Float_t         Jet_pt_jesTotalNoFlavorNoTimeUp[30];   //[nJet]
   Float_t         Jet_mass_jesTotalNoFlavorNoTimeUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesTotalNoFlavorNoTimeUp;
   Float_t         MET_T1_phi_jesTotalNoFlavorNoTimeUp;
   Float_t         MET_T1Smear_pt_jesTotalNoFlavorNoTimeUp;
   Float_t         MET_T1Smear_phi_jesTotalNoFlavorNoTimeUp;
   Float_t         Jet_pt_jesFlavorZJetUp[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorZJetUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesFlavorZJetUp;
   Float_t         MET_T1_phi_jesFlavorZJetUp;
   Float_t         MET_T1Smear_pt_jesFlavorZJetUp;
   Float_t         MET_T1Smear_phi_jesFlavorZJetUp;
   Float_t         Jet_pt_jesFlavorPhotonJetUp[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPhotonJetUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesFlavorPhotonJetUp;
   Float_t         MET_T1_phi_jesFlavorPhotonJetUp;
   Float_t         MET_T1Smear_pt_jesFlavorPhotonJetUp;
   Float_t         MET_T1Smear_phi_jesFlavorPhotonJetUp;
   Float_t         Jet_pt_jesFlavorPureGluonUp[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureGluonUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesFlavorPureGluonUp;
   Float_t         MET_T1_phi_jesFlavorPureGluonUp;
   Float_t         MET_T1Smear_pt_jesFlavorPureGluonUp;
   Float_t         MET_T1Smear_phi_jesFlavorPureGluonUp;
   Float_t         Jet_pt_jesFlavorPureQuarkUp[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureQuarkUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesFlavorPureQuarkUp;
   Float_t         MET_T1_phi_jesFlavorPureQuarkUp;
   Float_t         MET_T1Smear_pt_jesFlavorPureQuarkUp;
   Float_t         MET_T1Smear_phi_jesFlavorPureQuarkUp;
   Float_t         Jet_pt_jesFlavorPureCharmUp[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureCharmUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesFlavorPureCharmUp;
   Float_t         MET_T1_phi_jesFlavorPureCharmUp;
   Float_t         MET_T1Smear_pt_jesFlavorPureCharmUp;
   Float_t         MET_T1Smear_phi_jesFlavorPureCharmUp;
   Float_t         Jet_pt_jesFlavorPureBottomUp[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureBottomUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesFlavorPureBottomUp;
   Float_t         MET_T1_phi_jesFlavorPureBottomUp;
   Float_t         MET_T1Smear_pt_jesFlavorPureBottomUp;
   Float_t         MET_T1Smear_phi_jesFlavorPureBottomUp;
   Float_t         Jet_pt_jesTimeRunBUp[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunBUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesTimeRunBUp;
   Float_t         MET_T1_phi_jesTimeRunBUp;
   Float_t         MET_T1Smear_pt_jesTimeRunBUp;
   Float_t         MET_T1Smear_phi_jesTimeRunBUp;
   Float_t         Jet_pt_jesTimeRunCUp[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunCUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesTimeRunCUp;
   Float_t         MET_T1_phi_jesTimeRunCUp;
   Float_t         MET_T1Smear_pt_jesTimeRunCUp;
   Float_t         MET_T1Smear_phi_jesTimeRunCUp;
   Float_t         Jet_pt_jesTimeRunDEUp[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunDEUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesTimeRunDEUp;
   Float_t         MET_T1_phi_jesTimeRunDEUp;
   Float_t         MET_T1Smear_pt_jesTimeRunDEUp;
   Float_t         MET_T1Smear_phi_jesTimeRunDEUp;
   Float_t         Jet_pt_jesTimeRunFUp[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunFUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesTimeRunFUp;
   Float_t         MET_T1_phi_jesTimeRunFUp;
   Float_t         MET_T1Smear_pt_jesTimeRunFUp;
   Float_t         MET_T1Smear_phi_jesTimeRunFUp;
   Float_t         Jet_pt_jesCorrelationGroupMPFInSituUp[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupMPFInSituUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesCorrelationGroupMPFInSituUp;
   Float_t         MET_T1_phi_jesCorrelationGroupMPFInSituUp;
   Float_t         MET_T1Smear_pt_jesCorrelationGroupMPFInSituUp;
   Float_t         MET_T1Smear_phi_jesCorrelationGroupMPFInSituUp;
   Float_t         Jet_pt_jesCorrelationGroupIntercalibrationUp[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupIntercalibrationUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesCorrelationGroupIntercalibrationUp;
   Float_t         MET_T1_phi_jesCorrelationGroupIntercalibrationUp;
   Float_t         MET_T1Smear_pt_jesCorrelationGroupIntercalibrationUp;
   Float_t         MET_T1Smear_phi_jesCorrelationGroupIntercalibrationUp;
   Float_t         Jet_pt_jesCorrelationGroupbJESUp[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupbJESUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesCorrelationGroupbJESUp;
   Float_t         MET_T1_phi_jesCorrelationGroupbJESUp;
   Float_t         MET_T1Smear_pt_jesCorrelationGroupbJESUp;
   Float_t         MET_T1Smear_phi_jesCorrelationGroupbJESUp;
   Float_t         Jet_pt_jesCorrelationGroupFlavorUp[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupFlavorUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesCorrelationGroupFlavorUp;
   Float_t         MET_T1_phi_jesCorrelationGroupFlavorUp;
   Float_t         MET_T1Smear_pt_jesCorrelationGroupFlavorUp;
   Float_t         MET_T1Smear_phi_jesCorrelationGroupFlavorUp;
   Float_t         Jet_pt_jesCorrelationGroupUncorrelatedUp[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupUncorrelatedUp[30];   //[nJet]
   Float_t         MET_T1_pt_jesCorrelationGroupUncorrelatedUp;
   Float_t         MET_T1_phi_jesCorrelationGroupUncorrelatedUp;
   Float_t         MET_T1Smear_pt_jesCorrelationGroupUncorrelatedUp;
   Float_t         MET_T1Smear_phi_jesCorrelationGroupUncorrelatedUp;
   Float_t         MET_T1_pt_unclustEnUp;
   Float_t         MET_T1_phi_unclustEnUp;
   Float_t         MET_T1Smear_pt_unclustEnUp;
   Float_t         MET_T1Smear_phi_unclustEnUp;
   Float_t         Jet_pt_jer0Down[30];   //[nJet]
   Float_t         Jet_mass_jer0Down[30];   //[nJet]
   Float_t         MET_T1_pt_jer0Down;
   Float_t         MET_T1_phi_jer0Down;
   Float_t         MET_T1Smear_pt_jer0Down;
   Float_t         MET_T1Smear_phi_jer0Down;
   Float_t         Jet_pt_jer1Down[30];   //[nJet]
   Float_t         Jet_mass_jer1Down[30];   //[nJet]
   Float_t         MET_T1_pt_jer1Down;
   Float_t         MET_T1_phi_jer1Down;
   Float_t         MET_T1Smear_pt_jer1Down;
   Float_t         MET_T1Smear_phi_jer1Down;
   Float_t         Jet_pt_jer2Down[30];   //[nJet]
   Float_t         Jet_mass_jer2Down[30];   //[nJet]
   Float_t         MET_T1_pt_jer2Down;
   Float_t         MET_T1_phi_jer2Down;
   Float_t         MET_T1Smear_pt_jer2Down;
   Float_t         MET_T1Smear_phi_jer2Down;
   Float_t         Jet_pt_jer3Down[30];   //[nJet]
   Float_t         Jet_mass_jer3Down[30];   //[nJet]
   Float_t         MET_T1_pt_jer3Down;
   Float_t         MET_T1_phi_jer3Down;
   Float_t         MET_T1Smear_pt_jer3Down;
   Float_t         MET_T1Smear_phi_jer3Down;
   Float_t         Jet_pt_jer4Down[30];   //[nJet]
   Float_t         Jet_mass_jer4Down[30];   //[nJet]
   Float_t         MET_T1_pt_jer4Down;
   Float_t         MET_T1_phi_jer4Down;
   Float_t         MET_T1Smear_pt_jer4Down;
   Float_t         MET_T1Smear_phi_jer4Down;
   Float_t         Jet_pt_jer5Down[30];   //[nJet]
   Float_t         Jet_mass_jer5Down[30];   //[nJet]
   Float_t         MET_T1_pt_jer5Down;
   Float_t         MET_T1_phi_jer5Down;
   Float_t         MET_T1Smear_pt_jer5Down;
   Float_t         MET_T1Smear_phi_jer5Down;
   Float_t         Jet_pt_jesAbsoluteStatDown[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteStatDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesAbsoluteStatDown;
   Float_t         MET_T1_phi_jesAbsoluteStatDown;
   Float_t         MET_T1Smear_pt_jesAbsoluteStatDown;
   Float_t         MET_T1Smear_phi_jesAbsoluteStatDown;
   Float_t         Jet_pt_jesAbsoluteScaleDown[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteScaleDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesAbsoluteScaleDown;
   Float_t         MET_T1_phi_jesAbsoluteScaleDown;
   Float_t         MET_T1Smear_pt_jesAbsoluteScaleDown;
   Float_t         MET_T1Smear_phi_jesAbsoluteScaleDown;
   Float_t         Jet_pt_jesAbsoluteFlavMapDown[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteFlavMapDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesAbsoluteFlavMapDown;
   Float_t         MET_T1_phi_jesAbsoluteFlavMapDown;
   Float_t         MET_T1Smear_pt_jesAbsoluteFlavMapDown;
   Float_t         MET_T1Smear_phi_jesAbsoluteFlavMapDown;
   Float_t         Jet_pt_jesAbsoluteMPFBiasDown[30];   //[nJet]
   Float_t         Jet_mass_jesAbsoluteMPFBiasDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesAbsoluteMPFBiasDown;
   Float_t         MET_T1_phi_jesAbsoluteMPFBiasDown;
   Float_t         MET_T1Smear_pt_jesAbsoluteMPFBiasDown;
   Float_t         MET_T1Smear_phi_jesAbsoluteMPFBiasDown;
   Float_t         Jet_pt_jesFragmentationDown[30];   //[nJet]
   Float_t         Jet_mass_jesFragmentationDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesFragmentationDown;
   Float_t         MET_T1_phi_jesFragmentationDown;
   Float_t         MET_T1Smear_pt_jesFragmentationDown;
   Float_t         MET_T1Smear_phi_jesFragmentationDown;
   Float_t         Jet_pt_jesSinglePionECALDown[30];   //[nJet]
   Float_t         Jet_mass_jesSinglePionECALDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesSinglePionECALDown;
   Float_t         MET_T1_phi_jesSinglePionECALDown;
   Float_t         MET_T1Smear_pt_jesSinglePionECALDown;
   Float_t         MET_T1Smear_phi_jesSinglePionECALDown;
   Float_t         Jet_pt_jesSinglePionHCALDown[30];   //[nJet]
   Float_t         Jet_mass_jesSinglePionHCALDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesSinglePionHCALDown;
   Float_t         MET_T1_phi_jesSinglePionHCALDown;
   Float_t         MET_T1Smear_pt_jesSinglePionHCALDown;
   Float_t         MET_T1Smear_phi_jesSinglePionHCALDown;
   Float_t         Jet_pt_jesFlavorQCDDown[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorQCDDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesFlavorQCDDown;
   Float_t         MET_T1_phi_jesFlavorQCDDown;
   Float_t         MET_T1Smear_pt_jesFlavorQCDDown;
   Float_t         MET_T1Smear_phi_jesFlavorQCDDown;
   Float_t         Jet_pt_jesTimePtEtaDown[30];   //[nJet]
   Float_t         Jet_mass_jesTimePtEtaDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesTimePtEtaDown;
   Float_t         MET_T1_phi_jesTimePtEtaDown;
   Float_t         MET_T1Smear_pt_jesTimePtEtaDown;
   Float_t         MET_T1Smear_phi_jesTimePtEtaDown;
   Float_t         Jet_pt_jesRelativeJEREC1Down[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeJEREC1Down[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeJEREC1Down;
   Float_t         MET_T1_phi_jesRelativeJEREC1Down;
   Float_t         MET_T1Smear_pt_jesRelativeJEREC1Down;
   Float_t         MET_T1Smear_phi_jesRelativeJEREC1Down;
   Float_t         Jet_pt_jesRelativeJEREC2Down[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeJEREC2Down[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeJEREC2Down;
   Float_t         MET_T1_phi_jesRelativeJEREC2Down;
   Float_t         MET_T1Smear_pt_jesRelativeJEREC2Down;
   Float_t         MET_T1Smear_phi_jesRelativeJEREC2Down;
   Float_t         Jet_pt_jesRelativeJERHFDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeJERHFDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeJERHFDown;
   Float_t         MET_T1_phi_jesRelativeJERHFDown;
   Float_t         MET_T1Smear_pt_jesRelativeJERHFDown;
   Float_t         MET_T1Smear_phi_jesRelativeJERHFDown;
   Float_t         Jet_pt_jesRelativePtBBDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtBBDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativePtBBDown;
   Float_t         MET_T1_phi_jesRelativePtBBDown;
   Float_t         MET_T1Smear_pt_jesRelativePtBBDown;
   Float_t         MET_T1Smear_phi_jesRelativePtBBDown;
   Float_t         Jet_pt_jesRelativePtEC1Down[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtEC1Down[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativePtEC1Down;
   Float_t         MET_T1_phi_jesRelativePtEC1Down;
   Float_t         MET_T1Smear_pt_jesRelativePtEC1Down;
   Float_t         MET_T1Smear_phi_jesRelativePtEC1Down;
   Float_t         Jet_pt_jesRelativePtEC2Down[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtEC2Down[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativePtEC2Down;
   Float_t         MET_T1_phi_jesRelativePtEC2Down;
   Float_t         MET_T1Smear_pt_jesRelativePtEC2Down;
   Float_t         MET_T1Smear_phi_jesRelativePtEC2Down;
   Float_t         Jet_pt_jesRelativePtHFDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativePtHFDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativePtHFDown;
   Float_t         MET_T1_phi_jesRelativePtHFDown;
   Float_t         MET_T1Smear_pt_jesRelativePtHFDown;
   Float_t         MET_T1Smear_phi_jesRelativePtHFDown;
   Float_t         Jet_pt_jesRelativeBalDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeBalDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeBalDown;
   Float_t         MET_T1_phi_jesRelativeBalDown;
   Float_t         MET_T1Smear_pt_jesRelativeBalDown;
   Float_t         MET_T1Smear_phi_jesRelativeBalDown;
   Float_t         Jet_pt_jesRelativeSampleDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeSampleDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeSampleDown;
   Float_t         MET_T1_phi_jesRelativeSampleDown;
   Float_t         MET_T1Smear_pt_jesRelativeSampleDown;
   Float_t         MET_T1Smear_phi_jesRelativeSampleDown;
   Float_t         Jet_pt_jesRelativeFSRDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeFSRDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeFSRDown;
   Float_t         MET_T1_phi_jesRelativeFSRDown;
   Float_t         MET_T1Smear_pt_jesRelativeFSRDown;
   Float_t         MET_T1Smear_phi_jesRelativeFSRDown;
   Float_t         Jet_pt_jesRelativeStatFSRDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatFSRDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeStatFSRDown;
   Float_t         MET_T1_phi_jesRelativeStatFSRDown;
   Float_t         MET_T1Smear_pt_jesRelativeStatFSRDown;
   Float_t         MET_T1Smear_phi_jesRelativeStatFSRDown;
   Float_t         Jet_pt_jesRelativeStatECDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatECDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeStatECDown;
   Float_t         MET_T1_phi_jesRelativeStatECDown;
   Float_t         MET_T1Smear_pt_jesRelativeStatECDown;
   Float_t         MET_T1Smear_phi_jesRelativeStatECDown;
   Float_t         Jet_pt_jesRelativeStatHFDown[30];   //[nJet]
   Float_t         Jet_mass_jesRelativeStatHFDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesRelativeStatHFDown;
   Float_t         MET_T1_phi_jesRelativeStatHFDown;
   Float_t         MET_T1Smear_pt_jesRelativeStatHFDown;
   Float_t         MET_T1Smear_phi_jesRelativeStatHFDown;
   Float_t         Jet_pt_jesPileUpDataMCDown[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpDataMCDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpDataMCDown;
   Float_t         MET_T1_phi_jesPileUpDataMCDown;
   Float_t         MET_T1Smear_pt_jesPileUpDataMCDown;
   Float_t         MET_T1Smear_phi_jesPileUpDataMCDown;
   Float_t         Jet_pt_jesPileUpPtRefDown[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtRefDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpPtRefDown;
   Float_t         MET_T1_phi_jesPileUpPtRefDown;
   Float_t         MET_T1Smear_pt_jesPileUpPtRefDown;
   Float_t         MET_T1Smear_phi_jesPileUpPtRefDown;
   Float_t         Jet_pt_jesPileUpPtBBDown[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtBBDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpPtBBDown;
   Float_t         MET_T1_phi_jesPileUpPtBBDown;
   Float_t         MET_T1Smear_pt_jesPileUpPtBBDown;
   Float_t         MET_T1Smear_phi_jesPileUpPtBBDown;
   Float_t         Jet_pt_jesPileUpPtEC1Down[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtEC1Down[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpPtEC1Down;
   Float_t         MET_T1_phi_jesPileUpPtEC1Down;
   Float_t         MET_T1Smear_pt_jesPileUpPtEC1Down;
   Float_t         MET_T1Smear_phi_jesPileUpPtEC1Down;
   Float_t         Jet_pt_jesPileUpPtEC2Down[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtEC2Down[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpPtEC2Down;
   Float_t         MET_T1_phi_jesPileUpPtEC2Down;
   Float_t         MET_T1Smear_pt_jesPileUpPtEC2Down;
   Float_t         MET_T1Smear_phi_jesPileUpPtEC2Down;
   Float_t         Jet_pt_jesPileUpPtHFDown[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpPtHFDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpPtHFDown;
   Float_t         MET_T1_phi_jesPileUpPtHFDown;
   Float_t         MET_T1Smear_pt_jesPileUpPtHFDown;
   Float_t         MET_T1Smear_phi_jesPileUpPtHFDown;
   Float_t         Jet_pt_jesPileUpMuZeroDown[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpMuZeroDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpMuZeroDown;
   Float_t         MET_T1_phi_jesPileUpMuZeroDown;
   Float_t         MET_T1Smear_pt_jesPileUpMuZeroDown;
   Float_t         MET_T1Smear_phi_jesPileUpMuZeroDown;
   Float_t         Jet_pt_jesPileUpEnvelopeDown[30];   //[nJet]
   Float_t         Jet_mass_jesPileUpEnvelopeDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesPileUpEnvelopeDown;
   Float_t         MET_T1_phi_jesPileUpEnvelopeDown;
   Float_t         MET_T1Smear_pt_jesPileUpEnvelopeDown;
   Float_t         MET_T1Smear_phi_jesPileUpEnvelopeDown;
   Float_t         Jet_pt_jesSubTotalPileUpDown[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalPileUpDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesSubTotalPileUpDown;
   Float_t         MET_T1_phi_jesSubTotalPileUpDown;
   Float_t         MET_T1Smear_pt_jesSubTotalPileUpDown;
   Float_t         MET_T1Smear_phi_jesSubTotalPileUpDown;
   Float_t         Jet_pt_jesSubTotalRelativeDown[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalRelativeDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesSubTotalRelativeDown;
   Float_t         MET_T1_phi_jesSubTotalRelativeDown;
   Float_t         MET_T1Smear_pt_jesSubTotalRelativeDown;
   Float_t         MET_T1Smear_phi_jesSubTotalRelativeDown;
   Float_t         Jet_pt_jesSubTotalPtDown[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalPtDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesSubTotalPtDown;
   Float_t         MET_T1_phi_jesSubTotalPtDown;
   Float_t         MET_T1Smear_pt_jesSubTotalPtDown;
   Float_t         MET_T1Smear_phi_jesSubTotalPtDown;
   Float_t         Jet_pt_jesSubTotalScaleDown[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalScaleDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesSubTotalScaleDown;
   Float_t         MET_T1_phi_jesSubTotalScaleDown;
   Float_t         MET_T1Smear_pt_jesSubTotalScaleDown;
   Float_t         MET_T1Smear_phi_jesSubTotalScaleDown;
   Float_t         Jet_pt_jesSubTotalAbsoluteDown[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalAbsoluteDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesSubTotalAbsoluteDown;
   Float_t         MET_T1_phi_jesSubTotalAbsoluteDown;
   Float_t         MET_T1Smear_pt_jesSubTotalAbsoluteDown;
   Float_t         MET_T1Smear_phi_jesSubTotalAbsoluteDown;
   Float_t         Jet_pt_jesSubTotalMCDown[30];   //[nJet]
   Float_t         Jet_mass_jesSubTotalMCDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesSubTotalMCDown;
   Float_t         MET_T1_phi_jesSubTotalMCDown;
   Float_t         MET_T1Smear_pt_jesSubTotalMCDown;
   Float_t         MET_T1Smear_phi_jesSubTotalMCDown;
   Float_t         Jet_pt_jesTotalDown[30];   //[nJet]
   Float_t         Jet_mass_jesTotalDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesTotalDown;
   Float_t         MET_T1_phi_jesTotalDown;
   Float_t         MET_T1Smear_pt_jesTotalDown;
   Float_t         MET_T1Smear_phi_jesTotalDown;
   Float_t         Jet_pt_jesTotalNoFlavorDown[30];   //[nJet]
   Float_t         Jet_mass_jesTotalNoFlavorDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesTotalNoFlavorDown;
   Float_t         MET_T1_phi_jesTotalNoFlavorDown;
   Float_t         MET_T1Smear_pt_jesTotalNoFlavorDown;
   Float_t         MET_T1Smear_phi_jesTotalNoFlavorDown;
   Float_t         Jet_pt_jesTotalNoTimeDown[30];   //[nJet]
   Float_t         Jet_mass_jesTotalNoTimeDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesTotalNoTimeDown;
   Float_t         MET_T1_phi_jesTotalNoTimeDown;
   Float_t         MET_T1Smear_pt_jesTotalNoTimeDown;
   Float_t         MET_T1Smear_phi_jesTotalNoTimeDown;
   Float_t         Jet_pt_jesTotalNoFlavorNoTimeDown[30];   //[nJet]
   Float_t         Jet_mass_jesTotalNoFlavorNoTimeDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesTotalNoFlavorNoTimeDown;
   Float_t         MET_T1_phi_jesTotalNoFlavorNoTimeDown;
   Float_t         MET_T1Smear_pt_jesTotalNoFlavorNoTimeDown;
   Float_t         MET_T1Smear_phi_jesTotalNoFlavorNoTimeDown;
   Float_t         Jet_pt_jesFlavorZJetDown[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorZJetDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesFlavorZJetDown;
   Float_t         MET_T1_phi_jesFlavorZJetDown;
   Float_t         MET_T1Smear_pt_jesFlavorZJetDown;
   Float_t         MET_T1Smear_phi_jesFlavorZJetDown;
   Float_t         Jet_pt_jesFlavorPhotonJetDown[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPhotonJetDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesFlavorPhotonJetDown;
   Float_t         MET_T1_phi_jesFlavorPhotonJetDown;
   Float_t         MET_T1Smear_pt_jesFlavorPhotonJetDown;
   Float_t         MET_T1Smear_phi_jesFlavorPhotonJetDown;
   Float_t         Jet_pt_jesFlavorPureGluonDown[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureGluonDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesFlavorPureGluonDown;
   Float_t         MET_T1_phi_jesFlavorPureGluonDown;
   Float_t         MET_T1Smear_pt_jesFlavorPureGluonDown;
   Float_t         MET_T1Smear_phi_jesFlavorPureGluonDown;
   Float_t         Jet_pt_jesFlavorPureQuarkDown[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureQuarkDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesFlavorPureQuarkDown;
   Float_t         MET_T1_phi_jesFlavorPureQuarkDown;
   Float_t         MET_T1Smear_pt_jesFlavorPureQuarkDown;
   Float_t         MET_T1Smear_phi_jesFlavorPureQuarkDown;
   Float_t         Jet_pt_jesFlavorPureCharmDown[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureCharmDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesFlavorPureCharmDown;
   Float_t         MET_T1_phi_jesFlavorPureCharmDown;
   Float_t         MET_T1Smear_pt_jesFlavorPureCharmDown;
   Float_t         MET_T1Smear_phi_jesFlavorPureCharmDown;
   Float_t         Jet_pt_jesFlavorPureBottomDown[30];   //[nJet]
   Float_t         Jet_mass_jesFlavorPureBottomDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesFlavorPureBottomDown;
   Float_t         MET_T1_phi_jesFlavorPureBottomDown;
   Float_t         MET_T1Smear_pt_jesFlavorPureBottomDown;
   Float_t         MET_T1Smear_phi_jesFlavorPureBottomDown;
   Float_t         Jet_pt_jesTimeRunBDown[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunBDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesTimeRunBDown;
   Float_t         MET_T1_phi_jesTimeRunBDown;
   Float_t         MET_T1Smear_pt_jesTimeRunBDown;
   Float_t         MET_T1Smear_phi_jesTimeRunBDown;
   Float_t         Jet_pt_jesTimeRunCDown[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunCDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesTimeRunCDown;
   Float_t         MET_T1_phi_jesTimeRunCDown;
   Float_t         MET_T1Smear_pt_jesTimeRunCDown;
   Float_t         MET_T1Smear_phi_jesTimeRunCDown;
   Float_t         Jet_pt_jesTimeRunDEDown[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunDEDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesTimeRunDEDown;
   Float_t         MET_T1_phi_jesTimeRunDEDown;
   Float_t         MET_T1Smear_pt_jesTimeRunDEDown;
   Float_t         MET_T1Smear_phi_jesTimeRunDEDown;
   Float_t         Jet_pt_jesTimeRunFDown[30];   //[nJet]
   Float_t         Jet_mass_jesTimeRunFDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesTimeRunFDown;
   Float_t         MET_T1_phi_jesTimeRunFDown;
   Float_t         MET_T1Smear_pt_jesTimeRunFDown;
   Float_t         MET_T1Smear_phi_jesTimeRunFDown;
   Float_t         Jet_pt_jesCorrelationGroupMPFInSituDown[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupMPFInSituDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesCorrelationGroupMPFInSituDown;
   Float_t         MET_T1_phi_jesCorrelationGroupMPFInSituDown;
   Float_t         MET_T1Smear_pt_jesCorrelationGroupMPFInSituDown;
   Float_t         MET_T1Smear_phi_jesCorrelationGroupMPFInSituDown;
   Float_t         Jet_pt_jesCorrelationGroupIntercalibrationDown[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupIntercalibrationDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesCorrelationGroupIntercalibrationDown;
   Float_t         MET_T1_phi_jesCorrelationGroupIntercalibrationDown;
   Float_t         MET_T1Smear_pt_jesCorrelationGroupIntercalibrationDown;
   Float_t         MET_T1Smear_phi_jesCorrelationGroupIntercalibrationDown;
   Float_t         Jet_pt_jesCorrelationGroupbJESDown[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupbJESDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesCorrelationGroupbJESDown;
   Float_t         MET_T1_phi_jesCorrelationGroupbJESDown;
   Float_t         MET_T1Smear_pt_jesCorrelationGroupbJESDown;
   Float_t         MET_T1Smear_phi_jesCorrelationGroupbJESDown;
   Float_t         Jet_pt_jesCorrelationGroupFlavorDown[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupFlavorDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesCorrelationGroupFlavorDown;
   Float_t         MET_T1_phi_jesCorrelationGroupFlavorDown;
   Float_t         MET_T1Smear_pt_jesCorrelationGroupFlavorDown;
   Float_t         MET_T1Smear_phi_jesCorrelationGroupFlavorDown;
   Float_t         Jet_pt_jesCorrelationGroupUncorrelatedDown[30];   //[nJet]
   Float_t         Jet_mass_jesCorrelationGroupUncorrelatedDown[30];   //[nJet]
   Float_t         MET_T1_pt_jesCorrelationGroupUncorrelatedDown;
   Float_t         MET_T1_phi_jesCorrelationGroupUncorrelatedDown;
   Float_t         MET_T1Smear_pt_jesCorrelationGroupUncorrelatedDown;
   Float_t         MET_T1Smear_phi_jesCorrelationGroupUncorrelatedDown;
   Float_t         MET_T1_pt_unclustEnDown;
   Float_t         MET_T1_phi_unclustEnDown;
   Float_t         MET_T1Smear_pt_unclustEnDown;
   Float_t         MET_T1Smear_phi_unclustEnDown;
   Float_t         puWeight;
   Float_t         puWeightUp;
   Float_t         puWeightDown;
   Float_t         PrefireWeight;
   Float_t         PrefireWeight_Up;
   Float_t         PrefireWeight_Down;
   Float_t         Jet_btagSF_deepcsv_M_down[30];   //[nJet]
   Float_t         Jet_btagSF_deepcsv_M[30];   //[nJet]
   Float_t         Jet_btagSF_deepcsv_M_up[30];   //[nJet]
   Float_t         Jet_btagSF_deepcsv_L_down[30];   //[nJet]
   Float_t         Jet_btagSF_deepcsv_L[30];   //[nJet]
   Float_t         Jet_btagSF_deepcsv_L_up[30];   //[nJet]
   Float_t         Jet_btagSF_deepcsv_T_down[30];   //[nJet]
   Float_t         Jet_btagSF_deepcsv_T[30];   //[nJet]
   Float_t         Jet_btagSF_deepcsv_T_up[30];   //[nJet]
   Float_t         Jet_btagSF_deepjet_M_down[30];   //[nJet]
   Float_t         Jet_btagSF_deepjet_M[30];   //[nJet]
   Float_t         Jet_btagSF_deepjet_M_up[30];   //[nJet]
   Float_t         Jet_btagSF_deepjet_L_down[30];   //[nJet]
   Float_t         Jet_btagSF_deepjet_L[30];   //[nJet]
   Float_t         Jet_btagSF_deepjet_L_up[30];   //[nJet]
   Float_t         Jet_btagSF_deepjet_T_down[30];   //[nJet]
   Float_t         Jet_btagSF_deepjet_T[30];   //[nJet]
   Float_t         Jet_btagSF_deepjet_T_up[30];   //[nJet]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_HTXS_Higgs_pt;   //!
   TBranch        *b_HTXS_Higgs_y;   //!
   TBranch        *b_HTXS_stage1_1_cat_pTjet25GeV;   //!
   TBranch        *b_HTXS_stage1_1_cat_pTjet30GeV;   //!
   TBranch        *b_HTXS_stage1_1_fine_cat_pTjet25GeV;   //!
   TBranch        *b_HTXS_stage1_1_fine_cat_pTjet30GeV;   //!
   TBranch        *b_HTXS_stage1_2_cat_pTjet25GeV;   //!
   TBranch        *b_HTXS_stage1_2_cat_pTjet30GeV;   //!
   TBranch        *b_HTXS_stage1_2_fine_cat_pTjet25GeV;   //!
   TBranch        *b_HTXS_stage1_2_fine_cat_pTjet30GeV;   //!
   TBranch        *b_HTXS_stage_0;   //!
   TBranch        *b_HTXS_stage_1_pTjet25;   //!
   TBranch        *b_HTXS_stage_1_pTjet30;   //!
   TBranch        *b_HTXS_njets25;   //!
   TBranch        *b_HTXS_njets30;   //!
   TBranch        *b_nboostedTau;   //!
   TBranch        *b_boostedTau_chargedIso;   //!
   TBranch        *b_boostedTau_eta;   //!
   TBranch        *b_boostedTau_leadTkDeltaEta;   //!
   TBranch        *b_boostedTau_leadTkDeltaPhi;   //!
   TBranch        *b_boostedTau_leadTkPtOverTauPt;   //!
   TBranch        *b_boostedTau_mass;   //!
   TBranch        *b_boostedTau_neutralIso;   //!
   TBranch        *b_boostedTau_phi;   //!
   TBranch        *b_boostedTau_photonsOutsideSignalCone;   //!
   TBranch        *b_boostedTau_pt;   //!
   TBranch        *b_boostedTau_puCorr;   //!
   TBranch        *b_boostedTau_rawAntiEle2018;   //!
   TBranch        *b_boostedTau_rawIso;   //!
   TBranch        *b_boostedTau_rawIsodR03;   //!
   TBranch        *b_boostedTau_rawMVAnewDM2017v2;   //!
   TBranch        *b_boostedTau_rawMVAoldDM2017v2;   //!
   TBranch        *b_boostedTau_rawMVAoldDMdR032017v2;   //!
   TBranch        *b_boostedTau_charge;   //!
   TBranch        *b_boostedTau_decayMode;   //!
   TBranch        *b_boostedTau_jetIdx;   //!
   TBranch        *b_boostedTau_rawAntiEleCat2018;   //!
   TBranch        *b_boostedTau_idAntiEle2018;   //!
   TBranch        *b_boostedTau_idAntiMu;   //!
   TBranch        *b_boostedTau_idMVAnewDM2017v2;   //!
   TBranch        *b_boostedTau_idMVAoldDM2017v2;   //!
   TBranch        *b_boostedTau_idMVAoldDMdR032017v2;   //!
   TBranch        *b_btagWeight_CSVV2;   //!
   TBranch        *b_btagWeight_DeepCSVB;   //!
   TBranch        *b_CaloMET_phi;   //!
   TBranch        *b_CaloMET_pt;   //!
   TBranch        *b_CaloMET_sumEt;   //!
   TBranch        *b_ChsMET_phi;   //!
   TBranch        *b_ChsMET_pt;   //!
   TBranch        *b_ChsMET_sumEt;   //!
   TBranch        *b_nCorrT1METJet;   //!
   TBranch        *b_CorrT1METJet_area;   //!
   TBranch        *b_CorrT1METJet_eta;   //!
   TBranch        *b_CorrT1METJet_muonSubtrFactor;   //!
   TBranch        *b_CorrT1METJet_phi;   //!
   TBranch        *b_CorrT1METJet_rawPt;   //!
   TBranch        *b_DeepMETResolutionTune_phi;   //!
   TBranch        *b_DeepMETResolutionTune_pt;   //!
   TBranch        *b_DeepMETResponseTune_phi;   //!
   TBranch        *b_DeepMETResponseTune_pt;   //!
   TBranch        *b_nElectron;   //!
   TBranch        *b_Electron_dEscaleDown;   //!
   TBranch        *b_Electron_dEscaleUp;   //!
   TBranch        *b_Electron_dEsigmaDown;   //!
   TBranch        *b_Electron_dEsigmaUp;   //!
   TBranch        *b_Electron_deltaEtaSC;   //!
   TBranch        *b_Electron_dr03EcalRecHitSumEt;   //!
   TBranch        *b_Electron_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_Electron_dr03TkSumPt;   //!
   TBranch        *b_Electron_dr03TkSumPtHEEP;   //!
   TBranch        *b_Electron_dxy;   //!
   TBranch        *b_Electron_dxyErr;   //!
   TBranch        *b_Electron_dz;   //!
   TBranch        *b_Electron_dzErr;   //!
   TBranch        *b_Electron_eCorr;   //!
   TBranch        *b_Electron_eInvMinusPInv;   //!
   TBranch        *b_Electron_energyErr;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_hoe;   //!
   TBranch        *b_Electron_ip3d;   //!
   TBranch        *b_Electron_jetPtRelv2;   //!
   TBranch        *b_Electron_jetRelIso;   //!
   TBranch        *b_Electron_mass;   //!
   TBranch        *b_Electron_miniPFRelIso_all;   //!
   TBranch        *b_Electron_miniPFRelIso_chg;   //!
   TBranch        *b_Electron_mvaFall17V2Iso;   //!
   TBranch        *b_Electron_mvaFall17V2noIso;   //!
   TBranch        *b_Electron_pfRelIso03_all;   //!
   TBranch        *b_Electron_pfRelIso03_chg;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_r9;   //!
   TBranch        *b_Electron_scEtOverPt;   //!
   TBranch        *b_Electron_sieie;   //!
   TBranch        *b_Electron_sip3d;   //!
   TBranch        *b_Electron_mvaTTH;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_cutBased;   //!
   TBranch        *b_Electron_jetIdx;   //!
   TBranch        *b_Electron_pdgId;   //!
   TBranch        *b_Electron_photonIdx;   //!
   TBranch        *b_Electron_tightCharge;   //!
   TBranch        *b_Electron_vidNestedWPBitmap;   //!
   TBranch        *b_Electron_vidNestedWPBitmapHEEP;   //!
   TBranch        *b_Electron_convVeto;   //!
   TBranch        *b_Electron_cutBased_HEEP;   //!
   TBranch        *b_Electron_isPFcand;   //!
   TBranch        *b_Electron_jetNDauCharged;   //!
   TBranch        *b_Electron_lostHits;   //!
   TBranch        *b_Electron_mvaFall17V2Iso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V2Iso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V2Iso_WPL;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WPL;   //!
   TBranch        *b_Electron_seedGain;   //!
   TBranch        *b_nFatJet;   //!
   TBranch        *b_FatJet_area;   //!
   TBranch        *b_FatJet_btagCSVV2;   //!
   TBranch        *b_FatJet_btagDDBvLV2;   //!
   TBranch        *b_FatJet_btagDDCvBV2;   //!
   TBranch        *b_FatJet_btagDDCvLV2;   //!
   TBranch        *b_FatJet_btagDeepB;   //!
   TBranch        *b_FatJet_btagHbb;   //!
   TBranch        *b_FatJet_deepTagMD_H4qvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_HbbvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_TvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_WvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZHbbvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZHccvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZbbvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_bbvsLight;   //!
   TBranch        *b_FatJet_deepTagMD_ccvsLight;   //!
   TBranch        *b_FatJet_deepTag_H;   //!
   TBranch        *b_FatJet_deepTag_QCD;   //!
   TBranch        *b_FatJet_deepTag_QCDothers;   //!
   TBranch        *b_FatJet_deepTag_TvsQCD;   //!
   TBranch        *b_FatJet_deepTag_WvsQCD;   //!
   TBranch        *b_FatJet_deepTag_ZvsQCD;   //!
   TBranch        *b_FatJet_eta;   //!
   TBranch        *b_FatJet_mass;   //!
   TBranch        *b_FatJet_msoftdrop;   //!
   TBranch        *b_FatJet_n2b1;   //!
   TBranch        *b_FatJet_n3b1;   //!
   TBranch        *b_FatJet_particleNetMD_QCD;   //!
   TBranch        *b_FatJet_particleNetMD_Xbb;   //!
   TBranch        *b_FatJet_particleNetMD_Xcc;   //!
   TBranch        *b_FatJet_particleNetMD_Xqq;   //!
   TBranch        *b_FatJet_particleNet_H4qvsQCD;   //!
   TBranch        *b_FatJet_particleNet_HbbvsQCD;   //!
   TBranch        *b_FatJet_particleNet_HccvsQCD;   //!
   TBranch        *b_FatJet_particleNet_QCD;   //!
   TBranch        *b_FatJet_particleNet_TvsQCD;   //!
   TBranch        *b_FatJet_particleNet_WvsQCD;   //!
   TBranch        *b_FatJet_particleNet_ZvsQCD;   //!
   TBranch        *b_FatJet_particleNet_mass;   //!
   TBranch        *b_FatJet_phi;   //!
   TBranch        *b_FatJet_pt;   //!
   TBranch        *b_FatJet_rawFactor;   //!
   TBranch        *b_FatJet_tau1;   //!
   TBranch        *b_FatJet_tau2;   //!
   TBranch        *b_FatJet_tau3;   //!
   TBranch        *b_FatJet_tau4;   //!
   TBranch        *b_FatJet_lsf3;   //!
   TBranch        *b_FatJet_jetId;   //!
   TBranch        *b_FatJet_subJetIdx1;   //!
   TBranch        *b_FatJet_subJetIdx2;   //!
   TBranch        *b_FatJet_electronIdx3SJ;   //!
   TBranch        *b_FatJet_muonIdx3SJ;   //!
   TBranch        *b_FatJet_nConstituents;   //!
   TBranch        *b_nFsrPhoton;   //!
   TBranch        *b_FsrPhoton_dROverEt2;   //!
   TBranch        *b_FsrPhoton_eta;   //!
   TBranch        *b_FsrPhoton_phi;   //!
   TBranch        *b_FsrPhoton_pt;   //!
   TBranch        *b_FsrPhoton_relIso03;   //!
   TBranch        *b_FsrPhoton_muonIdx;   //!
   TBranch        *b_nGenJetAK8;   //!
   TBranch        *b_GenJetAK8_eta;   //!
   TBranch        *b_GenJetAK8_mass;   //!
   TBranch        *b_GenJetAK8_phi;   //!
   TBranch        *b_GenJetAK8_pt;   //!
   TBranch        *b_nGenJet;   //!
   TBranch        *b_GenJet_eta;   //!
   TBranch        *b_GenJet_mass;   //!
   TBranch        *b_GenJet_phi;   //!
   TBranch        *b_GenJet_pt;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_mass;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_genPartIdxMother;   //!
   TBranch        *b_GenPart_pdgId;   //!
   TBranch        *b_GenPart_status;   //!
   TBranch        *b_GenPart_statusFlags;   //!
   TBranch        *b_nSubGenJetAK8;   //!
   TBranch        *b_SubGenJetAK8_eta;   //!
   TBranch        *b_SubGenJetAK8_mass;   //!
   TBranch        *b_SubGenJetAK8_phi;   //!
   TBranch        *b_SubGenJetAK8_pt;   //!
   TBranch        *b_Generator_binvar;   //!
   TBranch        *b_Generator_scalePDF;   //!
   TBranch        *b_Generator_weight;   //!
   TBranch        *b_Generator_x1;   //!
   TBranch        *b_Generator_x2;   //!
   TBranch        *b_Generator_xpdf1;   //!
   TBranch        *b_Generator_xpdf2;   //!
   TBranch        *b_Generator_id1;   //!
   TBranch        *b_Generator_id2;   //!
   TBranch        *b_GenVtx_x;   //!
   TBranch        *b_GenVtx_y;   //!
   TBranch        *b_GenVtx_z;   //!
   TBranch        *b_nGenVisTau;   //!
   TBranch        *b_GenVisTau_eta;   //!
   TBranch        *b_GenVisTau_mass;   //!
   TBranch        *b_GenVisTau_phi;   //!
   TBranch        *b_GenVisTau_pt;   //!
   TBranch        *b_GenVisTau_charge;   //!
   TBranch        *b_GenVisTau_genPartIdxMother;   //!
   TBranch        *b_GenVisTau_status;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_nEFTfitCoefficients;   //!
   TBranch        *b_EFTfitCoefficients;   //!
   TBranch        *b_LHEWeight_originalXWGTUP;   //!
   TBranch        *b_nLHEPdfWeight;   //!
   TBranch        *b_LHEPdfWeight;   //!
   TBranch        *b_nLHEReweightingWeight;   //!
   TBranch        *b_LHEReweightingWeight;   //!
   TBranch        *b_nLHEScaleWeight;   //!
   TBranch        *b_LHEScaleWeight;   //!
   TBranch        *b_nPSWeight;   //!
   TBranch        *b_PSWeight;   //!
   TBranch        *b_nWCnames;   //!
   TBranch        *b_WCnames;   //!
   TBranch        *b_nIsoTrack;   //!
   TBranch        *b_IsoTrack_dxy;   //!
   TBranch        *b_IsoTrack_dz;   //!
   TBranch        *b_IsoTrack_eta;   //!
   TBranch        *b_IsoTrack_pfRelIso03_all;   //!
   TBranch        *b_IsoTrack_pfRelIso03_chg;   //!
   TBranch        *b_IsoTrack_phi;   //!
   TBranch        *b_IsoTrack_pt;   //!
   TBranch        *b_IsoTrack_miniPFRelIso_all;   //!
   TBranch        *b_IsoTrack_miniPFRelIso_chg;   //!
   TBranch        *b_IsoTrack_charge;   //!
   TBranch        *b_IsoTrack_fromPV;   //!
   TBranch        *b_IsoTrack_pdgId;   //!
   TBranch        *b_IsoTrack_isHighPurityTrack;   //!
   TBranch        *b_IsoTrack_isPFcand;   //!
   TBranch        *b_IsoTrack_isFromLostTrack;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_btagCSVV2;   //!
   TBranch        *b_Jet_btagDeepB;   //!
   TBranch        *b_Jet_btagDeepCvB;   //!
   TBranch        *b_Jet_btagDeepCvL;   //!
   TBranch        *b_Jet_btagDeepFlavB;   //!
   TBranch        *b_Jet_btagDeepFlavCvB;   //!
   TBranch        *b_Jet_btagDeepFlavCvL;   //!
   TBranch        *b_Jet_btagDeepFlavQG;   //!
   TBranch        *b_Jet_chEmEF;   //!
   TBranch        *b_Jet_chFPV0EF;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_hfsigmaEtaEta;   //!
   TBranch        *b_Jet_hfsigmaPhiPhi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_muEF;   //!
   TBranch        *b_Jet_muonSubtrFactor;   //!
   TBranch        *b_Jet_neEmEF;   //!
   TBranch        *b_Jet_neHEF;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_puIdDisc;   //!
   TBranch        *b_Jet_qgl;   //!
   TBranch        *b_Jet_rawFactor;   //!
   TBranch        *b_Jet_bRegCorr;   //!
   TBranch        *b_Jet_bRegRes;   //!
   TBranch        *b_Jet_cRegCorr;   //!
   TBranch        *b_Jet_cRegRes;   //!
   TBranch        *b_Jet_electronIdx1;   //!
   TBranch        *b_Jet_electronIdx2;   //!
   TBranch        *b_Jet_hfadjacentEtaStripsSize;   //!
   TBranch        *b_Jet_hfcentralEtaStripSize;   //!
   TBranch        *b_Jet_jetId;   //!
   TBranch        *b_Jet_muonIdx1;   //!
   TBranch        *b_Jet_muonIdx2;   //!
   TBranch        *b_Jet_nElectrons;   //!
   TBranch        *b_Jet_nMuons;   //!
   TBranch        *b_Jet_puId;   //!
   TBranch        *b_Jet_nConstituents;   //!
   TBranch        *b_L1PreFiringWeight_Dn;   //!
   TBranch        *b_L1PreFiringWeight_ECAL_Dn;   //!
   TBranch        *b_L1PreFiringWeight_ECAL_Nom;   //!
   TBranch        *b_L1PreFiringWeight_ECAL_Up;   //!
   TBranch        *b_L1PreFiringWeight_Muon_Nom;   //!
   TBranch        *b_L1PreFiringWeight_Muon_StatDn;   //!
   TBranch        *b_L1PreFiringWeight_Muon_StatUp;   //!
   TBranch        *b_L1PreFiringWeight_Muon_SystDn;   //!
   TBranch        *b_L1PreFiringWeight_Muon_SystUp;   //!
   TBranch        *b_L1PreFiringWeight_Nom;   //!
   TBranch        *b_L1PreFiringWeight_Up;   //!
   TBranch        *b_LHE_HT;   //!
   TBranch        *b_LHE_HTIncoming;   //!
   TBranch        *b_LHE_Vpt;   //!
   TBranch        *b_LHE_AlphaS;   //!
   TBranch        *b_LHE_Njets;   //!
   TBranch        *b_LHE_Nb;   //!
   TBranch        *b_LHE_Nc;   //!
   TBranch        *b_LHE_Nuds;   //!
   TBranch        *b_LHE_Nglu;   //!
   TBranch        *b_LHE_NpNLO;   //!
   TBranch        *b_LHE_NpLO;   //!
   TBranch        *b_nLHEPart;   //!
   TBranch        *b_LHEPart_pt;   //!
   TBranch        *b_LHEPart_eta;   //!
   TBranch        *b_LHEPart_phi;   //!
   TBranch        *b_LHEPart_mass;   //!
   TBranch        *b_LHEPart_incomingpz;   //!
   TBranch        *b_LHEPart_pdgId;   //!
   TBranch        *b_LHEPart_status;   //!
   TBranch        *b_LHEPart_spin;   //!
   TBranch        *b_nLowPtElectron;   //!
   TBranch        *b_LowPtElectron_ID;   //!
   TBranch        *b_LowPtElectron_convVtxRadius;   //!
   TBranch        *b_LowPtElectron_deltaEtaSC;   //!
   TBranch        *b_LowPtElectron_dxy;   //!
   TBranch        *b_LowPtElectron_dxyErr;   //!
   TBranch        *b_LowPtElectron_dz;   //!
   TBranch        *b_LowPtElectron_dzErr;   //!
   TBranch        *b_LowPtElectron_eInvMinusPInv;   //!
   TBranch        *b_LowPtElectron_embeddedID;   //!
   TBranch        *b_LowPtElectron_energyErr;   //!
   TBranch        *b_LowPtElectron_eta;   //!
   TBranch        *b_LowPtElectron_hoe;   //!
   TBranch        *b_LowPtElectron_mass;   //!
   TBranch        *b_LowPtElectron_miniPFRelIso_all;   //!
   TBranch        *b_LowPtElectron_miniPFRelIso_chg;   //!
   TBranch        *b_LowPtElectron_phi;   //!
   TBranch        *b_LowPtElectron_pt;   //!
   TBranch        *b_LowPtElectron_ptbiased;   //!
   TBranch        *b_LowPtElectron_r9;   //!
   TBranch        *b_LowPtElectron_scEtOverPt;   //!
   TBranch        *b_LowPtElectron_sieie;   //!
   TBranch        *b_LowPtElectron_unbiased;   //!
   TBranch        *b_LowPtElectron_charge;   //!
   TBranch        *b_LowPtElectron_convWP;   //!
   TBranch        *b_LowPtElectron_pdgId;   //!
   TBranch        *b_LowPtElectron_convVeto;   //!
   TBranch        *b_LowPtElectron_lostHits;   //!
   TBranch        *b_GenMET_phi;   //!
   TBranch        *b_GenMET_pt;   //!
   TBranch        *b_MET_MetUnclustEnUpDeltaX;   //!
   TBranch        *b_MET_MetUnclustEnUpDeltaY;   //!
   TBranch        *b_MET_covXX;   //!
   TBranch        *b_MET_covXY;   //!
   TBranch        *b_MET_covYY;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_pt;   //!
   TBranch        *b_MET_significance;   //!
   TBranch        *b_MET_sumEt;   //!
   TBranch        *b_MET_sumPtUnclustered;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_dxyErr;   //!
   TBranch        *b_Muon_dxybs;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_dzErr;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_ip3d;   //!
   TBranch        *b_Muon_jetPtRelv2;   //!
   TBranch        *b_Muon_jetRelIso;   //!
   TBranch        *b_Muon_mass;   //!
   TBranch        *b_Muon_miniPFRelIso_all;   //!
   TBranch        *b_Muon_miniPFRelIso_chg;   //!
   TBranch        *b_Muon_pfRelIso03_all;   //!
   TBranch        *b_Muon_pfRelIso03_chg;   //!
   TBranch        *b_Muon_pfRelIso04_all;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_ptErr;   //!
   TBranch        *b_Muon_segmentComp;   //!
   TBranch        *b_Muon_sip3d;   //!
   TBranch        *b_Muon_softMva;   //!
   TBranch        *b_Muon_tkRelIso;   //!
   TBranch        *b_Muon_tunepRelPt;   //!
   TBranch        *b_Muon_mvaLowPt;   //!
   TBranch        *b_Muon_mvaTTH;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_jetIdx;   //!
   TBranch        *b_Muon_nStations;   //!
   TBranch        *b_Muon_nTrackerLayers;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_tightCharge;   //!
   TBranch        *b_Muon_fsrPhotonIdx;   //!
   TBranch        *b_Muon_highPtId;   //!
   TBranch        *b_Muon_highPurity;   //!
   TBranch        *b_Muon_inTimeMuon;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_isPFcand;   //!
   TBranch        *b_Muon_isStandalone;   //!
   TBranch        *b_Muon_isTracker;   //!
   TBranch        *b_Muon_jetNDauCharged;   //!
   TBranch        *b_Muon_looseId;   //!
   TBranch        *b_Muon_mediumId;   //!
   TBranch        *b_Muon_mediumPromptId;   //!
   TBranch        *b_Muon_miniIsoId;   //!
   TBranch        *b_Muon_multiIsoId;   //!
   TBranch        *b_Muon_mvaId;   //!
   TBranch        *b_Muon_mvaLowPtId;   //!
   TBranch        *b_Muon_pfIsoId;   //!
   TBranch        *b_Muon_puppiIsoId;   //!
   TBranch        *b_Muon_softId;   //!
   TBranch        *b_Muon_softMvaId;   //!
   TBranch        *b_Muon_tightId;   //!
   TBranch        *b_Muon_tkIsoId;   //!
   TBranch        *b_Muon_triggerIdLoose;   //!
   TBranch        *b_nPhoton;   //!
   TBranch        *b_Photon_dEscaleDown;   //!
   TBranch        *b_Photon_dEscaleUp;   //!
   TBranch        *b_Photon_dEsigmaDown;   //!
   TBranch        *b_Photon_dEsigmaUp;   //!
   TBranch        *b_Photon_eCorr;   //!
   TBranch        *b_Photon_energyErr;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_hoe;   //!
   TBranch        *b_Photon_mass;   //!
   TBranch        *b_Photon_mvaID;   //!
   TBranch        *b_Photon_mvaID_Fall17V1p1;   //!
   TBranch        *b_Photon_pfRelIso03_all;   //!
   TBranch        *b_Photon_pfRelIso03_chg;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_r9;   //!
   TBranch        *b_Photon_sieie;   //!
   TBranch        *b_Photon_charge;   //!
   TBranch        *b_Photon_cutBased;   //!
   TBranch        *b_Photon_cutBased_Fall17V1Bitmap;   //!
   TBranch        *b_Photon_electronIdx;   //!
   TBranch        *b_Photon_jetIdx;   //!
   TBranch        *b_Photon_pdgId;   //!
   TBranch        *b_Photon_vidNestedWPBitmap;   //!
   TBranch        *b_Photon_electronVeto;   //!
   TBranch        *b_Photon_isScEtaEB;   //!
   TBranch        *b_Photon_isScEtaEE;   //!
   TBranch        *b_Photon_mvaID_WP80;   //!
   TBranch        *b_Photon_mvaID_WP90;   //!
   TBranch        *b_Photon_pixelSeed;   //!
   TBranch        *b_Photon_seedGain;   //!
   TBranch        *b_Pileup_nTrueInt;   //!
   TBranch        *b_Pileup_pudensity;   //!
   TBranch        *b_Pileup_gpudensity;   //!
   TBranch        *b_Pileup_nPU;   //!
   TBranch        *b_Pileup_sumEOOT;   //!
   TBranch        *b_Pileup_sumLOOT;   //!
   TBranch        *b_PuppiMET_phi;   //!
   TBranch        *b_PuppiMET_phiJERDown;   //!
   TBranch        *b_PuppiMET_phiJERUp;   //!
   TBranch        *b_PuppiMET_phiJESDown;   //!
   TBranch        *b_PuppiMET_phiJESUp;   //!
   TBranch        *b_PuppiMET_phiUnclusteredDown;   //!
   TBranch        *b_PuppiMET_phiUnclusteredUp;   //!
   TBranch        *b_PuppiMET_pt;   //!
   TBranch        *b_PuppiMET_ptJERDown;   //!
   TBranch        *b_PuppiMET_ptJERUp;   //!
   TBranch        *b_PuppiMET_ptJESDown;   //!
   TBranch        *b_PuppiMET_ptJESUp;   //!
   TBranch        *b_PuppiMET_ptUnclusteredDown;   //!
   TBranch        *b_PuppiMET_ptUnclusteredUp;   //!
   TBranch        *b_PuppiMET_sumEt;   //!
   TBranch        *b_RawMET_phi;   //!
   TBranch        *b_RawMET_pt;   //!
   TBranch        *b_RawMET_sumEt;   //!
   TBranch        *b_RawPuppiMET_phi;   //!
   TBranch        *b_RawPuppiMET_pt;   //!
   TBranch        *b_RawPuppiMET_sumEt;   //!
   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_fixedGridRhoFastjetCentral;   //!
   TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_nGenDressedLepton;   //!
   TBranch        *b_GenDressedLepton_eta;   //!
   TBranch        *b_GenDressedLepton_mass;   //!
   TBranch        *b_GenDressedLepton_phi;   //!
   TBranch        *b_GenDressedLepton_pt;   //!
   TBranch        *b_GenDressedLepton_pdgId;   //!
   TBranch        *b_GenDressedLepton_hasTauAnc;   //!
   TBranch        *b_nGenIsolatedPhoton;   //!
   TBranch        *b_GenIsolatedPhoton_eta;   //!
   TBranch        *b_GenIsolatedPhoton_mass;   //!
   TBranch        *b_GenIsolatedPhoton_phi;   //!
   TBranch        *b_GenIsolatedPhoton_pt;   //!
   TBranch        *b_nSoftActivityJet;   //!
   TBranch        *b_SoftActivityJet_eta;   //!
   TBranch        *b_SoftActivityJet_phi;   //!
   TBranch        *b_SoftActivityJet_pt;   //!
   TBranch        *b_SoftActivityJetHT;   //!
   TBranch        *b_SoftActivityJetHT10;   //!
   TBranch        *b_SoftActivityJetHT2;   //!
   TBranch        *b_SoftActivityJetHT5;   //!
   TBranch        *b_SoftActivityJetNjets10;   //!
   TBranch        *b_SoftActivityJetNjets2;   //!
   TBranch        *b_SoftActivityJetNjets5;   //!
   TBranch        *b_nSubJet;   //!
   TBranch        *b_SubJet_btagCSVV2;   //!
   TBranch        *b_SubJet_btagDeepB;   //!
   TBranch        *b_SubJet_eta;   //!
   TBranch        *b_SubJet_mass;   //!
   TBranch        *b_SubJet_n2b1;   //!
   TBranch        *b_SubJet_n3b1;   //!
   TBranch        *b_SubJet_phi;   //!
   TBranch        *b_SubJet_pt;   //!
   TBranch        *b_SubJet_rawFactor;   //!
   TBranch        *b_SubJet_tau1;   //!
   TBranch        *b_SubJet_tau2;   //!
   TBranch        *b_SubJet_tau3;   //!
   TBranch        *b_SubJet_tau4;   //!
   TBranch        *b_nTau;   //!
   TBranch        *b_Tau_chargedIso;   //!
   TBranch        *b_Tau_dxy;   //!
   TBranch        *b_Tau_dz;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_leadTkDeltaEta;   //!
   TBranch        *b_Tau_leadTkDeltaPhi;   //!
   TBranch        *b_Tau_leadTkPtOverTauPt;   //!
   TBranch        *b_Tau_mass;   //!
   TBranch        *b_Tau_neutralIso;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_photonsOutsideSignalCone;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_puCorr;   //!
   TBranch        *b_Tau_rawDeepTau2017v2p1VSe;   //!
   TBranch        *b_Tau_rawDeepTau2017v2p1VSjet;   //!
   TBranch        *b_Tau_rawDeepTau2017v2p1VSmu;   //!
   TBranch        *b_Tau_rawIso;   //!
   TBranch        *b_Tau_rawIsodR03;   //!
   TBranch        *b_Tau_charge;   //!
   TBranch        *b_Tau_decayMode;   //!
   TBranch        *b_Tau_jetIdx;   //!
   TBranch        *b_Tau_idAntiEleDeadECal;   //!
   TBranch        *b_Tau_idAntiMu;   //!
   TBranch        *b_Tau_idDecayModeOldDMs;   //!
   TBranch        *b_Tau_idDeepTau2017v2p1VSe;   //!
   TBranch        *b_Tau_idDeepTau2017v2p1VSjet;   //!
   TBranch        *b_Tau_idDeepTau2017v2p1VSmu;   //!
   TBranch        *b_TkMET_phi;   //!
   TBranch        *b_TkMET_pt;   //!
   TBranch        *b_TkMET_sumEt;   //!
   TBranch        *b_nTrigObj;   //!
   TBranch        *b_TrigObj_pt;   //!
   TBranch        *b_TrigObj_eta;   //!
   TBranch        *b_TrigObj_phi;   //!
   TBranch        *b_TrigObj_l1pt;   //!
   TBranch        *b_TrigObj_l1pt_2;   //!
   TBranch        *b_TrigObj_l2pt;   //!
   TBranch        *b_TrigObj_id;   //!
   TBranch        *b_TrigObj_l1iso;   //!
   TBranch        *b_TrigObj_l1charge;   //!
   TBranch        *b_TrigObj_filterBits;   //!
   TBranch        *b_genTtbarId;   //!
   TBranch        *b_nOtherPV;   //!
   TBranch        *b_OtherPV_z;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_nSV;   //!
   TBranch        *b_SV_dlen;   //!
   TBranch        *b_SV_dlenSig;   //!
   TBranch        *b_SV_dxy;   //!
   TBranch        *b_SV_dxySig;   //!
   TBranch        *b_SV_pAngle;   //!
   TBranch        *b_SV_charge;   //!
   TBranch        *b_boostedTau_genPartIdx;   //!
   TBranch        *b_boostedTau_genPartFlav;   //!
   TBranch        *b_Electron_genPartIdx;   //!
   TBranch        *b_Electron_genPartFlav;   //!
   TBranch        *b_FatJet_genJetAK8Idx;   //!
   TBranch        *b_FatJet_hadronFlavour;   //!
   TBranch        *b_FatJet_nBHadrons;   //!
   TBranch        *b_FatJet_nCHadrons;   //!
   TBranch        *b_GenJetAK8_partonFlavour;   //!
   TBranch        *b_GenJetAK8_hadronFlavour;   //!
   TBranch        *b_GenJet_partonFlavour;   //!
   TBranch        *b_GenJet_hadronFlavour;   //!
   TBranch        *b_GenVtx_t0;   //!
   TBranch        *b_Jet_genJetIdx;   //!
   TBranch        *b_Jet_hadronFlavour;   //!
   TBranch        *b_Jet_partonFlavour;   //!
   TBranch        *b_LowPtElectron_genPartIdx;   //!
   TBranch        *b_LowPtElectron_genPartFlav;   //!
   TBranch        *b_Muon_genPartIdx;   //!
   TBranch        *b_Muon_genPartFlav;   //!
   TBranch        *b_Photon_genPartIdx;   //!
   TBranch        *b_Photon_genPartFlav;   //!
   TBranch        *b_MET_fiducialGenPhi;   //!
   TBranch        *b_MET_fiducialGenPt;   //!
   TBranch        *b_Electron_cleanmask;   //!
   TBranch        *b_Jet_cleanmask;   //!
   TBranch        *b_Muon_cleanmask;   //!
   TBranch        *b_Photon_cleanmask;   //!
   TBranch        *b_Tau_cleanmask;   //!
   TBranch        *b_SubJet_hadronFlavour;   //!
   TBranch        *b_SubJet_nBHadrons;   //!
   TBranch        *b_SubJet_nCHadrons;   //!
   TBranch        *b_SV_chi2;   //!
   TBranch        *b_SV_eta;   //!
   TBranch        *b_SV_mass;   //!
   TBranch        *b_SV_ndof;   //!
   TBranch        *b_SV_phi;   //!
   TBranch        *b_SV_pt;   //!
   TBranch        *b_SV_x;   //!
   TBranch        *b_SV_y;   //!
   TBranch        *b_SV_z;   //!
   TBranch        *b_SV_ntracks;   //!
   TBranch        *b_Tau_genPartIdx;   //!
   TBranch        *b_Tau_genPartFlav;   //!
   TBranch        *b_L1_AlwaysTrue;   //!
   TBranch        *b_L1_BPTX_AND_Ref1_VME;   //!
   TBranch        *b_L1_BPTX_AND_Ref3_VME;   //!
   TBranch        *b_L1_BPTX_AND_Ref4_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_B1_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_B2_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_Ref1_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_Ref2_VME;   //!
   TBranch        *b_L1_BPTX_NotOR_VME;   //!
   TBranch        *b_L1_BPTX_OR_Ref3_VME;   //!
   TBranch        *b_L1_BPTX_OR_Ref4_VME;   //!
   TBranch        *b_L1_BPTX_RefAND_VME;   //!
   TBranch        *b_L1_BptxMinus;   //!
   TBranch        *b_L1_BptxOR;   //!
   TBranch        *b_L1_BptxPlus;   //!
   TBranch        *b_L1_BptxXOR;   //!
   TBranch        *b_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;   //!
   TBranch        *b_L1_DoubleEG6_HTT240er;   //!
   TBranch        *b_L1_DoubleEG6_HTT250er;   //!
   TBranch        *b_L1_DoubleEG6_HTT255er;   //!
   TBranch        *b_L1_DoubleEG6_HTT270er;   //!
   TBranch        *b_L1_DoubleEG6_HTT300er;   //!
   TBranch        *b_L1_DoubleEG8er2p6_HTT255er;   //!
   TBranch        *b_L1_DoubleEG8er2p6_HTT270er;   //!
   TBranch        *b_L1_DoubleEG8er2p6_HTT300er;   //!
   TBranch        *b_L1_DoubleEG_15_10;   //!
   TBranch        *b_L1_DoubleEG_18_17;   //!
   TBranch        *b_L1_DoubleEG_20_18;   //!
   TBranch        *b_L1_DoubleEG_22_10;   //!
   TBranch        *b_L1_DoubleEG_22_12;   //!
   TBranch        *b_L1_DoubleEG_22_15;   //!
   TBranch        *b_L1_DoubleEG_23_10;   //!
   TBranch        *b_L1_DoubleEG_24_17;   //!
   TBranch        *b_L1_DoubleEG_25_12;   //!
   TBranch        *b_L1_DoubleEG_25_13;   //!
   TBranch        *b_L1_DoubleEG_25_14;   //!
   TBranch        *b_L1_DoubleEG_LooseIso23_10;   //!
   TBranch        *b_L1_DoubleEG_LooseIso24_10;   //!
   TBranch        *b_L1_DoubleIsoTau28er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau30er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau32er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau33er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau34er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau35er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau36er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau38er2p1;   //!
   TBranch        *b_L1_DoubleJet100er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_DoubleJet100er2p7;   //!
   TBranch        *b_L1_DoubleJet112er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_DoubleJet112er2p7;   //!
   TBranch        *b_L1_DoubleJet120er2p7;   //!
   TBranch        *b_L1_DoubleJet150er2p7;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min300_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min320_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min340_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min360_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min380_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min400_Mu10;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min400_Mu6;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min400_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450;   //!
   TBranch        *b_L1_DoubleJet40er2p7;   //!
   TBranch        *b_L1_DoubleJet50er2p7;   //!
   TBranch        *b_L1_DoubleJet60er2p7;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM100;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM60;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM70;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM80;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM90;   //!
   TBranch        *b_L1_DoubleJet80er2p7;   //!
   TBranch        *b_L1_DoubleJet_100_30_DoubleJet30_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_100_35_DoubleJet35_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_110_40_DoubleJet40_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_115_35_DoubleJet35_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_90_30_DoubleJet30_Mass_Min620;   //!
   TBranch        *b_L1_DoubleLooseIsoEG22er2p1;   //!
   TBranch        *b_L1_DoubleLooseIsoEG24er2p1;   //!
   TBranch        *b_L1_DoubleMu0;   //!
   TBranch        *b_L1_DoubleMu0_ETM40;   //!
   TBranch        *b_L1_DoubleMu0_ETM55;   //!
   TBranch        *b_L1_DoubleMu0_ETM60;   //!
   TBranch        *b_L1_DoubleMu0_ETM65;   //!
   TBranch        *b_L1_DoubleMu0_ETM70;   //!
   TBranch        *b_L1_DoubleMu0_SQ;   //!
   TBranch        *b_L1_DoubleMu0_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu0er1p4_dEta_Max1p8_OS;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu0er2_SQ_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu18er2p1;   //!
   TBranch        *b_L1_DoubleMu22er2p1;   //!
   TBranch        *b_L1_DoubleMu3_OS_DoubleEG7p5Upsilon;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_HTT100er;   //!
   TBranch        *b_L1_DoubleMu3_SQ_HTT200er;   //!
   TBranch        *b_L1_DoubleMu3_SQ_HTT220er;   //!
   TBranch        *b_L1_DoubleMu3_SQ_HTT240er;   //!
   TBranch        *b_L1_DoubleMu4_OS_EG12;   //!
   TBranch        *b_L1_DoubleMu4_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu4_SQ_OS_dR_Max1p2;   //!
   TBranch        *b_L1_DoubleMu4p5_SQ;   //!
   TBranch        *b_L1_DoubleMu4p5_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu4p5_SQ_OS_dR_Max1p2;   //!
   TBranch        *b_L1_DoubleMu4p5er2p0_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18;   //!
   TBranch        *b_L1_DoubleMu5Upsilon_OS_DoubleEG3;   //!
   TBranch        *b_L1_DoubleMu5_OS_EG12;   //!
   TBranch        *b_L1_DoubleMu5_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu5_SQ_OS_Mass7to18;   //!
   TBranch        *b_L1_DoubleMu6_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu7_EG7;   //!
   TBranch        *b_L1_DoubleMu7_SQ_EG7;   //!
   TBranch        *b_L1_DoubleMu8_SQ;   //!
   TBranch        *b_L1_DoubleMu_10_0_dEta_Max1p8;   //!
   TBranch        *b_L1_DoubleMu_11_4;   //!
   TBranch        *b_L1_DoubleMu_12_5;   //!
   TBranch        *b_L1_DoubleMu_12_8;   //!
   TBranch        *b_L1_DoubleMu_13_6;   //!
   TBranch        *b_L1_DoubleMu_15_5;   //!
   TBranch        *b_L1_DoubleMu_15_5_SQ;   //!
   TBranch        *b_L1_DoubleMu_15_7;   //!
   TBranch        *b_L1_DoubleMu_15_7_SQ;   //!
   TBranch        *b_L1_DoubleMu_15_7_SQ_Mass_Min4;   //!
   TBranch        *b_L1_DoubleMu_20_2_SQ_Mass_Max20;   //!
   TBranch        *b_L1_DoubleTau50er2p1;   //!
   TBranch        *b_L1_DoubleTau70er2p1;   //!
   TBranch        *b_L1_EG25er2p1_HTT125er;   //!
   TBranch        *b_L1_EG27er2p1_HTT200er;   //!
   TBranch        *b_L1_ETM100;   //!
   TBranch        *b_L1_ETM100_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM105;   //!
   TBranch        *b_L1_ETM110;   //!
   TBranch        *b_L1_ETM110_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM115;   //!
   TBranch        *b_L1_ETM120;   //!
   TBranch        *b_L1_ETM150;   //!
   TBranch        *b_L1_ETM30;   //!
   TBranch        *b_L1_ETM40;   //!
   TBranch        *b_L1_ETM50;   //!
   TBranch        *b_L1_ETM60;   //!
   TBranch        *b_L1_ETM70;   //!
   TBranch        *b_L1_ETM75;   //!
   TBranch        *b_L1_ETM75_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM80;   //!
   TBranch        *b_L1_ETM80_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM85;   //!
   TBranch        *b_L1_ETM90;   //!
   TBranch        *b_L1_ETM90_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM95;   //!
   TBranch        *b_L1_ETMHF100;   //!
   TBranch        *b_L1_ETMHF100_HTT60er;   //!
   TBranch        *b_L1_ETMHF100_Jet60_OR_DiJet30woTT28;   //!
   TBranch        *b_L1_ETMHF100_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETMHF110;   //!
   TBranch        *b_L1_ETMHF110_HTT60er;   //!
   TBranch        *b_L1_ETMHF110_Jet60_OR_DiJet30woTT28;   //!
   TBranch        *b_L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETMHF120;   //!
   TBranch        *b_L1_ETMHF120_HTT60er;   //!
   TBranch        *b_L1_ETMHF120_Jet60_OR_DiJet30woTT28;   //!
   TBranch        *b_L1_ETMHF150;   //!
   TBranch        *b_L1_ETMHF70;   //!
   TBranch        *b_L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETMHF80;   //!
   TBranch        *b_L1_ETMHF80_HTT60er;   //!
   TBranch        *b_L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETMHF90;   //!
   TBranch        *b_L1_ETMHF90_HTT60er;   //!
   TBranch        *b_L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETT100_BptxAND;   //!
   TBranch        *b_L1_ETT110_BptxAND;   //!
   TBranch        *b_L1_ETT40_BptxAND;   //!
   TBranch        *b_L1_ETT50_BptxAND;   //!
   TBranch        *b_L1_ETT60_BptxAND;   //!
   TBranch        *b_L1_ETT70_BptxAND;   //!
   TBranch        *b_L1_ETT75_BptxAND;   //!
   TBranch        *b_L1_ETT80_BptxAND;   //!
   TBranch        *b_L1_ETT85_BptxAND;   //!
   TBranch        *b_L1_ETT90_BptxAND;   //!
   TBranch        *b_L1_ETT95_BptxAND;   //!
   TBranch        *b_L1_FirstBunchAfterTrain;   //!
   TBranch        *b_L1_FirstBunchInTrain;   //!
   TBranch        *b_L1_FirstCollisionInOrbit;   //!
   TBranch        *b_L1_FirstCollisionInTrain;   //!
   TBranch        *b_L1_HTT120er;   //!
   TBranch        *b_L1_HTT160er;   //!
   TBranch        *b_L1_HTT200er;   //!
   TBranch        *b_L1_HTT220er;   //!
   TBranch        *b_L1_HTT240er;   //!
   TBranch        *b_L1_HTT250er_QuadJet_70_55_40_35_er2p5;   //!
   TBranch        *b_L1_HTT255er;   //!
   TBranch        *b_L1_HTT270er;   //!
   TBranch        *b_L1_HTT280er;   //!
   TBranch        *b_L1_HTT280er_QuadJet_70_55_40_35_er2p5;   //!
   TBranch        *b_L1_HTT300er;   //!
   TBranch        *b_L1_HTT300er_QuadJet_70_55_40_35_er2p5;   //!
   TBranch        *b_L1_HTT320er;   //!
   TBranch        *b_L1_HTT320er_QuadJet_70_55_40_40_er2p4;   //!
   TBranch        *b_L1_HTT320er_QuadJet_70_55_40_40_er2p5;   //!
   TBranch        *b_L1_HTT320er_QuadJet_70_55_45_45_er2p5;   //!
   TBranch        *b_L1_HTT340er;   //!
   TBranch        *b_L1_HTT340er_QuadJet_70_55_40_40_er2p5;   //!
   TBranch        *b_L1_HTT340er_QuadJet_70_55_45_45_er2p5;   //!
   TBranch        *b_L1_HTT380er;   //!
   TBranch        *b_L1_HTT400er;   //!
   TBranch        *b_L1_HTT450er;   //!
   TBranch        *b_L1_HTT500er;   //!
   TBranch        *b_L1_IsoEG33_Mt40;   //!
   TBranch        *b_L1_IsoEG33_Mt44;   //!
   TBranch        *b_L1_IsoEG33_Mt48;   //!
   TBranch        *b_L1_IsoTau40er_ETM100;   //!
   TBranch        *b_L1_IsoTau40er_ETM105;   //!
   TBranch        *b_L1_IsoTau40er_ETM110;   //!
   TBranch        *b_L1_IsoTau40er_ETM115;   //!
   TBranch        *b_L1_IsoTau40er_ETM120;   //!
   TBranch        *b_L1_IsoTau40er_ETM80;   //!
   TBranch        *b_L1_IsoTau40er_ETM85;   //!
   TBranch        *b_L1_IsoTau40er_ETM90;   //!
   TBranch        *b_L1_IsoTau40er_ETM95;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF100;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF110;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF120;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF80;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF90;   //!
   TBranch        *b_L1_IsolatedBunch;   //!
   TBranch        *b_L1_LastCollisionInTrain;   //!
   TBranch        *b_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_HTT100er;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7;   //!
   TBranch        *b_L1_LooseIsoEG26er2p1_HTT100er;   //!
   TBranch        *b_L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG28er2p1_HTT100er;   //!
   TBranch        *b_L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3;   //!
   TBranch        *b_L1_MU20_EG15;   //!
   TBranch        *b_L1_MinimumBiasHF0_AND_BptxAND;   //!
   TBranch        *b_L1_MinimumBiasHF0_OR_BptxAND;   //!
   TBranch        *b_L1_Mu10er2p1_ETM30;   //!
   TBranch        *b_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_Mu12_EG10;   //!
   TBranch        *b_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_Mu14er2p1_ETM30;   //!
   TBranch        *b_L1_Mu15_HTT100er;   //!
   TBranch        *b_L1_Mu18_HTT100er;   //!
   TBranch        *b_L1_Mu18_Jet24er2p7;   //!
   TBranch        *b_L1_Mu18er2p1_IsoTau26er2p1;   //!
   TBranch        *b_L1_Mu18er2p1_Tau24er2p1;   //!
   TBranch        *b_L1_Mu20_EG10;   //!
   TBranch        *b_L1_Mu20_EG17;   //!
   TBranch        *b_L1_Mu20_LooseIsoEG6;   //!
   TBranch        *b_L1_Mu20er2p1_IsoTau26er2p1;   //!
   TBranch        *b_L1_Mu20er2p1_IsoTau27er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau28er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau30er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau32er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau33er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau34er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau35er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau36er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau38er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau40er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_Tau50er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_Tau70er2p1;   //!
   TBranch        *b_L1_Mu23_EG10;   //!
   TBranch        *b_L1_Mu23_LooseIsoEG10;   //!
   TBranch        *b_L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4;   //!
   TBranch        *b_L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4;   //!
   TBranch        *b_L1_Mu3_Jet30er2p5;   //!
   TBranch        *b_L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4;   //!
   TBranch        *b_L1_Mu5_EG15;   //!
   TBranch        *b_L1_Mu5_EG20;   //!
   TBranch        *b_L1_Mu5_EG23;   //!
   TBranch        *b_L1_Mu5_LooseIsoEG18;   //!
   TBranch        *b_L1_Mu5_LooseIsoEG20;   //!
   TBranch        *b_L1_Mu6_DoubleEG10;   //!
   TBranch        *b_L1_Mu6_DoubleEG17;   //!
   TBranch        *b_L1_Mu6_HTT200er;   //!
   TBranch        *b_L1_Mu6_HTT240er;   //!
   TBranch        *b_L1_Mu6_HTT250er;   //!
   TBranch        *b_L1_Mu7_EG23;   //!
   TBranch        *b_L1_Mu7_LooseIsoEG20;   //!
   TBranch        *b_L1_Mu7_LooseIsoEG23;   //!
   TBranch        *b_L1_Mu8_HTT150er;   //!
   TBranch        *b_L1_NotBptxOR;   //!
   TBranch        *b_L1_QuadJet36er2p7_IsoTau52er2p1;   //!
   TBranch        *b_L1_QuadJet36er2p7_Tau52;   //!
   TBranch        *b_L1_QuadJet40er2p7;   //!
   TBranch        *b_L1_QuadJet50er2p7;   //!
   TBranch        *b_L1_QuadJet60er2p7;   //!
   TBranch        *b_L1_QuadMu0;   //!
   TBranch        *b_L1_SingleEG10;   //!
   TBranch        *b_L1_SingleEG15;   //!
   TBranch        *b_L1_SingleEG18;   //!
   TBranch        *b_L1_SingleEG24;   //!
   TBranch        *b_L1_SingleEG26;   //!
   TBranch        *b_L1_SingleEG28;   //!
   TBranch        *b_L1_SingleEG2_BptxAND;   //!
   TBranch        *b_L1_SingleEG30;   //!
   TBranch        *b_L1_SingleEG32;   //!
   TBranch        *b_L1_SingleEG34;   //!
   TBranch        *b_L1_SingleEG34er2p1;   //!
   TBranch        *b_L1_SingleEG36;   //!
   TBranch        *b_L1_SingleEG36er2p1;   //!
   TBranch        *b_L1_SingleEG38;   //!
   TBranch        *b_L1_SingleEG38er2p1;   //!
   TBranch        *b_L1_SingleEG40;   //!
   TBranch        *b_L1_SingleEG42;   //!
   TBranch        *b_L1_SingleEG45;   //!
   TBranch        *b_L1_SingleEG5;   //!
   TBranch        *b_L1_SingleEG50;   //!
   TBranch        *b_L1_SingleIsoEG18;   //!
   TBranch        *b_L1_SingleIsoEG18er2p1;   //!
   TBranch        *b_L1_SingleIsoEG20;   //!
   TBranch        *b_L1_SingleIsoEG20er2p1;   //!
   TBranch        *b_L1_SingleIsoEG22;   //!
   TBranch        *b_L1_SingleIsoEG22er2p1;   //!
   TBranch        *b_L1_SingleIsoEG24;   //!
   TBranch        *b_L1_SingleIsoEG24er2p1;   //!
   TBranch        *b_L1_SingleIsoEG26;   //!
   TBranch        *b_L1_SingleIsoEG26er2p1;   //!
   TBranch        *b_L1_SingleIsoEG28;   //!
   TBranch        *b_L1_SingleIsoEG28er2p1;   //!
   TBranch        *b_L1_SingleIsoEG30;   //!
   TBranch        *b_L1_SingleIsoEG30er2p1;   //!
   TBranch        *b_L1_SingleIsoEG32;   //!
   TBranch        *b_L1_SingleIsoEG32er2p1;   //!
   TBranch        *b_L1_SingleIsoEG33er2p1;   //!
   TBranch        *b_L1_SingleIsoEG34;   //!
   TBranch        *b_L1_SingleIsoEG34er2p1;   //!
   TBranch        *b_L1_SingleIsoEG35;   //!
   TBranch        *b_L1_SingleIsoEG35er2p1;   //!
   TBranch        *b_L1_SingleIsoEG36;   //!
   TBranch        *b_L1_SingleIsoEG36er2p1;   //!
   TBranch        *b_L1_SingleIsoEG37;   //!
   TBranch        *b_L1_SingleIsoEG38;   //!
   TBranch        *b_L1_SingleIsoEG38er2p1;   //!
   TBranch        *b_L1_SingleIsoEG40;   //!
   TBranch        *b_L1_SingleIsoEG40er2p1;   //!
   TBranch        *b_L1_SingleJet120;   //!
   TBranch        *b_L1_SingleJet120_FWD;   //!
   TBranch        *b_L1_SingleJet12_BptxAND;   //!
   TBranch        *b_L1_SingleJet140;   //!
   TBranch        *b_L1_SingleJet150;   //!
   TBranch        *b_L1_SingleJet16;   //!
   TBranch        *b_L1_SingleJet160;   //!
   TBranch        *b_L1_SingleJet170;   //!
   TBranch        *b_L1_SingleJet180;   //!
   TBranch        *b_L1_SingleJet20;   //!
   TBranch        *b_L1_SingleJet200;   //!
   TBranch        *b_L1_SingleJet20er2p7_NotBptxOR;   //!
   TBranch        *b_L1_SingleJet20er2p7_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleJet35;   //!
   TBranch        *b_L1_SingleJet35_FWD;   //!
   TBranch        *b_L1_SingleJet35_HFm;   //!
   TBranch        *b_L1_SingleJet35_HFp;   //!
   TBranch        *b_L1_SingleJet43er2p7_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleJet46er2p7_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleJet60;   //!
   TBranch        *b_L1_SingleJet60_FWD;   //!
   TBranch        *b_L1_SingleJet60_HFm;   //!
   TBranch        *b_L1_SingleJet60_HFp;   //!
   TBranch        *b_L1_SingleJet90;   //!
   TBranch        *b_L1_SingleJet90_FWD;   //!
   TBranch        *b_L1_SingleMu0_BMTF;   //!
   TBranch        *b_L1_SingleMu0_EMTF;   //!
   TBranch        *b_L1_SingleMu0_OMTF;   //!
   TBranch        *b_L1_SingleMu10_LowQ;   //!
   TBranch        *b_L1_SingleMu11_LowQ;   //!
   TBranch        *b_L1_SingleMu12_LowQ_BMTF;   //!
   TBranch        *b_L1_SingleMu12_LowQ_EMTF;   //!
   TBranch        *b_L1_SingleMu12_LowQ_OMTF;   //!
   TBranch        *b_L1_SingleMu14er2p1;   //!
   TBranch        *b_L1_SingleMu16;   //!
   TBranch        *b_L1_SingleMu16er2p1;   //!
   TBranch        *b_L1_SingleMu18;   //!
   TBranch        *b_L1_SingleMu18er2p1;   //!
   TBranch        *b_L1_SingleMu20;   //!
   TBranch        *b_L1_SingleMu20er2p1;   //!
   TBranch        *b_L1_SingleMu22;   //!
   TBranch        *b_L1_SingleMu22_BMTF;   //!
   TBranch        *b_L1_SingleMu22_EMTF;   //!
   TBranch        *b_L1_SingleMu22_OMTF;   //!
   TBranch        *b_L1_SingleMu22er2p1;   //!
   TBranch        *b_L1_SingleMu25;   //!
   TBranch        *b_L1_SingleMu3;   //!
   TBranch        *b_L1_SingleMu30;   //!
   TBranch        *b_L1_SingleMu5;   //!
   TBranch        *b_L1_SingleMu7;   //!
   TBranch        *b_L1_SingleMuCosmics;   //!
   TBranch        *b_L1_SingleMuCosmics_BMTF;   //!
   TBranch        *b_L1_SingleMuCosmics_EMTF;   //!
   TBranch        *b_L1_SingleMuCosmics_OMTF;   //!
   TBranch        *b_L1_SingleMuOpen;   //!
   TBranch        *b_L1_SingleMuOpen_NotBptxOR;   //!
   TBranch        *b_L1_SingleMuOpen_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleTau100er2p1;   //!
   TBranch        *b_L1_SingleTau120er2p1;   //!
   TBranch        *b_L1_SingleTau130er2p1;   //!
   TBranch        *b_L1_SingleTau140er2p1;   //!
   TBranch        *b_L1_SingleTau20;   //!
   TBranch        *b_L1_SingleTau80er2p1;   //!
   TBranch        *b_L1_TripleEG_14_10_8;   //!
   TBranch        *b_L1_TripleEG_18_17_8;   //!
   TBranch        *b_L1_TripleEG_LooseIso20_10_5;   //!
   TBranch        *b_L1_TripleJet_100_85_72_VBF;   //!
   TBranch        *b_L1_TripleJet_105_85_76_VBF;   //!
   TBranch        *b_L1_TripleJet_84_68_48_VBF;   //!
   TBranch        *b_L1_TripleJet_88_72_56_VBF;   //!
   TBranch        *b_L1_TripleJet_92_76_64_VBF;   //!
   TBranch        *b_L1_TripleJet_98_83_71_VBF;   //!
   TBranch        *b_L1_TripleMu0;   //!
   TBranch        *b_L1_TripleMu0_OQ;   //!
   TBranch        *b_L1_TripleMu3;   //!
   TBranch        *b_L1_TripleMu3_SQ;   //!
   TBranch        *b_L1_TripleMu_4_4_4;   //!
   TBranch        *b_L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17;   //!
   TBranch        *b_L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0OQ;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9;   //!
   TBranch        *b_L1_TripleMu_5_0_0;   //!
   TBranch        *b_L1_TripleMu_5_3_3;   //!
   TBranch        *b_L1_TripleMu_5_3p5_2p5;   //!
   TBranch        *b_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17;   //!
   TBranch        *b_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17;   //!
   TBranch        *b_L1_TripleMu_5_5_3;   //!
   TBranch        *b_L1_UnpairedBunchBptxMinus;   //!
   TBranch        *b_L1_UnpairedBunchBptxPlus;   //!
   TBranch        *b_L1_ZeroBias;   //!
   TBranch        *b_L1_ZeroBias_copy;   //!
   TBranch        *b_L1_UnprefireableEvent;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_Flag_CSCTightHaloTrkMuUnvetoFilter;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_HcalStripHaloFilter;   //!
   TBranch        *b_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   TBranch        *b_Flag_trkPOGFilters;   //!
   TBranch        *b_Flag_chargedHadronTrackResolutionFilter;   //!
   TBranch        *b_Flag_muonBadTrackFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_BadPFMuonDzFilter;   //!
   TBranch        *b_Flag_hfNoisyHitsFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateSummer16Filter;   //!
   TBranch        *b_Flag_BadPFMuonSummer16Filter;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   TBranch        *b_Flag_METFilters;   //!
   TBranch        *b_L1Reco_step;   //!
   TBranch        *b_HLTriggerFirstPath;   //!
   TBranch        *b_HLT_AK8PFJet360_TrimMass30;   //!
   TBranch        *b_HLT_AK8PFJet380_TrimMass30;   //!
   TBranch        *b_HLT_AK8PFJet400_TrimMass30;   //!
   TBranch        *b_HLT_AK8PFJet420_TrimMass30;   //!
   TBranch        *b_HLT_AK8PFHT750_TrimMass50;   //!
   TBranch        *b_HLT_AK8PFHT800_TrimMass50;   //!
   TBranch        *b_HLT_AK8PFHT850_TrimMass50;   //!
   TBranch        *b_HLT_AK8PFHT900_TrimMass50;   //!
   TBranch        *b_HLT_CaloJet500_NoJetID;   //!
   TBranch        *b_HLT_CaloJet550_NoJetID;   //!
   TBranch        *b_HLT_Trimuon5_3p5_2_Upsilon_Muon;   //!
   TBranch        *b_HLT_DoubleEle25_CaloIdL_MW;   //!
   TBranch        *b_HLT_DoubleEle27_CaloIdL_MW;   //!
   TBranch        *b_HLT_DoubleEle33_CaloIdL_MW;   //!
   TBranch        *b_HLT_DoubleEle24_eta2p1_WPTight_Gsf;   //!
   TBranch        *b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;   //!
   TBranch        *b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350;   //!
   TBranch        *b_HLT_Ele27_Ele37_CaloIdL_MW;   //!
   TBranch        *b_HLT_Mu27_Ele37_CaloIdL_MW;   //!
   TBranch        *b_HLT_Mu37_Ele27_CaloIdL_MW;   //!
   TBranch        *b_HLT_Mu37_TkMu27;   //!
   TBranch        *b_HLT_DoubleMu4_3_Bs;   //!
   TBranch        *b_HLT_DoubleMu4_3_Jpsi_Displaced;   //!
   TBranch        *b_HLT_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced;   //!
   TBranch        *b_HLT_DoubleMu3_Trk_Tau3mu;   //!
   TBranch        *b_HLT_DoubleMu4_PsiPrimeTrk_Displaced;   //!
   TBranch        *b_HLT_DoubleMu4_Mass8_DZ_PFHT350;   //!
   TBranch        *b_HLT_DoubleMu8_Mass8_PFHT350;   //!
   TBranch        *b_HLT_Mu3_PFJet40;   //!
   TBranch        *b_HLT_Mu7p5_L2Mu2_Jpsi;   //!
   TBranch        *b_HLT_Mu7p5_L2Mu2_Upsilon;   //!
   TBranch        *b_HLT_Mu7p5_Track2_Jpsi;   //!
   TBranch        *b_HLT_Mu7p5_Track3p5_Jpsi;   //!
   TBranch        *b_HLT_Mu7p5_Track7_Jpsi;   //!
   TBranch        *b_HLT_Mu7p5_Track2_Upsilon;   //!
   TBranch        *b_HLT_Mu7p5_Track3p5_Upsilon;   //!
   TBranch        *b_HLT_Mu7p5_Track7_Upsilon;   //!
   TBranch        *b_HLT_DoublePhoton33_CaloIdL;   //!
   TBranch        *b_HLT_DoublePhoton70;   //!
   TBranch        *b_HLT_DoublePhoton85;   //!
   TBranch        *b_HLT_Ele20_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele20_WPLoose_Gsf;   //!
   TBranch        *b_HLT_Ele20_eta2p1_WPLoose_Gsf;   //!
   TBranch        *b_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;   //!
   TBranch        *b_HLT_Ele27_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele32_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele35_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele35_WPTight_Gsf_L1EGMT;   //!
   TBranch        *b_HLT_Ele38_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele40_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele32_WPTight_Gsf_L1DoubleEG;   //!
   TBranch        *b_HLT_HT450_Beamspot;   //!
   TBranch        *b_HLT_HT300_Beamspot;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1;   //!
   TBranch        *b_HLT_IsoMu20;   //!
   TBranch        *b_HLT_IsoMu24;   //!
   TBranch        *b_HLT_IsoTkMu24;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1;   //!
   TBranch        *b_HLT_IsoMu27;   //!
   TBranch        *b_HLT_IsoMu30;   //!
   TBranch        *b_HLT_UncorrectedJetE30_NoBPTX;   //!
   TBranch        *b_HLT_UncorrectedJetE30_NoBPTX3BX;   //!
   TBranch        *b_HLT_UncorrectedJetE60_NoBPTX3BX;   //!
   TBranch        *b_HLT_UncorrectedJetE70_NoBPTX3BX;   //!
   TBranch        *b_HLT_L1SingleMu18;   //!
   TBranch        *b_HLT_L1SingleMu25;   //!
   TBranch        *b_HLT_L2Mu10;   //!
   TBranch        *b_HLT_L2Mu10_NoVertex_NoBPTX3BX;   //!
   TBranch        *b_HLT_L2Mu10_NoVertex_NoBPTX;   //!
   TBranch        *b_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;   //!
   TBranch        *b_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;   //!
   TBranch        *b_HLT_L2Mu50;   //!
   TBranch        *b_HLT_DoubleL2Mu50;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_HLT_Mu25_TkMu0_Onia;   //!
   TBranch        *b_HLT_Mu30_TkMu0_Onia;   //!
   TBranch        *b_HLT_Mu20_TkMu0_Phi;   //!
   TBranch        *b_HLT_Mu25_TkMu0_Phi;   //!
   TBranch        *b_HLT_Mu20;   //!
   TBranch        *b_HLT_Mu27;   //!
   TBranch        *b_HLT_Mu50;   //!
   TBranch        *b_HLT_Mu55;   //!
   TBranch        *b_HLT_OldMu100;   //!
   TBranch        *b_HLT_TkMu100;   //!
   TBranch        *b_HLT_DiPFJet15_NoCaloMatched;   //!
   TBranch        *b_HLT_DiPFJet25_NoCaloMatched;   //!
   TBranch        *b_HLT_DiPFJet15_FBEta3_NoCaloMatched;   //!
   TBranch        *b_HLT_DiPFJet25_FBEta3_NoCaloMatched;   //!
   TBranch        *b_HLT_DiPFJetAve40;   //!
   TBranch        *b_HLT_DiPFJetAve60;   //!
   TBranch        *b_HLT_DiPFJetAve80;   //!
   TBranch        *b_HLT_DiPFJetAve140;   //!
   TBranch        *b_HLT_DiPFJetAve200;   //!
   TBranch        *b_HLT_DiPFJetAve260;   //!
   TBranch        *b_HLT_DiPFJetAve320;   //!
   TBranch        *b_HLT_DiPFJetAve400;   //!
   TBranch        *b_HLT_DiPFJetAve500;   //!
   TBranch        *b_HLT_DiPFJetAve15_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve25_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve35_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve60_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve80_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve100_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve160_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve220_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve300_HFJEC;   //!
   TBranch        *b_HLT_AK8PFJet40;   //!
   TBranch        *b_HLT_AK8PFJet60;   //!
   TBranch        *b_HLT_AK8PFJet80;   //!
   TBranch        *b_HLT_AK8PFJet140;   //!
   TBranch        *b_HLT_AK8PFJet200;   //!
   TBranch        *b_HLT_AK8PFJet260;   //!
   TBranch        *b_HLT_AK8PFJet320;   //!
   TBranch        *b_HLT_AK8PFJet400;   //!
   TBranch        *b_HLT_AK8PFJet450;   //!
   TBranch        *b_HLT_AK8PFJet500;   //!
   TBranch        *b_HLT_AK8PFJet550;   //!
   TBranch        *b_HLT_PFJet40;   //!
   TBranch        *b_HLT_PFJet60;   //!
   TBranch        *b_HLT_PFJet80;   //!
   TBranch        *b_HLT_PFJet140;   //!
   TBranch        *b_HLT_PFJet200;   //!
   TBranch        *b_HLT_PFJet260;   //!
   TBranch        *b_HLT_PFJet320;   //!
   TBranch        *b_HLT_PFJet400;   //!
   TBranch        *b_HLT_PFJet450;   //!
   TBranch        *b_HLT_PFJet500;   //!
   TBranch        *b_HLT_PFJet550;   //!
   TBranch        *b_HLT_PFJetFwd40;   //!
   TBranch        *b_HLT_PFJetFwd60;   //!
   TBranch        *b_HLT_PFJetFwd80;   //!
   TBranch        *b_HLT_PFJetFwd140;   //!
   TBranch        *b_HLT_PFJetFwd200;   //!
   TBranch        *b_HLT_PFJetFwd260;   //!
   TBranch        *b_HLT_PFJetFwd320;   //!
   TBranch        *b_HLT_PFJetFwd400;   //!
   TBranch        *b_HLT_PFJetFwd450;   //!
   TBranch        *b_HLT_PFJetFwd500;   //!
   TBranch        *b_HLT_AK8PFJetFwd40;   //!
   TBranch        *b_HLT_AK8PFJetFwd60;   //!
   TBranch        *b_HLT_AK8PFJetFwd80;   //!
   TBranch        *b_HLT_AK8PFJetFwd140;   //!
   TBranch        *b_HLT_AK8PFJetFwd200;   //!
   TBranch        *b_HLT_AK8PFJetFwd260;   //!
   TBranch        *b_HLT_AK8PFJetFwd320;   //!
   TBranch        *b_HLT_AK8PFJetFwd400;   //!
   TBranch        *b_HLT_AK8PFJetFwd450;   //!
   TBranch        *b_HLT_AK8PFJetFwd500;   //!
   TBranch        *b_HLT_PFHT180;   //!
   TBranch        *b_HLT_PFHT250;   //!
   TBranch        *b_HLT_PFHT370;   //!
   TBranch        *b_HLT_PFHT430;   //!
   TBranch        *b_HLT_PFHT510;   //!
   TBranch        *b_HLT_PFHT590;   //!
   TBranch        *b_HLT_PFHT680;   //!
   TBranch        *b_HLT_PFHT780;   //!
   TBranch        *b_HLT_PFHT890;   //!
   TBranch        *b_HLT_PFHT1050;   //!
   TBranch        *b_HLT_PFHT500_PFMET100_PFMHT100_IDTight;   //!
   TBranch        *b_HLT_PFHT500_PFMET110_PFMHT110_IDTight;   //!
   TBranch        *b_HLT_PFHT700_PFMET85_PFMHT85_IDTight;   //!
   TBranch        *b_HLT_PFHT700_PFMET95_PFMHT95_IDTight;   //!
   TBranch        *b_HLT_PFHT800_PFMET75_PFMHT75_IDTight;   //!
   TBranch        *b_HLT_PFHT800_PFMET85_PFMHT85_IDTight;   //!
   TBranch        *b_HLT_PFMET110_PFMHT110_IDTight;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight;   //!
   TBranch        *b_HLT_PFMET130_PFMHT130_IDTight;   //!
   TBranch        *b_HLT_PFMET140_PFMHT140_IDTight;   //!
   TBranch        *b_HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1;   //!
   TBranch        *b_HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1;   //!
   TBranch        *b_HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1;   //!
   TBranch        *b_HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight_PFHT60;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;   //!
   TBranch        *b_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60;   //!
   TBranch        *b_HLT_PFMETTypeOne110_PFMHT110_IDTight;   //!
   TBranch        *b_HLT_PFMETTypeOne120_PFMHT120_IDTight;   //!
   TBranch        *b_HLT_PFMETTypeOne130_PFMHT130_IDTight;   //!
   TBranch        *b_HLT_PFMETTypeOne140_PFMHT140_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;   //!
   TBranch        *b_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight;   //!
   TBranch        *b_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight;   //!
   TBranch        *b_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight;   //!
   TBranch        *b_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight;   //!
   TBranch        *b_HLT_L1ETMHadSeeds;   //!
   TBranch        *b_HLT_CaloMHT90;   //!
   TBranch        *b_HLT_CaloMET80_NotCleaned;   //!
   TBranch        *b_HLT_CaloMET90_NotCleaned;   //!
   TBranch        *b_HLT_CaloMET100_NotCleaned;   //!
   TBranch        *b_HLT_CaloMET110_NotCleaned;   //!
   TBranch        *b_HLT_CaloMET250_NotCleaned;   //!
   TBranch        *b_HLT_CaloMET70_HBHECleaned;   //!
   TBranch        *b_HLT_CaloMET80_HBHECleaned;   //!
   TBranch        *b_HLT_CaloMET90_HBHECleaned;   //!
   TBranch        *b_HLT_CaloMET100_HBHECleaned;   //!
   TBranch        *b_HLT_CaloMET250_HBHECleaned;   //!
   TBranch        *b_HLT_CaloMET300_HBHECleaned;   //!
   TBranch        *b_HLT_CaloMET350_HBHECleaned;   //!
   TBranch        *b_HLT_PFMET200_NotCleaned;   //!
   TBranch        *b_HLT_PFMET200_HBHECleaned;   //!
   TBranch        *b_HLT_PFMET250_HBHECleaned;   //!
   TBranch        *b_HLT_PFMET300_HBHECleaned;   //!
   TBranch        *b_HLT_PFMET200_HBHE_BeamHaloCleaned;   //!
   TBranch        *b_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned;   //!
   TBranch        *b_HLT_MET105_IsoTrk50;   //!
   TBranch        *b_HLT_MET120_IsoTrk50;   //!
   TBranch        *b_HLT_SingleJet30_Mu12_SinglePFJet40;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33;   //!
   TBranch        *b_HLT_DoublePFJets40_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_DoublePFJets100_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_DoublePFJets200_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_DoublePFJets350_CaloBTagCSV_p33;   //!
   TBranch        *b_HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33;   //!
   TBranch        *b_HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33;   //!
   TBranch        *b_HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33;   //!
   TBranch        *b_HLT_Photon300_NoHE;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;   //!
   TBranch        *b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;   //!
   TBranch        *b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet20_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet40_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet70_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet110_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet170_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4Jet300_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK8DiJet170_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK8Jet300_Mu5;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Mu12_DoublePhoton20;   //!
   TBranch        *b_HLT_TriplePhoton_20_20_20_CaloIdLV2;   //!
   TBranch        *b_HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL;   //!
   TBranch        *b_HLT_TriplePhoton_30_30_10_CaloIdLV2;   //!
   TBranch        *b_HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL;   //!
   TBranch        *b_HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL;   //!
   TBranch        *b_HLT_Photon25;   //!
   TBranch        *b_HLT_Photon33;   //!
   TBranch        *b_HLT_Photon50;   //!
   TBranch        *b_HLT_Photon75;   //!
   TBranch        *b_HLT_Photon90;   //!
   TBranch        *b_HLT_Photon120;   //!
   TBranch        *b_HLT_Photon150;   //!
   TBranch        *b_HLT_Photon175;   //!
   TBranch        *b_HLT_Photon200;   //!
   TBranch        *b_HLT_Photon50_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon75_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon90_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon120_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon165_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon90_CaloIdL_PFHT700;   //!
   TBranch        *b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;   //!
   TBranch        *b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95;   //!
   TBranch        *b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;   //!
   TBranch        *b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;   //!
   TBranch        *b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;   //!
   TBranch        *b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_L1_NoOS;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_NoVertexing_NoOS;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_NoVertexing;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi3p5_Muon2;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_4p5;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_5;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_4p5NoOS;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_4p5er2p0;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_4p5er2p0M;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_NoVertexing;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_5M;   //!
   TBranch        *b_HLT_Dimuon0_LowMass_L1_0er1p5R;   //!
   TBranch        *b_HLT_Dimuon0_LowMass_L1_0er1p5;   //!
   TBranch        *b_HLT_Dimuon0_LowMass;   //!
   TBranch        *b_HLT_Dimuon0_LowMass_L1_4;   //!
   TBranch        *b_HLT_Dimuon0_LowMass_L1_4R;   //!
   TBranch        *b_HLT_Dimuon0_LowMass_L1_TM530;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_L1_TM0;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_NoL1Mass;   //!
   TBranch        *b_HLT_TripleMu_5_3_3_Mass3p8to60_DZ;   //!
   TBranch        *b_HLT_TripleMu_10_5_5_DZ;   //!
   TBranch        *b_HLT_TripleMu_12_10_5;   //!
   TBranch        *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15;   //!
   TBranch        *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1;   //!
   TBranch        *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15;   //!
   TBranch        *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1;   //!
   TBranch        *b_HLT_DoubleMu3_DZ_PFMET50_PFMHT60;   //!
   TBranch        *b_HLT_DoubleMu3_DZ_PFMET70_PFMHT70;   //!
   TBranch        *b_HLT_DoubleMu3_DZ_PFMET90_PFMHT90;   //!
   TBranch        *b_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass;   //!
   TBranch        *b_HLT_DoubleMu4_Jpsi_Displaced;   //!
   TBranch        *b_HLT_DoubleMu4_Jpsi_NoVertexing;   //!
   TBranch        *b_HLT_DoubleMu4_JpsiTrkTrk_Displaced;   //!
   TBranch        *b_HLT_DoubleMu43NoFiltersNoVtx;   //!
   TBranch        *b_HLT_DoubleMu48NoFiltersNoVtx;   //!
   TBranch        *b_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL;   //!
   TBranch        *b_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL;   //!
   TBranch        *b_HLT_DoubleMu20_7_Mass0to30_L1_DM4;   //!
   TBranch        *b_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG;   //!
   TBranch        *b_HLT_HT425;   //!
   TBranch        *b_HLT_HT430_DisplacedDijet40_DisplacedTrack;   //!
   TBranch        *b_HLT_HT430_DisplacedDijet60_DisplacedTrack;   //!
   TBranch        *b_HLT_HT430_DisplacedDijet80_DisplacedTrack;   //!
   TBranch        *b_HLT_HT400_DisplacedDijet40_DisplacedTrack;   //!
   TBranch        *b_HLT_HT650_DisplacedDijet60_Inclusive;   //!
   TBranch        *b_HLT_HT550_DisplacedDijet80_Inclusive;   //!
   TBranch        *b_HLT_HT550_DisplacedDijet60_Inclusive;   //!
   TBranch        *b_HLT_HT650_DisplacedDijet80_Inclusive;   //!
   TBranch        *b_HLT_HT750_DisplacedDijet80_Inclusive;   //!
   TBranch        *b_HLT_DiJet110_35_Mjj650_PFMET110;   //!
   TBranch        *b_HLT_DiJet110_35_Mjj650_PFMET120;   //!
   TBranch        *b_HLT_DiJet110_35_Mjj650_PFMET130;   //!
   TBranch        *b_HLT_TripleJet110_35_35_Mjj650_PFMET110;   //!
   TBranch        *b_HLT_TripleJet110_35_35_Mjj650_PFMET120;   //!
   TBranch        *b_HLT_TripleJet110_35_35_Mjj650_PFMET130;   //!
   TBranch        *b_HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned;   //!
   TBranch        *b_HLT_Ele28_eta2p1_WPTight_Gsf_HT150;   //!
   TBranch        *b_HLT_Ele28_HighEta_SC20_Mass55;   //!
   TBranch        *b_HLT_DoubleMu20_7_Mass0to30_Photon23;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT450_PFMET50;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT450;   //!
   TBranch        *b_HLT_Ele50_IsoVVVL_PFHT450;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT600;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;   //!
   TBranch        *b_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT450_PFMET50;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT450;   //!
   TBranch        *b_HLT_Mu50_IsoVVVL_PFHT450;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT600;   //!
   TBranch        *b_HLT_Dimuon10_PsiPrime_Barrel_Seagulls;   //!
   TBranch        *b_HLT_Dimuon20_Jpsi_Barrel_Seagulls;   //!
   TBranch        *b_HLT_Dimuon10_Upsilon_Barrel_Seagulls;   //!
   TBranch        *b_HLT_Dimuon12_Upsilon_eta1p5;   //!
   TBranch        *b_HLT_Dimuon14_Phi_Barrel_Seagulls;   //!
   TBranch        *b_HLT_Dimuon18_PsiPrime;   //!
   TBranch        *b_HLT_Dimuon25_Jpsi;   //!
   TBranch        *b_HLT_Dimuon18_PsiPrime_noCorrL1;   //!
   TBranch        *b_HLT_Dimuon24_Upsilon_noCorrL1;   //!
   TBranch        *b_HLT_Dimuon24_Phi_noCorrL1;   //!
   TBranch        *b_HLT_Dimuon25_Jpsi_noCorrL1;   //!
   TBranch        *b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;   //!
   TBranch        *b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT_DoubleIsoMu20_eta2p1;   //!
   TBranch        *b_HLT_DoubleIsoMu24_eta2p1;   //!
   TBranch        *b_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx;   //!
   TBranch        *b_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx;   //!
   TBranch        *b_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx;   //!
   TBranch        *b_HLT_Mu8;   //!
   TBranch        *b_HLT_Mu17;   //!
   TBranch        *b_HLT_Mu19;   //!
   TBranch        *b_HLT_Mu17_Photon30_IsoCaloId;   //!
   TBranch        *b_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele17_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;   //!
   TBranch        *b_HLT_Ele115_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_Ele135_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_Ele145_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_Ele200_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_Ele250_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_Ele300_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_PFHT300PT30_QuadPFJet_75_60_45_40;   //!
   TBranch        *b_HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0;   //!
   TBranch        *b_HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2;   //!
   TBranch        *b_HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2;   //!
   TBranch        *b_HLT_PFHT380_SixPFJet32;   //!
   TBranch        *b_HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5;   //!
   TBranch        *b_HLT_PFHT430_SixPFJet40;   //!
   TBranch        *b_HLT_PFHT350;   //!
   TBranch        *b_HLT_PFHT350MinPFJet15;   //!
   TBranch        *b_HLT_Photon60_R9Id90_CaloIdL_IsoL;   //!
   TBranch        *b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL;   //!
   TBranch        *b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15;   //!
   TBranch        *b_HLT_FullTrack_Multiplicity85;   //!
   TBranch        *b_HLT_FullTrack_Multiplicity100;   //!
   TBranch        *b_HLT_FullTrack_Multiplicity130;   //!
   TBranch        *b_HLT_FullTrack_Multiplicity155;   //!
   TBranch        *b_HLT_ECALHT800;   //!
   TBranch        *b_HLT_DiSC30_18_EIso_AND_HE_Mass70;   //!
   TBranch        *b_HLT_Physics;   //!
   TBranch        *b_HLT_Physics_part0;   //!
   TBranch        *b_HLT_Physics_part1;   //!
   TBranch        *b_HLT_Physics_part2;   //!
   TBranch        *b_HLT_Physics_part3;   //!
   TBranch        *b_HLT_Physics_part4;   //!
   TBranch        *b_HLT_Physics_part5;   //!
   TBranch        *b_HLT_Physics_part6;   //!
   TBranch        *b_HLT_Physics_part7;   //!
   TBranch        *b_HLT_Random;   //!
   TBranch        *b_HLT_ZeroBias;   //!
   TBranch        *b_HLT_ZeroBias_part0;   //!
   TBranch        *b_HLT_ZeroBias_part1;   //!
   TBranch        *b_HLT_ZeroBias_part2;   //!
   TBranch        *b_HLT_ZeroBias_part3;   //!
   TBranch        *b_HLT_ZeroBias_part4;   //!
   TBranch        *b_HLT_ZeroBias_part5;   //!
   TBranch        *b_HLT_ZeroBias_part6;   //!
   TBranch        *b_HLT_ZeroBias_part7;   //!
   TBranch        *b_HLT_AK4CaloJet30;   //!
   TBranch        *b_HLT_AK4CaloJet40;   //!
   TBranch        *b_HLT_AK4CaloJet50;   //!
   TBranch        *b_HLT_AK4CaloJet80;   //!
   TBranch        *b_HLT_AK4CaloJet100;   //!
   TBranch        *b_HLT_AK4CaloJet120;   //!
   TBranch        *b_HLT_AK4PFJet30;   //!
   TBranch        *b_HLT_AK4PFJet50;   //!
   TBranch        *b_HLT_AK4PFJet80;   //!
   TBranch        *b_HLT_AK4PFJet100;   //!
   TBranch        *b_HLT_AK4PFJet120;   //!
   TBranch        *b_HLT_HISinglePhoton10_Eta3p1ForPPRef;   //!
   TBranch        *b_HLT_HISinglePhoton20_Eta3p1ForPPRef;   //!
   TBranch        *b_HLT_HISinglePhoton30_Eta3p1ForPPRef;   //!
   TBranch        *b_HLT_HISinglePhoton40_Eta3p1ForPPRef;   //!
   TBranch        *b_HLT_HISinglePhoton50_Eta3p1ForPPRef;   //!
   TBranch        *b_HLT_HISinglePhoton60_Eta3p1ForPPRef;   //!
   TBranch        *b_HLT_Photon20_HoverELoose;   //!
   TBranch        *b_HLT_Photon30_HoverELoose;   //!
   TBranch        *b_HLT_Photon40_HoverELoose;   //!
   TBranch        *b_HLT_Photon50_HoverELoose;   //!
   TBranch        *b_HLT_Photon60_HoverELoose;   //!
   TBranch        *b_HLT_EcalCalibration;   //!
   TBranch        *b_HLT_HcalCalibration;   //!
   TBranch        *b_HLT_L1UnpairedBunchBptxMinus;   //!
   TBranch        *b_HLT_L1UnpairedBunchBptxPlus;   //!
   TBranch        *b_HLT_L1NotBptxOR;   //!
   TBranch        *b_HLT_L1MinimumBiasHF_OR;   //!
   TBranch        *b_HLT_L1MinimumBiasHF0OR;   //!
   TBranch        *b_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;   //!
   TBranch        *b_HLT_HcalNZS;   //!
   TBranch        *b_HLT_HcalPhiSym;   //!
   TBranch        *b_HLT_HcalIsolatedbunch;   //!
   TBranch        *b_HLT_IsoTrackHB;   //!
   TBranch        *b_HLT_IsoTrackHE;   //!
   TBranch        *b_HLT_ZeroBias_FirstCollisionAfterAbortGap;   //!
   TBranch        *b_HLT_ZeroBias_IsolatedBunches;   //!
   TBranch        *b_HLT_ZeroBias_FirstCollisionInTrain;   //!
   TBranch        *b_HLT_ZeroBias_LastCollisionInTrain;   //!
   TBranch        *b_HLT_ZeroBias_FirstBXAfterTrain;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr;   //!
   TBranch        *b_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT_Rsq0p35;   //!
   TBranch        *b_HLT_Rsq0p40;   //!
   TBranch        *b_HLT_RsqMR300_Rsq0p09_MR200;   //!
   TBranch        *b_HLT_RsqMR320_Rsq0p09_MR200;   //!
   TBranch        *b_HLT_RsqMR300_Rsq0p09_MR200_4jet;   //!
   TBranch        *b_HLT_RsqMR320_Rsq0p09_MR200_4jet;   //!
   TBranch        *b_HLT_L1_DoubleJet30_Mass_Min400_Mu10;   //!
   TBranch        *b_HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1;   //!
   TBranch        *b_HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1;   //!
   TBranch        *b_HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1;   //!
   TBranch        *b_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50;   //!
   TBranch        *b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;   //!
   TBranch        *b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3;   //!
   TBranch        *b_HLT_PFMET100_PFMHT100_IDTight_PFHT60;   //!
   TBranch        *b_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;   //!
   TBranch        *b_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60;   //!
   TBranch        *b_HLT_Mu18_Mu9_SameSign;   //!
   TBranch        *b_HLT_Mu18_Mu9_SameSign_DZ;   //!
   TBranch        *b_HLT_Mu18_Mu9;   //!
   TBranch        *b_HLT_Mu18_Mu9_DZ;   //!
   TBranch        *b_HLT_Mu20_Mu10_SameSign;   //!
   TBranch        *b_HLT_Mu20_Mu10_SameSign_DZ;   //!
   TBranch        *b_HLT_Mu20_Mu10;   //!
   TBranch        *b_HLT_Mu20_Mu10_DZ;   //!
   TBranch        *b_HLT_Mu23_Mu12_SameSign;   //!
   TBranch        *b_HLT_Mu23_Mu12_SameSign_DZ;   //!
   TBranch        *b_HLT_Mu23_Mu12;   //!
   TBranch        *b_HLT_Mu23_Mu12_DZ;   //!
   TBranch        *b_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi;   //!
   TBranch        *b_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi;   //!
   TBranch        *b_HLT_DoubleMu3_DCA_PFMET50_PFMHT60;   //!
   TBranch        *b_HLT_TripleMu_5_3_3_Mass3p8to60_DCA;   //!
   TBranch        *b_HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1;   //!
   TBranch        *b_HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1;   //!
   TBranch        *b_HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1;   //!
   TBranch        *b_HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1;   //!
   TBranch        *b_HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2;   //!
   TBranch        *b_HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2;   //!
   TBranch        *b_HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2;   //!
   TBranch        *b_HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2;   //!
   TBranch        *b_HLT_QuadPFJet98_83_71_15;   //!
   TBranch        *b_HLT_QuadPFJet103_88_75_15;   //!
   TBranch        *b_HLT_QuadPFJet105_88_76_15;   //!
   TBranch        *b_HLT_QuadPFJet111_90_80_15;   //!
   TBranch        *b_HLT_AK8PFJet330_PFAK8BTagCSV_p17;   //!
   TBranch        *b_HLT_AK8PFJet330_PFAK8BTagCSV_p1;   //!
   TBranch        *b_HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;   //!
   TBranch        *b_HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   TBranch        *b_L1simulation_step;   //!
   TBranch        *b_Jet_pt_raw;   //!
   TBranch        *b_Jet_pt_nom;   //!
   TBranch        *b_Jet_mass_raw;   //!
   TBranch        *b_Jet_mass_nom;   //!
   TBranch        *b_Jet_corr_JEC;   //!
   TBranch        *b_Jet_corr_JER;   //!
   TBranch        *b_MET_T1_pt;   //!
   TBranch        *b_MET_T1_phi;   //!
   TBranch        *b_MET_T1Smear_pt;   //!
   TBranch        *b_MET_T1Smear_phi;   //!
   TBranch        *b_Jet_pt_jer0Up;   //!
   TBranch        *b_Jet_mass_jer0Up;   //!
   TBranch        *b_MET_T1_pt_jer0Up;   //!
   TBranch        *b_MET_T1_phi_jer0Up;   //!
   TBranch        *b_MET_T1Smear_pt_jer0Up;   //!
   TBranch        *b_MET_T1Smear_phi_jer0Up;   //!
   TBranch        *b_Jet_pt_jer1Up;   //!
   TBranch        *b_Jet_mass_jer1Up;   //!
   TBranch        *b_MET_T1_pt_jer1Up;   //!
   TBranch        *b_MET_T1_phi_jer1Up;   //!
   TBranch        *b_MET_T1Smear_pt_jer1Up;   //!
   TBranch        *b_MET_T1Smear_phi_jer1Up;   //!
   TBranch        *b_Jet_pt_jer2Up;   //!
   TBranch        *b_Jet_mass_jer2Up;   //!
   TBranch        *b_MET_T1_pt_jer2Up;   //!
   TBranch        *b_MET_T1_phi_jer2Up;   //!
   TBranch        *b_MET_T1Smear_pt_jer2Up;   //!
   TBranch        *b_MET_T1Smear_phi_jer2Up;   //!
   TBranch        *b_Jet_pt_jer3Up;   //!
   TBranch        *b_Jet_mass_jer3Up;   //!
   TBranch        *b_MET_T1_pt_jer3Up;   //!
   TBranch        *b_MET_T1_phi_jer3Up;   //!
   TBranch        *b_MET_T1Smear_pt_jer3Up;   //!
   TBranch        *b_MET_T1Smear_phi_jer3Up;   //!
   TBranch        *b_Jet_pt_jer4Up;   //!
   TBranch        *b_Jet_mass_jer4Up;   //!
   TBranch        *b_MET_T1_pt_jer4Up;   //!
   TBranch        *b_MET_T1_phi_jer4Up;   //!
   TBranch        *b_MET_T1Smear_pt_jer4Up;   //!
   TBranch        *b_MET_T1Smear_phi_jer4Up;   //!
   TBranch        *b_Jet_pt_jer5Up;   //!
   TBranch        *b_Jet_mass_jer5Up;   //!
   TBranch        *b_MET_T1_pt_jer5Up;   //!
   TBranch        *b_MET_T1_phi_jer5Up;   //!
   TBranch        *b_MET_T1Smear_pt_jer5Up;   //!
   TBranch        *b_MET_T1Smear_phi_jer5Up;   //!
   TBranch        *b_Jet_pt_jesAbsoluteStatUp;   //!
   TBranch        *b_Jet_mass_jesAbsoluteStatUp;   //!
   TBranch        *b_MET_T1_pt_jesAbsoluteStatUp;   //!
   TBranch        *b_MET_T1_phi_jesAbsoluteStatUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesAbsoluteStatUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesAbsoluteStatUp;   //!
   TBranch        *b_Jet_pt_jesAbsoluteScaleUp;   //!
   TBranch        *b_Jet_mass_jesAbsoluteScaleUp;   //!
   TBranch        *b_MET_T1_pt_jesAbsoluteScaleUp;   //!
   TBranch        *b_MET_T1_phi_jesAbsoluteScaleUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesAbsoluteScaleUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesAbsoluteScaleUp;   //!
   TBranch        *b_Jet_pt_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_Jet_mass_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_MET_T1_pt_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_MET_T1_phi_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesAbsoluteFlavMapUp;   //!
   TBranch        *b_Jet_pt_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_Jet_mass_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_MET_T1_pt_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_MET_T1_phi_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesAbsoluteMPFBiasUp;   //!
   TBranch        *b_Jet_pt_jesFragmentationUp;   //!
   TBranch        *b_Jet_mass_jesFragmentationUp;   //!
   TBranch        *b_MET_T1_pt_jesFragmentationUp;   //!
   TBranch        *b_MET_T1_phi_jesFragmentationUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesFragmentationUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesFragmentationUp;   //!
   TBranch        *b_Jet_pt_jesSinglePionECALUp;   //!
   TBranch        *b_Jet_mass_jesSinglePionECALUp;   //!
   TBranch        *b_MET_T1_pt_jesSinglePionECALUp;   //!
   TBranch        *b_MET_T1_phi_jesSinglePionECALUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesSinglePionECALUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesSinglePionECALUp;   //!
   TBranch        *b_Jet_pt_jesSinglePionHCALUp;   //!
   TBranch        *b_Jet_mass_jesSinglePionHCALUp;   //!
   TBranch        *b_MET_T1_pt_jesSinglePionHCALUp;   //!
   TBranch        *b_MET_T1_phi_jesSinglePionHCALUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesSinglePionHCALUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesSinglePionHCALUp;   //!
   TBranch        *b_Jet_pt_jesFlavorQCDUp;   //!
   TBranch        *b_Jet_mass_jesFlavorQCDUp;   //!
   TBranch        *b_MET_T1_pt_jesFlavorQCDUp;   //!
   TBranch        *b_MET_T1_phi_jesFlavorQCDUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesFlavorQCDUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesFlavorQCDUp;   //!
   TBranch        *b_Jet_pt_jesTimePtEtaUp;   //!
   TBranch        *b_Jet_mass_jesTimePtEtaUp;   //!
   TBranch        *b_MET_T1_pt_jesTimePtEtaUp;   //!
   TBranch        *b_MET_T1_phi_jesTimePtEtaUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesTimePtEtaUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesTimePtEtaUp;   //!
   TBranch        *b_Jet_pt_jesRelativeJEREC1Up;   //!
   TBranch        *b_Jet_mass_jesRelativeJEREC1Up;   //!
   TBranch        *b_MET_T1_pt_jesRelativeJEREC1Up;   //!
   TBranch        *b_MET_T1_phi_jesRelativeJEREC1Up;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeJEREC1Up;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeJEREC1Up;   //!
   TBranch        *b_Jet_pt_jesRelativeJEREC2Up;   //!
   TBranch        *b_Jet_mass_jesRelativeJEREC2Up;   //!
   TBranch        *b_MET_T1_pt_jesRelativeJEREC2Up;   //!
   TBranch        *b_MET_T1_phi_jesRelativeJEREC2Up;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeJEREC2Up;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeJEREC2Up;   //!
   TBranch        *b_Jet_pt_jesRelativeJERHFUp;   //!
   TBranch        *b_Jet_mass_jesRelativeJERHFUp;   //!
   TBranch        *b_MET_T1_pt_jesRelativeJERHFUp;   //!
   TBranch        *b_MET_T1_phi_jesRelativeJERHFUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeJERHFUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeJERHFUp;   //!
   TBranch        *b_Jet_pt_jesRelativePtBBUp;   //!
   TBranch        *b_Jet_mass_jesRelativePtBBUp;   //!
   TBranch        *b_MET_T1_pt_jesRelativePtBBUp;   //!
   TBranch        *b_MET_T1_phi_jesRelativePtBBUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativePtBBUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativePtBBUp;   //!
   TBranch        *b_Jet_pt_jesRelativePtEC1Up;   //!
   TBranch        *b_Jet_mass_jesRelativePtEC1Up;   //!
   TBranch        *b_MET_T1_pt_jesRelativePtEC1Up;   //!
   TBranch        *b_MET_T1_phi_jesRelativePtEC1Up;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativePtEC1Up;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativePtEC1Up;   //!
   TBranch        *b_Jet_pt_jesRelativePtEC2Up;   //!
   TBranch        *b_Jet_mass_jesRelativePtEC2Up;   //!
   TBranch        *b_MET_T1_pt_jesRelativePtEC2Up;   //!
   TBranch        *b_MET_T1_phi_jesRelativePtEC2Up;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativePtEC2Up;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativePtEC2Up;   //!
   TBranch        *b_Jet_pt_jesRelativePtHFUp;   //!
   TBranch        *b_Jet_mass_jesRelativePtHFUp;   //!
   TBranch        *b_MET_T1_pt_jesRelativePtHFUp;   //!
   TBranch        *b_MET_T1_phi_jesRelativePtHFUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativePtHFUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativePtHFUp;   //!
   TBranch        *b_Jet_pt_jesRelativeBalUp;   //!
   TBranch        *b_Jet_mass_jesRelativeBalUp;   //!
   TBranch        *b_MET_T1_pt_jesRelativeBalUp;   //!
   TBranch        *b_MET_T1_phi_jesRelativeBalUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeBalUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeBalUp;   //!
   TBranch        *b_Jet_pt_jesRelativeSampleUp;   //!
   TBranch        *b_Jet_mass_jesRelativeSampleUp;   //!
   TBranch        *b_MET_T1_pt_jesRelativeSampleUp;   //!
   TBranch        *b_MET_T1_phi_jesRelativeSampleUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeSampleUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeSampleUp;   //!
   TBranch        *b_Jet_pt_jesRelativeFSRUp;   //!
   TBranch        *b_Jet_mass_jesRelativeFSRUp;   //!
   TBranch        *b_MET_T1_pt_jesRelativeFSRUp;   //!
   TBranch        *b_MET_T1_phi_jesRelativeFSRUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeFSRUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeFSRUp;   //!
   TBranch        *b_Jet_pt_jesRelativeStatFSRUp;   //!
   TBranch        *b_Jet_mass_jesRelativeStatFSRUp;   //!
   TBranch        *b_MET_T1_pt_jesRelativeStatFSRUp;   //!
   TBranch        *b_MET_T1_phi_jesRelativeStatFSRUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeStatFSRUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeStatFSRUp;   //!
   TBranch        *b_Jet_pt_jesRelativeStatECUp;   //!
   TBranch        *b_Jet_mass_jesRelativeStatECUp;   //!
   TBranch        *b_MET_T1_pt_jesRelativeStatECUp;   //!
   TBranch        *b_MET_T1_phi_jesRelativeStatECUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeStatECUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeStatECUp;   //!
   TBranch        *b_Jet_pt_jesRelativeStatHFUp;   //!
   TBranch        *b_Jet_mass_jesRelativeStatHFUp;   //!
   TBranch        *b_MET_T1_pt_jesRelativeStatHFUp;   //!
   TBranch        *b_MET_T1_phi_jesRelativeStatHFUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeStatHFUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeStatHFUp;   //!
   TBranch        *b_Jet_pt_jesPileUpDataMCUp;   //!
   TBranch        *b_Jet_mass_jesPileUpDataMCUp;   //!
   TBranch        *b_MET_T1_pt_jesPileUpDataMCUp;   //!
   TBranch        *b_MET_T1_phi_jesPileUpDataMCUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpDataMCUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpDataMCUp;   //!
   TBranch        *b_Jet_pt_jesPileUpPtRefUp;   //!
   TBranch        *b_Jet_mass_jesPileUpPtRefUp;   //!
   TBranch        *b_MET_T1_pt_jesPileUpPtRefUp;   //!
   TBranch        *b_MET_T1_phi_jesPileUpPtRefUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpPtRefUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpPtRefUp;   //!
   TBranch        *b_Jet_pt_jesPileUpPtBBUp;   //!
   TBranch        *b_Jet_mass_jesPileUpPtBBUp;   //!
   TBranch        *b_MET_T1_pt_jesPileUpPtBBUp;   //!
   TBranch        *b_MET_T1_phi_jesPileUpPtBBUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpPtBBUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpPtBBUp;   //!
   TBranch        *b_Jet_pt_jesPileUpPtEC1Up;   //!
   TBranch        *b_Jet_mass_jesPileUpPtEC1Up;   //!
   TBranch        *b_MET_T1_pt_jesPileUpPtEC1Up;   //!
   TBranch        *b_MET_T1_phi_jesPileUpPtEC1Up;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpPtEC1Up;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpPtEC1Up;   //!
   TBranch        *b_Jet_pt_jesPileUpPtEC2Up;   //!
   TBranch        *b_Jet_mass_jesPileUpPtEC2Up;   //!
   TBranch        *b_MET_T1_pt_jesPileUpPtEC2Up;   //!
   TBranch        *b_MET_T1_phi_jesPileUpPtEC2Up;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpPtEC2Up;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpPtEC2Up;   //!
   TBranch        *b_Jet_pt_jesPileUpPtHFUp;   //!
   TBranch        *b_Jet_mass_jesPileUpPtHFUp;   //!
   TBranch        *b_MET_T1_pt_jesPileUpPtHFUp;   //!
   TBranch        *b_MET_T1_phi_jesPileUpPtHFUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpPtHFUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpPtHFUp;   //!
   TBranch        *b_Jet_pt_jesPileUpMuZeroUp;   //!
   TBranch        *b_Jet_mass_jesPileUpMuZeroUp;   //!
   TBranch        *b_MET_T1_pt_jesPileUpMuZeroUp;   //!
   TBranch        *b_MET_T1_phi_jesPileUpMuZeroUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpMuZeroUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpMuZeroUp;   //!
   TBranch        *b_Jet_pt_jesPileUpEnvelopeUp;   //!
   TBranch        *b_Jet_mass_jesPileUpEnvelopeUp;   //!
   TBranch        *b_MET_T1_pt_jesPileUpEnvelopeUp;   //!
   TBranch        *b_MET_T1_phi_jesPileUpEnvelopeUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpEnvelopeUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpEnvelopeUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalPileUpUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalPileUpUp;   //!
   TBranch        *b_MET_T1_pt_jesSubTotalPileUpUp;   //!
   TBranch        *b_MET_T1_phi_jesSubTotalPileUpUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesSubTotalPileUpUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesSubTotalPileUpUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalRelativeUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalRelativeUp;   //!
   TBranch        *b_MET_T1_pt_jesSubTotalRelativeUp;   //!
   TBranch        *b_MET_T1_phi_jesSubTotalRelativeUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesSubTotalRelativeUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesSubTotalRelativeUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalPtUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalPtUp;   //!
   TBranch        *b_MET_T1_pt_jesSubTotalPtUp;   //!
   TBranch        *b_MET_T1_phi_jesSubTotalPtUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesSubTotalPtUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesSubTotalPtUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalScaleUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalScaleUp;   //!
   TBranch        *b_MET_T1_pt_jesSubTotalScaleUp;   //!
   TBranch        *b_MET_T1_phi_jesSubTotalScaleUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesSubTotalScaleUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesSubTotalScaleUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_MET_T1_pt_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_MET_T1_phi_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesSubTotalAbsoluteUp;   //!
   TBranch        *b_Jet_pt_jesSubTotalMCUp;   //!
   TBranch        *b_Jet_mass_jesSubTotalMCUp;   //!
   TBranch        *b_MET_T1_pt_jesSubTotalMCUp;   //!
   TBranch        *b_MET_T1_phi_jesSubTotalMCUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesSubTotalMCUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesSubTotalMCUp;   //!
   TBranch        *b_Jet_pt_jesTotalUp;   //!
   TBranch        *b_Jet_mass_jesTotalUp;   //!
   TBranch        *b_MET_T1_pt_jesTotalUp;   //!
   TBranch        *b_MET_T1_phi_jesTotalUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesTotalUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesTotalUp;   //!
   TBranch        *b_Jet_pt_jesTotalNoFlavorUp;   //!
   TBranch        *b_Jet_mass_jesTotalNoFlavorUp;   //!
   TBranch        *b_MET_T1_pt_jesTotalNoFlavorUp;   //!
   TBranch        *b_MET_T1_phi_jesTotalNoFlavorUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesTotalNoFlavorUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesTotalNoFlavorUp;   //!
   TBranch        *b_Jet_pt_jesTotalNoTimeUp;   //!
   TBranch        *b_Jet_mass_jesTotalNoTimeUp;   //!
   TBranch        *b_MET_T1_pt_jesTotalNoTimeUp;   //!
   TBranch        *b_MET_T1_phi_jesTotalNoTimeUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesTotalNoTimeUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesTotalNoTimeUp;   //!
   TBranch        *b_Jet_pt_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_Jet_mass_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_MET_T1_pt_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_MET_T1_phi_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesTotalNoFlavorNoTimeUp;   //!
   TBranch        *b_Jet_pt_jesFlavorZJetUp;   //!
   TBranch        *b_Jet_mass_jesFlavorZJetUp;   //!
   TBranch        *b_MET_T1_pt_jesFlavorZJetUp;   //!
   TBranch        *b_MET_T1_phi_jesFlavorZJetUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesFlavorZJetUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesFlavorZJetUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPhotonJetUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPhotonJetUp;   //!
   TBranch        *b_MET_T1_pt_jesFlavorPhotonJetUp;   //!
   TBranch        *b_MET_T1_phi_jesFlavorPhotonJetUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesFlavorPhotonJetUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesFlavorPhotonJetUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPureGluonUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPureGluonUp;   //!
   TBranch        *b_MET_T1_pt_jesFlavorPureGluonUp;   //!
   TBranch        *b_MET_T1_phi_jesFlavorPureGluonUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesFlavorPureGluonUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesFlavorPureGluonUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPureQuarkUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPureQuarkUp;   //!
   TBranch        *b_MET_T1_pt_jesFlavorPureQuarkUp;   //!
   TBranch        *b_MET_T1_phi_jesFlavorPureQuarkUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesFlavorPureQuarkUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesFlavorPureQuarkUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPureCharmUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPureCharmUp;   //!
   TBranch        *b_MET_T1_pt_jesFlavorPureCharmUp;   //!
   TBranch        *b_MET_T1_phi_jesFlavorPureCharmUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesFlavorPureCharmUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesFlavorPureCharmUp;   //!
   TBranch        *b_Jet_pt_jesFlavorPureBottomUp;   //!
   TBranch        *b_Jet_mass_jesFlavorPureBottomUp;   //!
   TBranch        *b_MET_T1_pt_jesFlavorPureBottomUp;   //!
   TBranch        *b_MET_T1_phi_jesFlavorPureBottomUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesFlavorPureBottomUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesFlavorPureBottomUp;   //!
   TBranch        *b_Jet_pt_jesTimeRunBUp;   //!
   TBranch        *b_Jet_mass_jesTimeRunBUp;   //!
   TBranch        *b_MET_T1_pt_jesTimeRunBUp;   //!
   TBranch        *b_MET_T1_phi_jesTimeRunBUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesTimeRunBUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesTimeRunBUp;   //!
   TBranch        *b_Jet_pt_jesTimeRunCUp;   //!
   TBranch        *b_Jet_mass_jesTimeRunCUp;   //!
   TBranch        *b_MET_T1_pt_jesTimeRunCUp;   //!
   TBranch        *b_MET_T1_phi_jesTimeRunCUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesTimeRunCUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesTimeRunCUp;   //!
   TBranch        *b_Jet_pt_jesTimeRunDEUp;   //!
   TBranch        *b_Jet_mass_jesTimeRunDEUp;   //!
   TBranch        *b_MET_T1_pt_jesTimeRunDEUp;   //!
   TBranch        *b_MET_T1_phi_jesTimeRunDEUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesTimeRunDEUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesTimeRunDEUp;   //!
   TBranch        *b_Jet_pt_jesTimeRunFUp;   //!
   TBranch        *b_Jet_mass_jesTimeRunFUp;   //!
   TBranch        *b_MET_T1_pt_jesTimeRunFUp;   //!
   TBranch        *b_MET_T1_phi_jesTimeRunFUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesTimeRunFUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesTimeRunFUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_MET_T1_pt_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_MET_T1_phi_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesCorrelationGroupMPFInSituUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_MET_T1_pt_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_MET_T1_phi_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesCorrelationGroupIntercalibrationUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_MET_T1_pt_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_MET_T1_phi_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesCorrelationGroupbJESUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_MET_T1_pt_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_MET_T1_phi_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesCorrelationGroupFlavorUp;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_MET_T1_pt_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_MET_T1_phi_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_MET_T1Smear_pt_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_MET_T1Smear_phi_jesCorrelationGroupUncorrelatedUp;   //!
   TBranch        *b_MET_T1_pt_unclustEnUp;   //!
   TBranch        *b_MET_T1_phi_unclustEnUp;   //!
   TBranch        *b_MET_T1Smear_pt_unclustEnUp;   //!
   TBranch        *b_MET_T1Smear_phi_unclustEnUp;   //!
   TBranch        *b_Jet_pt_jer0Down;   //!
   TBranch        *b_Jet_mass_jer0Down;   //!
   TBranch        *b_MET_T1_pt_jer0Down;   //!
   TBranch        *b_MET_T1_phi_jer0Down;   //!
   TBranch        *b_MET_T1Smear_pt_jer0Down;   //!
   TBranch        *b_MET_T1Smear_phi_jer0Down;   //!
   TBranch        *b_Jet_pt_jer1Down;   //!
   TBranch        *b_Jet_mass_jer1Down;   //!
   TBranch        *b_MET_T1_pt_jer1Down;   //!
   TBranch        *b_MET_T1_phi_jer1Down;   //!
   TBranch        *b_MET_T1Smear_pt_jer1Down;   //!
   TBranch        *b_MET_T1Smear_phi_jer1Down;   //!
   TBranch        *b_Jet_pt_jer2Down;   //!
   TBranch        *b_Jet_mass_jer2Down;   //!
   TBranch        *b_MET_T1_pt_jer2Down;   //!
   TBranch        *b_MET_T1_phi_jer2Down;   //!
   TBranch        *b_MET_T1Smear_pt_jer2Down;   //!
   TBranch        *b_MET_T1Smear_phi_jer2Down;   //!
   TBranch        *b_Jet_pt_jer3Down;   //!
   TBranch        *b_Jet_mass_jer3Down;   //!
   TBranch        *b_MET_T1_pt_jer3Down;   //!
   TBranch        *b_MET_T1_phi_jer3Down;   //!
   TBranch        *b_MET_T1Smear_pt_jer3Down;   //!
   TBranch        *b_MET_T1Smear_phi_jer3Down;   //!
   TBranch        *b_Jet_pt_jer4Down;   //!
   TBranch        *b_Jet_mass_jer4Down;   //!
   TBranch        *b_MET_T1_pt_jer4Down;   //!
   TBranch        *b_MET_T1_phi_jer4Down;   //!
   TBranch        *b_MET_T1Smear_pt_jer4Down;   //!
   TBranch        *b_MET_T1Smear_phi_jer4Down;   //!
   TBranch        *b_Jet_pt_jer5Down;   //!
   TBranch        *b_Jet_mass_jer5Down;   //!
   TBranch        *b_MET_T1_pt_jer5Down;   //!
   TBranch        *b_MET_T1_phi_jer5Down;   //!
   TBranch        *b_MET_T1Smear_pt_jer5Down;   //!
   TBranch        *b_MET_T1Smear_phi_jer5Down;   //!
   TBranch        *b_Jet_pt_jesAbsoluteStatDown;   //!
   TBranch        *b_Jet_mass_jesAbsoluteStatDown;   //!
   TBranch        *b_MET_T1_pt_jesAbsoluteStatDown;   //!
   TBranch        *b_MET_T1_phi_jesAbsoluteStatDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesAbsoluteStatDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesAbsoluteStatDown;   //!
   TBranch        *b_Jet_pt_jesAbsoluteScaleDown;   //!
   TBranch        *b_Jet_mass_jesAbsoluteScaleDown;   //!
   TBranch        *b_MET_T1_pt_jesAbsoluteScaleDown;   //!
   TBranch        *b_MET_T1_phi_jesAbsoluteScaleDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesAbsoluteScaleDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesAbsoluteScaleDown;   //!
   TBranch        *b_Jet_pt_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_Jet_mass_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_MET_T1_pt_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_MET_T1_phi_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesAbsoluteFlavMapDown;   //!
   TBranch        *b_Jet_pt_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_Jet_mass_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_MET_T1_pt_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_MET_T1_phi_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesAbsoluteMPFBiasDown;   //!
   TBranch        *b_Jet_pt_jesFragmentationDown;   //!
   TBranch        *b_Jet_mass_jesFragmentationDown;   //!
   TBranch        *b_MET_T1_pt_jesFragmentationDown;   //!
   TBranch        *b_MET_T1_phi_jesFragmentationDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesFragmentationDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesFragmentationDown;   //!
   TBranch        *b_Jet_pt_jesSinglePionECALDown;   //!
   TBranch        *b_Jet_mass_jesSinglePionECALDown;   //!
   TBranch        *b_MET_T1_pt_jesSinglePionECALDown;   //!
   TBranch        *b_MET_T1_phi_jesSinglePionECALDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesSinglePionECALDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesSinglePionECALDown;   //!
   TBranch        *b_Jet_pt_jesSinglePionHCALDown;   //!
   TBranch        *b_Jet_mass_jesSinglePionHCALDown;   //!
   TBranch        *b_MET_T1_pt_jesSinglePionHCALDown;   //!
   TBranch        *b_MET_T1_phi_jesSinglePionHCALDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesSinglePionHCALDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesSinglePionHCALDown;   //!
   TBranch        *b_Jet_pt_jesFlavorQCDDown;   //!
   TBranch        *b_Jet_mass_jesFlavorQCDDown;   //!
   TBranch        *b_MET_T1_pt_jesFlavorQCDDown;   //!
   TBranch        *b_MET_T1_phi_jesFlavorQCDDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesFlavorQCDDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesFlavorQCDDown;   //!
   TBranch        *b_Jet_pt_jesTimePtEtaDown;   //!
   TBranch        *b_Jet_mass_jesTimePtEtaDown;   //!
   TBranch        *b_MET_T1_pt_jesTimePtEtaDown;   //!
   TBranch        *b_MET_T1_phi_jesTimePtEtaDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesTimePtEtaDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesTimePtEtaDown;   //!
   TBranch        *b_Jet_pt_jesRelativeJEREC1Down;   //!
   TBranch        *b_Jet_mass_jesRelativeJEREC1Down;   //!
   TBranch        *b_MET_T1_pt_jesRelativeJEREC1Down;   //!
   TBranch        *b_MET_T1_phi_jesRelativeJEREC1Down;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeJEREC1Down;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeJEREC1Down;   //!
   TBranch        *b_Jet_pt_jesRelativeJEREC2Down;   //!
   TBranch        *b_Jet_mass_jesRelativeJEREC2Down;   //!
   TBranch        *b_MET_T1_pt_jesRelativeJEREC2Down;   //!
   TBranch        *b_MET_T1_phi_jesRelativeJEREC2Down;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeJEREC2Down;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeJEREC2Down;   //!
   TBranch        *b_Jet_pt_jesRelativeJERHFDown;   //!
   TBranch        *b_Jet_mass_jesRelativeJERHFDown;   //!
   TBranch        *b_MET_T1_pt_jesRelativeJERHFDown;   //!
   TBranch        *b_MET_T1_phi_jesRelativeJERHFDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeJERHFDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeJERHFDown;   //!
   TBranch        *b_Jet_pt_jesRelativePtBBDown;   //!
   TBranch        *b_Jet_mass_jesRelativePtBBDown;   //!
   TBranch        *b_MET_T1_pt_jesRelativePtBBDown;   //!
   TBranch        *b_MET_T1_phi_jesRelativePtBBDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativePtBBDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativePtBBDown;   //!
   TBranch        *b_Jet_pt_jesRelativePtEC1Down;   //!
   TBranch        *b_Jet_mass_jesRelativePtEC1Down;   //!
   TBranch        *b_MET_T1_pt_jesRelativePtEC1Down;   //!
   TBranch        *b_MET_T1_phi_jesRelativePtEC1Down;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativePtEC1Down;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativePtEC1Down;   //!
   TBranch        *b_Jet_pt_jesRelativePtEC2Down;   //!
   TBranch        *b_Jet_mass_jesRelativePtEC2Down;   //!
   TBranch        *b_MET_T1_pt_jesRelativePtEC2Down;   //!
   TBranch        *b_MET_T1_phi_jesRelativePtEC2Down;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativePtEC2Down;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativePtEC2Down;   //!
   TBranch        *b_Jet_pt_jesRelativePtHFDown;   //!
   TBranch        *b_Jet_mass_jesRelativePtHFDown;   //!
   TBranch        *b_MET_T1_pt_jesRelativePtHFDown;   //!
   TBranch        *b_MET_T1_phi_jesRelativePtHFDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativePtHFDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativePtHFDown;   //!
   TBranch        *b_Jet_pt_jesRelativeBalDown;   //!
   TBranch        *b_Jet_mass_jesRelativeBalDown;   //!
   TBranch        *b_MET_T1_pt_jesRelativeBalDown;   //!
   TBranch        *b_MET_T1_phi_jesRelativeBalDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeBalDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeBalDown;   //!
   TBranch        *b_Jet_pt_jesRelativeSampleDown;   //!
   TBranch        *b_Jet_mass_jesRelativeSampleDown;   //!
   TBranch        *b_MET_T1_pt_jesRelativeSampleDown;   //!
   TBranch        *b_MET_T1_phi_jesRelativeSampleDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeSampleDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeSampleDown;   //!
   TBranch        *b_Jet_pt_jesRelativeFSRDown;   //!
   TBranch        *b_Jet_mass_jesRelativeFSRDown;   //!
   TBranch        *b_MET_T1_pt_jesRelativeFSRDown;   //!
   TBranch        *b_MET_T1_phi_jesRelativeFSRDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeFSRDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeFSRDown;   //!
   TBranch        *b_Jet_pt_jesRelativeStatFSRDown;   //!
   TBranch        *b_Jet_mass_jesRelativeStatFSRDown;   //!
   TBranch        *b_MET_T1_pt_jesRelativeStatFSRDown;   //!
   TBranch        *b_MET_T1_phi_jesRelativeStatFSRDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeStatFSRDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeStatFSRDown;   //!
   TBranch        *b_Jet_pt_jesRelativeStatECDown;   //!
   TBranch        *b_Jet_mass_jesRelativeStatECDown;   //!
   TBranch        *b_MET_T1_pt_jesRelativeStatECDown;   //!
   TBranch        *b_MET_T1_phi_jesRelativeStatECDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeStatECDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeStatECDown;   //!
   TBranch        *b_Jet_pt_jesRelativeStatHFDown;   //!
   TBranch        *b_Jet_mass_jesRelativeStatHFDown;   //!
   TBranch        *b_MET_T1_pt_jesRelativeStatHFDown;   //!
   TBranch        *b_MET_T1_phi_jesRelativeStatHFDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesRelativeStatHFDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesRelativeStatHFDown;   //!
   TBranch        *b_Jet_pt_jesPileUpDataMCDown;   //!
   TBranch        *b_Jet_mass_jesPileUpDataMCDown;   //!
   TBranch        *b_MET_T1_pt_jesPileUpDataMCDown;   //!
   TBranch        *b_MET_T1_phi_jesPileUpDataMCDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpDataMCDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpDataMCDown;   //!
   TBranch        *b_Jet_pt_jesPileUpPtRefDown;   //!
   TBranch        *b_Jet_mass_jesPileUpPtRefDown;   //!
   TBranch        *b_MET_T1_pt_jesPileUpPtRefDown;   //!
   TBranch        *b_MET_T1_phi_jesPileUpPtRefDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpPtRefDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpPtRefDown;   //!
   TBranch        *b_Jet_pt_jesPileUpPtBBDown;   //!
   TBranch        *b_Jet_mass_jesPileUpPtBBDown;   //!
   TBranch        *b_MET_T1_pt_jesPileUpPtBBDown;   //!
   TBranch        *b_MET_T1_phi_jesPileUpPtBBDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpPtBBDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpPtBBDown;   //!
   TBranch        *b_Jet_pt_jesPileUpPtEC1Down;   //!
   TBranch        *b_Jet_mass_jesPileUpPtEC1Down;   //!
   TBranch        *b_MET_T1_pt_jesPileUpPtEC1Down;   //!
   TBranch        *b_MET_T1_phi_jesPileUpPtEC1Down;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpPtEC1Down;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpPtEC1Down;   //!
   TBranch        *b_Jet_pt_jesPileUpPtEC2Down;   //!
   TBranch        *b_Jet_mass_jesPileUpPtEC2Down;   //!
   TBranch        *b_MET_T1_pt_jesPileUpPtEC2Down;   //!
   TBranch        *b_MET_T1_phi_jesPileUpPtEC2Down;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpPtEC2Down;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpPtEC2Down;   //!
   TBranch        *b_Jet_pt_jesPileUpPtHFDown;   //!
   TBranch        *b_Jet_mass_jesPileUpPtHFDown;   //!
   TBranch        *b_MET_T1_pt_jesPileUpPtHFDown;   //!
   TBranch        *b_MET_T1_phi_jesPileUpPtHFDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpPtHFDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpPtHFDown;   //!
   TBranch        *b_Jet_pt_jesPileUpMuZeroDown;   //!
   TBranch        *b_Jet_mass_jesPileUpMuZeroDown;   //!
   TBranch        *b_MET_T1_pt_jesPileUpMuZeroDown;   //!
   TBranch        *b_MET_T1_phi_jesPileUpMuZeroDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpMuZeroDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpMuZeroDown;   //!
   TBranch        *b_Jet_pt_jesPileUpEnvelopeDown;   //!
   TBranch        *b_Jet_mass_jesPileUpEnvelopeDown;   //!
   TBranch        *b_MET_T1_pt_jesPileUpEnvelopeDown;   //!
   TBranch        *b_MET_T1_phi_jesPileUpEnvelopeDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesPileUpEnvelopeDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesPileUpEnvelopeDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalPileUpDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalPileUpDown;   //!
   TBranch        *b_MET_T1_pt_jesSubTotalPileUpDown;   //!
   TBranch        *b_MET_T1_phi_jesSubTotalPileUpDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesSubTotalPileUpDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesSubTotalPileUpDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalRelativeDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalRelativeDown;   //!
   TBranch        *b_MET_T1_pt_jesSubTotalRelativeDown;   //!
   TBranch        *b_MET_T1_phi_jesSubTotalRelativeDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesSubTotalRelativeDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesSubTotalRelativeDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalPtDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalPtDown;   //!
   TBranch        *b_MET_T1_pt_jesSubTotalPtDown;   //!
   TBranch        *b_MET_T1_phi_jesSubTotalPtDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesSubTotalPtDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesSubTotalPtDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalScaleDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalScaleDown;   //!
   TBranch        *b_MET_T1_pt_jesSubTotalScaleDown;   //!
   TBranch        *b_MET_T1_phi_jesSubTotalScaleDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesSubTotalScaleDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesSubTotalScaleDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_MET_T1_pt_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_MET_T1_phi_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesSubTotalAbsoluteDown;   //!
   TBranch        *b_Jet_pt_jesSubTotalMCDown;   //!
   TBranch        *b_Jet_mass_jesSubTotalMCDown;   //!
   TBranch        *b_MET_T1_pt_jesSubTotalMCDown;   //!
   TBranch        *b_MET_T1_phi_jesSubTotalMCDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesSubTotalMCDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesSubTotalMCDown;   //!
   TBranch        *b_Jet_pt_jesTotalDown;   //!
   TBranch        *b_Jet_mass_jesTotalDown;   //!
   TBranch        *b_MET_T1_pt_jesTotalDown;   //!
   TBranch        *b_MET_T1_phi_jesTotalDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesTotalDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesTotalDown;   //!
   TBranch        *b_Jet_pt_jesTotalNoFlavorDown;   //!
   TBranch        *b_Jet_mass_jesTotalNoFlavorDown;   //!
   TBranch        *b_MET_T1_pt_jesTotalNoFlavorDown;   //!
   TBranch        *b_MET_T1_phi_jesTotalNoFlavorDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesTotalNoFlavorDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesTotalNoFlavorDown;   //!
   TBranch        *b_Jet_pt_jesTotalNoTimeDown;   //!
   TBranch        *b_Jet_mass_jesTotalNoTimeDown;   //!
   TBranch        *b_MET_T1_pt_jesTotalNoTimeDown;   //!
   TBranch        *b_MET_T1_phi_jesTotalNoTimeDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesTotalNoTimeDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesTotalNoTimeDown;   //!
   TBranch        *b_Jet_pt_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_Jet_mass_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_MET_T1_pt_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_MET_T1_phi_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesTotalNoFlavorNoTimeDown;   //!
   TBranch        *b_Jet_pt_jesFlavorZJetDown;   //!
   TBranch        *b_Jet_mass_jesFlavorZJetDown;   //!
   TBranch        *b_MET_T1_pt_jesFlavorZJetDown;   //!
   TBranch        *b_MET_T1_phi_jesFlavorZJetDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesFlavorZJetDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesFlavorZJetDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPhotonJetDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPhotonJetDown;   //!
   TBranch        *b_MET_T1_pt_jesFlavorPhotonJetDown;   //!
   TBranch        *b_MET_T1_phi_jesFlavorPhotonJetDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesFlavorPhotonJetDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesFlavorPhotonJetDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPureGluonDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPureGluonDown;   //!
   TBranch        *b_MET_T1_pt_jesFlavorPureGluonDown;   //!
   TBranch        *b_MET_T1_phi_jesFlavorPureGluonDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesFlavorPureGluonDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesFlavorPureGluonDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPureQuarkDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPureQuarkDown;   //!
   TBranch        *b_MET_T1_pt_jesFlavorPureQuarkDown;   //!
   TBranch        *b_MET_T1_phi_jesFlavorPureQuarkDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesFlavorPureQuarkDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesFlavorPureQuarkDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPureCharmDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPureCharmDown;   //!
   TBranch        *b_MET_T1_pt_jesFlavorPureCharmDown;   //!
   TBranch        *b_MET_T1_phi_jesFlavorPureCharmDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesFlavorPureCharmDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesFlavorPureCharmDown;   //!
   TBranch        *b_Jet_pt_jesFlavorPureBottomDown;   //!
   TBranch        *b_Jet_mass_jesFlavorPureBottomDown;   //!
   TBranch        *b_MET_T1_pt_jesFlavorPureBottomDown;   //!
   TBranch        *b_MET_T1_phi_jesFlavorPureBottomDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesFlavorPureBottomDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesFlavorPureBottomDown;   //!
   TBranch        *b_Jet_pt_jesTimeRunBDown;   //!
   TBranch        *b_Jet_mass_jesTimeRunBDown;   //!
   TBranch        *b_MET_T1_pt_jesTimeRunBDown;   //!
   TBranch        *b_MET_T1_phi_jesTimeRunBDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesTimeRunBDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesTimeRunBDown;   //!
   TBranch        *b_Jet_pt_jesTimeRunCDown;   //!
   TBranch        *b_Jet_mass_jesTimeRunCDown;   //!
   TBranch        *b_MET_T1_pt_jesTimeRunCDown;   //!
   TBranch        *b_MET_T1_phi_jesTimeRunCDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesTimeRunCDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesTimeRunCDown;   //!
   TBranch        *b_Jet_pt_jesTimeRunDEDown;   //!
   TBranch        *b_Jet_mass_jesTimeRunDEDown;   //!
   TBranch        *b_MET_T1_pt_jesTimeRunDEDown;   //!
   TBranch        *b_MET_T1_phi_jesTimeRunDEDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesTimeRunDEDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesTimeRunDEDown;   //!
   TBranch        *b_Jet_pt_jesTimeRunFDown;   //!
   TBranch        *b_Jet_mass_jesTimeRunFDown;   //!
   TBranch        *b_MET_T1_pt_jesTimeRunFDown;   //!
   TBranch        *b_MET_T1_phi_jesTimeRunFDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesTimeRunFDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesTimeRunFDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_MET_T1_pt_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_MET_T1_phi_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesCorrelationGroupMPFInSituDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_MET_T1_pt_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_MET_T1_phi_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesCorrelationGroupIntercalibrationDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_MET_T1_pt_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_MET_T1_phi_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesCorrelationGroupbJESDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_MET_T1_pt_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_MET_T1_phi_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesCorrelationGroupFlavorDown;   //!
   TBranch        *b_Jet_pt_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_Jet_mass_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_MET_T1_pt_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_MET_T1_phi_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_MET_T1Smear_pt_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_MET_T1Smear_phi_jesCorrelationGroupUncorrelatedDown;   //!
   TBranch        *b_MET_T1_pt_unclustEnDown;   //!
   TBranch        *b_MET_T1_phi_unclustEnDown;   //!
   TBranch        *b_MET_T1Smear_pt_unclustEnDown;   //!
   TBranch        *b_MET_T1Smear_phi_unclustEnDown;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_puWeightUp;   //!
   TBranch        *b_puWeightDown;   //!
   TBranch        *b_PrefireWeight;   //!
   TBranch        *b_PrefireWeight_Up;   //!
   TBranch        *b_PrefireWeight_Down;   //!
   TBranch        *b_Jet_btagSF_deepcsv_M_down;   //!
   TBranch        *b_Jet_btagSF_deepcsv_M;   //!
   TBranch        *b_Jet_btagSF_deepcsv_M_up;   //!
   TBranch        *b_Jet_btagSF_deepcsv_L_down;   //!
   TBranch        *b_Jet_btagSF_deepcsv_L;   //!
   TBranch        *b_Jet_btagSF_deepcsv_L_up;   //!
   TBranch        *b_Jet_btagSF_deepcsv_T_down;   //!
   TBranch        *b_Jet_btagSF_deepcsv_T;   //!
   TBranch        *b_Jet_btagSF_deepcsv_T_up;   //!
   TBranch        *b_Jet_btagSF_deepjet_M_down;   //!
   TBranch        *b_Jet_btagSF_deepjet_M;   //!
   TBranch        *b_Jet_btagSF_deepjet_M_up;   //!
   TBranch        *b_Jet_btagSF_deepjet_L_down;   //!
   TBranch        *b_Jet_btagSF_deepjet_L;   //!
   TBranch        *b_Jet_btagSF_deepjet_L_up;   //!
   TBranch        *b_Jet_btagSF_deepjet_T_down;   //!
   TBranch        *b_Jet_btagSF_deepjet_T;   //!
   TBranch        *b_Jet_btagSF_deepjet_T_up;   //!

   MyAnalysis(TTree *tree=0);
   virtual ~MyAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString fname, TString data, TString dataset ,string year, TString run, float xs, float lumi, float Nevent, int iseft, int nRuns);
//   virtual void     Loop(TString, TString, TString, TString, TString, TString, float,float,float);
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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_TTJets/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_UL17_TTJets/211106_221617/0000/tree_77.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_TTJets/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_UL17_TTJets/211106_221617/0000/tree_77.root");
      }
      f->GetObject("Events",tree);

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

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("HTXS_Higgs_pt", &HTXS_Higgs_pt, &b_HTXS_Higgs_pt);
   fChain->SetBranchAddress("HTXS_Higgs_y", &HTXS_Higgs_y, &b_HTXS_Higgs_y);
   fChain->SetBranchAddress("HTXS_stage1_1_cat_pTjet25GeV", &HTXS_stage1_1_cat_pTjet25GeV, &b_HTXS_stage1_1_cat_pTjet25GeV);
   fChain->SetBranchAddress("HTXS_stage1_1_cat_pTjet30GeV", &HTXS_stage1_1_cat_pTjet30GeV, &b_HTXS_stage1_1_cat_pTjet30GeV);
   fChain->SetBranchAddress("HTXS_stage1_1_fine_cat_pTjet25GeV", &HTXS_stage1_1_fine_cat_pTjet25GeV, &b_HTXS_stage1_1_fine_cat_pTjet25GeV);
   fChain->SetBranchAddress("HTXS_stage1_1_fine_cat_pTjet30GeV", &HTXS_stage1_1_fine_cat_pTjet30GeV, &b_HTXS_stage1_1_fine_cat_pTjet30GeV);
   fChain->SetBranchAddress("HTXS_stage1_2_cat_pTjet25GeV", &HTXS_stage1_2_cat_pTjet25GeV, &b_HTXS_stage1_2_cat_pTjet25GeV);
   fChain->SetBranchAddress("HTXS_stage1_2_cat_pTjet30GeV", &HTXS_stage1_2_cat_pTjet30GeV, &b_HTXS_stage1_2_cat_pTjet30GeV);
   fChain->SetBranchAddress("HTXS_stage1_2_fine_cat_pTjet25GeV", &HTXS_stage1_2_fine_cat_pTjet25GeV, &b_HTXS_stage1_2_fine_cat_pTjet25GeV);
   fChain->SetBranchAddress("HTXS_stage1_2_fine_cat_pTjet30GeV", &HTXS_stage1_2_fine_cat_pTjet30GeV, &b_HTXS_stage1_2_fine_cat_pTjet30GeV);
   fChain->SetBranchAddress("HTXS_stage_0", &HTXS_stage_0, &b_HTXS_stage_0);
   fChain->SetBranchAddress("HTXS_stage_1_pTjet25", &HTXS_stage_1_pTjet25, &b_HTXS_stage_1_pTjet25);
   fChain->SetBranchAddress("HTXS_stage_1_pTjet30", &HTXS_stage_1_pTjet30, &b_HTXS_stage_1_pTjet30);
   fChain->SetBranchAddress("HTXS_njets25", &HTXS_njets25, &b_HTXS_njets25);
   fChain->SetBranchAddress("HTXS_njets30", &HTXS_njets30, &b_HTXS_njets30);
   fChain->SetBranchAddress("nboostedTau", &nboostedTau, &b_nboostedTau);
   fChain->SetBranchAddress("boostedTau_chargedIso", boostedTau_chargedIso, &b_boostedTau_chargedIso);
   fChain->SetBranchAddress("boostedTau_eta", boostedTau_eta, &b_boostedTau_eta);
   fChain->SetBranchAddress("boostedTau_leadTkDeltaEta", boostedTau_leadTkDeltaEta, &b_boostedTau_leadTkDeltaEta);
   fChain->SetBranchAddress("boostedTau_leadTkDeltaPhi", boostedTau_leadTkDeltaPhi, &b_boostedTau_leadTkDeltaPhi);
   fChain->SetBranchAddress("boostedTau_leadTkPtOverTauPt", boostedTau_leadTkPtOverTauPt, &b_boostedTau_leadTkPtOverTauPt);
   fChain->SetBranchAddress("boostedTau_mass", boostedTau_mass, &b_boostedTau_mass);
   fChain->SetBranchAddress("boostedTau_neutralIso", boostedTau_neutralIso, &b_boostedTau_neutralIso);
   fChain->SetBranchAddress("boostedTau_phi", boostedTau_phi, &b_boostedTau_phi);
   fChain->SetBranchAddress("boostedTau_photonsOutsideSignalCone", boostedTau_photonsOutsideSignalCone, &b_boostedTau_photonsOutsideSignalCone);
   fChain->SetBranchAddress("boostedTau_pt", boostedTau_pt, &b_boostedTau_pt);
   fChain->SetBranchAddress("boostedTau_puCorr", boostedTau_puCorr, &b_boostedTau_puCorr);
   fChain->SetBranchAddress("boostedTau_rawAntiEle2018", boostedTau_rawAntiEle2018, &b_boostedTau_rawAntiEle2018);
   fChain->SetBranchAddress("boostedTau_rawIso", boostedTau_rawIso, &b_boostedTau_rawIso);
   fChain->SetBranchAddress("boostedTau_rawIsodR03", boostedTau_rawIsodR03, &b_boostedTau_rawIsodR03);
   fChain->SetBranchAddress("boostedTau_rawMVAnewDM2017v2", boostedTau_rawMVAnewDM2017v2, &b_boostedTau_rawMVAnewDM2017v2);
   fChain->SetBranchAddress("boostedTau_rawMVAoldDM2017v2", boostedTau_rawMVAoldDM2017v2, &b_boostedTau_rawMVAoldDM2017v2);
   fChain->SetBranchAddress("boostedTau_rawMVAoldDMdR032017v2", boostedTau_rawMVAoldDMdR032017v2, &b_boostedTau_rawMVAoldDMdR032017v2);
   fChain->SetBranchAddress("boostedTau_charge", boostedTau_charge, &b_boostedTau_charge);
   fChain->SetBranchAddress("boostedTau_decayMode", boostedTau_decayMode, &b_boostedTau_decayMode);
   fChain->SetBranchAddress("boostedTau_jetIdx", boostedTau_jetIdx, &b_boostedTau_jetIdx);
   fChain->SetBranchAddress("boostedTau_rawAntiEleCat2018", boostedTau_rawAntiEleCat2018, &b_boostedTau_rawAntiEleCat2018);
   fChain->SetBranchAddress("boostedTau_idAntiEle2018", boostedTau_idAntiEle2018, &b_boostedTau_idAntiEle2018);
   fChain->SetBranchAddress("boostedTau_idAntiMu", boostedTau_idAntiMu, &b_boostedTau_idAntiMu);
   fChain->SetBranchAddress("boostedTau_idMVAnewDM2017v2", boostedTau_idMVAnewDM2017v2, &b_boostedTau_idMVAnewDM2017v2);
   fChain->SetBranchAddress("boostedTau_idMVAoldDM2017v2", boostedTau_idMVAoldDM2017v2, &b_boostedTau_idMVAoldDM2017v2);
   fChain->SetBranchAddress("boostedTau_idMVAoldDMdR032017v2", boostedTau_idMVAoldDMdR032017v2, &b_boostedTau_idMVAoldDMdR032017v2);
   fChain->SetBranchAddress("btagWeight_CSVV2", &btagWeight_CSVV2, &b_btagWeight_CSVV2);
   fChain->SetBranchAddress("btagWeight_DeepCSVB", &btagWeight_DeepCSVB, &b_btagWeight_DeepCSVB);
   fChain->SetBranchAddress("CaloMET_phi", &CaloMET_phi, &b_CaloMET_phi);
   fChain->SetBranchAddress("CaloMET_pt", &CaloMET_pt, &b_CaloMET_pt);
   fChain->SetBranchAddress("CaloMET_sumEt", &CaloMET_sumEt, &b_CaloMET_sumEt);
   fChain->SetBranchAddress("ChsMET_phi", &ChsMET_phi, &b_ChsMET_phi);
   fChain->SetBranchAddress("ChsMET_pt", &ChsMET_pt, &b_ChsMET_pt);
   fChain->SetBranchAddress("ChsMET_sumEt", &ChsMET_sumEt, &b_ChsMET_sumEt);
   fChain->SetBranchAddress("nCorrT1METJet", &nCorrT1METJet, &b_nCorrT1METJet);
   fChain->SetBranchAddress("CorrT1METJet_area", CorrT1METJet_area, &b_CorrT1METJet_area);
   fChain->SetBranchAddress("CorrT1METJet_eta", CorrT1METJet_eta, &b_CorrT1METJet_eta);
   fChain->SetBranchAddress("CorrT1METJet_muonSubtrFactor", CorrT1METJet_muonSubtrFactor, &b_CorrT1METJet_muonSubtrFactor);
   fChain->SetBranchAddress("CorrT1METJet_phi", CorrT1METJet_phi, &b_CorrT1METJet_phi);
   fChain->SetBranchAddress("CorrT1METJet_rawPt", CorrT1METJet_rawPt, &b_CorrT1METJet_rawPt);
   fChain->SetBranchAddress("DeepMETResolutionTune_phi", &DeepMETResolutionTune_phi, &b_DeepMETResolutionTune_phi);
   fChain->SetBranchAddress("DeepMETResolutionTune_pt", &DeepMETResolutionTune_pt, &b_DeepMETResolutionTune_pt);
   fChain->SetBranchAddress("DeepMETResponseTune_phi", &DeepMETResponseTune_phi, &b_DeepMETResponseTune_phi);
   fChain->SetBranchAddress("DeepMETResponseTune_pt", &DeepMETResponseTune_pt, &b_DeepMETResponseTune_pt);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("Electron_dEscaleDown", Electron_dEscaleDown, &b_Electron_dEscaleDown);
   fChain->SetBranchAddress("Electron_dEscaleUp", Electron_dEscaleUp, &b_Electron_dEscaleUp);
   fChain->SetBranchAddress("Electron_dEsigmaDown", Electron_dEsigmaDown, &b_Electron_dEsigmaDown);
   fChain->SetBranchAddress("Electron_dEsigmaUp", Electron_dEsigmaUp, &b_Electron_dEsigmaUp);
   fChain->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC, &b_Electron_deltaEtaSC);
   fChain->SetBranchAddress("Electron_dr03EcalRecHitSumEt", Electron_dr03EcalRecHitSumEt, &b_Electron_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", Electron_dr03HcalDepth1TowerSumEt, &b_Electron_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Electron_dr03TkSumPt", Electron_dr03TkSumPt, &b_Electron_dr03TkSumPt);
   fChain->SetBranchAddress("Electron_dr03TkSumPtHEEP", Electron_dr03TkSumPtHEEP, &b_Electron_dr03TkSumPtHEEP);
   fChain->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dxyErr", Electron_dxyErr, &b_Electron_dxyErr);
   fChain->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron_dzErr", Electron_dzErr, &b_Electron_dzErr);
   fChain->SetBranchAddress("Electron_eCorr", Electron_eCorr, &b_Electron_eCorr);
   fChain->SetBranchAddress("Electron_eInvMinusPInv", Electron_eInvMinusPInv, &b_Electron_eInvMinusPInv);
   fChain->SetBranchAddress("Electron_energyErr", Electron_energyErr, &b_Electron_energyErr);
   fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_hoe", Electron_hoe, &b_Electron_hoe);
   fChain->SetBranchAddress("Electron_ip3d", Electron_ip3d, &b_Electron_ip3d);
   fChain->SetBranchAddress("Electron_jetPtRelv2", Electron_jetPtRelv2, &b_Electron_jetPtRelv2);
   fChain->SetBranchAddress("Electron_jetRelIso", Electron_jetRelIso, &b_Electron_jetRelIso);
   fChain->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
   fChain->SetBranchAddress("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all, &b_Electron_miniPFRelIso_all);
   fChain->SetBranchAddress("Electron_miniPFRelIso_chg", Electron_miniPFRelIso_chg, &b_Electron_miniPFRelIso_chg);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso", Electron_mvaFall17V2Iso, &b_Electron_mvaFall17V2Iso);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso", Electron_mvaFall17V2noIso, &b_Electron_mvaFall17V2noIso);
   fChain->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
   fChain->SetBranchAddress("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg, &b_Electron_pfRelIso03_chg);
   fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_r9", Electron_r9, &b_Electron_r9);
   fChain->SetBranchAddress("Electron_scEtOverPt", Electron_scEtOverPt, &b_Electron_scEtOverPt);
   fChain->SetBranchAddress("Electron_sieie", Electron_sieie, &b_Electron_sieie);
   fChain->SetBranchAddress("Electron_sip3d", Electron_sip3d, &b_Electron_sip3d);
   fChain->SetBranchAddress("Electron_mvaTTH", Electron_mvaTTH, &b_Electron_mvaTTH);
   fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_cutBased", Electron_cutBased, &b_Electron_cutBased);
   fChain->SetBranchAddress("Electron_jetIdx", Electron_jetIdx, &b_Electron_jetIdx);
   fChain->SetBranchAddress("Electron_pdgId", Electron_pdgId, &b_Electron_pdgId);
   fChain->SetBranchAddress("Electron_photonIdx", Electron_photonIdx, &b_Electron_photonIdx);
   fChain->SetBranchAddress("Electron_tightCharge", Electron_tightCharge, &b_Electron_tightCharge);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmap", Electron_vidNestedWPBitmap, &b_Electron_vidNestedWPBitmap);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmapHEEP", Electron_vidNestedWPBitmapHEEP, &b_Electron_vidNestedWPBitmapHEEP);
   fChain->SetBranchAddress("Electron_convVeto", Electron_convVeto, &b_Electron_convVeto);
   fChain->SetBranchAddress("Electron_cutBased_HEEP", Electron_cutBased_HEEP, &b_Electron_cutBased_HEEP);
   fChain->SetBranchAddress("Electron_isPFcand", Electron_isPFcand, &b_Electron_isPFcand);
   fChain->SetBranchAddress("Electron_jetNDauCharged", Electron_jetNDauCharged, &b_Electron_jetNDauCharged);
   fChain->SetBranchAddress("Electron_lostHits", Electron_lostHits, &b_Electron_lostHits);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WP80", Electron_mvaFall17V2Iso_WP80, &b_Electron_mvaFall17V2Iso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WP90", Electron_mvaFall17V2Iso_WP90, &b_Electron_mvaFall17V2Iso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WPL", Electron_mvaFall17V2Iso_WPL, &b_Electron_mvaFall17V2Iso_WPL);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WP80", Electron_mvaFall17V2noIso_WP80, &b_Electron_mvaFall17V2noIso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WP90", Electron_mvaFall17V2noIso_WP90, &b_Electron_mvaFall17V2noIso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WPL", Electron_mvaFall17V2noIso_WPL, &b_Electron_mvaFall17V2noIso_WPL);
   fChain->SetBranchAddress("Electron_seedGain", Electron_seedGain, &b_Electron_seedGain);
   fChain->SetBranchAddress("nFatJet", &nFatJet, &b_nFatJet);
   fChain->SetBranchAddress("FatJet_area", FatJet_area, &b_FatJet_area);
   fChain->SetBranchAddress("FatJet_btagCSVV2", FatJet_btagCSVV2, &b_FatJet_btagCSVV2);
   fChain->SetBranchAddress("FatJet_btagDDBvLV2", FatJet_btagDDBvLV2, &b_FatJet_btagDDBvLV2);
   fChain->SetBranchAddress("FatJet_btagDDCvBV2", FatJet_btagDDCvBV2, &b_FatJet_btagDDCvBV2);
   fChain->SetBranchAddress("FatJet_btagDDCvLV2", FatJet_btagDDCvLV2, &b_FatJet_btagDDCvLV2);
   fChain->SetBranchAddress("FatJet_btagDeepB", FatJet_btagDeepB, &b_FatJet_btagDeepB);
   fChain->SetBranchAddress("FatJet_btagHbb", FatJet_btagHbb, &b_FatJet_btagHbb);
   fChain->SetBranchAddress("FatJet_deepTagMD_H4qvsQCD", FatJet_deepTagMD_H4qvsQCD, &b_FatJet_deepTagMD_H4qvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_HbbvsQCD", FatJet_deepTagMD_HbbvsQCD, &b_FatJet_deepTagMD_HbbvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_TvsQCD", FatJet_deepTagMD_TvsQCD, &b_FatJet_deepTagMD_TvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_WvsQCD", FatJet_deepTagMD_WvsQCD, &b_FatJet_deepTagMD_WvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_ZHbbvsQCD", FatJet_deepTagMD_ZHbbvsQCD, &b_FatJet_deepTagMD_ZHbbvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_ZHccvsQCD", FatJet_deepTagMD_ZHccvsQCD, &b_FatJet_deepTagMD_ZHccvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_ZbbvsQCD", FatJet_deepTagMD_ZbbvsQCD, &b_FatJet_deepTagMD_ZbbvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_ZvsQCD", FatJet_deepTagMD_ZvsQCD, &b_FatJet_deepTagMD_ZvsQCD);
   fChain->SetBranchAddress("FatJet_deepTagMD_bbvsLight", FatJet_deepTagMD_bbvsLight, &b_FatJet_deepTagMD_bbvsLight);
   fChain->SetBranchAddress("FatJet_deepTagMD_ccvsLight", FatJet_deepTagMD_ccvsLight, &b_FatJet_deepTagMD_ccvsLight);
   fChain->SetBranchAddress("FatJet_deepTag_H", FatJet_deepTag_H, &b_FatJet_deepTag_H);
   fChain->SetBranchAddress("FatJet_deepTag_QCD", FatJet_deepTag_QCD, &b_FatJet_deepTag_QCD);
   fChain->SetBranchAddress("FatJet_deepTag_QCDothers", FatJet_deepTag_QCDothers, &b_FatJet_deepTag_QCDothers);
   fChain->SetBranchAddress("FatJet_deepTag_TvsQCD", FatJet_deepTag_TvsQCD, &b_FatJet_deepTag_TvsQCD);
   fChain->SetBranchAddress("FatJet_deepTag_WvsQCD", FatJet_deepTag_WvsQCD, &b_FatJet_deepTag_WvsQCD);
   fChain->SetBranchAddress("FatJet_deepTag_ZvsQCD", FatJet_deepTag_ZvsQCD, &b_FatJet_deepTag_ZvsQCD);
   fChain->SetBranchAddress("FatJet_eta", FatJet_eta, &b_FatJet_eta);
   fChain->SetBranchAddress("FatJet_mass", FatJet_mass, &b_FatJet_mass);
   fChain->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop, &b_FatJet_msoftdrop);
   fChain->SetBranchAddress("FatJet_n2b1", FatJet_n2b1, &b_FatJet_n2b1);
   fChain->SetBranchAddress("FatJet_n3b1", FatJet_n3b1, &b_FatJet_n3b1);
   fChain->SetBranchAddress("FatJet_particleNetMD_QCD", FatJet_particleNetMD_QCD, &b_FatJet_particleNetMD_QCD);
   fChain->SetBranchAddress("FatJet_particleNetMD_Xbb", FatJet_particleNetMD_Xbb, &b_FatJet_particleNetMD_Xbb);
   fChain->SetBranchAddress("FatJet_particleNetMD_Xcc", FatJet_particleNetMD_Xcc, &b_FatJet_particleNetMD_Xcc);
   fChain->SetBranchAddress("FatJet_particleNetMD_Xqq", FatJet_particleNetMD_Xqq, &b_FatJet_particleNetMD_Xqq);
   fChain->SetBranchAddress("FatJet_particleNet_H4qvsQCD", FatJet_particleNet_H4qvsQCD, &b_FatJet_particleNet_H4qvsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_HbbvsQCD", FatJet_particleNet_HbbvsQCD, &b_FatJet_particleNet_HbbvsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_HccvsQCD", FatJet_particleNet_HccvsQCD, &b_FatJet_particleNet_HccvsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_QCD", FatJet_particleNet_QCD, &b_FatJet_particleNet_QCD);
   fChain->SetBranchAddress("FatJet_particleNet_TvsQCD", FatJet_particleNet_TvsQCD, &b_FatJet_particleNet_TvsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_WvsQCD", FatJet_particleNet_WvsQCD, &b_FatJet_particleNet_WvsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_ZvsQCD", FatJet_particleNet_ZvsQCD, &b_FatJet_particleNet_ZvsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_mass", FatJet_particleNet_mass, &b_FatJet_particleNet_mass);
   fChain->SetBranchAddress("FatJet_phi", FatJet_phi, &b_FatJet_phi);
   fChain->SetBranchAddress("FatJet_pt", FatJet_pt, &b_FatJet_pt);
   fChain->SetBranchAddress("FatJet_rawFactor", FatJet_rawFactor, &b_FatJet_rawFactor);
   fChain->SetBranchAddress("FatJet_tau1", FatJet_tau1, &b_FatJet_tau1);
   fChain->SetBranchAddress("FatJet_tau2", FatJet_tau2, &b_FatJet_tau2);
   fChain->SetBranchAddress("FatJet_tau3", FatJet_tau3, &b_FatJet_tau3);
   fChain->SetBranchAddress("FatJet_tau4", FatJet_tau4, &b_FatJet_tau4);
   fChain->SetBranchAddress("FatJet_lsf3", FatJet_lsf3, &b_FatJet_lsf3);
   fChain->SetBranchAddress("FatJet_jetId", FatJet_jetId, &b_FatJet_jetId);
   fChain->SetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1, &b_FatJet_subJetIdx1);
   fChain->SetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2, &b_FatJet_subJetIdx2);
   fChain->SetBranchAddress("FatJet_electronIdx3SJ", FatJet_electronIdx3SJ, &b_FatJet_electronIdx3SJ);
   fChain->SetBranchAddress("FatJet_muonIdx3SJ", FatJet_muonIdx3SJ, &b_FatJet_muonIdx3SJ);
   fChain->SetBranchAddress("FatJet_nConstituents", FatJet_nConstituents, &b_FatJet_nConstituents);
   fChain->SetBranchAddress("nFsrPhoton", &nFsrPhoton, &b_nFsrPhoton);
   fChain->SetBranchAddress("FsrPhoton_dROverEt2", FsrPhoton_dROverEt2, &b_FsrPhoton_dROverEt2);
   fChain->SetBranchAddress("FsrPhoton_eta", FsrPhoton_eta, &b_FsrPhoton_eta);
   fChain->SetBranchAddress("FsrPhoton_phi", FsrPhoton_phi, &b_FsrPhoton_phi);
   fChain->SetBranchAddress("FsrPhoton_pt", FsrPhoton_pt, &b_FsrPhoton_pt);
   fChain->SetBranchAddress("FsrPhoton_relIso03", FsrPhoton_relIso03, &b_FsrPhoton_relIso03);
   fChain->SetBranchAddress("FsrPhoton_muonIdx", FsrPhoton_muonIdx, &b_FsrPhoton_muonIdx);
   fChain->SetBranchAddress("nGenJetAK8", &nGenJetAK8, &b_nGenJetAK8);
   fChain->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta, &b_GenJetAK8_eta);
   fChain->SetBranchAddress("GenJetAK8_mass", GenJetAK8_mass, &b_GenJetAK8_mass);
   fChain->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi, &b_GenJetAK8_phi);
   fChain->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt, &b_GenJetAK8_pt);
   fChain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
   fChain->SetBranchAddress("GenJet_eta", GenJet_eta, &b_GenJet_eta);
   fChain->SetBranchAddress("GenJet_mass", GenJet_mass, &b_GenJet_mass);
   fChain->SetBranchAddress("GenJet_phi", GenJet_phi, &b_GenJet_phi);
   fChain->SetBranchAddress("GenJet_pt", GenJet_pt, &b_GenJet_pt);
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
   fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
   fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
   fChain->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags, &b_GenPart_statusFlags);
   fChain->SetBranchAddress("nSubGenJetAK8", &nSubGenJetAK8, &b_nSubGenJetAK8);
   fChain->SetBranchAddress("SubGenJetAK8_eta", SubGenJetAK8_eta, &b_SubGenJetAK8_eta);
   fChain->SetBranchAddress("SubGenJetAK8_mass", SubGenJetAK8_mass, &b_SubGenJetAK8_mass);
   fChain->SetBranchAddress("SubGenJetAK8_phi", SubGenJetAK8_phi, &b_SubGenJetAK8_phi);
   fChain->SetBranchAddress("SubGenJetAK8_pt", SubGenJetAK8_pt, &b_SubGenJetAK8_pt);
   fChain->SetBranchAddress("Generator_binvar", &Generator_binvar, &b_Generator_binvar);
   fChain->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF, &b_Generator_scalePDF);
   fChain->SetBranchAddress("Generator_weight", &Generator_weight, &b_Generator_weight);
   fChain->SetBranchAddress("Generator_x1", &Generator_x1, &b_Generator_x1);
   fChain->SetBranchAddress("Generator_x2", &Generator_x2, &b_Generator_x2);
   fChain->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1, &b_Generator_xpdf1);
   fChain->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2, &b_Generator_xpdf2);
   fChain->SetBranchAddress("Generator_id1", &Generator_id1, &b_Generator_id1);
   fChain->SetBranchAddress("Generator_id2", &Generator_id2, &b_Generator_id2);
   fChain->SetBranchAddress("GenVtx_x", &GenVtx_x, &b_GenVtx_x);
   fChain->SetBranchAddress("GenVtx_y", &GenVtx_y, &b_GenVtx_y);
   fChain->SetBranchAddress("GenVtx_z", &GenVtx_z, &b_GenVtx_z);
   fChain->SetBranchAddress("nGenVisTau", &nGenVisTau, &b_nGenVisTau);
   fChain->SetBranchAddress("GenVisTau_eta", GenVisTau_eta, &b_GenVisTau_eta);
   fChain->SetBranchAddress("GenVisTau_mass", GenVisTau_mass, &b_GenVisTau_mass);
   fChain->SetBranchAddress("GenVisTau_phi", GenVisTau_phi, &b_GenVisTau_phi);
   fChain->SetBranchAddress("GenVisTau_pt", GenVisTau_pt, &b_GenVisTau_pt);
   fChain->SetBranchAddress("GenVisTau_charge", GenVisTau_charge, &b_GenVisTau_charge);
   fChain->SetBranchAddress("GenVisTau_genPartIdxMother", GenVisTau_genPartIdxMother, &b_GenVisTau_genPartIdxMother);
   fChain->SetBranchAddress("GenVisTau_status", GenVisTau_status, &b_GenVisTau_status);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("nEFTfitCoefficients", &nEFTfitCoefficients, &b_nEFTfitCoefficients);
   fChain->SetBranchAddress("EFTfitCoefficients", EFTfitCoefficients, &b_EFTfitCoefficients);
   fChain->SetBranchAddress("LHEWeight_originalXWGTUP", &LHEWeight_originalXWGTUP, &b_LHEWeight_originalXWGTUP);
   fChain->SetBranchAddress("nLHEPdfWeight", &nLHEPdfWeight, &b_nLHEPdfWeight);
   fChain->SetBranchAddress("LHEPdfWeight", LHEPdfWeight, &b_LHEPdfWeight);
   fChain->SetBranchAddress("nLHEReweightingWeight", &nLHEReweightingWeight, &b_nLHEReweightingWeight);
   fChain->SetBranchAddress("LHEReweightingWeight", &LHEReweightingWeight, &b_LHEReweightingWeight);
   fChain->SetBranchAddress("nLHEScaleWeight", &nLHEScaleWeight, &b_nLHEScaleWeight);
   fChain->SetBranchAddress("LHEScaleWeight", LHEScaleWeight, &b_LHEScaleWeight);
   fChain->SetBranchAddress("nPSWeight", &nPSWeight, &b_nPSWeight);
   fChain->SetBranchAddress("PSWeight", PSWeight, &b_PSWeight);
   fChain->SetBranchAddress("nWCnames", &nWCnames, &b_nWCnames);
   fChain->SetBranchAddress("WCnames", WCnames, &b_WCnames);
   fChain->SetBranchAddress("nIsoTrack", &nIsoTrack, &b_nIsoTrack);
   fChain->SetBranchAddress("IsoTrack_dxy", IsoTrack_dxy, &b_IsoTrack_dxy);
   fChain->SetBranchAddress("IsoTrack_dz", IsoTrack_dz, &b_IsoTrack_dz);
   fChain->SetBranchAddress("IsoTrack_eta", IsoTrack_eta, &b_IsoTrack_eta);
   fChain->SetBranchAddress("IsoTrack_pfRelIso03_all", IsoTrack_pfRelIso03_all, &b_IsoTrack_pfRelIso03_all);
   fChain->SetBranchAddress("IsoTrack_pfRelIso03_chg", IsoTrack_pfRelIso03_chg, &b_IsoTrack_pfRelIso03_chg);
   fChain->SetBranchAddress("IsoTrack_phi", IsoTrack_phi, &b_IsoTrack_phi);
   fChain->SetBranchAddress("IsoTrack_pt", IsoTrack_pt, &b_IsoTrack_pt);
   fChain->SetBranchAddress("IsoTrack_miniPFRelIso_all", IsoTrack_miniPFRelIso_all, &b_IsoTrack_miniPFRelIso_all);
   fChain->SetBranchAddress("IsoTrack_miniPFRelIso_chg", IsoTrack_miniPFRelIso_chg, &b_IsoTrack_miniPFRelIso_chg);
   fChain->SetBranchAddress("IsoTrack_charge", IsoTrack_charge, &b_IsoTrack_charge);
   fChain->SetBranchAddress("IsoTrack_fromPV", IsoTrack_fromPV, &b_IsoTrack_fromPV);
   fChain->SetBranchAddress("IsoTrack_pdgId", IsoTrack_pdgId, &b_IsoTrack_pdgId);
   fChain->SetBranchAddress("IsoTrack_isHighPurityTrack", IsoTrack_isHighPurityTrack, &b_IsoTrack_isHighPurityTrack);
   fChain->SetBranchAddress("IsoTrack_isPFcand", IsoTrack_isPFcand, &b_IsoTrack_isPFcand);
   fChain->SetBranchAddress("IsoTrack_isFromLostTrack", IsoTrack_isFromLostTrack, &b_IsoTrack_isFromLostTrack);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_btagCSVV2", Jet_btagCSVV2, &b_Jet_btagCSVV2);
   fChain->SetBranchAddress("Jet_btagDeepB", Jet_btagDeepB, &b_Jet_btagDeepB);
   fChain->SetBranchAddress("Jet_btagDeepCvB", Jet_btagDeepCvB, &b_Jet_btagDeepCvB);
   fChain->SetBranchAddress("Jet_btagDeepCvL", Jet_btagDeepCvL, &b_Jet_btagDeepCvL);
   fChain->SetBranchAddress("Jet_btagDeepFlavB", Jet_btagDeepFlavB, &b_Jet_btagDeepFlavB);
   fChain->SetBranchAddress("Jet_btagDeepFlavCvB", Jet_btagDeepFlavCvB, &b_Jet_btagDeepFlavCvB);
   fChain->SetBranchAddress("Jet_btagDeepFlavCvL", Jet_btagDeepFlavCvL, &b_Jet_btagDeepFlavCvL);
   fChain->SetBranchAddress("Jet_btagDeepFlavQG", Jet_btagDeepFlavQG, &b_Jet_btagDeepFlavQG);
   fChain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
   fChain->SetBranchAddress("Jet_chFPV0EF", Jet_chFPV0EF, &b_Jet_chFPV0EF);
   fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_hfsigmaEtaEta", Jet_hfsigmaEtaEta, &b_Jet_hfsigmaEtaEta);
   fChain->SetBranchAddress("Jet_hfsigmaPhiPhi", Jet_hfsigmaPhiPhi, &b_Jet_hfsigmaPhiPhi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_muEF", Jet_muEF, &b_Jet_muEF);
   fChain->SetBranchAddress("Jet_muonSubtrFactor", Jet_muonSubtrFactor, &b_Jet_muonSubtrFactor);
   fChain->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
   fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_puIdDisc", Jet_puIdDisc, &b_Jet_puIdDisc);
   fChain->SetBranchAddress("Jet_qgl", Jet_qgl, &b_Jet_qgl);
   fChain->SetBranchAddress("Jet_rawFactor", Jet_rawFactor, &b_Jet_rawFactor);
   fChain->SetBranchAddress("Jet_bRegCorr", Jet_bRegCorr, &b_Jet_bRegCorr);
   fChain->SetBranchAddress("Jet_bRegRes", Jet_bRegRes, &b_Jet_bRegRes);
   fChain->SetBranchAddress("Jet_cRegCorr", Jet_cRegCorr, &b_Jet_cRegCorr);
   fChain->SetBranchAddress("Jet_cRegRes", Jet_cRegRes, &b_Jet_cRegRes);
   fChain->SetBranchAddress("Jet_electronIdx1", Jet_electronIdx1, &b_Jet_electronIdx1);
   fChain->SetBranchAddress("Jet_electronIdx2", Jet_electronIdx2, &b_Jet_electronIdx2);
   fChain->SetBranchAddress("Jet_hfadjacentEtaStripsSize", Jet_hfadjacentEtaStripsSize, &b_Jet_hfadjacentEtaStripsSize);
   fChain->SetBranchAddress("Jet_hfcentralEtaStripSize", Jet_hfcentralEtaStripSize, &b_Jet_hfcentralEtaStripSize);
   fChain->SetBranchAddress("Jet_jetId", Jet_jetId, &b_Jet_jetId);
   fChain->SetBranchAddress("Jet_muonIdx1", Jet_muonIdx1, &b_Jet_muonIdx1);
   fChain->SetBranchAddress("Jet_muonIdx2", Jet_muonIdx2, &b_Jet_muonIdx2);
   fChain->SetBranchAddress("Jet_nElectrons", Jet_nElectrons, &b_Jet_nElectrons);
   fChain->SetBranchAddress("Jet_nMuons", Jet_nMuons, &b_Jet_nMuons);
   fChain->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);
   fChain->SetBranchAddress("Jet_nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
   fChain->SetBranchAddress("L1PreFiringWeight_Dn", &L1PreFiringWeight_Dn, &b_L1PreFiringWeight_Dn);
   fChain->SetBranchAddress("L1PreFiringWeight_ECAL_Dn", &L1PreFiringWeight_ECAL_Dn, &b_L1PreFiringWeight_ECAL_Dn);
   fChain->SetBranchAddress("L1PreFiringWeight_ECAL_Nom", &L1PreFiringWeight_ECAL_Nom, &b_L1PreFiringWeight_ECAL_Nom);
   fChain->SetBranchAddress("L1PreFiringWeight_ECAL_Up", &L1PreFiringWeight_ECAL_Up, &b_L1PreFiringWeight_ECAL_Up);
   fChain->SetBranchAddress("L1PreFiringWeight_Muon_Nom", &L1PreFiringWeight_Muon_Nom, &b_L1PreFiringWeight_Muon_Nom);
   fChain->SetBranchAddress("L1PreFiringWeight_Muon_StatDn", &L1PreFiringWeight_Muon_StatDn, &b_L1PreFiringWeight_Muon_StatDn);
   fChain->SetBranchAddress("L1PreFiringWeight_Muon_StatUp", &L1PreFiringWeight_Muon_StatUp, &b_L1PreFiringWeight_Muon_StatUp);
   fChain->SetBranchAddress("L1PreFiringWeight_Muon_SystDn", &L1PreFiringWeight_Muon_SystDn, &b_L1PreFiringWeight_Muon_SystDn);
   fChain->SetBranchAddress("L1PreFiringWeight_Muon_SystUp", &L1PreFiringWeight_Muon_SystUp, &b_L1PreFiringWeight_Muon_SystUp);
   fChain->SetBranchAddress("L1PreFiringWeight_Nom", &L1PreFiringWeight_Nom, &b_L1PreFiringWeight_Nom);
   fChain->SetBranchAddress("L1PreFiringWeight_Up", &L1PreFiringWeight_Up, &b_L1PreFiringWeight_Up);
   fChain->SetBranchAddress("LHE_HT", &LHE_HT, &b_LHE_HT);
   fChain->SetBranchAddress("LHE_HTIncoming", &LHE_HTIncoming, &b_LHE_HTIncoming);
   fChain->SetBranchAddress("LHE_Vpt", &LHE_Vpt, &b_LHE_Vpt);
   fChain->SetBranchAddress("LHE_AlphaS", &LHE_AlphaS, &b_LHE_AlphaS);
   fChain->SetBranchAddress("LHE_Njets", &LHE_Njets, &b_LHE_Njets);
   fChain->SetBranchAddress("LHE_Nb", &LHE_Nb, &b_LHE_Nb);
   fChain->SetBranchAddress("LHE_Nc", &LHE_Nc, &b_LHE_Nc);
   fChain->SetBranchAddress("LHE_Nuds", &LHE_Nuds, &b_LHE_Nuds);
   fChain->SetBranchAddress("LHE_Nglu", &LHE_Nglu, &b_LHE_Nglu);
   fChain->SetBranchAddress("LHE_NpNLO", &LHE_NpNLO, &b_LHE_NpNLO);
   fChain->SetBranchAddress("LHE_NpLO", &LHE_NpLO, &b_LHE_NpLO);
   fChain->SetBranchAddress("nLHEPart", &nLHEPart, &b_nLHEPart);
   fChain->SetBranchAddress("LHEPart_pt", LHEPart_pt, &b_LHEPart_pt);
   fChain->SetBranchAddress("LHEPart_eta", LHEPart_eta, &b_LHEPart_eta);
   fChain->SetBranchAddress("LHEPart_phi", LHEPart_phi, &b_LHEPart_phi);
   fChain->SetBranchAddress("LHEPart_mass", LHEPart_mass, &b_LHEPart_mass);
   fChain->SetBranchAddress("LHEPart_incomingpz", LHEPart_incomingpz, &b_LHEPart_incomingpz);
   fChain->SetBranchAddress("LHEPart_pdgId", LHEPart_pdgId, &b_LHEPart_pdgId);
   fChain->SetBranchAddress("LHEPart_status", LHEPart_status, &b_LHEPart_status);
   fChain->SetBranchAddress("LHEPart_spin", LHEPart_spin, &b_LHEPart_spin);
   fChain->SetBranchAddress("nLowPtElectron", &nLowPtElectron, &b_nLowPtElectron);
   fChain->SetBranchAddress("LowPtElectron_ID", LowPtElectron_ID, &b_LowPtElectron_ID);
   fChain->SetBranchAddress("LowPtElectron_convVtxRadius", LowPtElectron_convVtxRadius, &b_LowPtElectron_convVtxRadius);
   fChain->SetBranchAddress("LowPtElectron_deltaEtaSC", LowPtElectron_deltaEtaSC, &b_LowPtElectron_deltaEtaSC);
   fChain->SetBranchAddress("LowPtElectron_dxy", LowPtElectron_dxy, &b_LowPtElectron_dxy);
   fChain->SetBranchAddress("LowPtElectron_dxyErr", LowPtElectron_dxyErr, &b_LowPtElectron_dxyErr);
   fChain->SetBranchAddress("LowPtElectron_dz", LowPtElectron_dz, &b_LowPtElectron_dz);
   fChain->SetBranchAddress("LowPtElectron_dzErr", LowPtElectron_dzErr, &b_LowPtElectron_dzErr);
   fChain->SetBranchAddress("LowPtElectron_eInvMinusPInv", LowPtElectron_eInvMinusPInv, &b_LowPtElectron_eInvMinusPInv);
   fChain->SetBranchAddress("LowPtElectron_embeddedID", LowPtElectron_embeddedID, &b_LowPtElectron_embeddedID);
   fChain->SetBranchAddress("LowPtElectron_energyErr", LowPtElectron_energyErr, &b_LowPtElectron_energyErr);
   fChain->SetBranchAddress("LowPtElectron_eta", LowPtElectron_eta, &b_LowPtElectron_eta);
   fChain->SetBranchAddress("LowPtElectron_hoe", LowPtElectron_hoe, &b_LowPtElectron_hoe);
   fChain->SetBranchAddress("LowPtElectron_mass", LowPtElectron_mass, &b_LowPtElectron_mass);
   fChain->SetBranchAddress("LowPtElectron_miniPFRelIso_all", LowPtElectron_miniPFRelIso_all, &b_LowPtElectron_miniPFRelIso_all);
   fChain->SetBranchAddress("LowPtElectron_miniPFRelIso_chg", LowPtElectron_miniPFRelIso_chg, &b_LowPtElectron_miniPFRelIso_chg);
   fChain->SetBranchAddress("LowPtElectron_phi", LowPtElectron_phi, &b_LowPtElectron_phi);
   fChain->SetBranchAddress("LowPtElectron_pt", LowPtElectron_pt, &b_LowPtElectron_pt);
   fChain->SetBranchAddress("LowPtElectron_ptbiased", LowPtElectron_ptbiased, &b_LowPtElectron_ptbiased);
   fChain->SetBranchAddress("LowPtElectron_r9", LowPtElectron_r9, &b_LowPtElectron_r9);
   fChain->SetBranchAddress("LowPtElectron_scEtOverPt", LowPtElectron_scEtOverPt, &b_LowPtElectron_scEtOverPt);
   fChain->SetBranchAddress("LowPtElectron_sieie", LowPtElectron_sieie, &b_LowPtElectron_sieie);
   fChain->SetBranchAddress("LowPtElectron_unbiased", LowPtElectron_unbiased, &b_LowPtElectron_unbiased);
   fChain->SetBranchAddress("LowPtElectron_charge", LowPtElectron_charge, &b_LowPtElectron_charge);
   fChain->SetBranchAddress("LowPtElectron_convWP", LowPtElectron_convWP, &b_LowPtElectron_convWP);
   fChain->SetBranchAddress("LowPtElectron_pdgId", LowPtElectron_pdgId, &b_LowPtElectron_pdgId);
   fChain->SetBranchAddress("LowPtElectron_convVeto", LowPtElectron_convVeto, &b_LowPtElectron_convVeto);
   fChain->SetBranchAddress("LowPtElectron_lostHits", LowPtElectron_lostHits, &b_LowPtElectron_lostHits);
   fChain->SetBranchAddress("GenMET_phi", &GenMET_phi, &b_GenMET_phi);
   fChain->SetBranchAddress("GenMET_pt", &GenMET_pt, &b_GenMET_pt);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaX", &MET_MetUnclustEnUpDeltaX, &b_MET_MetUnclustEnUpDeltaX);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaY", &MET_MetUnclustEnUpDeltaY, &b_MET_MetUnclustEnUpDeltaY);
   fChain->SetBranchAddress("MET_covXX", &MET_covXX, &b_MET_covXX);
   fChain->SetBranchAddress("MET_covXY", &MET_covXY, &b_MET_covXY);
   fChain->SetBranchAddress("MET_covYY", &MET_covYY, &b_MET_covYY);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   fChain->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance);
   fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   fChain->SetBranchAddress("MET_sumPtUnclustered", &MET_sumPtUnclustered, &b_MET_sumPtUnclustered);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_dxy", Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dxyErr", Muon_dxyErr, &b_Muon_dxyErr);
   fChain->SetBranchAddress("Muon_dxybs", Muon_dxybs, &b_Muon_dxybs);
   fChain->SetBranchAddress("Muon_dz", Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_dzErr", Muon_dzErr, &b_Muon_dzErr);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_ip3d", Muon_ip3d, &b_Muon_ip3d);
   fChain->SetBranchAddress("Muon_jetPtRelv2", Muon_jetPtRelv2, &b_Muon_jetPtRelv2);
   fChain->SetBranchAddress("Muon_jetRelIso", Muon_jetRelIso, &b_Muon_jetRelIso);
   fChain->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all, &b_Muon_miniPFRelIso_all);
   fChain->SetBranchAddress("Muon_miniPFRelIso_chg", Muon_miniPFRelIso_chg, &b_Muon_miniPFRelIso_chg);
   fChain->SetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
   fChain->SetBranchAddress("Muon_pfRelIso03_chg", Muon_pfRelIso03_chg, &b_Muon_pfRelIso03_chg);
   fChain->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_ptErr", Muon_ptErr, &b_Muon_ptErr);
   fChain->SetBranchAddress("Muon_segmentComp", Muon_segmentComp, &b_Muon_segmentComp);
   fChain->SetBranchAddress("Muon_sip3d", Muon_sip3d, &b_Muon_sip3d);
   fChain->SetBranchAddress("Muon_softMva", Muon_softMva, &b_Muon_softMva);
   fChain->SetBranchAddress("Muon_tkRelIso", Muon_tkRelIso, &b_Muon_tkRelIso);
   fChain->SetBranchAddress("Muon_tunepRelPt", Muon_tunepRelPt, &b_Muon_tunepRelPt);
   fChain->SetBranchAddress("Muon_mvaLowPt", Muon_mvaLowPt, &b_Muon_mvaLowPt);
   fChain->SetBranchAddress("Muon_mvaTTH", Muon_mvaTTH, &b_Muon_mvaTTH);
   fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_jetIdx", Muon_jetIdx, &b_Muon_jetIdx);
   fChain->SetBranchAddress("Muon_nStations", Muon_nStations, &b_Muon_nStations);
   fChain->SetBranchAddress("Muon_nTrackerLayers", Muon_nTrackerLayers, &b_Muon_nTrackerLayers);
   fChain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_tightCharge", Muon_tightCharge, &b_Muon_tightCharge);
   fChain->SetBranchAddress("Muon_fsrPhotonIdx", Muon_fsrPhotonIdx, &b_Muon_fsrPhotonIdx);
   fChain->SetBranchAddress("Muon_highPtId", Muon_highPtId, &b_Muon_highPtId);
   fChain->SetBranchAddress("Muon_highPurity", Muon_highPurity, &b_Muon_highPurity);
   fChain->SetBranchAddress("Muon_inTimeMuon", Muon_inTimeMuon, &b_Muon_inTimeMuon);
   fChain->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_isPFcand", Muon_isPFcand, &b_Muon_isPFcand);
   fChain->SetBranchAddress("Muon_isStandalone", Muon_isStandalone, &b_Muon_isStandalone);
   fChain->SetBranchAddress("Muon_isTracker", Muon_isTracker, &b_Muon_isTracker);
   fChain->SetBranchAddress("Muon_jetNDauCharged", Muon_jetNDauCharged, &b_Muon_jetNDauCharged);
   fChain->SetBranchAddress("Muon_looseId", Muon_looseId, &b_Muon_looseId);
   fChain->SetBranchAddress("Muon_mediumId", Muon_mediumId, &b_Muon_mediumId);
   fChain->SetBranchAddress("Muon_mediumPromptId", Muon_mediumPromptId, &b_Muon_mediumPromptId);
   fChain->SetBranchAddress("Muon_miniIsoId", Muon_miniIsoId, &b_Muon_miniIsoId);
   fChain->SetBranchAddress("Muon_multiIsoId", Muon_multiIsoId, &b_Muon_multiIsoId);
   fChain->SetBranchAddress("Muon_mvaId", Muon_mvaId, &b_Muon_mvaId);
   fChain->SetBranchAddress("Muon_mvaLowPtId", Muon_mvaLowPtId, &b_Muon_mvaLowPtId);
   fChain->SetBranchAddress("Muon_pfIsoId", Muon_pfIsoId, &b_Muon_pfIsoId);
   fChain->SetBranchAddress("Muon_puppiIsoId", Muon_puppiIsoId, &b_Muon_puppiIsoId);
   fChain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
   fChain->SetBranchAddress("Muon_softMvaId", Muon_softMvaId, &b_Muon_softMvaId);
   fChain->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
   fChain->SetBranchAddress("Muon_tkIsoId", Muon_tkIsoId, &b_Muon_tkIsoId);
   fChain->SetBranchAddress("Muon_triggerIdLoose", Muon_triggerIdLoose, &b_Muon_triggerIdLoose);
   fChain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
   fChain->SetBranchAddress("Photon_dEscaleDown", Photon_dEscaleDown, &b_Photon_dEscaleDown);
   fChain->SetBranchAddress("Photon_dEscaleUp", Photon_dEscaleUp, &b_Photon_dEscaleUp);
   fChain->SetBranchAddress("Photon_dEsigmaDown", Photon_dEsigmaDown, &b_Photon_dEsigmaDown);
   fChain->SetBranchAddress("Photon_dEsigmaUp", Photon_dEsigmaUp, &b_Photon_dEsigmaUp);
   fChain->SetBranchAddress("Photon_eCorr", Photon_eCorr, &b_Photon_eCorr);
   fChain->SetBranchAddress("Photon_energyErr", Photon_energyErr, &b_Photon_energyErr);
   fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_hoe", Photon_hoe, &b_Photon_hoe);
   fChain->SetBranchAddress("Photon_mass", Photon_mass, &b_Photon_mass);
   fChain->SetBranchAddress("Photon_mvaID", Photon_mvaID, &b_Photon_mvaID);
   fChain->SetBranchAddress("Photon_mvaID_Fall17V1p1", Photon_mvaID_Fall17V1p1, &b_Photon_mvaID_Fall17V1p1);
   fChain->SetBranchAddress("Photon_pfRelIso03_all", Photon_pfRelIso03_all, &b_Photon_pfRelIso03_all);
   fChain->SetBranchAddress("Photon_pfRelIso03_chg", Photon_pfRelIso03_chg, &b_Photon_pfRelIso03_chg);
   fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_r9", Photon_r9, &b_Photon_r9);
   fChain->SetBranchAddress("Photon_sieie", Photon_sieie, &b_Photon_sieie);
   fChain->SetBranchAddress("Photon_charge", Photon_charge, &b_Photon_charge);
   fChain->SetBranchAddress("Photon_cutBased", Photon_cutBased, &b_Photon_cutBased);
   fChain->SetBranchAddress("Photon_cutBased_Fall17V1Bitmap", Photon_cutBased_Fall17V1Bitmap, &b_Photon_cutBased_Fall17V1Bitmap);
   fChain->SetBranchAddress("Photon_electronIdx", Photon_electronIdx, &b_Photon_electronIdx);
   fChain->SetBranchAddress("Photon_jetIdx", Photon_jetIdx, &b_Photon_jetIdx);
   fChain->SetBranchAddress("Photon_pdgId", Photon_pdgId, &b_Photon_pdgId);
   fChain->SetBranchAddress("Photon_vidNestedWPBitmap", Photon_vidNestedWPBitmap, &b_Photon_vidNestedWPBitmap);
   fChain->SetBranchAddress("Photon_electronVeto", Photon_electronVeto, &b_Photon_electronVeto);
   fChain->SetBranchAddress("Photon_isScEtaEB", Photon_isScEtaEB, &b_Photon_isScEtaEB);
   fChain->SetBranchAddress("Photon_isScEtaEE", Photon_isScEtaEE, &b_Photon_isScEtaEE);
   fChain->SetBranchAddress("Photon_mvaID_WP80", Photon_mvaID_WP80, &b_Photon_mvaID_WP80);
   fChain->SetBranchAddress("Photon_mvaID_WP90", Photon_mvaID_WP90, &b_Photon_mvaID_WP90);
   fChain->SetBranchAddress("Photon_pixelSeed", Photon_pixelSeed, &b_Photon_pixelSeed);
   fChain->SetBranchAddress("Photon_seedGain", Photon_seedGain, &b_Photon_seedGain);
   fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
   fChain->SetBranchAddress("Pileup_pudensity", &Pileup_pudensity, &b_Pileup_pudensity);
   fChain->SetBranchAddress("Pileup_gpudensity", &Pileup_gpudensity, &b_Pileup_gpudensity);
   fChain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
   fChain->SetBranchAddress("Pileup_sumEOOT", &Pileup_sumEOOT, &b_Pileup_sumEOOT);
   fChain->SetBranchAddress("Pileup_sumLOOT", &Pileup_sumLOOT, &b_Pileup_sumLOOT);
   fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
   fChain->SetBranchAddress("PuppiMET_phiJERDown", &PuppiMET_phiJERDown, &b_PuppiMET_phiJERDown);
   fChain->SetBranchAddress("PuppiMET_phiJERUp", &PuppiMET_phiJERUp, &b_PuppiMET_phiJERUp);
   fChain->SetBranchAddress("PuppiMET_phiJESDown", &PuppiMET_phiJESDown, &b_PuppiMET_phiJESDown);
   fChain->SetBranchAddress("PuppiMET_phiJESUp", &PuppiMET_phiJESUp, &b_PuppiMET_phiJESUp);
   fChain->SetBranchAddress("PuppiMET_phiUnclusteredDown", &PuppiMET_phiUnclusteredDown, &b_PuppiMET_phiUnclusteredDown);
   fChain->SetBranchAddress("PuppiMET_phiUnclusteredUp", &PuppiMET_phiUnclusteredUp, &b_PuppiMET_phiUnclusteredUp);
   fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
   fChain->SetBranchAddress("PuppiMET_ptJERDown", &PuppiMET_ptJERDown, &b_PuppiMET_ptJERDown);
   fChain->SetBranchAddress("PuppiMET_ptJERUp", &PuppiMET_ptJERUp, &b_PuppiMET_ptJERUp);
   fChain->SetBranchAddress("PuppiMET_ptJESDown", &PuppiMET_ptJESDown, &b_PuppiMET_ptJESDown);
   fChain->SetBranchAddress("PuppiMET_ptJESUp", &PuppiMET_ptJESUp, &b_PuppiMET_ptJESUp);
   fChain->SetBranchAddress("PuppiMET_ptUnclusteredDown", &PuppiMET_ptUnclusteredDown, &b_PuppiMET_ptUnclusteredDown);
   fChain->SetBranchAddress("PuppiMET_ptUnclusteredUp", &PuppiMET_ptUnclusteredUp, &b_PuppiMET_ptUnclusteredUp);
   fChain->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt, &b_PuppiMET_sumEt);
   fChain->SetBranchAddress("RawMET_phi", &RawMET_phi, &b_RawMET_phi);
   fChain->SetBranchAddress("RawMET_pt", &RawMET_pt, &b_RawMET_pt);
   fChain->SetBranchAddress("RawMET_sumEt", &RawMET_sumEt, &b_RawMET_sumEt);
   fChain->SetBranchAddress("RawPuppiMET_phi", &RawPuppiMET_phi, &b_RawPuppiMET_phi);
   fChain->SetBranchAddress("RawPuppiMET_pt", &RawPuppiMET_pt, &b_RawPuppiMET_pt);
   fChain->SetBranchAddress("RawPuppiMET_sumEt", &RawPuppiMET_sumEt, &b_RawPuppiMET_sumEt);
   fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentral", &fixedGridRhoFastjetCentral, &b_fixedGridRhoFastjetCentral);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, &b_fixedGridRhoFastjetCentralChargedPileUp);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("nGenDressedLepton", &nGenDressedLepton, &b_nGenDressedLepton);
   fChain->SetBranchAddress("GenDressedLepton_eta", GenDressedLepton_eta, &b_GenDressedLepton_eta);
   fChain->SetBranchAddress("GenDressedLepton_mass", GenDressedLepton_mass, &b_GenDressedLepton_mass);
   fChain->SetBranchAddress("GenDressedLepton_phi", GenDressedLepton_phi, &b_GenDressedLepton_phi);
   fChain->SetBranchAddress("GenDressedLepton_pt", GenDressedLepton_pt, &b_GenDressedLepton_pt);
   fChain->SetBranchAddress("GenDressedLepton_pdgId", GenDressedLepton_pdgId, &b_GenDressedLepton_pdgId);
   fChain->SetBranchAddress("GenDressedLepton_hasTauAnc", GenDressedLepton_hasTauAnc, &b_GenDressedLepton_hasTauAnc);
   fChain->SetBranchAddress("nGenIsolatedPhoton", &nGenIsolatedPhoton, &b_nGenIsolatedPhoton);
   fChain->SetBranchAddress("GenIsolatedPhoton_eta", GenIsolatedPhoton_eta, &b_GenIsolatedPhoton_eta);
   fChain->SetBranchAddress("GenIsolatedPhoton_mass", GenIsolatedPhoton_mass, &b_GenIsolatedPhoton_mass);
   fChain->SetBranchAddress("GenIsolatedPhoton_phi", GenIsolatedPhoton_phi, &b_GenIsolatedPhoton_phi);
   fChain->SetBranchAddress("GenIsolatedPhoton_pt", GenIsolatedPhoton_pt, &b_GenIsolatedPhoton_pt);
   fChain->SetBranchAddress("nSoftActivityJet", &nSoftActivityJet, &b_nSoftActivityJet);
   fChain->SetBranchAddress("SoftActivityJet_eta", SoftActivityJet_eta, &b_SoftActivityJet_eta);
   fChain->SetBranchAddress("SoftActivityJet_phi", SoftActivityJet_phi, &b_SoftActivityJet_phi);
   fChain->SetBranchAddress("SoftActivityJet_pt", SoftActivityJet_pt, &b_SoftActivityJet_pt);
   fChain->SetBranchAddress("SoftActivityJetHT", &SoftActivityJetHT, &b_SoftActivityJetHT);
   fChain->SetBranchAddress("SoftActivityJetHT10", &SoftActivityJetHT10, &b_SoftActivityJetHT10);
   fChain->SetBranchAddress("SoftActivityJetHT2", &SoftActivityJetHT2, &b_SoftActivityJetHT2);
   fChain->SetBranchAddress("SoftActivityJetHT5", &SoftActivityJetHT5, &b_SoftActivityJetHT5);
   fChain->SetBranchAddress("SoftActivityJetNjets10", &SoftActivityJetNjets10, &b_SoftActivityJetNjets10);
   fChain->SetBranchAddress("SoftActivityJetNjets2", &SoftActivityJetNjets2, &b_SoftActivityJetNjets2);
   fChain->SetBranchAddress("SoftActivityJetNjets5", &SoftActivityJetNjets5, &b_SoftActivityJetNjets5);
   fChain->SetBranchAddress("nSubJet", &nSubJet, &b_nSubJet);
   fChain->SetBranchAddress("SubJet_btagCSVV2", SubJet_btagCSVV2, &b_SubJet_btagCSVV2);
   fChain->SetBranchAddress("SubJet_btagDeepB", SubJet_btagDeepB, &b_SubJet_btagDeepB);
   fChain->SetBranchAddress("SubJet_eta", SubJet_eta, &b_SubJet_eta);
   fChain->SetBranchAddress("SubJet_mass", SubJet_mass, &b_SubJet_mass);
   fChain->SetBranchAddress("SubJet_n2b1", SubJet_n2b1, &b_SubJet_n2b1);
   fChain->SetBranchAddress("SubJet_n3b1", SubJet_n3b1, &b_SubJet_n3b1);
   fChain->SetBranchAddress("SubJet_phi", SubJet_phi, &b_SubJet_phi);
   fChain->SetBranchAddress("SubJet_pt", SubJet_pt, &b_SubJet_pt);
   fChain->SetBranchAddress("SubJet_rawFactor", SubJet_rawFactor, &b_SubJet_rawFactor);
   fChain->SetBranchAddress("SubJet_tau1", SubJet_tau1, &b_SubJet_tau1);
   fChain->SetBranchAddress("SubJet_tau2", SubJet_tau2, &b_SubJet_tau2);
   fChain->SetBranchAddress("SubJet_tau3", SubJet_tau3, &b_SubJet_tau3);
   fChain->SetBranchAddress("SubJet_tau4", SubJet_tau4, &b_SubJet_tau4);
   fChain->SetBranchAddress("nTau", &nTau, &b_nTau);
   fChain->SetBranchAddress("Tau_chargedIso", Tau_chargedIso, &b_Tau_chargedIso);
   fChain->SetBranchAddress("Tau_dxy", Tau_dxy, &b_Tau_dxy);
   fChain->SetBranchAddress("Tau_dz", Tau_dz, &b_Tau_dz);
   fChain->SetBranchAddress("Tau_eta", Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau_leadTkDeltaEta", Tau_leadTkDeltaEta, &b_Tau_leadTkDeltaEta);
   fChain->SetBranchAddress("Tau_leadTkDeltaPhi", Tau_leadTkDeltaPhi, &b_Tau_leadTkDeltaPhi);
   fChain->SetBranchAddress("Tau_leadTkPtOverTauPt", Tau_leadTkPtOverTauPt, &b_Tau_leadTkPtOverTauPt);
   fChain->SetBranchAddress("Tau_mass", Tau_mass, &b_Tau_mass);
   fChain->SetBranchAddress("Tau_neutralIso", Tau_neutralIso, &b_Tau_neutralIso);
   fChain->SetBranchAddress("Tau_phi", Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau_photonsOutsideSignalCone", Tau_photonsOutsideSignalCone, &b_Tau_photonsOutsideSignalCone);
   fChain->SetBranchAddress("Tau_pt", Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau_puCorr", Tau_puCorr, &b_Tau_puCorr);
   fChain->SetBranchAddress("Tau_rawDeepTau2017v2p1VSe", Tau_rawDeepTau2017v2p1VSe, &b_Tau_rawDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("Tau_rawDeepTau2017v2p1VSjet", Tau_rawDeepTau2017v2p1VSjet, &b_Tau_rawDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("Tau_rawDeepTau2017v2p1VSmu", Tau_rawDeepTau2017v2p1VSmu, &b_Tau_rawDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("Tau_rawIso", Tau_rawIso, &b_Tau_rawIso);
   fChain->SetBranchAddress("Tau_rawIsodR03", Tau_rawIsodR03, &b_Tau_rawIsodR03);
   fChain->SetBranchAddress("Tau_charge", Tau_charge, &b_Tau_charge);
   fChain->SetBranchAddress("Tau_decayMode", Tau_decayMode, &b_Tau_decayMode);
   fChain->SetBranchAddress("Tau_jetIdx", Tau_jetIdx, &b_Tau_jetIdx);
   fChain->SetBranchAddress("Tau_idAntiEleDeadECal", Tau_idAntiEleDeadECal, &b_Tau_idAntiEleDeadECal);
   fChain->SetBranchAddress("Tau_idAntiMu", Tau_idAntiMu, &b_Tau_idAntiMu);
   fChain->SetBranchAddress("Tau_idDecayModeOldDMs", Tau_idDecayModeOldDMs, &b_Tau_idDecayModeOldDMs);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSe", Tau_idDeepTau2017v2p1VSe, &b_Tau_idDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSjet", Tau_idDeepTau2017v2p1VSjet, &b_Tau_idDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSmu", Tau_idDeepTau2017v2p1VSmu, &b_Tau_idDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("TkMET_phi", &TkMET_phi, &b_TkMET_phi);
   fChain->SetBranchAddress("TkMET_pt", &TkMET_pt, &b_TkMET_pt);
   fChain->SetBranchAddress("TkMET_sumEt", &TkMET_sumEt, &b_TkMET_sumEt);
   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_l1pt", TrigObj_l1pt, &b_TrigObj_l1pt);
   fChain->SetBranchAddress("TrigObj_l1pt_2", TrigObj_l1pt_2, &b_TrigObj_l1pt_2);
   fChain->SetBranchAddress("TrigObj_l2pt", TrigObj_l2pt, &b_TrigObj_l2pt);
   fChain->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
   fChain->SetBranchAddress("TrigObj_l1iso", TrigObj_l1iso, &b_TrigObj_l1iso);
   fChain->SetBranchAddress("TrigObj_l1charge", TrigObj_l1charge, &b_TrigObj_l1charge);
   fChain->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits, &b_TrigObj_filterBits);
   fChain->SetBranchAddress("genTtbarId", &genTtbarId, &b_genTtbarId);
   fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   fChain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   fChain->SetBranchAddress("SV_dxy", SV_dxy, &b_SV_dxy);
   fChain->SetBranchAddress("SV_dxySig", SV_dxySig, &b_SV_dxySig);
   fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   fChain->SetBranchAddress("SV_charge", SV_charge, &b_SV_charge);
   fChain->SetBranchAddress("boostedTau_genPartIdx", boostedTau_genPartIdx, &b_boostedTau_genPartIdx);
   fChain->SetBranchAddress("boostedTau_genPartFlav", boostedTau_genPartFlav, &b_boostedTau_genPartFlav);
   fChain->SetBranchAddress("Electron_genPartIdx", Electron_genPartIdx, &b_Electron_genPartIdx);
   fChain->SetBranchAddress("Electron_genPartFlav", Electron_genPartFlav, &b_Electron_genPartFlav);
   fChain->SetBranchAddress("FatJet_genJetAK8Idx", FatJet_genJetAK8Idx, &b_FatJet_genJetAK8Idx);
   fChain->SetBranchAddress("FatJet_hadronFlavour", FatJet_hadronFlavour, &b_FatJet_hadronFlavour);
   fChain->SetBranchAddress("FatJet_nBHadrons", FatJet_nBHadrons, &b_FatJet_nBHadrons);
   fChain->SetBranchAddress("FatJet_nCHadrons", FatJet_nCHadrons, &b_FatJet_nCHadrons);
   fChain->SetBranchAddress("GenJetAK8_partonFlavour", GenJetAK8_partonFlavour, &b_GenJetAK8_partonFlavour);
   fChain->SetBranchAddress("GenJetAK8_hadronFlavour", GenJetAK8_hadronFlavour, &b_GenJetAK8_hadronFlavour);
   fChain->SetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour, &b_GenJet_partonFlavour);
   fChain->SetBranchAddress("GenJet_hadronFlavour", GenJet_hadronFlavour, &b_GenJet_hadronFlavour);
   fChain->SetBranchAddress("GenVtx_t0", &GenVtx_t0, &b_GenVtx_t0);
   fChain->SetBranchAddress("Jet_genJetIdx", Jet_genJetIdx, &b_Jet_genJetIdx);
   fChain->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
   fChain->SetBranchAddress("Jet_partonFlavour", Jet_partonFlavour, &b_Jet_partonFlavour);
   fChain->SetBranchAddress("LowPtElectron_genPartIdx", LowPtElectron_genPartIdx, &b_LowPtElectron_genPartIdx);
   fChain->SetBranchAddress("LowPtElectron_genPartFlav", LowPtElectron_genPartFlav, &b_LowPtElectron_genPartFlav);
   fChain->SetBranchAddress("Muon_genPartIdx", Muon_genPartIdx, &b_Muon_genPartIdx);
   fChain->SetBranchAddress("Muon_genPartFlav", Muon_genPartFlav, &b_Muon_genPartFlav);
   fChain->SetBranchAddress("Photon_genPartIdx", Photon_genPartIdx, &b_Photon_genPartIdx);
   fChain->SetBranchAddress("Photon_genPartFlav", Photon_genPartFlav, &b_Photon_genPartFlav);
   fChain->SetBranchAddress("MET_fiducialGenPhi", &MET_fiducialGenPhi, &b_MET_fiducialGenPhi);
   fChain->SetBranchAddress("MET_fiducialGenPt", &MET_fiducialGenPt, &b_MET_fiducialGenPt);
   fChain->SetBranchAddress("Electron_cleanmask", Electron_cleanmask, &b_Electron_cleanmask);
   fChain->SetBranchAddress("Jet_cleanmask", Jet_cleanmask, &b_Jet_cleanmask);
   fChain->SetBranchAddress("Muon_cleanmask", Muon_cleanmask, &b_Muon_cleanmask);
   fChain->SetBranchAddress("Photon_cleanmask", Photon_cleanmask, &b_Photon_cleanmask);
   fChain->SetBranchAddress("Tau_cleanmask", Tau_cleanmask, &b_Tau_cleanmask);
   fChain->SetBranchAddress("SubJet_hadronFlavour", SubJet_hadronFlavour, &b_SubJet_hadronFlavour);
   fChain->SetBranchAddress("SubJet_nBHadrons", SubJet_nBHadrons, &b_SubJet_nBHadrons);
   fChain->SetBranchAddress("SubJet_nCHadrons", SubJet_nCHadrons, &b_SubJet_nCHadrons);
   fChain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   fChain->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
   fChain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   fChain->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
   fChain->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
   fChain->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
   fChain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   fChain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   fChain->SetBranchAddress("SV_z", SV_z, &b_SV_z);
   fChain->SetBranchAddress("SV_ntracks", SV_ntracks, &b_SV_ntracks);
   fChain->SetBranchAddress("Tau_genPartIdx", Tau_genPartIdx, &b_Tau_genPartIdx);
   fChain->SetBranchAddress("Tau_genPartFlav", Tau_genPartFlav, &b_Tau_genPartFlav);
   fChain->SetBranchAddress("L1_AlwaysTrue", &L1_AlwaysTrue, &b_L1_AlwaysTrue);
   fChain->SetBranchAddress("L1_BPTX_AND_Ref1_VME", &L1_BPTX_AND_Ref1_VME, &b_L1_BPTX_AND_Ref1_VME);
   fChain->SetBranchAddress("L1_BPTX_AND_Ref3_VME", &L1_BPTX_AND_Ref3_VME, &b_L1_BPTX_AND_Ref3_VME);
   fChain->SetBranchAddress("L1_BPTX_AND_Ref4_VME", &L1_BPTX_AND_Ref4_VME, &b_L1_BPTX_AND_Ref4_VME);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_B1_VME", &L1_BPTX_BeamGas_B1_VME, &b_L1_BPTX_BeamGas_B1_VME);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_B2_VME", &L1_BPTX_BeamGas_B2_VME, &b_L1_BPTX_BeamGas_B2_VME);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_Ref1_VME", &L1_BPTX_BeamGas_Ref1_VME, &b_L1_BPTX_BeamGas_Ref1_VME);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_Ref2_VME", &L1_BPTX_BeamGas_Ref2_VME, &b_L1_BPTX_BeamGas_Ref2_VME);
   fChain->SetBranchAddress("L1_BPTX_NotOR_VME", &L1_BPTX_NotOR_VME, &b_L1_BPTX_NotOR_VME);
   fChain->SetBranchAddress("L1_BPTX_OR_Ref3_VME", &L1_BPTX_OR_Ref3_VME, &b_L1_BPTX_OR_Ref3_VME);
   fChain->SetBranchAddress("L1_BPTX_OR_Ref4_VME", &L1_BPTX_OR_Ref4_VME, &b_L1_BPTX_OR_Ref4_VME);
   fChain->SetBranchAddress("L1_BPTX_RefAND_VME", &L1_BPTX_RefAND_VME, &b_L1_BPTX_RefAND_VME);
   fChain->SetBranchAddress("L1_BptxMinus", &L1_BptxMinus, &b_L1_BptxMinus);
   fChain->SetBranchAddress("L1_BptxOR", &L1_BptxOR, &b_L1_BptxOR);
   fChain->SetBranchAddress("L1_BptxPlus", &L1_BptxPlus, &b_L1_BptxPlus);
   fChain->SetBranchAddress("L1_BptxXOR", &L1_BptxXOR, &b_L1_BptxXOR);
   fChain->SetBranchAddress("L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142, &b_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142);
   fChain->SetBranchAddress("L1_DoubleEG6_HTT240er", &L1_DoubleEG6_HTT240er, &b_L1_DoubleEG6_HTT240er);
   fChain->SetBranchAddress("L1_DoubleEG6_HTT250er", &L1_DoubleEG6_HTT250er, &b_L1_DoubleEG6_HTT250er);
   fChain->SetBranchAddress("L1_DoubleEG6_HTT255er", &L1_DoubleEG6_HTT255er, &b_L1_DoubleEG6_HTT255er);
   fChain->SetBranchAddress("L1_DoubleEG6_HTT270er", &L1_DoubleEG6_HTT270er, &b_L1_DoubleEG6_HTT270er);
   fChain->SetBranchAddress("L1_DoubleEG6_HTT300er", &L1_DoubleEG6_HTT300er, &b_L1_DoubleEG6_HTT300er);
   fChain->SetBranchAddress("L1_DoubleEG8er2p6_HTT255er", &L1_DoubleEG8er2p6_HTT255er, &b_L1_DoubleEG8er2p6_HTT255er);
   fChain->SetBranchAddress("L1_DoubleEG8er2p6_HTT270er", &L1_DoubleEG8er2p6_HTT270er, &b_L1_DoubleEG8er2p6_HTT270er);
   fChain->SetBranchAddress("L1_DoubleEG8er2p6_HTT300er", &L1_DoubleEG8er2p6_HTT300er, &b_L1_DoubleEG8er2p6_HTT300er);
   fChain->SetBranchAddress("L1_DoubleEG_15_10", &L1_DoubleEG_15_10, &b_L1_DoubleEG_15_10);
   fChain->SetBranchAddress("L1_DoubleEG_18_17", &L1_DoubleEG_18_17, &b_L1_DoubleEG_18_17);
   fChain->SetBranchAddress("L1_DoubleEG_20_18", &L1_DoubleEG_20_18, &b_L1_DoubleEG_20_18);
   fChain->SetBranchAddress("L1_DoubleEG_22_10", &L1_DoubleEG_22_10, &b_L1_DoubleEG_22_10);
   fChain->SetBranchAddress("L1_DoubleEG_22_12", &L1_DoubleEG_22_12, &b_L1_DoubleEG_22_12);
   fChain->SetBranchAddress("L1_DoubleEG_22_15", &L1_DoubleEG_22_15, &b_L1_DoubleEG_22_15);
   fChain->SetBranchAddress("L1_DoubleEG_23_10", &L1_DoubleEG_23_10, &b_L1_DoubleEG_23_10);
   fChain->SetBranchAddress("L1_DoubleEG_24_17", &L1_DoubleEG_24_17, &b_L1_DoubleEG_24_17);
   fChain->SetBranchAddress("L1_DoubleEG_25_12", &L1_DoubleEG_25_12, &b_L1_DoubleEG_25_12);
   fChain->SetBranchAddress("L1_DoubleEG_25_13", &L1_DoubleEG_25_13, &b_L1_DoubleEG_25_13);
   fChain->SetBranchAddress("L1_DoubleEG_25_14", &L1_DoubleEG_25_14, &b_L1_DoubleEG_25_14);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso23_10", &L1_DoubleEG_LooseIso23_10, &b_L1_DoubleEG_LooseIso23_10);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso24_10", &L1_DoubleEG_LooseIso24_10, &b_L1_DoubleEG_LooseIso24_10);
   fChain->SetBranchAddress("L1_DoubleIsoTau28er2p1", &L1_DoubleIsoTau28er2p1, &b_L1_DoubleIsoTau28er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau30er2p1", &L1_DoubleIsoTau30er2p1, &b_L1_DoubleIsoTau30er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau32er2p1", &L1_DoubleIsoTau32er2p1, &b_L1_DoubleIsoTau32er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau33er2p1", &L1_DoubleIsoTau33er2p1, &b_L1_DoubleIsoTau33er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau34er2p1", &L1_DoubleIsoTau34er2p1, &b_L1_DoubleIsoTau34er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau35er2p1", &L1_DoubleIsoTau35er2p1, &b_L1_DoubleIsoTau35er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau36er2p1", &L1_DoubleIsoTau36er2p1, &b_L1_DoubleIsoTau36er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau38er2p1", &L1_DoubleIsoTau38er2p1, &b_L1_DoubleIsoTau38er2p1);
   fChain->SetBranchAddress("L1_DoubleJet100er2p3_dEta_Max1p6", &L1_DoubleJet100er2p3_dEta_Max1p6, &b_L1_DoubleJet100er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("L1_DoubleJet100er2p7", &L1_DoubleJet100er2p7, &b_L1_DoubleJet100er2p7);
   fChain->SetBranchAddress("L1_DoubleJet112er2p3_dEta_Max1p6", &L1_DoubleJet112er2p3_dEta_Max1p6, &b_L1_DoubleJet112er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("L1_DoubleJet112er2p7", &L1_DoubleJet112er2p7, &b_L1_DoubleJet112er2p7);
   fChain->SetBranchAddress("L1_DoubleJet120er2p7", &L1_DoubleJet120er2p7, &b_L1_DoubleJet120er2p7);
   fChain->SetBranchAddress("L1_DoubleJet150er2p7", &L1_DoubleJet150er2p7, &b_L1_DoubleJet150er2p7);
   fChain->SetBranchAddress("L1_DoubleJet30_Mass_Min300_dEta_Max1p5", &L1_DoubleJet30_Mass_Min300_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min300_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet30_Mass_Min320_dEta_Max1p5", &L1_DoubleJet30_Mass_Min320_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min320_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet30_Mass_Min340_dEta_Max1p5", &L1_DoubleJet30_Mass_Min340_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min340_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet30_Mass_Min360_dEta_Max1p5", &L1_DoubleJet30_Mass_Min360_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min360_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet30_Mass_Min380_dEta_Max1p5", &L1_DoubleJet30_Mass_Min380_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min380_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet30_Mass_Min400_Mu10", &L1_DoubleJet30_Mass_Min400_Mu10, &b_L1_DoubleJet30_Mass_Min400_Mu10);
   fChain->SetBranchAddress("L1_DoubleJet30_Mass_Min400_Mu6", &L1_DoubleJet30_Mass_Min400_Mu6, &b_L1_DoubleJet30_Mass_Min400_Mu6);
   fChain->SetBranchAddress("L1_DoubleJet30_Mass_Min400_dEta_Max1p5", &L1_DoubleJet30_Mass_Min400_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min400_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450", &L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450, &b_L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450);
   fChain->SetBranchAddress("L1_DoubleJet40er2p7", &L1_DoubleJet40er2p7, &b_L1_DoubleJet40er2p7);
   fChain->SetBranchAddress("L1_DoubleJet50er2p7", &L1_DoubleJet50er2p7, &b_L1_DoubleJet50er2p7);
   fChain->SetBranchAddress("L1_DoubleJet60er2p7", &L1_DoubleJet60er2p7, &b_L1_DoubleJet60er2p7);
   fChain->SetBranchAddress("L1_DoubleJet60er2p7_ETM100", &L1_DoubleJet60er2p7_ETM100, &b_L1_DoubleJet60er2p7_ETM100);
   fChain->SetBranchAddress("L1_DoubleJet60er2p7_ETM60", &L1_DoubleJet60er2p7_ETM60, &b_L1_DoubleJet60er2p7_ETM60);
   fChain->SetBranchAddress("L1_DoubleJet60er2p7_ETM70", &L1_DoubleJet60er2p7_ETM70, &b_L1_DoubleJet60er2p7_ETM70);
   fChain->SetBranchAddress("L1_DoubleJet60er2p7_ETM80", &L1_DoubleJet60er2p7_ETM80, &b_L1_DoubleJet60er2p7_ETM80);
   fChain->SetBranchAddress("L1_DoubleJet60er2p7_ETM90", &L1_DoubleJet60er2p7_ETM90, &b_L1_DoubleJet60er2p7_ETM90);
   fChain->SetBranchAddress("L1_DoubleJet80er2p7", &L1_DoubleJet80er2p7, &b_L1_DoubleJet80er2p7);
   fChain->SetBranchAddress("L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_100_30_DoubleJet30_Mass_Min620, &b_L1_DoubleJet_100_30_DoubleJet30_Mass_Min620);
   fChain->SetBranchAddress("L1_DoubleJet_100_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_100_35_DoubleJet35_Mass_Min620, &b_L1_DoubleJet_100_35_DoubleJet35_Mass_Min620);
   fChain->SetBranchAddress("L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min620, &b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min620);
   fChain->SetBranchAddress("L1_DoubleJet_110_40_DoubleJet40_Mass_Min620", &L1_DoubleJet_110_40_DoubleJet40_Mass_Min620, &b_L1_DoubleJet_110_40_DoubleJet40_Mass_Min620);
   fChain->SetBranchAddress("L1_DoubleJet_115_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_115_35_DoubleJet35_Mass_Min620, &b_L1_DoubleJet_115_35_DoubleJet35_Mass_Min620);
   fChain->SetBranchAddress("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620, &b_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620);
   fChain->SetBranchAddress("L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_90_30_DoubleJet30_Mass_Min620, &b_L1_DoubleJet_90_30_DoubleJet30_Mass_Min620);
   fChain->SetBranchAddress("L1_DoubleLooseIsoEG22er2p1", &L1_DoubleLooseIsoEG22er2p1, &b_L1_DoubleLooseIsoEG22er2p1);
   fChain->SetBranchAddress("L1_DoubleLooseIsoEG24er2p1", &L1_DoubleLooseIsoEG24er2p1, &b_L1_DoubleLooseIsoEG24er2p1);
   fChain->SetBranchAddress("L1_DoubleMu0", &L1_DoubleMu0, &b_L1_DoubleMu0);
   fChain->SetBranchAddress("L1_DoubleMu0_ETM40", &L1_DoubleMu0_ETM40, &b_L1_DoubleMu0_ETM40);
   fChain->SetBranchAddress("L1_DoubleMu0_ETM55", &L1_DoubleMu0_ETM55, &b_L1_DoubleMu0_ETM55);
   fChain->SetBranchAddress("L1_DoubleMu0_ETM60", &L1_DoubleMu0_ETM60, &b_L1_DoubleMu0_ETM60);
   fChain->SetBranchAddress("L1_DoubleMu0_ETM65", &L1_DoubleMu0_ETM65, &b_L1_DoubleMu0_ETM65);
   fChain->SetBranchAddress("L1_DoubleMu0_ETM70", &L1_DoubleMu0_ETM70, &b_L1_DoubleMu0_ETM70);
   fChain->SetBranchAddress("L1_DoubleMu0_SQ", &L1_DoubleMu0_SQ, &b_L1_DoubleMu0_SQ);
   fChain->SetBranchAddress("L1_DoubleMu0_SQ_OS", &L1_DoubleMu0_SQ_OS, &b_L1_DoubleMu0_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4, &b_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4);
   fChain->SetBranchAddress("L1_DoubleMu0er1p4_dEta_Max1p8_OS", &L1_DoubleMu0er1p4_dEta_Max1p8_OS, &b_L1_DoubleMu0er1p4_dEta_Max1p8_OS);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS", &L1_DoubleMu0er1p5_SQ_OS, &b_L1_DoubleMu0er1p5_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4, &b_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_dR_Max1p4, &b_L1_DoubleMu0er1p5_SQ_dR_Max1p4);
   fChain->SetBranchAddress("L1_DoubleMu0er2_SQ_dR_Max1p4", &L1_DoubleMu0er2_SQ_dR_Max1p4, &b_L1_DoubleMu0er2_SQ_dR_Max1p4);
   fChain->SetBranchAddress("L1_DoubleMu18er2p1", &L1_DoubleMu18er2p1, &b_L1_DoubleMu18er2p1);
   fChain->SetBranchAddress("L1_DoubleMu22er2p1", &L1_DoubleMu22er2p1, &b_L1_DoubleMu22er2p1);
   fChain->SetBranchAddress("L1_DoubleMu3_OS_DoubleEG7p5Upsilon", &L1_DoubleMu3_OS_DoubleEG7p5Upsilon, &b_L1_DoubleMu3_OS_DoubleEG7p5Upsilon);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_HTT100er", &L1_DoubleMu3_SQ_HTT100er, &b_L1_DoubleMu3_SQ_HTT100er);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_HTT200er", &L1_DoubleMu3_SQ_HTT200er, &b_L1_DoubleMu3_SQ_HTT200er);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_HTT220er", &L1_DoubleMu3_SQ_HTT220er, &b_L1_DoubleMu3_SQ_HTT220er);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_HTT240er", &L1_DoubleMu3_SQ_HTT240er, &b_L1_DoubleMu3_SQ_HTT240er);
   fChain->SetBranchAddress("L1_DoubleMu4_OS_EG12", &L1_DoubleMu4_OS_EG12, &b_L1_DoubleMu4_OS_EG12);
   fChain->SetBranchAddress("L1_DoubleMu4_SQ_OS", &L1_DoubleMu4_SQ_OS, &b_L1_DoubleMu4_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu4_SQ_OS_dR_Max1p2", &L1_DoubleMu4_SQ_OS_dR_Max1p2, &b_L1_DoubleMu4_SQ_OS_dR_Max1p2);
   fChain->SetBranchAddress("L1_DoubleMu4p5_SQ", &L1_DoubleMu4p5_SQ, &b_L1_DoubleMu4p5_SQ);
   fChain->SetBranchAddress("L1_DoubleMu4p5_SQ_OS", &L1_DoubleMu4p5_SQ_OS, &b_L1_DoubleMu4p5_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2, &b_L1_DoubleMu4p5_SQ_OS_dR_Max1p2);
   fChain->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS", &L1_DoubleMu4p5er2p0_SQ_OS, &b_L1_DoubleMu4p5er2p0_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18", &L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18, &b_L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18);
   fChain->SetBranchAddress("L1_DoubleMu5Upsilon_OS_DoubleEG3", &L1_DoubleMu5Upsilon_OS_DoubleEG3, &b_L1_DoubleMu5Upsilon_OS_DoubleEG3);
   fChain->SetBranchAddress("L1_DoubleMu5_OS_EG12", &L1_DoubleMu5_OS_EG12, &b_L1_DoubleMu5_OS_EG12);
   fChain->SetBranchAddress("L1_DoubleMu5_SQ_OS", &L1_DoubleMu5_SQ_OS, &b_L1_DoubleMu5_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu5_SQ_OS_Mass7to18", &L1_DoubleMu5_SQ_OS_Mass7to18, &b_L1_DoubleMu5_SQ_OS_Mass7to18);
   fChain->SetBranchAddress("L1_DoubleMu6_SQ_OS", &L1_DoubleMu6_SQ_OS, &b_L1_DoubleMu6_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu7_EG7", &L1_DoubleMu7_EG7, &b_L1_DoubleMu7_EG7);
   fChain->SetBranchAddress("L1_DoubleMu7_SQ_EG7", &L1_DoubleMu7_SQ_EG7, &b_L1_DoubleMu7_SQ_EG7);
   fChain->SetBranchAddress("L1_DoubleMu8_SQ", &L1_DoubleMu8_SQ, &b_L1_DoubleMu8_SQ);
   fChain->SetBranchAddress("L1_DoubleMu_10_0_dEta_Max1p8", &L1_DoubleMu_10_0_dEta_Max1p8, &b_L1_DoubleMu_10_0_dEta_Max1p8);
   fChain->SetBranchAddress("L1_DoubleMu_11_4", &L1_DoubleMu_11_4, &b_L1_DoubleMu_11_4);
   fChain->SetBranchAddress("L1_DoubleMu_12_5", &L1_DoubleMu_12_5, &b_L1_DoubleMu_12_5);
   fChain->SetBranchAddress("L1_DoubleMu_12_8", &L1_DoubleMu_12_8, &b_L1_DoubleMu_12_8);
   fChain->SetBranchAddress("L1_DoubleMu_13_6", &L1_DoubleMu_13_6, &b_L1_DoubleMu_13_6);
   fChain->SetBranchAddress("L1_DoubleMu_15_5", &L1_DoubleMu_15_5, &b_L1_DoubleMu_15_5);
   fChain->SetBranchAddress("L1_DoubleMu_15_5_SQ", &L1_DoubleMu_15_5_SQ, &b_L1_DoubleMu_15_5_SQ);
   fChain->SetBranchAddress("L1_DoubleMu_15_7", &L1_DoubleMu_15_7, &b_L1_DoubleMu_15_7);
   fChain->SetBranchAddress("L1_DoubleMu_15_7_SQ", &L1_DoubleMu_15_7_SQ, &b_L1_DoubleMu_15_7_SQ);
   fChain->SetBranchAddress("L1_DoubleMu_15_7_SQ_Mass_Min4", &L1_DoubleMu_15_7_SQ_Mass_Min4, &b_L1_DoubleMu_15_7_SQ_Mass_Min4);
   fChain->SetBranchAddress("L1_DoubleMu_20_2_SQ_Mass_Max20", &L1_DoubleMu_20_2_SQ_Mass_Max20, &b_L1_DoubleMu_20_2_SQ_Mass_Max20);
   fChain->SetBranchAddress("L1_DoubleTau50er2p1", &L1_DoubleTau50er2p1, &b_L1_DoubleTau50er2p1);
   fChain->SetBranchAddress("L1_DoubleTau70er2p1", &L1_DoubleTau70er2p1, &b_L1_DoubleTau70er2p1);
   fChain->SetBranchAddress("L1_EG25er2p1_HTT125er", &L1_EG25er2p1_HTT125er, &b_L1_EG25er2p1_HTT125er);
   fChain->SetBranchAddress("L1_EG27er2p1_HTT200er", &L1_EG27er2p1_HTT200er, &b_L1_EG27er2p1_HTT200er);
   fChain->SetBranchAddress("L1_ETM100", &L1_ETM100, &b_L1_ETM100);
   fChain->SetBranchAddress("L1_ETM100_Jet60_dPhi_Min0p4", &L1_ETM100_Jet60_dPhi_Min0p4, &b_L1_ETM100_Jet60_dPhi_Min0p4);
   fChain->SetBranchAddress("L1_ETM105", &L1_ETM105, &b_L1_ETM105);
   fChain->SetBranchAddress("L1_ETM110", &L1_ETM110, &b_L1_ETM110);
   fChain->SetBranchAddress("L1_ETM110_Jet60_dPhi_Min0p4", &L1_ETM110_Jet60_dPhi_Min0p4, &b_L1_ETM110_Jet60_dPhi_Min0p4);
   fChain->SetBranchAddress("L1_ETM115", &L1_ETM115, &b_L1_ETM115);
   fChain->SetBranchAddress("L1_ETM120", &L1_ETM120, &b_L1_ETM120);
   fChain->SetBranchAddress("L1_ETM150", &L1_ETM150, &b_L1_ETM150);
   fChain->SetBranchAddress("L1_ETM30", &L1_ETM30, &b_L1_ETM30);
   fChain->SetBranchAddress("L1_ETM40", &L1_ETM40, &b_L1_ETM40);
   fChain->SetBranchAddress("L1_ETM50", &L1_ETM50, &b_L1_ETM50);
   fChain->SetBranchAddress("L1_ETM60", &L1_ETM60, &b_L1_ETM60);
   fChain->SetBranchAddress("L1_ETM70", &L1_ETM70, &b_L1_ETM70);
   fChain->SetBranchAddress("L1_ETM75", &L1_ETM75, &b_L1_ETM75);
   fChain->SetBranchAddress("L1_ETM75_Jet60_dPhi_Min0p4", &L1_ETM75_Jet60_dPhi_Min0p4, &b_L1_ETM75_Jet60_dPhi_Min0p4);
   fChain->SetBranchAddress("L1_ETM80", &L1_ETM80, &b_L1_ETM80);
   fChain->SetBranchAddress("L1_ETM80_Jet60_dPhi_Min0p4", &L1_ETM80_Jet60_dPhi_Min0p4, &b_L1_ETM80_Jet60_dPhi_Min0p4);
   fChain->SetBranchAddress("L1_ETM85", &L1_ETM85, &b_L1_ETM85);
   fChain->SetBranchAddress("L1_ETM90", &L1_ETM90, &b_L1_ETM90);
   fChain->SetBranchAddress("L1_ETM90_Jet60_dPhi_Min0p4", &L1_ETM90_Jet60_dPhi_Min0p4, &b_L1_ETM90_Jet60_dPhi_Min0p4);
   fChain->SetBranchAddress("L1_ETM95", &L1_ETM95, &b_L1_ETM95);
   fChain->SetBranchAddress("L1_ETMHF100", &L1_ETMHF100, &b_L1_ETMHF100);
   fChain->SetBranchAddress("L1_ETMHF100_HTT60er", &L1_ETMHF100_HTT60er, &b_L1_ETMHF100_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF100_Jet60_OR_DiJet30woTT28", &L1_ETMHF100_Jet60_OR_DiJet30woTT28, &b_L1_ETMHF100_Jet60_OR_DiJet30woTT28);
   fChain->SetBranchAddress("L1_ETMHF100_Jet60_OR_DoubleJet30", &L1_ETMHF100_Jet60_OR_DoubleJet30, &b_L1_ETMHF100_Jet60_OR_DoubleJet30);
   fChain->SetBranchAddress("L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30);
   fChain->SetBranchAddress("L1_ETMHF110", &L1_ETMHF110, &b_L1_ETMHF110);
   fChain->SetBranchAddress("L1_ETMHF110_HTT60er", &L1_ETMHF110_HTT60er, &b_L1_ETMHF110_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF110_Jet60_OR_DiJet30woTT28", &L1_ETMHF110_Jet60_OR_DiJet30woTT28, &b_L1_ETMHF110_Jet60_OR_DiJet30woTT28);
   fChain->SetBranchAddress("L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30);
   fChain->SetBranchAddress("L1_ETMHF120", &L1_ETMHF120, &b_L1_ETMHF120);
   fChain->SetBranchAddress("L1_ETMHF120_HTT60er", &L1_ETMHF120_HTT60er, &b_L1_ETMHF120_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF120_Jet60_OR_DiJet30woTT28", &L1_ETMHF120_Jet60_OR_DiJet30woTT28, &b_L1_ETMHF120_Jet60_OR_DiJet30woTT28);
   fChain->SetBranchAddress("L1_ETMHF150", &L1_ETMHF150, &b_L1_ETMHF150);
   fChain->SetBranchAddress("L1_ETMHF70", &L1_ETMHF70, &b_L1_ETMHF70);
   fChain->SetBranchAddress("L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30);
   fChain->SetBranchAddress("L1_ETMHF80", &L1_ETMHF80, &b_L1_ETMHF80);
   fChain->SetBranchAddress("L1_ETMHF80_HTT60er", &L1_ETMHF80_HTT60er, &b_L1_ETMHF80_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30);
   fChain->SetBranchAddress("L1_ETMHF90", &L1_ETMHF90, &b_L1_ETMHF90);
   fChain->SetBranchAddress("L1_ETMHF90_HTT60er", &L1_ETMHF90_HTT60er, &b_L1_ETMHF90_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30);
   fChain->SetBranchAddress("L1_ETT100_BptxAND", &L1_ETT100_BptxAND, &b_L1_ETT100_BptxAND);
   fChain->SetBranchAddress("L1_ETT110_BptxAND", &L1_ETT110_BptxAND, &b_L1_ETT110_BptxAND);
   fChain->SetBranchAddress("L1_ETT40_BptxAND", &L1_ETT40_BptxAND, &b_L1_ETT40_BptxAND);
   fChain->SetBranchAddress("L1_ETT50_BptxAND", &L1_ETT50_BptxAND, &b_L1_ETT50_BptxAND);
   fChain->SetBranchAddress("L1_ETT60_BptxAND", &L1_ETT60_BptxAND, &b_L1_ETT60_BptxAND);
   fChain->SetBranchAddress("L1_ETT70_BptxAND", &L1_ETT70_BptxAND, &b_L1_ETT70_BptxAND);
   fChain->SetBranchAddress("L1_ETT75_BptxAND", &L1_ETT75_BptxAND, &b_L1_ETT75_BptxAND);
   fChain->SetBranchAddress("L1_ETT80_BptxAND", &L1_ETT80_BptxAND, &b_L1_ETT80_BptxAND);
   fChain->SetBranchAddress("L1_ETT85_BptxAND", &L1_ETT85_BptxAND, &b_L1_ETT85_BptxAND);
   fChain->SetBranchAddress("L1_ETT90_BptxAND", &L1_ETT90_BptxAND, &b_L1_ETT90_BptxAND);
   fChain->SetBranchAddress("L1_ETT95_BptxAND", &L1_ETT95_BptxAND, &b_L1_ETT95_BptxAND);
   fChain->SetBranchAddress("L1_FirstBunchAfterTrain", &L1_FirstBunchAfterTrain, &b_L1_FirstBunchAfterTrain);
   fChain->SetBranchAddress("L1_FirstBunchInTrain", &L1_FirstBunchInTrain, &b_L1_FirstBunchInTrain);
   fChain->SetBranchAddress("L1_FirstCollisionInOrbit", &L1_FirstCollisionInOrbit, &b_L1_FirstCollisionInOrbit);
   fChain->SetBranchAddress("L1_FirstCollisionInTrain", &L1_FirstCollisionInTrain, &b_L1_FirstCollisionInTrain);
   fChain->SetBranchAddress("L1_HTT120er", &L1_HTT120er, &b_L1_HTT120er);
   fChain->SetBranchAddress("L1_HTT160er", &L1_HTT160er, &b_L1_HTT160er);
   fChain->SetBranchAddress("L1_HTT200er", &L1_HTT200er, &b_L1_HTT200er);
   fChain->SetBranchAddress("L1_HTT220er", &L1_HTT220er, &b_L1_HTT220er);
   fChain->SetBranchAddress("L1_HTT240er", &L1_HTT240er, &b_L1_HTT240er);
   fChain->SetBranchAddress("L1_HTT250er_QuadJet_70_55_40_35_er2p5", &L1_HTT250er_QuadJet_70_55_40_35_er2p5, &b_L1_HTT250er_QuadJet_70_55_40_35_er2p5);
   fChain->SetBranchAddress("L1_HTT255er", &L1_HTT255er, &b_L1_HTT255er);
   fChain->SetBranchAddress("L1_HTT270er", &L1_HTT270er, &b_L1_HTT270er);
   fChain->SetBranchAddress("L1_HTT280er", &L1_HTT280er, &b_L1_HTT280er);
   fChain->SetBranchAddress("L1_HTT280er_QuadJet_70_55_40_35_er2p5", &L1_HTT280er_QuadJet_70_55_40_35_er2p5, &b_L1_HTT280er_QuadJet_70_55_40_35_er2p5);
   fChain->SetBranchAddress("L1_HTT300er", &L1_HTT300er, &b_L1_HTT300er);
   fChain->SetBranchAddress("L1_HTT300er_QuadJet_70_55_40_35_er2p5", &L1_HTT300er_QuadJet_70_55_40_35_er2p5, &b_L1_HTT300er_QuadJet_70_55_40_35_er2p5);
   fChain->SetBranchAddress("L1_HTT320er", &L1_HTT320er, &b_L1_HTT320er);
   fChain->SetBranchAddress("L1_HTT320er_QuadJet_70_55_40_40_er2p4", &L1_HTT320er_QuadJet_70_55_40_40_er2p4, &b_L1_HTT320er_QuadJet_70_55_40_40_er2p4);
   fChain->SetBranchAddress("L1_HTT320er_QuadJet_70_55_40_40_er2p5", &L1_HTT320er_QuadJet_70_55_40_40_er2p5, &b_L1_HTT320er_QuadJet_70_55_40_40_er2p5);
   fChain->SetBranchAddress("L1_HTT320er_QuadJet_70_55_45_45_er2p5", &L1_HTT320er_QuadJet_70_55_45_45_er2p5, &b_L1_HTT320er_QuadJet_70_55_45_45_er2p5);
   fChain->SetBranchAddress("L1_HTT340er", &L1_HTT340er, &b_L1_HTT340er);
   fChain->SetBranchAddress("L1_HTT340er_QuadJet_70_55_40_40_er2p5", &L1_HTT340er_QuadJet_70_55_40_40_er2p5, &b_L1_HTT340er_QuadJet_70_55_40_40_er2p5);
   fChain->SetBranchAddress("L1_HTT340er_QuadJet_70_55_45_45_er2p5", &L1_HTT340er_QuadJet_70_55_45_45_er2p5, &b_L1_HTT340er_QuadJet_70_55_45_45_er2p5);
   fChain->SetBranchAddress("L1_HTT380er", &L1_HTT380er, &b_L1_HTT380er);
   fChain->SetBranchAddress("L1_HTT400er", &L1_HTT400er, &b_L1_HTT400er);
   fChain->SetBranchAddress("L1_HTT450er", &L1_HTT450er, &b_L1_HTT450er);
   fChain->SetBranchAddress("L1_HTT500er", &L1_HTT500er, &b_L1_HTT500er);
   fChain->SetBranchAddress("L1_IsoEG33_Mt40", &L1_IsoEG33_Mt40, &b_L1_IsoEG33_Mt40);
   fChain->SetBranchAddress("L1_IsoEG33_Mt44", &L1_IsoEG33_Mt44, &b_L1_IsoEG33_Mt44);
   fChain->SetBranchAddress("L1_IsoEG33_Mt48", &L1_IsoEG33_Mt48, &b_L1_IsoEG33_Mt48);
   fChain->SetBranchAddress("L1_IsoTau40er_ETM100", &L1_IsoTau40er_ETM100, &b_L1_IsoTau40er_ETM100);
   fChain->SetBranchAddress("L1_IsoTau40er_ETM105", &L1_IsoTau40er_ETM105, &b_L1_IsoTau40er_ETM105);
   fChain->SetBranchAddress("L1_IsoTau40er_ETM110", &L1_IsoTau40er_ETM110, &b_L1_IsoTau40er_ETM110);
   fChain->SetBranchAddress("L1_IsoTau40er_ETM115", &L1_IsoTau40er_ETM115, &b_L1_IsoTau40er_ETM115);
   fChain->SetBranchAddress("L1_IsoTau40er_ETM120", &L1_IsoTau40er_ETM120, &b_L1_IsoTau40er_ETM120);
   fChain->SetBranchAddress("L1_IsoTau40er_ETM80", &L1_IsoTau40er_ETM80, &b_L1_IsoTau40er_ETM80);
   fChain->SetBranchAddress("L1_IsoTau40er_ETM85", &L1_IsoTau40er_ETM85, &b_L1_IsoTau40er_ETM85);
   fChain->SetBranchAddress("L1_IsoTau40er_ETM90", &L1_IsoTau40er_ETM90, &b_L1_IsoTau40er_ETM90);
   fChain->SetBranchAddress("L1_IsoTau40er_ETM95", &L1_IsoTau40er_ETM95, &b_L1_IsoTau40er_ETM95);
   fChain->SetBranchAddress("L1_IsoTau40er_ETMHF100", &L1_IsoTau40er_ETMHF100, &b_L1_IsoTau40er_ETMHF100);
   fChain->SetBranchAddress("L1_IsoTau40er_ETMHF110", &L1_IsoTau40er_ETMHF110, &b_L1_IsoTau40er_ETMHF110);
   fChain->SetBranchAddress("L1_IsoTau40er_ETMHF120", &L1_IsoTau40er_ETMHF120, &b_L1_IsoTau40er_ETMHF120);
   fChain->SetBranchAddress("L1_IsoTau40er_ETMHF80", &L1_IsoTau40er_ETMHF80, &b_L1_IsoTau40er_ETMHF80);
   fChain->SetBranchAddress("L1_IsoTau40er_ETMHF90", &L1_IsoTau40er_ETMHF90, &b_L1_IsoTau40er_ETMHF90);
   fChain->SetBranchAddress("L1_IsolatedBunch", &L1_IsolatedBunch, &b_L1_IsolatedBunch);
   fChain->SetBranchAddress("L1_LastCollisionInTrain", &L1_LastCollisionInTrain, &b_L1_LastCollisionInTrain);
   fChain->SetBranchAddress("L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3, &b_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG24er2p1_HTT100er", &L1_LooseIsoEG24er2p1_HTT100er, &b_L1_LooseIsoEG24er2p1_HTT100er);
   fChain->SetBranchAddress("L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3, &b_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3", &L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3, &b_L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7", &L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7, &b_L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7);
   fChain->SetBranchAddress("L1_LooseIsoEG26er2p1_HTT100er", &L1_LooseIsoEG26er2p1_HTT100er, &b_L1_LooseIsoEG26er2p1_HTT100er);
   fChain->SetBranchAddress("L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3", &L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3, &b_L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG28er2p1_HTT100er", &L1_LooseIsoEG28er2p1_HTT100er, &b_L1_LooseIsoEG28er2p1_HTT100er);
   fChain->SetBranchAddress("L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3", &L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3, &b_L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3", &L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3, &b_L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3);
   fChain->SetBranchAddress("L1_MU20_EG15", &L1_MU20_EG15, &b_L1_MU20_EG15);
   fChain->SetBranchAddress("L1_MinimumBiasHF0_AND_BptxAND", &L1_MinimumBiasHF0_AND_BptxAND, &b_L1_MinimumBiasHF0_AND_BptxAND);
   fChain->SetBranchAddress("L1_MinimumBiasHF0_OR_BptxAND", &L1_MinimumBiasHF0_OR_BptxAND, &b_L1_MinimumBiasHF0_OR_BptxAND);
   fChain->SetBranchAddress("L1_Mu10er2p1_ETM30", &L1_Mu10er2p1_ETM30, &b_L1_Mu10er2p1_ETM30);
   fChain->SetBranchAddress("L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6, &b_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("L1_Mu12_EG10", &L1_Mu12_EG10, &b_L1_Mu12_EG10);
   fChain->SetBranchAddress("L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6, &b_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("L1_Mu14er2p1_ETM30", &L1_Mu14er2p1_ETM30, &b_L1_Mu14er2p1_ETM30);
   fChain->SetBranchAddress("L1_Mu15_HTT100er", &L1_Mu15_HTT100er, &b_L1_Mu15_HTT100er);
   fChain->SetBranchAddress("L1_Mu18_HTT100er", &L1_Mu18_HTT100er, &b_L1_Mu18_HTT100er);
   fChain->SetBranchAddress("L1_Mu18_Jet24er2p7", &L1_Mu18_Jet24er2p7, &b_L1_Mu18_Jet24er2p7);
   fChain->SetBranchAddress("L1_Mu18er2p1_IsoTau26er2p1", &L1_Mu18er2p1_IsoTau26er2p1, &b_L1_Mu18er2p1_IsoTau26er2p1);
   fChain->SetBranchAddress("L1_Mu18er2p1_Tau24er2p1", &L1_Mu18er2p1_Tau24er2p1, &b_L1_Mu18er2p1_Tau24er2p1);
   fChain->SetBranchAddress("L1_Mu20_EG10", &L1_Mu20_EG10, &b_L1_Mu20_EG10);
   fChain->SetBranchAddress("L1_Mu20_EG17", &L1_Mu20_EG17, &b_L1_Mu20_EG17);
   fChain->SetBranchAddress("L1_Mu20_LooseIsoEG6", &L1_Mu20_LooseIsoEG6, &b_L1_Mu20_LooseIsoEG6);
   fChain->SetBranchAddress("L1_Mu20er2p1_IsoTau26er2p1", &L1_Mu20er2p1_IsoTau26er2p1, &b_L1_Mu20er2p1_IsoTau26er2p1);
   fChain->SetBranchAddress("L1_Mu20er2p1_IsoTau27er2p1", &L1_Mu20er2p1_IsoTau27er2p1, &b_L1_Mu20er2p1_IsoTau27er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau28er2p1", &L1_Mu22er2p1_IsoTau28er2p1, &b_L1_Mu22er2p1_IsoTau28er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau30er2p1", &L1_Mu22er2p1_IsoTau30er2p1, &b_L1_Mu22er2p1_IsoTau30er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau32er2p1", &L1_Mu22er2p1_IsoTau32er2p1, &b_L1_Mu22er2p1_IsoTau32er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau33er2p1", &L1_Mu22er2p1_IsoTau33er2p1, &b_L1_Mu22er2p1_IsoTau33er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau34er2p1", &L1_Mu22er2p1_IsoTau34er2p1, &b_L1_Mu22er2p1_IsoTau34er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau35er2p1", &L1_Mu22er2p1_IsoTau35er2p1, &b_L1_Mu22er2p1_IsoTau35er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau36er2p1", &L1_Mu22er2p1_IsoTau36er2p1, &b_L1_Mu22er2p1_IsoTau36er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau38er2p1", &L1_Mu22er2p1_IsoTau38er2p1, &b_L1_Mu22er2p1_IsoTau38er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau40er2p1", &L1_Mu22er2p1_IsoTau40er2p1, &b_L1_Mu22er2p1_IsoTau40er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_Tau50er2p1", &L1_Mu22er2p1_Tau50er2p1, &b_L1_Mu22er2p1_Tau50er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_Tau70er2p1", &L1_Mu22er2p1_Tau70er2p1, &b_L1_Mu22er2p1_Tau70er2p1);
   fChain->SetBranchAddress("L1_Mu23_EG10", &L1_Mu23_EG10, &b_L1_Mu23_EG10);
   fChain->SetBranchAddress("L1_Mu23_LooseIsoEG10", &L1_Mu23_LooseIsoEG10, &b_L1_Mu23_LooseIsoEG10);
   fChain->SetBranchAddress("L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4", &L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4, &b_L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4);
   fChain->SetBranchAddress("L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4", &L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4, &b_L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4);
   fChain->SetBranchAddress("L1_Mu3_Jet30er2p5", &L1_Mu3_Jet30er2p5, &b_L1_Mu3_Jet30er2p5);
   fChain->SetBranchAddress("L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4", &L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4, &b_L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4);
   fChain->SetBranchAddress("L1_Mu5_EG15", &L1_Mu5_EG15, &b_L1_Mu5_EG15);
   fChain->SetBranchAddress("L1_Mu5_EG20", &L1_Mu5_EG20, &b_L1_Mu5_EG20);
   fChain->SetBranchAddress("L1_Mu5_EG23", &L1_Mu5_EG23, &b_L1_Mu5_EG23);
   fChain->SetBranchAddress("L1_Mu5_LooseIsoEG18", &L1_Mu5_LooseIsoEG18, &b_L1_Mu5_LooseIsoEG18);
   fChain->SetBranchAddress("L1_Mu5_LooseIsoEG20", &L1_Mu5_LooseIsoEG20, &b_L1_Mu5_LooseIsoEG20);
   fChain->SetBranchAddress("L1_Mu6_DoubleEG10", &L1_Mu6_DoubleEG10, &b_L1_Mu6_DoubleEG10);
   fChain->SetBranchAddress("L1_Mu6_DoubleEG17", &L1_Mu6_DoubleEG17, &b_L1_Mu6_DoubleEG17);
   fChain->SetBranchAddress("L1_Mu6_HTT200er", &L1_Mu6_HTT200er, &b_L1_Mu6_HTT200er);
   fChain->SetBranchAddress("L1_Mu6_HTT240er", &L1_Mu6_HTT240er, &b_L1_Mu6_HTT240er);
   fChain->SetBranchAddress("L1_Mu6_HTT250er", &L1_Mu6_HTT250er, &b_L1_Mu6_HTT250er);
   fChain->SetBranchAddress("L1_Mu7_EG23", &L1_Mu7_EG23, &b_L1_Mu7_EG23);
   fChain->SetBranchAddress("L1_Mu7_LooseIsoEG20", &L1_Mu7_LooseIsoEG20, &b_L1_Mu7_LooseIsoEG20);
   fChain->SetBranchAddress("L1_Mu7_LooseIsoEG23", &L1_Mu7_LooseIsoEG23, &b_L1_Mu7_LooseIsoEG23);
   fChain->SetBranchAddress("L1_Mu8_HTT150er", &L1_Mu8_HTT150er, &b_L1_Mu8_HTT150er);
   fChain->SetBranchAddress("L1_NotBptxOR", &L1_NotBptxOR, &b_L1_NotBptxOR);
   fChain->SetBranchAddress("L1_QuadJet36er2p7_IsoTau52er2p1", &L1_QuadJet36er2p7_IsoTau52er2p1, &b_L1_QuadJet36er2p7_IsoTau52er2p1);
   fChain->SetBranchAddress("L1_QuadJet36er2p7_Tau52", &L1_QuadJet36er2p7_Tau52, &b_L1_QuadJet36er2p7_Tau52);
   fChain->SetBranchAddress("L1_QuadJet40er2p7", &L1_QuadJet40er2p7, &b_L1_QuadJet40er2p7);
   fChain->SetBranchAddress("L1_QuadJet50er2p7", &L1_QuadJet50er2p7, &b_L1_QuadJet50er2p7);
   fChain->SetBranchAddress("L1_QuadJet60er2p7", &L1_QuadJet60er2p7, &b_L1_QuadJet60er2p7);
   fChain->SetBranchAddress("L1_QuadMu0", &L1_QuadMu0, &b_L1_QuadMu0);
   fChain->SetBranchAddress("L1_SingleEG10", &L1_SingleEG10, &b_L1_SingleEG10);
   fChain->SetBranchAddress("L1_SingleEG15", &L1_SingleEG15, &b_L1_SingleEG15);
   fChain->SetBranchAddress("L1_SingleEG18", &L1_SingleEG18, &b_L1_SingleEG18);
   fChain->SetBranchAddress("L1_SingleEG24", &L1_SingleEG24, &b_L1_SingleEG24);
   fChain->SetBranchAddress("L1_SingleEG26", &L1_SingleEG26, &b_L1_SingleEG26);
   fChain->SetBranchAddress("L1_SingleEG28", &L1_SingleEG28, &b_L1_SingleEG28);
   fChain->SetBranchAddress("L1_SingleEG2_BptxAND", &L1_SingleEG2_BptxAND, &b_L1_SingleEG2_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG30", &L1_SingleEG30, &b_L1_SingleEG30);
   fChain->SetBranchAddress("L1_SingleEG32", &L1_SingleEG32, &b_L1_SingleEG32);
   fChain->SetBranchAddress("L1_SingleEG34", &L1_SingleEG34, &b_L1_SingleEG34);
   fChain->SetBranchAddress("L1_SingleEG34er2p1", &L1_SingleEG34er2p1, &b_L1_SingleEG34er2p1);
   fChain->SetBranchAddress("L1_SingleEG36", &L1_SingleEG36, &b_L1_SingleEG36);
   fChain->SetBranchAddress("L1_SingleEG36er2p1", &L1_SingleEG36er2p1, &b_L1_SingleEG36er2p1);
   fChain->SetBranchAddress("L1_SingleEG38", &L1_SingleEG38, &b_L1_SingleEG38);
   fChain->SetBranchAddress("L1_SingleEG38er2p1", &L1_SingleEG38er2p1, &b_L1_SingleEG38er2p1);
   fChain->SetBranchAddress("L1_SingleEG40", &L1_SingleEG40, &b_L1_SingleEG40);
   fChain->SetBranchAddress("L1_SingleEG42", &L1_SingleEG42, &b_L1_SingleEG42);
   fChain->SetBranchAddress("L1_SingleEG45", &L1_SingleEG45, &b_L1_SingleEG45);
   fChain->SetBranchAddress("L1_SingleEG5", &L1_SingleEG5, &b_L1_SingleEG5);
   fChain->SetBranchAddress("L1_SingleEG50", &L1_SingleEG50, &b_L1_SingleEG50);
   fChain->SetBranchAddress("L1_SingleIsoEG18", &L1_SingleIsoEG18, &b_L1_SingleIsoEG18);
   fChain->SetBranchAddress("L1_SingleIsoEG18er2p1", &L1_SingleIsoEG18er2p1, &b_L1_SingleIsoEG18er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG20", &L1_SingleIsoEG20, &b_L1_SingleIsoEG20);
   fChain->SetBranchAddress("L1_SingleIsoEG20er2p1", &L1_SingleIsoEG20er2p1, &b_L1_SingleIsoEG20er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG22", &L1_SingleIsoEG22, &b_L1_SingleIsoEG22);
   fChain->SetBranchAddress("L1_SingleIsoEG22er2p1", &L1_SingleIsoEG22er2p1, &b_L1_SingleIsoEG22er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG24", &L1_SingleIsoEG24, &b_L1_SingleIsoEG24);
   fChain->SetBranchAddress("L1_SingleIsoEG24er2p1", &L1_SingleIsoEG24er2p1, &b_L1_SingleIsoEG24er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG26", &L1_SingleIsoEG26, &b_L1_SingleIsoEG26);
   fChain->SetBranchAddress("L1_SingleIsoEG26er2p1", &L1_SingleIsoEG26er2p1, &b_L1_SingleIsoEG26er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG28", &L1_SingleIsoEG28, &b_L1_SingleIsoEG28);
   fChain->SetBranchAddress("L1_SingleIsoEG28er2p1", &L1_SingleIsoEG28er2p1, &b_L1_SingleIsoEG28er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG30", &L1_SingleIsoEG30, &b_L1_SingleIsoEG30);
   fChain->SetBranchAddress("L1_SingleIsoEG30er2p1", &L1_SingleIsoEG30er2p1, &b_L1_SingleIsoEG30er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG32", &L1_SingleIsoEG32, &b_L1_SingleIsoEG32);
   fChain->SetBranchAddress("L1_SingleIsoEG32er2p1", &L1_SingleIsoEG32er2p1, &b_L1_SingleIsoEG32er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG33er2p1", &L1_SingleIsoEG33er2p1, &b_L1_SingleIsoEG33er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG34", &L1_SingleIsoEG34, &b_L1_SingleIsoEG34);
   fChain->SetBranchAddress("L1_SingleIsoEG34er2p1", &L1_SingleIsoEG34er2p1, &b_L1_SingleIsoEG34er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG35", &L1_SingleIsoEG35, &b_L1_SingleIsoEG35);
   fChain->SetBranchAddress("L1_SingleIsoEG35er2p1", &L1_SingleIsoEG35er2p1, &b_L1_SingleIsoEG35er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG36", &L1_SingleIsoEG36, &b_L1_SingleIsoEG36);
   fChain->SetBranchAddress("L1_SingleIsoEG36er2p1", &L1_SingleIsoEG36er2p1, &b_L1_SingleIsoEG36er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG37", &L1_SingleIsoEG37, &b_L1_SingleIsoEG37);
   fChain->SetBranchAddress("L1_SingleIsoEG38", &L1_SingleIsoEG38, &b_L1_SingleIsoEG38);
   fChain->SetBranchAddress("L1_SingleIsoEG38er2p1", &L1_SingleIsoEG38er2p1, &b_L1_SingleIsoEG38er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG40", &L1_SingleIsoEG40, &b_L1_SingleIsoEG40);
   fChain->SetBranchAddress("L1_SingleIsoEG40er2p1", &L1_SingleIsoEG40er2p1, &b_L1_SingleIsoEG40er2p1);
   fChain->SetBranchAddress("L1_SingleJet120", &L1_SingleJet120, &b_L1_SingleJet120);
   fChain->SetBranchAddress("L1_SingleJet120_FWD", &L1_SingleJet120_FWD, &b_L1_SingleJet120_FWD);
   fChain->SetBranchAddress("L1_SingleJet12_BptxAND", &L1_SingleJet12_BptxAND, &b_L1_SingleJet12_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet140", &L1_SingleJet140, &b_L1_SingleJet140);
   fChain->SetBranchAddress("L1_SingleJet150", &L1_SingleJet150, &b_L1_SingleJet150);
   fChain->SetBranchAddress("L1_SingleJet16", &L1_SingleJet16, &b_L1_SingleJet16);
   fChain->SetBranchAddress("L1_SingleJet160", &L1_SingleJet160, &b_L1_SingleJet160);
   fChain->SetBranchAddress("L1_SingleJet170", &L1_SingleJet170, &b_L1_SingleJet170);
   fChain->SetBranchAddress("L1_SingleJet180", &L1_SingleJet180, &b_L1_SingleJet180);
   fChain->SetBranchAddress("L1_SingleJet20", &L1_SingleJet20, &b_L1_SingleJet20);
   fChain->SetBranchAddress("L1_SingleJet200", &L1_SingleJet200, &b_L1_SingleJet200);
   fChain->SetBranchAddress("L1_SingleJet20er2p7_NotBptxOR", &L1_SingleJet20er2p7_NotBptxOR, &b_L1_SingleJet20er2p7_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleJet20er2p7_NotBptxOR_3BX", &L1_SingleJet20er2p7_NotBptxOR_3BX, &b_L1_SingleJet20er2p7_NotBptxOR_3BX);
   fChain->SetBranchAddress("L1_SingleJet35", &L1_SingleJet35, &b_L1_SingleJet35);
   fChain->SetBranchAddress("L1_SingleJet35_FWD", &L1_SingleJet35_FWD, &b_L1_SingleJet35_FWD);
   fChain->SetBranchAddress("L1_SingleJet35_HFm", &L1_SingleJet35_HFm, &b_L1_SingleJet35_HFm);
   fChain->SetBranchAddress("L1_SingleJet35_HFp", &L1_SingleJet35_HFp, &b_L1_SingleJet35_HFp);
   fChain->SetBranchAddress("L1_SingleJet43er2p7_NotBptxOR_3BX", &L1_SingleJet43er2p7_NotBptxOR_3BX, &b_L1_SingleJet43er2p7_NotBptxOR_3BX);
   fChain->SetBranchAddress("L1_SingleJet46er2p7_NotBptxOR_3BX", &L1_SingleJet46er2p7_NotBptxOR_3BX, &b_L1_SingleJet46er2p7_NotBptxOR_3BX);
   fChain->SetBranchAddress("L1_SingleJet60", &L1_SingleJet60, &b_L1_SingleJet60);
   fChain->SetBranchAddress("L1_SingleJet60_FWD", &L1_SingleJet60_FWD, &b_L1_SingleJet60_FWD);
   fChain->SetBranchAddress("L1_SingleJet60_HFm", &L1_SingleJet60_HFm, &b_L1_SingleJet60_HFm);
   fChain->SetBranchAddress("L1_SingleJet60_HFp", &L1_SingleJet60_HFp, &b_L1_SingleJet60_HFp);
   fChain->SetBranchAddress("L1_SingleJet90", &L1_SingleJet90, &b_L1_SingleJet90);
   fChain->SetBranchAddress("L1_SingleJet90_FWD", &L1_SingleJet90_FWD, &b_L1_SingleJet90_FWD);
   fChain->SetBranchAddress("L1_SingleMu0_BMTF", &L1_SingleMu0_BMTF, &b_L1_SingleMu0_BMTF);
   fChain->SetBranchAddress("L1_SingleMu0_EMTF", &L1_SingleMu0_EMTF, &b_L1_SingleMu0_EMTF);
   fChain->SetBranchAddress("L1_SingleMu0_OMTF", &L1_SingleMu0_OMTF, &b_L1_SingleMu0_OMTF);
   fChain->SetBranchAddress("L1_SingleMu10_LowQ", &L1_SingleMu10_LowQ, &b_L1_SingleMu10_LowQ);
   fChain->SetBranchAddress("L1_SingleMu11_LowQ", &L1_SingleMu11_LowQ, &b_L1_SingleMu11_LowQ);
   fChain->SetBranchAddress("L1_SingleMu12_LowQ_BMTF", &L1_SingleMu12_LowQ_BMTF, &b_L1_SingleMu12_LowQ_BMTF);
   fChain->SetBranchAddress("L1_SingleMu12_LowQ_EMTF", &L1_SingleMu12_LowQ_EMTF, &b_L1_SingleMu12_LowQ_EMTF);
   fChain->SetBranchAddress("L1_SingleMu12_LowQ_OMTF", &L1_SingleMu12_LowQ_OMTF, &b_L1_SingleMu12_LowQ_OMTF);
   fChain->SetBranchAddress("L1_SingleMu14er2p1", &L1_SingleMu14er2p1, &b_L1_SingleMu14er2p1);
   fChain->SetBranchAddress("L1_SingleMu16", &L1_SingleMu16, &b_L1_SingleMu16);
   fChain->SetBranchAddress("L1_SingleMu16er2p1", &L1_SingleMu16er2p1, &b_L1_SingleMu16er2p1);
   fChain->SetBranchAddress("L1_SingleMu18", &L1_SingleMu18, &b_L1_SingleMu18);
   fChain->SetBranchAddress("L1_SingleMu18er2p1", &L1_SingleMu18er2p1, &b_L1_SingleMu18er2p1);
   fChain->SetBranchAddress("L1_SingleMu20", &L1_SingleMu20, &b_L1_SingleMu20);
   fChain->SetBranchAddress("L1_SingleMu20er2p1", &L1_SingleMu20er2p1, &b_L1_SingleMu20er2p1);
   fChain->SetBranchAddress("L1_SingleMu22", &L1_SingleMu22, &b_L1_SingleMu22);
   fChain->SetBranchAddress("L1_SingleMu22_BMTF", &L1_SingleMu22_BMTF, &b_L1_SingleMu22_BMTF);
   fChain->SetBranchAddress("L1_SingleMu22_EMTF", &L1_SingleMu22_EMTF, &b_L1_SingleMu22_EMTF);
   fChain->SetBranchAddress("L1_SingleMu22_OMTF", &L1_SingleMu22_OMTF, &b_L1_SingleMu22_OMTF);
   fChain->SetBranchAddress("L1_SingleMu22er2p1", &L1_SingleMu22er2p1, &b_L1_SingleMu22er2p1);
   fChain->SetBranchAddress("L1_SingleMu25", &L1_SingleMu25, &b_L1_SingleMu25);
   fChain->SetBranchAddress("L1_SingleMu3", &L1_SingleMu3, &b_L1_SingleMu3);
   fChain->SetBranchAddress("L1_SingleMu30", &L1_SingleMu30, &b_L1_SingleMu30);
   fChain->SetBranchAddress("L1_SingleMu5", &L1_SingleMu5, &b_L1_SingleMu5);
   fChain->SetBranchAddress("L1_SingleMu7", &L1_SingleMu7, &b_L1_SingleMu7);
   fChain->SetBranchAddress("L1_SingleMuCosmics", &L1_SingleMuCosmics, &b_L1_SingleMuCosmics);
   fChain->SetBranchAddress("L1_SingleMuCosmics_BMTF", &L1_SingleMuCosmics_BMTF, &b_L1_SingleMuCosmics_BMTF);
   fChain->SetBranchAddress("L1_SingleMuCosmics_EMTF", &L1_SingleMuCosmics_EMTF, &b_L1_SingleMuCosmics_EMTF);
   fChain->SetBranchAddress("L1_SingleMuCosmics_OMTF", &L1_SingleMuCosmics_OMTF, &b_L1_SingleMuCosmics_OMTF);
   fChain->SetBranchAddress("L1_SingleMuOpen", &L1_SingleMuOpen, &b_L1_SingleMuOpen);
   fChain->SetBranchAddress("L1_SingleMuOpen_NotBptxOR", &L1_SingleMuOpen_NotBptxOR, &b_L1_SingleMuOpen_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleMuOpen_NotBptxOR_3BX", &L1_SingleMuOpen_NotBptxOR_3BX, &b_L1_SingleMuOpen_NotBptxOR_3BX);
   fChain->SetBranchAddress("L1_SingleTau100er2p1", &L1_SingleTau100er2p1, &b_L1_SingleTau100er2p1);
   fChain->SetBranchAddress("L1_SingleTau120er2p1", &L1_SingleTau120er2p1, &b_L1_SingleTau120er2p1);
   fChain->SetBranchAddress("L1_SingleTau130er2p1", &L1_SingleTau130er2p1, &b_L1_SingleTau130er2p1);
   fChain->SetBranchAddress("L1_SingleTau140er2p1", &L1_SingleTau140er2p1, &b_L1_SingleTau140er2p1);
   fChain->SetBranchAddress("L1_SingleTau20", &L1_SingleTau20, &b_L1_SingleTau20);
   fChain->SetBranchAddress("L1_SingleTau80er2p1", &L1_SingleTau80er2p1, &b_L1_SingleTau80er2p1);
   fChain->SetBranchAddress("L1_TripleEG_14_10_8", &L1_TripleEG_14_10_8, &b_L1_TripleEG_14_10_8);
   fChain->SetBranchAddress("L1_TripleEG_18_17_8", &L1_TripleEG_18_17_8, &b_L1_TripleEG_18_17_8);
   fChain->SetBranchAddress("L1_TripleEG_LooseIso20_10_5", &L1_TripleEG_LooseIso20_10_5, &b_L1_TripleEG_LooseIso20_10_5);
   fChain->SetBranchAddress("L1_TripleJet_100_85_72_VBF", &L1_TripleJet_100_85_72_VBF, &b_L1_TripleJet_100_85_72_VBF);
   fChain->SetBranchAddress("L1_TripleJet_105_85_76_VBF", &L1_TripleJet_105_85_76_VBF, &b_L1_TripleJet_105_85_76_VBF);
   fChain->SetBranchAddress("L1_TripleJet_84_68_48_VBF", &L1_TripleJet_84_68_48_VBF, &b_L1_TripleJet_84_68_48_VBF);
   fChain->SetBranchAddress("L1_TripleJet_88_72_56_VBF", &L1_TripleJet_88_72_56_VBF, &b_L1_TripleJet_88_72_56_VBF);
   fChain->SetBranchAddress("L1_TripleJet_92_76_64_VBF", &L1_TripleJet_92_76_64_VBF, &b_L1_TripleJet_92_76_64_VBF);
   fChain->SetBranchAddress("L1_TripleJet_98_83_71_VBF", &L1_TripleJet_98_83_71_VBF, &b_L1_TripleJet_98_83_71_VBF);
   fChain->SetBranchAddress("L1_TripleMu0", &L1_TripleMu0, &b_L1_TripleMu0);
   fChain->SetBranchAddress("L1_TripleMu0_OQ", &L1_TripleMu0_OQ, &b_L1_TripleMu0_OQ);
   fChain->SetBranchAddress("L1_TripleMu3", &L1_TripleMu3, &b_L1_TripleMu3);
   fChain->SetBranchAddress("L1_TripleMu3_SQ", &L1_TripleMu3_SQ, &b_L1_TripleMu3_SQ);
   fChain->SetBranchAddress("L1_TripleMu_4_4_4", &L1_TripleMu_4_4_4, &b_L1_TripleMu_4_4_4);
   fChain->SetBranchAddress("L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", &L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17, &b_L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17);
   fChain->SetBranchAddress("L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14", &L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14, &b_L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14);
   fChain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0OQ", &L1_TripleMu_5SQ_3SQ_0OQ, &b_L1_TripleMu_5SQ_3SQ_0OQ);
   fChain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9, &b_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9);
   fChain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9, &b_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9);
   fChain->SetBranchAddress("L1_TripleMu_5_0_0", &L1_TripleMu_5_0_0, &b_L1_TripleMu_5_0_0);
   fChain->SetBranchAddress("L1_TripleMu_5_3_3", &L1_TripleMu_5_3_3, &b_L1_TripleMu_5_3_3);
   fChain->SetBranchAddress("L1_TripleMu_5_3p5_2p5", &L1_TripleMu_5_3p5_2p5, &b_L1_TripleMu_5_3p5_2p5);
   fChain->SetBranchAddress("L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17, &b_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17);
   fChain->SetBranchAddress("L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17, &b_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17);
   fChain->SetBranchAddress("L1_TripleMu_5_5_3", &L1_TripleMu_5_5_3, &b_L1_TripleMu_5_5_3);
   fChain->SetBranchAddress("L1_UnpairedBunchBptxMinus", &L1_UnpairedBunchBptxMinus, &b_L1_UnpairedBunchBptxMinus);
   fChain->SetBranchAddress("L1_UnpairedBunchBptxPlus", &L1_UnpairedBunchBptxPlus, &b_L1_UnpairedBunchBptxPlus);
   fChain->SetBranchAddress("L1_ZeroBias", &L1_ZeroBias, &b_L1_ZeroBias);
   fChain->SetBranchAddress("L1_ZeroBias_copy", &L1_ZeroBias_copy, &b_L1_ZeroBias_copy);
   fChain->SetBranchAddress("L1_UnprefireableEvent", &L1_UnprefireableEvent, &b_L1_UnprefireableEvent);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter, &b_Flag_BadPFMuonDzFilter);
   fChain->SetBranchAddress("Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter, &b_Flag_hfNoisyHitsFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter, &b_Flag_BadChargedCandidateSummer16Filter);
   fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter, &b_Flag_BadPFMuonSummer16Filter);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
   fChain->SetBranchAddress("L1Reco_step", &L1Reco_step, &b_L1Reco_step);
   fChain->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath, &b_HLTriggerFirstPath);
   fChain->SetBranchAddress("HLT_AK8PFJet360_TrimMass30", &HLT_AK8PFJet360_TrimMass30, &b_HLT_AK8PFJet360_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet380_TrimMass30", &HLT_AK8PFJet380_TrimMass30, &b_HLT_AK8PFJet380_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet400_TrimMass30", &HLT_AK8PFJet400_TrimMass30, &b_HLT_AK8PFJet400_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet420_TrimMass30", &HLT_AK8PFJet420_TrimMass30, &b_HLT_AK8PFJet420_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFHT750_TrimMass50", &HLT_AK8PFHT750_TrimMass50, &b_HLT_AK8PFHT750_TrimMass50);
   fChain->SetBranchAddress("HLT_AK8PFHT800_TrimMass50", &HLT_AK8PFHT800_TrimMass50, &b_HLT_AK8PFHT800_TrimMass50);
   fChain->SetBranchAddress("HLT_AK8PFHT850_TrimMass50", &HLT_AK8PFHT850_TrimMass50, &b_HLT_AK8PFHT850_TrimMass50);
   fChain->SetBranchAddress("HLT_AK8PFHT900_TrimMass50", &HLT_AK8PFHT900_TrimMass50, &b_HLT_AK8PFHT900_TrimMass50);
   fChain->SetBranchAddress("HLT_CaloJet500_NoJetID", &HLT_CaloJet500_NoJetID, &b_HLT_CaloJet500_NoJetID);
   fChain->SetBranchAddress("HLT_CaloJet550_NoJetID", &HLT_CaloJet550_NoJetID, &b_HLT_CaloJet550_NoJetID);
   fChain->SetBranchAddress("HLT_Trimuon5_3p5_2_Upsilon_Muon", &HLT_Trimuon5_3p5_2_Upsilon_Muon, &b_HLT_Trimuon5_3p5_2_Upsilon_Muon);
   fChain->SetBranchAddress("HLT_DoubleEle25_CaloIdL_MW", &HLT_DoubleEle25_CaloIdL_MW, &b_HLT_DoubleEle25_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle27_CaloIdL_MW", &HLT_DoubleEle27_CaloIdL_MW, &b_HLT_DoubleEle27_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL_MW", &HLT_DoubleEle33_CaloIdL_MW, &b_HLT_DoubleEle33_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle24_eta2p1_WPTight_Gsf", &HLT_DoubleEle24_eta2p1_WPTight_Gsf, &b_HLT_DoubleEle24_eta2p1_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350, &b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350);
   fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350, &b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350);
   fChain->SetBranchAddress("HLT_Ele27_Ele37_CaloIdL_MW", &HLT_Ele27_Ele37_CaloIdL_MW, &b_HLT_Ele27_Ele37_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_Mu27_Ele37_CaloIdL_MW", &HLT_Mu27_Ele37_CaloIdL_MW, &b_HLT_Mu27_Ele37_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_Mu37_Ele27_CaloIdL_MW", &HLT_Mu37_Ele27_CaloIdL_MW, &b_HLT_Mu37_Ele27_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_Mu37_TkMu27", &HLT_Mu37_TkMu27, &b_HLT_Mu37_TkMu27);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Bs", &HLT_DoubleMu4_3_Bs, &b_HLT_DoubleMu4_3_Bs);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Jpsi_Displaced", &HLT_DoubleMu4_3_Jpsi_Displaced, &b_HLT_DoubleMu4_3_Jpsi_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Displaced", &HLT_DoubleMu4_JpsiTrk_Displaced, &b_HLT_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced", &HLT_DoubleMu4_LowMassNonResonantTrk_Displaced, &b_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu", &HLT_DoubleMu3_Trk_Tau3mu, &b_HLT_DoubleMu3_Trk_Tau3mu);
   fChain->SetBranchAddress("HLT_DoubleMu4_PsiPrimeTrk_Displaced", &HLT_DoubleMu4_PsiPrimeTrk_Displaced, &b_HLT_DoubleMu4_PsiPrimeTrk_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_Mass8_DZ_PFHT350", &HLT_DoubleMu4_Mass8_DZ_PFHT350, &b_HLT_DoubleMu4_Mass8_DZ_PFHT350);
   fChain->SetBranchAddress("HLT_DoubleMu8_Mass8_PFHT350", &HLT_DoubleMu8_Mass8_PFHT350, &b_HLT_DoubleMu8_Mass8_PFHT350);
   fChain->SetBranchAddress("HLT_Mu3_PFJet40", &HLT_Mu3_PFJet40, &b_HLT_Mu3_PFJet40);
   fChain->SetBranchAddress("HLT_Mu7p5_L2Mu2_Jpsi", &HLT_Mu7p5_L2Mu2_Jpsi, &b_HLT_Mu7p5_L2Mu2_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_L2Mu2_Upsilon", &HLT_Mu7p5_L2Mu2_Upsilon, &b_HLT_Mu7p5_L2Mu2_Upsilon);
   fChain->SetBranchAddress("HLT_Mu7p5_Track2_Jpsi", &HLT_Mu7p5_Track2_Jpsi, &b_HLT_Mu7p5_Track2_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_Track3p5_Jpsi", &HLT_Mu7p5_Track3p5_Jpsi, &b_HLT_Mu7p5_Track3p5_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_Track7_Jpsi", &HLT_Mu7p5_Track7_Jpsi, &b_HLT_Mu7p5_Track7_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_Track2_Upsilon", &HLT_Mu7p5_Track2_Upsilon, &b_HLT_Mu7p5_Track2_Upsilon);
   fChain->SetBranchAddress("HLT_Mu7p5_Track3p5_Upsilon", &HLT_Mu7p5_Track3p5_Upsilon, &b_HLT_Mu7p5_Track3p5_Upsilon);
   fChain->SetBranchAddress("HLT_Mu7p5_Track7_Upsilon", &HLT_Mu7p5_Track7_Upsilon, &b_HLT_Mu7p5_Track7_Upsilon);
   fChain->SetBranchAddress("HLT_DoublePhoton33_CaloIdL", &HLT_DoublePhoton33_CaloIdL, &b_HLT_DoublePhoton33_CaloIdL);
   fChain->SetBranchAddress("HLT_DoublePhoton70", &HLT_DoublePhoton70, &b_HLT_DoublePhoton70);
   fChain->SetBranchAddress("HLT_DoublePhoton85", &HLT_DoublePhoton85, &b_HLT_DoublePhoton85);
   fChain->SetBranchAddress("HLT_Ele20_WPTight_Gsf", &HLT_Ele20_WPTight_Gsf, &b_HLT_Ele20_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele20_WPLoose_Gsf", &HLT_Ele20_WPLoose_Gsf, &b_HLT_Ele20_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele20_eta2p1_WPLoose_Gsf", &HLT_Ele20_eta2p1_WPLoose_Gsf, &b_HLT_Ele20_eta2p1_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG", &HLT_DiEle27_WPTightCaloOnly_L1DoubleEG, &b_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG);
   fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf, &b_HLT_Ele27_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf, &b_HLT_Ele32_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf, &b_HLT_Ele35_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf_L1EGMT", &HLT_Ele35_WPTight_Gsf_L1EGMT, &b_HLT_Ele35_WPTight_Gsf_L1EGMT);
   fChain->SetBranchAddress("HLT_Ele38_WPTight_Gsf", &HLT_Ele38_WPTight_Gsf, &b_HLT_Ele38_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele40_WPTight_Gsf", &HLT_Ele40_WPTight_Gsf, &b_HLT_Ele40_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG", &HLT_Ele32_WPTight_Gsf_L1DoubleEG, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG);
   fChain->SetBranchAddress("HLT_HT450_Beamspot", &HLT_HT450_Beamspot, &b_HLT_HT450_Beamspot);
   fChain->SetBranchAddress("HLT_HT300_Beamspot", &HLT_HT300_Beamspot, &b_HLT_HT300_Beamspot);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1, &b_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1, &b_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1, &b_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1, &b_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1, &b_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1", &HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1, &b_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1", &HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1", &HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1, &b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20, &b_HLT_IsoMu20);
   fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
   fChain->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu24, &b_HLT_IsoTkMu24);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1", &HLT_IsoMu24_eta2p1, &b_HLT_IsoMu24_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
   fChain->SetBranchAddress("HLT_IsoMu30", &HLT_IsoMu30, &b_HLT_IsoMu30);
   fChain->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX", &HLT_UncorrectedJetE30_NoBPTX, &b_HLT_UncorrectedJetE30_NoBPTX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX3BX", &HLT_UncorrectedJetE30_NoBPTX3BX, &b_HLT_UncorrectedJetE30_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE60_NoBPTX3BX", &HLT_UncorrectedJetE60_NoBPTX3BX, &b_HLT_UncorrectedJetE60_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE70_NoBPTX3BX", &HLT_UncorrectedJetE70_NoBPTX3BX, &b_HLT_UncorrectedJetE70_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L1SingleMu18", &HLT_L1SingleMu18, &b_HLT_L1SingleMu18);
   fChain->SetBranchAddress("HLT_L1SingleMu25", &HLT_L1SingleMu25, &b_HLT_L1SingleMu25);
   fChain->SetBranchAddress("HLT_L2Mu10", &HLT_L2Mu10, &b_HLT_L2Mu10);
   fChain->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX3BX", &HLT_L2Mu10_NoVertex_NoBPTX3BX, &b_HLT_L2Mu10_NoVertex_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX", &HLT_L2Mu10_NoVertex_NoBPTX, &b_HLT_L2Mu10_NoVertex_NoBPTX);
   fChain->SetBranchAddress("HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX, &b_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX, &b_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L2Mu50", &HLT_L2Mu50, &b_HLT_L2Mu50);
   fChain->SetBranchAddress("HLT_DoubleL2Mu50", &HLT_DoubleL2Mu50, &b_HLT_DoubleL2Mu50);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Mu25_TkMu0_Onia", &HLT_Mu25_TkMu0_Onia, &b_HLT_Mu25_TkMu0_Onia);
   fChain->SetBranchAddress("HLT_Mu30_TkMu0_Onia", &HLT_Mu30_TkMu0_Onia, &b_HLT_Mu30_TkMu0_Onia);
   fChain->SetBranchAddress("HLT_Mu20_TkMu0_Phi", &HLT_Mu20_TkMu0_Phi, &b_HLT_Mu20_TkMu0_Phi);
   fChain->SetBranchAddress("HLT_Mu25_TkMu0_Phi", &HLT_Mu25_TkMu0_Phi, &b_HLT_Mu25_TkMu0_Phi);
   fChain->SetBranchAddress("HLT_Mu20", &HLT_Mu20, &b_HLT_Mu20);
   fChain->SetBranchAddress("HLT_Mu27", &HLT_Mu27, &b_HLT_Mu27);
   fChain->SetBranchAddress("HLT_Mu50", &HLT_Mu50, &b_HLT_Mu50);
   fChain->SetBranchAddress("HLT_Mu55", &HLT_Mu55, &b_HLT_Mu55);
   fChain->SetBranchAddress("HLT_OldMu100", &HLT_OldMu100, &b_HLT_OldMu100);
   fChain->SetBranchAddress("HLT_TkMu100", &HLT_TkMu100, &b_HLT_TkMu100);
   fChain->SetBranchAddress("HLT_DiPFJet15_NoCaloMatched", &HLT_DiPFJet15_NoCaloMatched, &b_HLT_DiPFJet15_NoCaloMatched);
   fChain->SetBranchAddress("HLT_DiPFJet25_NoCaloMatched", &HLT_DiPFJet25_NoCaloMatched, &b_HLT_DiPFJet25_NoCaloMatched);
   fChain->SetBranchAddress("HLT_DiPFJet15_FBEta3_NoCaloMatched", &HLT_DiPFJet15_FBEta3_NoCaloMatched, &b_HLT_DiPFJet15_FBEta3_NoCaloMatched);
   fChain->SetBranchAddress("HLT_DiPFJet25_FBEta3_NoCaloMatched", &HLT_DiPFJet25_FBEta3_NoCaloMatched, &b_HLT_DiPFJet25_FBEta3_NoCaloMatched);
   fChain->SetBranchAddress("HLT_DiPFJetAve40", &HLT_DiPFJetAve40, &b_HLT_DiPFJetAve40);
   fChain->SetBranchAddress("HLT_DiPFJetAve60", &HLT_DiPFJetAve60, &b_HLT_DiPFJetAve60);
   fChain->SetBranchAddress("HLT_DiPFJetAve80", &HLT_DiPFJetAve80, &b_HLT_DiPFJetAve80);
   fChain->SetBranchAddress("HLT_DiPFJetAve140", &HLT_DiPFJetAve140, &b_HLT_DiPFJetAve140);
   fChain->SetBranchAddress("HLT_DiPFJetAve200", &HLT_DiPFJetAve200, &b_HLT_DiPFJetAve200);
   fChain->SetBranchAddress("HLT_DiPFJetAve260", &HLT_DiPFJetAve260, &b_HLT_DiPFJetAve260);
   fChain->SetBranchAddress("HLT_DiPFJetAve320", &HLT_DiPFJetAve320, &b_HLT_DiPFJetAve320);
   fChain->SetBranchAddress("HLT_DiPFJetAve400", &HLT_DiPFJetAve400, &b_HLT_DiPFJetAve400);
   fChain->SetBranchAddress("HLT_DiPFJetAve500", &HLT_DiPFJetAve500, &b_HLT_DiPFJetAve500);
   fChain->SetBranchAddress("HLT_DiPFJetAve15_HFJEC", &HLT_DiPFJetAve15_HFJEC, &b_HLT_DiPFJetAve15_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve25_HFJEC", &HLT_DiPFJetAve25_HFJEC, &b_HLT_DiPFJetAve25_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve35_HFJEC", &HLT_DiPFJetAve35_HFJEC, &b_HLT_DiPFJetAve35_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve60_HFJEC", &HLT_DiPFJetAve60_HFJEC, &b_HLT_DiPFJetAve60_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve80_HFJEC", &HLT_DiPFJetAve80_HFJEC, &b_HLT_DiPFJetAve80_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve100_HFJEC", &HLT_DiPFJetAve100_HFJEC, &b_HLT_DiPFJetAve100_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve160_HFJEC", &HLT_DiPFJetAve160_HFJEC, &b_HLT_DiPFJetAve160_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve220_HFJEC", &HLT_DiPFJetAve220_HFJEC, &b_HLT_DiPFJetAve220_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve300_HFJEC", &HLT_DiPFJetAve300_HFJEC, &b_HLT_DiPFJetAve300_HFJEC);
   fChain->SetBranchAddress("HLT_AK8PFJet40", &HLT_AK8PFJet40, &b_HLT_AK8PFJet40);
   fChain->SetBranchAddress("HLT_AK8PFJet60", &HLT_AK8PFJet60, &b_HLT_AK8PFJet60);
   fChain->SetBranchAddress("HLT_AK8PFJet80", &HLT_AK8PFJet80, &b_HLT_AK8PFJet80);
   fChain->SetBranchAddress("HLT_AK8PFJet140", &HLT_AK8PFJet140, &b_HLT_AK8PFJet140);
   fChain->SetBranchAddress("HLT_AK8PFJet200", &HLT_AK8PFJet200, &b_HLT_AK8PFJet200);
   fChain->SetBranchAddress("HLT_AK8PFJet260", &HLT_AK8PFJet260, &b_HLT_AK8PFJet260);
   fChain->SetBranchAddress("HLT_AK8PFJet320", &HLT_AK8PFJet320, &b_HLT_AK8PFJet320);
   fChain->SetBranchAddress("HLT_AK8PFJet400", &HLT_AK8PFJet400, &b_HLT_AK8PFJet400);
   fChain->SetBranchAddress("HLT_AK8PFJet450", &HLT_AK8PFJet450, &b_HLT_AK8PFJet450);
   fChain->SetBranchAddress("HLT_AK8PFJet500", &HLT_AK8PFJet500, &b_HLT_AK8PFJet500);
   fChain->SetBranchAddress("HLT_AK8PFJet550", &HLT_AK8PFJet550, &b_HLT_AK8PFJet550);
   fChain->SetBranchAddress("HLT_PFJet40", &HLT_PFJet40, &b_HLT_PFJet40);
   fChain->SetBranchAddress("HLT_PFJet60", &HLT_PFJet60, &b_HLT_PFJet60);
   fChain->SetBranchAddress("HLT_PFJet80", &HLT_PFJet80, &b_HLT_PFJet80);
   fChain->SetBranchAddress("HLT_PFJet140", &HLT_PFJet140, &b_HLT_PFJet140);
   fChain->SetBranchAddress("HLT_PFJet200", &HLT_PFJet200, &b_HLT_PFJet200);
   fChain->SetBranchAddress("HLT_PFJet260", &HLT_PFJet260, &b_HLT_PFJet260);
   fChain->SetBranchAddress("HLT_PFJet320", &HLT_PFJet320, &b_HLT_PFJet320);
   fChain->SetBranchAddress("HLT_PFJet400", &HLT_PFJet400, &b_HLT_PFJet400);
   fChain->SetBranchAddress("HLT_PFJet450", &HLT_PFJet450, &b_HLT_PFJet450);
   fChain->SetBranchAddress("HLT_PFJet500", &HLT_PFJet500, &b_HLT_PFJet500);
   fChain->SetBranchAddress("HLT_PFJet550", &HLT_PFJet550, &b_HLT_PFJet550);
   fChain->SetBranchAddress("HLT_PFJetFwd40", &HLT_PFJetFwd40, &b_HLT_PFJetFwd40);
   fChain->SetBranchAddress("HLT_PFJetFwd60", &HLT_PFJetFwd60, &b_HLT_PFJetFwd60);
   fChain->SetBranchAddress("HLT_PFJetFwd80", &HLT_PFJetFwd80, &b_HLT_PFJetFwd80);
   fChain->SetBranchAddress("HLT_PFJetFwd140", &HLT_PFJetFwd140, &b_HLT_PFJetFwd140);
   fChain->SetBranchAddress("HLT_PFJetFwd200", &HLT_PFJetFwd200, &b_HLT_PFJetFwd200);
   fChain->SetBranchAddress("HLT_PFJetFwd260", &HLT_PFJetFwd260, &b_HLT_PFJetFwd260);
   fChain->SetBranchAddress("HLT_PFJetFwd320", &HLT_PFJetFwd320, &b_HLT_PFJetFwd320);
   fChain->SetBranchAddress("HLT_PFJetFwd400", &HLT_PFJetFwd400, &b_HLT_PFJetFwd400);
   fChain->SetBranchAddress("HLT_PFJetFwd450", &HLT_PFJetFwd450, &b_HLT_PFJetFwd450);
   fChain->SetBranchAddress("HLT_PFJetFwd500", &HLT_PFJetFwd500, &b_HLT_PFJetFwd500);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd40", &HLT_AK8PFJetFwd40, &b_HLT_AK8PFJetFwd40);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd60", &HLT_AK8PFJetFwd60, &b_HLT_AK8PFJetFwd60);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd80", &HLT_AK8PFJetFwd80, &b_HLT_AK8PFJetFwd80);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd140", &HLT_AK8PFJetFwd140, &b_HLT_AK8PFJetFwd140);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd200", &HLT_AK8PFJetFwd200, &b_HLT_AK8PFJetFwd200);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd260", &HLT_AK8PFJetFwd260, &b_HLT_AK8PFJetFwd260);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd320", &HLT_AK8PFJetFwd320, &b_HLT_AK8PFJetFwd320);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd400", &HLT_AK8PFJetFwd400, &b_HLT_AK8PFJetFwd400);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd450", &HLT_AK8PFJetFwd450, &b_HLT_AK8PFJetFwd450);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd500", &HLT_AK8PFJetFwd500, &b_HLT_AK8PFJetFwd500);
   fChain->SetBranchAddress("HLT_PFHT180", &HLT_PFHT180, &b_HLT_PFHT180);
   fChain->SetBranchAddress("HLT_PFHT250", &HLT_PFHT250, &b_HLT_PFHT250);
   fChain->SetBranchAddress("HLT_PFHT370", &HLT_PFHT370, &b_HLT_PFHT370);
   fChain->SetBranchAddress("HLT_PFHT430", &HLT_PFHT430, &b_HLT_PFHT430);
   fChain->SetBranchAddress("HLT_PFHT510", &HLT_PFHT510, &b_HLT_PFHT510);
   fChain->SetBranchAddress("HLT_PFHT590", &HLT_PFHT590, &b_HLT_PFHT590);
   fChain->SetBranchAddress("HLT_PFHT680", &HLT_PFHT680, &b_HLT_PFHT680);
   fChain->SetBranchAddress("HLT_PFHT780", &HLT_PFHT780, &b_HLT_PFHT780);
   fChain->SetBranchAddress("HLT_PFHT890", &HLT_PFHT890, &b_HLT_PFHT890);
   fChain->SetBranchAddress("HLT_PFHT1050", &HLT_PFHT1050, &b_HLT_PFHT1050);
   fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight", &HLT_PFHT500_PFMET100_PFMHT100_IDTight, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight);
   fChain->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight", &HLT_PFHT500_PFMET110_PFMHT110_IDTight, &b_HLT_PFHT500_PFMET110_PFMHT110_IDTight);
   fChain->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight", &HLT_PFHT700_PFMET85_PFMHT85_IDTight, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight);
   fChain->SetBranchAddress("HLT_PFHT700_PFMET95_PFMHT95_IDTight", &HLT_PFHT700_PFMET95_PFMHT95_IDTight, &b_HLT_PFHT700_PFMET95_PFMHT95_IDTight);
   fChain->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight", &HLT_PFHT800_PFMET75_PFMHT75_IDTight, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight);
   fChain->SetBranchAddress("HLT_PFHT800_PFMET85_PFMHT85_IDTight", &HLT_PFHT800_PFMET85_PFMHT85_IDTight, &b_HLT_PFHT800_PFMET85_PFMHT85_IDTight);
   fChain->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight", &HLT_PFMET110_PFMHT110_IDTight, &b_HLT_PFMET110_PFMHT110_IDTight);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight, &b_HLT_PFMET120_PFMHT120_IDTight);
   fChain->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight", &HLT_PFMET130_PFMHT130_IDTight, &b_HLT_PFMET130_PFMHT130_IDTight);
   fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight", &HLT_PFMET140_PFMHT140_IDTight, &b_HLT_PFMET140_PFMHT140_IDTight);
   fChain->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1", &HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1);
   fChain->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1", &HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1", &HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1);
   fChain->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1", &HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1);
   fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1", &HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1, &b_HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60", &HLT_PFMET120_PFMHT120_IDTight_PFHT60, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60", &HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60, &b_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETTypeOne110_PFMHT110_IDTight", &HLT_PFMETTypeOne110_PFMHT110_IDTight, &b_HLT_PFMETTypeOne110_PFMHT110_IDTight);
   fChain->SetBranchAddress("HLT_PFMETTypeOne120_PFMHT120_IDTight", &HLT_PFMETTypeOne120_PFMHT120_IDTight, &b_HLT_PFMETTypeOne120_PFMHT120_IDTight);
   fChain->SetBranchAddress("HLT_PFMETTypeOne130_PFMHT130_IDTight", &HLT_PFMETTypeOne130_PFMHT130_IDTight, &b_HLT_PFMETTypeOne130_PFMHT130_IDTight);
   fChain->SetBranchAddress("HLT_PFMETTypeOne140_PFMHT140_IDTight", &HLT_PFMETTypeOne140_PFMHT140_IDTight, &b_HLT_PFMETTypeOne140_PFMHT140_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight, &b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight);
   fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight);
   fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight);
   fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight);
   fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight);
   fChain->SetBranchAddress("HLT_L1ETMHadSeeds", &HLT_L1ETMHadSeeds, &b_HLT_L1ETMHadSeeds);
   fChain->SetBranchAddress("HLT_CaloMHT90", &HLT_CaloMHT90, &b_HLT_CaloMHT90);
   fChain->SetBranchAddress("HLT_CaloMET80_NotCleaned", &HLT_CaloMET80_NotCleaned, &b_HLT_CaloMET80_NotCleaned);
   fChain->SetBranchAddress("HLT_CaloMET90_NotCleaned", &HLT_CaloMET90_NotCleaned, &b_HLT_CaloMET90_NotCleaned);
   fChain->SetBranchAddress("HLT_CaloMET100_NotCleaned", &HLT_CaloMET100_NotCleaned, &b_HLT_CaloMET100_NotCleaned);
   fChain->SetBranchAddress("HLT_CaloMET110_NotCleaned", &HLT_CaloMET110_NotCleaned, &b_HLT_CaloMET110_NotCleaned);
   fChain->SetBranchAddress("HLT_CaloMET250_NotCleaned", &HLT_CaloMET250_NotCleaned, &b_HLT_CaloMET250_NotCleaned);
   fChain->SetBranchAddress("HLT_CaloMET70_HBHECleaned", &HLT_CaloMET70_HBHECleaned, &b_HLT_CaloMET70_HBHECleaned);
   fChain->SetBranchAddress("HLT_CaloMET80_HBHECleaned", &HLT_CaloMET80_HBHECleaned, &b_HLT_CaloMET80_HBHECleaned);
   fChain->SetBranchAddress("HLT_CaloMET90_HBHECleaned", &HLT_CaloMET90_HBHECleaned, &b_HLT_CaloMET90_HBHECleaned);
   fChain->SetBranchAddress("HLT_CaloMET100_HBHECleaned", &HLT_CaloMET100_HBHECleaned, &b_HLT_CaloMET100_HBHECleaned);
   fChain->SetBranchAddress("HLT_CaloMET250_HBHECleaned", &HLT_CaloMET250_HBHECleaned, &b_HLT_CaloMET250_HBHECleaned);
   fChain->SetBranchAddress("HLT_CaloMET300_HBHECleaned", &HLT_CaloMET300_HBHECleaned, &b_HLT_CaloMET300_HBHECleaned);
   fChain->SetBranchAddress("HLT_CaloMET350_HBHECleaned", &HLT_CaloMET350_HBHECleaned, &b_HLT_CaloMET350_HBHECleaned);
   fChain->SetBranchAddress("HLT_PFMET200_NotCleaned", &HLT_PFMET200_NotCleaned, &b_HLT_PFMET200_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET200_HBHECleaned", &HLT_PFMET200_HBHECleaned, &b_HLT_PFMET200_HBHECleaned);
   fChain->SetBranchAddress("HLT_PFMET250_HBHECleaned", &HLT_PFMET250_HBHECleaned, &b_HLT_PFMET250_HBHECleaned);
   fChain->SetBranchAddress("HLT_PFMET300_HBHECleaned", &HLT_PFMET300_HBHECleaned, &b_HLT_PFMET300_HBHECleaned);
   fChain->SetBranchAddress("HLT_PFMET200_HBHE_BeamHaloCleaned", &HLT_PFMET200_HBHE_BeamHaloCleaned, &b_HLT_PFMET200_HBHE_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned", &HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned, &b_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_MET105_IsoTrk50", &HLT_MET105_IsoTrk50, &b_HLT_MET105_IsoTrk50);
   fChain->SetBranchAddress("HLT_MET120_IsoTrk50", &HLT_MET120_IsoTrk50, &b_HLT_MET120_IsoTrk50);
   fChain->SetBranchAddress("HLT_SingleJet30_Mu12_SinglePFJet40", &HLT_SingleJet30_Mu12_SinglePFJet40, &b_HLT_SingleJet30_Mu12_SinglePFJet40);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33", &HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33", &HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33", &HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33", &HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_DoublePFJets40_CaloBTagCSV_p33", &HLT_DoublePFJets40_CaloBTagCSV_p33, &b_HLT_DoublePFJets40_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_DoublePFJets100_CaloBTagCSV_p33", &HLT_DoublePFJets100_CaloBTagCSV_p33, &b_HLT_DoublePFJets100_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_DoublePFJets200_CaloBTagCSV_p33", &HLT_DoublePFJets200_CaloBTagCSV_p33, &b_HLT_DoublePFJets200_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_DoublePFJets350_CaloBTagCSV_p33", &HLT_DoublePFJets350_CaloBTagCSV_p33, &b_HLT_DoublePFJets350_CaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33", &HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33, &b_HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33);
   fChain->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE, &b_HLT_Photon300_NoHE);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL", &HLT_Mu8_TrkIsoVVL, &b_HLT_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ, &b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ);
   fChain->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL, &b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ);
   fChain->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet20_Mu5", &HLT_BTagMu_AK4DiJet20_Mu5, &b_HLT_BTagMu_AK4DiJet20_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet40_Mu5", &HLT_BTagMu_AK4DiJet40_Mu5, &b_HLT_BTagMu_AK4DiJet40_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet70_Mu5", &HLT_BTagMu_AK4DiJet70_Mu5, &b_HLT_BTagMu_AK4DiJet70_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet110_Mu5", &HLT_BTagMu_AK4DiJet110_Mu5, &b_HLT_BTagMu_AK4DiJet110_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet170_Mu5", &HLT_BTagMu_AK4DiJet170_Mu5, &b_HLT_BTagMu_AK4DiJet170_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4Jet300_Mu5", &HLT_BTagMu_AK4Jet300_Mu5, &b_HLT_BTagMu_AK4Jet300_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK8DiJet170_Mu5", &HLT_BTagMu_AK8DiJet170_Mu5, &b_HLT_BTagMu_AK8DiJet170_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK8Jet300_Mu5", &HLT_BTagMu_AK8Jet300_Mu5, &b_HLT_BTagMu_AK8Jet300_Mu5);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu12_DoublePhoton20", &HLT_Mu12_DoublePhoton20, &b_HLT_Mu12_DoublePhoton20);
   fChain->SetBranchAddress("HLT_TriplePhoton_20_20_20_CaloIdLV2", &HLT_TriplePhoton_20_20_20_CaloIdLV2, &b_HLT_TriplePhoton_20_20_20_CaloIdLV2);
   fChain->SetBranchAddress("HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL, &b_HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL);
   fChain->SetBranchAddress("HLT_TriplePhoton_30_30_10_CaloIdLV2", &HLT_TriplePhoton_30_30_10_CaloIdLV2, &b_HLT_TriplePhoton_30_30_10_CaloIdLV2);
   fChain->SetBranchAddress("HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL, &b_HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL);
   fChain->SetBranchAddress("HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL", &HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL, &b_HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL);
   fChain->SetBranchAddress("HLT_Photon25", &HLT_Photon25, &b_HLT_Photon25);
   fChain->SetBranchAddress("HLT_Photon33", &HLT_Photon33, &b_HLT_Photon33);
   fChain->SetBranchAddress("HLT_Photon50", &HLT_Photon50, &b_HLT_Photon50);
   fChain->SetBranchAddress("HLT_Photon75", &HLT_Photon75, &b_HLT_Photon75);
   fChain->SetBranchAddress("HLT_Photon90", &HLT_Photon90, &b_HLT_Photon90);
   fChain->SetBranchAddress("HLT_Photon120", &HLT_Photon120, &b_HLT_Photon120);
   fChain->SetBranchAddress("HLT_Photon150", &HLT_Photon150, &b_HLT_Photon150);
   fChain->SetBranchAddress("HLT_Photon175", &HLT_Photon175, &b_HLT_Photon175);
   fChain->SetBranchAddress("HLT_Photon200", &HLT_Photon200, &b_HLT_Photon200);
   fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM, &b_HLT_Photon50_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM, &b_HLT_Photon75_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM, &b_HLT_Photon90_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM, &b_HLT_Photon120_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM, &b_HLT_Photon165_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon90_CaloIdL_PFHT700", &HLT_Photon90_CaloIdL_PFHT700, &b_HLT_Photon90_CaloIdL_PFHT700);
   fChain->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90, &b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90);
   fChain->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95, &b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95);
   fChain->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55, &b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55);
   fChain->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55, &b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55);
   fChain->SetBranchAddress("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", &HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55, &b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55);
   fChain->SetBranchAddress("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", &HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55, &b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_L1_NoOS", &HLT_Dimuon0_Jpsi_L1_NoOS, &b_HLT_Dimuon0_Jpsi_L1_NoOS);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_NoOS", &HLT_Dimuon0_Jpsi_NoVertexing_NoOS, &b_HLT_Dimuon0_Jpsi_NoVertexing_NoOS);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi", &HLT_Dimuon0_Jpsi, &b_HLT_Dimuon0_Jpsi);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing", &HLT_Dimuon0_Jpsi_NoVertexing, &b_HLT_Dimuon0_Jpsi_NoVertexing);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R", &HLT_Dimuon0_Jpsi_L1_4R_0er1p5R, &b_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R", &HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R, &b_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi3p5_Muon2", &HLT_Dimuon0_Jpsi3p5_Muon2, &b_HLT_Dimuon0_Jpsi3p5_Muon2);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5", &HLT_Dimuon0_Upsilon_L1_4p5, &b_HLT_Dimuon0_Upsilon_L1_4p5);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_5", &HLT_Dimuon0_Upsilon_L1_5, &b_HLT_Dimuon0_Upsilon_L1_5);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5NoOS", &HLT_Dimuon0_Upsilon_L1_4p5NoOS, &b_HLT_Dimuon0_Upsilon_L1_4p5NoOS);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5er2p0", &HLT_Dimuon0_Upsilon_L1_4p5er2p0, &b_HLT_Dimuon0_Upsilon_L1_4p5er2p0);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5er2p0M", &HLT_Dimuon0_Upsilon_L1_4p5er2p0M, &b_HLT_Dimuon0_Upsilon_L1_4p5er2p0M);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_NoVertexing", &HLT_Dimuon0_Upsilon_NoVertexing, &b_HLT_Dimuon0_Upsilon_NoVertexing);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_5M", &HLT_Dimuon0_Upsilon_L1_5M, &b_HLT_Dimuon0_Upsilon_L1_5M);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_0er1p5R", &HLT_Dimuon0_LowMass_L1_0er1p5R, &b_HLT_Dimuon0_LowMass_L1_0er1p5R);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_0er1p5", &HLT_Dimuon0_LowMass_L1_0er1p5, &b_HLT_Dimuon0_LowMass_L1_0er1p5);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass", &HLT_Dimuon0_LowMass, &b_HLT_Dimuon0_LowMass);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_4", &HLT_Dimuon0_LowMass_L1_4, &b_HLT_Dimuon0_LowMass_L1_4);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_4R", &HLT_Dimuon0_LowMass_L1_4R, &b_HLT_Dimuon0_LowMass_L1_4R);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_TM530", &HLT_Dimuon0_LowMass_L1_TM530, &b_HLT_Dimuon0_LowMass_L1_TM530);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_L1_TM0", &HLT_Dimuon0_Upsilon_Muon_L1_TM0, &b_HLT_Dimuon0_Upsilon_Muon_L1_TM0);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_NoL1Mass", &HLT_Dimuon0_Upsilon_Muon_NoL1Mass, &b_HLT_Dimuon0_Upsilon_Muon_NoL1Mass);
   fChain->SetBranchAddress("HLT_TripleMu_5_3_3_Mass3p8to60_DZ", &HLT_TripleMu_5_3_3_Mass3p8to60_DZ, &b_HLT_TripleMu_5_3_3_Mass3p8to60_DZ);
   fChain->SetBranchAddress("HLT_TripleMu_10_5_5_DZ", &HLT_TripleMu_10_5_5_DZ, &b_HLT_TripleMu_10_5_5_DZ);
   fChain->SetBranchAddress("HLT_TripleMu_12_10_5", &HLT_TripleMu_12_10_5, &b_HLT_TripleMu_12_10_5);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1);
   fChain->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET50_PFMHT60", &HLT_DoubleMu3_DZ_PFMET50_PFMHT60, &b_HLT_DoubleMu3_DZ_PFMET50_PFMHT60);
   fChain->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET70_PFMHT70", &HLT_DoubleMu3_DZ_PFMET70_PFMHT70, &b_HLT_DoubleMu3_DZ_PFMET70_PFMHT70);
   fChain->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET90_PFMHT90", &HLT_DoubleMu3_DZ_PFMET90_PFMHT90, &b_HLT_DoubleMu3_DZ_PFMET90_PFMHT90);
   fChain->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass", &HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass, &b_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass);
   fChain->SetBranchAddress("HLT_DoubleMu4_Jpsi_Displaced", &HLT_DoubleMu4_Jpsi_Displaced, &b_HLT_DoubleMu4_Jpsi_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_Jpsi_NoVertexing", &HLT_DoubleMu4_Jpsi_NoVertexing, &b_HLT_DoubleMu4_Jpsi_NoVertexing);
   fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrkTrk_Displaced", &HLT_DoubleMu4_JpsiTrkTrk_Displaced, &b_HLT_DoubleMu4_JpsiTrkTrk_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu43NoFiltersNoVtx", &HLT_DoubleMu43NoFiltersNoVtx, &b_HLT_DoubleMu43NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_DoubleMu48NoFiltersNoVtx", &HLT_DoubleMu48NoFiltersNoVtx, &b_HLT_DoubleMu48NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL", &HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL, &b_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL", &HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL, &b_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL);
   fChain->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_L1_DM4", &HLT_DoubleMu20_7_Mass0to30_L1_DM4, &b_HLT_DoubleMu20_7_Mass0to30_L1_DM4);
   fChain->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_L1_DM4EG", &HLT_DoubleMu20_7_Mass0to30_L1_DM4EG, &b_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG);
   fChain->SetBranchAddress("HLT_HT425", &HLT_HT425, &b_HLT_HT425);
   fChain->SetBranchAddress("HLT_HT430_DisplacedDijet40_DisplacedTrack", &HLT_HT430_DisplacedDijet40_DisplacedTrack, &b_HLT_HT430_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT430_DisplacedDijet60_DisplacedTrack", &HLT_HT430_DisplacedDijet60_DisplacedTrack, &b_HLT_HT430_DisplacedDijet60_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT430_DisplacedDijet80_DisplacedTrack", &HLT_HT430_DisplacedDijet80_DisplacedTrack, &b_HLT_HT430_DisplacedDijet80_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT400_DisplacedDijet40_DisplacedTrack", &HLT_HT400_DisplacedDijet40_DisplacedTrack, &b_HLT_HT400_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT650_DisplacedDijet60_Inclusive", &HLT_HT650_DisplacedDijet60_Inclusive, &b_HLT_HT650_DisplacedDijet60_Inclusive);
   fChain->SetBranchAddress("HLT_HT550_DisplacedDijet80_Inclusive", &HLT_HT550_DisplacedDijet80_Inclusive, &b_HLT_HT550_DisplacedDijet80_Inclusive);
   fChain->SetBranchAddress("HLT_HT550_DisplacedDijet60_Inclusive", &HLT_HT550_DisplacedDijet60_Inclusive, &b_HLT_HT550_DisplacedDijet60_Inclusive);
   fChain->SetBranchAddress("HLT_HT650_DisplacedDijet80_Inclusive", &HLT_HT650_DisplacedDijet80_Inclusive, &b_HLT_HT650_DisplacedDijet80_Inclusive);
   fChain->SetBranchAddress("HLT_HT750_DisplacedDijet80_Inclusive", &HLT_HT750_DisplacedDijet80_Inclusive, &b_HLT_HT750_DisplacedDijet80_Inclusive);
   fChain->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET110", &HLT_DiJet110_35_Mjj650_PFMET110, &b_HLT_DiJet110_35_Mjj650_PFMET110);
   fChain->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET120", &HLT_DiJet110_35_Mjj650_PFMET120, &b_HLT_DiJet110_35_Mjj650_PFMET120);
   fChain->SetBranchAddress("HLT_DiJet110_35_Mjj650_PFMET130", &HLT_DiJet110_35_Mjj650_PFMET130, &b_HLT_DiJet110_35_Mjj650_PFMET130);
   fChain->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET110", &HLT_TripleJet110_35_35_Mjj650_PFMET110, &b_HLT_TripleJet110_35_35_Mjj650_PFMET110);
   fChain->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET120", &HLT_TripleJet110_35_35_Mjj650_PFMET120, &b_HLT_TripleJet110_35_35_Mjj650_PFMET120);
   fChain->SetBranchAddress("HLT_TripleJet110_35_35_Mjj650_PFMET130", &HLT_TripleJet110_35_35_Mjj650_PFMET130, &b_HLT_TripleJet110_35_35_Mjj650_PFMET130);
   fChain->SetBranchAddress("HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg", &HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg, &b_HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg", &HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg, &b_HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg", &HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg, &b_HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned", &HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned, &b_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned);
   fChain->SetBranchAddress("HLT_Ele28_eta2p1_WPTight_Gsf_HT150", &HLT_Ele28_eta2p1_WPTight_Gsf_HT150, &b_HLT_Ele28_eta2p1_WPTight_Gsf_HT150);
   fChain->SetBranchAddress("HLT_Ele28_HighEta_SC20_Mass55", &HLT_Ele28_HighEta_SC20_Mass55, &b_HLT_Ele28_HighEta_SC20_Mass55);
   fChain->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_Photon23", &HLT_DoubleMu20_7_Mass0to30_Photon23, &b_HLT_DoubleMu20_7_Mass0to30_Photon23);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5", &HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5, &b_HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450_PFMET50", &HLT_Ele15_IsoVVVL_PFHT450_PFMET50, &b_HLT_Ele15_IsoVVVL_PFHT450_PFMET50);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450", &HLT_Ele15_IsoVVVL_PFHT450, &b_HLT_Ele15_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Ele50_IsoVVVL_PFHT450", &HLT_Ele50_IsoVVVL_PFHT450, &b_HLT_Ele50_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT600", &HLT_Ele15_IsoVVVL_PFHT600, &b_HLT_Ele15_IsoVVVL_PFHT600);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60, &b_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60);
   fChain->SetBranchAddress("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60", &HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60, &b_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5", &HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5, &b_HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450_PFMET50", &HLT_Mu15_IsoVVVL_PFHT450_PFMET50, &b_HLT_Mu15_IsoVVVL_PFHT450_PFMET50);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450", &HLT_Mu15_IsoVVVL_PFHT450, &b_HLT_Mu15_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Mu50_IsoVVVL_PFHT450", &HLT_Mu50_IsoVVVL_PFHT450, &b_HLT_Mu50_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT600", &HLT_Mu15_IsoVVVL_PFHT600, &b_HLT_Mu15_IsoVVVL_PFHT600);
   fChain->SetBranchAddress("HLT_Dimuon10_PsiPrime_Barrel_Seagulls", &HLT_Dimuon10_PsiPrime_Barrel_Seagulls, &b_HLT_Dimuon10_PsiPrime_Barrel_Seagulls);
   fChain->SetBranchAddress("HLT_Dimuon20_Jpsi_Barrel_Seagulls", &HLT_Dimuon20_Jpsi_Barrel_Seagulls, &b_HLT_Dimuon20_Jpsi_Barrel_Seagulls);
   fChain->SetBranchAddress("HLT_Dimuon10_Upsilon_Barrel_Seagulls", &HLT_Dimuon10_Upsilon_Barrel_Seagulls, &b_HLT_Dimuon10_Upsilon_Barrel_Seagulls);
   fChain->SetBranchAddress("HLT_Dimuon12_Upsilon_eta1p5", &HLT_Dimuon12_Upsilon_eta1p5, &b_HLT_Dimuon12_Upsilon_eta1p5);
   fChain->SetBranchAddress("HLT_Dimuon14_Phi_Barrel_Seagulls", &HLT_Dimuon14_Phi_Barrel_Seagulls, &b_HLT_Dimuon14_Phi_Barrel_Seagulls);
   fChain->SetBranchAddress("HLT_Dimuon18_PsiPrime", &HLT_Dimuon18_PsiPrime, &b_HLT_Dimuon18_PsiPrime);
   fChain->SetBranchAddress("HLT_Dimuon25_Jpsi", &HLT_Dimuon25_Jpsi, &b_HLT_Dimuon25_Jpsi);
   fChain->SetBranchAddress("HLT_Dimuon18_PsiPrime_noCorrL1", &HLT_Dimuon18_PsiPrime_noCorrL1, &b_HLT_Dimuon18_PsiPrime_noCorrL1);
   fChain->SetBranchAddress("HLT_Dimuon24_Upsilon_noCorrL1", &HLT_Dimuon24_Upsilon_noCorrL1, &b_HLT_Dimuon24_Upsilon_noCorrL1);
   fChain->SetBranchAddress("HLT_Dimuon24_Phi_noCorrL1", &HLT_Dimuon24_Phi_noCorrL1, &b_HLT_Dimuon24_Phi_noCorrL1);
   fChain->SetBranchAddress("HLT_Dimuon25_Jpsi_noCorrL1", &HLT_Dimuon25_Jpsi_noCorrL1, &b_HLT_Dimuon25_Jpsi_noCorrL1);
   fChain->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ, &b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ);
   fChain->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL, &b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_DoubleIsoMu20_eta2p1", &HLT_DoubleIsoMu20_eta2p1, &b_HLT_DoubleIsoMu20_eta2p1);
   fChain->SetBranchAddress("HLT_DoubleIsoMu24_eta2p1", &HLT_DoubleIsoMu24_eta2p1, &b_HLT_DoubleIsoMu24_eta2p1);
   fChain->SetBranchAddress("HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx", &HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx, &b_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx", &HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx, &b_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx", &HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx, &b_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_Mu8", &HLT_Mu8, &b_HLT_Mu8);
   fChain->SetBranchAddress("HLT_Mu17", &HLT_Mu17, &b_HLT_Mu17);
   fChain->SetBranchAddress("HLT_Mu19", &HLT_Mu19, &b_HLT_Mu19);
   fChain->SetBranchAddress("HLT_Mu17_Photon30_IsoCaloId", &HLT_Mu17_Photon30_IsoCaloId, &b_HLT_Mu17_Photon30_IsoCaloId);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &HLT_Ele8_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &HLT_Ele17_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele17_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &HLT_Ele23_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
   fChain->SetBranchAddress("HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_Ele115_CaloIdVT_GsfTrkIdT, &b_HLT_Ele115_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele135_CaloIdVT_GsfTrkIdT", &HLT_Ele135_CaloIdVT_GsfTrkIdT, &b_HLT_Ele135_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele145_CaloIdVT_GsfTrkIdT", &HLT_Ele145_CaloIdVT_GsfTrkIdT, &b_HLT_Ele145_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele200_CaloIdVT_GsfTrkIdT", &HLT_Ele200_CaloIdVT_GsfTrkIdT, &b_HLT_Ele200_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele250_CaloIdVT_GsfTrkIdT", &HLT_Ele250_CaloIdVT_GsfTrkIdT, &b_HLT_Ele250_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele300_CaloIdVT_GsfTrkIdT", &HLT_Ele300_CaloIdVT_GsfTrkIdT, &b_HLT_Ele300_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_PFHT300PT30_QuadPFJet_75_60_45_40", &HLT_PFHT300PT30_QuadPFJet_75_60_45_40, &b_HLT_PFHT300PT30_QuadPFJet_75_60_45_40);
   fChain->SetBranchAddress("HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0", &HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0, &b_HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0);
   fChain->SetBranchAddress("HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2", &HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2, &b_HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2);
   fChain->SetBranchAddress("HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2", &HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2, &b_HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2);
   fChain->SetBranchAddress("HLT_PFHT380_SixPFJet32", &HLT_PFHT380_SixPFJet32, &b_HLT_PFHT380_SixPFJet32);
   fChain->SetBranchAddress("HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5", &HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5, &b_HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5);
   fChain->SetBranchAddress("HLT_PFHT430_SixPFJet40", &HLT_PFHT430_SixPFJet40, &b_HLT_PFHT430_SixPFJet40);
   fChain->SetBranchAddress("HLT_PFHT350", &HLT_PFHT350, &b_HLT_PFHT350);
   fChain->SetBranchAddress("HLT_PFHT350MinPFJet15", &HLT_PFHT350MinPFJet15, &b_HLT_PFHT350MinPFJet15);
   fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL", &HLT_Photon60_R9Id90_CaloIdL_IsoL, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL);
   fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL);
   fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15);
   fChain->SetBranchAddress("HLT_FullTrack_Multiplicity85", &HLT_FullTrack_Multiplicity85, &b_HLT_FullTrack_Multiplicity85);
   fChain->SetBranchAddress("HLT_FullTrack_Multiplicity100", &HLT_FullTrack_Multiplicity100, &b_HLT_FullTrack_Multiplicity100);
   fChain->SetBranchAddress("HLT_FullTrack_Multiplicity130", &HLT_FullTrack_Multiplicity130, &b_HLT_FullTrack_Multiplicity130);
   fChain->SetBranchAddress("HLT_FullTrack_Multiplicity155", &HLT_FullTrack_Multiplicity155, &b_HLT_FullTrack_Multiplicity155);
   fChain->SetBranchAddress("HLT_ECALHT800", &HLT_ECALHT800, &b_HLT_ECALHT800);
   fChain->SetBranchAddress("HLT_DiSC30_18_EIso_AND_HE_Mass70", &HLT_DiSC30_18_EIso_AND_HE_Mass70, &b_HLT_DiSC30_18_EIso_AND_HE_Mass70);
   fChain->SetBranchAddress("HLT_Physics", &HLT_Physics, &b_HLT_Physics);
   fChain->SetBranchAddress("HLT_Physics_part0", &HLT_Physics_part0, &b_HLT_Physics_part0);
   fChain->SetBranchAddress("HLT_Physics_part1", &HLT_Physics_part1, &b_HLT_Physics_part1);
   fChain->SetBranchAddress("HLT_Physics_part2", &HLT_Physics_part2, &b_HLT_Physics_part2);
   fChain->SetBranchAddress("HLT_Physics_part3", &HLT_Physics_part3, &b_HLT_Physics_part3);
   fChain->SetBranchAddress("HLT_Physics_part4", &HLT_Physics_part4, &b_HLT_Physics_part4);
   fChain->SetBranchAddress("HLT_Physics_part5", &HLT_Physics_part5, &b_HLT_Physics_part5);
   fChain->SetBranchAddress("HLT_Physics_part6", &HLT_Physics_part6, &b_HLT_Physics_part6);
   fChain->SetBranchAddress("HLT_Physics_part7", &HLT_Physics_part7, &b_HLT_Physics_part7);
   fChain->SetBranchAddress("HLT_Random", &HLT_Random, &b_HLT_Random);
   fChain->SetBranchAddress("HLT_ZeroBias", &HLT_ZeroBias, &b_HLT_ZeroBias);
   fChain->SetBranchAddress("HLT_ZeroBias_part0", &HLT_ZeroBias_part0, &b_HLT_ZeroBias_part0);
   fChain->SetBranchAddress("HLT_ZeroBias_part1", &HLT_ZeroBias_part1, &b_HLT_ZeroBias_part1);
   fChain->SetBranchAddress("HLT_ZeroBias_part2", &HLT_ZeroBias_part2, &b_HLT_ZeroBias_part2);
   fChain->SetBranchAddress("HLT_ZeroBias_part3", &HLT_ZeroBias_part3, &b_HLT_ZeroBias_part3);
   fChain->SetBranchAddress("HLT_ZeroBias_part4", &HLT_ZeroBias_part4, &b_HLT_ZeroBias_part4);
   fChain->SetBranchAddress("HLT_ZeroBias_part5", &HLT_ZeroBias_part5, &b_HLT_ZeroBias_part5);
   fChain->SetBranchAddress("HLT_ZeroBias_part6", &HLT_ZeroBias_part6, &b_HLT_ZeroBias_part6);
   fChain->SetBranchAddress("HLT_ZeroBias_part7", &HLT_ZeroBias_part7, &b_HLT_ZeroBias_part7);
   fChain->SetBranchAddress("HLT_AK4CaloJet30", &HLT_AK4CaloJet30, &b_HLT_AK4CaloJet30);
   fChain->SetBranchAddress("HLT_AK4CaloJet40", &HLT_AK4CaloJet40, &b_HLT_AK4CaloJet40);
   fChain->SetBranchAddress("HLT_AK4CaloJet50", &HLT_AK4CaloJet50, &b_HLT_AK4CaloJet50);
   fChain->SetBranchAddress("HLT_AK4CaloJet80", &HLT_AK4CaloJet80, &b_HLT_AK4CaloJet80);
   fChain->SetBranchAddress("HLT_AK4CaloJet100", &HLT_AK4CaloJet100, &b_HLT_AK4CaloJet100);
   fChain->SetBranchAddress("HLT_AK4CaloJet120", &HLT_AK4CaloJet120, &b_HLT_AK4CaloJet120);
   fChain->SetBranchAddress("HLT_AK4PFJet30", &HLT_AK4PFJet30, &b_HLT_AK4PFJet30);
   fChain->SetBranchAddress("HLT_AK4PFJet50", &HLT_AK4PFJet50, &b_HLT_AK4PFJet50);
   fChain->SetBranchAddress("HLT_AK4PFJet80", &HLT_AK4PFJet80, &b_HLT_AK4PFJet80);
   fChain->SetBranchAddress("HLT_AK4PFJet100", &HLT_AK4PFJet100, &b_HLT_AK4PFJet100);
   fChain->SetBranchAddress("HLT_AK4PFJet120", &HLT_AK4PFJet120, &b_HLT_AK4PFJet120);
   fChain->SetBranchAddress("HLT_HISinglePhoton10_Eta3p1ForPPRef", &HLT_HISinglePhoton10_Eta3p1ForPPRef, &b_HLT_HISinglePhoton10_Eta3p1ForPPRef);
   fChain->SetBranchAddress("HLT_HISinglePhoton20_Eta3p1ForPPRef", &HLT_HISinglePhoton20_Eta3p1ForPPRef, &b_HLT_HISinglePhoton20_Eta3p1ForPPRef);
   fChain->SetBranchAddress("HLT_HISinglePhoton30_Eta3p1ForPPRef", &HLT_HISinglePhoton30_Eta3p1ForPPRef, &b_HLT_HISinglePhoton30_Eta3p1ForPPRef);
   fChain->SetBranchAddress("HLT_HISinglePhoton40_Eta3p1ForPPRef", &HLT_HISinglePhoton40_Eta3p1ForPPRef, &b_HLT_HISinglePhoton40_Eta3p1ForPPRef);
   fChain->SetBranchAddress("HLT_HISinglePhoton50_Eta3p1ForPPRef", &HLT_HISinglePhoton50_Eta3p1ForPPRef, &b_HLT_HISinglePhoton50_Eta3p1ForPPRef);
   fChain->SetBranchAddress("HLT_HISinglePhoton60_Eta3p1ForPPRef", &HLT_HISinglePhoton60_Eta3p1ForPPRef, &b_HLT_HISinglePhoton60_Eta3p1ForPPRef);
   fChain->SetBranchAddress("HLT_Photon20_HoverELoose", &HLT_Photon20_HoverELoose, &b_HLT_Photon20_HoverELoose);
   fChain->SetBranchAddress("HLT_Photon30_HoverELoose", &HLT_Photon30_HoverELoose, &b_HLT_Photon30_HoverELoose);
   fChain->SetBranchAddress("HLT_Photon40_HoverELoose", &HLT_Photon40_HoverELoose, &b_HLT_Photon40_HoverELoose);
   fChain->SetBranchAddress("HLT_Photon50_HoverELoose", &HLT_Photon50_HoverELoose, &b_HLT_Photon50_HoverELoose);
   fChain->SetBranchAddress("HLT_Photon60_HoverELoose", &HLT_Photon60_HoverELoose, &b_HLT_Photon60_HoverELoose);
   fChain->SetBranchAddress("HLT_EcalCalibration", &HLT_EcalCalibration, &b_HLT_EcalCalibration);
   fChain->SetBranchAddress("HLT_HcalCalibration", &HLT_HcalCalibration, &b_HLT_HcalCalibration);
   fChain->SetBranchAddress("HLT_L1UnpairedBunchBptxMinus", &HLT_L1UnpairedBunchBptxMinus, &b_HLT_L1UnpairedBunchBptxMinus);
   fChain->SetBranchAddress("HLT_L1UnpairedBunchBptxPlus", &HLT_L1UnpairedBunchBptxPlus, &b_HLT_L1UnpairedBunchBptxPlus);
   fChain->SetBranchAddress("HLT_L1NotBptxOR", &HLT_L1NotBptxOR, &b_HLT_L1NotBptxOR);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF_OR", &HLT_L1MinimumBiasHF_OR, &b_HLT_L1MinimumBiasHF_OR);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF0OR", &HLT_L1MinimumBiasHF0OR, &b_HLT_L1MinimumBiasHF0OR);
   fChain->SetBranchAddress("HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142, &b_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142);
   fChain->SetBranchAddress("HLT_HcalNZS", &HLT_HcalNZS, &b_HLT_HcalNZS);
   fChain->SetBranchAddress("HLT_HcalPhiSym", &HLT_HcalPhiSym, &b_HLT_HcalPhiSym);
   fChain->SetBranchAddress("HLT_HcalIsolatedbunch", &HLT_HcalIsolatedbunch, &b_HLT_HcalIsolatedbunch);
   fChain->SetBranchAddress("HLT_IsoTrackHB", &HLT_IsoTrackHB, &b_HLT_IsoTrackHB);
   fChain->SetBranchAddress("HLT_IsoTrackHE", &HLT_IsoTrackHE, &b_HLT_IsoTrackHE);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap", &HLT_ZeroBias_FirstCollisionAfterAbortGap, &b_HLT_ZeroBias_FirstCollisionAfterAbortGap);
   fChain->SetBranchAddress("HLT_ZeroBias_IsolatedBunches", &HLT_ZeroBias_IsolatedBunches, &b_HLT_ZeroBias_IsolatedBunches);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionInTrain", &HLT_ZeroBias_FirstCollisionInTrain, &b_HLT_ZeroBias_FirstCollisionInTrain);
   fChain->SetBranchAddress("HLT_ZeroBias_LastCollisionInTrain", &HLT_ZeroBias_LastCollisionInTrain, &b_HLT_ZeroBias_LastCollisionInTrain);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstBXAfterTrain", &HLT_ZeroBias_FirstBXAfterTrain, &b_HLT_ZeroBias_FirstBXAfterTrain);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg, &b_HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg, &b_HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg, &b_HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg, &b_HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg, &b_HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg, &b_HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg", &HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg", &HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg", &HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg", &HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg, &b_HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", &HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr, &b_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr", &HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr, &b_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr);
   fChain->SetBranchAddress("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1", &HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1, &b_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1", &HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1, &b_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL, &b_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_Rsq0p35", &HLT_Rsq0p35, &b_HLT_Rsq0p35);
   fChain->SetBranchAddress("HLT_Rsq0p40", &HLT_Rsq0p40, &b_HLT_Rsq0p40);
   fChain->SetBranchAddress("HLT_RsqMR300_Rsq0p09_MR200", &HLT_RsqMR300_Rsq0p09_MR200, &b_HLT_RsqMR300_Rsq0p09_MR200);
   fChain->SetBranchAddress("HLT_RsqMR320_Rsq0p09_MR200", &HLT_RsqMR320_Rsq0p09_MR200, &b_HLT_RsqMR320_Rsq0p09_MR200);
   fChain->SetBranchAddress("HLT_RsqMR300_Rsq0p09_MR200_4jet", &HLT_RsqMR300_Rsq0p09_MR200_4jet, &b_HLT_RsqMR300_Rsq0p09_MR200_4jet);
   fChain->SetBranchAddress("HLT_RsqMR320_Rsq0p09_MR200_4jet", &HLT_RsqMR320_Rsq0p09_MR200_4jet, &b_HLT_RsqMR320_Rsq0p09_MR200_4jet);
   fChain->SetBranchAddress("HLT_L1_DoubleJet30_Mass_Min400_Mu10", &HLT_L1_DoubleJet30_Mass_Min400_Mu10, &b_HLT_L1_DoubleJet30_Mass_Min400_Mu10);
   fChain->SetBranchAddress("HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1", &HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1", &HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1", &HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1, &b_HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50", &HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50, &b_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3, &b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3, &b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3);
   fChain->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_PFHT60", &HLT_PFMET100_PFMHT100_IDTight_PFHT60, &b_HLT_PFMET100_PFMHT100_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", &HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60, &b_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60", &HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60, &b_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_Mu18_Mu9_SameSign", &HLT_Mu18_Mu9_SameSign, &b_HLT_Mu18_Mu9_SameSign);
   fChain->SetBranchAddress("HLT_Mu18_Mu9_SameSign_DZ", &HLT_Mu18_Mu9_SameSign_DZ, &b_HLT_Mu18_Mu9_SameSign_DZ);
   fChain->SetBranchAddress("HLT_Mu18_Mu9", &HLT_Mu18_Mu9, &b_HLT_Mu18_Mu9);
   fChain->SetBranchAddress("HLT_Mu18_Mu9_DZ", &HLT_Mu18_Mu9_DZ, &b_HLT_Mu18_Mu9_DZ);
   fChain->SetBranchAddress("HLT_Mu20_Mu10_SameSign", &HLT_Mu20_Mu10_SameSign, &b_HLT_Mu20_Mu10_SameSign);
   fChain->SetBranchAddress("HLT_Mu20_Mu10_SameSign_DZ", &HLT_Mu20_Mu10_SameSign_DZ, &b_HLT_Mu20_Mu10_SameSign_DZ);
   fChain->SetBranchAddress("HLT_Mu20_Mu10", &HLT_Mu20_Mu10, &b_HLT_Mu20_Mu10);
   fChain->SetBranchAddress("HLT_Mu20_Mu10_DZ", &HLT_Mu20_Mu10_DZ, &b_HLT_Mu20_Mu10_DZ);
   fChain->SetBranchAddress("HLT_Mu23_Mu12_SameSign", &HLT_Mu23_Mu12_SameSign, &b_HLT_Mu23_Mu12_SameSign);
   fChain->SetBranchAddress("HLT_Mu23_Mu12_SameSign_DZ", &HLT_Mu23_Mu12_SameSign_DZ, &b_HLT_Mu23_Mu12_SameSign_DZ);
   fChain->SetBranchAddress("HLT_Mu23_Mu12", &HLT_Mu23_Mu12, &b_HLT_Mu23_Mu12);
   fChain->SetBranchAddress("HLT_Mu23_Mu12_DZ", &HLT_Mu23_Mu12_DZ, &b_HLT_Mu23_Mu12_DZ);
   fChain->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi", &HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi, &b_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi);
   fChain->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", &HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi, &b_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi);
   fChain->SetBranchAddress("HLT_DoubleMu3_DCA_PFMET50_PFMHT60", &HLT_DoubleMu3_DCA_PFMET50_PFMHT60, &b_HLT_DoubleMu3_DCA_PFMET50_PFMHT60);
   fChain->SetBranchAddress("HLT_TripleMu_5_3_3_Mass3p8to60_DCA", &HLT_TripleMu_5_3_3_Mass3p8to60_DCA, &b_HLT_TripleMu_5_3_3_Mass3p8to60_DCA);
   fChain->SetBranchAddress("HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1", &HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1, &b_HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1", &HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1, &b_HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1", &HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1, &b_HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1", &HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1, &b_HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2", &HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2, &b_HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2", &HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2, &b_HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2", &HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2, &b_HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2", &HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2, &b_HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet98_83_71_15", &HLT_QuadPFJet98_83_71_15, &b_HLT_QuadPFJet98_83_71_15);
   fChain->SetBranchAddress("HLT_QuadPFJet103_88_75_15", &HLT_QuadPFJet103_88_75_15, &b_HLT_QuadPFJet103_88_75_15);
   fChain->SetBranchAddress("HLT_QuadPFJet105_88_76_15", &HLT_QuadPFJet105_88_76_15, &b_HLT_QuadPFJet105_88_76_15);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_15", &HLT_QuadPFJet111_90_80_15, &b_HLT_QuadPFJet111_90_80_15);
   fChain->SetBranchAddress("HLT_AK8PFJet330_PFAK8BTagCSV_p17", &HLT_AK8PFJet330_PFAK8BTagCSV_p17, &b_HLT_AK8PFJet330_PFAK8BTagCSV_p17);
   fChain->SetBranchAddress("HLT_AK8PFJet330_PFAK8BTagCSV_p1", &HLT_AK8PFJet330_PFAK8BTagCSV_p1, &b_HLT_AK8PFJet330_PFAK8BTagCSV_p1);
   fChain->SetBranchAddress("HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", &HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55, &b_HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55);
   fChain->SetBranchAddress("HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", &HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55, &b_HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55);
   fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   fChain->SetBranchAddress("L1simulation_step", &L1simulation_step, &b_L1simulation_step);
   fChain->SetBranchAddress("Jet_pt_raw", Jet_pt_raw, &b_Jet_pt_raw);
   fChain->SetBranchAddress("Jet_pt_nom", Jet_pt_nom, &b_Jet_pt_nom);
   fChain->SetBranchAddress("Jet_mass_raw", Jet_mass_raw, &b_Jet_mass_raw);
   fChain->SetBranchAddress("Jet_mass_nom", Jet_mass_nom, &b_Jet_mass_nom);
   fChain->SetBranchAddress("Jet_corr_JEC", Jet_corr_JEC, &b_Jet_corr_JEC);
   fChain->SetBranchAddress("Jet_corr_JER", Jet_corr_JER, &b_Jet_corr_JER);
   fChain->SetBranchAddress("MET_T1_pt", &MET_T1_pt, &b_MET_T1_pt);
   fChain->SetBranchAddress("MET_T1_phi", &MET_T1_phi, &b_MET_T1_phi);
   fChain->SetBranchAddress("MET_T1Smear_pt", &MET_T1Smear_pt, &b_MET_T1Smear_pt);
   fChain->SetBranchAddress("MET_T1Smear_phi", &MET_T1Smear_phi, &b_MET_T1Smear_phi);
   fChain->SetBranchAddress("Jet_pt_jer0Up", Jet_pt_jer0Up, &b_Jet_pt_jer0Up);
   fChain->SetBranchAddress("Jet_mass_jer0Up", Jet_mass_jer0Up, &b_Jet_mass_jer0Up);
   fChain->SetBranchAddress("MET_T1_pt_jer0Up", &MET_T1_pt_jer0Up, &b_MET_T1_pt_jer0Up);
   fChain->SetBranchAddress("MET_T1_phi_jer0Up", &MET_T1_phi_jer0Up, &b_MET_T1_phi_jer0Up);
   fChain->SetBranchAddress("MET_T1Smear_pt_jer0Up", &MET_T1Smear_pt_jer0Up, &b_MET_T1Smear_pt_jer0Up);
   fChain->SetBranchAddress("MET_T1Smear_phi_jer0Up", &MET_T1Smear_phi_jer0Up, &b_MET_T1Smear_phi_jer0Up);
   fChain->SetBranchAddress("Jet_pt_jer1Up", Jet_pt_jer1Up, &b_Jet_pt_jer1Up);
   fChain->SetBranchAddress("Jet_mass_jer1Up", Jet_mass_jer1Up, &b_Jet_mass_jer1Up);
   fChain->SetBranchAddress("MET_T1_pt_jer1Up", &MET_T1_pt_jer1Up, &b_MET_T1_pt_jer1Up);
   fChain->SetBranchAddress("MET_T1_phi_jer1Up", &MET_T1_phi_jer1Up, &b_MET_T1_phi_jer1Up);
   fChain->SetBranchAddress("MET_T1Smear_pt_jer1Up", &MET_T1Smear_pt_jer1Up, &b_MET_T1Smear_pt_jer1Up);
   fChain->SetBranchAddress("MET_T1Smear_phi_jer1Up", &MET_T1Smear_phi_jer1Up, &b_MET_T1Smear_phi_jer1Up);
   fChain->SetBranchAddress("Jet_pt_jer2Up", Jet_pt_jer2Up, &b_Jet_pt_jer2Up);
   fChain->SetBranchAddress("Jet_mass_jer2Up", Jet_mass_jer2Up, &b_Jet_mass_jer2Up);
   fChain->SetBranchAddress("MET_T1_pt_jer2Up", &MET_T1_pt_jer2Up, &b_MET_T1_pt_jer2Up);
   fChain->SetBranchAddress("MET_T1_phi_jer2Up", &MET_T1_phi_jer2Up, &b_MET_T1_phi_jer2Up);
   fChain->SetBranchAddress("MET_T1Smear_pt_jer2Up", &MET_T1Smear_pt_jer2Up, &b_MET_T1Smear_pt_jer2Up);
   fChain->SetBranchAddress("MET_T1Smear_phi_jer2Up", &MET_T1Smear_phi_jer2Up, &b_MET_T1Smear_phi_jer2Up);
   fChain->SetBranchAddress("Jet_pt_jer3Up", Jet_pt_jer3Up, &b_Jet_pt_jer3Up);
   fChain->SetBranchAddress("Jet_mass_jer3Up", Jet_mass_jer3Up, &b_Jet_mass_jer3Up);
   fChain->SetBranchAddress("MET_T1_pt_jer3Up", &MET_T1_pt_jer3Up, &b_MET_T1_pt_jer3Up);
   fChain->SetBranchAddress("MET_T1_phi_jer3Up", &MET_T1_phi_jer3Up, &b_MET_T1_phi_jer3Up);
   fChain->SetBranchAddress("MET_T1Smear_pt_jer3Up", &MET_T1Smear_pt_jer3Up, &b_MET_T1Smear_pt_jer3Up);
   fChain->SetBranchAddress("MET_T1Smear_phi_jer3Up", &MET_T1Smear_phi_jer3Up, &b_MET_T1Smear_phi_jer3Up);
   fChain->SetBranchAddress("Jet_pt_jer4Up", Jet_pt_jer4Up, &b_Jet_pt_jer4Up);
   fChain->SetBranchAddress("Jet_mass_jer4Up", Jet_mass_jer4Up, &b_Jet_mass_jer4Up);
   fChain->SetBranchAddress("MET_T1_pt_jer4Up", &MET_T1_pt_jer4Up, &b_MET_T1_pt_jer4Up);
   fChain->SetBranchAddress("MET_T1_phi_jer4Up", &MET_T1_phi_jer4Up, &b_MET_T1_phi_jer4Up);
   fChain->SetBranchAddress("MET_T1Smear_pt_jer4Up", &MET_T1Smear_pt_jer4Up, &b_MET_T1Smear_pt_jer4Up);
   fChain->SetBranchAddress("MET_T1Smear_phi_jer4Up", &MET_T1Smear_phi_jer4Up, &b_MET_T1Smear_phi_jer4Up);
   fChain->SetBranchAddress("Jet_pt_jer5Up", Jet_pt_jer5Up, &b_Jet_pt_jer5Up);
   fChain->SetBranchAddress("Jet_mass_jer5Up", Jet_mass_jer5Up, &b_Jet_mass_jer5Up);
   fChain->SetBranchAddress("MET_T1_pt_jer5Up", &MET_T1_pt_jer5Up, &b_MET_T1_pt_jer5Up);
   fChain->SetBranchAddress("MET_T1_phi_jer5Up", &MET_T1_phi_jer5Up, &b_MET_T1_phi_jer5Up);
   fChain->SetBranchAddress("MET_T1Smear_pt_jer5Up", &MET_T1Smear_pt_jer5Up, &b_MET_T1Smear_pt_jer5Up);
   fChain->SetBranchAddress("MET_T1Smear_phi_jer5Up", &MET_T1Smear_phi_jer5Up, &b_MET_T1Smear_phi_jer5Up);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteStatUp", Jet_pt_jesAbsoluteStatUp, &b_Jet_pt_jesAbsoluteStatUp);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteStatUp", Jet_mass_jesAbsoluteStatUp, &b_Jet_mass_jesAbsoluteStatUp);
   fChain->SetBranchAddress("MET_T1_pt_jesAbsoluteStatUp", &MET_T1_pt_jesAbsoluteStatUp, &b_MET_T1_pt_jesAbsoluteStatUp);
   fChain->SetBranchAddress("MET_T1_phi_jesAbsoluteStatUp", &MET_T1_phi_jesAbsoluteStatUp, &b_MET_T1_phi_jesAbsoluteStatUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesAbsoluteStatUp", &MET_T1Smear_pt_jesAbsoluteStatUp, &b_MET_T1Smear_pt_jesAbsoluteStatUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesAbsoluteStatUp", &MET_T1Smear_phi_jesAbsoluteStatUp, &b_MET_T1Smear_phi_jesAbsoluteStatUp);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteScaleUp", Jet_pt_jesAbsoluteScaleUp, &b_Jet_pt_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteScaleUp", Jet_mass_jesAbsoluteScaleUp, &b_Jet_mass_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("MET_T1_pt_jesAbsoluteScaleUp", &MET_T1_pt_jesAbsoluteScaleUp, &b_MET_T1_pt_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("MET_T1_phi_jesAbsoluteScaleUp", &MET_T1_phi_jesAbsoluteScaleUp, &b_MET_T1_phi_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesAbsoluteScaleUp", &MET_T1Smear_pt_jesAbsoluteScaleUp, &b_MET_T1Smear_pt_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesAbsoluteScaleUp", &MET_T1Smear_phi_jesAbsoluteScaleUp, &b_MET_T1Smear_phi_jesAbsoluteScaleUp);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteFlavMapUp", Jet_pt_jesAbsoluteFlavMapUp, &b_Jet_pt_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteFlavMapUp", Jet_mass_jesAbsoluteFlavMapUp, &b_Jet_mass_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("MET_T1_pt_jesAbsoluteFlavMapUp", &MET_T1_pt_jesAbsoluteFlavMapUp, &b_MET_T1_pt_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("MET_T1_phi_jesAbsoluteFlavMapUp", &MET_T1_phi_jesAbsoluteFlavMapUp, &b_MET_T1_phi_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesAbsoluteFlavMapUp", &MET_T1Smear_pt_jesAbsoluteFlavMapUp, &b_MET_T1Smear_pt_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesAbsoluteFlavMapUp", &MET_T1Smear_phi_jesAbsoluteFlavMapUp, &b_MET_T1Smear_phi_jesAbsoluteFlavMapUp);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteMPFBiasUp", Jet_pt_jesAbsoluteMPFBiasUp, &b_Jet_pt_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteMPFBiasUp", Jet_mass_jesAbsoluteMPFBiasUp, &b_Jet_mass_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("MET_T1_pt_jesAbsoluteMPFBiasUp", &MET_T1_pt_jesAbsoluteMPFBiasUp, &b_MET_T1_pt_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("MET_T1_phi_jesAbsoluteMPFBiasUp", &MET_T1_phi_jesAbsoluteMPFBiasUp, &b_MET_T1_phi_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesAbsoluteMPFBiasUp", &MET_T1Smear_pt_jesAbsoluteMPFBiasUp, &b_MET_T1Smear_pt_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesAbsoluteMPFBiasUp", &MET_T1Smear_phi_jesAbsoluteMPFBiasUp, &b_MET_T1Smear_phi_jesAbsoluteMPFBiasUp);
   fChain->SetBranchAddress("Jet_pt_jesFragmentationUp", Jet_pt_jesFragmentationUp, &b_Jet_pt_jesFragmentationUp);
   fChain->SetBranchAddress("Jet_mass_jesFragmentationUp", Jet_mass_jesFragmentationUp, &b_Jet_mass_jesFragmentationUp);
   fChain->SetBranchAddress("MET_T1_pt_jesFragmentationUp", &MET_T1_pt_jesFragmentationUp, &b_MET_T1_pt_jesFragmentationUp);
   fChain->SetBranchAddress("MET_T1_phi_jesFragmentationUp", &MET_T1_phi_jesFragmentationUp, &b_MET_T1_phi_jesFragmentationUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFragmentationUp", &MET_T1Smear_pt_jesFragmentationUp, &b_MET_T1Smear_pt_jesFragmentationUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFragmentationUp", &MET_T1Smear_phi_jesFragmentationUp, &b_MET_T1Smear_phi_jesFragmentationUp);
   fChain->SetBranchAddress("Jet_pt_jesSinglePionECALUp", Jet_pt_jesSinglePionECALUp, &b_Jet_pt_jesSinglePionECALUp);
   fChain->SetBranchAddress("Jet_mass_jesSinglePionECALUp", Jet_mass_jesSinglePionECALUp, &b_Jet_mass_jesSinglePionECALUp);
   fChain->SetBranchAddress("MET_T1_pt_jesSinglePionECALUp", &MET_T1_pt_jesSinglePionECALUp, &b_MET_T1_pt_jesSinglePionECALUp);
   fChain->SetBranchAddress("MET_T1_phi_jesSinglePionECALUp", &MET_T1_phi_jesSinglePionECALUp, &b_MET_T1_phi_jesSinglePionECALUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSinglePionECALUp", &MET_T1Smear_pt_jesSinglePionECALUp, &b_MET_T1Smear_pt_jesSinglePionECALUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSinglePionECALUp", &MET_T1Smear_phi_jesSinglePionECALUp, &b_MET_T1Smear_phi_jesSinglePionECALUp);
   fChain->SetBranchAddress("Jet_pt_jesSinglePionHCALUp", Jet_pt_jesSinglePionHCALUp, &b_Jet_pt_jesSinglePionHCALUp);
   fChain->SetBranchAddress("Jet_mass_jesSinglePionHCALUp", Jet_mass_jesSinglePionHCALUp, &b_Jet_mass_jesSinglePionHCALUp);
   fChain->SetBranchAddress("MET_T1_pt_jesSinglePionHCALUp", &MET_T1_pt_jesSinglePionHCALUp, &b_MET_T1_pt_jesSinglePionHCALUp);
   fChain->SetBranchAddress("MET_T1_phi_jesSinglePionHCALUp", &MET_T1_phi_jesSinglePionHCALUp, &b_MET_T1_phi_jesSinglePionHCALUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSinglePionHCALUp", &MET_T1Smear_pt_jesSinglePionHCALUp, &b_MET_T1Smear_pt_jesSinglePionHCALUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSinglePionHCALUp", &MET_T1Smear_phi_jesSinglePionHCALUp, &b_MET_T1Smear_phi_jesSinglePionHCALUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorQCDUp", Jet_pt_jesFlavorQCDUp, &b_Jet_pt_jesFlavorQCDUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorQCDUp", Jet_mass_jesFlavorQCDUp, &b_Jet_mass_jesFlavorQCDUp);
   fChain->SetBranchAddress("MET_T1_pt_jesFlavorQCDUp", &MET_T1_pt_jesFlavorQCDUp, &b_MET_T1_pt_jesFlavorQCDUp);
   fChain->SetBranchAddress("MET_T1_phi_jesFlavorQCDUp", &MET_T1_phi_jesFlavorQCDUp, &b_MET_T1_phi_jesFlavorQCDUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFlavorQCDUp", &MET_T1Smear_pt_jesFlavorQCDUp, &b_MET_T1Smear_pt_jesFlavorQCDUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFlavorQCDUp", &MET_T1Smear_phi_jesFlavorQCDUp, &b_MET_T1Smear_phi_jesFlavorQCDUp);
   fChain->SetBranchAddress("Jet_pt_jesTimePtEtaUp", Jet_pt_jesTimePtEtaUp, &b_Jet_pt_jesTimePtEtaUp);
   fChain->SetBranchAddress("Jet_mass_jesTimePtEtaUp", Jet_mass_jesTimePtEtaUp, &b_Jet_mass_jesTimePtEtaUp);
   fChain->SetBranchAddress("MET_T1_pt_jesTimePtEtaUp", &MET_T1_pt_jesTimePtEtaUp, &b_MET_T1_pt_jesTimePtEtaUp);
   fChain->SetBranchAddress("MET_T1_phi_jesTimePtEtaUp", &MET_T1_phi_jesTimePtEtaUp, &b_MET_T1_phi_jesTimePtEtaUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTimePtEtaUp", &MET_T1Smear_pt_jesTimePtEtaUp, &b_MET_T1Smear_pt_jesTimePtEtaUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTimePtEtaUp", &MET_T1Smear_phi_jesTimePtEtaUp, &b_MET_T1Smear_phi_jesTimePtEtaUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJEREC1Up", Jet_pt_jesRelativeJEREC1Up, &b_Jet_pt_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJEREC1Up", Jet_mass_jesRelativeJEREC1Up, &b_Jet_mass_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeJEREC1Up", &MET_T1_pt_jesRelativeJEREC1Up, &b_MET_T1_pt_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeJEREC1Up", &MET_T1_phi_jesRelativeJEREC1Up, &b_MET_T1_phi_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeJEREC1Up", &MET_T1Smear_pt_jesRelativeJEREC1Up, &b_MET_T1Smear_pt_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeJEREC1Up", &MET_T1Smear_phi_jesRelativeJEREC1Up, &b_MET_T1Smear_phi_jesRelativeJEREC1Up);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJEREC2Up", Jet_pt_jesRelativeJEREC2Up, &b_Jet_pt_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJEREC2Up", Jet_mass_jesRelativeJEREC2Up, &b_Jet_mass_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeJEREC2Up", &MET_T1_pt_jesRelativeJEREC2Up, &b_MET_T1_pt_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeJEREC2Up", &MET_T1_phi_jesRelativeJEREC2Up, &b_MET_T1_phi_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeJEREC2Up", &MET_T1Smear_pt_jesRelativeJEREC2Up, &b_MET_T1Smear_pt_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeJEREC2Up", &MET_T1Smear_phi_jesRelativeJEREC2Up, &b_MET_T1Smear_phi_jesRelativeJEREC2Up);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJERHFUp", Jet_pt_jesRelativeJERHFUp, &b_Jet_pt_jesRelativeJERHFUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJERHFUp", Jet_mass_jesRelativeJERHFUp, &b_Jet_mass_jesRelativeJERHFUp);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeJERHFUp", &MET_T1_pt_jesRelativeJERHFUp, &b_MET_T1_pt_jesRelativeJERHFUp);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeJERHFUp", &MET_T1_phi_jesRelativeJERHFUp, &b_MET_T1_phi_jesRelativeJERHFUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeJERHFUp", &MET_T1Smear_pt_jesRelativeJERHFUp, &b_MET_T1Smear_pt_jesRelativeJERHFUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeJERHFUp", &MET_T1Smear_phi_jesRelativeJERHFUp, &b_MET_T1Smear_phi_jesRelativeJERHFUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtBBUp", Jet_pt_jesRelativePtBBUp, &b_Jet_pt_jesRelativePtBBUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtBBUp", Jet_mass_jesRelativePtBBUp, &b_Jet_mass_jesRelativePtBBUp);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativePtBBUp", &MET_T1_pt_jesRelativePtBBUp, &b_MET_T1_pt_jesRelativePtBBUp);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativePtBBUp", &MET_T1_phi_jesRelativePtBBUp, &b_MET_T1_phi_jesRelativePtBBUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativePtBBUp", &MET_T1Smear_pt_jesRelativePtBBUp, &b_MET_T1Smear_pt_jesRelativePtBBUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativePtBBUp", &MET_T1Smear_phi_jesRelativePtBBUp, &b_MET_T1Smear_phi_jesRelativePtBBUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtEC1Up", Jet_pt_jesRelativePtEC1Up, &b_Jet_pt_jesRelativePtEC1Up);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtEC1Up", Jet_mass_jesRelativePtEC1Up, &b_Jet_mass_jesRelativePtEC1Up);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativePtEC1Up", &MET_T1_pt_jesRelativePtEC1Up, &b_MET_T1_pt_jesRelativePtEC1Up);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativePtEC1Up", &MET_T1_phi_jesRelativePtEC1Up, &b_MET_T1_phi_jesRelativePtEC1Up);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativePtEC1Up", &MET_T1Smear_pt_jesRelativePtEC1Up, &b_MET_T1Smear_pt_jesRelativePtEC1Up);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativePtEC1Up", &MET_T1Smear_phi_jesRelativePtEC1Up, &b_MET_T1Smear_phi_jesRelativePtEC1Up);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtEC2Up", Jet_pt_jesRelativePtEC2Up, &b_Jet_pt_jesRelativePtEC2Up);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtEC2Up", Jet_mass_jesRelativePtEC2Up, &b_Jet_mass_jesRelativePtEC2Up);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativePtEC2Up", &MET_T1_pt_jesRelativePtEC2Up, &b_MET_T1_pt_jesRelativePtEC2Up);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativePtEC2Up", &MET_T1_phi_jesRelativePtEC2Up, &b_MET_T1_phi_jesRelativePtEC2Up);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativePtEC2Up", &MET_T1Smear_pt_jesRelativePtEC2Up, &b_MET_T1Smear_pt_jesRelativePtEC2Up);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativePtEC2Up", &MET_T1Smear_phi_jesRelativePtEC2Up, &b_MET_T1Smear_phi_jesRelativePtEC2Up);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtHFUp", Jet_pt_jesRelativePtHFUp, &b_Jet_pt_jesRelativePtHFUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtHFUp", Jet_mass_jesRelativePtHFUp, &b_Jet_mass_jesRelativePtHFUp);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativePtHFUp", &MET_T1_pt_jesRelativePtHFUp, &b_MET_T1_pt_jesRelativePtHFUp);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativePtHFUp", &MET_T1_phi_jesRelativePtHFUp, &b_MET_T1_phi_jesRelativePtHFUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativePtHFUp", &MET_T1Smear_pt_jesRelativePtHFUp, &b_MET_T1Smear_pt_jesRelativePtHFUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativePtHFUp", &MET_T1Smear_phi_jesRelativePtHFUp, &b_MET_T1Smear_phi_jesRelativePtHFUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeBalUp", Jet_pt_jesRelativeBalUp, &b_Jet_pt_jesRelativeBalUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeBalUp", Jet_mass_jesRelativeBalUp, &b_Jet_mass_jesRelativeBalUp);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeBalUp", &MET_T1_pt_jesRelativeBalUp, &b_MET_T1_pt_jesRelativeBalUp);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeBalUp", &MET_T1_phi_jesRelativeBalUp, &b_MET_T1_phi_jesRelativeBalUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeBalUp", &MET_T1Smear_pt_jesRelativeBalUp, &b_MET_T1Smear_pt_jesRelativeBalUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeBalUp", &MET_T1Smear_phi_jesRelativeBalUp, &b_MET_T1Smear_phi_jesRelativeBalUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeSampleUp", Jet_pt_jesRelativeSampleUp, &b_Jet_pt_jesRelativeSampleUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeSampleUp", Jet_mass_jesRelativeSampleUp, &b_Jet_mass_jesRelativeSampleUp);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeSampleUp", &MET_T1_pt_jesRelativeSampleUp, &b_MET_T1_pt_jesRelativeSampleUp);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeSampleUp", &MET_T1_phi_jesRelativeSampleUp, &b_MET_T1_phi_jesRelativeSampleUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeSampleUp", &MET_T1Smear_pt_jesRelativeSampleUp, &b_MET_T1Smear_pt_jesRelativeSampleUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeSampleUp", &MET_T1Smear_phi_jesRelativeSampleUp, &b_MET_T1Smear_phi_jesRelativeSampleUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeFSRUp", Jet_pt_jesRelativeFSRUp, &b_Jet_pt_jesRelativeFSRUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeFSRUp", Jet_mass_jesRelativeFSRUp, &b_Jet_mass_jesRelativeFSRUp);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeFSRUp", &MET_T1_pt_jesRelativeFSRUp, &b_MET_T1_pt_jesRelativeFSRUp);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeFSRUp", &MET_T1_phi_jesRelativeFSRUp, &b_MET_T1_phi_jesRelativeFSRUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeFSRUp", &MET_T1Smear_pt_jesRelativeFSRUp, &b_MET_T1Smear_pt_jesRelativeFSRUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeFSRUp", &MET_T1Smear_phi_jesRelativeFSRUp, &b_MET_T1Smear_phi_jesRelativeFSRUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatFSRUp", Jet_pt_jesRelativeStatFSRUp, &b_Jet_pt_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatFSRUp", Jet_mass_jesRelativeStatFSRUp, &b_Jet_mass_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeStatFSRUp", &MET_T1_pt_jesRelativeStatFSRUp, &b_MET_T1_pt_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeStatFSRUp", &MET_T1_phi_jesRelativeStatFSRUp, &b_MET_T1_phi_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeStatFSRUp", &MET_T1Smear_pt_jesRelativeStatFSRUp, &b_MET_T1Smear_pt_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeStatFSRUp", &MET_T1Smear_phi_jesRelativeStatFSRUp, &b_MET_T1Smear_phi_jesRelativeStatFSRUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatECUp", Jet_pt_jesRelativeStatECUp, &b_Jet_pt_jesRelativeStatECUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatECUp", Jet_mass_jesRelativeStatECUp, &b_Jet_mass_jesRelativeStatECUp);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeStatECUp", &MET_T1_pt_jesRelativeStatECUp, &b_MET_T1_pt_jesRelativeStatECUp);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeStatECUp", &MET_T1_phi_jesRelativeStatECUp, &b_MET_T1_phi_jesRelativeStatECUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeStatECUp", &MET_T1Smear_pt_jesRelativeStatECUp, &b_MET_T1Smear_pt_jesRelativeStatECUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeStatECUp", &MET_T1Smear_phi_jesRelativeStatECUp, &b_MET_T1Smear_phi_jesRelativeStatECUp);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatHFUp", Jet_pt_jesRelativeStatHFUp, &b_Jet_pt_jesRelativeStatHFUp);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatHFUp", Jet_mass_jesRelativeStatHFUp, &b_Jet_mass_jesRelativeStatHFUp);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeStatHFUp", &MET_T1_pt_jesRelativeStatHFUp, &b_MET_T1_pt_jesRelativeStatHFUp);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeStatHFUp", &MET_T1_phi_jesRelativeStatHFUp, &b_MET_T1_phi_jesRelativeStatHFUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeStatHFUp", &MET_T1Smear_pt_jesRelativeStatHFUp, &b_MET_T1Smear_pt_jesRelativeStatHFUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeStatHFUp", &MET_T1Smear_phi_jesRelativeStatHFUp, &b_MET_T1Smear_phi_jesRelativeStatHFUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpDataMCUp", Jet_pt_jesPileUpDataMCUp, &b_Jet_pt_jesPileUpDataMCUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpDataMCUp", Jet_mass_jesPileUpDataMCUp, &b_Jet_mass_jesPileUpDataMCUp);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpDataMCUp", &MET_T1_pt_jesPileUpDataMCUp, &b_MET_T1_pt_jesPileUpDataMCUp);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpDataMCUp", &MET_T1_phi_jesPileUpDataMCUp, &b_MET_T1_phi_jesPileUpDataMCUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpDataMCUp", &MET_T1Smear_pt_jesPileUpDataMCUp, &b_MET_T1Smear_pt_jesPileUpDataMCUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpDataMCUp", &MET_T1Smear_phi_jesPileUpDataMCUp, &b_MET_T1Smear_phi_jesPileUpDataMCUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtRefUp", Jet_pt_jesPileUpPtRefUp, &b_Jet_pt_jesPileUpPtRefUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtRefUp", Jet_mass_jesPileUpPtRefUp, &b_Jet_mass_jesPileUpPtRefUp);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpPtRefUp", &MET_T1_pt_jesPileUpPtRefUp, &b_MET_T1_pt_jesPileUpPtRefUp);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpPtRefUp", &MET_T1_phi_jesPileUpPtRefUp, &b_MET_T1_phi_jesPileUpPtRefUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpPtRefUp", &MET_T1Smear_pt_jesPileUpPtRefUp, &b_MET_T1Smear_pt_jesPileUpPtRefUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpPtRefUp", &MET_T1Smear_phi_jesPileUpPtRefUp, &b_MET_T1Smear_phi_jesPileUpPtRefUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtBBUp", Jet_pt_jesPileUpPtBBUp, &b_Jet_pt_jesPileUpPtBBUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtBBUp", Jet_mass_jesPileUpPtBBUp, &b_Jet_mass_jesPileUpPtBBUp);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpPtBBUp", &MET_T1_pt_jesPileUpPtBBUp, &b_MET_T1_pt_jesPileUpPtBBUp);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpPtBBUp", &MET_T1_phi_jesPileUpPtBBUp, &b_MET_T1_phi_jesPileUpPtBBUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpPtBBUp", &MET_T1Smear_pt_jesPileUpPtBBUp, &b_MET_T1Smear_pt_jesPileUpPtBBUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpPtBBUp", &MET_T1Smear_phi_jesPileUpPtBBUp, &b_MET_T1Smear_phi_jesPileUpPtBBUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtEC1Up", Jet_pt_jesPileUpPtEC1Up, &b_Jet_pt_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtEC1Up", Jet_mass_jesPileUpPtEC1Up, &b_Jet_mass_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpPtEC1Up", &MET_T1_pt_jesPileUpPtEC1Up, &b_MET_T1_pt_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpPtEC1Up", &MET_T1_phi_jesPileUpPtEC1Up, &b_MET_T1_phi_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpPtEC1Up", &MET_T1Smear_pt_jesPileUpPtEC1Up, &b_MET_T1Smear_pt_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpPtEC1Up", &MET_T1Smear_phi_jesPileUpPtEC1Up, &b_MET_T1Smear_phi_jesPileUpPtEC1Up);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtEC2Up", Jet_pt_jesPileUpPtEC2Up, &b_Jet_pt_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtEC2Up", Jet_mass_jesPileUpPtEC2Up, &b_Jet_mass_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpPtEC2Up", &MET_T1_pt_jesPileUpPtEC2Up, &b_MET_T1_pt_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpPtEC2Up", &MET_T1_phi_jesPileUpPtEC2Up, &b_MET_T1_phi_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpPtEC2Up", &MET_T1Smear_pt_jesPileUpPtEC2Up, &b_MET_T1Smear_pt_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpPtEC2Up", &MET_T1Smear_phi_jesPileUpPtEC2Up, &b_MET_T1Smear_phi_jesPileUpPtEC2Up);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtHFUp", Jet_pt_jesPileUpPtHFUp, &b_Jet_pt_jesPileUpPtHFUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtHFUp", Jet_mass_jesPileUpPtHFUp, &b_Jet_mass_jesPileUpPtHFUp);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpPtHFUp", &MET_T1_pt_jesPileUpPtHFUp, &b_MET_T1_pt_jesPileUpPtHFUp);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpPtHFUp", &MET_T1_phi_jesPileUpPtHFUp, &b_MET_T1_phi_jesPileUpPtHFUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpPtHFUp", &MET_T1Smear_pt_jesPileUpPtHFUp, &b_MET_T1Smear_pt_jesPileUpPtHFUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpPtHFUp", &MET_T1Smear_phi_jesPileUpPtHFUp, &b_MET_T1Smear_phi_jesPileUpPtHFUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpMuZeroUp", Jet_pt_jesPileUpMuZeroUp, &b_Jet_pt_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpMuZeroUp", Jet_mass_jesPileUpMuZeroUp, &b_Jet_mass_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpMuZeroUp", &MET_T1_pt_jesPileUpMuZeroUp, &b_MET_T1_pt_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpMuZeroUp", &MET_T1_phi_jesPileUpMuZeroUp, &b_MET_T1_phi_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpMuZeroUp", &MET_T1Smear_pt_jesPileUpMuZeroUp, &b_MET_T1Smear_pt_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpMuZeroUp", &MET_T1Smear_phi_jesPileUpMuZeroUp, &b_MET_T1Smear_phi_jesPileUpMuZeroUp);
   fChain->SetBranchAddress("Jet_pt_jesPileUpEnvelopeUp", Jet_pt_jesPileUpEnvelopeUp, &b_Jet_pt_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("Jet_mass_jesPileUpEnvelopeUp", Jet_mass_jesPileUpEnvelopeUp, &b_Jet_mass_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpEnvelopeUp", &MET_T1_pt_jesPileUpEnvelopeUp, &b_MET_T1_pt_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpEnvelopeUp", &MET_T1_phi_jesPileUpEnvelopeUp, &b_MET_T1_phi_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpEnvelopeUp", &MET_T1Smear_pt_jesPileUpEnvelopeUp, &b_MET_T1Smear_pt_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpEnvelopeUp", &MET_T1Smear_phi_jesPileUpEnvelopeUp, &b_MET_T1Smear_phi_jesPileUpEnvelopeUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalPileUpUp", Jet_pt_jesSubTotalPileUpUp, &b_Jet_pt_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalPileUpUp", Jet_mass_jesSubTotalPileUpUp, &b_Jet_mass_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("MET_T1_pt_jesSubTotalPileUpUp", &MET_T1_pt_jesSubTotalPileUpUp, &b_MET_T1_pt_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("MET_T1_phi_jesSubTotalPileUpUp", &MET_T1_phi_jesSubTotalPileUpUp, &b_MET_T1_phi_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSubTotalPileUpUp", &MET_T1Smear_pt_jesSubTotalPileUpUp, &b_MET_T1Smear_pt_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSubTotalPileUpUp", &MET_T1Smear_phi_jesSubTotalPileUpUp, &b_MET_T1Smear_phi_jesSubTotalPileUpUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalRelativeUp", Jet_pt_jesSubTotalRelativeUp, &b_Jet_pt_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalRelativeUp", Jet_mass_jesSubTotalRelativeUp, &b_Jet_mass_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("MET_T1_pt_jesSubTotalRelativeUp", &MET_T1_pt_jesSubTotalRelativeUp, &b_MET_T1_pt_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("MET_T1_phi_jesSubTotalRelativeUp", &MET_T1_phi_jesSubTotalRelativeUp, &b_MET_T1_phi_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSubTotalRelativeUp", &MET_T1Smear_pt_jesSubTotalRelativeUp, &b_MET_T1Smear_pt_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSubTotalRelativeUp", &MET_T1Smear_phi_jesSubTotalRelativeUp, &b_MET_T1Smear_phi_jesSubTotalRelativeUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalPtUp", Jet_pt_jesSubTotalPtUp, &b_Jet_pt_jesSubTotalPtUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalPtUp", Jet_mass_jesSubTotalPtUp, &b_Jet_mass_jesSubTotalPtUp);
   fChain->SetBranchAddress("MET_T1_pt_jesSubTotalPtUp", &MET_T1_pt_jesSubTotalPtUp, &b_MET_T1_pt_jesSubTotalPtUp);
   fChain->SetBranchAddress("MET_T1_phi_jesSubTotalPtUp", &MET_T1_phi_jesSubTotalPtUp, &b_MET_T1_phi_jesSubTotalPtUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSubTotalPtUp", &MET_T1Smear_pt_jesSubTotalPtUp, &b_MET_T1Smear_pt_jesSubTotalPtUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSubTotalPtUp", &MET_T1Smear_phi_jesSubTotalPtUp, &b_MET_T1Smear_phi_jesSubTotalPtUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalScaleUp", Jet_pt_jesSubTotalScaleUp, &b_Jet_pt_jesSubTotalScaleUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalScaleUp", Jet_mass_jesSubTotalScaleUp, &b_Jet_mass_jesSubTotalScaleUp);
   fChain->SetBranchAddress("MET_T1_pt_jesSubTotalScaleUp", &MET_T1_pt_jesSubTotalScaleUp, &b_MET_T1_pt_jesSubTotalScaleUp);
   fChain->SetBranchAddress("MET_T1_phi_jesSubTotalScaleUp", &MET_T1_phi_jesSubTotalScaleUp, &b_MET_T1_phi_jesSubTotalScaleUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSubTotalScaleUp", &MET_T1Smear_pt_jesSubTotalScaleUp, &b_MET_T1Smear_pt_jesSubTotalScaleUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSubTotalScaleUp", &MET_T1Smear_phi_jesSubTotalScaleUp, &b_MET_T1Smear_phi_jesSubTotalScaleUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalAbsoluteUp", Jet_pt_jesSubTotalAbsoluteUp, &b_Jet_pt_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalAbsoluteUp", Jet_mass_jesSubTotalAbsoluteUp, &b_Jet_mass_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("MET_T1_pt_jesSubTotalAbsoluteUp", &MET_T1_pt_jesSubTotalAbsoluteUp, &b_MET_T1_pt_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("MET_T1_phi_jesSubTotalAbsoluteUp", &MET_T1_phi_jesSubTotalAbsoluteUp, &b_MET_T1_phi_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSubTotalAbsoluteUp", &MET_T1Smear_pt_jesSubTotalAbsoluteUp, &b_MET_T1Smear_pt_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSubTotalAbsoluteUp", &MET_T1Smear_phi_jesSubTotalAbsoluteUp, &b_MET_T1Smear_phi_jesSubTotalAbsoluteUp);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalMCUp", Jet_pt_jesSubTotalMCUp, &b_Jet_pt_jesSubTotalMCUp);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalMCUp", Jet_mass_jesSubTotalMCUp, &b_Jet_mass_jesSubTotalMCUp);
   fChain->SetBranchAddress("MET_T1_pt_jesSubTotalMCUp", &MET_T1_pt_jesSubTotalMCUp, &b_MET_T1_pt_jesSubTotalMCUp);
   fChain->SetBranchAddress("MET_T1_phi_jesSubTotalMCUp", &MET_T1_phi_jesSubTotalMCUp, &b_MET_T1_phi_jesSubTotalMCUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSubTotalMCUp", &MET_T1Smear_pt_jesSubTotalMCUp, &b_MET_T1Smear_pt_jesSubTotalMCUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSubTotalMCUp", &MET_T1Smear_phi_jesSubTotalMCUp, &b_MET_T1Smear_phi_jesSubTotalMCUp);
   fChain->SetBranchAddress("Jet_pt_jesTotalUp", Jet_pt_jesTotalUp, &b_Jet_pt_jesTotalUp);
   fChain->SetBranchAddress("Jet_mass_jesTotalUp", Jet_mass_jesTotalUp, &b_Jet_mass_jesTotalUp);
   fChain->SetBranchAddress("MET_T1_pt_jesTotalUp", &MET_T1_pt_jesTotalUp, &b_MET_T1_pt_jesTotalUp);
   fChain->SetBranchAddress("MET_T1_phi_jesTotalUp", &MET_T1_phi_jesTotalUp, &b_MET_T1_phi_jesTotalUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTotalUp", &MET_T1Smear_pt_jesTotalUp, &b_MET_T1Smear_pt_jesTotalUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTotalUp", &MET_T1Smear_phi_jesTotalUp, &b_MET_T1Smear_phi_jesTotalUp);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoFlavorUp", Jet_pt_jesTotalNoFlavorUp, &b_Jet_pt_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoFlavorUp", Jet_mass_jesTotalNoFlavorUp, &b_Jet_mass_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("MET_T1_pt_jesTotalNoFlavorUp", &MET_T1_pt_jesTotalNoFlavorUp, &b_MET_T1_pt_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("MET_T1_phi_jesTotalNoFlavorUp", &MET_T1_phi_jesTotalNoFlavorUp, &b_MET_T1_phi_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTotalNoFlavorUp", &MET_T1Smear_pt_jesTotalNoFlavorUp, &b_MET_T1Smear_pt_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTotalNoFlavorUp", &MET_T1Smear_phi_jesTotalNoFlavorUp, &b_MET_T1Smear_phi_jesTotalNoFlavorUp);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoTimeUp", Jet_pt_jesTotalNoTimeUp, &b_Jet_pt_jesTotalNoTimeUp);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoTimeUp", Jet_mass_jesTotalNoTimeUp, &b_Jet_mass_jesTotalNoTimeUp);
   fChain->SetBranchAddress("MET_T1_pt_jesTotalNoTimeUp", &MET_T1_pt_jesTotalNoTimeUp, &b_MET_T1_pt_jesTotalNoTimeUp);
   fChain->SetBranchAddress("MET_T1_phi_jesTotalNoTimeUp", &MET_T1_phi_jesTotalNoTimeUp, &b_MET_T1_phi_jesTotalNoTimeUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTotalNoTimeUp", &MET_T1Smear_pt_jesTotalNoTimeUp, &b_MET_T1Smear_pt_jesTotalNoTimeUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTotalNoTimeUp", &MET_T1Smear_phi_jesTotalNoTimeUp, &b_MET_T1Smear_phi_jesTotalNoTimeUp);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoFlavorNoTimeUp", Jet_pt_jesTotalNoFlavorNoTimeUp, &b_Jet_pt_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoFlavorNoTimeUp", Jet_mass_jesTotalNoFlavorNoTimeUp, &b_Jet_mass_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("MET_T1_pt_jesTotalNoFlavorNoTimeUp", &MET_T1_pt_jesTotalNoFlavorNoTimeUp, &b_MET_T1_pt_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("MET_T1_phi_jesTotalNoFlavorNoTimeUp", &MET_T1_phi_jesTotalNoFlavorNoTimeUp, &b_MET_T1_phi_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTotalNoFlavorNoTimeUp", &MET_T1Smear_pt_jesTotalNoFlavorNoTimeUp, &b_MET_T1Smear_pt_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTotalNoFlavorNoTimeUp", &MET_T1Smear_phi_jesTotalNoFlavorNoTimeUp, &b_MET_T1Smear_phi_jesTotalNoFlavorNoTimeUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorZJetUp", Jet_pt_jesFlavorZJetUp, &b_Jet_pt_jesFlavorZJetUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorZJetUp", Jet_mass_jesFlavorZJetUp, &b_Jet_mass_jesFlavorZJetUp);
   fChain->SetBranchAddress("MET_T1_pt_jesFlavorZJetUp", &MET_T1_pt_jesFlavorZJetUp, &b_MET_T1_pt_jesFlavorZJetUp);
   fChain->SetBranchAddress("MET_T1_phi_jesFlavorZJetUp", &MET_T1_phi_jesFlavorZJetUp, &b_MET_T1_phi_jesFlavorZJetUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFlavorZJetUp", &MET_T1Smear_pt_jesFlavorZJetUp, &b_MET_T1Smear_pt_jesFlavorZJetUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFlavorZJetUp", &MET_T1Smear_phi_jesFlavorZJetUp, &b_MET_T1Smear_phi_jesFlavorZJetUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPhotonJetUp", Jet_pt_jesFlavorPhotonJetUp, &b_Jet_pt_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPhotonJetUp", Jet_mass_jesFlavorPhotonJetUp, &b_Jet_mass_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("MET_T1_pt_jesFlavorPhotonJetUp", &MET_T1_pt_jesFlavorPhotonJetUp, &b_MET_T1_pt_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("MET_T1_phi_jesFlavorPhotonJetUp", &MET_T1_phi_jesFlavorPhotonJetUp, &b_MET_T1_phi_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFlavorPhotonJetUp", &MET_T1Smear_pt_jesFlavorPhotonJetUp, &b_MET_T1Smear_pt_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFlavorPhotonJetUp", &MET_T1Smear_phi_jesFlavorPhotonJetUp, &b_MET_T1Smear_phi_jesFlavorPhotonJetUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureGluonUp", Jet_pt_jesFlavorPureGluonUp, &b_Jet_pt_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureGluonUp", Jet_mass_jesFlavorPureGluonUp, &b_Jet_mass_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("MET_T1_pt_jesFlavorPureGluonUp", &MET_T1_pt_jesFlavorPureGluonUp, &b_MET_T1_pt_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("MET_T1_phi_jesFlavorPureGluonUp", &MET_T1_phi_jesFlavorPureGluonUp, &b_MET_T1_phi_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFlavorPureGluonUp", &MET_T1Smear_pt_jesFlavorPureGluonUp, &b_MET_T1Smear_pt_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFlavorPureGluonUp", &MET_T1Smear_phi_jesFlavorPureGluonUp, &b_MET_T1Smear_phi_jesFlavorPureGluonUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureQuarkUp", Jet_pt_jesFlavorPureQuarkUp, &b_Jet_pt_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureQuarkUp", Jet_mass_jesFlavorPureQuarkUp, &b_Jet_mass_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("MET_T1_pt_jesFlavorPureQuarkUp", &MET_T1_pt_jesFlavorPureQuarkUp, &b_MET_T1_pt_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("MET_T1_phi_jesFlavorPureQuarkUp", &MET_T1_phi_jesFlavorPureQuarkUp, &b_MET_T1_phi_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFlavorPureQuarkUp", &MET_T1Smear_pt_jesFlavorPureQuarkUp, &b_MET_T1Smear_pt_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFlavorPureQuarkUp", &MET_T1Smear_phi_jesFlavorPureQuarkUp, &b_MET_T1Smear_phi_jesFlavorPureQuarkUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureCharmUp", Jet_pt_jesFlavorPureCharmUp, &b_Jet_pt_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureCharmUp", Jet_mass_jesFlavorPureCharmUp, &b_Jet_mass_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("MET_T1_pt_jesFlavorPureCharmUp", &MET_T1_pt_jesFlavorPureCharmUp, &b_MET_T1_pt_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("MET_T1_phi_jesFlavorPureCharmUp", &MET_T1_phi_jesFlavorPureCharmUp, &b_MET_T1_phi_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFlavorPureCharmUp", &MET_T1Smear_pt_jesFlavorPureCharmUp, &b_MET_T1Smear_pt_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFlavorPureCharmUp", &MET_T1Smear_phi_jesFlavorPureCharmUp, &b_MET_T1Smear_phi_jesFlavorPureCharmUp);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureBottomUp", Jet_pt_jesFlavorPureBottomUp, &b_Jet_pt_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureBottomUp", Jet_mass_jesFlavorPureBottomUp, &b_Jet_mass_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("MET_T1_pt_jesFlavorPureBottomUp", &MET_T1_pt_jesFlavorPureBottomUp, &b_MET_T1_pt_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("MET_T1_phi_jesFlavorPureBottomUp", &MET_T1_phi_jesFlavorPureBottomUp, &b_MET_T1_phi_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFlavorPureBottomUp", &MET_T1Smear_pt_jesFlavorPureBottomUp, &b_MET_T1Smear_pt_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFlavorPureBottomUp", &MET_T1Smear_phi_jesFlavorPureBottomUp, &b_MET_T1Smear_phi_jesFlavorPureBottomUp);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunBUp", Jet_pt_jesTimeRunBUp, &b_Jet_pt_jesTimeRunBUp);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunBUp", Jet_mass_jesTimeRunBUp, &b_Jet_mass_jesTimeRunBUp);
   fChain->SetBranchAddress("MET_T1_pt_jesTimeRunBUp", &MET_T1_pt_jesTimeRunBUp, &b_MET_T1_pt_jesTimeRunBUp);
   fChain->SetBranchAddress("MET_T1_phi_jesTimeRunBUp", &MET_T1_phi_jesTimeRunBUp, &b_MET_T1_phi_jesTimeRunBUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTimeRunBUp", &MET_T1Smear_pt_jesTimeRunBUp, &b_MET_T1Smear_pt_jesTimeRunBUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTimeRunBUp", &MET_T1Smear_phi_jesTimeRunBUp, &b_MET_T1Smear_phi_jesTimeRunBUp);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunCUp", Jet_pt_jesTimeRunCUp, &b_Jet_pt_jesTimeRunCUp);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunCUp", Jet_mass_jesTimeRunCUp, &b_Jet_mass_jesTimeRunCUp);
   fChain->SetBranchAddress("MET_T1_pt_jesTimeRunCUp", &MET_T1_pt_jesTimeRunCUp, &b_MET_T1_pt_jesTimeRunCUp);
   fChain->SetBranchAddress("MET_T1_phi_jesTimeRunCUp", &MET_T1_phi_jesTimeRunCUp, &b_MET_T1_phi_jesTimeRunCUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTimeRunCUp", &MET_T1Smear_pt_jesTimeRunCUp, &b_MET_T1Smear_pt_jesTimeRunCUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTimeRunCUp", &MET_T1Smear_phi_jesTimeRunCUp, &b_MET_T1Smear_phi_jesTimeRunCUp);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunDEUp", Jet_pt_jesTimeRunDEUp, &b_Jet_pt_jesTimeRunDEUp);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunDEUp", Jet_mass_jesTimeRunDEUp, &b_Jet_mass_jesTimeRunDEUp);
   fChain->SetBranchAddress("MET_T1_pt_jesTimeRunDEUp", &MET_T1_pt_jesTimeRunDEUp, &b_MET_T1_pt_jesTimeRunDEUp);
   fChain->SetBranchAddress("MET_T1_phi_jesTimeRunDEUp", &MET_T1_phi_jesTimeRunDEUp, &b_MET_T1_phi_jesTimeRunDEUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTimeRunDEUp", &MET_T1Smear_pt_jesTimeRunDEUp, &b_MET_T1Smear_pt_jesTimeRunDEUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTimeRunDEUp", &MET_T1Smear_phi_jesTimeRunDEUp, &b_MET_T1Smear_phi_jesTimeRunDEUp);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunFUp", Jet_pt_jesTimeRunFUp, &b_Jet_pt_jesTimeRunFUp);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunFUp", Jet_mass_jesTimeRunFUp, &b_Jet_mass_jesTimeRunFUp);
   fChain->SetBranchAddress("MET_T1_pt_jesTimeRunFUp", &MET_T1_pt_jesTimeRunFUp, &b_MET_T1_pt_jesTimeRunFUp);
   fChain->SetBranchAddress("MET_T1_phi_jesTimeRunFUp", &MET_T1_phi_jesTimeRunFUp, &b_MET_T1_phi_jesTimeRunFUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTimeRunFUp", &MET_T1Smear_pt_jesTimeRunFUp, &b_MET_T1Smear_pt_jesTimeRunFUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTimeRunFUp", &MET_T1Smear_phi_jesTimeRunFUp, &b_MET_T1Smear_phi_jesTimeRunFUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupMPFInSituUp", Jet_pt_jesCorrelationGroupMPFInSituUp, &b_Jet_pt_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupMPFInSituUp", Jet_mass_jesCorrelationGroupMPFInSituUp, &b_Jet_mass_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("MET_T1_pt_jesCorrelationGroupMPFInSituUp", &MET_T1_pt_jesCorrelationGroupMPFInSituUp, &b_MET_T1_pt_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("MET_T1_phi_jesCorrelationGroupMPFInSituUp", &MET_T1_phi_jesCorrelationGroupMPFInSituUp, &b_MET_T1_phi_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesCorrelationGroupMPFInSituUp", &MET_T1Smear_pt_jesCorrelationGroupMPFInSituUp, &b_MET_T1Smear_pt_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesCorrelationGroupMPFInSituUp", &MET_T1Smear_phi_jesCorrelationGroupMPFInSituUp, &b_MET_T1Smear_phi_jesCorrelationGroupMPFInSituUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupIntercalibrationUp", Jet_pt_jesCorrelationGroupIntercalibrationUp, &b_Jet_pt_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupIntercalibrationUp", Jet_mass_jesCorrelationGroupIntercalibrationUp, &b_Jet_mass_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("MET_T1_pt_jesCorrelationGroupIntercalibrationUp", &MET_T1_pt_jesCorrelationGroupIntercalibrationUp, &b_MET_T1_pt_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("MET_T1_phi_jesCorrelationGroupIntercalibrationUp", &MET_T1_phi_jesCorrelationGroupIntercalibrationUp, &b_MET_T1_phi_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesCorrelationGroupIntercalibrationUp", &MET_T1Smear_pt_jesCorrelationGroupIntercalibrationUp, &b_MET_T1Smear_pt_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesCorrelationGroupIntercalibrationUp", &MET_T1Smear_phi_jesCorrelationGroupIntercalibrationUp, &b_MET_T1Smear_phi_jesCorrelationGroupIntercalibrationUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupbJESUp", Jet_pt_jesCorrelationGroupbJESUp, &b_Jet_pt_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupbJESUp", Jet_mass_jesCorrelationGroupbJESUp, &b_Jet_mass_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("MET_T1_pt_jesCorrelationGroupbJESUp", &MET_T1_pt_jesCorrelationGroupbJESUp, &b_MET_T1_pt_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("MET_T1_phi_jesCorrelationGroupbJESUp", &MET_T1_phi_jesCorrelationGroupbJESUp, &b_MET_T1_phi_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesCorrelationGroupbJESUp", &MET_T1Smear_pt_jesCorrelationGroupbJESUp, &b_MET_T1Smear_pt_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesCorrelationGroupbJESUp", &MET_T1Smear_phi_jesCorrelationGroupbJESUp, &b_MET_T1Smear_phi_jesCorrelationGroupbJESUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupFlavorUp", Jet_pt_jesCorrelationGroupFlavorUp, &b_Jet_pt_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupFlavorUp", Jet_mass_jesCorrelationGroupFlavorUp, &b_Jet_mass_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("MET_T1_pt_jesCorrelationGroupFlavorUp", &MET_T1_pt_jesCorrelationGroupFlavorUp, &b_MET_T1_pt_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("MET_T1_phi_jesCorrelationGroupFlavorUp", &MET_T1_phi_jesCorrelationGroupFlavorUp, &b_MET_T1_phi_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesCorrelationGroupFlavorUp", &MET_T1Smear_pt_jesCorrelationGroupFlavorUp, &b_MET_T1Smear_pt_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesCorrelationGroupFlavorUp", &MET_T1Smear_phi_jesCorrelationGroupFlavorUp, &b_MET_T1Smear_phi_jesCorrelationGroupFlavorUp);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupUncorrelatedUp", Jet_pt_jesCorrelationGroupUncorrelatedUp, &b_Jet_pt_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupUncorrelatedUp", Jet_mass_jesCorrelationGroupUncorrelatedUp, &b_Jet_mass_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("MET_T1_pt_jesCorrelationGroupUncorrelatedUp", &MET_T1_pt_jesCorrelationGroupUncorrelatedUp, &b_MET_T1_pt_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("MET_T1_phi_jesCorrelationGroupUncorrelatedUp", &MET_T1_phi_jesCorrelationGroupUncorrelatedUp, &b_MET_T1_phi_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesCorrelationGroupUncorrelatedUp", &MET_T1Smear_pt_jesCorrelationGroupUncorrelatedUp, &b_MET_T1Smear_pt_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesCorrelationGroupUncorrelatedUp", &MET_T1Smear_phi_jesCorrelationGroupUncorrelatedUp, &b_MET_T1Smear_phi_jesCorrelationGroupUncorrelatedUp);
   fChain->SetBranchAddress("MET_T1_pt_unclustEnUp", &MET_T1_pt_unclustEnUp, &b_MET_T1_pt_unclustEnUp);
   fChain->SetBranchAddress("MET_T1_phi_unclustEnUp", &MET_T1_phi_unclustEnUp, &b_MET_T1_phi_unclustEnUp);
   fChain->SetBranchAddress("MET_T1Smear_pt_unclustEnUp", &MET_T1Smear_pt_unclustEnUp, &b_MET_T1Smear_pt_unclustEnUp);
   fChain->SetBranchAddress("MET_T1Smear_phi_unclustEnUp", &MET_T1Smear_phi_unclustEnUp, &b_MET_T1Smear_phi_unclustEnUp);
   fChain->SetBranchAddress("Jet_pt_jer0Down", Jet_pt_jer0Down, &b_Jet_pt_jer0Down);
   fChain->SetBranchAddress("Jet_mass_jer0Down", Jet_mass_jer0Down, &b_Jet_mass_jer0Down);
   fChain->SetBranchAddress("MET_T1_pt_jer0Down", &MET_T1_pt_jer0Down, &b_MET_T1_pt_jer0Down);
   fChain->SetBranchAddress("MET_T1_phi_jer0Down", &MET_T1_phi_jer0Down, &b_MET_T1_phi_jer0Down);
   fChain->SetBranchAddress("MET_T1Smear_pt_jer0Down", &MET_T1Smear_pt_jer0Down, &b_MET_T1Smear_pt_jer0Down);
   fChain->SetBranchAddress("MET_T1Smear_phi_jer0Down", &MET_T1Smear_phi_jer0Down, &b_MET_T1Smear_phi_jer0Down);
   fChain->SetBranchAddress("Jet_pt_jer1Down", Jet_pt_jer1Down, &b_Jet_pt_jer1Down);
   fChain->SetBranchAddress("Jet_mass_jer1Down", Jet_mass_jer1Down, &b_Jet_mass_jer1Down);
   fChain->SetBranchAddress("MET_T1_pt_jer1Down", &MET_T1_pt_jer1Down, &b_MET_T1_pt_jer1Down);
   fChain->SetBranchAddress("MET_T1_phi_jer1Down", &MET_T1_phi_jer1Down, &b_MET_T1_phi_jer1Down);
   fChain->SetBranchAddress("MET_T1Smear_pt_jer1Down", &MET_T1Smear_pt_jer1Down, &b_MET_T1Smear_pt_jer1Down);
   fChain->SetBranchAddress("MET_T1Smear_phi_jer1Down", &MET_T1Smear_phi_jer1Down, &b_MET_T1Smear_phi_jer1Down);
   fChain->SetBranchAddress("Jet_pt_jer2Down", Jet_pt_jer2Down, &b_Jet_pt_jer2Down);
   fChain->SetBranchAddress("Jet_mass_jer2Down", Jet_mass_jer2Down, &b_Jet_mass_jer2Down);
   fChain->SetBranchAddress("MET_T1_pt_jer2Down", &MET_T1_pt_jer2Down, &b_MET_T1_pt_jer2Down);
   fChain->SetBranchAddress("MET_T1_phi_jer2Down", &MET_T1_phi_jer2Down, &b_MET_T1_phi_jer2Down);
   fChain->SetBranchAddress("MET_T1Smear_pt_jer2Down", &MET_T1Smear_pt_jer2Down, &b_MET_T1Smear_pt_jer2Down);
   fChain->SetBranchAddress("MET_T1Smear_phi_jer2Down", &MET_T1Smear_phi_jer2Down, &b_MET_T1Smear_phi_jer2Down);
   fChain->SetBranchAddress("Jet_pt_jer3Down", Jet_pt_jer3Down, &b_Jet_pt_jer3Down);
   fChain->SetBranchAddress("Jet_mass_jer3Down", Jet_mass_jer3Down, &b_Jet_mass_jer3Down);
   fChain->SetBranchAddress("MET_T1_pt_jer3Down", &MET_T1_pt_jer3Down, &b_MET_T1_pt_jer3Down);
   fChain->SetBranchAddress("MET_T1_phi_jer3Down", &MET_T1_phi_jer3Down, &b_MET_T1_phi_jer3Down);
   fChain->SetBranchAddress("MET_T1Smear_pt_jer3Down", &MET_T1Smear_pt_jer3Down, &b_MET_T1Smear_pt_jer3Down);
   fChain->SetBranchAddress("MET_T1Smear_phi_jer3Down", &MET_T1Smear_phi_jer3Down, &b_MET_T1Smear_phi_jer3Down);
   fChain->SetBranchAddress("Jet_pt_jer4Down", Jet_pt_jer4Down, &b_Jet_pt_jer4Down);
   fChain->SetBranchAddress("Jet_mass_jer4Down", Jet_mass_jer4Down, &b_Jet_mass_jer4Down);
   fChain->SetBranchAddress("MET_T1_pt_jer4Down", &MET_T1_pt_jer4Down, &b_MET_T1_pt_jer4Down);
   fChain->SetBranchAddress("MET_T1_phi_jer4Down", &MET_T1_phi_jer4Down, &b_MET_T1_phi_jer4Down);
   fChain->SetBranchAddress("MET_T1Smear_pt_jer4Down", &MET_T1Smear_pt_jer4Down, &b_MET_T1Smear_pt_jer4Down);
   fChain->SetBranchAddress("MET_T1Smear_phi_jer4Down", &MET_T1Smear_phi_jer4Down, &b_MET_T1Smear_phi_jer4Down);
   fChain->SetBranchAddress("Jet_pt_jer5Down", Jet_pt_jer5Down, &b_Jet_pt_jer5Down);
   fChain->SetBranchAddress("Jet_mass_jer5Down", Jet_mass_jer5Down, &b_Jet_mass_jer5Down);
   fChain->SetBranchAddress("MET_T1_pt_jer5Down", &MET_T1_pt_jer5Down, &b_MET_T1_pt_jer5Down);
   fChain->SetBranchAddress("MET_T1_phi_jer5Down", &MET_T1_phi_jer5Down, &b_MET_T1_phi_jer5Down);
   fChain->SetBranchAddress("MET_T1Smear_pt_jer5Down", &MET_T1Smear_pt_jer5Down, &b_MET_T1Smear_pt_jer5Down);
   fChain->SetBranchAddress("MET_T1Smear_phi_jer5Down", &MET_T1Smear_phi_jer5Down, &b_MET_T1Smear_phi_jer5Down);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteStatDown", Jet_pt_jesAbsoluteStatDown, &b_Jet_pt_jesAbsoluteStatDown);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteStatDown", Jet_mass_jesAbsoluteStatDown, &b_Jet_mass_jesAbsoluteStatDown);
   fChain->SetBranchAddress("MET_T1_pt_jesAbsoluteStatDown", &MET_T1_pt_jesAbsoluteStatDown, &b_MET_T1_pt_jesAbsoluteStatDown);
   fChain->SetBranchAddress("MET_T1_phi_jesAbsoluteStatDown", &MET_T1_phi_jesAbsoluteStatDown, &b_MET_T1_phi_jesAbsoluteStatDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesAbsoluteStatDown", &MET_T1Smear_pt_jesAbsoluteStatDown, &b_MET_T1Smear_pt_jesAbsoluteStatDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesAbsoluteStatDown", &MET_T1Smear_phi_jesAbsoluteStatDown, &b_MET_T1Smear_phi_jesAbsoluteStatDown);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteScaleDown", Jet_pt_jesAbsoluteScaleDown, &b_Jet_pt_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteScaleDown", Jet_mass_jesAbsoluteScaleDown, &b_Jet_mass_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("MET_T1_pt_jesAbsoluteScaleDown", &MET_T1_pt_jesAbsoluteScaleDown, &b_MET_T1_pt_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("MET_T1_phi_jesAbsoluteScaleDown", &MET_T1_phi_jesAbsoluteScaleDown, &b_MET_T1_phi_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesAbsoluteScaleDown", &MET_T1Smear_pt_jesAbsoluteScaleDown, &b_MET_T1Smear_pt_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesAbsoluteScaleDown", &MET_T1Smear_phi_jesAbsoluteScaleDown, &b_MET_T1Smear_phi_jesAbsoluteScaleDown);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteFlavMapDown", Jet_pt_jesAbsoluteFlavMapDown, &b_Jet_pt_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteFlavMapDown", Jet_mass_jesAbsoluteFlavMapDown, &b_Jet_mass_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("MET_T1_pt_jesAbsoluteFlavMapDown", &MET_T1_pt_jesAbsoluteFlavMapDown, &b_MET_T1_pt_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("MET_T1_phi_jesAbsoluteFlavMapDown", &MET_T1_phi_jesAbsoluteFlavMapDown, &b_MET_T1_phi_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesAbsoluteFlavMapDown", &MET_T1Smear_pt_jesAbsoluteFlavMapDown, &b_MET_T1Smear_pt_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesAbsoluteFlavMapDown", &MET_T1Smear_phi_jesAbsoluteFlavMapDown, &b_MET_T1Smear_phi_jesAbsoluteFlavMapDown);
   fChain->SetBranchAddress("Jet_pt_jesAbsoluteMPFBiasDown", Jet_pt_jesAbsoluteMPFBiasDown, &b_Jet_pt_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("Jet_mass_jesAbsoluteMPFBiasDown", Jet_mass_jesAbsoluteMPFBiasDown, &b_Jet_mass_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("MET_T1_pt_jesAbsoluteMPFBiasDown", &MET_T1_pt_jesAbsoluteMPFBiasDown, &b_MET_T1_pt_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("MET_T1_phi_jesAbsoluteMPFBiasDown", &MET_T1_phi_jesAbsoluteMPFBiasDown, &b_MET_T1_phi_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesAbsoluteMPFBiasDown", &MET_T1Smear_pt_jesAbsoluteMPFBiasDown, &b_MET_T1Smear_pt_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesAbsoluteMPFBiasDown", &MET_T1Smear_phi_jesAbsoluteMPFBiasDown, &b_MET_T1Smear_phi_jesAbsoluteMPFBiasDown);
   fChain->SetBranchAddress("Jet_pt_jesFragmentationDown", Jet_pt_jesFragmentationDown, &b_Jet_pt_jesFragmentationDown);
   fChain->SetBranchAddress("Jet_mass_jesFragmentationDown", Jet_mass_jesFragmentationDown, &b_Jet_mass_jesFragmentationDown);
   fChain->SetBranchAddress("MET_T1_pt_jesFragmentationDown", &MET_T1_pt_jesFragmentationDown, &b_MET_T1_pt_jesFragmentationDown);
   fChain->SetBranchAddress("MET_T1_phi_jesFragmentationDown", &MET_T1_phi_jesFragmentationDown, &b_MET_T1_phi_jesFragmentationDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFragmentationDown", &MET_T1Smear_pt_jesFragmentationDown, &b_MET_T1Smear_pt_jesFragmentationDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFragmentationDown", &MET_T1Smear_phi_jesFragmentationDown, &b_MET_T1Smear_phi_jesFragmentationDown);
   fChain->SetBranchAddress("Jet_pt_jesSinglePionECALDown", Jet_pt_jesSinglePionECALDown, &b_Jet_pt_jesSinglePionECALDown);
   fChain->SetBranchAddress("Jet_mass_jesSinglePionECALDown", Jet_mass_jesSinglePionECALDown, &b_Jet_mass_jesSinglePionECALDown);
   fChain->SetBranchAddress("MET_T1_pt_jesSinglePionECALDown", &MET_T1_pt_jesSinglePionECALDown, &b_MET_T1_pt_jesSinglePionECALDown);
   fChain->SetBranchAddress("MET_T1_phi_jesSinglePionECALDown", &MET_T1_phi_jesSinglePionECALDown, &b_MET_T1_phi_jesSinglePionECALDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSinglePionECALDown", &MET_T1Smear_pt_jesSinglePionECALDown, &b_MET_T1Smear_pt_jesSinglePionECALDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSinglePionECALDown", &MET_T1Smear_phi_jesSinglePionECALDown, &b_MET_T1Smear_phi_jesSinglePionECALDown);
   fChain->SetBranchAddress("Jet_pt_jesSinglePionHCALDown", Jet_pt_jesSinglePionHCALDown, &b_Jet_pt_jesSinglePionHCALDown);
   fChain->SetBranchAddress("Jet_mass_jesSinglePionHCALDown", Jet_mass_jesSinglePionHCALDown, &b_Jet_mass_jesSinglePionHCALDown);
   fChain->SetBranchAddress("MET_T1_pt_jesSinglePionHCALDown", &MET_T1_pt_jesSinglePionHCALDown, &b_MET_T1_pt_jesSinglePionHCALDown);
   fChain->SetBranchAddress("MET_T1_phi_jesSinglePionHCALDown", &MET_T1_phi_jesSinglePionHCALDown, &b_MET_T1_phi_jesSinglePionHCALDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSinglePionHCALDown", &MET_T1Smear_pt_jesSinglePionHCALDown, &b_MET_T1Smear_pt_jesSinglePionHCALDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSinglePionHCALDown", &MET_T1Smear_phi_jesSinglePionHCALDown, &b_MET_T1Smear_phi_jesSinglePionHCALDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorQCDDown", Jet_pt_jesFlavorQCDDown, &b_Jet_pt_jesFlavorQCDDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorQCDDown", Jet_mass_jesFlavorQCDDown, &b_Jet_mass_jesFlavorQCDDown);
   fChain->SetBranchAddress("MET_T1_pt_jesFlavorQCDDown", &MET_T1_pt_jesFlavorQCDDown, &b_MET_T1_pt_jesFlavorQCDDown);
   fChain->SetBranchAddress("MET_T1_phi_jesFlavorQCDDown", &MET_T1_phi_jesFlavorQCDDown, &b_MET_T1_phi_jesFlavorQCDDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFlavorQCDDown", &MET_T1Smear_pt_jesFlavorQCDDown, &b_MET_T1Smear_pt_jesFlavorQCDDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFlavorQCDDown", &MET_T1Smear_phi_jesFlavorQCDDown, &b_MET_T1Smear_phi_jesFlavorQCDDown);
   fChain->SetBranchAddress("Jet_pt_jesTimePtEtaDown", Jet_pt_jesTimePtEtaDown, &b_Jet_pt_jesTimePtEtaDown);
   fChain->SetBranchAddress("Jet_mass_jesTimePtEtaDown", Jet_mass_jesTimePtEtaDown, &b_Jet_mass_jesTimePtEtaDown);
   fChain->SetBranchAddress("MET_T1_pt_jesTimePtEtaDown", &MET_T1_pt_jesTimePtEtaDown, &b_MET_T1_pt_jesTimePtEtaDown);
   fChain->SetBranchAddress("MET_T1_phi_jesTimePtEtaDown", &MET_T1_phi_jesTimePtEtaDown, &b_MET_T1_phi_jesTimePtEtaDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTimePtEtaDown", &MET_T1Smear_pt_jesTimePtEtaDown, &b_MET_T1Smear_pt_jesTimePtEtaDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTimePtEtaDown", &MET_T1Smear_phi_jesTimePtEtaDown, &b_MET_T1Smear_phi_jesTimePtEtaDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJEREC1Down", Jet_pt_jesRelativeJEREC1Down, &b_Jet_pt_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJEREC1Down", Jet_mass_jesRelativeJEREC1Down, &b_Jet_mass_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeJEREC1Down", &MET_T1_pt_jesRelativeJEREC1Down, &b_MET_T1_pt_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeJEREC1Down", &MET_T1_phi_jesRelativeJEREC1Down, &b_MET_T1_phi_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeJEREC1Down", &MET_T1Smear_pt_jesRelativeJEREC1Down, &b_MET_T1Smear_pt_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeJEREC1Down", &MET_T1Smear_phi_jesRelativeJEREC1Down, &b_MET_T1Smear_phi_jesRelativeJEREC1Down);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJEREC2Down", Jet_pt_jesRelativeJEREC2Down, &b_Jet_pt_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJEREC2Down", Jet_mass_jesRelativeJEREC2Down, &b_Jet_mass_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeJEREC2Down", &MET_T1_pt_jesRelativeJEREC2Down, &b_MET_T1_pt_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeJEREC2Down", &MET_T1_phi_jesRelativeJEREC2Down, &b_MET_T1_phi_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeJEREC2Down", &MET_T1Smear_pt_jesRelativeJEREC2Down, &b_MET_T1Smear_pt_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeJEREC2Down", &MET_T1Smear_phi_jesRelativeJEREC2Down, &b_MET_T1Smear_phi_jesRelativeJEREC2Down);
   fChain->SetBranchAddress("Jet_pt_jesRelativeJERHFDown", Jet_pt_jesRelativeJERHFDown, &b_Jet_pt_jesRelativeJERHFDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeJERHFDown", Jet_mass_jesRelativeJERHFDown, &b_Jet_mass_jesRelativeJERHFDown);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeJERHFDown", &MET_T1_pt_jesRelativeJERHFDown, &b_MET_T1_pt_jesRelativeJERHFDown);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeJERHFDown", &MET_T1_phi_jesRelativeJERHFDown, &b_MET_T1_phi_jesRelativeJERHFDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeJERHFDown", &MET_T1Smear_pt_jesRelativeJERHFDown, &b_MET_T1Smear_pt_jesRelativeJERHFDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeJERHFDown", &MET_T1Smear_phi_jesRelativeJERHFDown, &b_MET_T1Smear_phi_jesRelativeJERHFDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtBBDown", Jet_pt_jesRelativePtBBDown, &b_Jet_pt_jesRelativePtBBDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtBBDown", Jet_mass_jesRelativePtBBDown, &b_Jet_mass_jesRelativePtBBDown);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativePtBBDown", &MET_T1_pt_jesRelativePtBBDown, &b_MET_T1_pt_jesRelativePtBBDown);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativePtBBDown", &MET_T1_phi_jesRelativePtBBDown, &b_MET_T1_phi_jesRelativePtBBDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativePtBBDown", &MET_T1Smear_pt_jesRelativePtBBDown, &b_MET_T1Smear_pt_jesRelativePtBBDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativePtBBDown", &MET_T1Smear_phi_jesRelativePtBBDown, &b_MET_T1Smear_phi_jesRelativePtBBDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtEC1Down", Jet_pt_jesRelativePtEC1Down, &b_Jet_pt_jesRelativePtEC1Down);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtEC1Down", Jet_mass_jesRelativePtEC1Down, &b_Jet_mass_jesRelativePtEC1Down);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativePtEC1Down", &MET_T1_pt_jesRelativePtEC1Down, &b_MET_T1_pt_jesRelativePtEC1Down);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativePtEC1Down", &MET_T1_phi_jesRelativePtEC1Down, &b_MET_T1_phi_jesRelativePtEC1Down);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativePtEC1Down", &MET_T1Smear_pt_jesRelativePtEC1Down, &b_MET_T1Smear_pt_jesRelativePtEC1Down);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativePtEC1Down", &MET_T1Smear_phi_jesRelativePtEC1Down, &b_MET_T1Smear_phi_jesRelativePtEC1Down);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtEC2Down", Jet_pt_jesRelativePtEC2Down, &b_Jet_pt_jesRelativePtEC2Down);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtEC2Down", Jet_mass_jesRelativePtEC2Down, &b_Jet_mass_jesRelativePtEC2Down);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativePtEC2Down", &MET_T1_pt_jesRelativePtEC2Down, &b_MET_T1_pt_jesRelativePtEC2Down);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativePtEC2Down", &MET_T1_phi_jesRelativePtEC2Down, &b_MET_T1_phi_jesRelativePtEC2Down);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativePtEC2Down", &MET_T1Smear_pt_jesRelativePtEC2Down, &b_MET_T1Smear_pt_jesRelativePtEC2Down);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativePtEC2Down", &MET_T1Smear_phi_jesRelativePtEC2Down, &b_MET_T1Smear_phi_jesRelativePtEC2Down);
   fChain->SetBranchAddress("Jet_pt_jesRelativePtHFDown", Jet_pt_jesRelativePtHFDown, &b_Jet_pt_jesRelativePtHFDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativePtHFDown", Jet_mass_jesRelativePtHFDown, &b_Jet_mass_jesRelativePtHFDown);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativePtHFDown", &MET_T1_pt_jesRelativePtHFDown, &b_MET_T1_pt_jesRelativePtHFDown);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativePtHFDown", &MET_T1_phi_jesRelativePtHFDown, &b_MET_T1_phi_jesRelativePtHFDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativePtHFDown", &MET_T1Smear_pt_jesRelativePtHFDown, &b_MET_T1Smear_pt_jesRelativePtHFDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativePtHFDown", &MET_T1Smear_phi_jesRelativePtHFDown, &b_MET_T1Smear_phi_jesRelativePtHFDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeBalDown", Jet_pt_jesRelativeBalDown, &b_Jet_pt_jesRelativeBalDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeBalDown", Jet_mass_jesRelativeBalDown, &b_Jet_mass_jesRelativeBalDown);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeBalDown", &MET_T1_pt_jesRelativeBalDown, &b_MET_T1_pt_jesRelativeBalDown);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeBalDown", &MET_T1_phi_jesRelativeBalDown, &b_MET_T1_phi_jesRelativeBalDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeBalDown", &MET_T1Smear_pt_jesRelativeBalDown, &b_MET_T1Smear_pt_jesRelativeBalDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeBalDown", &MET_T1Smear_phi_jesRelativeBalDown, &b_MET_T1Smear_phi_jesRelativeBalDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeSampleDown", Jet_pt_jesRelativeSampleDown, &b_Jet_pt_jesRelativeSampleDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeSampleDown", Jet_mass_jesRelativeSampleDown, &b_Jet_mass_jesRelativeSampleDown);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeSampleDown", &MET_T1_pt_jesRelativeSampleDown, &b_MET_T1_pt_jesRelativeSampleDown);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeSampleDown", &MET_T1_phi_jesRelativeSampleDown, &b_MET_T1_phi_jesRelativeSampleDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeSampleDown", &MET_T1Smear_pt_jesRelativeSampleDown, &b_MET_T1Smear_pt_jesRelativeSampleDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeSampleDown", &MET_T1Smear_phi_jesRelativeSampleDown, &b_MET_T1Smear_phi_jesRelativeSampleDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeFSRDown", Jet_pt_jesRelativeFSRDown, &b_Jet_pt_jesRelativeFSRDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeFSRDown", Jet_mass_jesRelativeFSRDown, &b_Jet_mass_jesRelativeFSRDown);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeFSRDown", &MET_T1_pt_jesRelativeFSRDown, &b_MET_T1_pt_jesRelativeFSRDown);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeFSRDown", &MET_T1_phi_jesRelativeFSRDown, &b_MET_T1_phi_jesRelativeFSRDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeFSRDown", &MET_T1Smear_pt_jesRelativeFSRDown, &b_MET_T1Smear_pt_jesRelativeFSRDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeFSRDown", &MET_T1Smear_phi_jesRelativeFSRDown, &b_MET_T1Smear_phi_jesRelativeFSRDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatFSRDown", Jet_pt_jesRelativeStatFSRDown, &b_Jet_pt_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatFSRDown", Jet_mass_jesRelativeStatFSRDown, &b_Jet_mass_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeStatFSRDown", &MET_T1_pt_jesRelativeStatFSRDown, &b_MET_T1_pt_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeStatFSRDown", &MET_T1_phi_jesRelativeStatFSRDown, &b_MET_T1_phi_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeStatFSRDown", &MET_T1Smear_pt_jesRelativeStatFSRDown, &b_MET_T1Smear_pt_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeStatFSRDown", &MET_T1Smear_phi_jesRelativeStatFSRDown, &b_MET_T1Smear_phi_jesRelativeStatFSRDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatECDown", Jet_pt_jesRelativeStatECDown, &b_Jet_pt_jesRelativeStatECDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatECDown", Jet_mass_jesRelativeStatECDown, &b_Jet_mass_jesRelativeStatECDown);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeStatECDown", &MET_T1_pt_jesRelativeStatECDown, &b_MET_T1_pt_jesRelativeStatECDown);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeStatECDown", &MET_T1_phi_jesRelativeStatECDown, &b_MET_T1_phi_jesRelativeStatECDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeStatECDown", &MET_T1Smear_pt_jesRelativeStatECDown, &b_MET_T1Smear_pt_jesRelativeStatECDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeStatECDown", &MET_T1Smear_phi_jesRelativeStatECDown, &b_MET_T1Smear_phi_jesRelativeStatECDown);
   fChain->SetBranchAddress("Jet_pt_jesRelativeStatHFDown", Jet_pt_jesRelativeStatHFDown, &b_Jet_pt_jesRelativeStatHFDown);
   fChain->SetBranchAddress("Jet_mass_jesRelativeStatHFDown", Jet_mass_jesRelativeStatHFDown, &b_Jet_mass_jesRelativeStatHFDown);
   fChain->SetBranchAddress("MET_T1_pt_jesRelativeStatHFDown", &MET_T1_pt_jesRelativeStatHFDown, &b_MET_T1_pt_jesRelativeStatHFDown);
   fChain->SetBranchAddress("MET_T1_phi_jesRelativeStatHFDown", &MET_T1_phi_jesRelativeStatHFDown, &b_MET_T1_phi_jesRelativeStatHFDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesRelativeStatHFDown", &MET_T1Smear_pt_jesRelativeStatHFDown, &b_MET_T1Smear_pt_jesRelativeStatHFDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesRelativeStatHFDown", &MET_T1Smear_phi_jesRelativeStatHFDown, &b_MET_T1Smear_phi_jesRelativeStatHFDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpDataMCDown", Jet_pt_jesPileUpDataMCDown, &b_Jet_pt_jesPileUpDataMCDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpDataMCDown", Jet_mass_jesPileUpDataMCDown, &b_Jet_mass_jesPileUpDataMCDown);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpDataMCDown", &MET_T1_pt_jesPileUpDataMCDown, &b_MET_T1_pt_jesPileUpDataMCDown);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpDataMCDown", &MET_T1_phi_jesPileUpDataMCDown, &b_MET_T1_phi_jesPileUpDataMCDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpDataMCDown", &MET_T1Smear_pt_jesPileUpDataMCDown, &b_MET_T1Smear_pt_jesPileUpDataMCDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpDataMCDown", &MET_T1Smear_phi_jesPileUpDataMCDown, &b_MET_T1Smear_phi_jesPileUpDataMCDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtRefDown", Jet_pt_jesPileUpPtRefDown, &b_Jet_pt_jesPileUpPtRefDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtRefDown", Jet_mass_jesPileUpPtRefDown, &b_Jet_mass_jesPileUpPtRefDown);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpPtRefDown", &MET_T1_pt_jesPileUpPtRefDown, &b_MET_T1_pt_jesPileUpPtRefDown);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpPtRefDown", &MET_T1_phi_jesPileUpPtRefDown, &b_MET_T1_phi_jesPileUpPtRefDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpPtRefDown", &MET_T1Smear_pt_jesPileUpPtRefDown, &b_MET_T1Smear_pt_jesPileUpPtRefDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpPtRefDown", &MET_T1Smear_phi_jesPileUpPtRefDown, &b_MET_T1Smear_phi_jesPileUpPtRefDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtBBDown", Jet_pt_jesPileUpPtBBDown, &b_Jet_pt_jesPileUpPtBBDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtBBDown", Jet_mass_jesPileUpPtBBDown, &b_Jet_mass_jesPileUpPtBBDown);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpPtBBDown", &MET_T1_pt_jesPileUpPtBBDown, &b_MET_T1_pt_jesPileUpPtBBDown);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpPtBBDown", &MET_T1_phi_jesPileUpPtBBDown, &b_MET_T1_phi_jesPileUpPtBBDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpPtBBDown", &MET_T1Smear_pt_jesPileUpPtBBDown, &b_MET_T1Smear_pt_jesPileUpPtBBDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpPtBBDown", &MET_T1Smear_phi_jesPileUpPtBBDown, &b_MET_T1Smear_phi_jesPileUpPtBBDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtEC1Down", Jet_pt_jesPileUpPtEC1Down, &b_Jet_pt_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtEC1Down", Jet_mass_jesPileUpPtEC1Down, &b_Jet_mass_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpPtEC1Down", &MET_T1_pt_jesPileUpPtEC1Down, &b_MET_T1_pt_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpPtEC1Down", &MET_T1_phi_jesPileUpPtEC1Down, &b_MET_T1_phi_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpPtEC1Down", &MET_T1Smear_pt_jesPileUpPtEC1Down, &b_MET_T1Smear_pt_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpPtEC1Down", &MET_T1Smear_phi_jesPileUpPtEC1Down, &b_MET_T1Smear_phi_jesPileUpPtEC1Down);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtEC2Down", Jet_pt_jesPileUpPtEC2Down, &b_Jet_pt_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtEC2Down", Jet_mass_jesPileUpPtEC2Down, &b_Jet_mass_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpPtEC2Down", &MET_T1_pt_jesPileUpPtEC2Down, &b_MET_T1_pt_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpPtEC2Down", &MET_T1_phi_jesPileUpPtEC2Down, &b_MET_T1_phi_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpPtEC2Down", &MET_T1Smear_pt_jesPileUpPtEC2Down, &b_MET_T1Smear_pt_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpPtEC2Down", &MET_T1Smear_phi_jesPileUpPtEC2Down, &b_MET_T1Smear_phi_jesPileUpPtEC2Down);
   fChain->SetBranchAddress("Jet_pt_jesPileUpPtHFDown", Jet_pt_jesPileUpPtHFDown, &b_Jet_pt_jesPileUpPtHFDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpPtHFDown", Jet_mass_jesPileUpPtHFDown, &b_Jet_mass_jesPileUpPtHFDown);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpPtHFDown", &MET_T1_pt_jesPileUpPtHFDown, &b_MET_T1_pt_jesPileUpPtHFDown);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpPtHFDown", &MET_T1_phi_jesPileUpPtHFDown, &b_MET_T1_phi_jesPileUpPtHFDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpPtHFDown", &MET_T1Smear_pt_jesPileUpPtHFDown, &b_MET_T1Smear_pt_jesPileUpPtHFDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpPtHFDown", &MET_T1Smear_phi_jesPileUpPtHFDown, &b_MET_T1Smear_phi_jesPileUpPtHFDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpMuZeroDown", Jet_pt_jesPileUpMuZeroDown, &b_Jet_pt_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpMuZeroDown", Jet_mass_jesPileUpMuZeroDown, &b_Jet_mass_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpMuZeroDown", &MET_T1_pt_jesPileUpMuZeroDown, &b_MET_T1_pt_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpMuZeroDown", &MET_T1_phi_jesPileUpMuZeroDown, &b_MET_T1_phi_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpMuZeroDown", &MET_T1Smear_pt_jesPileUpMuZeroDown, &b_MET_T1Smear_pt_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpMuZeroDown", &MET_T1Smear_phi_jesPileUpMuZeroDown, &b_MET_T1Smear_phi_jesPileUpMuZeroDown);
   fChain->SetBranchAddress("Jet_pt_jesPileUpEnvelopeDown", Jet_pt_jesPileUpEnvelopeDown, &b_Jet_pt_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("Jet_mass_jesPileUpEnvelopeDown", Jet_mass_jesPileUpEnvelopeDown, &b_Jet_mass_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("MET_T1_pt_jesPileUpEnvelopeDown", &MET_T1_pt_jesPileUpEnvelopeDown, &b_MET_T1_pt_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("MET_T1_phi_jesPileUpEnvelopeDown", &MET_T1_phi_jesPileUpEnvelopeDown, &b_MET_T1_phi_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesPileUpEnvelopeDown", &MET_T1Smear_pt_jesPileUpEnvelopeDown, &b_MET_T1Smear_pt_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesPileUpEnvelopeDown", &MET_T1Smear_phi_jesPileUpEnvelopeDown, &b_MET_T1Smear_phi_jesPileUpEnvelopeDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalPileUpDown", Jet_pt_jesSubTotalPileUpDown, &b_Jet_pt_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalPileUpDown", Jet_mass_jesSubTotalPileUpDown, &b_Jet_mass_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("MET_T1_pt_jesSubTotalPileUpDown", &MET_T1_pt_jesSubTotalPileUpDown, &b_MET_T1_pt_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("MET_T1_phi_jesSubTotalPileUpDown", &MET_T1_phi_jesSubTotalPileUpDown, &b_MET_T1_phi_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSubTotalPileUpDown", &MET_T1Smear_pt_jesSubTotalPileUpDown, &b_MET_T1Smear_pt_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSubTotalPileUpDown", &MET_T1Smear_phi_jesSubTotalPileUpDown, &b_MET_T1Smear_phi_jesSubTotalPileUpDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalRelativeDown", Jet_pt_jesSubTotalRelativeDown, &b_Jet_pt_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalRelativeDown", Jet_mass_jesSubTotalRelativeDown, &b_Jet_mass_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("MET_T1_pt_jesSubTotalRelativeDown", &MET_T1_pt_jesSubTotalRelativeDown, &b_MET_T1_pt_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("MET_T1_phi_jesSubTotalRelativeDown", &MET_T1_phi_jesSubTotalRelativeDown, &b_MET_T1_phi_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSubTotalRelativeDown", &MET_T1Smear_pt_jesSubTotalRelativeDown, &b_MET_T1Smear_pt_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSubTotalRelativeDown", &MET_T1Smear_phi_jesSubTotalRelativeDown, &b_MET_T1Smear_phi_jesSubTotalRelativeDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalPtDown", Jet_pt_jesSubTotalPtDown, &b_Jet_pt_jesSubTotalPtDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalPtDown", Jet_mass_jesSubTotalPtDown, &b_Jet_mass_jesSubTotalPtDown);
   fChain->SetBranchAddress("MET_T1_pt_jesSubTotalPtDown", &MET_T1_pt_jesSubTotalPtDown, &b_MET_T1_pt_jesSubTotalPtDown);
   fChain->SetBranchAddress("MET_T1_phi_jesSubTotalPtDown", &MET_T1_phi_jesSubTotalPtDown, &b_MET_T1_phi_jesSubTotalPtDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSubTotalPtDown", &MET_T1Smear_pt_jesSubTotalPtDown, &b_MET_T1Smear_pt_jesSubTotalPtDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSubTotalPtDown", &MET_T1Smear_phi_jesSubTotalPtDown, &b_MET_T1Smear_phi_jesSubTotalPtDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalScaleDown", Jet_pt_jesSubTotalScaleDown, &b_Jet_pt_jesSubTotalScaleDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalScaleDown", Jet_mass_jesSubTotalScaleDown, &b_Jet_mass_jesSubTotalScaleDown);
   fChain->SetBranchAddress("MET_T1_pt_jesSubTotalScaleDown", &MET_T1_pt_jesSubTotalScaleDown, &b_MET_T1_pt_jesSubTotalScaleDown);
   fChain->SetBranchAddress("MET_T1_phi_jesSubTotalScaleDown", &MET_T1_phi_jesSubTotalScaleDown, &b_MET_T1_phi_jesSubTotalScaleDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSubTotalScaleDown", &MET_T1Smear_pt_jesSubTotalScaleDown, &b_MET_T1Smear_pt_jesSubTotalScaleDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSubTotalScaleDown", &MET_T1Smear_phi_jesSubTotalScaleDown, &b_MET_T1Smear_phi_jesSubTotalScaleDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalAbsoluteDown", Jet_pt_jesSubTotalAbsoluteDown, &b_Jet_pt_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalAbsoluteDown", Jet_mass_jesSubTotalAbsoluteDown, &b_Jet_mass_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("MET_T1_pt_jesSubTotalAbsoluteDown", &MET_T1_pt_jesSubTotalAbsoluteDown, &b_MET_T1_pt_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("MET_T1_phi_jesSubTotalAbsoluteDown", &MET_T1_phi_jesSubTotalAbsoluteDown, &b_MET_T1_phi_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSubTotalAbsoluteDown", &MET_T1Smear_pt_jesSubTotalAbsoluteDown, &b_MET_T1Smear_pt_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSubTotalAbsoluteDown", &MET_T1Smear_phi_jesSubTotalAbsoluteDown, &b_MET_T1Smear_phi_jesSubTotalAbsoluteDown);
   fChain->SetBranchAddress("Jet_pt_jesSubTotalMCDown", Jet_pt_jesSubTotalMCDown, &b_Jet_pt_jesSubTotalMCDown);
   fChain->SetBranchAddress("Jet_mass_jesSubTotalMCDown", Jet_mass_jesSubTotalMCDown, &b_Jet_mass_jesSubTotalMCDown);
   fChain->SetBranchAddress("MET_T1_pt_jesSubTotalMCDown", &MET_T1_pt_jesSubTotalMCDown, &b_MET_T1_pt_jesSubTotalMCDown);
   fChain->SetBranchAddress("MET_T1_phi_jesSubTotalMCDown", &MET_T1_phi_jesSubTotalMCDown, &b_MET_T1_phi_jesSubTotalMCDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesSubTotalMCDown", &MET_T1Smear_pt_jesSubTotalMCDown, &b_MET_T1Smear_pt_jesSubTotalMCDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesSubTotalMCDown", &MET_T1Smear_phi_jesSubTotalMCDown, &b_MET_T1Smear_phi_jesSubTotalMCDown);
   fChain->SetBranchAddress("Jet_pt_jesTotalDown", Jet_pt_jesTotalDown, &b_Jet_pt_jesTotalDown);
   fChain->SetBranchAddress("Jet_mass_jesTotalDown", Jet_mass_jesTotalDown, &b_Jet_mass_jesTotalDown);
   fChain->SetBranchAddress("MET_T1_pt_jesTotalDown", &MET_T1_pt_jesTotalDown, &b_MET_T1_pt_jesTotalDown);
   fChain->SetBranchAddress("MET_T1_phi_jesTotalDown", &MET_T1_phi_jesTotalDown, &b_MET_T1_phi_jesTotalDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTotalDown", &MET_T1Smear_pt_jesTotalDown, &b_MET_T1Smear_pt_jesTotalDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTotalDown", &MET_T1Smear_phi_jesTotalDown, &b_MET_T1Smear_phi_jesTotalDown);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoFlavorDown", Jet_pt_jesTotalNoFlavorDown, &b_Jet_pt_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoFlavorDown", Jet_mass_jesTotalNoFlavorDown, &b_Jet_mass_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("MET_T1_pt_jesTotalNoFlavorDown", &MET_T1_pt_jesTotalNoFlavorDown, &b_MET_T1_pt_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("MET_T1_phi_jesTotalNoFlavorDown", &MET_T1_phi_jesTotalNoFlavorDown, &b_MET_T1_phi_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTotalNoFlavorDown", &MET_T1Smear_pt_jesTotalNoFlavorDown, &b_MET_T1Smear_pt_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTotalNoFlavorDown", &MET_T1Smear_phi_jesTotalNoFlavorDown, &b_MET_T1Smear_phi_jesTotalNoFlavorDown);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoTimeDown", Jet_pt_jesTotalNoTimeDown, &b_Jet_pt_jesTotalNoTimeDown);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoTimeDown", Jet_mass_jesTotalNoTimeDown, &b_Jet_mass_jesTotalNoTimeDown);
   fChain->SetBranchAddress("MET_T1_pt_jesTotalNoTimeDown", &MET_T1_pt_jesTotalNoTimeDown, &b_MET_T1_pt_jesTotalNoTimeDown);
   fChain->SetBranchAddress("MET_T1_phi_jesTotalNoTimeDown", &MET_T1_phi_jesTotalNoTimeDown, &b_MET_T1_phi_jesTotalNoTimeDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTotalNoTimeDown", &MET_T1Smear_pt_jesTotalNoTimeDown, &b_MET_T1Smear_pt_jesTotalNoTimeDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTotalNoTimeDown", &MET_T1Smear_phi_jesTotalNoTimeDown, &b_MET_T1Smear_phi_jesTotalNoTimeDown);
   fChain->SetBranchAddress("Jet_pt_jesTotalNoFlavorNoTimeDown", Jet_pt_jesTotalNoFlavorNoTimeDown, &b_Jet_pt_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("Jet_mass_jesTotalNoFlavorNoTimeDown", Jet_mass_jesTotalNoFlavorNoTimeDown, &b_Jet_mass_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("MET_T1_pt_jesTotalNoFlavorNoTimeDown", &MET_T1_pt_jesTotalNoFlavorNoTimeDown, &b_MET_T1_pt_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("MET_T1_phi_jesTotalNoFlavorNoTimeDown", &MET_T1_phi_jesTotalNoFlavorNoTimeDown, &b_MET_T1_phi_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTotalNoFlavorNoTimeDown", &MET_T1Smear_pt_jesTotalNoFlavorNoTimeDown, &b_MET_T1Smear_pt_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTotalNoFlavorNoTimeDown", &MET_T1Smear_phi_jesTotalNoFlavorNoTimeDown, &b_MET_T1Smear_phi_jesTotalNoFlavorNoTimeDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorZJetDown", Jet_pt_jesFlavorZJetDown, &b_Jet_pt_jesFlavorZJetDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorZJetDown", Jet_mass_jesFlavorZJetDown, &b_Jet_mass_jesFlavorZJetDown);
   fChain->SetBranchAddress("MET_T1_pt_jesFlavorZJetDown", &MET_T1_pt_jesFlavorZJetDown, &b_MET_T1_pt_jesFlavorZJetDown);
   fChain->SetBranchAddress("MET_T1_phi_jesFlavorZJetDown", &MET_T1_phi_jesFlavorZJetDown, &b_MET_T1_phi_jesFlavorZJetDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFlavorZJetDown", &MET_T1Smear_pt_jesFlavorZJetDown, &b_MET_T1Smear_pt_jesFlavorZJetDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFlavorZJetDown", &MET_T1Smear_phi_jesFlavorZJetDown, &b_MET_T1Smear_phi_jesFlavorZJetDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPhotonJetDown", Jet_pt_jesFlavorPhotonJetDown, &b_Jet_pt_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPhotonJetDown", Jet_mass_jesFlavorPhotonJetDown, &b_Jet_mass_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("MET_T1_pt_jesFlavorPhotonJetDown", &MET_T1_pt_jesFlavorPhotonJetDown, &b_MET_T1_pt_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("MET_T1_phi_jesFlavorPhotonJetDown", &MET_T1_phi_jesFlavorPhotonJetDown, &b_MET_T1_phi_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFlavorPhotonJetDown", &MET_T1Smear_pt_jesFlavorPhotonJetDown, &b_MET_T1Smear_pt_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFlavorPhotonJetDown", &MET_T1Smear_phi_jesFlavorPhotonJetDown, &b_MET_T1Smear_phi_jesFlavorPhotonJetDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureGluonDown", Jet_pt_jesFlavorPureGluonDown, &b_Jet_pt_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureGluonDown", Jet_mass_jesFlavorPureGluonDown, &b_Jet_mass_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("MET_T1_pt_jesFlavorPureGluonDown", &MET_T1_pt_jesFlavorPureGluonDown, &b_MET_T1_pt_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("MET_T1_phi_jesFlavorPureGluonDown", &MET_T1_phi_jesFlavorPureGluonDown, &b_MET_T1_phi_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFlavorPureGluonDown", &MET_T1Smear_pt_jesFlavorPureGluonDown, &b_MET_T1Smear_pt_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFlavorPureGluonDown", &MET_T1Smear_phi_jesFlavorPureGluonDown, &b_MET_T1Smear_phi_jesFlavorPureGluonDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureQuarkDown", Jet_pt_jesFlavorPureQuarkDown, &b_Jet_pt_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureQuarkDown", Jet_mass_jesFlavorPureQuarkDown, &b_Jet_mass_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("MET_T1_pt_jesFlavorPureQuarkDown", &MET_T1_pt_jesFlavorPureQuarkDown, &b_MET_T1_pt_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("MET_T1_phi_jesFlavorPureQuarkDown", &MET_T1_phi_jesFlavorPureQuarkDown, &b_MET_T1_phi_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFlavorPureQuarkDown", &MET_T1Smear_pt_jesFlavorPureQuarkDown, &b_MET_T1Smear_pt_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFlavorPureQuarkDown", &MET_T1Smear_phi_jesFlavorPureQuarkDown, &b_MET_T1Smear_phi_jesFlavorPureQuarkDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureCharmDown", Jet_pt_jesFlavorPureCharmDown, &b_Jet_pt_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureCharmDown", Jet_mass_jesFlavorPureCharmDown, &b_Jet_mass_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("MET_T1_pt_jesFlavorPureCharmDown", &MET_T1_pt_jesFlavorPureCharmDown, &b_MET_T1_pt_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("MET_T1_phi_jesFlavorPureCharmDown", &MET_T1_phi_jesFlavorPureCharmDown, &b_MET_T1_phi_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFlavorPureCharmDown", &MET_T1Smear_pt_jesFlavorPureCharmDown, &b_MET_T1Smear_pt_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFlavorPureCharmDown", &MET_T1Smear_phi_jesFlavorPureCharmDown, &b_MET_T1Smear_phi_jesFlavorPureCharmDown);
   fChain->SetBranchAddress("Jet_pt_jesFlavorPureBottomDown", Jet_pt_jesFlavorPureBottomDown, &b_Jet_pt_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("Jet_mass_jesFlavorPureBottomDown", Jet_mass_jesFlavorPureBottomDown, &b_Jet_mass_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("MET_T1_pt_jesFlavorPureBottomDown", &MET_T1_pt_jesFlavorPureBottomDown, &b_MET_T1_pt_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("MET_T1_phi_jesFlavorPureBottomDown", &MET_T1_phi_jesFlavorPureBottomDown, &b_MET_T1_phi_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesFlavorPureBottomDown", &MET_T1Smear_pt_jesFlavorPureBottomDown, &b_MET_T1Smear_pt_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesFlavorPureBottomDown", &MET_T1Smear_phi_jesFlavorPureBottomDown, &b_MET_T1Smear_phi_jesFlavorPureBottomDown);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunBDown", Jet_pt_jesTimeRunBDown, &b_Jet_pt_jesTimeRunBDown);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunBDown", Jet_mass_jesTimeRunBDown, &b_Jet_mass_jesTimeRunBDown);
   fChain->SetBranchAddress("MET_T1_pt_jesTimeRunBDown", &MET_T1_pt_jesTimeRunBDown, &b_MET_T1_pt_jesTimeRunBDown);
   fChain->SetBranchAddress("MET_T1_phi_jesTimeRunBDown", &MET_T1_phi_jesTimeRunBDown, &b_MET_T1_phi_jesTimeRunBDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTimeRunBDown", &MET_T1Smear_pt_jesTimeRunBDown, &b_MET_T1Smear_pt_jesTimeRunBDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTimeRunBDown", &MET_T1Smear_phi_jesTimeRunBDown, &b_MET_T1Smear_phi_jesTimeRunBDown);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunCDown", Jet_pt_jesTimeRunCDown, &b_Jet_pt_jesTimeRunCDown);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunCDown", Jet_mass_jesTimeRunCDown, &b_Jet_mass_jesTimeRunCDown);
   fChain->SetBranchAddress("MET_T1_pt_jesTimeRunCDown", &MET_T1_pt_jesTimeRunCDown, &b_MET_T1_pt_jesTimeRunCDown);
   fChain->SetBranchAddress("MET_T1_phi_jesTimeRunCDown", &MET_T1_phi_jesTimeRunCDown, &b_MET_T1_phi_jesTimeRunCDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTimeRunCDown", &MET_T1Smear_pt_jesTimeRunCDown, &b_MET_T1Smear_pt_jesTimeRunCDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTimeRunCDown", &MET_T1Smear_phi_jesTimeRunCDown, &b_MET_T1Smear_phi_jesTimeRunCDown);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunDEDown", Jet_pt_jesTimeRunDEDown, &b_Jet_pt_jesTimeRunDEDown);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunDEDown", Jet_mass_jesTimeRunDEDown, &b_Jet_mass_jesTimeRunDEDown);
   fChain->SetBranchAddress("MET_T1_pt_jesTimeRunDEDown", &MET_T1_pt_jesTimeRunDEDown, &b_MET_T1_pt_jesTimeRunDEDown);
   fChain->SetBranchAddress("MET_T1_phi_jesTimeRunDEDown", &MET_T1_phi_jesTimeRunDEDown, &b_MET_T1_phi_jesTimeRunDEDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTimeRunDEDown", &MET_T1Smear_pt_jesTimeRunDEDown, &b_MET_T1Smear_pt_jesTimeRunDEDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTimeRunDEDown", &MET_T1Smear_phi_jesTimeRunDEDown, &b_MET_T1Smear_phi_jesTimeRunDEDown);
   fChain->SetBranchAddress("Jet_pt_jesTimeRunFDown", Jet_pt_jesTimeRunFDown, &b_Jet_pt_jesTimeRunFDown);
   fChain->SetBranchAddress("Jet_mass_jesTimeRunFDown", Jet_mass_jesTimeRunFDown, &b_Jet_mass_jesTimeRunFDown);
   fChain->SetBranchAddress("MET_T1_pt_jesTimeRunFDown", &MET_T1_pt_jesTimeRunFDown, &b_MET_T1_pt_jesTimeRunFDown);
   fChain->SetBranchAddress("MET_T1_phi_jesTimeRunFDown", &MET_T1_phi_jesTimeRunFDown, &b_MET_T1_phi_jesTimeRunFDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesTimeRunFDown", &MET_T1Smear_pt_jesTimeRunFDown, &b_MET_T1Smear_pt_jesTimeRunFDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesTimeRunFDown", &MET_T1Smear_phi_jesTimeRunFDown, &b_MET_T1Smear_phi_jesTimeRunFDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupMPFInSituDown", Jet_pt_jesCorrelationGroupMPFInSituDown, &b_Jet_pt_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupMPFInSituDown", Jet_mass_jesCorrelationGroupMPFInSituDown, &b_Jet_mass_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("MET_T1_pt_jesCorrelationGroupMPFInSituDown", &MET_T1_pt_jesCorrelationGroupMPFInSituDown, &b_MET_T1_pt_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("MET_T1_phi_jesCorrelationGroupMPFInSituDown", &MET_T1_phi_jesCorrelationGroupMPFInSituDown, &b_MET_T1_phi_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesCorrelationGroupMPFInSituDown", &MET_T1Smear_pt_jesCorrelationGroupMPFInSituDown, &b_MET_T1Smear_pt_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesCorrelationGroupMPFInSituDown", &MET_T1Smear_phi_jesCorrelationGroupMPFInSituDown, &b_MET_T1Smear_phi_jesCorrelationGroupMPFInSituDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupIntercalibrationDown", Jet_pt_jesCorrelationGroupIntercalibrationDown, &b_Jet_pt_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupIntercalibrationDown", Jet_mass_jesCorrelationGroupIntercalibrationDown, &b_Jet_mass_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("MET_T1_pt_jesCorrelationGroupIntercalibrationDown", &MET_T1_pt_jesCorrelationGroupIntercalibrationDown, &b_MET_T1_pt_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("MET_T1_phi_jesCorrelationGroupIntercalibrationDown", &MET_T1_phi_jesCorrelationGroupIntercalibrationDown, &b_MET_T1_phi_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesCorrelationGroupIntercalibrationDown", &MET_T1Smear_pt_jesCorrelationGroupIntercalibrationDown, &b_MET_T1Smear_pt_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesCorrelationGroupIntercalibrationDown", &MET_T1Smear_phi_jesCorrelationGroupIntercalibrationDown, &b_MET_T1Smear_phi_jesCorrelationGroupIntercalibrationDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupbJESDown", Jet_pt_jesCorrelationGroupbJESDown, &b_Jet_pt_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupbJESDown", Jet_mass_jesCorrelationGroupbJESDown, &b_Jet_mass_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("MET_T1_pt_jesCorrelationGroupbJESDown", &MET_T1_pt_jesCorrelationGroupbJESDown, &b_MET_T1_pt_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("MET_T1_phi_jesCorrelationGroupbJESDown", &MET_T1_phi_jesCorrelationGroupbJESDown, &b_MET_T1_phi_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesCorrelationGroupbJESDown", &MET_T1Smear_pt_jesCorrelationGroupbJESDown, &b_MET_T1Smear_pt_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesCorrelationGroupbJESDown", &MET_T1Smear_phi_jesCorrelationGroupbJESDown, &b_MET_T1Smear_phi_jesCorrelationGroupbJESDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupFlavorDown", Jet_pt_jesCorrelationGroupFlavorDown, &b_Jet_pt_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupFlavorDown", Jet_mass_jesCorrelationGroupFlavorDown, &b_Jet_mass_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("MET_T1_pt_jesCorrelationGroupFlavorDown", &MET_T1_pt_jesCorrelationGroupFlavorDown, &b_MET_T1_pt_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("MET_T1_phi_jesCorrelationGroupFlavorDown", &MET_T1_phi_jesCorrelationGroupFlavorDown, &b_MET_T1_phi_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesCorrelationGroupFlavorDown", &MET_T1Smear_pt_jesCorrelationGroupFlavorDown, &b_MET_T1Smear_pt_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesCorrelationGroupFlavorDown", &MET_T1Smear_phi_jesCorrelationGroupFlavorDown, &b_MET_T1Smear_phi_jesCorrelationGroupFlavorDown);
   fChain->SetBranchAddress("Jet_pt_jesCorrelationGroupUncorrelatedDown", Jet_pt_jesCorrelationGroupUncorrelatedDown, &b_Jet_pt_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("Jet_mass_jesCorrelationGroupUncorrelatedDown", Jet_mass_jesCorrelationGroupUncorrelatedDown, &b_Jet_mass_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("MET_T1_pt_jesCorrelationGroupUncorrelatedDown", &MET_T1_pt_jesCorrelationGroupUncorrelatedDown, &b_MET_T1_pt_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("MET_T1_phi_jesCorrelationGroupUncorrelatedDown", &MET_T1_phi_jesCorrelationGroupUncorrelatedDown, &b_MET_T1_phi_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_jesCorrelationGroupUncorrelatedDown", &MET_T1Smear_pt_jesCorrelationGroupUncorrelatedDown, &b_MET_T1Smear_pt_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_jesCorrelationGroupUncorrelatedDown", &MET_T1Smear_phi_jesCorrelationGroupUncorrelatedDown, &b_MET_T1Smear_phi_jesCorrelationGroupUncorrelatedDown);
   fChain->SetBranchAddress("MET_T1_pt_unclustEnDown", &MET_T1_pt_unclustEnDown, &b_MET_T1_pt_unclustEnDown);
   fChain->SetBranchAddress("MET_T1_phi_unclustEnDown", &MET_T1_phi_unclustEnDown, &b_MET_T1_phi_unclustEnDown);
   fChain->SetBranchAddress("MET_T1Smear_pt_unclustEnDown", &MET_T1Smear_pt_unclustEnDown, &b_MET_T1Smear_pt_unclustEnDown);
   fChain->SetBranchAddress("MET_T1Smear_phi_unclustEnDown", &MET_T1Smear_phi_unclustEnDown, &b_MET_T1Smear_phi_unclustEnDown);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("puWeightUp", &puWeightUp, &b_puWeightUp);
   fChain->SetBranchAddress("puWeightDown", &puWeightDown, &b_puWeightDown);
   fChain->SetBranchAddress("PrefireWeight", &PrefireWeight, &b_PrefireWeight);
   fChain->SetBranchAddress("PrefireWeight_Up", &PrefireWeight_Up, &b_PrefireWeight_Up);
   fChain->SetBranchAddress("PrefireWeight_Down", &PrefireWeight_Down, &b_PrefireWeight_Down);
   fChain->SetBranchAddress("Jet_btagSF_deepcsv_M_down", Jet_btagSF_deepcsv_M_down, &b_Jet_btagSF_deepcsv_M_down);
   fChain->SetBranchAddress("Jet_btagSF_deepcsv_M", Jet_btagSF_deepcsv_M, &b_Jet_btagSF_deepcsv_M);
   fChain->SetBranchAddress("Jet_btagSF_deepcsv_M_up", Jet_btagSF_deepcsv_M_up, &b_Jet_btagSF_deepcsv_M_up);
   fChain->SetBranchAddress("Jet_btagSF_deepcsv_L_down", Jet_btagSF_deepcsv_L_down, &b_Jet_btagSF_deepcsv_L_down);
   fChain->SetBranchAddress("Jet_btagSF_deepcsv_L", Jet_btagSF_deepcsv_L, &b_Jet_btagSF_deepcsv_L);
   fChain->SetBranchAddress("Jet_btagSF_deepcsv_L_up", Jet_btagSF_deepcsv_L_up, &b_Jet_btagSF_deepcsv_L_up);
   fChain->SetBranchAddress("Jet_btagSF_deepcsv_T_down", Jet_btagSF_deepcsv_T_down, &b_Jet_btagSF_deepcsv_T_down);
   fChain->SetBranchAddress("Jet_btagSF_deepcsv_T", Jet_btagSF_deepcsv_T, &b_Jet_btagSF_deepcsv_T);
   fChain->SetBranchAddress("Jet_btagSF_deepcsv_T_up", Jet_btagSF_deepcsv_T_up, &b_Jet_btagSF_deepcsv_T_up);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_M_down", Jet_btagSF_deepjet_M_down, &b_Jet_btagSF_deepjet_M_down);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_M", Jet_btagSF_deepjet_M, &b_Jet_btagSF_deepjet_M);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_M_up", Jet_btagSF_deepjet_M_up, &b_Jet_btagSF_deepjet_M_up);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_L_down", Jet_btagSF_deepjet_L_down, &b_Jet_btagSF_deepjet_L_down);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_L", Jet_btagSF_deepjet_L, &b_Jet_btagSF_deepjet_L);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_L_up", Jet_btagSF_deepjet_L_up, &b_Jet_btagSF_deepjet_L_up);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_T_down", Jet_btagSF_deepjet_T_down, &b_Jet_btagSF_deepjet_T_down);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_T", Jet_btagSF_deepjet_T, &b_Jet_btagSF_deepjet_T);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_T_up", Jet_btagSF_deepjet_T_up, &b_Jet_btagSF_deepjet_T_up);
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
