#include "../include/MyAnalysis.h"
int main(){
    TChain* ch    = new TChain("Events") ;
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_tW/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/crab_UL17_tW/211106_221743/0000/tree_27.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_E_MuonEG/MuonEG/crab_data_UL17_E_MuonEG/211110_104657/0000/tree_15.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_ST_TDCE/tree_3.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_tuFCNC_tllProduction/tree_81.root");
//    ch ->Add("/afs/crc.nd.edu/user/r/rgoldouz/MakeLobsterJobs/UL/mgprod/lobster_workflow/ul_cfgs/NAOD-00000.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_ST_TDUE/tree_5.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL16postVFP/v1/data_UL16postVFP_H_SingleMuon/SingleMuon/crab_data_UL16postVFP_H_SingleMuon/220101_070435/0000/tree_31.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_D_SingleMuon/SingleMuon/crab_data_UL17_D_SingleMuon/211109_183635/0000/tree_28.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL16preVFP/v1/UL16preVFP_TTTo2L2Nu/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_UL16preVFP_TTTo2L2Nu/220104_170332/0000/tree_6.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_TTTo2L2Nu/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_UL17_TTTo2L2Nu/211106_220814/0000/tree_28.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_DY50/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_UL17_DY50/211203_184109/0000/tree_99.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_ST_TSUE_DAS/*");
    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL16preVFP/v2/UL16preVFP_TTTo2L2Nu/tree_2104.root");    
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v2/UL17_TTBNV_TBUMu/tree_584.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v2/UL17_TTBNV_TBUE/tree_7.root");
//    ch ->Add("root://ndcms.crc.nd.edu//store/mc/RunIISummer20UL17NanoAODv9/ST_TuneCP5_BNV_TDUMu_13TeV-madgraph-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/2540000/556711E3-EA23-F44A-BDC1-F3F3E846EC7D.root");
//    ch ->Add("root://ndcms.crc.nd.edu//store/mc/RunIISummer20UL17NanoAODv9/ST_TuneCP5_BNV_TDUE_13TeV-madgraph-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/2550000/283C5366-FD22-C34C-BBB0-D64B3C8E85E9.root");
//    ch ->Add("/hadoop/store/user/rgoldouz/NanoAodPostProcessingUL/UL17/v2/UL17_STBNV_TDUE/tree_204.root");
//    ch ->Add("/afs/crc.nd.edu/user/r/rgoldouz/Ntupleproducer/MyNano/CMSSW_10_2_24/src/PhysicsTools/NanoAODTools/lobster/Skimmedtree_3.root");
//    ch ->Add("/afs/crc.nd.edu/user/r/rgoldouz/Ntupleproducer/MyNano/CMSSW_10_2_24/src/PhysicsTools/NanoAODTools/lobster/tree.root");
    MyAnalysis * t1 = new MyAnalysis(ch);

//    ch ->Add("/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runD/191203_070404/0000/outfile_85.root");
//    MyAnalysis t1(ch);
//    t1.Loop("2016_D_MuonEG", "2016_D_MuonEG_0_0.root", "data" , "MuonEG" , "2016" , "D" , 1 , 1 , 1);
//    t1->Loop("BNV.root",  "mc" , "" , "2016preVFP" , "" , 87.31 , 35.92 , 67312164,1,102);
t1->Loop("TTTo2L2Nu", "mc","SingleMuon", "2017", "D", 1, 41.48, 494000, 0, 38);
//t1->Loop("BV", "data", "MuonEG", "2016preVFP", "Bv2", 87.31, 19.52, 32368794.5651, 0, 1);

//    ch ->Add("/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_ST_clequ1_emutu/outfile_12096.root");
//    ch ->Add("/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_ST_clequ1_emutu/outfile_12146.root");
//    MyAnalysis t1(ch);
//    t1.Loop("2017_LFVStScalarU_0_0.root","2017_LFVStScalarU", "mc" , "" , "2017" , "" , 0.102 , 41.53 , 500000);
delete t1;
}
