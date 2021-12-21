import sys 
import os 
import subprocess 
import readline 
import string 

UL17={
"UL17_BNV_ST_TDUE":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_ST_TDUE'], 'mc', 'none', '2017', 'none', '1', '41.53', '300000', '1', '300'],
"UL17_ZZTo2L2Nu":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_ZZTo2L2Nu/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/crab_UL17_ZZTo2L2Nu/211106_220010/0000'], 'mc', 'none', '2017', 'none', '0.564', '41.53', '39846157.8106', '0', '1'],
"UL17_DY50":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_DY50/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_UL17_DY50/211203_184109/0000'], 'mc', 'none', '2017', 'none', '6077.22', '41.53', '131552424.895', '0', '1'],
"UL17_WWZ_4F":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_WWZ_4F/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/crab_UL17_WWZ_4F/211106_220254/0000'], 'mc', 'none', '2017', 'none', '0.1651', '41.53', '143284.002169', '0', '1'],
"UL17_ZZTo4L":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_ZZTo4L/ZZTo4L_13TeV_powheg_pythia8/crab_UL17_ZZTo4L/211106_220413/0000'], 'mc', 'none', '2017', 'none', '1.256 ', '41.53', '3791240.98217', '0', '1'],
"UL17_TTZToLLNuNu_M_10":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_TTZToLLNuNu_M_10/TTZToLLNuNu_M-10_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/crab_UL17_TTZToLLNuNu_M_10/211106_220534/0000'], 'mc', 'none', '2017', 'none', '0.2529', '41.53', '5148568.50021', '0', '1'],
"UL17_TTTo2L2Nu":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_TTTo2L2Nu/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_UL17_TTTo2L2Nu/211106_220814/0000'], 'mc', 'none', '2017', 'none', '87.31', '41.53', '96019993.7435', '0', '1'],
"UL17_DY10to50":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_DY10to50/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/crab_UL17_DY10to50/211106_220933/0000'], 'mc', 'none', '2017', 'none', '18610', '41.53', '65673537.0', '0', '1'],
"UL17_BNV_TT_TBUE":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_TT_TBUE'], 'mc', 'none', '2017', 'none', '1', '41.53', '300000', '1', '300'],
"UL17_ZZZ":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_ZZZ/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/crab_UL17_ZZZ/211106_221340/0000'], 'mc', 'none', '2017', 'none', '0.01398', '41.53', '150627.999747', '0', '1'],
"UL17_TTW":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_TTW/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/crab_UL17_TTWJetsToLNu/211106_220653/0000'], 'mc', 'none', '2017', 'none', '0.2043', '41.53', '2683394.35363', '0', '1'],
"UL17_WZTo2L2Q":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_WZTo2L2Q/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_UL17_WZTo2L2Q/211106_221459/0000'], 'mc', 'none', '2017', 'none', '5.595', '41.53', '15660416.1584', '0', '1'],
"UL17_BNV_TT_TSUE":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_TT_TSUE'], 'mc', 'none', '2017', 'none', '1', '41.53', '300000', '1', '300'],
"UL17_BNV_ST_TBUE":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_ST_TBUE'], 'mc', 'none', '2017', 'none', '1', '41.53', '300000', '1', '300'],
"UL17_tW":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_tW/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/crab_UL17_tW/211106_221743/0000'], 'mc', 'none', '2017', 'none', '35.85', '41.53', '5668711.86458', '0', '1'],
"UL17_WJetsToLNu":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_WJetsToLNu/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/crab_UL17_WJetsToLNu/211106_222020/0000'], 'mc', 'none', '2017', 'none', '61526.7', '41.53', '94710339.5406', '0', '1'],
"UL17_BNV_ST_TDCE":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_ST_TDCE'], 'mc', 'none', '2017', 'none', '1', '41.53', '300000', '1', '300'],
"UL17_tbarW":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_tbarW/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/crab_UL17_tbarW/211106_215449/0000'], 'mc', 'none', '2017', 'none', '35.85', '41.53', '5517710.19165', '0', '1'],
"UL17_WWTo2L2Nu":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_WWTo2L2Nu/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_UL17_WWTo2L2Nu/211106_215608/0000'], 'mc', 'none', '2017', 'none', '12.178', '41.53', '6808356.09056', '0', '1'],
"UL17_WWW_4F":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_WWW_4F/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/crab_UL17_WWW_4F/211106_215733/0000'], 'mc', 'none', '2017', 'none', '0.2086', '41.53', '145774.003086', '0', '1'],
"UL17_tuFCNC_tHProduction":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_tuFCNC_tHProduction'], 'mc', 'none', '2017', 'none', '1', '41.53', '997500', '1', '1995'],
"UL17_BNV_ST_TSUE":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_ST_TSUE'], 'mc', 'none', '2017', 'none', '1', '41.53', '300000', '1', '300'],
"UL17_WZTo3LNu":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_WZTo3LNu/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_UL17_WZTo3LNu/211106_215852/0000'], 'mc', 'none', '2017', 'none', '4.43', '41.53', '6826898.12612', '0', '1'],
"UL17_tuFCNC_tllProduction":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_tuFCNC_tllProduction'], 'mc', 'none', '2017', 'none', '1', '41.53', '4986500', '1', '9973'],
"UL17_tuFCNC_ullDecay":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_tuFCNC_ullDecay'], 'mc', 'none', '2017', 'none', '1', '41.53', '1994000', '1', '3988'],
"UL17_BNV_TT_TDCE":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_TT_TDCE'], 'mc', 'none', '2017', 'none', '1', '41.53', '299000', '1', '299'],
"UL17_BNV_ST_TBCE":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_ST_TBCE'], 'mc', 'none', '2017', 'none', '1', '41.53', '300000', '1', '300'],
"UL17_BNV_ST_TSCE":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_ST_TSCE'], 'mc', 'none', '2017', 'none', '1', '41.53', '300000', '1', '300'],
"UL17_BNV_TT_TSCE":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_TT_TSCE'], 'mc', 'none', '2017', 'none', '1', '41.53', '300000', '1', '300'],
"UL17_BNV_TT_TDUE":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_TT_TDUE'], 'mc', 'none', '2017', 'none', '1', '41.53', '300000', '1', '300'],
"UL17_tuFCNC_uHDecay":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_tuFCNC_uHDecay'], 'mc', 'none', '2017', 'none', '1', '41.53', '498000', '1', '996'],
"UL17_BNV_TT_TBCE":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/UL17_BNV_TT_TBCE'], 'mc', 'none', '2017', 'none', '1', '41.53', '300000', '1', '300'],

 
"data_UL17_B_SingleMuon":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_B_SingleMuon/SingleMuon/crab_data_UL17_B_SingleMuon/211109_183753/0000'], 'data', 'SingleMuon', '2017', 'B', '1', '41.53', '1', '0', '1'],
"data_UL17_B_MuonEG":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_B_MuonEG/MuonEG/crab_data_UL17_B_MuonEG/211109_183121/0000'], 'data', 'MuonEG', '2017', 'B', '1', '41.53', '1', '0', '1'],
"data_UL17_C_DoubleEG":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_C_DoubleEG/DoubleEG/crab_data_UL17_C_DoubleEG/211109_183239/0000'], 'data', 'DoubleEG', '2017', 'C', '1', '41.53', '1', '0', '1'],
"data_UL17_F_DoubleEG":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_F_DoubleEG/DoubleEG/crab_data_UL17_F_DoubleEG/211109_184307/0000'], 'data', 'DoubleEG', '2017', 'F', '1', '41.53', '1', '0', '1'],
"data_UL17_C_DoubleMuon":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_C_DoubleMuon/DoubleMuon/crab_data_UL17_C_DoubleMuon/211217_175254/0000'], 'data', 'DoubleMuon', '2017', 'C', '1', '41.53', '1', '0', '1'],
"data_UL17_D_DoubleEG":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_D_DoubleEG/DoubleEG/crab_data_UL17_D_DoubleEG/211109_183003/0000'], 'data', 'DoubleEG', '2017', 'D', '1', '41.53', '1', '0', '1'],
"data_UL17_F_DoubleMuon":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_F_DoubleMuon/DoubleMuon/crab_data_UL17_F_DoubleMuon/211109_184820/0000'], 'data', 'DoubleMuon', '2017', 'F', '1', '41.53', '1', '0', '1'],
"data_UL17_B_SingleElectron":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_B_SingleElectron/SingleElectron/crab_data_UL17_B_SingleElectron/211110_104300/0000'], 'data', 'SingleElectron', '2017', 'B', '1', '41.53', '1', '0', '1'],
"data_UL17_D_SingleMuon":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_D_SingleMuon/SingleMuon/crab_data_UL17_D_SingleMuon/211109_183635/0000'], 'data', 'SingleMuon', '2017', 'D', '1', '41.53', '1', '0', '1'],
"data_UL17_F_MuonEG":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_F_MuonEG/MuonEG/crab_data_UL17_F_MuonEG/211217_175526/0000'], 'data', 'MuonEG', '2017', 'F', '1', '41.53', '1', '0', '1'],
"data_UL17_F_SingleMuon":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_F_SingleMuon/SingleMuon/crab_data_UL17_F_SingleMuon/211217_175410/0000'], 'data', 'SingleMuon', '2017', 'F', '1', '41.53', '1', '0', '1'],
"data_UL17_E_SingleMuon":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_E_SingleMuon/SingleMuon/crab_data_UL17_E_SingleMuon/211109_184030/0000'], 'data', 'SingleMuon', '2017', 'E', '1', '41.53', '1', '0', '1'],
"data_UL17_B_DoubleMuon":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_B_DoubleMuon/DoubleMuon/crab_data_UL17_B_DoubleMuon/211109_184149/0000'], 'data', 'DoubleMuon', '2017', 'B', '1', '41.53', '1', '0', '1'],
"data_UL17_E_SingleElectron":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_E_SingleElectron/SingleElectron/crab_data_UL17_E_SingleElectron/211109_184938/0000'], 'data', 'SingleElectron', '2017', 'E', '1', '41.53', '1', '0', '1'],
"data_UL17_E_MuonEG":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_E_MuonEG/MuonEG/crab_data_UL17_E_MuonEG/211217_175641/0000'], 'data', 'MuonEG', '2017', 'E', '1', '41.53', '1', '0', '1'],
"data_UL17_D_DoubleMuon":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_D_DoubleMuon/DoubleMuon/crab_data_UL17_D_DoubleMuon/211109_184425/0000'], 'data', 'DoubleMuon', '2017', 'D', '1', '41.53', '1', '0', '1'],
"data_UL17_E_DoubleEG":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_E_DoubleEG/DoubleEG/crab_data_UL17_E_DoubleEG/211109_185056/0000'], 'data', 'DoubleEG', '2017', 'E', '1', '41.53', '1', '0', '1'],
"data_UL17_D_SingleElectron":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_D_SingleElectron/SingleElectron/crab_data_UL17_D_SingleElectron/211109_185214/0000'], 'data', 'SingleElectron', '2017', 'D', '1', '41.53', '1', '0', '1'],
"data_UL17_C_MuonEG":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_C_MuonEG/MuonEG/crab_data_UL17_C_MuonEG/211109_183516/0000'], 'data', 'MuonEG', '2017', 'C', '1', '41.53', '1', '0', '1'],
"data_UL17_C_SingleMuon":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_C_SingleMuon/SingleMuon/crab_data_UL17_C_SingleMuon/211217_175756/0000'], 'data', 'SingleMuon', '2017', 'C', '1', '41.53', '1', '0', '1'],
"data_UL17_D_MuonEG":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_D_MuonEG/MuonEG/crab_data_UL17_D_MuonEG/211109_184702/0000'], 'data', 'MuonEG', '2017', 'D', '1', '41.53', '1', '0', '1'],
"data_UL17_B_DoubleEG":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_B_DoubleEG/DoubleEG/crab_data_UL17_B_DoubleEG/211109_185412/0000'], 'data', 'DoubleEG', '2017', 'B', '1', '41.53', '1', '0', '1'],
"data_UL17_E_DoubleMuon":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_E_DoubleMuon/DoubleMuon/crab_data_UL17_E_DoubleMuon/211110_104539/0000'], 'data', 'DoubleMuon', '2017', 'E', '1', '41.53', '1', '0', '1'],
"data_UL17_F_SingleElectron":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_F_SingleElectron/SingleElectron/crab_data_UL17_F_SingleElectron/211217_175912/0000'], 'data', 'SingleElectron', '2017', 'F', '1', '41.53', '1', '0', '1'],
"data_UL17_C_SingleElectron":[['rgoldouz/NanoAodPostProcessingUL/UL17/v1/data_UL17_C_SingleElectron/SingleElectron/crab_data_UL17_C_SingleElectron/211109_184543/0000'], 'data', 'SingleElectron', '2017', 'C', '1', '41.53', '1', '0', '1'],
}