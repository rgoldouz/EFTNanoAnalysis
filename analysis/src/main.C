#include "MyAnalysis.h"
int main(){
    TChain* ch    = new TChain("IIHEAnalysis") ;
    ch ->Add("/hadoop/store/user/rgoldouz/FullProduction/BNV/ntuple_BNV_ST_TDUE/*");
//    ch ->Add("/hadoop/store/user/rgoldouz/FullProduction/FCNC/ntuple_tuFCNC_ullDecay_noH/outfile_328.root");    
    TChain ch1("meta") ;
//    ch1.Add("/hadoop/store/user/rgoldouz/FullProduction/FCNC/ntuple_tuFCNC_ullDecay_noH/outfile_328.root");
    ch1.Add("/hadoop/store/user/rgoldouz/FullProduction/BNV/ntuple_BNV_ST_TDUE/outfile_198.root");    
    ch ->AddFriend("meta");
    MyAnalysis t1(ch);
//void MyAnalysis::Loop(TString fname, TString data, TString dataset ,TString year, TString run, float xs, float lumi, float Nevent, int iseft, int nRuns)
    t1.Loop("test.root", "mc" , "" , "2018" , "" , 2.06 ,1 , 3000,1,40);
}
