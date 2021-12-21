#include "jet_candidate.h"

jet_candidate::jet_candidate(float pt_in, float eta_in, float phi_in, float E_in, float btag_in, TString year, int flavor_in, int ind_in, int NtopObj_in){
  pt_ = pt_in;
  eta_ = eta_in;
  phi_ = phi_in;
  btag_ = isb(btag_in,year);
  p4_.SetPtEtaPhiE(pt_, eta_, phi_, E_in) ;
  flavor_ = flavor_in;
  indice_ = ind_in;
  NtopObj_ = NtopObj_in;
}


int jet_candidate::isb(float btag_in , TString year){
  int R = 0;
  if (year == "2016" && btag_in > 0.6321) R=1;
  if (year == "2017" && btag_in > 0.4941) R=1;
  if (year == "2018" && btag_in > 0.4184) R=1;
  return R;
}
  
jet_candidate::~jet_candidate(){}


