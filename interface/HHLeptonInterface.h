#ifndef HHLeptonInterface_h
#define HHLeptonInterface_h

// -------------------------------------------------------------------------------------------------------------- //
//                                                                                                                //
//   class HHLeptonInterface                                                                                      //
//                                                                                                                //
// -------------------------------------------------------------------------------------------------------------- //

// Standard libraries
#include <vector>
#include <string>
#include <cmath>

// ROOT libraries
#include <TLorentzVector.h>
#include <ROOT/RVec.hxx>

// CMSSW
#include "DataFormats/Math/interface/deltaPhi.h"

typedef ROOT::VecOps::RVec<float> fRVec;
typedef ROOT::VecOps::RVec<bool> bRVec;
typedef ROOT::VecOps::RVec<int> iRVec;

struct tau_pair {
  int index1;
  float iso1;
  float pt1;
  int index2;
  float iso2;
  float pt2;
  int isTauTauJetTrigger;
  int isVBFtrigger;
};

struct trig_req {
  bool pass;
  float pt1;
  float eta1;
  float pt2;
  float eta2;
  std::vector<std::vector<int>> bits;  
};

struct lepton_output {
    int pairType;
    int dau1_index;
    int dau2_index;
    int isTauTauJetTrigger;
    int isVBFtrigger;
    int isOS;

    float dau1_eta;
    float dau1_phi;
    float dau1_iso;
    int dau1_decayMode;
    int dau1_tauIdVSe;
    int dau1_tauIdVSmu;
    int dau1_tauIdVSjet;

    float dau2_eta;
    float dau2_phi;
    int dau2_decayMode;
    int dau2_tauIdVSe;
    int dau2_tauIdVSmu;
    int dau2_tauIdVSjet;
};

bool pairSort (const tau_pair& pA, const tau_pair& pB)
{
  // first leg 1 iso
  float isoA = pA.iso1;
  float isoB = pB.iso1;
  if (isoA > isoB) return true; // NB: MVA iso ! Most iso --> highest MVA score
  else if (isoA < isoB) return false;

  // then leg 1 pt
  float ptA = pA.pt1;
  float ptB = pB.pt1;
  if (ptA > ptB) return true;
  else if (ptA < ptB) return false;

  // then leg 2 iso
  isoA = pA.iso2;
  isoB = pB.iso2;
  if (isoA > isoB) return true;
  else if (isoA < isoB) return false;

  // then leg 2 pt
  ptA = pA.pt2;
  ptB = pB.pt2;
  if (ptA > ptB) return true;
  else if (ptA < ptB) return false;

  // should be never here..
  return false;
}

// HHLeptonInterface class
class HHLeptonInterface {

  public:
    HHLeptonInterface (int vvvl_vsjet, int vl_vse, int vvl_vse, int t_vsmu, int vl_vsmu);
    ~HHLeptonInterface ();
    lepton_output get_dau_indexes(
      fRVec Muon_pt, fRVec Muon_eta, fRVec Muon_phi, fRVec Muon_mass,
      fRVec Muon_pfRelIso04_all, fRVec Muon_dxy, fRVec Muon_dz,
      bRVec Muon_mediumId, bRVec Muon_tightId, iRVec Muon_charge,
      fRVec Electron_pt, fRVec Electron_eta, fRVec Electron_phi, fRVec Electron_mass,
      bRVec Electron_mvaFall17V2Iso_WP80, bRVec Electron_mvaFall17V2noIso_WP90,
      bRVec Electron_mvaFall17V2Iso_WP90, fRVec Electron_pfRelIso03_all,
      fRVec Electron_dxy, fRVec Electron_dz, iRVec Electron_charge,
      fRVec Tau_pt, fRVec Tau_eta, fRVec Tau_phi, fRVec Tau_mass,
      iRVec Tau_tauIdVSmu, iRVec Tau_tauIdVSe,
      iRVec Tau_tauIdVSjet, fRVec Tau_rawTauIdVSjet,
      fRVec Tau_dz, iRVec Tau_decayMode, iRVec Tau_charge,
      iRVec TrigObj_id, iRVec TrigObj_filterBits, fRVec TrigObj_eta, fRVec TrigObj_phi,
      std::vector<trig_req> mutau_triggers, std::vector<trig_req> etau_triggers,
      std::vector<trig_req> tautau_triggers, std::vector<trig_req> tautaujet_triggers, 
      std::vector<trig_req> vbf_triggers
    );

    bool lepton_veto(int muon_index, int electron_index,
      fRVec Muon_pt, fRVec Muon_eta, fRVec Muon_dz, fRVec Muon_dxy,
      fRVec Muon_pfRelIso04_all, bRVec Muon_mediumId, bRVec Muon_tightId,
      fRVec Electron_pt, fRVec Electron_eta, fRVec Electron_dz, fRVec Electron_dxy,
      bRVec Electron_mvaFall17V2noIso_WP90, bRVec Electron_mvaFall17V2Iso_WP90,
      fRVec Electron_pfRelIso03_all);

    bool pass_trigger(
      float off_pt1, float off_eta1, float off_phi1, int obj_id1,
      float off_pt2, float off_eta2, float off_phi2, int obj_id2,
      std::vector<trig_req> triggers, 
      iRVec TrigObj_id, iRVec TrigObj_filterBits, fRVec TrigObj_eta, fRVec TrigObj_phi);
    
    bool match_trigger_object(float off_eta, float off_phi, int obj_id,
      iRVec TrigObj_id, iRVec TrigObj_filterBits, fRVec TrigObj_eta, fRVec TrigObj_phi,
      std::vector<int> bits);

  private:
    int vvvl_vsjet_;
    int vl_vse_;
    int vvl_vse_;
    int t_vsmu_;
    int vl_vsmu_;
};

#endif // HHLeptonInterface_h