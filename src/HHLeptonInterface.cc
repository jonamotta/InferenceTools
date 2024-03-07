#include "Tools/Tools/interface/HHLeptonInterface.h"

// Constructors

HHLeptonInterface::HHLeptonInterface (
    int vvvl_vsjet, int vl_vse, int vvl_vse, int t_vsmu, int vl_vsmu) {
  vvvl_vsjet_ = vvvl_vsjet;
  vl_vse_ = vl_vse;
  vvl_vse_ = vvl_vse;
  t_vsmu_ = t_vsmu;
  vl_vsmu_ = vl_vsmu;
};

// Destructor
HHLeptonInterface::~HHLeptonInterface() {}

lepton_output HHLeptonInterface::get_dau_indexes(
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
  )
{
  // mutau channel
  std::vector<int> goodmuons;
  for (size_t imuon = 0; imuon < Muon_pt.size(); imuon ++) {
    if (fabs(Muon_eta[imuon]) > 2.1 || Muon_pfRelIso04_all[imuon] > 0.15
        || fabs(Muon_dxy[imuon]) > 0.045 || fabs(Muon_dz[imuon]) > 0.2
        || !Muon_tightId[imuon])
      continue;
    goodmuons.push_back(imuon);
  }  // loop over muons
  if (goodmuons.size() >= 1) {
    std::vector<int> goodtaus;
    for (size_t itau = 0; itau < Tau_pt.size(); itau ++) {
      if (Tau_tauIdVSmu[itau] < t_vsmu_ || Tau_tauIdVSe[itau] < vl_vse_
          || Tau_tauIdVSjet[itau] < vvvl_vsjet_)
        continue;
      if (fabs(Tau_dz[itau]) > 0.2)
        continue;
      if (Tau_decayMode[itau] != 0 && Tau_decayMode[itau] != 1
          && Tau_decayMode[itau] != 10 && Tau_decayMode[itau] != 11)
        continue;
      // common tau pt req for both single and cross triggers
      if (Tau_pt[itau] <= 20.)
        continue;
      goodtaus.push_back(itau);
    } // loop over taus

    std::vector<tau_pair> tau_pairs;

    for (auto & imuon: goodmuons) {
      auto muon_tlv = TLorentzVector();
      muon_tlv.SetPtEtaPhiM(Muon_pt[imuon], Muon_eta[imuon], Muon_phi[imuon], Muon_mass[imuon]);
      for (auto & itau: goodtaus) {
        auto tau_tlv = TLorentzVector();
        tau_tlv.SetPtEtaPhiM(Tau_pt[itau], Tau_eta[itau], Tau_phi[itau], Tau_mass[itau]);
        if (tau_tlv.DeltaR(muon_tlv) < 0.5)
          continue;

        if (!pass_trigger(
            muon_tlv.Pt(), muon_tlv.Eta(), muon_tlv.Phi(), 13,
            tau_tlv.Pt(), tau_tlv.Eta(), tau_tlv.Phi(), 15,
            mutau_triggers, TrigObj_id, TrigObj_filterBits, TrigObj_eta, TrigObj_phi))
          continue;

        tau_pairs.push_back(tau_pair({imuon, Muon_pfRelIso04_all[imuon], Muon_pt[imuon],
          itau, Tau_rawTauIdVSjet[itau], Tau_pt[itau], 0, 0}));

      } // loop over goodtaus
    } // loop over goodmuons
    if (tau_pairs.size() > 0) {
      std::stable_sort(tau_pairs.begin(), tau_pairs.end(), pairSort);

      if (lepton_veto(tau_pairs[0].index1, -1,
          Muon_pt, Muon_eta, Muon_dz, Muon_dxy,
          Muon_pfRelIso04_all, Muon_mediumId, Muon_tightId,
          Electron_pt, Electron_eta, Electron_dz, Electron_dxy,
          Electron_mvaFall17V2noIso_WP90, Electron_mvaFall17V2Iso_WP90,
          Electron_pfRelIso03_all))
        return lepton_output({
          -1, -1, -1, -1, -1, -1,
          -1., -1., -1., -1, -1, -1, -1,
          -1., -1., -1, -1, -1, -1});
      int ind1 = tau_pairs[0].index1;
      int ind2 = tau_pairs[0].index2;
      int isOS = (int) (Muon_charge[ind1] != Tau_charge[ind2]);

      return lepton_output({0, ind1, ind2, 0, 0, isOS,
        Muon_eta[ind1], Muon_phi[ind1], Muon_pfRelIso04_all[ind1], -1, -1, -1, -1,
        Tau_eta[ind2], Tau_phi[ind2], Tau_decayMode[ind2],
        Tau_tauIdVSe[ind2], Tau_tauIdVSmu[ind2],
        Tau_tauIdVSjet[ind2]});
    }
  }  // goodmuons stuff


  // etau channel
  std::vector<int> goodelectrons;
  for (size_t iele = 0; iele < Electron_pt.size(); iele ++) {
    if (!Electron_mvaFall17V2Iso_WP80[iele]
        || fabs(Electron_dxy[iele]) > 0.045 || fabs(Electron_dz[iele]) > 0.2)
      continue;
    goodelectrons.push_back(iele);
  }  // loop over electrons
  if (goodelectrons.size() >= 1) {
    std::vector<int> goodtaus;
    for (size_t itau = 0; itau < Tau_pt.size(); itau ++) {
      if (Tau_tauIdVSmu[itau] < t_vsmu_ || Tau_tauIdVSe[itau] < vl_vse_
          || Tau_tauIdVSjet[itau] < vvvl_vsjet_)
        continue;
      if (fabs(Tau_dz[itau]) > 0.2)
        continue;
      if (Tau_decayMode[itau] != 0 && Tau_decayMode[itau] != 1
          && Tau_decayMode[itau] != 10 && Tau_decayMode[itau] != 11)
        continue;
      // common tau pt req for both single and cross triggers
      if (Tau_pt[itau] <= 20.)
      continue;
      goodtaus.push_back(itau);
    } // loop over taus

    std::vector<tau_pair> tau_pairs;

    for (auto & iele: goodelectrons) {
      auto electron_tlv = TLorentzVector();
      electron_tlv.SetPtEtaPhiM(Electron_pt[iele], Electron_eta[iele],
        Electron_phi[iele], Electron_mass[iele]);
      for (auto & itau: goodtaus) {
        auto tau_tlv = TLorentzVector();
        tau_tlv.SetPtEtaPhiM(Tau_pt[itau], Tau_eta[itau], Tau_phi[itau], Tau_mass[itau]);
        if (tau_tlv.DeltaR(electron_tlv) < 0.5)
          continue;
        if (!pass_trigger(
            electron_tlv.Pt(), electron_tlv.Eta(), electron_tlv.Phi(), 11,
            tau_tlv.Pt(), tau_tlv.Eta(), tau_tlv.Phi(), 15,
            etau_triggers, TrigObj_id, TrigObj_filterBits, TrigObj_eta, TrigObj_phi)) {
          continue;
        }
        tau_pairs.push_back(tau_pair({iele, Electron_pfRelIso03_all[iele], Electron_pt[iele],
          itau, Tau_rawTauIdVSjet[itau], Tau_pt[itau], 0, 0}));

      } // loop over goodtaus
    } // loop over goodelectrons
    if (tau_pairs.size() > 0) {
      std::stable_sort(tau_pairs.begin(), tau_pairs.end(), pairSort);

      if (lepton_veto(-1, tau_pairs[0].index1,
          Muon_pt, Muon_eta, Muon_dz, Muon_dxy,
          Muon_pfRelIso04_all, Muon_mediumId, Muon_tightId,
          Electron_pt, Electron_eta, Electron_dz, Electron_dxy,
          Electron_mvaFall17V2noIso_WP90, Electron_mvaFall17V2Iso_WP90,
          Electron_pfRelIso03_all))
        return lepton_output({
          -1, -1, -1, -1, -1, -1,
          -1., -1., -1., -1, -1, -1, -1,
          -1., -1., -1, -1, -1, -1});

      int ind1 = tau_pairs[0].index1;
      int ind2 = tau_pairs[0].index2;
      int isOS = (int) (Electron_charge[ind1] != Tau_charge[ind2]);

      return lepton_output({1, ind1, ind2, 0, 0, isOS,
        Electron_eta[ind1], Electron_phi[ind1], Electron_pfRelIso03_all[ind1], -1,
        -1, -1, -1,
        Tau_eta[ind2], Tau_phi[ind2], Tau_decayMode[ind2],
        Tau_tauIdVSe[ind2], Tau_tauIdVSmu[ind2],
        Tau_tauIdVSjet[ind2]});
    }
  }  // goodmuons stuff

  std::vector<int> goodtaus;
  for (size_t itau = 0; itau < Tau_pt.size(); itau ++) {
    if (Tau_tauIdVSmu[itau] < vl_vsmu_ || Tau_tauIdVSe[itau] < vvl_vse_
        || Tau_tauIdVSjet[itau] < vvvl_vsjet_)
      continue;
    if (fabs(Tau_dz[itau]) > 0.2)
      continue;
    if (Tau_decayMode[itau] != 0 && Tau_decayMode[itau] != 1
        && Tau_decayMode[itau] != 10 && Tau_decayMode[itau] != 11)
      continue;
    // common tau pt req for both single and cross triggers
    if (Tau_pt[itau] <= 20.)
      continue;
    goodtaus.push_back(itau);
  } // loop over taus

  if (goodtaus.size() >= 2) {
    std::vector<tau_pair> tau_pairs;
    for (auto & itau1 : goodtaus) {
      auto tau1_tlv = TLorentzVector();
      tau1_tlv.SetPtEtaPhiM(Tau_pt[itau1], Tau_eta[itau1], Tau_phi[itau1], Tau_mass[itau1]);
      for (auto & itau2 : goodtaus) {
        if (itau1 == itau2)
          continue;
        auto tau2_tlv = TLorentzVector();
        tau2_tlv.SetPtEtaPhiM(Tau_pt[itau2], Tau_eta[itau2], Tau_phi[itau2], Tau_mass[itau2]);
        if (tau1_tlv.DeltaR(tau2_tlv) < 0.5)
          continue;

        int pass_tautaujet = 0;
        int pass_vbf = 0;
        if (!pass_trigger(
            tau1_tlv.Pt(), tau1_tlv.Eta(), tau1_tlv.Phi(), 15,
            tau2_tlv.Pt(), tau2_tlv.Eta(), tau2_tlv.Phi(), 15,
            tautau_triggers, TrigObj_id, TrigObj_filterBits, TrigObj_eta, TrigObj_phi)) {
          if (pass_trigger(
              tau1_tlv.Pt(), tau1_tlv.Eta(), tau1_tlv.Phi(), 15,
              tau2_tlv.Pt(), tau2_tlv.Eta(), tau2_tlv.Phi(), 15,
              tautaujet_triggers, TrigObj_id, TrigObj_filterBits, TrigObj_eta, TrigObj_phi))
            pass_tautaujet = 1;
          else if (pass_trigger(
              tau1_tlv.Pt(), tau1_tlv.Eta(), tau1_tlv.Phi(), 15,
              tau2_tlv.Pt(), tau2_tlv.Eta(), tau2_tlv.Phi(), 15,
              vbf_triggers, TrigObj_id, TrigObj_filterBits, TrigObj_eta, TrigObj_phi))
            pass_vbf = 1;
          else continue;
        }
        tau_pairs.push_back(tau_pair({itau1, Tau_rawTauIdVSjet[itau1], Tau_pt[itau1],
          itau2, Tau_rawTauIdVSjet[itau2], Tau_pt[itau2], pass_tautaujet, pass_vbf}));
      }
    }
    if (tau_pairs.size() > 0) {
      std::stable_sort(tau_pairs.begin(), tau_pairs.end(), pairSort);

      if (lepton_veto(-1, -1,
          Muon_pt, Muon_eta, Muon_dz, Muon_dxy,
          Muon_pfRelIso04_all, Muon_mediumId, Muon_tightId,
          Electron_pt, Electron_eta, Electron_dz, Electron_dxy,
          Electron_mvaFall17V2noIso_WP90, Electron_mvaFall17V2Iso_WP90,
          Electron_pfRelIso03_all))
        return lepton_output({
          -1, -1, -1, -1, -1, -1,
          -1., -1., -1., -1, -1, -1, -1,
          -1., -1., -1, -1, -1, -1});

      int ind1 = tau_pairs[0].index1;
      int ind2 = tau_pairs[0].index2;
      int isOS = (int) (Tau_charge[ind1] != Tau_charge[ind2]);

      return lepton_output({2, ind1, ind2,
        tau_pairs[0].isTauTauJetTrigger, tau_pairs[0].isVBFtrigger, isOS,
        Tau_eta[ind1], Tau_phi[ind1], -1., Tau_decayMode[ind1],
        Tau_tauIdVSe[ind1], Tau_tauIdVSmu[ind1],
        Tau_tauIdVSjet[ind1],
        Tau_eta[ind2], Tau_phi[ind2], Tau_decayMode[ind2],
        Tau_tauIdVSe[ind2], Tau_tauIdVSmu[ind2],
        Tau_tauIdVSjet[ind2]});
    }
  }

  return lepton_output({
    -1, -1, -1, -1, -1, -1,
    -1., -1., -1., -1, -1, -1, -1,
    -1., -1., -1, -1, -1, -1});
}


bool HHLeptonInterface::lepton_veto(int muon_index, int electron_index,
    fRVec Muon_pt, fRVec Muon_eta, fRVec Muon_dz, fRVec Muon_dxy,
    fRVec Muon_pfRelIso04_all, bRVec Muon_mediumId, bRVec Muon_tightId,
    fRVec Electron_pt, fRVec Electron_eta, fRVec Electron_dz, fRVec Electron_dxy,
    bRVec Electron_mvaFall17V2noIso_WP90, bRVec Electron_mvaFall17V2Iso_WP90,
    fRVec Electron_pfRelIso03_all)
{
  size_t nleps = 0;
  for (size_t imuon = 0; imuon < Muon_pt.size(); imuon++) {
    if ((int) imuon == muon_index)  // same muon than the chosen in the pair
      continue;
    if (fabs(Muon_eta[imuon]) > 2.4 || Muon_pt[imuon] < 10 || fabs(Muon_dz[imuon]) > 0.2
        || fabs(Muon_dxy[imuon]) > 0.045 || Muon_pfRelIso04_all[imuon] > 0.3)
      continue;
    if (!Muon_mediumId[imuon] and !Muon_tightId[imuon])
      continue;
    nleps++;
  }
  for (size_t iele = 0; iele < Electron_pt.size(); iele++) {
    if ((int) iele == electron_index)  // same electron than the chosen in the pair
      continue;
    if (fabs(Electron_eta[iele]) > 2.5 || Electron_pt[iele] < 10 || fabs(Electron_dz[iele]) > 0.2
        || fabs(Electron_dxy[iele]) > 0.045)
      continue;
    if (!((Electron_pfRelIso03_all[iele] < 0.3 && Electron_mvaFall17V2noIso_WP90[iele])
        || Electron_mvaFall17V2Iso_WP90[iele]))
      continue;
    nleps++;
  }
  return (nleps > 0);
}

bool HHLeptonInterface::pass_trigger(
    float off_pt1, float off_eta1, float off_phi1, int obj_id1,
    float off_pt2, float off_eta2, float off_phi2, int obj_id2,
    std::vector<trig_req> triggers, 
    iRVec TrigObj_id, iRVec TrigObj_filterBits, fRVec TrigObj_eta, fRVec TrigObj_phi)
{
  for (auto &trigger: triggers) {
    if (!trigger.pass)
      continue;
    if (off_pt1 < trigger.pt1 || off_pt2 < trigger.pt2
        || abs(off_eta1) > trigger.eta1 || abs(off_eta2) > trigger.eta2)
      continue;
    if (trigger.bits[0].size() > 0) {
      if (!match_trigger_object(off_eta1, off_phi1, obj_id1,
          TrigObj_id, TrigObj_filterBits, TrigObj_eta, TrigObj_phi, trigger.bits[0]))
        continue;
    }
    if (trigger.bits[1].size() > 0) {
      if (!match_trigger_object(off_eta2, off_phi2, obj_id2,
          TrigObj_id, TrigObj_filterBits, TrigObj_eta, TrigObj_phi, trigger.bits[1]))
        continue;
    }
    return true;
  }
  return false;
}

bool HHLeptonInterface::match_trigger_object(float off_eta, float off_phi, int obj_id,
    iRVec TrigObj_id, iRVec TrigObj_filterBits, fRVec TrigObj_eta, fRVec TrigObj_phi,
    std::vector<int> bits)
{
  for (size_t iobj = 0; iobj < TrigObj_id.size(); iobj++) {
    if (TrigObj_id[iobj] != obj_id) continue;
    auto const dPhi(std::abs(reco::deltaPhi(off_phi, TrigObj_phi[iobj])));
    auto const dEta(std::abs(off_eta - TrigObj_eta[iobj]));
    auto const delR2(dPhi * dPhi + dEta * dEta);
    if (delR2 > 0.5 * 0.5)
      continue;
    bool matched_bits = true;
    for (auto & bit : bits) {
      if ((TrigObj_filterBits[iobj] & bit) == 0) {
        matched_bits = false;
        break;
      }
    }
    if (!matched_bits)
      continue;
    return true;
  }
  return false;
}
