#include "Tools/Tools/interface/HHJetsInterface.h"

// Constructor
HHJetsInterface::HHJetsInterface (std::string model_0, std::string model_1, int year, bool extEtaAcc):
  HHbtagger_(std::array<std::string, 2> { {model_0, model_1} })
{
  year_ = year;
  if (extEtaAcc) max_bjet_eta = 2.5;
}


// Destructor
HHJetsInterface::~HHJetsInterface() {}



output HHJetsInterface::GetHHJets(
    unsigned long long int event, int pairType,
    fRVec Jet_pt, fRVec Jet_eta, fRVec Jet_phi, fRVec Jet_mass,
    iRVec Jet_puId, fRVec Jet_jetId, fRVec Jet_btag,
    fRVec SubJet_pt, fRVec SubJet_eta, fRVec SubJet_phi, fRVec SubJet_mass,
    fRVec FatJet_msoftdrop, iRVec FatJet_subJetIdx1, iRVec FatJet_subJetIdx2,
    float dau1_pt, float dau1_eta, float dau1_phi, float dau1_mass,
    float dau2_pt, float dau2_eta, float dau2_phi, float dau2_mass,
    float met_pt, float met_phi)
{
  std::vector<float> all_HHbtag_scores;
  for (size_t ijet = 0; ijet < Jet_pt.size(); ijet++) {
    all_HHbtag_scores.push_back(-999.);
  }

  int bjet1_idx = -1;
  int bjet2_idx = -1;
  int vbfjet1_idx = -1;
  int vbfjet2_idx = -1;
  int isBoosted_ = 0;
  std::vector <int> ctjet_indexes, fwjet_indexes;

  if (pairType < 0) {
    return output({all_HHbtag_scores, bjet1_idx, bjet2_idx, vbfjet1_idx, vbfjet2_idx,
      ctjet_indexes, fwjet_indexes, isBoosted_});
  }

  auto dau1_tlv = TLorentzVector();
  auto dau2_tlv = TLorentzVector();
  auto met_tlv = TLorentzVector();
  dau1_tlv.SetPtEtaPhiM(dau1_pt, dau1_eta, dau1_phi, dau1_mass);
  dau2_tlv.SetPtEtaPhiM(dau2_pt, dau2_eta, dau2_phi, dau2_mass);
  met_tlv.SetPxPyPzE(met_pt * cos(met_phi), met_pt * sin(met_phi), 0, met_pt);

  std::vector <jet_idx_btag> jet_indexes;
  std::vector <int> all_jet_indexes;
  for (size_t ijet = 0; ijet < Jet_pt.size(); ijet++) {
    if ((Jet_puId[ijet] < 1 && Jet_pt[ijet] <= 50) || Jet_jetId[ijet] < 2)
      continue;
    auto jet_tlv = TLorentzVector();
    jet_tlv.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet], Jet_mass[ijet]);
    if (jet_tlv.DeltaR(dau1_tlv) < 0.5 || jet_tlv.DeltaR(dau2_tlv) < 0.5)
      continue;
    if (Jet_pt[ijet] > 20 && fabs(Jet_eta[ijet]) < max_bjet_eta)
      jet_indexes.push_back(jet_idx_btag({(int) ijet, Jet_btag[ijet]}));
    if (Jet_pt[ijet] > 20 && fabs(Jet_eta[ijet]) < 4.7)
      all_jet_indexes.push_back(ijet);
  }
  if (jet_indexes.size() >= 2) {
    std::stable_sort(jet_indexes.begin(), jet_indexes.end(), jetSort);

    auto htt_tlv = dau1_tlv + dau2_tlv;
    std::vector <float> HHbtag_jet_pt_, HHbtag_jet_eta_, HHbtag_rel_jet_M_pt_, HHbtag_rel_jet_E_pt_;
    std::vector <float> HHbtag_jet_htt_deta_, HHbtag_jet_htt_dphi_, HHbtag_jet_btagOutput_;

    for (auto &jet : jet_indexes) {
      auto jet_tlv = TLorentzVector();
      jet_tlv.SetPtEtaPhiM(Jet_pt[jet.idx], Jet_eta[jet.idx], Jet_phi[jet.idx], Jet_mass[jet.idx]);
      HHbtag_jet_pt_.push_back(jet_tlv.Pt());
      HHbtag_jet_eta_.push_back(jet_tlv.Eta());
      HHbtag_rel_jet_M_pt_.push_back(jet_tlv.M() / jet_tlv.Pt());
      HHbtag_rel_jet_E_pt_.push_back(jet_tlv.E() / jet_tlv.Pt());
      HHbtag_jet_htt_deta_.push_back(htt_tlv.Eta() - jet_tlv.Eta());
      HHbtag_jet_htt_dphi_.push_back(ROOT::Math::VectorUtil::DeltaPhi(htt_tlv, jet_tlv));
      HHbtag_jet_btagOutput_.push_back(jet.btag);
    }

    auto HHbtag_htt_met_dphi_ = (float) ROOT::Math::VectorUtil::DeltaPhi(htt_tlv, met_tlv);
    auto HHbtag_htt_scalar_pt_ = (float) (dau1_tlv.Pt() + dau2_tlv.Pt());
    auto HHbtag_rel_met_pt_htt_pt_ = (float) met_tlv.Pt() / HHbtag_htt_scalar_pt_;
    auto HHbtag_htt_pt_ = (float) htt_tlv.Pt();
    auto HHbtag_htt_eta_ = (float) htt_tlv.Eta();
    int HHbtag_channel_ = -1;

    if (pairType == 0)
      HHbtag_channel_ = 1;
    else if (pairType == 1)
      HHbtag_channel_ = 0;
    else if (pairType == 2)
      HHbtag_channel_ = 1;

    auto HHbtag_scores = HHbtagger_.GetScore(HHbtag_jet_pt_, HHbtag_jet_eta_,
      HHbtag_rel_jet_M_pt_, HHbtag_rel_jet_E_pt_, HHbtag_jet_htt_deta_,
      HHbtag_jet_btagOutput_, HHbtag_jet_htt_dphi_, year_, HHbtag_channel_,
      HHbtag_htt_pt_, HHbtag_htt_eta_, HHbtag_htt_met_dphi_,
      HHbtag_rel_met_pt_htt_pt_, HHbtag_htt_scalar_pt_, event);

    for (size_t ijet = 0; ijet < jet_indexes.size(); ijet++) {
      all_HHbtag_scores[jet_indexes[ijet].idx] = HHbtag_scores[ijet];
    }
    std::vector <jet_idx_btag> jet_indexes_hhbtag;
    for (size_t ijet = 0; ijet < jet_indexes.size(); ijet++) {
      jet_indexes_hhbtag.push_back(jet_idx_btag({jet_indexes[ijet].idx, HHbtag_scores[ijet]}));
    }
    std::stable_sort(jet_indexes_hhbtag.begin(), jet_indexes_hhbtag.end(), jetSort);
    bjet1_idx = jet_indexes_hhbtag[0].idx;
    bjet2_idx = jet_indexes_hhbtag[1].idx;

    if (Jet_pt[bjet1_idx] < Jet_pt[bjet2_idx]) {
      auto aux = bjet1_idx;
      bjet1_idx = bjet2_idx;
      bjet2_idx = aux;
    }

    if (all_jet_indexes.size() >= 4) { // 2 bjets + 2 vbf jets
      std::vector <jet_pair_mass> vbfjet_indexes;
      for (size_t ijet = 0; ijet < all_jet_indexes.size(); ijet++) {
        auto jet1_index = all_jet_indexes[ijet];
        if (jet1_index == (int) bjet1_idx || jet1_index == (int) bjet2_idx)
          continue;
        for (size_t jjet = ijet + 1; jjet < all_jet_indexes.size(); jjet++) {
          auto jet2_index = all_jet_indexes[jjet];
          if (jet2_index == (int) bjet1_idx || jet2_index == (int) bjet2_idx)
            continue;
          if (Jet_pt[jet1_index] < 30 || Jet_pt[jet2_index] < 30)
            continue;
          auto jet1_tlv = TLorentzVector();
          auto jet2_tlv = TLorentzVector();
          jet1_tlv.SetPtEtaPhiM(Jet_pt[jet1_index], Jet_eta[jet1_index],
            Jet_phi[jet1_index], Jet_mass[jet1_index]);
          jet2_tlv.SetPtEtaPhiM(Jet_pt[jet2_index], Jet_eta[jet2_index],
            Jet_phi[jet2_index], Jet_mass[jet2_index]);
          auto jj_tlv = jet1_tlv + jet2_tlv;
          vbfjet_indexes.push_back(jet_pair_mass({jet1_index, jet2_index, (float) jj_tlv.M()}));
        }
      }
      if (vbfjet_indexes.size() > 0) {
        std::stable_sort(vbfjet_indexes.begin(), vbfjet_indexes.end(), jetPairSort);
        vbfjet1_idx = vbfjet_indexes[0].idx1;
        vbfjet2_idx = vbfjet_indexes[0].idx2;

        if (Jet_pt[vbfjet1_idx] < Jet_pt[vbfjet2_idx]) {
          auto aux = vbfjet1_idx;
          vbfjet1_idx = vbfjet2_idx;
          vbfjet2_idx = aux;
        }
      }
    }

    // additional central and forward jets
    for (auto & ijet : all_jet_indexes) {
      if ((int) ijet == bjet1_idx || (int) ijet == bjet2_idx
          || (int) ijet == vbfjet1_idx || (int) ijet == vbfjet2_idx)
        continue;
      if (fabs(Jet_eta[ijet]) < max_bjet_eta)
        ctjet_indexes.push_back(ijet);
      else if (fabs(Jet_eta[ijet]) < 4.7 && Jet_pt[ijet] > 30)
        fwjet_indexes.push_back(ijet);
    }

    // is the event boosted?
    // we loop over the fat AK8 jets, apply a mass cut and verify that its subjets match
    // the jets we selected before.
    auto bjet1_tlv = TLorentzVector();
    auto bjet2_tlv = TLorentzVector();
    bjet1_tlv.SetPtEtaPhiM(Jet_pt[bjet1_idx], Jet_eta[bjet1_idx],
      Jet_phi[bjet1_idx], Jet_mass[bjet1_idx]);
    bjet2_tlv.SetPtEtaPhiM(Jet_pt[bjet2_idx], Jet_eta[bjet2_idx],
      Jet_phi[bjet2_idx], Jet_mass[bjet2_idx]);
    for (size_t ifatjet = 0; ifatjet < FatJet_msoftdrop.size(); ifatjet++) {
      if (FatJet_msoftdrop[ifatjet] < 30)
        continue;
      if (FatJet_subJetIdx1[ifatjet] == -1 || FatJet_subJetIdx2[ifatjet] == -1)
        continue;
      auto subidx1 = FatJet_subJetIdx1[ifatjet];
      auto subidx2 = FatJet_subJetIdx2[ifatjet];
      auto subj1_tlv = TLorentzVector();
      auto subj2_tlv = TLorentzVector();
      subj1_tlv.SetPtEtaPhiM(SubJet_pt[subidx1], SubJet_eta[subidx1],
        SubJet_phi[subidx1], SubJet_mass[subidx1]);
      subj2_tlv.SetPtEtaPhiM(SubJet_pt[subidx2], SubJet_eta[subidx2],
        SubJet_phi[subidx2], SubJet_mass[subidx2]);
      if ((fabs(bjet1_tlv.DeltaR(subj1_tlv)) > 0.4 || fabs(bjet2_tlv.DeltaR(subj2_tlv)) > 0.4) &&
          (fabs(bjet1_tlv.DeltaR(subj2_tlv)) > 0.4 || fabs(bjet2_tlv.DeltaR(subj1_tlv)) > 0.4))
        continue;
      // setBoosted(1);
      isBoosted_ = 1;
    }

  }
  return output({all_HHbtag_scores, bjet1_idx, bjet2_idx, vbfjet1_idx, vbfjet2_idx,
    ctjet_indexes, fwjet_indexes, isBoosted_});
}


// GetScore
std::vector<float> HHJetsInterface::GetScore(
  std::vector<float> HHbtag_jet_pt_, std::vector<float> HHbtag_jet_eta_, std::vector<float> HHbtag_rel_jet_M_pt_,
  std::vector<float> HHbtag_rel_jet_E_pt_, std::vector<float> HHbtag_jet_htt_deta_, std::vector<float> HHbtag_jet_btagOutput_,
  std::vector<float> HHbtag_jet_htt_dphi_, int HHbtag_year_, int HHbtag_channel_, float HHbtag_tauH_pt_, float HHbtag_tauH_eta_,
  float HHbtag_htt_met_dphi_, float HHbtag_rel_met_pt_htt_pt_, float HHbtag_htt_scalar_pt_, unsigned long long int HHbtag_evt_)
{
  // Get HHbtag score
  auto HHbtag_scores = HHbtagger_.GetScore(HHbtag_jet_pt_, HHbtag_jet_eta_, HHbtag_rel_jet_M_pt_,
      HHbtag_rel_jet_E_pt_, HHbtag_jet_htt_deta_, HHbtag_jet_btagOutput_, HHbtag_jet_htt_dphi_,
      HHbtag_year_, HHbtag_channel_, HHbtag_tauH_pt_, HHbtag_tauH_eta_, HHbtag_htt_met_dphi_,
      HHbtag_rel_met_pt_htt_pt_, HHbtag_htt_scalar_pt_, HHbtag_evt_);

  // Store HHbtag scores in a map<jet_idx,HHbtag_score>
  std::vector<float> jets_and_HHbtag;
  for (unsigned int i = 0; i < HHbtag_jet_pt_.size(); i++)
  {
    jets_and_HHbtag.push_back(HHbtag_scores.at(i));
  }

  return jets_and_HHbtag;
}
