import os
from array import array

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from analysis_tools.utils import import_root
from Tools.Tools.jet_utils import JetPair
from Base.Modules.baseModules import JetLepMetModule, JetLepMetSyst

ROOT = import_root()

class HHJetsProducer(JetLepMetModule):
    def __init__(self, *args, **kwargs):
        extEtaAcc = kwargs.pop("extEtaAcc")
        self.btag_algo = kwargs.pop("btag_algo")
        super(HHJetsProducer, self).__init__(self, *args, **kwargs)

        if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("libToolsTools.so")

        if os.path.expandvars("$CMT_SCRAM_ARCH") == "slc7_amd64_gcc10":
            ROOT.gROOT.ProcessLine(".include /cvmfs/cms.cern.ch/slc7_amd64_gcc10/"
                "external/eigen/d812f411c3f9-cms/include/")
            ROOT.gROOT.ProcessLine(".include /cvmfs/cms.cern.ch/slc7_amd64_gcc10/"
                "external/tensorflow/2.5.0/include/")
        elif os.path.expandvars("$CMT_SCRAM_ARCH") == "slc7_amd64_gcc820":
            ROOT.gROOT.ProcessLine(".include /cvmfs/cms.cern.ch/slc7_amd64_gcc820/"
                "external/eigen/d812f411c3f9-bcolbf/include/eigen3")
            ROOT.gROOT.ProcessLine(".include /cvmfs/cms.cern.ch/slc7_amd64_gcc820/"
                "external/tensorflow/2.1.0-bcolbf/include")
        else:
            raise ValueError("Architecture not considered")

        base = "{}/{}/src/Tools/Tools".format(
            os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))

        ROOT.gROOT.ProcessLine(".L {}/interface/HHJetsInterface.h".format(base))

        self.year = kwargs.pop("year")
        base_hhbtag = "{}/{}/src/HHTools/HHbtag".format(
            os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
        models = [base_hhbtag + "/models/HHbtag_v1_par_%i" % i for i in range(2)]

        self.HHJets = ROOT.HHJetsInterface(models[0], models[1], self.year, extEtaAcc)

        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch('bjet1_JetIdx', 'I')
        self.out.branch('bjet2_JetIdx', 'I')
        self.out.branch('VBFjet1_JetIdx', 'I')
        self.out.branch('VBFjet2_JetIdx', 'I')

        self.out.branch('Jet_HHbtag', "F", lenVar='nJet')

        self.out.branch('isBoosted', 'I')
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        jets = Collection(event, "Jet")
        fatjets = Collection(event, "FatJet")
        subjets = Collection(event, "SubJet")
        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        taus = Collection(event, "Tau")

        dau1, dau2, dau1_tlv, dau2_tlv = self.get_daus(event, muons, electrons, taus)
        met, met_tlv = self.get_met(event)

        bjets = []
        all_jet_indexes = []
        for ijet, jet in enumerate(jets):
            if (jet.puId < 4 and eval("jet.pt%s" % self.jet_syst) <= 50) or jet.jetId < 2:
                continue
            jet_tlv = ROOT.TLorentzVector()
            jet_tlv.SetPtEtaPhiM(
                eval("jet.pt%s" % self.jet_syst),
                jet.eta,
                jet.phi,
                eval("jet.mass%s" % self.jet_syst)
            )
            if abs(jet_tlv.DeltaR(dau1_tlv)) < 0.5 or abs(jet.DeltaR(dau2_tlv)) < 0.5:
                continue
            if eval("jet.pt%s" % self.jet_syst) > 20 and abs(jet.eta) < 2.4:
                bjets.append((ijet, jet))
            # store also jets w/ eta < 4.7 for vbf analysis
            if eval("jet.pt%s" % self.jet_syst) > 20 and abs(jet.eta) < 4.7:
                all_jet_indexes.append(ijet)

        if len(bjets) < 2:
            return False

        bjets.sort(key=lambda x: x[1]["btag"+self.btag_algo], reverse=True)
        htt_tlv = dau1_tlv + dau2_tlv

        HHbtag_jet_pt_ = ROOT.vector(float)()
        HHbtag_jet_eta_ = ROOT.vector(float)()
        HHbtag_rel_jet_M_pt_ = ROOT.vector(float)()
        HHbtag_rel_jet_E_pt_ = ROOT.vector(float)()
        HHbtag_jet_htt_deta_ = ROOT.vector(float)()
        HHbtag_jet_htt_dphi_ = ROOT.vector(float)()
        HHbtag_jet_btagOutput_ = ROOT.vector(float)()

        for jet in bjets:
            jet_tlv = ROOT.TLorentzVector()
            jet_tlv.SetPtEtaPhiM(
                eval("jet[1].pt%s" % self.jet_syst),
                jet[1].eta,
                jet[1].phi,
                eval("jet[1].mass%s" % self.jet_syst)
            )

            HHbtag_jet_pt_.push_back(jet_tlv.Pt())
            HHbtag_jet_eta_.push_back(jet_tlv.Eta())
            HHbtag_rel_jet_M_pt_.push_back(jet_tlv.M() / jet_tlv.Pt())
            HHbtag_rel_jet_E_pt_.push_back(jet_tlv.E() / jet_tlv.Pt())
            HHbtag_jet_htt_deta_.push_back(htt_tlv.Eta() - jet_tlv.Eta())
            HHbtag_jet_htt_dphi_.push_back(ROOT.Math.VectorUtil.DeltaPhi(htt_tlv, jet_tlv))
            HHbtag_jet_btagOutput_.push_back(jet[1]["btag"+self.btag_algo])

        HHbtag_htt_met_dphi_ = ROOT.Math.VectorUtil.DeltaPhi(htt_tlv, met_tlv)
        HHbtag_htt_scalar_pt_ = dau1.pt + dau2.pt
        HHbtag_rel_met_pt_htt_pt_ = met_tlv.Pt() / HHbtag_htt_scalar_pt_
        HHbtag_htt_pt_ = htt_tlv.Pt()
        HHbtag_htt_eta_ = htt_tlv.Eta()

        HHbtag_evt_ = event.event
        HHbtag_year_ = self.year
        if event.pairType == 0:
            HHbtag_channel_ = 1
        elif event.pairType == 1:
            HHbtag_channel_ = 0
        elif event.pairType == 2:
            HHbtag_channel_ = 1
        else:
            raise ValueError("Pairtype {} is not supported for HHbtag computation".format(
                event.pairType))

        HHbtag_scores = self.HHJets.GetScore(HHbtag_jet_pt_, HHbtag_jet_eta_,
            HHbtag_rel_jet_M_pt_, HHbtag_rel_jet_E_pt_, HHbtag_jet_htt_deta_,
            HHbtag_jet_btagOutput_, HHbtag_jet_htt_dphi_, HHbtag_year_, HHbtag_channel_,
            HHbtag_htt_pt_, HHbtag_htt_eta_, HHbtag_htt_met_dphi_,
            HHbtag_rel_met_pt_htt_pt_, HHbtag_htt_scalar_pt_, HHbtag_evt_)

        HHbtag_scores = list(zip([bjet[0] for bjet in bjets], HHbtag_scores))
        HHbtag_scores.sort(key=lambda x:x[1], reverse=True)  # sort by the obtained HHbtag score

        # 2 "bjets" with the higher HHbtag score are the selected H(bb) candidates
        bjet1_JetIdx = HHbtag_scores[0][0]
        bjet2_JetIdx = HHbtag_scores[1][0]

        if (eval("jets[bjet1_JetIdx].pt%s" % self.jet_syst) <
                eval("jets[bjet2_JetIdx].pt%s" % self.jet_syst)):
            bjet1_JetIdx = HHbtag_scores[1][0]
            bjet2_JetIdx = HHbtag_scores[0][0]

        # let's get the H(bb) object and tlv for the boosted analysis (later)
        bjet1_obj = jets[bjet1_JetIdx]
        bjet2_obj = jets[bjet2_JetIdx]
        bjet1_tlv = ROOT.TLorentzVector()
        bjet2_tlv = ROOT.TLorentzVector()
        bjet1_tlv.SetPtEtaPhiM(eval("bjet1_obj.pt%s" % self.jet_syst), bjet1_obj.eta,
            bjet1_obj.phi, eval("bjet1_obj.mass%s" % self.jet_syst))
        bjet2_tlv.SetPtEtaPhiM(eval("bjet2_obj.pt%s" % self.jet_syst), bjet2_obj.eta,
            bjet2_obj.phi, eval("bjet2_obj.mass%s" % self.jet_syst))

        vbf_jet_pairs = []
        if len(all_jet_indexes) >= 4:
            for i in range(len(all_jet_indexes)):
                if all_jet_indexes[i] in [bjet1_JetIdx, bjet2_JetIdx]:
                    continue
                for j in range(i + 1, len(all_jet_indexes)):
                    if all_jet_indexes[j] in [bjet1_JetIdx, bjet2_JetIdx]:
                        continue

                    jet1_idx = all_jet_indexes[i]
                    jet2_idx = all_jet_indexes[j]
                    if (eval("jets[jet1_idx].pt%s" % self.jet_syst) < 30
                            or eval("jets[jet2_idx].pt%s" % self.jet_syst) < 30
                            or abs(jets[jet1_idx].eta) > 4.7 or abs(jets[jet2_idx].eta) > 4.7):
                        continue
                    vbf_jet_pairs.append(JetPair(jets[jet1_idx], jets[jet2_idx], self.jet_syst,
                        index1=jet1_idx, index2=jet2_idx))
                    # print eval("vbf_jet_pairs[-1].obj1.pt%s" % self.jet_syst), eval("vbf_jet_pairs[-1].obj2.pt%s" % self.jet_syst), vbf_jet_pairs[-1].inv_mass

            if vbf_jet_pairs:
                vbf_pair = max(vbf_jet_pairs)

        self.out.fillBranch("bjet1_JetIdx", bjet1_JetIdx)
        self.out.fillBranch("bjet2_JetIdx", bjet2_JetIdx)

        if vbf_jet_pairs:
            if vbf_pair.obj1.pt > vbf_pair.obj2.pt:
                self.out.fillBranch("VBFjet1_JetIdx", vbf_pair.obj1_index)
                self.out.fillBranch("VBFjet2_JetIdx", vbf_pair.obj2_index)
            else:
                self.out.fillBranch("VBFjet1_JetIdx", vbf_pair.obj2_index)
                self.out.fillBranch("VBFjet2_JetIdx", vbf_pair.obj1_index)
        else:
            self.out.fillBranch("VBFjet1_JetIdx", -1)
            self.out.fillBranch("VBFjet2_JetIdx", -1)

        Jet_HHbtag = []
        HHbtag_scores = dict(HHbtag_scores)  # so it's easier to check the index from the good jets
        for i in range(event.nJet):
            Jet_HHbtag.append(HHbtag_scores[i] if i in HHbtag_scores.keys() else -999.)

        self.out.fillBranch("Jet_HHbtag", Jet_HHbtag)

        # is the event boosted?
        # we loop over the fat AK8 jets, apply a mass cut and verify that its subjets match
        # the jets we selected before.
        is_boosted = 0
        for ifatjet, fatjet in enumerate(fatjets):
            if fatjet.msoftdrop < 30:
                continue
            if fatjet.subJetIdx1 == -1 or fatjet.subJetIdx2 == -1:
                continue
            subj1 = subjets[fatjet.subJetIdx1]
            subj2 = subjets[fatjet.subJetIdx2]
            subj1_tlv = ROOT.TLorentzVector()
            subj2_tlv = ROOT.TLorentzVector()
            subj1_tlv.SetPtEtaPhiM(subj1.pt, subj1.eta, subj1.phi, subj1.mass)
            subj2_tlv.SetPtEtaPhiM(subj2.pt, subj2.eta, subj2.phi, subj2.mass)
            if ((abs(bjet1_tlv.DeltaR(subj1_tlv)) > 0.4
                    or abs(bjet2_tlv.DeltaR(subj2_tlv)) > 0.4)
                and
                (abs(bjet1_tlv.DeltaR(subj2_tlv)) > 0.4
                    or abs(bjet2_tlv.DeltaR(subj1_tlv)) > 0.4)):
                continue
            is_boosted = 1
        self.out.fillBranch("isBoosted", is_boosted)
        return True


def HHJets(**kwargs):
    return lambda: HHJetsProducer(**kwargs)


class HHJetsRDFProducer(JetLepMetSyst):
    def __init__(self, df_filter, *args, **kwargs):
        extEtaAcc = kwargs.pop("extEtaAcc")
        self.btag_algo = kwargs.pop("btag_algo")
        super(HHJetsRDFProducer, self).__init__(self, *args, **kwargs)

        self.df_filter = df_filter

        if not os.getenv("_HHJets"):
            os.environ["_HHJets"] = "HHJets"

            if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libToolsTools.so")
            if os.path.expandvars("$CMT_SCRAM_ARCH") == "slc7_amd64_gcc10":
                ROOT.gROOT.ProcessLine(".include /cvmfs/cms.cern.ch/slc7_amd64_gcc10/"
                    "external/eigen/d812f411c3f9-cms/include/")
                ROOT.gROOT.ProcessLine(".include /cvmfs/cms.cern.ch/slc7_amd64_gcc10/external/"
                    "tensorflow/2.5.0/include/")
            elif s.path.expandvars("$CMT_SCRAM_ARCH") == "slc7_amd64_gcc820":
                ROOT.gROOT.ProcessLine(".include /cvmfs/cms.cern.ch/slc7_amd64_gcc820/"
                    "external/eigen/d812f411c3f9-bcolbf/include/eigen3")
                ROOT.gROOT.ProcessLine(".include /cvmfs/cms.cern.ch/slc7_amd64_gcc820/"
                    "external/tensorflow/2.1.0-bcolbf/include")
            else:
                raise ValueError("Architecture not considered")

            base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))

            ROOT.gROOT.ProcessLine(".L {}/interface/HHJetsInterface.h".format(base))

            self.year = kwargs.pop("year")
            base_hhbtag = "{}/{}/src/HHTools/HHbtag".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            models = [base_hhbtag + "/models/HHbtag_v1_par_%i" % i for i in range(2)]

            ROOT.gInterpreter.Declare("""
                auto HHJets = HHJetsInterface("%s", "%s", %s, %i);
            """ % (models[0], models[1], int(self.year), extEtaAcc))

            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                using VInt = const ROOT::RVec<int>&;
                output get_hh_jets (
                    unsigned long long int event,
                    Vfloat Jet_pt, Vfloat Jet_eta, Vfloat Jet_phi, Vfloat Jet_mass,
                    VInt Jet_puId, Vfloat Jet_jetId, Vfloat Jet_btag%s,
                    Vfloat SubJet_pt, Vfloat SubJet_eta, Vfloat SubJet_phi, Vfloat SubJet_mass,
                    Vfloat FatJet_msoftdrop, VInt FatJet_subJetIdx1, VInt FatJet_subJetIdx2,
                    int pairType, int dau1_index, int dau2_index,
                    Vfloat muon_pt, Vfloat muon_eta, Vfloat muon_phi, Vfloat muon_mass,
                    Vfloat electron_pt, Vfloat electron_eta, Vfloat electron_phi, Vfloat electron_mass,
                    Vfloat tau_pt, Vfloat tau_eta, Vfloat tau_phi, Vfloat tau_mass,
                    float met_pt, float met_phi
                )
                {
                    float dau1_pt, dau1_eta, dau1_phi, dau1_mass, dau2_pt, dau2_eta, dau2_phi, dau2_mass;
                    if (pairType == 0) {
                        dau1_pt = muon_pt.at(dau1_index);
                        dau1_eta = muon_eta.at(dau1_index);
                        dau1_phi = muon_phi.at(dau1_index);
                        dau1_mass = muon_mass.at(dau1_index);
                    } else if (pairType == 1) {
                        dau1_pt = electron_pt.at(dau1_index);
                        dau1_eta = electron_eta.at(dau1_index);
                        dau1_phi = electron_phi.at(dau1_index);
                        dau1_mass = electron_mass.at(dau1_index);
                    } else if (pairType == 2) {
                        dau1_pt = tau_pt.at(dau1_index);
                        dau1_eta = tau_eta.at(dau1_index);
                        dau1_phi = tau_phi.at(dau1_index);
                        dau1_mass = tau_mass.at(dau1_index);
                    } else {
                        dau1_pt = -999.;
                        dau1_eta = -999.;
                        dau1_phi = -999.;
                        dau1_mass = -999.;
                        dau2_pt = -999.;
                        dau2_eta = -999.;
                        dau2_phi = -999.;
                        dau2_mass = -999.;
                    }
                    if (pairType >= 0) {
                        dau2_pt = tau_pt.at(dau2_index);
                        dau2_eta = tau_eta.at(dau2_index);
                        dau2_phi = tau_phi.at(dau2_index);
                        dau2_mass = tau_mass.at(dau2_index);
                    }

                    return HHJets.GetHHJets(event, pairType,
                        Jet_pt, Jet_eta, Jet_phi, Jet_mass,
                        Jet_puId, Jet_jetId, Jet_btag%s,
                        SubJet_pt, SubJet_eta, SubJet_phi, SubJet_mass,
                        FatJet_msoftdrop, FatJet_subJetIdx1, FatJet_subJetIdx2,
                        dau1_pt, dau1_eta, dau1_phi, dau1_mass,
                        dau2_pt, dau2_eta, dau2_phi, dau2_mass,
                        met_pt, met_phi);
                }
            """ % (self.btag_algo, self.btag_algo))

    def run(self, df):
        df = df.Define("HHJets", "get_hh_jets(event, "
            "Jet_pt{5}, Jet_eta, Jet_phi, Jet_mass{5}, "
            "Jet_puId, Jet_jetId, Jet_btag{6}, "
            "SubJet_pt, SubJet_eta, SubJet_phi, SubJet_mass, "
            "FatJet_msoftdrop, FatJet_subJetIdx1, FatJet_subJetIdx2, "
            "pairType, dau1_index, dau2_index, "
            "Muon_pt{0}, Muon_eta, Muon_phi, Muon_mass{0}, "
            "Electron_pt{1}, Electron_eta, Electron_phi, Electron_mass{1}, "
            "Tau_pt{2}, Tau_eta, Tau_phi, Tau_mass{2}, "
            "MET{4}_pt{3}, MET{4}_phi{3})".format(
                self.muon_syst, self.electron_syst, self.tau_syst, self.met_syst,
                self.met_smear_tag, self.jet_syst, self.btag_algo))

        df = df.Define("Jet_HHbtag", "HHJets.hhbtag")
        df = df.Define("bjet1_JetIdx", "HHJets.bjet_idx1")
        df = df.Define("bjet2_JetIdx", "HHJets.bjet_idx2")

        df = df.Define("VBFjet1_JetIdx", "HHJets.vbfjet_idx1")
        df = df.Define("VBFjet2_JetIdx", "HHJets.vbfjet_idx2")

        df = df.Define("ctjet_indexes", "HHJets.ctjet_indexes")
        df = df.Define("fwjet_indexes", "HHJets.fwjet_indexes")

        df = df.Define("isBoosted", "HHJets.isBoosted")

        if self.df_filter:
            df = df.Filter("bjet1_JetIdx >= 0")
        return df, ["Jet_HHbtag", "bjet1_JetIdx", "bjet2_JetIdx",
            "VBFjet1_JetIdx", "VBFjet2_JetIdx", "ctjet_indexes", "fwjet_indexes", "isBoosted"]


def HHJetsRDF(**kwargs):
    """
    Returns the HHbtag output, the indexes from the 2 bjets and 2 vbfjets (if existing),
    the indexes of the additional central and forward jets (if existing) and if the
    event has a boosted topology.

    Lepton and jet systematics (used for pt and mass variables) can be modified using the parameters
    from :ref:`BaseModules_JetLepMetSyst`.

    :param bpair_filter: whether to bpair_filter out output events if they don't have 2 bjet candidates
    :type bpair_filter: bool

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HHJetsRDF
            path: Tools.Tools.HHJets
            parameters:
                year: self.config.year
                isMC: self.dataset.process.isMC
                extEtaAcc: True/False
                bpair_filter: True/False

    """
    df_filter = kwargs.pop("bpair_filter")
    return lambda: HHJetsRDFProducer(df_filter=df_filter, **kwargs)
