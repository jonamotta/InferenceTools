import os

from analysis_tools.utils import import_root

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from Tools.Tools.tau_utils import LeptonTauPair, TriggerChecker, lepton_veto
from Base.Modules.baseModules import JetLepMetSyst, JetLepMetModule

ROOT = import_root()

class HHLeptonProducer(JetLepMetModule):
    def __init__(self, isMC, year, runPeriod, *args, **kwargs):
        super(HHLeptonProducer, self).__init__(*args, **kwargs)
        self.isMC = isMC
        self.year = year
        self.runPeriod = runPeriod
        self.trigger_checker = TriggerChecker(year)

        if self.year == 2016:
            self.trigger_checker.mutau_triggers = ["HLT_IsoMu22", "HLT_IsoMu22_eta2p1",
                "HLT_IsoTkMu22", "HLT_IsoTkMu22_eta2p1"]
            self.trigger_checker.mutau_crosstriggers = ["HLT_IsoMu19_eta2p1_LooseIsoPFTau20",
                "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1"]
            self.trigger_checker.etau_triggers = ["HLT_Ele25_eta2p1_WPTight_Gsf"]
            self.trigger_checker.etau_crosstriggers = []
            if not self.isMC:
                if self.runPeriod != "H":
                    self.trigger_checker.tautau_triggers = [
                        "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg"]
                else:
                    self.trigger_checker.tautau_triggers = [
                        "HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg"]
            else:
                self.trigger_checker.tautau_triggers = [
                    "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg",
                    "HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg"]
            self.trigger_checker.vbf_triggers = []

        elif self.year == 2017:
            self.trigger_checker.mutau_triggers = ["HLT_IsoMu24", "HLT_IsoMu27"]
            self.trigger_checker.mutau_crosstriggers = [
                "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1"]
            self.trigger_checker.etau_triggers = ["HLT_Ele32_WPTight_Gsf_L1DoubleEG",
                "HLT_Ele32_WPTight_Gsf", "HLT_Ele35_WPTight_Gsf"]
            self.trigger_checker.etau_crosstriggers = [
                "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1"]
            self.trigger_checker.tautau_triggers = [
                "HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
                "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
                "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg"]
            if self.runPeriod in ["D", "E", "F"] or self.isMC:
                self.trigger_checker.vbf_triggers = [
                    "HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg"]
            else:
                self.trigger_checker.vbf_triggers = []

        elif self.year == 2018:
            self.trigger_checker.mutau_triggers = ["HLT_IsoMu24", "HLT_IsoMu27"]
            # lista = lambda e: ([1, 1] if e == 0 else [2, 2])
            self.trigger_checker.mutau_crosstriggers = lambda e: (
                ["HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1"]
                    if (e.run < 317509 and not self.isMC)
                    else ["HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1"])
            self.trigger_checker.etau_triggers = ["HLT_Ele32_WPTight_Gsf", "HLT_Ele35_WPTight_Gsf"]
            self.trigger_checker.etau_crosstriggers = lambda e: (
                ["HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1"]
                    if (e.run < 317509 and not self.isMC)
                    else ["HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1"])
            self.trigger_checker.tautau_triggers = lambda e: (
                ["HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
                "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
                "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg"]
                    if (e.run < 317509 and not self.isMC)
                    else ["HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg"])
            self.trigger_checker.vbf_triggers = lambda e: (
                ["HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1"]
                    if (e.run < 317509 and not self.isMC)
                    else ["HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1"])
                
        pass

    #def beginJob(self):
    #    pass

    #def endJob(self):
    #    pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch('pairType', 'I')
        self.out.branch('dau1_index', 'I')
        self.out.branch('dau2_index', 'I')
        self.out.branch('isVBFtrigger', 'I')
        self.out.branch('isOS', 'I')

        self.out.branch('dau1_eta', 'F')
        self.out.branch('dau1_phi', 'F')
        self.out.branch('dau1_dxy', 'F')
        self.out.branch('dau1_dz', 'F')
        self.out.branch('dau1_q', 'I')
        self.out.branch('dau1_iso', 'F')
        self.out.branch('dau1_decayMode', 'I')
        self.out.branch('dau1_idDeepTau2017v2p1VSe', 'I')
        self.out.branch('dau1_idDeepTau2017v2p1VSmu', 'I')
        self.out.branch('dau1_idDeepTau2017v2p1VSjet', 'I')

        self.out.branch('dau2_eta', 'F')
        self.out.branch('dau2_phi', 'F')
        self.out.branch('dau2_dxy', 'F')
        self.out.branch('dau2_dz', 'F')
        self.out.branch('dau2_q', 'I')
        self.out.branch('dau2_iso', 'F')
        self.out.branch('dau2_decayMode', 'I')
        self.out.branch('dau2_idDeepTau2017v2p1VSe', 'I')
        self.out.branch('dau2_idDeepTau2017v2p1VSmu', 'I')
        self.out.branch('dau2_idDeepTau2017v2p1VSjet', 'I')
        
        self.histo = ROOT.TH1D("InsideHHLepton", "", 21, -1, 20)
        bins = [
            "all", "goodmuon events", "muon-tau pairs", "deltaR > 0.5",
            "pass muon trigger", "pass lepton veto",
            "goodele events", "ele-tau pairs", "deltaR > 0.5",
            "pass ele trigger", "pass lepton veto",
            "goodtau events", "tau-tau pairs",
            "pass tau trigger", "pass lepton veto",
        ]
        for ibin, binname in enumerate(bins):
            self.histo.GetXaxis().SetBinLabel(ibin + 1, binname)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        prevdir = ROOT.gDirectory
        outputFile.cd()
        if "histos" not in [key.GetName() for key in outputFile.GetListOfKeys()]:
            outputFile.mkdir("histos")
        outputFile.cd("histos")
        self.histo.Write()
        prevdir.cd()
    
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        taus = Collection(event, "Tau")
        self.histo.Fill(-1)
        # muon-tau channels
        goodmuons = []
        for imuon, muon in enumerate(muons):
            if (abs(muon.eta) > 2.1 or muon.pfRelIso04_all > 0.15 or abs(muon.dxy) > 0.045
                    or abs(muon.dz) > 0.2 or not muon.tightId):
                continue
            goodmuons.append((imuon, muon))
        if goodmuons:
            self.histo.Fill(0)
            goodtaus = []
            for itau, tau in enumerate(taus):
                if (tau.idDeepTau2017v2p1VSmu < 15 or tau.idDeepTau2017v2p1VSe < 7
                        or tau.idDeepTau2017v2p1VSjet < 1):
                    continue
                if abs(tau.dz) > 0.2:
                    continue
                if tau.decayMode not in [0, 1, 10, 11]:
                    continue
                # common tau pt req for both single and cross triggers
                if eval("tau.pt%s" % self.tau_syst) <= 20.:
                    continue
                goodtaus.append((itau, tau))
            muontaupairs = []
            for (imuon, muon) in goodmuons:
                for (itau, tau) in goodtaus:
                    self.histo.Fill(1)
                    if tau.DeltaR(muon) < 0.5: continue
                    self.histo.Fill(2)
                    if not self.trigger_checker.check_mutau(event,
                            eval("muon.pt%s" % self.muon_syst), muon.eta, muon.phi,
                            eval("tau.pt%s" % self.tau_syst), tau.eta, tau.phi, th1=1, th2=5):
                        continue
                    self.histo.Fill(3)
                    muontaupair = LeptonTauPair(
                        muon, eval("muon.pt%s" % self.muon_syst), muon.pfRelIso04_all,
                        tau, eval("tau.pt%s" % self.tau_syst), tau.rawDeepTau2017v2p1VSjet)
                    # if muontaupair.check_charge():
                    muontaupairs.append((imuon, itau, muontaupair))

            if len(muontaupairs) != 0:
                muontaupairs.sort(key=lambda x: x[2], reverse=True)
                muon, tau = muontaupairs[0][2].pair

                fail_lepton_veto, _ = lepton_veto(electrons, muons, taus, muon)
                if fail_lepton_veto:
                    return False
                self.histo.Fill(4)

                self.out.fillBranch("pairType", 0)
                self.out.fillBranch("isVBFtrigger", 0)
                self.out.fillBranch("isOS", int(muontaupairs[0][2].check_charge()))

                self.out.fillBranch("dau1_index", muontaupairs[0][0])
                self.out.fillBranch("dau1_eta", muon.eta)
                self.out.fillBranch("dau1_phi", muon.phi)
                self.out.fillBranch("dau1_dxy", muon.dxy)
                self.out.fillBranch("dau1_dz", muon.dz)
                self.out.fillBranch("dau1_q", muon.charge)
                self.out.fillBranch("dau1_iso", muon.pfRelIso04_all)
                self.out.fillBranch("dau1_decayMode", -1)
                self.out.fillBranch("dau1_idDeepTau2017v2p1VSe", -1)
                self.out.fillBranch("dau1_idDeepTau2017v2p1VSmu", -1)
                self.out.fillBranch("dau1_idDeepTau2017v2p1VSjet", -1)

                self.out.fillBranch("dau2_index", muontaupairs[0][1])
                self.out.fillBranch("dau2_eta", tau.eta)
                self.out.fillBranch("dau2_phi", tau.phi)
                self.out.fillBranch("dau2_dxy", tau.dxy)
                self.out.fillBranch("dau2_dz", tau.dz)
                self.out.fillBranch("dau2_q", tau.charge)
                self.out.fillBranch("dau2_iso", tau.rawIso)
                self.out.fillBranch("dau2_decayMode", tau.decayMode)
                self.out.fillBranch("dau2_idDeepTau2017v2p1VSe", tau.idDeepTau2017v2p1VSe)
                self.out.fillBranch("dau2_idDeepTau2017v2p1VSmu", tau.idDeepTau2017v2p1VSmu)
                self.out.fillBranch("dau2_idDeepTau2017v2p1VSjet", tau.idDeepTau2017v2p1VSjet)
                return True

        # electron-tau channels
        goodelectrons = []
        for ielectron, electron in enumerate(electrons):
            # if (not (electron.mvaFall17V2Iso_WP80 or electron.mvaFall17V2noIso_WP80)
            if ((not electron.mvaFall17V2Iso_WP80)
                    or abs(electron.dxy) > 0.045 or abs(electron.dz) > 0.2):
                continue
            goodelectrons.append((ielectron, electron))
        if goodelectrons:
            self.histo.Fill(5)
            goodtaus = []
            for itau, tau in enumerate(taus):
                if (tau.idDeepTau2017v2p1VSmu < 15 or tau.idDeepTau2017v2p1VSe < 7
                        or tau.idDeepTau2017v2p1VSjet < 1):
                    continue
                if abs(tau.dz) > 0.2:
                    continue
                if tau.decayMode not in [0, 1, 10, 11]:
                    continue
                # common tau pt req for both single and cross triggers
                if eval("tau.pt%s" % self.tau_syst) <= 20.:
                    continue
                goodtaus.append((itau, tau))

            electrontaupairs = []
            for (ielectron, electron) in goodelectrons:
                for (itau, tau) in goodtaus:
                    self.histo.Fill(6)
                    if tau.DeltaR(electron) < 0.5: continue
                    self.histo.Fill(7)
                    if not self.trigger_checker.check_etau(event,
                            eval("electron.pt%s" % self.electron_syst), electron.eta, electron.phi,
                            eval("tau.pt%s" % self.tau_syst), tau.eta, tau.phi, th1=1, th2=5):
                        continue
                    self.histo.Fill(8)
                    electrontaupair = LeptonTauPair(
                        electron, eval("electron.pt%s" % self.electron_syst), electron.pfRelIso03_all,
                        tau, eval("tau.pt%s" % self.tau_syst), tau.rawDeepTau2017v2p1VSjet)
                    # if electrontaupair.check_charge():
                    electrontaupairs.append((ielectron, itau, electrontaupair))

            if len(electrontaupairs) != 0:
                electrontaupairs.sort(key=lambda x: x[2], reverse=True)
                electron, tau = electrontaupairs[0][2].pair

                fail_lepton_veto, _ = lepton_veto(electrons, muons, taus, electron)
                if fail_lepton_veto:
                    return False
                self.histo.Fill(9)
                self.out.fillBranch("pairType", 1)
                self.out.fillBranch("isVBFtrigger", 0)
                self.out.fillBranch("isOS", int(electrontaupairs[0][2].check_charge()))

                self.out.fillBranch("dau1_index", electrontaupairs[0][0])
                self.out.fillBranch("dau1_eta", electron.eta)
                self.out.fillBranch("dau1_phi", electron.phi)
                self.out.fillBranch("dau1_dxy", electron.dxy)
                self.out.fillBranch("dau1_dz", electron.dz)
                self.out.fillBranch("dau1_q", electron.charge)
                self.out.fillBranch("dau1_iso", electron.pfRelIso03_all)
                self.out.fillBranch("dau1_decayMode", -1)
                self.out.fillBranch("dau1_idDeepTau2017v2p1VSe", -1)
                self.out.fillBranch("dau1_idDeepTau2017v2p1VSmu", -1)
                self.out.fillBranch("dau1_idDeepTau2017v2p1VSjet", -1)

                self.out.fillBranch("dau2_index", electrontaupairs[0][1])
                self.out.fillBranch("dau2_eta", tau.eta)
                self.out.fillBranch("dau2_phi", tau.phi)
                self.out.fillBranch("dau2_dxy", tau.dxy)
                self.out.fillBranch("dau2_dz", tau.dz)
                self.out.fillBranch("dau2_q", tau.charge)
                self.out.fillBranch("dau2_iso", tau.rawIso)
                self.out.fillBranch("dau2_decayMode", tau.decayMode)
                self.out.fillBranch("dau2_idDeepTau2017v2p1VSe", tau.idDeepTau2017v2p1VSe)
                self.out.fillBranch("dau2_idDeepTau2017v2p1VSmu", tau.idDeepTau2017v2p1VSmu)
                self.out.fillBranch("dau2_idDeepTau2017v2p1VSjet", tau.idDeepTau2017v2p1VSjet)

                return True

        goodtaus = []
        for itau, tau in enumerate(taus):
            if (tau.idDeepTau2017v2p1VSmu < 1 or tau.idDeepTau2017v2p1VSe < 3
                    or tau.idDeepTau2017v2p1VSjet < 1):
                continue
            if abs(tau.dz) > 0.2:
                continue
            if tau.decayMode not in [0, 1, 10, 11]:
                continue
            goodtaus.append((itau, tau))

        if goodtaus:
            self.histo.Fill(10)
        tautaupairs = []
        for i in range(len(goodtaus)):
            for j in range(len(goodtaus)):
                if i == j:
                    continue
                self.histo.Fill(11)
                tau1_index = goodtaus[i][0]
                tau1 = goodtaus[i][1]
                tau2_index = goodtaus[j][0]
                tau2 = goodtaus[j][1]

                if tau1.DeltaR(tau2) < 0.5: continue

                pass_ditau = self.trigger_checker.check_tautau(event,
                    eval("tau1.pt%s" % self.tau_syst), tau1.eta, tau1.phi,
                    eval("tau2.pt%s" % self.tau_syst), tau2.eta, tau2.phi, abs_th1=40, abs_th2=40)
                # passing vbf trigger ONLY
                pass_vbf = (not pass_ditau) and self.trigger_checker.check_vbftautau(event,
                    eval("tau1.pt%s" % self.tau_syst), tau1.eta, tau1.phi,
                    eval("tau2.pt%s" % self.tau_syst), tau2.eta, tau2.phi, abs_th1=25, abs_th2=25)

                if not (pass_ditau or pass_vbf):
                    continue
                self.histo.Fill(12)
                pass_vbf = int(pass_vbf)
                tautaupair = LeptonTauPair(
                    tau1, eval("tau1.pt%s" % self.tau_syst), tau1.rawDeepTau2017v2p1VSjet,
                    tau2, eval("tau2.pt%s" % self.tau_syst), tau2.rawDeepTau2017v2p1VSjet)
                # if tautaupair.check_charge():
                tautaupairs.append((tau1_index, tau2_index, tautaupair, pass_vbf))

        if len(tautaupairs) != 0:
            tautaupairs.sort(key=lambda x: x[2], reverse=True)
            tau1, tau2 = tautaupairs[0][2].pair

            fail_lepton_veto, _ = lepton_veto(electrons, muons, taus)
            if fail_lepton_veto:
                return False
            self.histo.Fill(13)

            self.out.fillBranch("pairType", 2)
            self.out.fillBranch("isVBFtrigger", tautaupairs[0][3])
            self.out.fillBranch("isOS", int(tautaupairs[0][2].check_charge()))

            self.out.fillBranch("dau1_index", tautaupairs[0][0])
            self.out.fillBranch("dau1_eta", tau1.eta)
            self.out.fillBranch("dau1_phi", tau1.phi)
            self.out.fillBranch("dau1_dxy", tau1.dxy)
            self.out.fillBranch("dau1_dz", tau1.dz)
            self.out.fillBranch("dau1_q", tau1.charge)
            self.out.fillBranch("dau1_iso", tau1.rawIso)
            self.out.fillBranch("dau1_decayMode", tau1.decayMode)
            self.out.fillBranch("dau1_idDeepTau2017v2p1VSe", tau1.idDeepTau2017v2p1VSe)
            self.out.fillBranch("dau1_idDeepTau2017v2p1VSmu", tau1.idDeepTau2017v2p1VSmu)
            self.out.fillBranch("dau1_idDeepTau2017v2p1VSjet", tau1.idDeepTau2017v2p1VSjet)

            self.out.fillBranch("dau2_index", tautaupairs[0][1])
            self.out.fillBranch("dau2_eta", tau2.eta)
            self.out.fillBranch("dau2_phi", tau2.phi)
            self.out.fillBranch("dau2_dxy", tau2.dxy)
            self.out.fillBranch("dau2_dz", tau2.dz)
            self.out.fillBranch("dau2_q", tau2.charge)
            self.out.fillBranch("dau2_iso", tau2.rawIso)
            self.out.fillBranch("dau2_decayMode", tau2.decayMode)
            self.out.fillBranch("dau2_idDeepTau2017v2p1VSe", tau2.idDeepTau2017v2p1VSe)
            self.out.fillBranch("dau2_idDeepTau2017v2p1VSmu", tau2.idDeepTau2017v2p1VSmu)
            self.out.fillBranch("dau2_idDeepTau2017v2p1VSjet", tau2.idDeepTau2017v2p1VSjet)

            return True
        return False


class HHLeptonVariableProducer(JetLepMetModule):
    def __init__(self, isMC, *args, **kwargs):
        super(HHLeptonVariableProducer, self).__init__(*args, **kwargs)
        self.isMC = isMC

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch('dau1_pt%s' % self.lep_syst, 'F')
        self.out.branch('dau1_mass%s' % self.lep_syst, 'F')
        self.out.branch('dau2_pt%s' % self.lep_syst, 'F')
        self.out.branch('dau2_mass%s' % self.lep_syst, 'F')
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        taus = Collection(event, "Tau")

        _, _, dau1_tlv, dau2_tlv = self.get_daus(event, muons, electrons, taus)
        self.out.fillBranch("dau1_pt%s" % self.lep_syst, dau1_tlv.Pt())
        self.out.fillBranch("dau1_mass%s" % self.lep_syst, dau1_tlv.M())
        self.out.fillBranch("dau2_pt%s" % self.lep_syst, dau2_tlv.Pt())
        self.out.fillBranch("dau2_mass%s" % self.lep_syst, dau2_tlv.M())
        return True


def HHLepton(**kwargs):
    return lambda: HHLeptonProducer(**kwargs)


def HHLeptonVariable(**kwargs):
    return lambda: HHLeptonVariableProducer(**kwargs)


class HHLeptonRDFProducer(JetLepMetSyst):
    def __init__(self, isMC, year, runPeriod, df_filter, *args, **kwargs):
        super(HHLeptonRDFProducer, self).__init__(*args, **kwargs)
        self.isMC = isMC
        self.year = year
        self.runPeriod = runPeriod
        self.df_filter = df_filter
        vvvl_vsjet = kwargs.pop("vvvl_vsjet")
        vl_vse = kwargs.pop("vl_vse")
        vvl_vse = kwargs.pop("vvl_vse")
        t_vsmu = kwargs.pop("t_vsmu")
        vl_vsmu = kwargs.pop("vl_vsmu")

        if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("libToolsTools.so")

        base = "{}/{}/src/Tools/Tools".format(
            os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
        ROOT.gROOT.ProcessLine(".L {}/interface/HHLeptonInterface.h".format(base))
        ROOT.gInterpreter.Declare("""
            auto HHLepton = HHLeptonInterface(%s, %s, %s, %s, %s);
        """ % (vvvl_vsjet, vl_vse, vvl_vse, t_vsmu, vl_vsmu))

        self.mutau_triggers = ["HLT_IsoMu22", "HLT_IsoMu22_eta2p1",
            "HLT_IsoTkMu22", "HLT_IsoTkMu22_eta2p1", "HLT_IsoMu24", "HLT_IsoMu27",
            "HLT_IsoMu19_eta2p1_LooseIsoPFTau20", "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1",
            "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1",
            "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1"]
        self.etau_triggers = ["HLT_Ele25_eta2p1_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
            "HLT_Ele32_WPTight_Gsf", "HLT_Ele35_WPTight_Gsf",
            "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1",
            "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1"]
        self.tautau_triggers = ["HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg",
            "HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg",
            "HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
            "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
            "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg",
            "HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg"]
        self.vbf_triggers = ["HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg",
            "HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1",
            "HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1"]

        if self.year == 2018:
            ROOT.gInterpreter.Declare("""
                using Vbool = const ROOT::RVec<Bool_t>&;  
                std::vector<trig_req> get_mutau_triggers(
                        Vbool triggers, bool isMC, int run, int runPeriod) {
                    std::vector<trig_req> trigger_reqs;
                    trigger_reqs.push_back(trig_req({triggers[4], 25, 2.1, 20, 2.3, {{2, 8}, {}}}));
                    trigger_reqs.push_back(trig_req({triggers[5], 28, 2.1, 20, 2.3, {{2, 8}, {}}}));
                    if (!isMC && run < 317509)
                        trigger_reqs.push_back(trig_req({triggers[8], 21, 2.1, 32, 2.1, {{64}, {1, 256}}}));
                    else
                        trigger_reqs.push_back(trig_req({triggers[9], 21, 2.1, 32, 2.1, {{4}, {1, 16}}}));
                    return trigger_reqs;
                }
                std::vector<trig_req> get_etau_triggers(
                        Vbool triggers, bool isMC, int run, int runPeriod) {
                    std::vector<trig_req> trigger_reqs;
                    trigger_reqs.push_back(trig_req({triggers[2], 33, 2.1, 20, 2.3, {{2}, {}}}));
                    trigger_reqs.push_back(trig_req({triggers[3], 36, 2.1, 20, 2.3, {{2}, {}}}));
                    if (!isMC && run < 317509)
                        trigger_reqs.push_back(trig_req({triggers[4], 25, 2.1, 35, 2.1, {{64}, {1, 128}}}));
                    else
                        trigger_reqs.push_back(trig_req({triggers[5], 25, 2.1, 35, 2.1, {{8}, {1, 16}}}));
                    return trigger_reqs;
                }
                std::vector<trig_req> get_tautau_triggers(
                        Vbool triggers, bool isMC, int run, int runPeriod) {
                    std::vector<trig_req> trigger_reqs;
                    trigger_reqs.push_back(trig_req({triggers[2], 40, 2.1, 40, 2.1, {{64, 4, 8}, {64, 4, 8}}}));
                    trigger_reqs.push_back(trig_req({triggers[3], 40, 2.1, 40, 2.1, {{64, 4, 8}, {64, 4, 8}}}));
                    if (!isMC && run < 317509)
                        trigger_reqs.push_back(trig_req({triggers[4], 40, 2.1, 40, 2.1, {{64, 4}, {64, 4}}}));
                    else
                        trigger_reqs.push_back(trig_req({triggers[5], 40, 2.1, 40, 2.1, {{2, 16}, {2, 16}}}));
                    return trigger_reqs;
                }
                std::vector<trig_req> get_vbf_triggers(
                        Vbool triggers, bool isMC, int run, int runPeriod) {
                    std::vector<trig_req> trigger_reqs;
                    if (!isMC && run < 317509)
                        trigger_reqs.push_back(trig_req({triggers[1], 25, 2.1, 25, 2.1, {{64, 4, 8}, {64, 4, 8}}}));
                    else
                        trigger_reqs.push_back(trig_req({triggers[2], 25, 2.1, 25, 2.1, {{512, 1, 16}, {512, 1, 16}}}));
                    return trigger_reqs;
                }
            """)

    def run(self, df):
        branches = ["pairType", "dau1_index", "dau2_index", "isVBFtrigger", "isOS",
            "dau1_eta", "dau1_phi", "dau1_iso", "dau1_decayMode",
            "dau1_idDeepTau2017v2p1VSe", "dau1_idDeepTau2017v2p1VSmu",
            "dau1_idDeepTau2017v2p1VSjet",
            "dau2_eta", "dau2_phi", "dau2_decayMode",
            "dau2_idDeepTau2017v2p1VSe", "dau2_idDeepTau2017v2p1VSmu",
            "dau2_idDeepTau2017v2p1VSjet"
        ]

        all_branches = df.GetColumnNames()
        for ib, branch in enumerate(self.mutau_triggers):
            if branch not in all_branches:
                self.mutau_triggers[ib] = "false"
        for ib, branch in enumerate(self.etau_triggers):
            if branch not in all_branches:
                self.etau_triggers[ib] = "false"
        for ib, branch in enumerate(self.tautau_triggers):
            if branch not in all_branches:
                self.tautau_triggers[ib] = "false"
        for ib, branch in enumerate(self.vbf_triggers):
            if branch not in all_branches:
                self.vbf_triggers[ib] = "false"

        runPeriods = ["dum", "A", "B", "C", "D", "E", "F", "G", "H"]
        runPeriod = None
        for irun, runPeriod in enumerate(runPeriods):
            if self.runPeriod == runPeriod:
                runPeriod = irun
                break
        assert runPeriod != None

        df = df.Define("mutau_triggers", "get_mutau_triggers({%s}, %s, run, %s)" % (
            ", ".join(self.mutau_triggers), ("true" if self.isMC else "false"), runPeriod))
        df = df.Define("etau_triggers", "get_etau_triggers({%s}, %s, run, %s)" % (
            ", ".join(self.etau_triggers), ("true" if self.isMC else "false"), runPeriod))
        df = df.Define("tautau_triggers", "get_tautau_triggers({%s}, %s, run, %s)" % (
            ", ".join(self.tautau_triggers), ("true" if self.isMC else "false"), runPeriod))
        df = df.Define("vbf_triggers", "get_vbf_triggers({%s}, %s, run, %s)" % (
            ", ".join(self.vbf_triggers), ("true" if self.isMC else "false"), runPeriod))

        df = df.Define("hh_lepton_results", "HHLepton.get_dau_indexes("
            "Muon_pt{0}, Muon_eta, Muon_phi, Muon_mass{0}, "
            "Muon_pfRelIso04_all, Muon_dxy, Muon_dz, Muon_mediumId, Muon_tightId, Muon_charge, "
            "Electron_pt{1}, Electron_eta, Electron_phi, Electron_mass{1}, "
            "Electron_mvaFall17V2Iso_WP80, Electron_mvaFall17V2noIso_WP90, "
            "Electron_mvaFall17V2Iso_WP90, Electron_pfRelIso03_all, "
            "Electron_dxy, Electron_dz, Electron_charge, "
            "Tau_pt{2}, Tau_eta, Tau_phi, Tau_mass{2}, "
            "Tau_idDeepTau2017v2p1VSmu, Tau_idDeepTau2017v2p1VSe, "
            "Tau_idDeepTau2017v2p1VSjet, Tau_rawDeepTau2017v2p1VSjet, "
            "Tau_dz, Tau_decayMode, Tau_charge, "
            "TrigObj_id, TrigObj_filterBits, TrigObj_eta, TrigObj_phi, "
            "mutau_triggers, etau_triggers, tautau_triggers, vbf_triggers"
        ")".format(self.muon_syst, self.electron_syst, self.tau_syst))

        for branch in branches:
            df = df.Define(branch, "hh_lepton_results.%s" % branch)

        if self.df_filter:
            df = df.Filter("pairType >= 0")

        return df, branches
        # return df, []


def HHLeptonRDF(**kwargs):
    """
    Returns the index of the two selected taus + several of their variables not affected by
    systematics.

    Lepton systematics (used for pt and mass variables) can be modified using the parameters from 
    :ref:`BaseModules_JetLepMetSyst`.

    :param runPeriod: run period in caps (data only)
    :type runPeriod: str

    :param vvvl_vsjet: VVVLoose DeepTauVSjet WP value
    :type vvvl_vsjet: int

    :param vl_vse: VLoose DeepTauVSe WP value
    :type vl_vse: int

    :param vvl_vse: VVLoose DeepTaVSe WP value
    :type vvl_vse: int

    :param t_vsmu: Tight DeepTauVSmu WP value
    :type t_vsmu: int

    :param vl_vsmu: VLoose DeepTauVSmu WP value
    :type vl_vsmu: int

    :param filter: whether to filter out output events if they don't have 2 lepton candidates
    :type filter: bool

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HHLeptonRDF
            path: Tools.Tools.HHLepton
            parameters:
                isMC: self.dataset.process.isMC
                year: self.config.year
                runPeriod: self.dataset.runPeriod
                vvvl_vsjet: self.config.deeptau.vsjet.VVVLoose
                vl_vse: self.config.deeptau.vse.VLoose
                vvl_vse: self.config.deeptau.vse.VVLoose
                t_vsmu: self.config.deeptau.vsmu.Tight
                vl_vsmu: self.config.deeptau.vsmu.VLoose
                filter: True

    """
    df_filter = kwargs.pop("filter")
    return lambda: HHLeptonRDFProducer(df_filter=df_filter, **kwargs)


class HHLeptonVarRDFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        super(HHLeptonVarRDFProducer, self).__init__(*args, **kwargs)
        if not os.getenv("DAU_VAR"):
            os.environ["DAU_VAR"] = "true"
            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                using VInt = const ROOT::RVec<int>&;
                std::vector<float> get_lepton_values (
                    int pairType, int dau1_index, int dau2_index,
                    Vfloat muon_pt, Vfloat muon_mass,
                    Vfloat electron_pt, Vfloat electron_mass,
                    Vfloat tau_pt, Vfloat tau_mass
                )
                {
                    float dau1_pt, dau1_mass, dau2_pt, dau2_mass;
                    if (pairType == 0) {
                        dau1_pt = muon_pt.at(dau1_index);
                        dau1_mass = muon_mass.at(dau1_index);
                    } else if (pairType == 1) {
                        dau1_pt = electron_pt.at(dau1_index);
                        dau1_mass = electron_mass.at(dau1_index);
                    } else if (pairType == 2) {
                        dau1_pt = tau_pt.at(dau1_index);
                        dau1_mass = tau_mass.at(dau1_index);
                    } else {
                        return {-999., -999., -999., -999.};
                    }
                    dau2_pt = tau_pt.at(dau2_index);
                    dau2_mass = tau_mass.at(dau2_index);
                    return {dau1_pt, dau1_mass, dau2_pt, dau2_mass};
                }
            """)


    def run(self, df):
        branches = "dau1_pt{0}, dau1_mass{0}, dau2_pt{0}, dau2_mass{0}".format(self.lep_syst)
        branches = branches.split(", ")

        df = df.Define("lepton_values%s" % self.lep_syst, "get_lepton_values("
            "pairType, dau1_index, dau2_index, "
            "Muon_pt{0}, Muon_mass{0}, Electron_pt{1}, Electron_mass{1}, "
            "Tau_pt{2}, Tau_mass{2})".format(self.muon_syst, self.electron_syst, self.tau_syst))

        for ib, branch in enumerate(branches):
            df = df.Define(branch, "lepton_values%s[%s]" % (self.lep_syst, ib))

        return df, branches


def HHLeptonVarRDF(**kwargs):
    """
    Returns the pt and mass of the selected taus possibly including systematics

    Lepton systematics (used for pt and mass variables) can be modified using the parameters from 
    :ref:`BaseModules_JetLepMetSyst`.

    Required RDFModules: :ref:`HHLepton_HHLeptonRDF`.

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HHLeptonVarRDF
            path: Tools.Tools.HHLepton
            parameters:
                isMC: self.dataset.process.isMC

    """
    return lambda: HHLeptonVarRDFProducer(**kwargs)
