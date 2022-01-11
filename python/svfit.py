import os
from array import array
import math

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from Base.Modules.baseModules import JetLepMetSyst, JetLepMetModule
from analysis_tools.utils import import_root


ROOT = import_root()

class SVFitProducer(JetLepMetModule):
    def __init__(self, *args, **kwargs):
        super(SVFitProducer, self).__init__(*args, **kwargs)
        if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("libToolsTools.so")

        base = "{}/{}/src/Tools/Tools".format(
            os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
        ROOT.gROOT.ProcessLine(".L {}/interface/SVfitinterface.h".format(base))
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch('Htt_svfit_pt%s' % self.systs, 'F')
        self.out.branch('Htt_svfit_eta%s' % self.systs, 'F')
        self.out.branch('Htt_svfit_phi%s' % self.systs, 'F')
        self.out.branch('Htt_svfit_mass%s' % self.systs, 'F')
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        taus = Collection(event, "Tau")

        dau1, dau2, dau1_tlv, dau2_tlv = self.get_daus(event, muons, electrons, taus)
        met, met_tlv = self.get_met(event)

        decayMode1 = (dau1.decayMode if dau1 in taus else -1)
        decayMode2 = (dau2.decayMode if dau2 in taus else -1)

        svfit = ROOT.SVfitinterface(
            0, event.pairType, decayMode1, decayMode2,
            dau1_tlv.Pt(), dau1_tlv.Eta(), dau1_tlv.Phi(), dau1_tlv.M(),
            dau2_tlv.Pt(), dau2_tlv.Eta(), dau2_tlv.Phi(), dau2_tlv.M(),
            met_tlv.Pt(), met_tlv.Phi(), met.covXX, met.covXY, met.covYY
        )

        result = svfit.FitAndGetResult()

        self.out.fillBranch("Htt_svfit_pt%s" % self.systs, result[0])
        self.out.fillBranch("Htt_svfit_eta%s" % self.systs, result[1])
        self.out.fillBranch("Htt_svfit_phi%s" % self.systs, result[2])
        self.out.fillBranch("Htt_svfit_mass%s" % self.systs, result[3])
        return True


class SVFitRDFProducer(JetLepMetSyst):
    def run(self, df):
        if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("libToolsTools.so")
        base = "{}/{}/src/Tools/Tools".format(
            os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
        ROOT.gROOT.ProcessLine(".L {}/interface/SVfitinterface.h".format(base))
        ROOT.gInterpreter.Declare("""
            // #include "Tools/Tools/interface/SVfitinterface.h"
            using Vfloat = const ROOT::RVec<float>&;
            ROOT::RVec<double> compute_svfit(
                    int pairType, int dau1_index, int dau2_index, int DM1, int DM2,
                    Vfloat muon_pt, Vfloat muon_eta, Vfloat muon_phi, Vfloat muon_mass,
                    Vfloat electron_pt, Vfloat electron_eta, Vfloat electron_phi, Vfloat electron_mass,
                    Vfloat tau_pt, Vfloat tau_eta, Vfloat tau_phi, Vfloat tau_mass,
                    float met_pt, float met_phi, float met_covXX, float met_covXY, float met_covYY) {
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
                }
                dau2_pt = tau_pt.at(dau2_index);
                dau2_eta = tau_eta.at(dau2_index);
                dau2_phi = tau_phi.at(dau2_index);
                dau2_mass = tau_mass.at(dau2_index);

                auto svfit = SVfitinterface();
                svfit.SetInputs(0, pairType, DM1, DM2,
                    dau1_pt, dau1_eta, dau1_phi, dau1_mass,
                    dau2_pt, dau2_eta, dau2_phi, dau2_mass,
                    met_pt, met_phi, met_covXX, met_covXY, met_covYY);

                return svfit.FitAndGetResult();
            }
        """)

        df = df.Define("svfit_result",
            "compute_svfit(pairType, dau1_index, dau2_index, dau1_decayMode, dau2_decayMode, "
                "Muon_pt{0}, Muon_eta, Muon_phi, Muon_mass{0}, "
                "Electron_pt{1}, Electron_eta, Electron_phi, Electron_mass{1}, "
                "Tau_pt{2}, Tau_eta, Tau_phi, Tau_mass{2}, "
                "MET{4}_pt{3}, MET{4}_phi{3}, MET_covXX, MET_covXY, MET_covYY)".format(
                    self.muon_syst, self.electron_syst, self.tau_syst, self.met_syst,
                    self.met_smear_tag)).Define(
            "Htt_svfit_pt%s" % self.systs, "svfit_result[0]").Define(
            "Htt_svfit_eta%s" % self.systs, "svfit_result[1]").Define(
            "Htt_svfit_phi%s" % self.systs, "svfit_result[2]").Define(
            "Htt_svfit_mass%s" % self.systs, "svfit_result[3]")
        return df, ["Htt_svfit_pt%s" % self.systs, "Htt_svfit_eta%s" % self.systs,
            "Htt_svfit_phi%s" % self.systs, "Htt_svfit_mass%s" % self.systs]


def SVFit(**kwargs):
    return lambda: SVFitProducer(**kwargs)
