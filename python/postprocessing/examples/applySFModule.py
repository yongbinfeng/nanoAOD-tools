from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import os


class applySFProducer(Module):
    def __init__(self):
        if "/AppEffSF.so" not in ROOT.gSystem.GetLibraries():
            print("Load C++ Worker")
            ROOT.gROOT.ProcessLine(
                ".L %s/src/PhysicsTools/NanoAODTools/python/postprocessing/helpers/AppEffSF.cc+"
                % os.environ['CMSSW_BASE'])
        
        self.ELE_MASS = 0.000511
        self.MU_MASS = 0.1057

    def beginJob(self):
        self.effs_mu = ROOT.AppEffSF('/afs/cern.ch/work/a/arapyan/public/lowpu_data/Efficiency/lowpu_13TeV/results/Zmm/')
        self.effs_mu.loadHLT("MuHLTEff_aMCxPythia","Positive","Negative");
        self.effs_mu.loadSel("MuSITEff_aMCxPythia","Combined","Combined");
        self.effs_mu.loadSta("MuStaEff_aMCxPythia","Combined","Combined");
        self.effs_mu.loadUncSel("/afs/cern.ch/work/a/arapyan/public/lowpu_data/Efficiency/lowpu_13TeV/Systematics/SysUnc_MuSITEff.root");
        self.effs_mu.loadUncSta("/afs/cern.ch/work/a/arapyan/public/lowpu_data/Efficiency/lowpu_13TeV/Systematics/SysUnc_MuStaEff.root");

        self.effs_ele = ROOT.AppEffSF('/afs/cern.ch/work/a/arapyan/public/lowpu_data/Efficiency/lowpu_13TeV/results/Zee/')
        self.effs_ele.loadHLT("EleHLTEff_aMCxPythia","Positive","Negative");
        self.effs_ele.loadSel("EleGSFSelEff_aMCxPythia","Combined","Combined");
        self.effs_ele.loadUncSel("/afs/cern.ch/work/a/arapyan/public/lowpu_data/Efficiency/lowpu_13TeV/Systematics/SysUnc_EleGSFSelEff.root");

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("lepsfweight", "F")
        self.out.branch("lepsfweight_hlt", "F")
        self.out.branch("lepsfweight_fsr", "F")
        self.out.branch("lepsfweight_mc",  "F")
        self.out.branch("lepsfweight_tag", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        wmain = 1.0
        w_hlt = 1.0
        wfsr = 1.0
        wmc = 1.0
        wbkg = 1.0
        wtag = 1.0

        if not event.isData:
            l1 = ROOT.TLorentzVector()
            l1.SetPtEtaPhiM(event.lep1_corrected_pt, event.eval('lep1.Eta()'), event.eval('lep1.Phi()'), self.MU_MASS if event.isMuChannel else self.ELE_MASS)
            l2 = ROOT.TLorentzVector()
            l2.SetPtEtaPhiM(event.lep2_corrected_pt, event.eval('lep2.Eta()'), event.eval('lep2.Phi()'), self.MU_MASS if event.isMuChannel else self.ELE_MASS)
            uncs_sta = [1.0, 1.0, 1.0, 1.0]

            if event.isMuChannel:
                corr = self.effs_mu.fullEfficiencies(l1, event.q1, l2, event.q2)
                uncs_sta = self.effs_mu.getUncSta(   l1, event.q1, l2, event.q2)
                uncs_sit = self.effs_mu.getUncSel(   l1, event.q1, l2, event.q2)
                w_hlt *= self.effs_mu.computeHLTSF(l1,event.q1, l2, event.q2)
            else:
                corr = self.effs_ele.fullEfficiencies(l1, event.q1, l2, event.q2)
                uncs_sit = self.effs_ele.getUncSel(   l1, event.q1, l2, event.q2)
                w_hlt *= self.effs_ele.computeHLTSF(l1,event.q1, l2, event.q2)

            wmain *= corr
            wfsr *= uncs_sta[0]*uncs_sit[0]*w_hlt
            wmc  *= uncs_sta[1]*uncs_sit[1]*w_hlt
            wbkg *= uncs_sta[2]*uncs_sit[2]*w_hlt
            wtag *= uncs_sta[3]*uncs_sit[3]*w_hlt

        self.out.fillBranch("lepsfweight", wmain)
        self.out.fillBranch("lepsfweight_hlt", w_hlt)
        self.out.fillBranch("lepsfweight_fsr", wfsr)
        self.out.fillBranch("lepsfweight_mc", wmc)
        self.out.fillBranch("lepsfweight_tag", wtag)
        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

applySFModuleConstr = lambda: applySFProducer()
