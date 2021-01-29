from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import os


class recoilCorrectionProducer(Module):
    def __init__(self):
        if "/AppEffSF.so" not in ROOT.gSystem.GetLibraries():
            print("Load C++ Worker")
            ROOT.gROOT.ProcessLine(
                ".L %s/src/PhysicsTools/NanoAODTools/python/postprocessing/helpers/RecoilCorrector.cc+"
                % os.environ['CMSSW_BASE'])

    def beginJob(self):
        self.rcMain = ROOT.RecoilCorrector("", "")
        self.rcMain.loadRooWorkspacesMCtoCorrect("/afs/cern.ch/work/a/arapyan/public/lowpu_data/Recoil/ZmmMC_PF_13TeV_2G/")
        self.rcMain.loadRooWorkspacesData("/afs/cern.ch/work/a/arapyan/public/lowpu_data/Recoil/ZmmData_PF_13TeV_2G_bkg_fixRoch/")
        self.rcMain.loadRooWorkspacesMC("/afs/cern.ch/work/a/arapyan/public/lowpu_data/Recoil/ZmmMC_PF_13TeV_2G/")

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("met_raw_pt",  "F")
        self.out.branch("met_raw_phi", "F")

        self.out.branch("met_wlepcorr_pt",  "F")
        self.out.branch("met_wlepcorr_phi", "F")

        self.out.branch("met_corrected_pt", "F")
        self.out.branch("met_corrected_phi","F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        lep1_pt = event.lep1_corrected_pt
        lep1_phi = event.eval('lep1.Phi()')
        lep2_pt = event.lep2_corrected_pt
        lep2_phi = event.eval('lep2.Phi()')
        if event.isMuChannel:
            lep1_rawpt = event.eval('lep1.Pt()')
            lep1_rawphi = lep1_phi
            lep2_rawpt = event.eval('lep2.Pt()')
            lep2_rawphi = lep2_phi
        else:
            lep1_rawpt = event.eval('lep1_raw.Pt()')
            lep1_rawphi = event.eval('lep1_raw.Phi()')
            lep2_rawpt = event.eval('lep2_raw.Pt()')
            lep2_rawphi = event.eval('lep2_raw.Phi()')

        vmet_wlepcorr = self.rcMain.CorrectfromLepPt( event.met, event.metPhi, lep1_pt, lep1_phi, lep1_rawpt, lep1_rawphi, lep2_pt, lep2_phi, lep2_rawpt, lep2_rawphi )

        # met with lepton pt corrections
        met_wlepcorr = vmet_wlepcorr.Mod()
        metPhi_wlepcorr = vmet_wlepcorr.Phi()

        if not event.isData:
            # corrected lepton pts
            vlep1 = ROOT.TVector2()
            vlep1.SetMagPhi(lep1_pt, lep1_phi)
            vlep2 = ROOT.TVector2()
            vlep2.SetMagPhi(lep2_pt, lep2_phi)
            vV = vlep1 + vlep2

            vmet_recoilcorrected = self.rcMain.CorrectInvCdf(met_wlepcorr, metPhi_wlepcorr, event.genVPt, event.genVPhi, vV.Mod(), vV.Phi(), 0, 0, 0, False, False)

        else:
            vmet_recoilcorrected = vmet_wlepcorr

        self.out.fillBranch("met_raw_pt",  event.met)
        self.out.fillBranch("met_raw_phi", event.metPhi)

        self.out.fillBranch("met_wlepcorr_pt",  met_wlepcorr)
        self.out.fillBranch("met_wlepcorr_phi", ROOT.TVector2.Phi_mpi_pi(metPhi_wlepcorr))

        self.out.fillBranch("met_corrected_pt", vmet_recoilcorrected.Mod())
        self.out.fillBranch("met_corrected_phi", ROOT.TVector2.Phi_mpi_pi(vmet_recoilcorrected.Phi()))

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

recoilCorrectionModuleConstr = lambda: recoilCorrectionProducer()
