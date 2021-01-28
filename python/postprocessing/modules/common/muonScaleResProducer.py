from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
import ROOT
import os
import random
ROOT.PyConfig.IgnoreCommandLineOptions = True


def mk_safe(fct, *args):
    try:
        return fct(*args)
    except Exception as e:
        if any('Error in function boost::math::erf_inv' in arg
               for arg in e.args):
            print(
                'WARNING: catching exception and returning -1. Exception arguments: %s'
                % e.args)
            return -1.
        else:
            raise e


class muonScaleResProducer(Module):
    def __init__(self, rc_dir, rc_corrections):
        print("load muon module")
        p_postproc = '%s/src/PhysicsTools/NanoAODTools/python/postprocessing' % os.environ[
            'CMSSW_BASE']
        p_roccor = p_postproc + '/data/' + rc_dir
        if "/RoccoR_cc.so" not in ROOT.gSystem.GetLibraries():
            p_helper = '%s/RoccoR.cc' % p_roccor
            print('Loading C++ helper from ' + p_helper)
            ROOT.gROOT.ProcessLine('.L ' + p_helper)
        self._roccor = ROOT.RoccoR(p_roccor + '/' + rc_corrections)

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("isData", "O")
        self.out.branch("isMuChannel", "O")
        self.out.branch("lep1_corrected_pt",     "F")
        self.out.branch("lep1_correctedUp_pt",   "F")
        self.out.branch("lep1_correctedDown_pt", "F")
        self.out.branch("lep2_corrected_pt",     "F")
        self.out.branch("lep2_correctedUp_pt",   "F")
        self.out.branch("lep2_correctedDown_pt", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        isData = (event.runNum != 1)

        # check if the lepton is muon by checking the lepton mass
        isMuon = (abs(event.eval('lep1.M()') - 0.1056583) < 0.01)
        if not isMuon:
            self.out.fillBranch("isData", isData)
            self.out.fillBranch("isMuChannel", isMuon)
            lep1_pt = event.eval('lep1.Pt()')
            lep2_pt = event.eval('lep2.Pt()')
            self.out.fillBranch("lep1_corrected_pt",     lep1_pt)
            self.out.fillBranch("lep2_corrected_pt",     lep2_pt)
            self.out.fillBranch("lep1_correctedUp_pt",   lep1_pt)
            self.out.fillBranch("lep1_correctedDown_pt", lep1_pt)
            self.out.fillBranch("lep2_correctedUp_pt",   lep2_pt)
            self.out.fillBranch("lep2_correctedDown_pt", lep2_pt)

            return True
            
        mu1 = Object(event, "mu1")
        mu1.pt  = event.eval('lep1.Pt()')
        mu1.eta = event.eval('lep1.Eta()')
        mu1.phi = event.eval('lep1.Phi()')
        mu1.charge = event.q1
        mu1.nTkLayers = event.nTkLayers1

        mu2 = Object(event, "mu2")
        mu2.pt  = event.eval('lep2.Pt()')
        mu2.eta = event.eval('lep2.Eta()')
        mu2.phi = event.eval('lep2.Phi()')
        mu2.charge = event.q2
        mu2.nTkLayers = event.nTkLayers2

        genmu1_pt = event.eval('genlep1.Pt()')
        genmu2_pt = event.eval('genlep2.Pt()')

        roccor = self._roccor
        pt_corr = []
        pt_err = []
        if not isData:
            for mu, genmu_pt in zip([mu1, mu2], [genmu1_pt, genmu2_pt]):
                if genmu_pt>0:
                    pt_corr.append(mu.pt *
                                   mk_safe(roccor.kSpreadMC, mu.charge, mu.pt,
                                           mu.eta, mu.phi, genmu_pt))
                    pt_err.append(mu.pt *
                                  mk_safe(roccor.kSpreadMCerror, mu.charge,
                                          mu.pt, mu.eta, mu.phi, genmu_pt))
                else:
                    u1 = random.uniform(0.0, 1.0)
                    pt_corr.append(
                        mu.pt * mk_safe(roccor.kSmearMC, mu.charge, mu.pt,
                                        mu.eta, mu.phi, mu.nTkLayers, u1))
                    pt_err.append(
                        mu.pt * mk_safe(roccor.kSmearMCerror, mu.charge, mu.pt,
                                        mu.eta, mu.phi, mu.nTkLayers, u1))

        else:
            for mu in [mu1, mu2]:
                pt_corr.append( mu.pt *
                                mk_safe(roccor.kScaleDT, mu.charge, mu.pt, mu.eta, mu.phi) )
                pt_err.append( mu.pt *
                                mk_safe(roccor.kScaleDTerror, mu.charge, mu.pt, mu.eta, mu.phi) )

        self.out.fillBranch("isData", isData)
        self.out.fillBranch("isMuChannel", isMuon)
        self.out.fillBranch("lep1_corrected_pt", pt_corr[0])
        self.out.fillBranch("lep2_corrected_pt", pt_corr[1])
        self.out.fillBranch("lep1_correctedUp_pt",   pt_corr[0] + pt_err[0])
        self.out.fillBranch("lep1_correctedDown_pt", pt_corr[0] - pt_err[0])
        self.out.fillBranch("lep2_correctedUp_pt",   pt_corr[1] + pt_err[1])
        self.out.fillBranch("lep2_correctedDown_pt", pt_corr[1] - pt_err[1])

        return True


#muonScaleRes2016 = lambda: muonScaleResProducer('roccor.Run2.v3',
#                                                'RoccoR2016.txt', 2016)
#muonScaleRes2017 = lambda: muonScaleResProducer('roccor.Run2.v3',
#                                                'RoccoR2017.txt', 2017)
#muonScaleRes2018 = lambda: muonScaleResProducer('roccor.Run2.v3',
#                                                'RoccoR2018.txt', 2018)
muonScaleRes = lambda: muonScaleResProducer('rocco', 'RoccoR2017.txt')
