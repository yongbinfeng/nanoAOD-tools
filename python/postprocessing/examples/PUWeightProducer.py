import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class PUWeightProducer(Module):
    """ decide if the input file is data or MC"""
    def __init__(self):
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("PUWeight",  "F");
        fdata = ROOT.TFile("/afs/cern.ch/work/y/yofeng/public/PFCand/CMSSW_10_2_10/src/PhysicsTools/NanoAODTools/python/postprocessing/data/pileup/mine/MyDataPileupHistogram.root")
        self.pudata = fdata.Get("pileup").Clone()
        self.pudata.SetDirectory(0)
        fmc = ROOT.TFile("/afs/cern.ch/work/y/yofeng/public/PFCand/CMSSW_10_2_10/src/PhysicsTools/NanoAODTools/python/postprocessing/data/pileup/mine/pileupWeight_DYJets.root")
        self.pumc = fmc.Get("pile_truthIn").Clone()
        self.pumc.SetDirectory(0)
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        if event.isData:
            PUWeight = 1.
        else:
            npu = event.Pileup_nTrueInt
            ntotal_data = self.pudata.Integral()
            ntotal_mc = self.pumc.Integral()
            den = self.pumc.GetBinContent( self.pumc.FindBin(npu) ) / ntotal_mc
            num = self.pudata.GetBinContent( self.pudata.FindBin(npu) ) / ntotal_data
            PUWeight = ( num / den )
            #print npu, num, den, PUWeight
        self.out.fillBranch("PUWeight",PUWeight )
        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

PUWeightModuleConstr = lambda : PUWeightProducer() 
 
