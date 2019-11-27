import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class isDataProducer(Module):
    """ decide if the input file is data or MC"""
    def __init__(self):
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("isData",  "O");
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        if event.run > 1: 
            isData = True
        else:
            isData = False
        print event.run, "isData", isData
        self.out.fillBranch("isData",isData )
        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

isDataModuleConstr = lambda : isDataProducer() 
 
