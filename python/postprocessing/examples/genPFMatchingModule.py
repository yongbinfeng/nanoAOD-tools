import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR,matchObjectCollection 

class genPFMatchingProducer(Module):
    """ match final state Gen Particles (from PackedGenParticles) with reco particles (from PF candidates) """
    def __init__(self):
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        self.out = wrappedOutputTree
        
        self.out.branch("packedGenPart_toPF",      "O", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_toPFIndex", "I", lenVar="npackedGenPart")

        self.out.branch("PF_toGen",                "O", lenVar="nPF")
        self.out.branch("PF_toGenIndex",           "I", lenVar="nPF") # matched Gen index
        self.out.branch("PF_toGenrpt",             "F", lenVar="nPF") # pT_reco/ pT_matched_Gen
        self.out.branch("PF_toGendR",              "F", lenVar="nPF") # delta R between reco and matched Gen


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        pfCands  = Collection(event, "PF")
        packedGenParts = Collection(event, "packedGenPart") # gen particles from packedGenParticles (status==1)

        pfCands_selected = []
        index = 0
        for p in pfCands:
            p.toGen = False
            p.toGenIndex = -999
            p.toGenrpt = 0.
            p.toGendR = -999
            p.index = index
            if p.pt > 0.4 and abs( p.eta ) < 3.1:
                pfCands_selected.append( p )
            index += 1

        packedGenParts_selected = []
        index = 0 
        for gp in packedGenParts:
            gp.toPF = False
            gp.toPFIndex = -999
            gp.index = index
            if gp.pt > 0.4 and abs( gp.eta ) < 3.1:
                packedGenParts_selected.append( gp ) 
            index += 1

        pairs = matchObjectCollection( packedGenParts_selected, pfCands_selected, dRmax=0.1 )

        for gp, p in pairs.iteritems():
            if p is not None:
            # found matched particles
                gp.toPF = True
                gp.toPFIndex = p.index

                p.toGen = True
                p.toGenIndex = gp.index
                p.toGenrpt = p.pt / (gp.pt + 1e-6)
                p.toGendR = deltaR( gp.eta, gp.phi, p.eta, p.phi )

        # prepare the list and write to output 
        outputs = {}
        outputs["packedGenPart_toPF"]      = [ gp.toPF      for gp in packedGenParts ]
        outputs["packedGenPart_toPFIndex"] = [ gp.toPFIndex for gp in packedGenParts ]        
        outputs["PF_toGen"]                = [ p.toGen      for p  in pfCands        ]
        outputs["PF_toGenIndex"]           = [ p.toGenIndex for p  in pfCands        ]        
        outputs["PF_toGenrpt"]             = [ p.toGenrpt   for p  in pfCands        ]
        outputs["PF_toGendR"]              = [ p.toGendR    for p  in pfCands        ]

        for ikey, ioutput in outputs.iteritems():
            self.out.fillBranch( ikey,  ioutput )

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

genPFMatchingModuleConstr = lambda : genPFMatchingProducer()
 
