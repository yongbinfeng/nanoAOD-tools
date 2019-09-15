import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR,matchObjectCollection

import numpy as np
import time

class pfGenMatchingProducer(Module):
    """
    chg refers to charge particles (not just charged hadrons)
    neu refers to neutral hadrons
    pho refers to photons
    """
    def __init__(self):
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        self.out = wrappedOutputTree

        # match PF to Gen
        self.out.branch("packedGenPart_toPFmindRIndex",   "I", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_toPFmindR",        "F", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_toPFN",            "I", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_ptype",            "I", lenVar="npackedGenPart")

        self.out.branch("PF_toGendR",              "F", lenVar="nPF") # smallest delta R between reco and Gen 
        self.out.branch("PF_toGenIndex",           "I", lenVar="nPF") # index of the Gen with mindR to PF
        self.out.branch("PF_ptype",                "I", lenVar="nPF")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def timing(f):
        def wrap(*args):
            time1 = time.time()
            ret = f(*args)
            time2 = time.time()
            print '%s function took %0.4f s' % (f.func_name, (time2-time1))
            return ret
        return wrap 

    def getGenType(self, gp):
        """ return the type of gen particles """
        ptype = 0
        if gp.charge!=0 and abs(gp.eta)<3.0:
            # charged particle
            # there are still some charged pf outside tracker (2.5). 
            # so try matching all charged gen within 3.0 first.
            # If not matched, change the gen type from charged to neutral
            ptype = 1
        elif gp.charge==0 and abs(gp.pdgId) not in [12, 14, 16, 22]:
            # neutral hadrons
            ptype = 2
        elif gp.charge==0 and abs(gp.pdgId)==22:
            # photon
            ptype= 3
        return ptype

    def getPFType(self, p):
        """ return the type for pf particles """
        ptype = 0
        if p.charge!=0 and abs(p.eta)<3.0:
            ptype = 1
        elif p.charge==0 and abs(p.pdgId) == 130:
            ptype = 2
        elif p.charge==0 and abs(p.pdgId) == 22:
            ptype = 3
        return ptype


    def mapPFToGen(self, GenColl, PFColl, selGen=lambda x: True, selPF=lambda x:True):
        # map charged
        for p in PFColl:
            if not selPF(p): continue
            dRmin = p.toGendR
            matchedIndex = -999

            for igp in xrange(len(GenColl)):
                gp = GenColl[igp]
                if not selGen(gp): continue
                dR = deltaR(gp.eta, gp.phi, p.eta, p.phi)
                if dR < dRmin:
                    dRmin = dR
                    matchedIndex = igp

            if matchedIndex>=0:
                gp = GenColl[matchedIndex]
                p.toGenIndex = gp.index
                p.toGendR = dRmin
                gp.toPFN += 1
                if dRmin < gp.toPF_mindR:
                    gp.toPF_mindR = dRmin
                    gp.toPF_mindRIndex = p.index


    def writeToOutput(self, packedGenParts, pfCands):
        # prepare the list and write to output
        outputs = {}
        outputs["packedGenPart_toPFmindRIndex"]   = [ gp.toPF_mindRIndex   for gp in packedGenParts ]
        outputs["packedGenPart_toPFmindR"]        = [ gp.toPF_mindR        for gp in packedGenParts ]
        outputs["packedGenPart_toPFN"]            = [ gp.toPFN             for gp in packedGenParts ]
        outputs["packedGenPart_ptype"]            = [ gp.ptype             for gp in packedGenParts ]

        outputs["PF_toGendR"]          = [ p.toGendR                      for p in pfCands ]
        outputs["PF_toGenIndex"]       = [ p.toGenIndex                   for p in pfCands ]
        outputs["PF_ptype"]            = [ p.ptype                        for p in pfCands ]

        for ikey, ioutput in outputs.iteritems():
            self.out.fillBranch( ikey,  ioutput )

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        pfCands  = Collection(event, "PF")
        packedGenParts = Collection(event, "packedGenPart") # gen particles from packedGenParticles (status==1)

        pfCands_selected = []
        index = 0
        # match PF candidates and gen with |eta|<3.0
        for p in pfCands:
            p.toGendR = 999.0
            p.toGenIndex = -999
            p.index = index
            index += 1
            p.ptype = 0
            if abs( p.eta ) > 3.0: continue
            p.ptype = self.getPFType( p )
            pfCands_selected.append( p )

        packedGenParts_selected = []
        index = 0 
        for gp in packedGenParts:
            gp.toPF_mindRIndex = -999
            gp.toPF_mindR  = 999.0
            gp.toPFN = 0
            gp.index = index
            index += 1
            gp.ptype = 0
            if abs(gp.pdgId) in [12, 14, 16]: continue
            if abs(gp.eta)>3.0: continue
            gp.ptype = self.getGenType( gp ) 
            packedGenParts_selected.append( gp ) 

        # map charged PF to charged Gen
        self.mapPFToGen( packedGenParts_selected, pfCands_selected, selGen=lambda gp: gp.ptype==1, selPF=lambda p: p.ptype==1 )
        
        # for gen chg, if outside 2.5 and not matched to PF, treat it as neutral hadrons or photons
        for gp in packedGenParts_selected:
            if gp.ptype==1 and abs(gp.eta)>2.5 and gp.toPF_mindR>0.04:
                # electrons outside the tracker will be photons
                # otherwise neutral hadrons
                gp.ptype = 3 if abs(gp.pdgId)==11 else 2
        
        # map PF photons to Gen photons
        self.mapPFToGen( packedGenParts_selected, pfCands_selected, selGen=lambda gp: gp.ptype==3, selPF=lambda p: p.ptype==3 )

        # map PF neutral hadrons and PF photons unmatched to Gen photons and electrons to neutral hadrons
        self.mapPFToGen( packedGenParts_selected, pfCands_selected, selGen=lambda gp: gp.ptype==2, selPF=lambda p: (p.ptype==2 or (p.ptype==3 and p.toGendR>0.04)) )

        #  write the new properties of packedGenParts and PF cands to output 
        self.writeToOutput( packedGenParts, pfCands )

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

pfGenMatchingModuleConstr = lambda : pfGenMatchingProducer()
 
