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

        # match gen to PF
        self.out.branch("packedGenPart_toPFIndex", "I", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_toPFdR",    "F", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_toPFrpt",   "F", lenVar="npackedGenPart")
        # gen isolation ( please note there is a pt > 0.4 cut for the gen particles. So the isolation may not be correct.)
        self.out.branch("packedGenPart_chgIso3",   "F", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_neuIso3",   "F", lenVar="npackedGenPart")        
        self.out.branch("packedGenPart_phoIso3",   "F", lenVar="npackedGenPart")

        #self.out.branch("PF_toGen",                "O", lenVar="nPF")
        #self.out.branch("PF_toGenIndex",           "I", lenVar="nPF") # closest Gen index
        #self.out.branch("PF_toGenrpt",             "F", lenVar="nPF") # pT_reco/ pT_closest_Gen
        #self.out.branch("PF_toGendR",              "F", lenVar="nPF") # delta R between reco and closest Gen
        #self.out.branch("PF_neuIso1",              "O", lenVar="nPF")
        #self.out.branch("PF_neuIso3",              "O", lenVar="nPF")
        #self.out.branch("PF_neuIso4",              "O", lenVar="nPF")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def getIsolations(self, particleColl):
        """ 
        calculate isolatios for particles in particleColl
        """
        for ip1 in xrange(len(particleColl)):
            vchg3 = ROOT.TLorentzVector()
            vneu3 = ROOT.TLorentzVector()
            vpho3 = ROOT.TLorentzVector()
            for ip2 in xrange(len(particleColl)):
                if ip2 == ip1:
                    continue
                dR = deltaR( particleColl[ip1].eta, particleColl[ip1].phi, particleColl[ip2].eta, particleColl[ip2].phi )
                if dR < 0.3:
                    p2 = particleColl[ip2]
                    plorenz = ROOT.TLorentzVector()
                    plorenz.SetPtEtaPhiM( p2.pt, p2.eta, p2.phi, p2.mass )
                    
                    vchg3 += plorenz * ( p2.charge!=0 )
                    vneu3 += plorenz * ( p2.charge==0 ) * ( abs(p2.pdgId) not in [12, 14, 16, 22] )
                    vpho3 += plorenz * ( abs(p2.pdgId)==22 )

            p1 = particleColl[ip1]
            p1.chgIso3 = vchg3.Pt()
            p1.neuIso3 = vneu3.Pt()
            p1.phoIso3 = vpho3.Pt()
            

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        pfCands  = Collection(event, "PF")
        packedGenParts = Collection(event, "packedGenPart") # gen particles from packedGenParticles (status==1)

        pfCands_selected = []
        index = 0
        # match PF candidates and gen with pt > 0.5 and |eta|<3.0
        # cut at 0.4 and 3.1 to leave some space for the resolution
        for p in pfCands:
            #p.toGenIndex = -999
            #p.toGenrpt = 0.
            #p.toGendR = -999
            #p.neuIso1 = True
            #p.neuIso3 = True
            #p.neuIso4 = True
            p.index = index
            if p.pt > 0.4 and abs( p.eta ) < 3.1:
                pfCands_selected.append( p )
            index += 1

        packedGenParts_selected = []
        index = 0 
        for gp in packedGenParts:
            gp.toPFIndex = -999
            gp.toPFrpt = -999.0
            gp.toPFdR  = -999.0
            gp.chgIso3 = -999.0
            gp.neuIso3 = -999.0
            gp.phoIso3 = -999.0
            gp.index = index
            if gp.pt > 0.4 and abs( gp.eta ) < 3.1 and abs( gp.pdgId ) not in [12, 14, 16]:
                packedGenParts_selected.append( gp ) 
            index += 1

        # calculate isolations for selected packedGenParts
        self.getIsolations( packedGenParts_selected )

        # match gen with the closest PF candidates
        pairs = matchObjectCollection( packedGenParts_selected, pfCands_selected, dRmax=1000. )

        for gp, p in pairs.iteritems():
            if p is not None:
            # found closest particles
                gp.toPFIndex = p.index
                gp.toPFdR = deltaR( gp.eta, gp.phi, p.eta, p.phi )
                gp.toPFrpt = p.pt / ( gp.pt + 1e-6 )

                #p.toGenIndex = gp.index
                #p.toGenrpt = p.pt / (gp.pt + 1e-6)
                #p.toGendR = gp.toPFdR 

        # prepare the list and write to output 
        outputs = {}
        outputs["packedGenPart_toPFIndex"] = [ gp.toPFIndex for gp in packedGenParts ]        
        outputs["packedGenPart_toPFdR"]    = [ gp.toPFdR    for gp in packedGenParts ]
        outputs["packedGenPart_toPFrpt"]   = [ gp.toPFrpt   for gp in packedGenParts ]
        outputs["packedGenPart_chgIso3"]   = [ gp.chgIso3   for gp in packedGenParts ]
        outputs["packedGenPart_neuIso3"]   = [ gp.neuIso3   for gp in packedGenParts ]
        outputs["packedGenPart_phoIso3"]   = [ gp.phoIso3   for gp in packedGenParts ]

        #outputs["PF_toGenIndex"]           = [ p.toGenIndex for p  in pfCands        ]        
        #outputs["PF_toGenrpt"]             = [ p.toGenrpt   for p  in pfCands        ]
        #outputs["PF_toGendR"]              = [ p.toGendR    for p  in pfCands        ]
        #outputs["PF_neuIso1"]              = [ p.neuIso1    for p  in pfCands        ]
        #outputs["PF_neuIso3"]              = [ p.neuIso3    for p  in pfCands        ]
        #outputs["PF_neuIso4"]              = [ p.neuIso4    for p  in pfCands        ]

        for ikey, ioutput in outputs.iteritems():
            self.out.fillBranch( ikey,  ioutput )

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

genPFMatchingModuleConstr = lambda : genPFMatchingProducer()
 
