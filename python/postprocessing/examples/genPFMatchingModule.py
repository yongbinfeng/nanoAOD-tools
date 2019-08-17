import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR,matchObjectCollection

import numpy as np
import time

class genPFMatchingProducer(Module):
    """ 
    match final state Gen Particles (from PackedGenParticles) with reco particles (from PF candidates).
    Currently using deltaR matching: for one gp, loop over pf cands and find the closest p.
    If more than one gp is matched to the p, the rpt and dR for p is calculated based on p and sum(gp_matched)

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

        # match gen to PF
        self.out.branch("packedGenPart_toPFIndex",   "I", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_toPFdR",      "F", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_toPFrpt",     "F", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_toSTPFIndex", "I", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_toSTPFdR",    "F", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_toSTPFrpt",   "F", lenVar="npackedGenPart")
        # gen isolation ( please note there is a pt > 0.4 cut for the gen particles. So the isolation may not be correct.)
        self.out.branch("packedGenPart_chgIso1",     "F", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_neuIso1",     "F", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_phoIso1",     "F", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_chgPFIso1",   "F", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_neuPFIso1",   "F", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_phoPFIso1",   "F", lenVar="npackedGenPart")
        self.out.branch("packedGenPart_ptype",       "I", lenVar="npackedGenPart")

        self.out.branch("PF_toGenN",               "I", lenVar="nPF") # number of gen particles closest to this PF
        self.out.branch("PF_toGenrpt",             "F", lenVar="nPF") # pT_reco/ pT_sum_Gen
        self.out.branch("PF_toGendR",              "F", lenVar="nPF") # delta R between reco and the direction of sum Gen
        self.out.branch("PF_toGenmindR",           "F", lenVar="nPF") # min delta R between reco PF and the gens that are closest to this PF
        self.out.branch("PF_toGenmindR_rpt",       "F", lenVar="nPF") # pT_reco / pT_Gen of the above Gen
        self.out.branch("PF_toSTGenN",             "I", lenVar="nPF") # number of gen particles closest to this PF
        self.out.branch("PF_toSTGenrpt",           "F", lenVar="nPF") # pT_reco/ pT_sum_Gen
        self.out.branch("PF_toSTGendR",            "F", lenVar="nPF") # delta R between reco and the direction of sum Gen
        self.out.branch("PF_toSTGenmindR",         "F", lenVar="nPF") # min delta R between reco PF and the gens that are closest to this PF
        self.out.branch("PF_toSTGenmindR_rpt",     "F", lenVar="nPF") # pT_reco / pT_Gen of the above Gen
        self.out.branch("PF_ptype",                "I", lenVar="nPF")
        #self.out.branch("PF_neuIso1",              "O", lenVar="nPF")
        #self.out.branch("PF_neuIso3",              "O", lenVar="nPF")
        #self.out.branch("PF_neuIso4",              "O", lenVar="nPF")


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
        if gp.charge!=0 and abs(gp.eta)<2.5:
            # charged particle
            ptype = 1
        elif ( ( gp.charge==0 and abs(gp.pdgId) not in [12, 14, 16, 22] ) or ( gp.charge!=0 and abs(gp.eta)>2.5 ) ):
            # neutral hadrons
            ptype = 2
        if gp.charge==0 and abs(gp.pdgId)==22:
            # photon
            ptype= 3
        return ptype

    def getPFType(self, p):
        """ return the type for pf particles """
        ptype = 0
        if p.charge!=0 and abs(p.eta)<2.5:
            ptype = 1
        if p.charge==0 and abs(p.pdgId) == 130:
            ptype = 2
        if p.charge==0 and abs(p.pdgId) in [1,2,22]:# treat HF hadron and em as pho fraction for now
            ptype = 3
        return ptype

    def getIsolations(self, particleColl):
        """ 
        calculate isolatios for particles in particleColl
        """
        for ip1 in xrange(len(particleColl)):
            vchg1 = ROOT.TLorentzVector()
            vneu1 = ROOT.TLorentzVector()
            vpho1 = ROOT.TLorentzVector()
            for ip2 in xrange(len(particleColl)):
                if ip2 == ip1:
                    continue
                dR = deltaR( particleColl[ip1].eta, particleColl[ip1].phi, particleColl[ip2].eta, particleColl[ip2].phi )
                if dR < 0.1:
                    p2 = particleColl[ip2]
                    plorenz = p2.p4()
                    vchg1 += plorenz * ( p2.ptype==1 )
                    vneu1 += plorenz * ( p2.ptype==2 )
                    vpho1 += plorenz * ( p2.ptype==3 )

            p1 = particleColl[ip1]
            p1.chgIso1 = vchg1.Pt()
            p1.neuIso1 = vneu1.Pt()
            p1.phoIso1 = vpho1.Pt()

    def getRecoIsolationForGen(self, GenColl, PFColl):
        """ calculate isolation using pf candiates for gen particles """
        for igp in xrange(len(GenColl)):
            vchg1 = ROOT.TLorentzVector()
            vneu1 = ROOT.TLorentzVector()
            vpho1 = ROOT.TLorentzVector()

            for ip in np.argwhere(self.AdrGenPF[igp,:]<0.1).flatten():
                p = PFColl[ip]
                plorenz = p.p4()
                vchg1 += plorenz * ( p.ptype==1 )
                vneu1 += plorenz * ( p.ptype==2 )
                vpho1 += plorenz * ( p.ptype==3 )

            gp = GenColl[igp]
            gp.chgPFIso1 = vchg1.Pt()
            gp.neuPFIso1 = vneu1.Pt()
            gp.phoPFIso1 = vpho1.Pt()

    def makeDRMapGenPF(self, GenColl, PFColl):
        """ save the deltaR values betwee Gen and PF into numpy arrays """

        self.AdrGenPF   = np.full( (len(GenColl), len(PFColl)), np.inf )
        # map from Gen to PF of the same type (chg to chg, neu to neu, and pho to pho)
        self.AdrGenPFST = np.full( (len(GenColl), len(PFColl)), np.inf )

        for igp in xrange(len(GenColl)):
            for ip in xrange(len(PFColl)):
                dR = deltaR( GenColl[igp].eta, GenColl[igp].phi, PFColl[ip].eta, PFColl[ip].phi )
                ## only save dR<0.1 values now, to save time in changing np array values
                #if dR<0.1:
                self.AdrGenPF[igp, ip] = dR
                if GenColl[igp].ptype == PFColl[ip].ptype:
                    self.AdrGenPFST[igp, ip] = dR


    def mapGenPF(self, GenColl, PFColl):
        """ map gen particles to PF candidates """
        # inclusive mapping: loop over all gen particles and pf candidates in the GenColl and PFColl collection
        ApfIndices   = np.argmin(self.AdrGenPF, axis=1)
        ApfIndicesST = np.argmin(self.AdrGenPFST, axis=1)

        for igp, ip in enumerate(ApfIndices):
            gp = GenColl[igp]
            p  = PFColl[ip]
            dR = self.AdrGenPF[igp, ip]

            gp.toPFIndex = p.index
            gp.toPFdR = dR            
            gp.toPFrpt = p.pt / ( gp.pt + 1e-6 )

            p.toGenN += 1
            p.toGen_sumpt += gp.p4()
            if dR< p.toGen_mindR:
                p.toGen_mindR = dR
                p.toGen_mindR_rpt = p.pt / ( gp.pt + 1e-6 )

        for igp, ip in enumerate(ApfIndicesST):
            gp = GenColl[igp]
            p  = PFColl[ip]
            dR = self.AdrGenPF[igp, ip]

            gp.toSTPFIndex = p.index
            gp.toSTPFdR = dR
            gp.toSTPFrpt = p.pt / ( gp.pt + 1e-6 )

            p.toSTGenN += 1
            p.toSTGen_sumpt += gp.p4()
            if dR< p.toSTGen_mindR:
                p.toSTGen_mindR = dR
                p.toSTGen_mindR_rpt = p.pt / ( gp.pt + 1e-6 )

    def writeToOutput(self, packedGenParts, pfCands):
        # prepare the list and write to output
        outputs = {}
        outputs["packedGenPart_toPFIndex"]   = [ gp.toPFIndex   for gp in packedGenParts ]
        outputs["packedGenPart_toPFdR"]      = [ gp.toPFdR      for gp in packedGenParts ]
        outputs["packedGenPart_toPFrpt"]     = [ gp.toPFrpt     for gp in packedGenParts ]
        outputs["packedGenPart_toSTPFIndex"] = [ gp.toSTPFIndex for gp in packedGenParts ]
        outputs["packedGenPart_toSTPFdR"]    = [ gp.toSTPFdR    for gp in packedGenParts ]
        outputs["packedGenPart_toSTPFrpt"]   = [ gp.toSTPFrpt   for gp in packedGenParts ]
        outputs["packedGenPart_chgIso1"]     = [ gp.chgIso1     for gp in packedGenParts ]
        outputs["packedGenPart_neuIso1"]     = [ gp.neuIso1     for gp in packedGenParts ]
        outputs["packedGenPart_phoIso1"]     = [ gp.phoIso1     for gp in packedGenParts ]
        outputs["packedGenPart_chgPFIso1"]   = [ gp.chgPFIso1   for gp in packedGenParts ]
        outputs["packedGenPart_neuPFIso1"]   = [ gp.neuPFIso1   for gp in packedGenParts ]
        outputs["packedGenPart_phoPFIso1"]   = [ gp.phoPFIso1   for gp in packedGenParts ]
        outputs["packedGenPart_ptype"]       = [ gp.ptype       for gp in packedGenParts ]

        outputs["PF_toGenN"]           = [ p.toGenN                           for p in pfCands ]
        outputs["PF_toGenmindR"]       = [ p.toGen_mindR                      for p in pfCands ]
        outputs["PF_toGenmindR_rpt"]   = [ p.toGen_mindR_rpt                  for p in pfCands ]
        outputs["PF_toGenrpt"]         = [ p.pt/(p.toGen_sumpt.Pt()+1e-6)     for p in pfCands ]
        outputs["PF_toGendR"]          = [ p.p4().DeltaR(p.toGen_sumpt)       for p in pfCands ]
        outputs["PF_toSTGenN"]         = [ p.toSTGenN                         for p in pfCands ]
        outputs["PF_toSTGenmindR"]     = [ p.toSTGen_mindR                    for p in pfCands ]
        outputs["PF_toSTGenmindR_rpt"] = [ p.toSTGen_mindR_rpt                for p in pfCands ]
        outputs["PF_toSTGenrpt"]       = [ p.pt/(p.toSTGen_sumpt.Pt()+1e-6)   for p in pfCands ]
        outputs["PF_toSTGendR"]        = [ p.p4().DeltaR(p.toSTGen_sumpt)     for p in pfCands ]
        outputs["PF_ptype"]            = [ p.ptype                            for p in pfCands ]
        #outputs["PF_neuIso1"]              = [ p.neuIso1    for p  in pfCands        ]
        #outputs["PF_neuIso3"]              = [ p.neuIso3    for p  in pfCands        ]
        #outputs["PF_neuIso4"]              = [ p.neuIso4    for p  in pfCands        ]

        for ikey, ioutput in outputs.iteritems():
            self.out.fillBranch( ikey,  ioutput )

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        pfCands  = Collection(event, "PF")
        packedGenParts = Collection(event, "packedGenPart") # gen particles from packedGenParticles (status==1)

        pfCands_selected = []
        index = 0
        # match PF candidates and gen with pt > 0.5 and |eta|<3.0
        # cut at 0.4 and 3.1 to leave some space for the resolution
        for p in pfCands:
            p.toGenN = 0
            p.toGen_sumpt = ROOT.TLorentzVector()
            p.toGen_mindR = 999.0
            p.toGen_mindR_rpt = -999.0
            p.toSTGenN = 0
            p.toSTGen_sumpt = ROOT.TLorentzVector()
            p.toSTGen_mindR = 999.0
            p.toSTGen_mindR_rpt = -999.0
            p.index = index
            p.ptype = 0
            if p.pt > 0.4 and abs( p.eta ) < 3.1:
                p.ptype = self.getPFType( p )
                pfCands_selected.append( p )
            index += 1

        packedGenParts_selected = []
        index = 0 
        for gp in packedGenParts:
            gp.toPFIndex = -999
            gp.toPFrpt = -999.0
            gp.toPFdR  = -999.0
            gp.toSTPFIndex = -999
            gp.toSTPFrpt = -999.0
            gp.toSTPFdR  = -999.0
            gp.chgIso1 = -999.0
            gp.neuIso1 = -999.0
            gp.phoIso1 = -999.0
            gp.chgPFIso1 = -999.0
            gp.neuPFIso1 = -999.0
            gp.phoPFIso1 = -999.0
            gp.index = index
            gp.ptype = 0
            if gp.pt > 0.4 and abs( gp.eta ) < 3.1 and abs( gp.pdgId ) not in [12, 14, 16]:
                gp.ptype = self.getGenType( gp ) 
                packedGenParts_selected.append( gp ) 
            index += 1

        # calculate deltaR between gen and pf and save to numpy array
        self.makeDRMapGenPF( packedGenParts_selected, pfCands_selected )

        # calculate isolations for selected packedGenParts
        #self.getIsolations( packedGenParts_selected )
        #self.getRecoIsolationForGen( packedGenParts_selected, pfCands_selected )

        ## match gen with the closest PF candidates
        self.mapGenPF( packedGenParts_selected, pfCands_selected )

        #  write the new properties of packedGenParts and PF cands to output 
        self.writeToOutput( packedGenParts, pfCands )

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

genPFMatchingModuleConstr = lambda : genPFMatchingProducer()
 
