import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR

class recoilProducer(Module):
    """ produce some flags and variables for the truth and reco recoil """
    def __init__(self):
        self.recoils = ["tk", "tk_pt500", 
                        "tk_eta2p5", "tk_pt500_eta2p5", 
                        "puppi", "puppi_pt500",
                        "puppi_eta2p5",  "puppi_pt500_eta2p5", 
                        "puppi_eta3p0", "puppi_pt500_eta3p0", 
                        "puppi_eta5p0", "puppi_pt500_eta5p0"]
        self.genrecoils = ["tkGen", "tkGen_pt500", 
                           "tkGen_eta2p5", "tkGen_pt500_eta2p5", 
                           "Gen", "Gen_pt500",
                           "Gen_eta2p5", "Gen_pt500_eta2p5",
                           "Gen_eta3p0", "Gen_pt500_eta3p0",
                           "Gen_eta5p0", "Gen_pt500_eta5p0",
                           ]
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        self.out = wrappedOutputTree
        self.out.branch("foundGenW",     "O");
        self.out.branch("pT_W",          "F");
        self.out.branch("eta_W",         "F");
        self.out.branch("phi_W",         "F");

        self.out.branch("hasGoodMuon",   "O");

        for u in self.recoils + self.genrecoils:
            self.out.branch("u_%s_pt"  %u,    "F");
            self.out.branch("u_%s_eta" %u,    "F");
            self.out.branch("u_%s_phi" %u,    "F");
            self.out.branch("u_%s_mass"%u,    "F");


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def getGenW(self, genParts):
        """ to find the Gen W boson and its momentum from genParts """
        foundGenW = False
        pT_W  = -999.0
        phi_W = -999.0
        eta_W = -999.0

        for p in genParts:
            if abs(p.pdgId) != 24: continue
            if p.status != 62: continue
            foundGenW = True
            pT_W = p.pt
            eta_W = p.eta
            phi_W = p.phi
            break

        return (foundGenW, pT_W, eta_W, phi_W)

    def vetoPF(self, pfcand, vetoCands):
        """ given input pfcand and vetoCands, decide if the PF cand
            should be vetoed or not. Currently only apply DeltaR cut """

        for ic in xrange(0, len(vetoCands)):
            if deltaR(pfcand.eta, pfcand.phi, ic.eta, ic.phi) < 0.01:
                return True
        return False


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        muons = Collection(event, "Muon")
        genParts = Collection(event, "GenPart")        
        pfCands  = Collection(event, "PF")
        packedGenParts = Collection(event, "packedGenPart") # gen particles from packedGenParticles (status==1)

        # get the Gen W
        foundGenW, pT_W, eta_W, phi_W = self.getGenW( genParts )
        self.out.fillBranch("foundGenW", foundGenW )
        self.out.fillBranch("pT_W",      pT_W      )
        self.out.fillBranch("eta_W",     eta_W     )
        self.out.fillBranch("phi_W",     phi_W     )

        # check if the event has a good muon
        hasGoodMuon = ( muons[0].pt<20.0 and abs(muons[0].eta)>2.4 and muons[0].Muon_tightId )
        self.out.fillBranch("hasGoodMuon",  hasGoodMuon )
        vetoCands = []
        if hasGoodMuon:
            vetoCands.append( muons[0] )

        # calculate different kinds of reco recoil
        ulorenzs = {}
        for u in self.recoils:
            ulorenzs[u] = ROOT.TLorentzVector()

        for p in pfCands:
            if self.vetoPF(p, vetoCands ): continue

            plorenz = ROOT.TLorentzVector()
            plorenz.SetPtEtaPhiM( p.pt, p.eta, p.phi, p.mass )
            isGoodTk = ( p.charge!=0 and p.fromPV >= 2 )

            ulorenzs["puppi"]              += plorenz * p.puppiWeightNoLep
            ulorenzs["puppi_pt500"]        += plorenz * p.puppiWeightNoLep * ( p.pt> 0.5 )
            ulorenzs["puppi_eta2p5"]       += plorenz * p.puppiWeightNoLep * ( abs(p.eta) < 2.5 )
            ulorenzs["puppi_pt500_eta2p5"] += plorenz * p.puppiWeightNoLep * ( p.pt> 0.5 ) * ( abs(p.eta) < 2.5 )
            ulorenzs["puppi_eta3p0"]       += plorenz * p.puppiWeightNoLep * ( abs(p.eta) < 3.0 )
            ulorenzs["puppi_pt500_eta3p0"] += plorenz * p.puppiWeightNoLep * ( p.pt> 0.5 ) * ( abs(p.eta) < 3.0 )
            ulorenzs["puppi_eta5p0"]       += plorenz * p.puppiWeightNoLep * ( abs(p.eta) < 5.0 )
            ulorenzs["puppi_pt500_eta5p0"] += plorenz * p.puppiWeightNoLep * ( p.pt> 0.5 ) * ( abs(p.eta) < 5.0 )

            ulorenzs["tk"]                 += plorenz * isGoodTk
            ulorenzs["tk_pt500"]           += plorenz * isGoodTk * ( p.pt> 0.5 )
            ulorenzs["tk_eta2p5"]          += plorenz * isGoodTk * ( abs(p.eta) < 2.5 )
            ulorenzs["tk_pt500_eta2p5"]    += plorenz * isGoodTk * ( p.pt> 0.5 ) * ( abs(p.eta) < 2.5 )

        for uname, ulorenz in ulorenzs.iteritems():
            self.out.fillBranch( "u_%s_pt"  %uname,    ulorenz.Pt()  )
            self.out.fillBranch( "u_%s_eta" %uname,    ulorenz.Eta() ) 
            self.out.fillBranch( "u_%s_phi" %uname,    ulorenz.Phi() )
            self.out.fillBranch( "u_%s_mass"%uname,    ulorenz.M()   )

        # calculate the Gen recoils
        ugenlorenzs = {}
        for ugen in self.genrecoils:
            ugenlorenzs[ugen] = ROOT.TLorentzVector()

        nlep = 0
        for p in packedGenParts:
            if abs(p.pdgId) in [12, 14, 16]: continue # veto neutrino
            if abs(p.pdgId) in [11, 13] and nlep < 1:
                # veto the first lepton
                nlep += 1
                continue

            plorenz = ROOT.TLorentzVector()
            plorenz.SetPtEtaPhiM( p.pt, p.eta, p.phi, p.mass )

            ugenlorenzs["Gen"]              += plorenz
            ugenlorenzs["Gen_pt500"]        += plorenz * ( p.pt > 0.5 )
            ugenlorenzs["Gen_eta2p5"]       += plorenz * ( abs(p.eta)<2.5 )
            ugenlorenzs["Gen_pt500_eta2p5"] += plorenz * ( p.pt > 0.5 ) * (abs(p.eta)<2.5 )
            ugenlorenzs["Gen_eta3p0"]       += plorenz * ( abs(p.eta)<3.0 )
            ugenlorenzs["Gen_pt500_eta3p0"] += plorenz * ( p.pt > 0.5 ) * (abs(p.eta)<3.0 )
            ugenlorenzs["Gen_eta5p0"]       += plorenz * ( abs(p.eta)<5.0 )
            ugenlorenzs["Gen_pt500_eta5p0"] += plorenz * ( p.pt > 0.5 ) * (abs(p.eta)<5.0 )

            ugenlorenzs["tkGen"]              += plorenz * ( p.charge!=0 )
            ugenlorenzs["tkGen_pt500"]        += plorenz * ( p.charge!=0 ) * ( p.pt > 0.5 )
            ugenlorenzs["tkGen_eta2p5"]       += plorenz * ( p.charge!=0 ) * ( abs(p.eta)<2.5 )
            ugenlorenzs["tkGen_pt500_eta2p5"] += plorenz * ( p.charge!=0 ) * ( p.pt > 0.5 ) * (abs(p.eta)<2.5 )

        for ugenname, ugenlorenz in ugenlorenzs.iteritems():
            self.out.fillBranch( "u_%s_pt"  %ugenname,    ugenlorenz.Pt()  )
            self.out.fillBranch( "u_%s_eta" %ugenname,    ugenlorenz.Eta() )
            self.out.fillBranch( "u_%s_phi" %ugenname,    ugenlorenz.Phi() )
            self.out.fillBranch( "u_%s_mass"%ugenname,    ugenlorenz.M()   )


        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

recoilModuleConstr = lambda : recoilProducer()
 
