import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR,matchObjectCollection

import numpy as np

class HCALSegmentation(object):
    """
    document the hcal segment information
    numbers taken from http://cds.cern.ch/record/896897/files/note05_016.pdf
    """
    def __init__(self):
        self.etaPosSeg = [
                            0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 
                            0.696, 0.783, 0.870, 0.957, 1.044, 1.131, 1.218, 
                            1.305, 1.392, 1.479, 1.566, 1.653, 1.740, 1.830, 
                            1.930, 2.043, 2.172, 2.322, 2.500, 2.650, 3.00001
                          ]
        self.etaSeg = [ 0. ]
        for eta in self.etaPosSeg:
            self.etaSeg.insert(0, -eta)
            self.etaSeg.append(eta)

        # below this etaCut, the phi is 5 degrees per cell. 
        # above it the phi segmentation is 10 degrees
        self.etaCut = 1.740
        self.phiHB = np.pi/36.0
        self.phiHE = np.pi/18.0

    def getEtaIndex(self, eta):
        assert eta>-3.00001 and eta<3.00001, "eta %f goes out of the range"%eta
        return np.searchsorted( self.etaSeg, eta )-1

    # return the middle of the eta cell given the index
    def getEta(self, etaIndex):
        return ( self.etaSeg[etaIndex] + self.etaSeg[etaIndex+1] )/2.0
    
    # index of phi depends on the index of eta
    def getPhiIndex(self, eta, phi):
        if phi > np.pi:
            phi = phi - 2*np.pi
        elif phi< -np.pi:
            phi = phi + 2*np.pi
        assert phi >= -np.pi and phi<= np.pi, "phi %f goes out of the range"%phi
        cell = self.phiHB if abs(eta)<self.etaCut else self.phiHE
        return int((phi + np.pi )/cell)

    def getPhi(self, etaIndex, phiIndex):
        eta = self.getEta(etaIndex)
        cell = self.phiHB if abs(eta)<self.etaCut else self.phiHE 
        return ( -np.pi +cell * (phiIndex+0.5) )

    def getEtaPhiIndex(self, eta, phi):
        return self.getEtaIndex(eta), self.getPhiIndex(eta, phi)

    def getEtaPhi(self, etaIndex, phiIndex):
        return self.getEta(etaIndex), self.getPhi(etaIndex, phiIndex)


class neuToNVPProducer(Module):
    """ 
    convert the neutral particles ( neutral hadrons and photons ) into pixels.
    Each pixel will be treated as a 'virtual' particle in the training later
    """
    def __init__(self):
        self.segment = HCALSegmentation()
    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for itype in ["PF", "Gen"]:
            self.out.branch("nNVP_%s"%itype, "I")
            self.out.branch("NVP_%s_pt"%itype,  "F", lenVar="nNVP_%s"%itype)
            self.out.branch("NVP_%s_eta"%itype, "F", lenVar="nNVP_%s"%itype)
            self.out.branch("NVP_%s_phi"%itype, "F", lenVar="nNVP_%s"%itype)
            self.out.branch("NVP_%s_photonFraction"%itype, "F", lenVar="nNVP_%s"%itype)
        self.out.branch("NVP_PF_PVFraction", "F", lenVar="nNVP_PF")
        self.out.branch("NVP_Gen_recoFraction", "F", lenVar="nNVP_Gen")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def Neutrals2Pixels(self, NeuColl):
        # 56 cells between eta[-3.0, 3.0]
        # 72/36 cells in phi HB/HE
        # save the total energy of different cells
        neuPixels = np.zeros((56, 72))
        phoPixels = np.zeros((56, 72))
        for p in NeuColl:
            ieta, iphi = self.segment.getEtaPhiIndex( p.eta, p.phi )
            ptot = p.pt * np.cosh(p.eta)
            neuPixels[ieta][iphi] += np.hypot( ptot, p.mass )
            phoPixels[ieta][iphi] += ptot * (abs(p.pdgId)==22) 
        
        phoFracPixels = phoPixels / (neuPixels + 1e-8)

        return neuPixels, phoFracPixels


    def Pixels2NVPs(self, neuPixels, phoFracPixels):
        """ 
        convert the energy deposit in cells to neutral 'virtual particles' 
        """
        NVP_pt       = []
        NVP_eta      = []
        NVP_phi      = []
        NVP_photonFraction  = []
        for ieta, iphi in zip(*np.where(neuPixels>0)):
            eta, phi = self.segment.getEtaPhi( ieta, iphi )
            pt = neuPixels[ieta, iphi] / np.cosh(eta)
            phofrac = phoFracPixels[ieta, iphi]

            NVP_pt .append( pt  )
            NVP_eta.append( eta )
            NVP_phi.append( phi )
            NVP_photonFraction.append( phofrac ) 

        return NVP_pt, NVP_eta, NVP_phi, NVP_photonFraction

    def mapPFGen(self, neuPixels_PF, neuPixels_Gen):
        """
        calculate the PV pt fractions comparing PF with Gen
        """
        NVP_PF_PVFraction = []
        NVP_Gen_recoFraction = []

        pvFractions = neuPixels_Gen / ( neuPixels_PF + 1e-8)
        for ieta, iphi in zip(*np.where(neuPixels_PF>0)):
            NVP_PF_PVFraction.append( pvFractions[ieta, iphi] )
        for ieta, iphi in zip(*np.where(neuPixels_Gen>0)):
            NVP_Gen_recoFraction.append( 1.0/pvFractions[ieta, iphi] )

        return NVP_PF_PVFraction, NVP_Gen_recoFraction



    def vetoCandidate(self, pcand, vetoCands):
        """ given input pcand and vetoCands, decide if the pcand
            should be vetoed or not. Currently only apply DeltaR cut """

        for vc in vetoCands:
            if deltaR(pcand.eta, pcand.phi, vc.eta, vc.phi) < 0.01:
                return True
        return False

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""        
        muons = Collection(event, "Muon")
        pfCands  = Collection(event, "PF")
        packedGenParts = Collection(event, "packedGenPart") # gen particles from packedGenParticles (status==1)

        neuPFCands = []
        # check if the event has a good muon
        hasGoodMuon = ( muons[0].pt>20.0 and abs(muons[0].eta)<2.4 and muons[0].tightId )
        vetoCands = []
        if hasGoodMuon:
            vetoCands.append( muons[0] )
        for p in pfCands:
            if abs(p.eta)>3.0: continue
            if p.charge!=0 and abs(p.eta)<2.5: continue
            if self.vetoCandidate(p, vetoCands ): continue
            neuPFCands.append( p )

        neuPixels_PF, phoFracPixels_PF = self.Neutrals2Pixels( neuPFCands )
        NVP_PF_pt, NVP_PF_eta, NVP_PF_phi, NVP_PF_photonFraction = self.Pixels2NVPs( neuPixels_PF, phoFracPixels_PF )

        self.out.fillBranch("nNVP_PF",    len(NVP_PF_pt))
        self.out.fillBranch("NVP_PF_pt",  NVP_PF_pt )
        self.out.fillBranch("NVP_PF_eta", NVP_PF_eta)
        self.out.fillBranch("NVP_PF_phi", NVP_PF_phi)
        self.out.fillBranch("NVP_PF_photonFraction", NVP_PF_photonFraction)

        vetoGens = []
        for gp in packedGenParts:
            if abs(gp.pdgId) in [11, 13] and gp.pt > 20.0 and abs(gp.eta)<2.4 :
                vetoGens.append( gp )
                break

        neuGenCands = []
        for gp in packedGenParts:
            if abs(gp.eta)>3.0: continue
            if gp.charge!=0 and abs(gp.eta)<2.5: continue
            if abs(gp.pdgId) in [12, 14, 16]: continue # veto neutrino
            if abs(gp.pdgId)==22 and self.vetoCandidate(gp, vetoGens): continue
            neuGenCands.append( gp )

        neuPixels_Gen, phoFracPixels_Gen = self.Neutrals2Pixels( neuGenCands )
        NVP_Gen_pt, NVP_Gen_eta, NVP_Gen_phi, NVP_Gen_photonFraction = self.Pixels2NVPs( neuPixels_Gen, phoFracPixels_Gen )

        self.out.fillBranch("nNVP_Gen",    len(NVP_Gen_pt))
        self.out.fillBranch("NVP_Gen_pt",  NVP_Gen_pt )
        self.out.fillBranch("NVP_Gen_eta", NVP_Gen_eta)
        self.out.fillBranch("NVP_Gen_phi", NVP_Gen_phi)
        self.out.fillBranch("NVP_Gen_photonFraction", NVP_Gen_photonFraction)

        # build the map between PF and Gen
        NVP_PF_PVFraction, NVP_Gen_recoFraction = self.mapPFGen( neuPixels_PF, neuPixels_Gen )
        self.out.fillBranch("NVP_PF_PVFraction",    NVP_PF_PVFraction   )
        self.out.fillBranch("NVP_Gen_recoFraction", NVP_Gen_recoFraction )

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

neuToNVPModuleConstr = lambda : neuToNVPProducer()
