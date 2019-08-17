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
        self.nPhiHB = 72
        self.nPhiHE = 36

    def getIEta(self, eta):
        assert eta>-3.00001 and eta<3.00001, "eta %f goes out of the range"%eta
        return np.searchsorted( self.etaSeg, eta )-1

    # return the middle of the eta cell given the index
    def getEta(self, ieta):
        return ( self.etaSeg[ieta] + self.etaSeg[ieta+1] )/2.0
    
    # index of phi depends on the index of eta
    def getIPhi(self, eta, phi):
        if phi >= np.pi:
            phi = phi - 2*np.pi
        elif phi< -np.pi:
            phi = phi + 2*np.pi
        assert phi >= -np.pi and phi< np.pi, "phi %f goes out of the range"%phi
        cell = 2*np.pi/self.nPhiHB if abs(eta)<self.etaCut else 2*np.pi/self.nPhiHE
        return int((phi + np.pi )/cell)

    def getPhi(self, ieta, iphi):
        eta = self.getEta(ieta)
        nPhi = self.nPhiHB if abs(eta)<self.etaCut else self.nPhiHE
        cell = 2*np.pi/nPhi
        return ( -np.pi +cell * (iphi+0.5) )

    def getIEtaIPhi(self, eta, phi):
        return self.getIEta(eta), self.getIPhi(eta, phi)

    def getEtaPhi(self, ieta, iphi):
        return self.getEta(ieta), self.getPhi(ieta, iphi)

    def getConnectedIEtaIPhis(self, ieta, iphi):
        # return the ieta, iphi of the cells connected the to the (ieta, iphi) cell
        ietaiphicands = [(ieta-1, iphi), (ieta+1, iphi), (ieta, iphi-1), (ieta, iphi+1)]
        ietaiphis = []
        for ieta, iphi in ietaiphicands:
            if ieta < 0 or ieta >= len(self.etaSeg)-1:
                continue

            eta = self.getEta(ieta)
            nPhi = self.nPhiHB if abs(eta)<self.etaCut else self.nPhiHE
            if iphi<0:
                iphi = iphi + nPhi
            elif iphi>=nPhi:
                iphi = iphi - nPhi
            ietaiphis.append( (ieta, iphi) )

        return ietaiphis

    def getAvgEta(self, ieta1, ieta2):
        # return the average eta of ieta1 and ieta2
        return ( self.getEta(ieta1) + self.getEta(ieta2) )/2.0

    def getAvgPhi(self, ieta1, iphi1, ieta2, iphi2):
        phi1 = self.getPhi(ieta1, iphi1)
        phi2 = self.getPhi(ieta2, iphi2)
        avg = (phi1+phi2)/2.0
        if abs(phi1-phi2) < np.pi:
            return avg
        elif avg<0:
            return (avg+np.pi)
        else:
            return (avg-np.pi)


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
            self.out.branch("NVP_%s_isMerged"%itype, "O", lenVar="nNVP_%s"%itype)
        self.out.branch("NVP_PF_toGenIndex", "I", lenVar="nNVP_PF")
        self.out.branch("NVP_Gen_toPFIndex", "I", lenVar="nNVP_Gen")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def Neutrals2Pixels(self, NeuColl):
        # 56 cells between eta[-3.0, 3.0]
        # 72/36 cells in phi HB/HE
        # save the total energy of different cells
        neuPixels = np.zeros((56, 72))
        phoPixels = np.zeros((56, 72))
        for p in NeuColl:
            ieta, iphi = self.segment.getIEtaIPhi( p.eta, p.phi )
            ptot = p.pt * np.cosh(p.eta)
            neuPixels[ieta][iphi] += np.hypot( ptot, p.mass )
            phoPixels[ieta][iphi] += ptot * (abs(p.pdgId)==22) 
        
        return neuPixels, phoPixels

    def MergePixels(self):
        """
        Merge the pixels for Gen and PF if one is PF > Gen and the other is PF < Gen
        """
        for ieta, iphi in zip(*np.where(self.neuPixels_Gen>0)):

            E_Gen = self.neuPixels_Gen[ieta, iphi]
            if E_Gen/np.cosh(self.segment.getEta(ieta)) < 1.0:
                # only merge pT > 1.0 GeV cells
                continue
            if self.neuPixels_PF[ieta, iphi] > 0.8 * E_Gen:
                continue 

            r_best = 1e8
            ietaiphim = None
            for ietac, iphic in self.segment.getConnectedIEtaIPhis(ieta, iphi):
                r_best = min(abs(E_Gen/(self.neuPixels_PF[ieta, iphi]+1e-8)-1.0), abs(self.neuPixels_Gen[ietac, iphic]/(self.neuPixels_PF[ietac, iphic]+1e-8)-1.0), r_best)
                r_merged = abs((E_Gen+self.neuPixels_Gen[ietac, iphic])/(self.neuPixels_PF[ieta, iphi]+self.neuPixels_PF[ietac,iphic]+1e-8)-1.0)
                if r_best> r_merged:
                    r_best = r_merged
                    ietaiphim = (ietac, iphic)

            if ietaiphim:
                ietam = ietaiphim[0]
                iphim = ietaiphim[1]

                E_Gen_merged        = E_Gen                          + self.neuPixels_Gen[ietam, iphim]
                E_pho_Gen_merged    = self.phoPixels_Gen[ieta, iphi] + self.phoPixels_Gen[ietam, iphim]
                E_PF_merged         = self.neuPixels_PF [ieta, iphi] + self.neuPixels_PF [ietam, iphim]
                E_pho_PF_merged     = self.phoPixels_PF [ieta, iphi] + self.phoPixels_PF [ietam, iphim]
                eta_merged = self.segment.getAvgEta(ieta, ietam)
                phi_merged = self.segment.getAvgPhi(ieta, iphi, ietam, iphim)

                self.neuPixels_Gen_merged.append((E_Gen_merged, E_pho_Gen_merged, eta_merged, phi_merged))
                self.neuPixels_PF_merged .append((E_PF_merged,  E_pho_PF_merged,  eta_merged, phi_merged))

                self.neuPixels_Gen[ieta,  iphi ]= 0
                self.neuPixels_Gen[ietam, iphim]=0
                self.neuPixels_PF[ieta,  iphi ] = 0
                self.neuPixels_PF[ietam, iphim]=0


    def Pixels2NVPs(self, neuPixels, phoPixels, neuPixels_merged = []):
        """ 
        convert the energy deposit in cells to neutral 'virtual particles' 
        """
        NVP_pt       = []
        NVP_eta      = []
        NVP_phi      = []
        NVP_photonFraction  = []
        NVP_isMerged = []
        for ieta, iphi in zip(*np.where(neuPixels>0)):
            eta, phi = self.segment.getEtaPhi( ieta, iphi )
            pt = neuPixels[ieta, iphi] / np.cosh(eta)
            phofrac = phoPixels[ieta, iphi] / (neuPixels[ieta, iphi] + 1e-8)

            NVP_pt .append( pt  )
            NVP_eta.append( eta )
            NVP_phi.append( phi )
            NVP_photonFraction.append( phofrac ) 
            NVP_isMerged.append( 0 )
        
        for E, Epho, eta, phi in neuPixels_merged:
            NVP_pt .append( E/np.cosh(eta) )
            NVP_eta.append( eta )
            NVP_phi.append( phi )
            NVP_photonFraction.append( Epho/(E+1e-8) )
            NVP_isMerged.append( 1 )

        return NVP_pt, NVP_eta, NVP_phi, NVP_photonFraction, NVP_isMerged

    def linkPFGen(self):
        """
        return the index of the PF and Gen with the same ieta iphi
        """
        NVP_PF_toGenIndex = []
        NVP_Gen_toPFIndex = []

        PF_index = -1
        Gen_index = -1

        for ieta, iphi in zip(*np.where( (self.neuPixels_PF>0)|(self.neuPixels_Gen>0) )):
            if self.neuPixels_PF[ieta, iphi]>0 and self.neuPixels_Gen[ieta, iphi]>0:
                PF_index += 1
                Gen_index += 1
                NVP_PF_toGenIndex.append( Gen_index )
                NVP_Gen_toPFIndex.append( PF_index  )
            elif self.neuPixels_PF[ieta, iphi]>0 and self.neuPixels_Gen[ieta, iphi]<=0:
                PF_index += 1
                NVP_PF_toGenIndex.append( -1 )
            elif self.neuPixels_PF[ieta, iphi]<=0 and self.neuPixels_Gen[ieta, iphi]>0:
                Gen_index += 1
                NVP_Gen_toPFIndex.append( -1 )

        for iGenm, iPFm in zip(self.neuPixels_Gen_merged, self.neuPixels_PF_merged):
            PF_index += 1
            Gen_index += 1
            NVP_PF_toGenIndex.append( Gen_index )
            NVP_Gen_toPFIndex.append( PF_index )

        return NVP_PF_toGenIndex, NVP_Gen_toPFIndex
                

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

        # PF
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

        self.neuPixels_PF, self.phoPixels_PF = self.Neutrals2Pixels( neuPFCands )

        # Gen
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

        self.neuPixels_Gen, self.phoPixels_Gen = self.Neutrals2Pixels( neuGenCands )

        self.neuPixels_Gen_merged = []
        self.neuPixels_PF_merged  = []

        self.MergePixels()

        NVP_PF_pt, NVP_PF_eta, NVP_PF_phi, NVP_PF_photonFraction, NVP_PF_isMerged = self.Pixels2NVPs( self.neuPixels_PF, self.phoPixels_PF, self.neuPixels_PF_merged)
        NVP_Gen_pt, NVP_Gen_eta, NVP_Gen_phi, NVP_Gen_photonFraction, NVP_Gen_isMerged = self.Pixels2NVPs( self.neuPixels_Gen, self.phoPixels_Gen, self.neuPixels_Gen_merged)

        self.out.fillBranch("nNVP_PF",    len(NVP_PF_pt))
        self.out.fillBranch("NVP_PF_pt",  NVP_PF_pt )
        self.out.fillBranch("NVP_PF_eta", NVP_PF_eta)
        self.out.fillBranch("NVP_PF_phi", NVP_PF_phi)
        self.out.fillBranch("NVP_PF_photonFraction", NVP_PF_photonFraction)
        self.out.fillBranch("NVP_PF_isMerged", NVP_PF_isMerged)

        self.out.fillBranch("nNVP_Gen",    len(NVP_Gen_pt))
        self.out.fillBranch("NVP_Gen_pt",  NVP_Gen_pt )
        self.out.fillBranch("NVP_Gen_eta", NVP_Gen_eta)
        self.out.fillBranch("NVP_Gen_phi", NVP_Gen_phi)
        self.out.fillBranch("NVP_Gen_photonFraction", NVP_Gen_photonFraction)
        self.out.fillBranch("NVP_Gen_isMerged", NVP_Gen_isMerged)

        # build the map between PF and Gen
        NVP_PF_toGenIndex, NVP_Gen_toPFIndex = self.linkPFGen()
        self.out.fillBranch("NVP_PF_toGenIndex",    NVP_PF_toGenIndex  )
        self.out.fillBranch("NVP_Gen_toPFIndex",    NVP_Gen_toPFIndex  )

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

neuToNVPModuleConstr = lambda : neuToNVPProducer()
