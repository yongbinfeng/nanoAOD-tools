import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR,matchObjectCollection

import numpy as np

"""
neuToNVPProducer has to be run after the genPartMatching,
because it has used Gen and PF particle types
"""

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
                            1.930, 2.043, 2.172, 2.322, 2.500, 2.650, 2.868,
                            3.00001
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

    def getNEta(self):
        # return the number of eta towers
        return len(self.etaSeg)-1

    def getNPhi(self):
        return max(self.nPhiHB, self.nPhiHE)

    def getIEta(self, eta):
        assert eta>-3.00001 and eta<3.00001, "eta %f goes out of the range"%eta
        return np.searchsorted( self.etaSeg, eta )-1

    def getEta(self, ieta):
        # return the middle of the eta cell given the index
        return ( self.etaSeg[ieta] + self.etaSeg[ieta+1] )/2.0
    
    def getIPhi(self, eta, phi):
        # index of phi depends on the index of eta
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

    def formatIEtaIPhi(self, ietaiphicands):
        # given ieta, iphi, check if it makes sense:
        # e.g., ieta > nEtas, 
        # and use the periodic property to regulate iphi in [0, 2*nPhi)
        ietaiphis = []
        for ieta, iphi in ietaiphicands:
            if ieta<0 or ieta >= self.getNEta():
                continue
            eta = self.getEta(ieta)
            nPhi = self.nPhiHB if abs(eta)<self.etaCut else self.nPhiHE
            if iphi<0:
                iphi = iphi + nPhi
            elif iphi>=nPhi:
                iphi = iphi - nPhi
            ietaiphis.append((ieta, iphi))

        return ietaiphis

    def getConnectedIEtaIPhis(self, ieta, iphi):
        # return the ieta, iphi of the cells connected the to the (ieta, iphi) cell
        ietaiphicands = [(ieta-1, iphi), (ieta+1, iphi), (ieta, iphi-1), (ieta, iphi+1)]
        return self.formatIEtaIPhi( ietaiphicands )

    def getSurroundingIEtaIPhis(self, ieta, iphi):
        # the connected ieta, iphis + the 4 at the coner
        ietaiphicands = [(ieta-1, iphi-1), (ieta-1, iphi), (ieta-1, iphi+1),
                         (ieta,   iphi-1),                 (ieta,   iphi+1),
                         (ieta+1, iphi-1), (ieta+1, iphi), (ieta+1, iphi+1),
                         ]
        return self.formatIEtaIPhi( ietaiphicands )
        

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

class HCALTower(object):
    """
    A small struct to save the pt, eta, phi, nparticles,
    and puppiEnergy (puppiWeight * total energy) information, 
    from the particles in the tower region
    """
    def __init__(self, ieta, iphi, pt, eta, phi, mass, pindex, photonFraction = 0., puppiWeight=1.0, isMerged = False):
        # ieta, iphi is for the HCAL
        # pt, eta, phi, mass, photonFraction, puppiWeight is from the neutral particle in this tower
        self.ieta = ieta
        self.iphi = iphi
        self.segment = HCALSegmentation()
        self.eta, self.phi = self.segment.getEtaPhi(ieta, iphi)

        self.energy = np.hypot( pt * np.cosh(eta), mass ) 
        self.photonEnergy = photonFraction * self.energy
        self.puppiEnergy = puppiWeight * self.energy

        self.nPart = 1
        self.pIndices = [ pindex ]
        self.index = 0
        self.isMerged = False

    def __add__(self, other):
        assert self.ieta==other.ieta or self.iphi==other.iphi, "HCAL Towers with different eta and phis can not be added!"
        self.energy += other.energy
        self.photonEnergy += other.photonEnergy
        self.puppiEnergy += other.puppiEnergy
        self.nPart += other.nPart
        self.pIndices += other.pIndices
        return self

    def pt(self):
        return self.energy / np.cosh(self.eta)

    def photonFraction(self):
        return self.photonEnergy / (self.energy + 1e-8)

    def puppiWeight(self):
        return self.puppiEnergy / ( self.energy + 1e-8)

    def firstPIndex(self):
        # return the index of the fist particle in this tower
        return self.pIndices[0]


class neuToNVPProducer(Module):
    """ 
    convert the neutral particles ( neutral hadrons and photons ) into towers.
    Each tower will be treated as a 'virtual' particle in the training later
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
            self.out.branch("NVP_%s_nPart"%itype,    "I", lenVar="nNVP_%s"%itype)
            self.out.branch("NVP_%s_isMerged"%itype, "O", lenVar="nNVP_%s"%itype)
            self.out.branch("NVP_%s_firstPIndex"%itype, "I", lenVar="nNVP_%s"%itype)
        self.out.branch("NVP_PF_puppiWeightNoLep",   "F", lenVar="nNVP_PF")
        self.out.branch("NVP_PF_toGenIndex", "I", lenVar="nNVP_PF")
        self.out.branch("NVP_Gen_toPFIndex", "I", lenVar="nNVP_Gen")
        self.out.branch("packedGenPart_TowerIndex", "I", lenVar="npackedGenPart")
        self.out.branch("PF_TowerIndex",            "I", lenVar="nPF")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def Neutrals2Towers(self, NeuColl, isPF = False):
        # 58 cells between eta[-3.0, 3.0]
        # 72/36 cells in phi HB/HE
        # save the tower object
        neuTowers = np.ndarray((self.segment.getNEta(), self.segment.getNPhi()), dtype=np.object)
        for p in NeuColl:
            ieta, iphi = self.segment.getIEtaIPhi( p.eta, p.phi )
            puppiWeight = p.puppiWeightNoLep if isPF else 1.0
            ptower = HCALTower( ieta, iphi, p.pt, p.eta, p.phi, p.mass, p.index, (abs(p.pdgId)==22), puppiWeight )
            if neuTowers[ieta, iphi]==None:
                neuTowers[ieta, iphi] = ptower
            else:
                neuTowers[ieta, iphi] += ptower

        # index the tower based on ieta, iphi
        index = 0
        for ieta, iphi in zip(*np.where(neuTowers[:,:])):
            tower = neuTowers[ieta, iphi]
            tower.index = index
            for ip in tower.pIndices:
                if isPF:
                    self.pfCands[ip].TIndex = index
                else:
                    self.packedGenParts[ip].TIndex = index
            index += 1

        return neuTowers

    def MergeTowers(self):
        """
        Merge the towers for Gen and PF if one is Gen > 1.25 * PF, and
        if merging the connected cell could make the merged one smaller
        """
        for ieta, iphi in zip(*np.where(self.neuTowers_Gen[:,:,0]>0)):

            E_Gen = self.neuTowers_Gen[ieta, iphi, 0]
            if E_Gen/np.cosh(self.segment.getEta(ieta)) < 0.5:
                # only merge pT > 0.5 GeV cells
                continue
            if self.neuTowers_PF[ieta, iphi, 0] > 0.8 * E_Gen:
                continue 

            r_best = 1e8
            ietaiphim = None
            for ietac, iphic in self.segment.getConnectedIEtaIPhis(ieta, iphi):
                r_merged = (E_Gen+self.neuTowers_Gen[ietac, iphic, 0])/(self.neuTowers_PF[ieta, iphi, 0]+self.neuTowers_PF[ietac,iphic,0]+1e-8)
                r_org  = E_Gen / ( self.neuTowers_PF[ieta, iphi, 0] + 1e-8)
                r_orgc = self.neuTowers_Gen[ietac, iphic, 0] / ( self.neuTowers_PF[ietac,iphic,0]+1e-8 )
                if r_merged < r_org and r_merged > r_orgc and abs(r_merged-1.0) < abs(r_best-1.0):
                    r_best = r_merged
                    ietaiphim = (ietac, iphic)

            if ietaiphim:
                ietam = ietaiphim[0]
                iphim = ietaiphim[1]

                neuGen_merged = self.neuTowers_Gen[ieta, iphi] + self.neuTowers_Gen[ietam, iphim]
                neuPF_merged  = self.neuTowers_PF [ieta, iphi] + self.neuTowers_PF [ietam, iphim]
                eta_merged = self.segment.getAvgEta(ieta, ietam)
                phi_merged = self.segment.getAvgPhi(ieta, iphi, ietam, iphim)

                self.neuTowers_Gen_merged.append((list(neuGen_merged), eta_merged, phi_merged))
                self.neuTowers_PF_merged .append((list(neuPF_merged),  eta_merged, phi_merged))


                # clear the merged towers info
                self.neuTowers_Gen[ieta,  iphi ]= 0
                self.neuTowers_Gen[ietam, iphim]=0
                self.neuTowers_PF[ieta,  iphi ] = 0
                self.neuTowers_PF[ietam, iphim]=0


    def Towers2NVPs(self, neuTowers, neuTowers_merged = []):
        """ 
        convert the energy deposit in cells to neutral 'virtual particles' 
        """
        NVP_pt       = []
        NVP_eta      = []
        NVP_phi      = []
        NVP_photonFraction  = []
        NVP_nPart    = []
        NVP_puppiWeight = []
        NVP_isMerged = []
        NVP_firstPIndex = []
        for ieta, iphi in zip(*np.where(neuTowers[:,:])):
            eta, phi = self.segment.getEtaPhi( ieta, iphi )
            tower = neuTowers[ieta, iphi]

            NVP_pt .append( tower.pt()  )
            NVP_eta.append( tower.eta )
            NVP_phi.append( tower.phi )
            NVP_photonFraction.append( tower.photonFraction() ) 
            NVP_nPart.append( tower.nPart )
            NVP_puppiWeight.append( tower.puppiWeight() )
            NVP_isMerged.append( 0 )
            NVP_firstPIndex.append( tower.firstPIndex() )
        
        for neumerged, eta, phi in neuTowers_merged:
            NVP_pt .append( neumerged[0]/np.cosh(eta) )
            NVP_eta.append( eta )
            NVP_phi.append( phi )
            NVP_photonFraction.append( neumerged[1]/(neumerged[0]+1e-8) )
            NVP_nPart.append( int(neumerged[2]) )
            NVP_pt_puppi.append( neumerged[3]/np.cosh(eta) )
            NVP_isMerged.append( 1 )

        return NVP_pt, NVP_eta, NVP_phi, NVP_photonFraction, NVP_nPart, NVP_puppiWeight, NVP_isMerged, NVP_firstPIndex

    def linkPFGenTowers(self, GenColl, PFColl):
        """
        return the index of the PF and Gen with the same ieta iphi
        """
        NVP_PF_toGenIndex = []
        NVP_Gen_toPFIndex = []

        PF_index = -1
        Gen_index = -1

        for ieta, iphi in zip(*np.where( (self.neuTowers_PF[:,:]>0)|(self.neuTowers_Gen[:,:]>0) )):
            if self.neuTowers_PF[ieta, iphi] and self.neuTowers_Gen[ieta, iphi]:
                PF_index = self.neuTowers_PF[ieta, iphi].index
                Gen_index = self.neuTowers_Gen[ieta, iphi].index
                NVP_PF_toGenIndex.append( Gen_index )
                NVP_Gen_toPFIndex.append( PF_index  )
            elif self.neuTowers_PF[ieta, iphi] and not self.neuTowers_Gen[ieta, iphi]:
                NVP_PF_toGenIndex.append( -1 )
            elif ( not self.neuTowers_PF[ieta, iphi] ) and self.neuTowers_Gen[ieta, iphi]:
                NVP_Gen_toPFIndex.append( -1 )

        for iGenm, iPFm in zip(self.neuTowers_Gen_merged, self.neuTowers_PF_merged):
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
        self.pfCands  = Collection(event, "PF")
        self.packedGenParts = Collection(event, "packedGenPart") # gen particles from packedGenParticles (status==1)

        # check if the event has a good muon
        hasGoodMuon = ( muons[0].pt>20.0 and abs(muons[0].eta)<2.4 and muons[0].tightId )
        vetoCands = []
        if hasGoodMuon:
            vetoCands.append( muons[0] )

        # PF
        neuPFCands = []
        index = 0
        for p in self.pfCands:
            p.TIndex = -999
            p.index = index
            index += 1
            if abs(p.eta)>3.0: continue
            # remove charged, HF, and photon candidates
            #if p.charge!=0 and abs(p.eta)<2.5: continue
            if p.charge!=0: continue
            if p.pdgId in [1,2,22]: continue
            if p.ptype!=2: continue
            if self.vetoCandidate(p, vetoCands ): continue
            neuPFCands.append( p )

        self.neuTowers_PF = self.Neutrals2Towers( neuPFCands, isPF = True )

        # Gen
        vetoGens = []
        for gp in self.packedGenParts:
            if abs(gp.pdgId) in [11, 13] and gp.pt > 20.0 and abs(gp.eta)<2.4 :
                vetoGens.append( gp )
                break

        neuGenCands = []
        index = 0
        for gp in self.packedGenParts:
            gp.TIndex = -999
            gp.index = index
            index += 1
            if abs(gp.eta)>3.0: continue
            if gp.charge!=0 and abs(gp.eta)<2.5: continue
            if abs(gp.pdgId) in [12, 14, 16]: continue # veto neutrino
            if abs(gp.pdgId)==22: continue # veto photon
            if gp.ptype!=2: continue
            #if abs(gp.pdgId)==22 and self.vetoCandidate(gp, vetoGens): continue
            neuGenCands.append( gp )

        self.neuTowers_Gen = self.Neutrals2Towers( neuGenCands )

        self.neuTowers_Gen_merged = []
        self.neuTowers_PF_merged  = []

        #self.MergeTowers()

        NVP_PF_pt, NVP_PF_eta, NVP_PF_phi, NVP_PF_photonFraction, NVP_PF_nPart, NVP_PF_puppiWeight, NVP_PF_isMerged, NVP_PF_firstPIndex = self.Towers2NVPs( self.neuTowers_PF, self.neuTowers_PF_merged)
        NVP_Gen_pt, NVP_Gen_eta, NVP_Gen_phi, NVP_Gen_photonFraction, NVP_Gen_nPart, NVP_Gen_puppiWeight, NVP_Gen_isMerged, NVP_Gen_firstPIndex = self.Towers2NVPs( self.neuTowers_Gen, self.neuTowers_Gen_merged)

        self.out.fillBranch("nNVP_PF",    len(NVP_PF_pt))
        self.out.fillBranch("NVP_PF_pt",  NVP_PF_pt )
        self.out.fillBranch("NVP_PF_eta", NVP_PF_eta)
        self.out.fillBranch("NVP_PF_phi", NVP_PF_phi)
        self.out.fillBranch("NVP_PF_photonFraction", NVP_PF_photonFraction)
        self.out.fillBranch("NVP_PF_nPart",    NVP_PF_nPart)
        self.out.fillBranch("NVP_PF_puppiWeightNoLep", NVP_PF_puppiWeight)
        self.out.fillBranch("NVP_PF_isMerged", NVP_PF_isMerged)
        self.out.fillBranch("NVP_PF_firstPIndex", NVP_PF_firstPIndex)

        self.out.fillBranch("nNVP_Gen",    len(NVP_Gen_pt))
        self.out.fillBranch("NVP_Gen_pt",  NVP_Gen_pt )
        self.out.fillBranch("NVP_Gen_eta", NVP_Gen_eta)
        self.out.fillBranch("NVP_Gen_phi", NVP_Gen_phi)
        self.out.fillBranch("NVP_Gen_photonFraction", NVP_Gen_photonFraction)
        self.out.fillBranch("NVP_Gen_nPart",    NVP_Gen_nPart)
        self.out.fillBranch("NVP_Gen_isMerged", NVP_Gen_isMerged)
        self.out.fillBranch("NVP_Gen_firstPIndex", NVP_Gen_firstPIndex)

        # build the map between PF and Gen
        NVP_PF_toGenIndex, NVP_Gen_toPFIndex = self.linkPFGenTowers( self.packedGenParts, self.pfCands )
        self.out.fillBranch("NVP_PF_toGenIndex",    NVP_PF_toGenIndex  )
        self.out.fillBranch("NVP_Gen_toPFIndex",    NVP_Gen_toPFIndex  )

        # save the tower index of PF and Gen
        TIndices_Gen = [ gp.TIndex for gp in self.packedGenParts ]
        TIndices_PF  = [ p.TIndex  for p  in self.pfCands ]

        self.out.fillBranch("packedGenPart_TowerIndex", TIndices_Gen) 
        self.out.fillBranch("PF_TowerIndex",            TIndices_PF)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

neuToNVPModuleConstr = lambda : neuToNVPProducer()
