import ROOT
import os, array, copy
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi
from math import *

class SimpleVBoson:
    def __init__(self,legs,pdgId=24):
        self.pdgId=pdgId
        self.legs = legs
        if len(legs)<2:
            print "ERROR: making a VBoson w/ < 2 legs!"
        self.pt1 = legs[0].Pt()
        self.pt2 = legs[1].Pt()
        self.dphi = self.legs[0].DeltaPhi(self.legs[1])        
        self.deta = self.legs[0].Eta()-self.legs[1].Eta()
        self.px1 = legs[0].Px(); self.py1 = legs[0].Py();
        self.px2 = legs[1].Px(); self.py2 = legs[1].Py();
    def pt(self):
        return (self.legs[0]+self.legs[1]).Pt()
    def phi(self):
        return (self.legs[0]+self.legs[1]).Phi()
    def y(self):
        return (self.legs[0]+self.legs[1]).Rapidity()
    def mt(self):
        return sqrt(2*self.pt1*self.pt2*(1-cos(self.dphi)))
    def ux(self):
        return (-self.px1-self.px2)
    def uy(self):
        return (-self.py1-self.py2)    
    def mll(self):
        return sqrt(2*self.pt1*self.pt2*(cosh(self.deta)-cos(self.dphi)))


class GenQEDJetProducer(Module):
    def __init__(self,deltaR,beamEn=7000.):
        self.beamEn=beamEn
        self.deltaR = deltaR
        self.vars = ("pt","eta","phi","mass","pdgId","fsrDR","fsrPt")
        self.genwvars = ("charge","pt","mass","y","decayId")
        ## if "genQEDJetHelper_cc.so" not in ROOT.gSystem.GetLibraries():
        ##     print "Load C++ Worker"
        ##     ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/WMass/python/postprocessing/helpers/genQEDJetHelper.cc+" % os.environ['CMSSW_BASE'])
        ## else:
        ##     print "genQEDJetHelper_cc.so found in ROOT libraries"
        ## self._worker = ROOT.GenQEDJetHelper(deltaR)
        self.massWeights = range(-20,21) ## order them by integer range(80300, 80505, 5) #masses in MeV
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.initReaders(inputTree) # initReaders must be called in beginFile
        self.out = wrappedOutputTree
        self.out.branch("GenPart_isPromptFinalState",        "O", lenVar="nGenPart")
        self.out.branch("GenPart_isPromptDecayed",           "O", lenVar="nGenPart")
        self.out.branch("GenPart_fromHardProcessFinalState", "O", lenVar="nGenPart")
        #self.out.branch("GenPart_fromHardProcessDecayed",    "O", lenVar="nGenPart")
        self.out.branch("pdfCentralWgt"  , "F")
        self.out.branch("partonId1"      , "I")
        self.out.branch("partonId2"      , "I")
        self.out.branch("nGenLepDressed" , "I")
        self.out.branch("nGenLepBare"    , "I")
        self.out.branch("nGenPromptNu"   , "I")
        for V in self.vars:
            self.out.branch("GenLepDressed_"+V, "F", lenVar="nGenLepDressed")
            self.out.branch("GenLepPreFSR_" +V, "F", lenVar="nGenLepPreFSR")
            self.out.branch("GenLepBare_"   +V, "F", lenVar="nGenLepBare")
            self.out.branch("GenPromptNu_"  +V, "F", lenVar="nGenPromptNu")
        for V in self.genwvars:
            self.out.branch("genw_"   +V, "F")
            self.out.branch("prefsrw_"+V, "F")
        for imass in self.massWeights:
            masssign = 'm' if imass < 0 else 'p' if imass > 0 else ''
            self.out.branch("mass_{s}{mass}".format(s=masssign,mass=abs(imass)), "F")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def initReaders(self,tree): # this function gets the pointers to Value and ArrayReaders and sets them in the C++ worker class
        ## marc try:
        ## marc     self.nGenPart = tree.valueReader("nGenPart")
        ## marc     for B in ("pt","eta","phi","mass","pdgId","isPromptHard","motherId","status") : setattr(self,"GenPart_"+B, tree.arrayReader("GenPart_"+B))
        ## marc     self._worker.setGenParticles(self.nGenPart,self.GenPart_pt,self.GenPart_eta,self.GenPart_phi,self.GenPart_mass,self.GenPart_pdgId,self.GenPart_isPromptHard,self.GenPart_motherId,self.GenPart_status)
        ## marc     #self.nLHEweight = tree.valueReader("nLHEweight")
        ## marc     #self.LHEweight_wgt = tree.arrayReader("LHEweight_wgt")
        ## marc except:
        ## marc     print '[genFriendProducer][Warning] Unable to attach to generator-level particles (data only?). No info will be produced'
        self._ttreereaderversion = tree._ttreereaderversion # self._ttreereaderversion must be set AFTER all calls to tree.valueReader or tree.arrayReader

    def bwWeight(self,genMass,imass, isW):
        # default mass calculated from MG5 inputs
        # width calculated with MG5_aMC_v2_6_3_2 loop_sm-ckm_no_b_mass for w+ > all all --> 2.05 +/- 7.65e-06 (GeV)
        (m0,gamma) = (80419.,2050.0) if isW else (91187.6, 2495.2) # MeV . the Z values are from the PDG, not from the MC production!!!
        newmass = m0 + imass*5.
        s_hat = pow(genMass,2)
        return (pow(s_hat - m0*m0,2) + pow(gamma*m0,2)) / (pow(s_hat - newmass*newmass,2) + pow(gamma*newmass,2))

    def getNeutrino(self):
        nus = []
        for p in self.genParts:
            if abs(p.pdgId) in [12, 14, 16] : ##and p.isPromptFinalState > 0: ## requiring isPromptFinalState results in 12% of events without neutrino
                nus.append(p)

        ##nus = sorted(nus, key = lambda x: x.pt)
        nus.sort(key = lambda x: x.pt, reverse=True)

        if len(nus) < 1:
            return 0, 0

        nu = nus[0]
        nuvec = ROOT.TLorentzVector()
        nuvec.SetPtEtaPhiM(nu.pt, nu.eta, nu.phi, nu.mass if nu.mass > 0. else 0.)

        return nuvec, nu.pdgId

    def getBareLeptons(self, strictlyPrompt=True):
        leptons = []
        for p in self.genParts:
            if not abs(p.pdgId) in [11,13,15] : continue
            if strictlyPrompt:
                if abs(p.pdgId) in [11, 13] and not p.isPromptFinalState : continue
                if abs(p.pdgId) == 15       and not p.isPromptDecayed : continue
            lepton = ROOT.TLorentzVector()
            lepton.SetPtEtaPhiM(p.pt, p.eta, p.phi, p.mass if p.mass >= 0. else 0.)
            leptons.append( [lepton, p.pdgId] )

        #leptons = sorted(leptons, key = lambda x: x[0].Pt() )
        leptons.sort(key = lambda x: x[0].Pt(), reverse=True )
        return leptons

    def getDressedLeptons(self, bareLeptons, strictlyPrompt=True, cone=0.1):

        if len(bareLeptons) ==0: 
            return bareLeptons
            ##print 'ERROR: DID NOT FIND A LEPTON FOR DRESSING !!!!!'

        leptons = copy.deepcopy(bareLeptons)
        for p in self.genParts:
            if not abs(p.pdgId)==22: continue
            if strictlyPrompt and not p.isPromptFinalState: continue
            tmp_photon = ROOT.TLorentzVector()
            tmp_photon.SetPtEtaPhiM(p.pt, p.eta, p.phi, p.mass if p.mass >= 0. else 0.)

            #first match in cone
            for il in xrange(0,len(leptons)):
                if leptons[il][0].DeltaR(tmp_photon) > cone: continue
                leptons[il][0] = leptons[il][0]+tmp_photon
                break

        ##print 'this is dressed leptons', leptons

        return leptons
    
    def getPreFSRLepton(self):
        
        lepInds = []
        for ip,p in enumerate(self.genParts):
            if abs(p.pdgId) in [11, 13] and p.fromHardProcessFinalState:
                lepInds.append( (p, ip) )
            if abs(p.pdgId) == 15       and p.fromHardProcessDecayed:
                lepInds.append( (p, ip) )
            ##if abs(p.pdgId) == 22:
            ##    if p.genPartIdxMother >=0:
            ##        print 'photon with mother', self.genParts[int(p.genPartIdxMother)].pdgId, 'and status', p.status

        if len(lepInds) < 1:
            return []

        retLeps = []

        for (lep, ind) in lepInds:
            finallep = lep

            ## see if there's a mother with same pdgId and lastBeforeFSR
            while (self.genParts[int(lep.genPartIdxMother)].pdgId == lep.pdgId):
                #if self.genParts[int(lep.genPartIdxMother)].lastBeforeFSR:
                if ( self.genParts[int(lep.genPartIdxMother)].statusFlags & 16384 ):
                    finallep = self.genParts[int(lep.genPartIdxMother)]
                    break
                ind = int(lep.genPartIdxMother)
                lep = self.genParts[ind]
                if lep.genPartIdxMother < 0:
                    break

            lepton = ROOT.TLorentzVector()
            lepton.SetPtEtaPhiM(finallep.pt, finallep.eta, finallep.phi, finallep.mass)

            retLeps.append( (lepton, finallep.pdgId) )

        ## if len(retLeps) == 4:
        ##     print '============'
        ##     for (lep, ind) in lepInds:
        ##         print 'this is pt {pt:4.2f} pdg id {pdg} motherpdgid {mid}'.format(pt=lep.pt,pdg= lep.pdgId, mid= self.genParts[int(lep.genPartIdxMother)].pdgId)

        ##     print 'below the original leptons'
        ##     for ip,p in enumerate(self.genParts):
        ##         if abs(p.pdgId) in [11, 13]:
        ##             print 'pt {pt:4.2f} eta {eta:4.2f} phi {phi:4.2f} status {st}'.format(pt= p.pt, eta=p.eta, phi=p.phi, st=p.status)
        ##             lepInds.append( (p, ip) )
        ##     print '================='

        ##     print 'invariant mass of the first two leptons:', (retLeps[0][0]+retLeps[1][0]).M()
        ##     print 'invariant mass of the last  two leptons:', (retLeps[2][0]+retLeps[3][0]).M()
            

        return retLeps

    def getFSRPhotons(self, lepton, leptonPdgId, cone=10):

        photons = []
        for p in self.genParts:
            if abs(p.pdgId) != 22 : continue
            if int(p.genPartIdxMother) < 0 : continue
            if self.genParts[int(p.genPartIdxMother)].pdgId != leptonPdgId: continue 
            if not p.isPromptFinalState : continue
            photon = ROOT.TLorentzVector()
            photon.SetPtEtaPhiM(p.pt, p.eta, p.phi, p.mass if p.mass >= 0. else 0.)
            if photon.DeltaR(lepton) < cone: 
                photons.append(photon)

        photons.sort(key = lambda x: x.Pt(), reverse=True )

        return photons

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        #nothing to do if this is data
        if event.isData: return True

        self.genParts = Collection(event, 'GenPart')

        # add the isPromptFinalState and isPromptDecayed flags
        genParts_isPromptFinalState = []
        genParts_isPromptDecayed = []
        genParts_fromHardProcessFinalState = []
        genParts_fromHardProcessDecayed = []
        for p in self.genParts:  
            p.isPromptFinalState = ( ( p.statusFlags & 1 ) and ( p.status==1 ) ) # isPrompt && status==1
            p.isPromptDecayed    = ( ( p.statusFlags & 2 ) and ( p.statusFlags & 1 ) ) # isDecayedLeptonHadron && isPrompt
            p.fromHardProcessFinalState = ( (p.statusFlags & 256 ) and ( p.status==1 ) ) # fromHardProcess && status==1
            p.fromHardProcessDecayed    = ( (p.statusFlags & 2 )   and ( p.statusFlags & 256 ) ) # isDecayedLeptonHadron && fromHardProcess
            
            genParts_isPromptFinalState.append( p.isPromptFinalState )
            genParts_isPromptDecayed.append( p.isPromptDecayed )
            genParts_fromHardProcessFinalState.append( p.fromHardProcessFinalState )
            genParts_fromHardProcessDecayed.append( p.fromHardProcessDecayed )

        self.out.fillBranch('GenPart_isPromptFinalState', genParts_isPromptFinalState)
        self.out.fillBranch('GenPart_isPromptDecayed',    genParts_isPromptDecayed)
        self.out.fillBranch('GenPart_fromHardProcessFinalState', genParts_fromHardProcessFinalState)
        #self.out.fillBranch('GenPart_fromHardProcessDecayed',    genParts_fromHardProcessDecayed)

        #this is a stupid check if variables are present. do this on the first gen particle
        makeBosons = False
        for ip,p in enumerate(self.genParts):
            if ip:
                break
            if hasattr(p, 'isPromptFinalState'):
                makeBosons = True
            

        if event._tree._ttreereaderversion > self._ttreereaderversion: # do this check at every event, as other modules might have read further branches
            self.initReaders(event._tree)

        ## check if the proper Ws can be made
        bareLeptonCollection                = self.getBareLeptons(strictlyPrompt=makeBosons)
        dressedLeptonCollection             = self.getDressedLeptons(bareLeptonCollection,strictlyPrompt=makeBosons)
        dressedLepton, dressedLeptonPdgId   = dressedLeptonCollection[0] if len(dressedLeptonCollection) else (0,0)
        bareLepton,    bareLeptonPdgId      = bareLeptonCollection[0]    if len(bareLeptonCollection)    else (0,0)
        #(preFSRLepton , preFSRLeptonPdgId ) = self.getPreFSRLepton()     if makeBosons else (0,0)
        (neutrino     , neutrinoPdgId     ) = self.getNeutrino()         if makeBosons else (0,0)
        ## new to work also for Z
        preFSRLeptonsAndPdgIds = self.getPreFSRLepton() if makeBosons else []

        ## take the two hardest preFSR leptons if there are 4, which happens with photos
        if len(preFSRLeptonsAndPdgIds) == 4:
            preFSRLeptonsAndPdgIds = preFSRLeptonsAndPdgIds[:2]
        ## same for the dressed leptons
        if len(dressedLeptonCollection) == 4:
            dressedLeptonCollection = dressedLeptonCollection[:2]
        ## =================================


        #always produce the dressed lepton collection
        if len(dressedLeptonCollection):
            self.out.fillBranch("nGenLepDressed", len(dressedLeptonCollection))
            retL={}
            #leptonsToTake = range(0,1 if makeBosons else len(dressedLeptonCollection))
            retL["pt"]    = [leppy[0].Pt()  for leppy in dressedLeptonCollection ]
            retL["eta"]   = [leppy[0].Eta() for leppy in dressedLeptonCollection ]
            retL["phi"]   = [leppy[0].Phi() for leppy in dressedLeptonCollection ]
            retL["pdgId"] = [leppy[1]       for leppy in dressedLeptonCollection ]
            retL["mass"]  = [leppy[0].M() if leppy[0].M() >= 0. else 0. ]

            retL["fsrDR"] = []
            retL["fsrPt"] = []
            for  leppy in dressedLeptonCollection:
                photonsFSR = self.getFSRPhotons(leppy[0],leppy[1])
                retL["fsrDR"].append(photonsFSR[0].DeltaR(leppy[0]) if len(photonsFSR)>0 else -1)
                retL["fsrPt"].append(photonsFSR[0].Pt() if len(photonsFSR)>0 else -1)

            for V in self.vars:
                self.out.fillBranch("GenLepDressed_"+V, retL[V])

        #always produce the bare lepton collection
        if len(bareLeptonCollection)>0:
            retB={}
            leptonsToTake = range(0,1 if makeBosons else len(bareLeptonCollection))
            retB["pt"]    = [bareLeptonCollection[i][0].Pt() for i in leptonsToTake ]
            retB["eta"]   = [bareLeptonCollection[i][0].Eta() for i in leptonsToTake ]
            retB["phi"]   = [bareLeptonCollection[i][0].Phi() for i in leptonsToTake ]
            retB["mass"]  = [bareLeptonCollection[i][0].M() if dressedLeptonCollection[i][0].M() >= 0. else 0. for i in leptonsToTake]
            retB["pdgId"] = [bareLeptonCollection[i][1] for i in leptonsToTake ]
            retB["fsrDR"] = []
            retB["fsrPt"] = []
            for  i in leptonsToTake:
                photonsFSR = self.getFSRPhotons(bareLeptonCollection[i][0],bareLeptonCollection[i][1])
                retB["fsrDR"].append(photonsFSR[0].DeltaR(bareLeptonCollection[i][0]) if len(photonsFSR)>0 else -1)
                retB["fsrPt"].append(photonsFSR[0].Pt() if len(photonsFSR)>0 else -1)

            self.out.fillBranch("nGenLepBare", leptonsToTake[-1])
            for V in self.vars:
                self.out.fillBranch("GenLepBare_"+V, retB[V])

        #W-specific
        if neutrino:
            retN={}
            retN["pt"]    = [neutrino.Pt()  ]
            retN["eta"]   = [neutrino.Eta() ]
            retN["phi"]   = [neutrino.Phi() ]
            retN["mass"]  = [neutrino.M() if neutrino.M() >= 0. else 0.  ] ## weird protection... also already checked. but better safe than sorry
            retN["pdgId"] = [neutrinoPdgId  ]
            retN["fsrDR"] = [-1 ]
            retN["fsrPt"] = [-1 ]
            self.out.fillBranch("nGenPromptNu", 1)
            for V in self.vars:
                self.out.fillBranch("GenPromptNu_"+V, retN[V])
            #self.out.fillBranch("GenPromptNu_pdgId", [pdgId for pdgId in nuPdgIds])            

        if len(preFSRLeptonsAndPdgIds):
            self.out.fillBranch("nGenLepPreFSR", len(preFSRLeptonsAndPdgIds))
            retP={}
            retP["pt"]    = [leppy[0].Pt()  for leppy in preFSRLeptonsAndPdgIds]
            retP["eta"]   = [leppy[0].Eta() for leppy in preFSRLeptonsAndPdgIds]
            retP["phi"]   = [leppy[0].Phi() for leppy in preFSRLeptonsAndPdgIds]
            retP["pdgId"] = [leppy[1]       for leppy in preFSRLeptonsAndPdgIds]
            retP["mass"]  = [leppy[0].M()   if leppy[0].M() >= 0. else 0. for leppy in preFSRLeptonsAndPdgIds]
            retP["fsrDR"] = []
            retP["fsrPt"] = []
            for leppy in preFSRLeptonsAndPdgIds:
                photonsFSR = self.getFSRPhotons(leppy[0],leppy[1])
                retP["fsrDR"].append(photonsFSR[0].DeltaR(leppy[0]) if len(photonsFSR)>0 else -1 )
                retP["fsrPt"].append(photonsFSR[0].Pt() if len(photonsFSR)>0 else -1 )
        
            for V in self.vars:
                self.out.fillBranch("GenLepPreFSR_"+V, retP[V])


        if len(preFSRLeptonsAndPdgIds) == 2 or (len(preFSRLeptonsAndPdgIds) == 1 and neutrino):

            if len(preFSRLeptonsAndPdgIds) == 2: ## these are Zs
                isWBoson = False
                prefsrw = preFSRLeptonsAndPdgIds[0][0] + preFSRLeptonsAndPdgIds[1][0]
                self.out.fillBranch("prefsrw_charge" , 0)
                l1, l1pdg, l2 = preFSRLeptonsAndPdgIds[0][0], preFSRLeptonsAndPdgIds[0][1], preFSRLeptonsAndPdgIds[1][0]
                # convention for phiCS: use l- direction for W-, use neutrino for W+
                (lplus,lminus) = (l2,l1) if l1pdg<0 else (l1,l2)
                self.out.fillBranch("prefsrw_decayId", abs(l1pdg))

            else: ## these are Ws
                isWBoson = True
                prefsrw = preFSRLeptonsAndPdgIds[0][0] + neutrino
                self.out.fillBranch("prefsrw_charge" , float(-1*np.sign(preFSRLeptonsAndPdgIds[0][1])))
                l1, l1pdg = preFSRLeptonsAndPdgIds[0][0], preFSRLeptonsAndPdgIds[0][1]
                # convention for phiCS: use l- direction for W-, use neutrino for W+
                (lplus,lminus) = (neutrino,l1) if l1pdg<0 else (l1,neutrino)
                self.out.fillBranch("prefsrw_decayId", abs(neutrinoPdgId))

            self.out.fillBranch("prefsrw_pt"     , prefsrw.Pt())
            self.out.fillBranch("prefsrw_y"      , prefsrw.Rapidity())
            self.out.fillBranch("prefsrw_mass"   , prefsrw.M())
            for imass in self.massWeights:
                masssign = 'm' if imass < 0 else 'p' if imass > 0 else ''
                self.out.fillBranch("mass_{s}{mass}".format(s=masssign,mass=abs(imass)), self.bwWeight(genMass=prefsrw.M()*1000,imass=imass,isW=isWBoson))

        if len(dressedLeptonCollection) == 2 or (len(dressedLeptonCollection) == 1 and neutrino):
            if len(dressedLeptonCollection) == 2: ## these are Zs
                isWBoson = False
                genw = dressedLeptonCollection[0][0] + dressedLeptonCollection[1][0]
                self.out.fillBranch("genw_charge" , 0)
                l1, l1pdg, l2 = dressedLeptonCollection[0][0], dressedLeptonCollection[0][1], dressedLeptonCollection[1][0]
                # convention for phiCS: use l- direction for W-, use neutrino for W+
                (lplus,lminus) = (l2,l1) if l1pdg<0 else (l1,l2)
                self.out.fillBranch("genw_decayId", abs(l1pdg))

            else: ## these are Ws
                isWBoson = True
                genw = dressedLeptonCollection[0][0] + neutrino
                self.out.fillBranch("genw_charge" , float(-1*np.sign(dressedLeptonCollection[0][1])))
                l1, l1pdg = dressedLeptonCollection[0][0], dressedLeptonCollection[0][1]
                # convention for phiCS: use l- direction for W-, use neutrino for W+
                (lplus,lminus) = (neutrino,l1) if l1pdg<0 else (l1,neutrino)
                self.out.fillBranch("genw_decayId", abs(neutrinoPdgId))

            self.out.fillBranch("genw_pt"     , genw.Pt())
            self.out.fillBranch("genw_y"      , genw.Rapidity())
            self.out.fillBranch("genw_mass"   , genw.M())
            for imass in self.massWeights:
                masssign = 'm' if imass < 0 else 'p' if imass > 0 else ''
                self.out.fillBranch("mass_{s}{mass}".format(s=masssign,mass=abs(imass)), self.bwWeight(genMass=genw.M()*1000,imass=imass,isW=isWBoson))

            ## remove these variables to save space
            #self.out.fillBranch("genw_mt"   , sqrt(2*lplus.Pt()*lminus.Pt()*(1.-cos(deltaPhi(lplus.Phi(),lminus.Phi())) )))
            #self.out.fillBranch("genw_eta",genw.Eta())
            #self.out.fillBranch("genw_phi",genw.Phi())
        else: ## no neutrino found!
            ##if not len(dressedLeptons): 
            ##    print '================================'
            ##    print 'no dressed leptons found!'
            ##    print 'no dressed leptons fround :lumi:evt: {a}:{b}:{c}'.format(a=getattr(event, "run"),b=getattr(event, "lumi"),c=getattr(event, "evt"))
            ##else:
            ##    print '================================'
            ##    print 'no neutrinos found, in run:lumi:evt: {a}:{b}:{c}'.format(a=getattr(event, "run"),b=getattr(event, "lumi"),c=getattr(event, "evt"))
            for V in self.genwvars:
                self.out.fillBranch("genw_"+V   , -999)
                self.out.fillBranch("prefsrw_"+V, -999)
            for imass in self.massWeights:
                masssign = 'm' if imass < 0 else 'p' if imass > 0 else ''
                self.out.fillBranch("mass_{s}{mass}".format(s=masssign,mass=abs(imass)), 1.)

        if hasattr(event,"genWeight"):
            self.out.fillBranch("partonId1" , getattr(event, "id1") if hasattr(event, 'id1') else -999)
            self.out.fillBranch("partonId2" , getattr(event, "id2") if hasattr(event, 'id2') else -999)
        else:
            self.out.fillBranch("partonId1" , -999 )
            self.out.fillBranch("partonId2" , -999 )

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

genQEDJets14TeV = lambda : GenQEDJetProducer(deltaR=0.1,beamEn=7000.)
genQEDJets = lambda : GenQEDJetProducer(deltaR=0.1,beamEn=6500.)
