# output ntuple for H->tautau analysis for CMSSW_10_2_X

from ROOT import TLorentzVector
from math import sqrt, sin, cos, pi
import tauFun 
import ROOT
import os

class outTuple() :
    
    def __init__(self,fileName):
        from array import array
        from ROOT import TFile, TTree

        # Tau Decay types
        self.kUndefinedDecayType, self.kTauToHadDecay,  self.kTauToElecDecay, self.kTauToMuDecay = 0, 1, 2, 3    
        ROOT.gInterpreter.ProcessLine(".include .")
        for baseName in ['MeasuredTauLepton','svFitAuxFunctions','FastMTT'] : 
            if os.path.isfile("{0:s}_cc.so".format(baseName)) :
                ROOT.gInterpreter.ProcessLine(".L {0:s}_cc.so".format(baseName))
            else :
                ROOT.gInterpreter.ProcessLine(".L {0:s}.cc++".format(baseName))   # .L is not just for .so files, also .cc
     
        self.f = TFile( fileName, 'recreate' )
        self.t = TTree( 'Events', 'Output tree' )

        self.entries = 0 
        self.run = array('l',[0])
        self.lumi = array('l',[0])
        self.evt = array('l',[0])
        self.weight = array('f',[0])
        
        self.pt_1 = array('f',[0])
        self.phi_1 = array('f',[0])
        self.eta_1 = array('f',[0])
        self.m_1 = array('f',[0])
        self.q_1 = array('f',[0])
        self.d0_1 = array('f',[0])
        self.dZ_1 = array('f',[0])
        self.mt_1 = array('f',[0])
        self.pfmt_1 = array('f',[0])
        self.puppimt_1 = array('f',[0])
        self.iso_1 = array('f',[0])
        self.gen_match_1 = array('l',[0])
        self.againstElectronLooseMVA6_1 = array('f',[0])
        self.againstElectronMediumMVA6_1 = array('f',[0])
        self.againstElectronTightMVA6_1 = array('f',[0])
        self.againstElectronVLooseMVA6_1 = array('f',[0])
        self.againstElectronVTightMVA6_1 = array('f',[0])
        self.againstMuonLoose3_1 = array('f',[0])
        self.againstMuonTight3_1 = array('f',[0])
        self.byIsolationMVA3oldDMwLTraw_1 = array('f',[0])
        self.trigweight_1 = array('f',[0])
        self.idisoweight_1 = array('f',[0])

        self.pt_2 = array('f',[0])
        self.phi_2 = array('f',[0])
        self.eta_2 = array('f',[0])
        self.m_2 = array('f',[0])
        self.q_2 = array('f',[0])
        self.d0_2 = array('f',[0])
        self.dZ_2 = array('f',[0])
        self.mt_2 = array('f',[0])
        self.pfmt_2 = array('f',[0])
        self.puppimt_2 = array('f',[0])
        self.iso_2 = array('f',[0])
        self.gen_match_2 = array('l',[0])
        self.againstElectronLooseMVA6_2 = array('f',[0])
        self.againstElectronMediumMVA6_2 = array('f',[0])
        self.againstElectronTightMVA6_2 = array('f',[0])
        self.againstElectronVLooseMVA6_2 = array('f',[0])
        self.againstElectronVTightMVA6_2 = array('f',[0])
        self.againstMuonLoose3_2 = array('f',[0])
        self.againstMuonTight3_2 = array('f',[0])
        self.byIsolationMVA3oldDMwLTraw_2 = array('f',[0])
        self.trigweight_2 = array('f',[0])
        self.idisoweight_2 = array('f',[0])

        # di-tau variables
        self.pt_tt  = array('f',[0])
        self.mt_tot = array('f',[0])
        self.m_vis  = array('f',[0])
        self.m_sv   = array('f',[0])
        self.mt_sv  = array('f',[0])

        # di-lepton variables.   _p and _m refer to plus and minus charge
        # ll_lmass is mass of decay lepton 
        self.ll_lmass = array('f',[0])     
        self.mll       = array('f',[0])
        self.ll_pt_p   = array('f',[0])
        self.ll_phi_p  = array('f',[0])
        self.ll_eta_p  = array('f',[0])
        self.ll_pt_m   = array('f',[0])
        self.ll_phi_m  = array('f',[0])
        self.ll_eta_m  = array('f',[0])
        
        # MET variables
        self.met = array('f',[0])
        self.metphi = array('f',[0])
        self.puppimet = array('f',[0])
        self.puppimetphi = array('f',[0])
        self.metcov00 = array('f',[0])
        self.metcov01 = array('f',[0])
        self.metcov10 = array('f',[0])
        self.metcov11 = array('f',[0])

        # jet variables
        self.njetspt20 = array('f',[0])
        self.njets = array('f',[0])
        self.nbtag = array('f',[0])

        self.jpt_1  = array('f',[0])
        self.jeta_1  = array('f',[0])
        self.jphi_1  = array('f',[0])
        self.jcsv_1 = array('f',[0])
        self.jpt_2  = array('f',[0])
        self.jeta_2  = array('f',[0])
        self.jphi_2  = array('f',[0])
        self.jcsv_2 = array('f',[0])

        self.bpt_1  = array('f',[0])
        self.beta_1  = array('f',[0])
        self.bphi_1  = array('f',[0])
        self.bcsv_1 = array('f',[0])
        self.bpt_2  = array('f',[0])
        self.beta_2  = array('f',[0])
        self.bphi_2  = array('f',[0])
        self.bcsv_2 = array('f',[0])
      
        self.t.Branch('run',  self.run,    'run/l' )
        self.t.Branch('lumi', self.lumi,   'lumi/I' )
        self.t.Branch('evt',  self.evt,  'evt/I' )
        self.t.Branch('weight',  self.weight,  'weight/F' )

        self.t.Branch('pt_1', self.pt_1, 'pt_1/F')
        self.t.Branch('phi_1', self.phi_1, 'phi_1/F')
        self.t.Branch('eta_1', self.eta_1, 'eta_1/F')
        self.t.Branch('m_1', self.m_1, 'm_1/F')
        self.t.Branch('q_1', self.q_1, 'q_1/F')
        self.t.Branch('d0_1', self.d0_1, 'd0_1/F')
        self.t.Branch('dZ_1', self.dZ_1, 'dZ_1/F')
        self.t.Branch('mt_1', self.mt_1, 'mt_1/F')
        self.t.Branch('pfmt_1', self.pfmt_1, 'pfmt_1/F')
        self.t.Branch('puppimt_1', self.puppimt_1, 'puppimt_1/F')
        self.t.Branch('iso_1', self.iso_1, 'iso_1/F')
        self.t.Branch('gen_match_1', self.gen_match_1, 'gen_match_1/l')
        self.t.Branch('againstElectronLooseMVA6_1', self.againstElectronLooseMVA6_1, 'againstElectronLooseMVA6_1/F')
        self.t.Branch('againstElectronMediumMVA6_1', self.againstElectronMediumMVA6_1, 'againstElectronMediumMVA6_1/F')
        self.t.Branch('againstElectronTightMVA6_1', self.againstElectronTightMVA6_1, 'againstElectronTightMVA6_1/F')
        self.t.Branch('againstElectronVLooseMVA6_1', self.againstElectronVLooseMVA6_1, 'againstElectronVLooseMVA6_1/F')
        self.t.Branch('againstElectronVTightMVA6_1', self.againstElectronVTightMVA6_1, 'againstElectronVTightMVA6_1/F')
        self.t.Branch('againstMuonLoose3_1', self.againstMuonLoose3_1, 'againstMuonLoose3_1/F')
        self.t.Branch('againstMuonTight3_1', self.againstMuonTight3_1, 'againstMuonTight3_1/F')
        self.t.Branch('byIsolationMVA3oldDMwLTraw_1', self.byIsolationMVA3oldDMwLTraw_1, 'byIsolationMVA3oldDMwLTraw_1/F')
        self.t.Branch('trigweight_1', self.trigweight_1, 'trigweight_1/F')
        self.t.Branch('idisoweight_1', self.idisoweight_1, 'idisoweight_1/F')

        self.t.Branch('pt_2', self.pt_2, 'pt_2/F')
        self.t.Branch('phi_2', self.phi_2, 'phi_2/F')
        self.t.Branch('eta_2', self.eta_2, 'eta_2/F')
        self.t.Branch('m_2', self.m_2, 'm_2/F')
        self.t.Branch('q_2', self.q_2, 'q_2/F')
        self.t.Branch('d0_2', self.d0_2, 'd0_2/F')
        self.t.Branch('dZ_2', self.dZ_2, 'dZ_2/F')
        self.t.Branch('mt_2', self.mt_2, 'mt_2/F')
        self.t.Branch('pfmt_2', self.pfmt_2, 'pfmt_2/F')
        self.t.Branch('puppimt_2', self.puppimt_2, 'puppimt_2/F')
        self.t.Branch('iso_2', self.iso_2, 'iso_2/F')
        self.t.Branch('gen_match_2', self.gen_match_2, 'gen_match_2/l')
        self.t.Branch('againstElectronLooseMVA6_2', self.againstElectronLooseMVA6_2, 'againstElectronLooseMVA6_2/F')
        self.t.Branch('againstElectronMediumMVA6_2', self.againstElectronMediumMVA6_2, 'againstElectronMediumMVA6_2/F')
        self.t.Branch('againstElectronTightMVA6_2', self.againstElectronTightMVA6_2, 'againstElectronTightMVA6_2/F')
        self.t.Branch('againstElectronVLooseMVA6_2', self.againstElectronVLooseMVA6_2, 'againstElectronVLooseMVA6_2/F')
        self.t.Branch('againstElectronVTightMVA6_2', self.againstElectronVTightMVA6_2, 'againstElectronVTightMVA6_2/F')
        self.t.Branch('againstMuonLoose3_2', self.againstMuonLoose3_2, 'againstMuonLoose3_2/F')
        self.t.Branch('againstMuonTight3_2', self.againstMuonTight3_2, 'againstMuonTight3_2/F')
        self.t.Branch('byIsolationMVA3oldDMwLTraw_2', self.byIsolationMVA3oldDMwLTraw_2, 'byIsolationMVA3oldDMwLTraw_2/F')
        self.t.Branch('trigweight_2', self.trigweight_2, 'trigweight_2/F')
        self.t.Branch('idisoweight_2', self.idisoweight_2, 'idisoweight_2/F')

        # di-tau variables
        self.t.Branch('pt_tt', self.pt_tt, 'pt_tt/F')
        self.t.Branch('mt_tot', self.mt_tot, 'mt_tot/F')
        self.t.Branch('m_vis', self.m_vis, 'm_vis/F')
        self.t.Branch('m_sv', self.m_sv, 'm_sv/F')
        self.t.Branch('mt_sv', self.mt_sv, 'mt_sv/F') 

        # di-lepton variables. 
        self.t.Branch('ll_lmass',  self.ll_lmass,  'll_lmass/F')
        self.t.Branch('mll',       self.mll,       'mll/F')   
        self.t.Branch('ll_pt_p',   self.ll_pt_p,   'll_pt_p/F')   
        self.t.Branch('ll_phi_p',  self.ll_phi_p,  'll_phi_p/F')  
        self.t.Branch('ll_eta_p',  self.ll_eta_p,  'll_eta_p/F')    
        self.t.Branch('ll_pt_m',   self.ll_pt_m,   'll_pt_m/F')      
        self.t.Branch('ll_phi_m',  self.ll_phi_m,  'll_phi_m/F')    
        self.t.Branch('ll_eta_m',  self.ll_eta_m,  'll_eta_m/F')      
        
        # MET variables
        self.t.Branch('met', self.met, 'met/F')
        self.t.Branch('metphi', self.metphi, 'metphi/F')
        self.t.Branch('puppimet', self.puppimet, 'puppimet/F')
        self.t.Branch('puppimetphi', self.puppimetphi, 'puppimetphi/F')
        self.t.Branch('metcov00', self.metcov00, 'metcov00/F')
        self.t.Branch('metcov01', self.metcov01, 'metcov01/F')
        self.t.Branch('metcov10', self.metcov10, 'metcov10/F')
        self.t.Branch('metcov11', self.metcov11, 'metcov11/F')

        # jet variables
        self.t.Branch('njetspt20', self.njetspt20, 'njetspt20/F') 
        self.t.Branch('njets', self.njets, 'njets/F')
        self.t.Branch('nbtag', self.nbtag, 'nbtag/F')

        self.t.Branch('jpt_1', self.jpt_1, 'jpt_1/F' )
        self.t.Branch('jeta_1', self.jeta_1, 'jeta_1/F' ) 
        self.t.Branch('jphi_1', self.jphi_1, 'jphi_1/F' )
        self.t.Branch('jcsv_1', self.jcsv_1, 'jcsv_1/F' )
        self.t.Branch('jpt_2', self.jpt_2, 'jpt_2/F' )
        self.t.Branch('jeta_2', self.jeta_2, 'jeta_2/F' ) 
        self.t.Branch('jphi_2', self.jphi_2, 'jphi_2/F' )
        self.t.Branch('jcsv_2', self.jcsv_2, 'jcsv_2/F' )

        self.t.Branch('bpt_1', self.bpt_1, 'bpt_1/F' )
        self.t.Branch('beta_1', self.beta_1, 'beta_1/F' ) 
        self.t.Branch('bphi_1', self.bphi_1, 'bphi_1/F' )
        self.t.Branch('bcsv_1', self.bcsv_1, 'bcsv_1/F' )
        self.t.Branch('bpt_2', self.bpt_2, 'bpt_2/F' )
        self.t.Branch('beta_2', self.beta_2, 'beta_2/F' ) 
        self.t.Branch('bphi_2', self.bphi_2, 'bphi_2/F' )
        self.t.Branch('bcsv_2', self.bcsv_2, 'bcsv_2/F' )

    def getAntiEle(self,entry,j,bitPos) :
        if ord(entry.Tau_idAntiEle[j]) & bitPos > 0 : return 1.
        else : return 0.

    def getAntiMu(self,entry,j,bitPos) :
        if ord(entry.Tau_idAntiMu[j]) & bitPos > 0 : return 1.
        else : return 0.

    def get_mt(self,METtype,entry,tau) :
        if METtype == 'MVAMet' :
            # temporary choice 
            dphi = tau.Phi() - entry.MET_phi
            return sqrt(2.*tau.Pt()*entry.MET_pt*(1. - cos(dphi)))
        elif METtype == 'PFMet' :
            dphi = tau.Phi() - entry.MET_phi
            return sqrt(2.*tau.Pt()*entry.MET_pt*(1. - cos(dphi)))
        elif METtype == 'PUPPIMet' :
            dphi = tau.Phi() - entry.PuppiMET_phi
            return sqrt(2.*tau.Pt()*entry.PuppiMET_pt*(1. - cos(dphi)))
        else :
            print("Invalid METtype={0:s} in outTuple.get_mt().   Exiting".format(METtype))

    def getPt_tt(self,entry,tau1,tau2) :
        ptMiss = TLorentzVector() 
        ptMiss.SetPtEtaPhiM(entry.MET_pt,0.,entry.MET_phi,0.)
        return (tau1+tau2+ptMiss).Pt()

    def getMt_tot(self,entry,tau1,tau2) :
        pt1, pt2, met = tau1.Pt(), tau2.Pt(), entry.MET_pt
        phi1, phi2, metphi = tau1.Phi(), tau2.Phi(), entry.MET_phi
        arg = 2.*(pt1*met*(1. - cos(phi1-metphi)) + pt2*met*(1. - cos(phi2-metphi)) + pt1*pt2*(1. - cos(phi2-phi1)))
        return sqrt(arg)

    def getM_vis(self,entry,tau1,tau2) :
        return (tau1+tau2).M()

    def getJets(self,entry,tau1,tau2) :
        nJet30, jetList, bJetList = 0, [], []
        phi2_1, eta2_1 = tau1.Phi(), tau1.Eta() 
        phi2_2, eta2_2 = tau2.Phi(), tau2.Eta() 
        for j in range(entry.nJet) :
            if entry.Jet_pt[j] < 20. : break
            if abs(entry.Jet_eta[j]) > 4.7 : continue
            phi1, eta1 = entry.Jet_phi[j], entry.Jet_eta[j]
            dPhi = min(abs(phi2_1-phi1),2.*pi-abs(phi2_1-phi1))
            DR = sqrt(dPhi**2 + (eta2_1-eta1)**2)
            dPhi = min(abs(phi2_2-phi1),2.*pi-abs(phi2_2-phi1))
            DR = min(DR,sqrt(dPhi**2 + (eta2_2-eta1)**2))
            if DR < 0.5 : continue
            if True  and abs(entry.Jet_eta[j]) < 2.5 and entry.Jet_btagDeepB[j] > 0.4941 : bJetList.append(j)
            #if True and abs(entry.Jet_eta[j]) < 2.4 and entry.Jet_btagCSVV2[j] > 0.800 and entry.Jet_pt[j] > 30. : bJetList.append(j)
            if entry.Jet_jetId[j] & 2 == 0 : continue
            if entry.Jet_pt[j] < 30. : continue
            nJet30 += 1
            jetList.append(j) 

        return nJet30, jetList, bJetList 

    def runSVFit(self, entry, channel, jt1, jt2, tau1, tau2 ) :
                      
        measuredMETx = entry.MET_pt*cos(entry.MET_phi)
        measuredMETy = entry.MET_pt*sin(entry.MET_phi)

        #define MET covariance
        covMET = ROOT.TMatrixD(2,2)
        #covMET[0][0] = entry.MET_covXX
        #covMET[1][0] = entry.MET_covXY
        #covMET[0][1] = entry.MET_covXY
        #covMET[1][1] = entry.MET_covYY
        covMET[0][0] = 772.
        covMET[1][0] = 0.
        covMET[0][1] = 0.
        covMET[1][1] = 772.

        #self.kUndefinedDecayType, self.kTauToHadDecay,  self.kTauToElecDecay, self.kTauToMuDecay = 0, 1, 2, 3

        if channel == 'et' :
            measTau1 = ROOT.MeasuredTauLepton(self.kTauToElecDecay, tau1.Pt(), tau1.Eta(), tau1.Phi(), 0.000511) 
        elif channel == 'mt' :
            measTau1 = ROOT.MeasuredTauLepton(self.kTauToMuDecay, tau1.Pt(), tau1.Eta(), tau1.Phi(), 0.106) 
        elif channel == 'tt' :
            measTau1 = ROOT.MeasuredTauLepton(self.kTauToHadDecay, tau1.Pt(), tau1.Eta(), tau1.Phi(), entry.Tau_mass[jt1])
                        
        measTau2 = ROOT.MeasuredTauLepton(self.kTauToHadDecay, tau2.Pt(), tau2.Eta(), tau2.Phi(), entry.Tau_mass[jt2])

        VectorOfTaus = ROOT.std.vector('MeasuredTauLepton')
        instance = VectorOfTaus()
        instance.push_back(measTau1)
        instance.push_back(measTau2)

        FMTT = ROOT.FastMTT()
        FMTT.run(instance, measuredMETx, measuredMETy, covMET)
        ttP4 = FMTT.getBestP4()
        return ttP4.M(), ttP4.Mt() 
    
        
    def Fill(self,entry,SVFit,channel,jt1,jt2,LepP,LepM) :

        # jt1 and jt2 point to the selected tau candidates according to the table below.
        # if e.g., channel = 'et', the jt1 points to the electron list and jt2 points to the tau list.
        # LepP and LepM are TLorentz vectors for the positive and negative members of the dilepton pair
        
        self.entries += 1

        self.run[0]  = entry.run
        self.lumi[0] = entry.luminosityBlock 
        self.evt[0]  = entry.event

        self.againstElectronVLooseMVA6_1[0] = -1.
        self.againstElectronLooseMVA6_1[0]  = -1.
        self.againstElectronMediumMVA6_1[0] = -1.
        self.againstElectronTightMVA6_1[0]  = -1.
        self.againstElectronVTightMVA6_1[0] = -1.
        self.againstMuonLoose3_1[0] = -1.
        self.againstMuonTight3_1[0] = -1. 
        self.byIsolationMVA3oldDMwLTraw_1[0] = -1.

        tauMass = 1.7768 
        tau1, tau2 = TLorentzVector(), TLorentzVector()
        if channel == 'et' :
            self.pt_1[0] = entry.Electron_pt[jt1]
            self.phi_1[0] = entry.Electron_phi[jt1]
            self.eta_1[0] = entry.Electron_eta[jt1]
            self.m_1[0] = entry.Electron_mass[jt1]
            self.q_1[0] = entry.Electron_charge[jt1]
            self.d0_1[0] = entry.Electron_dxy[jt1]
            self.dZ_1[0] = entry.Electron_dz[jt1]
            #self.iso_1[0] = entry.Electron_mvaFall17noIso[jt1]
            self.iso_1[0] = entry.Electron_mvaFall17V2noIso[jt1]
            tau1.SetPtEtaPhiM(entry.Electron_pt[jt1],entry.Electron_eta[jt1], entry.Electron_phi[jt1], tauMass)
            tau2.SetPtEtaPhiM(entry.Tau_pt[jt2],entry.Tau_eta[jt2],entry.Tau_phi[jt2],tauMass) 
        elif channel == 'mt' :
            self.pt_1[0] = entry.Muon_pt[jt1]
            self.phi_1[0] = entry.Muon_phi[jt1]
            self.eta_1[0] = entry.Muon_eta[jt1]
            self.m_1[0] = entry.Muon_mass[jt1]
            self.q_1[0] = entry.Muon_charge[jt1]
            self.d0_1[0] = entry.Muon_dxy[jt1]
            self.dZ_1[0] = entry.Muon_dz[jt1]
            self.iso_1[0] = entry.Muon_pfRelIso04_all[jt1]
            tau1.SetPtEtaPhiM(entry.Muon_pt[jt1],entry.Muon_eta[jt1], entry.Muon_phi[jt1], tauMass)
            tau2.SetPtEtaPhiM(entry.Tau_pt[jt2],entry.Tau_eta[jt2],entry.Tau_phi[jt2],tauMass) 
        elif channel == 'tt' :
            #self.gen_match_1[0] = entry.Tau_genPartIdx[jt1]
            self.gen_match_1[0] = ord(entry.Tau_genPartFlav[jt1])   
	    self.pt_1[0] = entry.Tau_pt[jt1]
            self.phi_1[0] = entry.Tau_phi[jt1]
            self.eta_1[0] = entry.Tau_eta[jt1]
            self.m_1[0] = entry.Tau_mass[jt1]
            self.q_1[0] = entry.Tau_charge[jt1]
            self.d0_1[0] = entry.Tau_dxy[jt1]
            self.dZ_1[0] = entry.Tau_dz[jt1]
            self.iso_1[0] = entry.Tau_rawMVAoldDM2017v2[jt1]
            self.againstElectronVLooseMVA6_1[0] = self.getAntiEle(entry,jt1,1)
            self.againstElectronLooseMVA6_1[0]  = self.getAntiEle(entry,jt1,2)
            self.againstElectronMediumMVA6_1[0] = self.getAntiEle(entry,jt1,4)
            self.againstElectronTightMVA6_1[0]  = self.getAntiEle(entry,jt1,8)
            self.againstElectronVTightMVA6_1[0] = self.getAntiEle(entry,jt1,16)
            self.againstMuonLoose3_1[0] = self.getAntiMu(entry,jt1,1)
            self.againstMuonTight3_1[0] = self.getAntiMu(entry,jt1,2)                            
            #self.byIsolationMVA3oldDMwLTraw_1[0] = float(ord(entry.Tau_idMVAoldDMdR032017v2[jt1]))  # check this
            self.byIsolationMVA3oldDMwLTraw_1[0] = 0.    # does not seem to exist in nanoAOD
            tau1.SetPtEtaPhiM(entry.Tau_pt[jt1],entry.Tau_eta[jt1], entry.Tau_phi[jt1], tauMass)
            tau2.SetPtEtaPhiM(entry.Tau_pt[jt2],entry.Tau_eta[jt2],entry.Tau_phi[jt2],tauMass) 
        else :
            print("Invalid channel={0:s} in outTuple(). Exiting.".format(channel))
            exit()

        self.mt_1[0] = self.get_mt('MVAMet',entry,tau1)
        self.pfmt_1[0] = self.get_mt('PFMet',entry,tau1)
        self.puppimt_1[0] = self.get_mt('PUPPIMet',entry,tau1)

        self.trigweight_1[0] =  -999.   # requires sf need help from Sam on these
        self.idisoweight_1[0] = -999.   # requires sf need help from Sam on these 

	self.pt_2[0] = entry.Tau_pt[jt2]
        self.phi_2[0] = entry.Tau_phi[jt2]
        self.eta_2[0] = entry.Tau_eta[jt2]
        self.m_2[0] = entry.Tau_mass[jt2]
        self.q_2[0] = entry.Tau_charge[jt2]
        self.d0_2[0] = entry.Tau_dxy[jt2]
        self.dZ_2[0] = entry.Tau_dz[jt2]
        phi, pt  = entry.Tau_phi[jt2], entry.Tau_pt[jt2]
        self.mt_2[0] = self.get_mt('MVAMet',entry,tau2) 
        self.pfmt_2[0] = self.get_mt('PFMet',entry,tau2)
        self.puppimt_2[0] = self.get_mt('PUPPIMet',entry,tau2) 
        self.iso_2[0] = entry.Tau_rawMVAoldDM2017v2[jt2]     # check this
        #self.iso_2[0] = entry.Tau_rawMVAoldDM2017v1[jt2]
        self.gen_match_2[0] = entry.Tau_genPartIdx[jt2]   

        self.againstElectronVLooseMVA6_2[0] = self.getAntiEle(entry,jt2,1)
        self.againstElectronLooseMVA6_2[0]  = self.getAntiEle(entry,jt2,2)
        self.againstElectronMediumMVA6_2[0] = self.getAntiEle(entry,jt2,4)
        self.againstElectronTightMVA6_2[0]  = self.getAntiEle(entry,jt2,8)
        self.againstElectronVTightMVA6_2[0] = self.getAntiEle(entry,jt2,16)
        self.againstMuonLoose3_2[0] = self.getAntiMu(entry,jt2,1)
        self.againstMuonTight3_2[0] = self.getAntiMu(entry,jt2,2)

        self.byIsolationMVA3oldDMwLTraw_2[0] = float(ord(entry.Tau_idMVAoldDMdR032017v2[jt2]))  # check this	
        self.trigweight_2[0] = -999.    # requires sf need help from Sam on these
        self.idisoweight_2[0] = -999.   # requires sf need help from Sam on these

        # di-tau variables
        self.pt_tt[0] = self.getPt_tt(entry,tau1,tau2)
        self.mt_tot[0] = self.getMt_tot(entry,tau1,tau2)
        self.m_vis[0] = self.getM_vis(entry,tau1,tau2)
        if SVFit :
            fastMTTmass, fastMTTtransverseMass = self.runSVFit(entry, channel, jt1, jt2, tau1, tau2) 
        else :
            fastMTTmass, fastMTTtransverseMass = -999., -999.
            
        self.m_sv[0] = fastMTTmass 
        self.mt_sv[0] = fastMTTtransverseMass  

        # di-lepton variables.   _p and _m refer to plus and minus charge

        self.ll_lmass[0]  = LepP.M() 
        self.mll[0]       = (LepP + LepM).M() 
        self.ll_pt_p[0]   = LepP.Pt()
        self.ll_phi_p[0]  = LepP.Phi()
        self.ll_eta_p[0]  = LepP.Eta()
        self.ll_pt_m[0]   = LepM.Pt()
        self.ll_phi_m[0]  = LepM.Phi()
        self.ll_eta_m[0]  = LepM.Eta()
        
        # MET variables
        self.met[0] = entry.MET_pt    
        self.metphi[0] = entry.MET_phi
        self.puppimet[0] = entry.PuppiMET_pt
        self.puppimetphi[0] = entry.PuppiMET_phi
        
        #self.metcov00[0] = entry.MET_covXX
        #self.metcov01[0] = entry.MET_covXY
        #self.metcov10[0] = entry.MET_covXY	
        #self.metcov11[0] = entry.MET_covYY

        self.metcov00[0] = 772.
        self.metcov01[0] = 0. 
        self.metcov10[0] = 0. 
        self.metcov11[0] = 772.

        # jet variables
        nJet30, jetList, bJetList = self.getJets(entry,tau1,tau2) 
        self.njetspt20[0] = len(jetList)
        self.njets[0] = nJet30
        self.nbtag[0] = len(bJetList)
        
        self.jpt_1[0], self.jeta_1[0], self.jphi_1[0], self.jcsv_1[0] = -9.99, -9.99, -9.99, -9.99 
        if len(jetList) > 0 :
            jj1 = jetList[0]
            self.jpt_1[0] = entry.Jet_pt[jj1]
            self.jeta_1[0] = entry.Jet_eta[jj1]
            self.jphi_1[0] = entry.Jet_phi[jj1]
            self.jcsv_1[0] = entry.Jet_btagDeepB[jj1]

        self.jpt_2[0], self.jeta_2[0], self.jphi_2[0], self.jcsv_2[0] = -9.99, -9.99, -9.99, -9.99 
        if len(jetList) > 1 :
            jj2 = jetList[1] 
            self.jpt_2[0] = entry.Jet_pt[jj2]
            self.jeta_2[0] = entry.Jet_eta[jj2]
            self.jphi_2[0] = entry.Jet_phi[jj2]
            self.jcsv_2[0] = entry.Jet_btagDeepB[jj2]

        self.bpt_1[0], self.beta_1[0], self.bphi_1[0], self.bcsv_1[0] = -9.99, -9.99, -9.99, -9.99
        if len(bJetList) > 0 :
            jbj1 = bJetList[0]
            self.bpt_1[0] = entry.Jet_pt[jbj1]
            self.beta_1[0] = entry.Jet_eta[jbj1]
            self.bphi_1[0] = entry.Jet_phi[jbj1]
            self.bcsv_1[0] = entry.Jet_btagDeepB[jbj1] 

        self.bpt_2[0], self.beta_2[0], self.bphi_2[0], self.bcsv_2[0] = -9.99, -9.99, -9.99, -9.99
        if len(bJetList) > 1 :
            jbj2 = bJetList[1] 
            self.bpt_2[0] = entry.Jet_pt[jbj2]
            self.beta_2[0] = entry.Jet_eta[jbj2]
            self.bphi_2[0] = entry.Jet_phi[jbj2]
            self.bcsv_2[0] = entry.Jet_btagDeepB[jbj2]

        self.t.Fill()
        #self.weight[0] = 1.
        return

    def setWeight(self,weight) :
        self.weight[0] = weight
        print("outTuple.setWeight() weight={0:f}".format(weight))
        return

    def writeTree(self) :
        print("In outTuple.writeTree() entries={0:d}".format(self.entries)) 
        self.f.Write()
        self.f.Close()
        return

    
