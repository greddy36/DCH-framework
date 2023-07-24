# !/usr/bin/env python


import io
import yaml
import subprocess
from ROOT import TLorentzVector
from math import sqrt, sin, cos, pi
import os
import os.path
import sys

__author__ = "Alexis Kalogeropoulos"

basedir=os.getenv("CMSSW_BASE")
print  'THIS IS THE BASEDIR OF CMSSW', basedir

#sys.path.append('./TauPOG')
sys.path.append(str(basedir)+'/src/TauPOG/')
#sys.path.append('./TauPOG')
#sys.path.append('TauPOG')


from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool
from TauPOG.TauIDSFs.TauIDSFTool import TauESTool
from TauPOG.TauIDSFs.TauIDSFTool import TauFESTool

'''
except ImportError : 
    sys.path.append('../TauPOG')
    from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool
    from TauPOG.TauIDSFs.TauIDSFTool import TauESTool
    from TauPOG.TauIDSFs.TauIDSFTool import TauFESTool
'''

class Weights() :
    def __init__(self,year):
        self.campaign = {2016:'2016Legacy', 2017:'2017ReReco', 2018:'2018ReReco'}

	self.weights_muToTauFR={''}
	self.weights_elToTauFR={''}
	self.weights_muTotauES={''}
	self.weights_elTotauES={''}
	self.weights_tauTotauES={''}
	self.weights_muES={''}
	self.weights_electronES={''}
	self.weights_electronES={''}
        self.LeptV = TLorentzVector()
        self.MetV = TLorentzVector()

	self.weights_muES = {'eta0to1p2' : 0.4, 'eta1p2to2p1' : 0.9, 'etagt2p1' : 2.7 }
	self.weights_electronES = {'eta0to1p2' : 1, 'eta1p2to2p1' : 1, 'etagt2p1' : 2 }
	#self.weights_electronES = {'0to1p2' : 2, '1p2to2p1' : 0.9, '2p1to2p4' : 2.7 }

	if year == 2016 : 

	    self.weights_muToTauFR = {'DM1' : 1.38, 'lt0p4' : 0.80, '0p4to0p8' : 0.81, '0p8to1p2' : 0.79, '1p2to1p7' : 0.68, '1p7to2p3' : 1.34 }
	    self.weights_elToTauFR = {'lt1p479_DM0' : 0.80, 'gt1p479_DM0' : 0.72, 'lt1p479_DM1' : 1.14, 'gt1p479_DM1' : 0.64 }


	    self.weights_tauTotauES = {'DM0' :-0.9, 'DM1' : -0.1, 'DM10' : 0.3, 'DM11' : -0.2}
	    self.weights_elTotauES = {'DM0' :-0.5, 'DM1' : 6, 'DM10' : 0, 'DM11' :0}
	    self.weights_muTotauES = {'DM0' :0., 'DM1' : -0.5, 'DM10' : 0, 'DM11' :0}
            self.TESSF={'dir' : '../TauPOG/TauIDSFs/data/', 'fileTES' : 'TauES_eta-dm_DeepTau2017v2p1VSe_2016Legacy.root'}
            self.FESSF={'dir' : '../TauPOG/TauIDSFs/data/', 'fileFES' : 'TauFES_eta-dm_DeepTau2017v2p1VSe_2016Legacy.root'}


	if year == 2017 : 

	    self.weights_muToTauFR = {'DM1' : 0.69, 'lt0p4' : 1.14, '0p4to0p8' : 1., '0p8to1p2' : 0.87, '1p2to1p7' : 0.52, '1p7to2p3' : 1.47 }
	    self.weights_elToTauFR = {'lt1p479_DM0' : 0.98, 'gt1p479_DM0' : 0.80, 'lt1p479_DM1' : 1.07, 'gt1p479_DM1' : 0.64 }


	    self.weights_tauTotauES = {'DM0' :0.4, 'DM1' : 0.2, 'DM10' : 0.1, 'DM11' : -1.3}
	    self.weights_elTotauES = {'DM0' :0.3, 'DM1' : 3.6, 'DM10' : 0, 'DM11' :0}
	    self.weights_muTotauES = {'DM0' :-0.2, 'DM1' : -0.8, 'DM10' : 0, 'DM11' :0}  ##that was different from Cecile
            self.TESSF={'dir' : '../TauPOG/TauIDSFs/data/', 'fileTES' : 'TauES_dm_DeepTau2017v2p1VSjet_2017ReReco.root'}
            self.FESSF={'dir' : '../TauPOG/TauIDSFs/data/', 'fileFES' : 'TauFES_eta-dm_DeepTau2017v2p1VSjet_2017ReReco.root'}

	if year == 2018 : 

	    self.weights_muToTauFR = {'DM1' : 0.55, 'lt0p4' : 1.08, '0p4to0p8' : 0.78, '0p8to1p2' : 0.77, '1p2to1p7' : 0.75, '1p7to2p3' : 2.02 }
	    self.weights_elToTauFR = {'lt1p479_DM0' : 1.09, 'gt1p479_DM0' : 0.80, 'lt1p479_DM1' : 0.85, 'gt1p479_DM1' : 0.49 }

	    self.weights_tauTotauES = {'DM0' :-1.6, 'DM1' : -0.4, 'DM10' : -1.2, 'DM11' : -0.4}
	    self.weights_elTotauES = {'DM0' :-3.2, 'DM1' : 2.6, 'DM10' : 0, 'DM11' :0}
	    self.weights_muTotauES = {'DM0' :-0.2, 'DM1' : -1., 'DM10' : 0, 'DM11' :0}
            self.TESSF={'dir' : '../TauPOG/TauIDSFs/data/', 'fileTES' : 'TauES_dm_DeepTau2017v2p1VSjet_2018ReReco.root'}
            self.FESSF={'dir' : '../TauPOG/TauIDSFs/data/', 'fileFES' : 'TauFES_eta-dm_DeepTau2017v2p1VSe_2018ReReco.root'}


	''' 
	for reco tauh matched to gen tauh at gen level in the format (dm0, dm1, dm10, dm11): for 2016 (-0.9%, -0.1%, +0.3%, -0.2%), for 2017 (+0.4%, +0.2%, +0.1%, -1.3%), for 2018 (-1.6%, -0.4%, -1.2%, -0.4%)

	for reco tauh matched to electrons at gen level in the format (dm0, dm1): for 2016 (-0.5%, +6.0%), for 2017 (+0.3%, +3.6%), for 2018 (-3.2%, +2.6%)
	for reco tauh matched to muons at gen level in the format (dm0, dm1): for 2016 (+0.0%, -0.5%), for 2017 (+0.0%, +0.0%), for 2018 (-0.2%, -1.0%)
	'''
        wpp = 'Medium'
	self.tauSFTool = TauIDSFTool(self.campaign[year],'DeepTau2017v2p1VSjet',wpp)
	self.testool = TauESTool(self.campaign[year],'DeepTau2017v2p1VSjet',self.TESSF['dir'])
	self.festool = TauFESTool(self.campaign[year],'DeepTau2017v2p1VSe',self.FESSF['dir'])



    def transverseVEC(self,vec):
	rvec = TLorentzVector()
	rvec.SetPtEtaPhiM(vec.Pt(), 0, vec.Phi(), 0)
	return rvec
        
    def correctallMET(self,entry, allmets, year, CorList, LeptList):
        
        list_of_arrays=[]
        list_of_arraysPhi=[]
        tryV=TLorentzVector()
        leptV=TLorentzVector()
	for i, v in enumerate(allmets) :

	    if str(year)=='2017' :
		#i_ should be the righ-hand of the branch and should retain the METFixEE2017 if y=2017
		#iMET should appear always at the branch name...
		v = v.replace('MET','METFixEE2017')
	    iMET= v.replace('METFixEE2017','MET')

	    try : 
                j = getattr(entry, "{0:s}".format(str(v)))
                vphi= v.replace('_pt','_phi')
                jphi = getattr(entry, "{0:s}".format(str(vphi)))
                tryV.SetPx(j*cos(jphi))
                tryV.SetPy(j*sin(jphi))
                #if i==0 : print 'Correction LISTTTTTTTTTTTTTTT', CorList, LeptList
                for j,ic in enumerate(LeptList) :
                    mass = entry.Tau_mass[ic]
                    gen_match = ord(entry.Tau_genPartFlav[ic])
                    if entry.Tau_decayMode[ic] == 0 and gen_match==5 : mass = 0.13960
                    leptV.SetPtEtaPhiM(entry.Tau_pt[ic], entry.Tau_eta[ic], entry.Tau_phi[ic], mass)
                    
                    tryV += (leptV-leptV*CorList[j])
                    #if i==0 : print 'correction for lepton', ic,  CorList[j], entry.nTau

                #if entry.event==8904 and 'pt_jerUp' in v : print 'correct pt', tryV.Pt(), entry.event, v, 'from ttree', j
                list_of_arrays.append(tryV.Pt())
                list_of_arraysPhi.append(tryV.Phi())
                #if i==0 or i==24 or i==25 : print 'inside all MET',  i, v, 'nomMET', entry.MET_T1_pt, 'sysMETuncor', j, 'sysMET_cor', tryV.Pt(), entry.event, entry.luminosityBlock, entry.run
            
	    except AttributeError : 
                setattr(entry, "{0:s}".format(str(v)), j)

        return list_of_arrays, list_of_arraysPhi

    def applyES(self,entry, year, systematic,  metPtPhi, allMETs, printOn=False) :
    #def applyES(self,entry, year, systematic,  metPtPhi,  printOn=False) :


        self.LeptV.SetPtEtaPhiM(0,0,0,0)
        self.MetV.SetPtEtaPhiM(0,0,0,0)
        metpt = float(metPtPhi[0])
        metphi = float(metPtPhi[1])
        self.MetV.SetPx(metpt * cos (metphi))
        self.MetV.SetPy(metpt * sin (metphi))
        cor=1.

        sign = 1.

        corList =[]
        leptList=[]
        metlist=[]
        philist=[]
        if 'Down' in systematic : sign = -1.


        # tauES should be applied by default, except if you plan to apply the tauES +/- systematics - make sure you don't apply the tauES twice...
        #if 'Central' in systematic and 'scale_e' not in systematic and 'scale_m_' not in systematic: 
        #if systematic =='Central' : 
        if True : 

	    for j in range(entry.nTau):    
	       
		dm = entry.Tau_decayMode[j]
		dmm='DM{0:d}'.format(dm)  
		pt= entry.Tau_pt[j]
		eta= entry.Tau_eta[j]
		gen_match = ord(entry.Tau_genPartFlav[j]) 


                isDM0 =  '1prong' in systematic and 'zero' not in systematic and dm == 0
                isDM1 =  '1prong' in systematic and 'zero' in systematic and dm == 1
                isDM10 =  '3prong' in systematic and 'zero' not in systematic and dm == 10
                isDM11 =  '3prong' in systematic and 'zero' in systematic and dm == 11

		if dm != 0 and dm != 1 and dm!=10 and dm!=11 : continue

                
		self.LeptV.SetPtEtaPhiM(entry.Tau_pt[j], entry.Tau_eta[j], entry.Tau_phi[j], entry.Tau_mass[j])

		if gen_match == 5 :

		    tes = self.testool.getTES(pt,dm,gen_match)
                    #if systematic == 'Central' :
                    #    metlist ,philist = self.correctallMET(entry,allMETs, self.LeptV,year,tes)
                    corList.append(tes)
                    leptList.append(j)
                    dirr=''
		    if 'Up'  in systematic : dirr = 'Up'
		    if 'Down' in systematic : dirr = 'Down'

                    if isDM0 or isDM1 or isDM10 or isDM11 : tes = self.testool.getTES(pt,dm,gen_match, unc=dirr)

                    self.MetV +=( self.LeptV - self.LeptV*tes)

                    self.LeptV *= tes          
                    #if 'prong' in systematic : print '----------------------------- inside for tau with gen_match=5', j, tes , oldMET, self.MetV.Pt(), oldphi, self.MetV.Phi(), systematic, entry.event, oldpt, self.LeptV.Pt(), oldphii, self.LeptV.Phi()

                    if printOn :print '----------------------------- inside for tau with gen_match=5', tes , self.MetV.Pt(), systematic

                    #print 'taus with gen_match=5', entry.Tau_pt[j], self.LeptV.Pt(),  'tesUp/Down',tes, j, entry.nTau, systematic, entry.event, 'corrected MetV.Pt()', self.MetV.Pt(), 'MEtV.Pt() before this correction', met_uncor, 'fed in', metpt, 'METPx diff', self.MetV.Px(), self.LeptV.Px()

		    if dm == 0 :
			entry.Tau_mass[j] = 0.13960

		#if 'prong' not in systematic :
		if gen_match == 2 or gen_match == 4 :

                    cor=1+self.weights_muTotauES[dmm]*0.01
                    #if systematic == 'Central' :
                    #    metlist,philist = self.correctallMET(entry,allMETs, self.LeptV,year,cor)
                    corList.append(cor)
                    leptList.append(j)

		    self.MetV += ( self.LeptV - self.LeptV *(1 + self.weights_muTotauES[dmm]*0.01))

		    self.LeptV *=  (1 + self.weights_muTotauES[dmm]*0.01)
		    if printOn : print 'will correct for  muon_faking_tau ES', self.weights_muTotauES[dmm]*0.01,  metPtPhi[0], '---->', self.MetV.Pt()
		
		    # leptons faking taus // electron->tau
		if gen_match == 1 or gen_match == 3 :

		    fes = self.festool.getFES(eta,dm,gen_match)
                    #if systematic == 'Central' :
                    #    metlist,philist = self.correctallMET(entry,allMETs, self.LeptV,year,fes)
                    corList.append(fes)
                    leptList.append(j)
		   
		    self.MetV +=( self.LeptV - self.LeptV*fes)

		    self.LeptV *= fes          
		    if printOn : print 'will correct for  electron_faking_tau ES', cor,  metPtPhi[0], '---->',self.MetV.Pt()

		entry.Tau_mass[j] = self.LeptV.M()
		entry.Tau_pt[j] = self.LeptV.Pt()
		entry.Tau_phi[j] = self.LeptV.Phi()
                #entry.Tau_phi[j] = self.LeptV.Phi()

		#print 'returnin LeptV ---------------->pt ', self.LeptV.Pt(),  'met' , self.MetV.Pt(), systematic
	    #print 'returning---------------->', entry.METFixEE2017_pt, self.MetV.Pt(), self.MetV.Phi(), self.LeptV.M(), entry.Tau_mass[j], self.LeptV.Pt(), 
	    #print 'returning---------------->pt ', self.MetV.Pt(), 'phi ', self.MetV.Phi(), 'uncorrected values', entry.MET_pt, entry.MET_phi, entry.MET_T1_pt, entry.MET_T1_phi, systematic



        #### muon_energy_scale , in three eta bins
        if 'scale_m' in systematic : 


	    #self.weights_muES = {'eta0to1p2' : 0.4, 'eta1p2to2p1' : 0.9, 'etagt2p1' : 2.7 }
	    for j in range(entry.nMuon):    

                fact = 1.
		self.LeptV.SetPtEtaPhiM(entry.Muon_pt[j], entry.Muon_eta[j], entry.Muon_phi[j], entry.Muon_mass[j])

                if  'scale_m_etalt1p2' in systematic and abs(entry.Muon_eta[j]) < 1.2 : fact = 1 + sign*self.weights_muES['eta0to1p2']*0.01  
                if  'scale_m_eta1p2to2p1' in systematic and  abs(entry.Muon_eta[j]) > 1.2 and abs(entry.Muon_eta[j]) < 2.1 : fact = 1+ sign * self.weights_muES['eta1p2to2p1']*0.01 
                if  'scale_m_etagt2p1' in systematic and abs(entry.Muon_eta[j]) > 2.1 : fact = 1+ sign * self.weights_muES['etagt2p1']*0.01 


                self.MetV += (self.LeptV - self.LeptV*fact)
                self.LeptV *= fact
                entry.Muon_pt[j] = self.LeptV.Pt()
                entry.Muon_phi[j] = self.LeptV.Phi()
                entry.Muon_mass[j] = self.LeptV.M()
                #entry.Muon_phi[j] = self.LeptV.Phi()
                #cor=fact


        #### electron_energy_scale
        if 'scale_e' in systematic : 
            
            if printOn :
		for j in range(entry.nTau):    
		    print 'taus in electron_e tauPt', entry.Tau_pt[j], j, entry.nTau, 'systematic', systematic, entry.event, 'met in applyES', self.MetV.Pt(), 'met fed in', metpt

	    #self.weights_muES = {'eta0to1p2' : 0.4, 'eta1p2to2p1' : 0.9, 'etagt2p1' : 2.7 }
	    for j in range(entry.nElectron):    


                fact=1.
		self.LeptV.SetPtEtaPhiM(entry.Electron_pt[j], entry.Electron_eta[j], entry.Electron_phi[j], entry.Electron_mass[j])

                if  abs(entry.Electron_eta[j]) < 1.2 : fact = 1 + sign*self.weights_electronES['eta0to1p2']*0.01 
                if  abs(entry.Electron_eta[j]) > 1.2 and abs(entry.Electron_eta[j]) < 2.1 :    fact =   1+ sign * self.weights_electronES['eta1p2to2p1']*0.01 
                if  abs(entry.Electron_eta[j]) > 2.1 : fact = 1+ sign * self.weights_electronES['etagt2p1']*0.01 

                #uncormet = self.MetV.Pt()
                self.MetV += (self.LeptV - self.LeptV*fact)
                self.LeptV *= fact
                entry.Electron_pt[j] = self.LeptV.Pt()
                entry.Electron_phi[j] = self.LeptV.Phi()
                entry.Electron_mass[j] = self.LeptV.M()
                #entry.Electron_phi[j] = self.LeptV.Phi()
                #if printOn : print 'scale_e systematic ', systematic, entry.event, 'uncorrected met', uncormet, 'corrected met', self.MetV.Pt(),  'met fed in', metpt


	if systematic == 'Central' :
	    metlist ,philist = self.correctallMET(entry,allMETs, year,corList, leptList)

	#if systematic == 'Central' : print 'returning---------------->pt ', self.MetV.Pt(), 'phi ', self.MetV.Phi(), 'uncorrected values', entry.MET_T1_pt, entry.MET_T1_phi, systematic,  entry.event
	return self.MetV.Pt(), self.MetV.Phi(), metlist , philist
	#return self.MetV.Pt(), self.MetV.Phi()

