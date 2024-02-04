# !/usr/bin/env python

""" tauFun.py: apply selection sequence to four-lepton final state """

import io
import yaml
import subprocess
from ROOT import TLorentzVector
from math import sqrt, sin, cos, pi
from itertools import combinations
import os
import os.path
import sys
import generalFunctions as GF
sys.path.append('../TauPOG')


# get selections from configZH.yaml:
with io.open('cuts.yaml', 'r') as stream:
    selections = yaml.load(stream)
print "Using selections:\n", selections



def goodTrigger(e, year):
    trig = selections['trig']
    if not (trig['singleLepton'] or trig['doubleLepton']) : return True
    #single mu 2016: HLT IsoMu22 v, HLT IsoMu22 eta2p1 v, HLT IsoTkMu22 v, HLT IsoTkMu22 eta2p1 v and cut pt(mu)>23, eta(mu)<2.1
    #single ele 2016: HLT Ele25 eta2p1 WPTight Gsf v and cut pt(ele)>26, eta(ele)<2.1
    #single mu 2017: HLT IsoMu24 v, HLT IsoMu27 v and cut pt(mu)>25, eta(mu)<2.4
    #single ele 2017: HLT Ele27 WPTight Gsf v, HLT Ele32 WPTight Gsf v, HLT Ele35 WPTight Gsf v and cut pt(ele)>28, eta(ele)<2.1
    #single mu 2018: HLT IsoMu24 v, HLT IsoMu27 v and cut pt(mu)>25, eta(mu)<2.4
    #single ele 2018:  HLT Ele32 WPTight Gsf v, HLT Ele35 WPTight Gsf v and cut pt(ele)>33, eta(ele)<2.1
   
    
    if year == 2016 :
        goodSingle = (e.HLT_IsoMu22 or e.HLT_IsoMu22_eta2p1 or e.HLT_IsoTkMu22 or e.HLT_IsoTkMu22_eta2p1 or e.HLT_Ele25_eta2p1_WPTight_Gsf or e.HLT_Ele27_eta2p1_WPTight_Gsf or e.HLT_IsoMu24 or e.HLT_IsoTkMu24 or e.HLT_IsoMu27)

        goodDouble = (e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ or e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or e.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ )
    elif (year == 2017 or year == 2018) :
        goodSingle = (e.HLT_Ele27_WPTight_Gsf or e.HLT_Ele35_WPTight_Gsf or e.HLT_Ele32_WPTight_Gsf or e.HLT_IsoMu24 or e.HLT_IsoMu27)

        goodDouble = (e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL or e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ  or e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 or e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)
    

    else :
        print("Invalid year={0:d} in goodTrigger()".format(year))
        return False
    
    return (trig['singleLepton'] and goodSingle) or (trig['doubleLepton'] and goodDouble)  

def goodPhotonTrigger(e, year):
   
    goodPhoton =  e.HLT_Photon50_R9Id90_HE10_IsoM or e.HLT_Photon75_R9Id90_HE10_IsoM or e.HLT_Photon90_R9Id90_HE10_IsoM  or e.HLT_Photon120_R9Id90_HE10_IsoM  or e.HLT_Photon165_R9Id90_HE10_IsoM

    return goodPhoton

    

def getTauList(dch, entry, pairList=[],printOn=False, isDCH2=False, signC=0) :
    """ tauFun.getTauList(): return a list of taus that 
                             pass the basic selection cuts               
    """

    if not dch in ['et','mt','tt'] :
        print("Warning: invalid dch={0:s} in tauFun.getTauList()".format(dch))
        exit()

    if printOn : print ' getTauList : will be checking nTau', entry.nTau

    if entry.nTau == 0: 
        if printOn : print ' failed nTau', entry.nTau
        return []

    tauList = []
    tt = selections['tt'] # selections for H->tau(h)+tau(h)
    for j in range(entry.nTau):   
        #print entry.Tau_pt[j] 
        # apply tau(h) selections 
        if printOn : print 'looking for Tau j', j, 'Q', entry.Tau_charge[j]

        if signC !=0 and entry.Tau_charge[j] == signC:
            if printOn : print("        fail same Tau Q {0:d} as lepton Q".format(entry.Tau_charge[j]))
            continue

        if entry.Tau_pt[j] < tt['tau_pt']: 
            if printOn : print("        fail Tau pT {0:f}".format(entry.Tau_pt[j]))
            continue
        if abs(entry.Tau_eta[j]) > tt['tau_eta']: 
       
            if printOn : print("        fail Tau eta {0:f}".format(abs(entry.Tau_eta[j])))
            continue
        if abs(entry.Tau_dz[j]) > tt['tau_dz']: 
            if printOn : print("        fail Tau dZ {0:f}".format(entry.Tau_dz[j]))
            continue
        #if not entry.Tau_idDecayModeNewDMs[j]: 
        #    if printOn : print("        fail Tau decayModeNewDM {0:f}".format(entry.Tau_idDecayModeNewDMs[j]))
        #    continue
	if  entry.Tau_decayMode[j] == 5 or entry.Tau_decayMode[j] == 6 : 
            if printOn : print("        fail Tau decayMode {0:f}".format(entry.Tau_decayMode[j]))
            continue
        if abs(entry.Tau_charge[j]) != 1: 
            if printOn : print("        fail Tau Q {0:d}".format(entry.Tau_charge[j]))
            continue

        
        if tt['tau_vJet'] > 0  and not ord(entry.Tau_idDeepTau2017v2p1VSjet[j]) & tt['tau_vJet'] > 0 :
            if printOn : print("        fail DeepTau vs. Jet={0:d}".format(ord(entry.Tau_idDeepTau2017v2p1VSjet[j])))
            continue
	if tt['tau_vEle'] > 0 and not ord(entry.Tau_idDeepTau2017v2p1VSe[j]) & tt['tau_vEle'] > 0 :
            if printOn : print("        fail DeepTau vs. ele={0:d}".format(ord(entry.Tau_idDeepTau2017v2p1VSe[j])))
            continue
        if tt['tau_vMu'] > 0 and not ord(entry.Tau_idDeepTau2017v2p1VSmu[j]) & tt['tau_vMu'] > 0 :
            if printOn : print("        fail DeepTau vs.  mu={0:d}".format(ord(entry.Tau_idDeepTau2017v2p1VSmu[j])))
            continue

        if isDCH2 :#Set false for dch1 and true for dch2
            eta, phi = entry.Tau_eta[j], entry.Tau_phi[j]
            DR0, DR1 =  lTauDR(eta,phi, pairList[0]), lTauDR(eta,phi,pairList[1])
	    if DR0 < tt['tt_DR'] or DR1 < tt['tt_DR']: 
		if printOn : print("        fail DR0={0:f} or DR1={1:f} ". format(DR0, DR1))
		continue

        if printOn: print "charge of taus ->sign of opposite one ", signC, "charge of Tau", entry.Tau_charge[j], "idx", j
        #print ord(entry.Tau_idDeepTau2017v2p1VSmu[j]), ord(entry.Tau_idDeepTau2017v2p1VSe[j]), ord(entry.Tau_idDeepTau2017v2p1VSjet[j])
       
        tauList.append(j)
    
    if printOn  : print 'returning with tauList from getTauList', tauList
    return tauList

def makeGoodTauList(entry, cutflow, printOn=False) :
    """ tauFun.makeGoodTauList(): return a list of taus that 
                             pass the basic selection cuts               
    """

    if entry.nTau == 0: return []

    tauList = []
    tt = selections['tt'] # selections for H->tau(h)+tau(h)

    #for j in range(entry.nTau):   
    '''
    for reco tauh matched to electrons at gen level in the format (dm0, dm1): for 2016 (-0.5%, +6.0%), for 2017 (+0.3%, +3.6%), for 2018 (-3.2%, +2.6%)
    for reco tauh matched to muons at gen level in the format (dm0, dm1): for 2016 (+0.0%, -0.5%), for 2017 (+0.0%, +0.0%), for 2018 (-0.2%, -1.0%)
    '''
    for j in range(entry.nTau):    
        # apply tau(h) selections 
        cutflow.count('cut0')
        if entry.Tau_pt[j] < tt['tau_pt']:
            cutflow.count('cut1')
            continue
        if abs(entry.Tau_eta[j]) > tt['tau_eta']:
            cutflow.count('cut2')
            continue
        if abs(entry.Tau_dz[j]) > tt['tau_dz']:
            cutflow.count('cut3')
            continue
        #if not entry.Tau_idDecayModeNewDMs[j]: continue
	if  entry.Tau_decayMode[j] == 5 or entry.Tau_decayMode[j] == 6 :
            cutflow.count('cut4')
            continue
        if abs(entry.Tau_charge[j]) != 1:
            cutflow.count('cut5')
            continue

        if tt['tau_vJet'] > 0  and not ord(entry.Tau_idDeepTau2017v2p1VSjet[j]) & tt['tau_vJet'] > 0 :
            cutflow.count('cut6')
            if printOn : print("        fail DeepTau vs. Jet={0:d}".format(ord(entry.Tau_idDeepTau2017v2p1VSjet[j])))
            continue
        if tt['tau_vEle'] > 0 and not ord(entry.Tau_idDeepTau2017v2p1VSe[j]) & tt['tau_vEle'] > 0 :
            cutflow.count('cut7')
            if printOn : print("        fail DeepTau vs. ele={0:d}".format(ord(entry.Tau_idDeepTau2017v2p1VSe[j])))
            continue
        if tt['tau_vMu'] > 0 and not ord(entry.Tau_idDeepTau2017v2p1VSmu[j]) & tt['tau_vMu'] > 0 :
            cutflow.count('cut8')
            if printOn : print("        fail DeepTau vs.  mu={0:d}".format(ord(entry.Tau_idDeepTau2017v2p1VSmu[j])))
            continue
        tauList.append(j)
    return tauList

def makeGoodTauListWjets(channel, entry, muIndex,  printOn=False) :
    """ tauFun.getTauList(): return a list of taus that 
                             pass the basic selection cuts               
    """

    if entry.nTau == 0: return []

    tauList = []

    phi2,eta=-100,-100
    if channel=='mnu' : 
        phi2, eta2 = entry.Muon_phi[muIndex], entry.Muon_eta[muIndex]
    if channel=='enu' : 
        phi2, eta2 = entry.Electron_phi[muIndex], entry.Electron_eta[muIndex]
    for j in range(entry.nTau):    
        # apply tau(h) selections 
        if channel=='mnu' and llDR(entry.Tau_eta[j], entry.Tau_phi[j], eta2, phi2) > 0.2 or int(ord(entry.Tau_idDeepTau2017v2p1VSmu[j])) > 0 :  
            #print 'failed dR or antiMu', llDR(entry.Tau_eta[j], entry.Tau_phi[j], eta2, phi2),  entry.Tau_pt[j], 'antiMu', int(ord(entry.Tau_idDeepTau2017v2p1VSmu[j]))
            tauList.append(j)
        if channel=='enu' and llDR(entry.Tau_eta[j], entry.Tau_phi[j], eta2, phi2) > 0.2 or int(ord(entry.Tau_idDeepTau2017v2p1VSe[j])) > 0 :  
            #print 'failed dR or antiMu', llDR(entry.Tau_eta[j], entry.Tau_phi[j], eta2, phi2),  entry.Tau_pt[j], 'antiMu', int(ord(entry.Tau_idDeepTau2017v2p1VSmu[j]))
            tauList.append(j)

        #if tt['tau_vJet'] > 0  and not ord(entry.Tau_idDeepTau2017v2p1VSjet[j]) & tt['tau_vJet'] > 0 :
        #    if printOn : print("        fail DeepTau vs. Jet={0:d}".format(ord(entry.Tau_idDeepTau2017v2p1VSjet[j])))
        #    continue
	#if tt['tau_vEle'] > 0 and not ord(entry.Tau_idDeepTau2017v2p1VSe[j]) & tt['tau_vEle'] > 0 :
        #    if printOn : print("        fail DeepTau vs. ele={0:d}".format(ord(entry.Tau_idDeepTau2017v2p1VSe[j])))
        #    continue
    return tauList


def tauDR(entry, j1,j2) :
    if j1 == j2 : return 0. 
    phi1, eta1, phi2, eta2 = entry.Tau_phi[j1], entry.Tau_eta[j1], entry.Tau_phi[j2], entry.Tau_eta[j2]
    return sqrt( (phi2-phi1)**2 + (eta2-eta1)**2 )

def llDR(eta1,phi1,eta2,phi2) :
    dPhi = min(abs(phi2-phi1),2.*pi-abs(phi2-phi1))
    return sqrt(dPhi**2 + (eta2-eta1)**2)


def lTauDR(eta1,phi1,Lep) :
    phi2, eta2 = Lep.Phi(), Lep.Eta()
    dPhi = min(abs(phi2-phi1),2.*pi-abs(phi2-phi1))
    return sqrt(dPhi**2 + (eta2-eta1)**2)


def DRobj(eta1,phi1,eta2,phi2) :
    dPhi = min(abs(phi2-phi1),2.*pi-abs(phi2-phi1))
    return sqrt(dPhi**2 + (eta2-eta1)**2)

def mllCut(mll) :
    mllcuts = selections['mll']
    if mll < mllcuts['mll_low'] or mll > mllcuts['mll_high'] : return False
    return True

def getTauPointer(entry, eta1, phi1) :
    # find the j value that most closely matches the specified eta or phi value
    bestMatch, jBest = 999., -1
    for j in range(entry.nTau) :
        eta2, phi2 = entry.Tau_eta[j], entry.Tau_phi[j]
        dPhi = min(abs(phi2-phi1),2.*pi-abs(phi2-phi1))
        DR = sqrt(dPhi**2 + (eta2-eta1)**2)
        if DR < bestMatch : bestMatch, jBest = DR, j
    if bestMatch > 0.1 :
        jBest = -1 
        print("Error in getTauPointer():   No match found eta={0:.3f} phi={1:.3f}".format(eta1,phi1))
    return jBest

 

def comparePairvspT(entry, tauPairList, printOn=False) :
    """ comparevsPt : return the index of the pair with the highest scalar sum 
    """
    SumList=[]
    for i in range(0,len(tauPairList)) :
         
        j1, j2 = tauPairList[i][0], tauPairList[i][1] # look at leading pt tau in each pair
        if printOn : print 'appending now', entry.Tau_pt[j1] + entry.Tau_pt[j2], j1, j2, tauPairList[i]
        SumList.append(entry.Tau_pt[j1] + entry.Tau_pt[j2])

    
    maxI=SumList.index(max(SumList))

    return maxI
       
def comparePair(entry, pair1, pair2) :
    """ tauFun.comparePair.py: return true if pair2 is 
                               better than pair1 
    """
    
    j1, j2 = pair1[0], pair2[0] # look at leading pt tau in each pair
    j3, j4 = pair1[1], pair2[1] # look at leading pt tau in each pair
    if entry.Tau_rawDeepTau2017v2p1VSjet[j2] > entry.Tau_rawDeepTau2017v2p1VSjet[j1] :
        return True
    if  entry.Tau_rawDeepTau2017v2p1VSjet[j2] < entry.Tau_rawDeepTau2017v2p1VSjet[j1] :
        return False
    if  entry.Tau_rawDeepTau2017v2p1VSjet[j2] == entry.Tau_rawDeepTau2017v2p1VSjet[j1] :
        if entry.Tau_pt[j2] > entry.Tau_pt[j1] :
            return True
        if entry.Tau_pt[j2] < entry.Tau_pt[j1] :
            return False
        if entry.Tau_pt[j2] == entry.Tau_pt[j1] :
            if  entry.Tau_rawDeepTau2017v2p1VSjet[j4] > entry.Tau_rawDeepTau2017v2p1VSjet[j3] : 
                return True
            if  entry.Tau_rawDeepTau2017v2p1VSjet[j4] < entry.Tau_rawDeepTau2017v2p1VSjet[j3] : 
                return False
            if  entry.Tau_rawDeepTau2017v2p1VSjet[j4] == entry.Tau_rawDeepTau2017v2p1VSjet[j3] : 
                if entry.Tau_pt[j4] > entry.Tau_pt[j3] : return True
                if entry.Tau_pt[j4] < entry.Tau_pt[j3] : return False
                if entry.Tau_pt[j4] == entry.Tau_pt[j3] : return False

    # do it once more swapping the tau pairs
    j1, j2 = pair1[1], pair2[1] # look at leading pt tau in each pair
    j3, j4 = pair1[0], pair2[0] # look at leading pt tau in each pair
    if entry.Tau_rawDeepTau2017v2p1VSjet[j2] > entry.Tau_rawDeepTau2017v2p1VSjet[j1] :
        return True
    if  entry.Tau_rawDeepTau2017v2p1VSjet[j2] < entry.Tau_rawDeepTau2017v2p1VSjet[j1] :
        return False
    if  entry.Tau_rawDeepTau2017v2p1VSjet[j2] == entry.Tau_rawDeepTau2017v2p1VSjet[j1] :
        if entry.Tau_pt[j2] > entry.Tau_pt[j1] :
            return True
        if entry.Tau_pt[j2] < entry.Tau_pt[j1] :
            return False
        if entry.Tau_pt[j2] == entry.Tau_pt[j1] :
            if  entry.Tau_rawDeepTau2017v2p1VSjet[j4] > entry.Tau_rawDeepTau2017v2p1VSjet[j3] : 
                return True
            if  entry.Tau_rawDeepTau2017v2p1VSjet[j4] < entry.Tau_rawDeepTau2017v2p1VSjet[j3] : 
                return False
            if  entry.Tau_rawDeepTau2017v2p1VSjet[j4] == entry.Tau_rawDeepTau2017v2p1VSjet[j3] : 
                if entry.Tau_pt[j4] > entry.Tau_pt[j3] : return True
                if entry.Tau_pt[j4] < entry.Tau_pt[j3] : return False
                if entry.Tau_pt[j4] == entry.Tau_pt[j3] : return False




def getBestTauPair(dch, entry, tauList,printOn=False) :
    """ tauFun.getBestTauPair(): return two taus that 
                                 best represent H->tt
    """ 

    if dch != 'tt' : 
        if printOn : print("Invalid dch={0:s} in tauFun.getBestTauPair()".format(dch))
        exit()

    if len(tauList) < 2: 
        if printOn : print 'failed to find good TauTau Pair'
        if printOn : print("Entering getTauPairs failing nTauList={0:s}".format(tauList))
        return [] 
    
    # form all possible pairs that satisfy DR requirement
    tauPairList = []
    tt = selections['tt'] # selections for H->(tau_h)(tau_h)
    for i in range(len(tauList)) :
        idx_tau1 = tauList[i]
        for j in range(len(tauList)) :
            if i <= j: continue
            idx_tau2 = tauList[j]
            if not entry.Tau_charge[idx_tau1] == entry.Tau_charge[idx_tau2]: continue
            if tauDR(entry, idx_tau1, idx_tau2) < tt['tt_DR'] : 
                tdR= tauDR(entry, idx_tau1, idx_tau2) 
                if printOn : print "failed tDR=",tdR, "tt_DR =", tt['tt_DR']
                continue
            tauPairList.append([idx_tau1, idx_tau2])
            #if idx_tau1 not in tauPairList : tauPairList.append(idx_tau1)
            #if idx_tau2 not in tauPairList : tauPairList.append(idx_tau2)


    if len(tauPairList) == 0 : 
        if printOn : print 'fail tauPairList =0'
        return []
    if len(tauPairList) >= 1 :
        '''for i in range(len(tauPairList)-1, 0, -1) :
            if comparePair(entry, tauPairList[i], tauPairList[i-1]) :
               tauPairList[i-1], tauPairList[i] = tauPairList[i], tauPairList[i-1]
        bestPair = tauPairList[0]#best pair is the first element
        '''
        bestPair=tauPairList[-1]
        for i in range(0,len(tauPairList)-1):
            if not comparePair(entry, tauPairList[i], bestPair) :
               bestPair=tauPairList[i]

    '''
    c = set(combinations(tauPairList, 2))
    if printOn : print 'these are the combinations', c, 'from', tauPairList
    tauPairList=list(c)

    if printOn : print 'Do I need to do the comparePair ? ', len(tauPairList),  tauPairList
    
    maxI= comparePairvspT(entry, tauPairList, printOn)
    tauPairList=tauPairList[maxI]
    if printOn : print 'this is the list of sums', maxI, tauPairList, tauPairList[0], tauPairList[1]

    idx_tau1, idx_tau2 = tauPairList[0], tauPairList[1]
    if entry.Tau_pt[idx_tau2] > entry.Tau_pt[idx_tau1] : 
        if printOn : print 'Sorting tt tauPair', entry.Tau_pt[idx_tau2], entry.Tau_pt[idx_tau1], idx_tau2, idx_tau1
        tauPairList = {}
        tauPairList[0] = idx_tau2
        tauPairList[1] = idx_tau1    
    if printOn : print 'returning tt tauPairList', tauPairList
    '''
    return bestPair

def checkOverlapMuon(entry,i,j,checkDR=False,printOn=False) :

    overlapM = False
    
    mt = selections['mmoverlap'] 
    if mt['mu_type'] and (entry.Muon_isGlobal[i] or entry.Muon_isTracker[i]) : 
	if mt['mu_ID'] and entry.Muon_tightId[i] : 
	    if abs(entry.Muon_dxy[i]) < mt['mu_dxy'] :
		if abs(entry.Muon_dz[i]) < mt['mu_dz']:
		    if entry.Muon_pt[i] > mt['mu_pt']:
			if abs(entry.Muon_eta[i]) < mt['mu_eta']:
			    if  mt['mu_iso_f'] and entry.Muon_pfRelIso04_all[i] < mt['mu_iso']:

				if mt['mu_type'] and (entry.Muon_isGlobal[j] or entry.Muon_isTracker[j]) : 
				    if mt['mu_ID'] and  entry.Muon_tightId[j] : 
					if abs(entry.Muon_dxy[j]) < mt['mu_dxy'] :
					    if abs(entry.Muon_dz[j]) < mt['mu_dz']:
						if entry.Muon_pt[j] > mt['mu_pt']:
						    if abs(entry.Muon_eta[j]) < mt['mu_eta']:
							if  mt['mu_iso_f'] and entry.Muon_pfRelIso04_all[j] < mt['mu_iso']:
							    if printOn: print 'checking overlap....', i, j, DRobj(entry.Muon_eta[i],entry.Muon_phi[i], entry.Muon_eta[j],entry.Muon_phi[j]), entry.Muon_pt[i], entry.Muon_eta[i],entry.Muon_phi[i], 'j-->', entry.Muon_pt[j], entry.Muon_eta[j],entry.Muon_phi[j]
							    if checkDR  :
                                                                if DRobj(entry.Muon_eta[i],entry.Muon_phi[i], entry.Muon_eta[j],entry.Muon_phi[j])  < mt['ll_DR'] : overlapM = True
                                                            else : overlapM = True
                                                            
    return overlapM

def getMuTauPairs(entry,dch='mt',pairList=[],printOn=False,isDCH2=False,signC=0) :
    """  tauFun.getMuTauPairs.py: return list of acceptable pairs
                                 of muons and taus 
    """
   

    if entry.nMuon < dch.count('m') or entry.nTau < dch.count('t'):
        #if printOn : print("Entering getMuTauPairs failing nMuon={0:d} nTau={1:d} lumi={2:s} run={3:s} event={4:s}".format(entry.nMuon,entry.nTau, str(entry.luminosityBlock), str(entry.run), str(entry.event)))
        return []

    muTauPairs = []
    mt = selections['mt'] # H->tau(mu)+tau(h) selections
    #if printOn : print("Entering tauFun.getMuTauPairs() nMuon={0:d} nTau={1:d}".format(entry.nMuon,entry.nTau))
    #printOn=True
    if printOn : print("Entering getMuTauPairs some info nMuon={0:d} nTau={1:d} lumi={2:s} run={3:s} event={4:s}".format(entry.nMuon,entry.nTau, str(entry.luminosityBlock), str(entry.run), str(entry.event)))
    leptOverlap = False
    if dch==0:
    #if entry.nMuon>1 and len(pairList)>0:#dropping too many events
    #if (dch=='mmmt' and entry.nMuon>3) or (dch=='eemt' and entry.nMuon>1):
	for i in range(entry.nMuon):
	
	    if printOn : print 'checking pT with findZ nMuon', entry.nMuon, 'i', i, entry.run, entry.luminosityBlock, entry.event, 'ZpT', pairList[0].Pt(), pairList[1].Pt(), 'vs', entry.Muon_pt[i] ,  pairList[0].Pt()-entry.Muon_pt[i], pairList[1].Pt()-entry.Muon_pt[i], DRobj(entry.Muon_eta[i],entry.Muon_phi[i], pairList[0].Eta(), pairList[0].Phi()), DRobj(entry.Muon_eta[i],entry.Muon_phi[i], pairList[1].Eta(), pairList[1].Phi())

	    if DRobj(entry.Muon_eta[i],entry.Muon_phi[i], pairList[0].Eta(), pairList[0].Phi())<0.1 or DRobj(entry.Muon_eta[i],entry.Muon_phi[i], pairList[1].Eta(), pairList[1].Phi())<0.1 :
                #print 'found a match on findZ', i, entry.Muon_pt[i]

                continue# make sure that you don't consider the findZ leptons
	    for j in range(i+1, entry.nMuon):
	        #if printOn : print 'checking in loop pT with findZ nMuon', entry.nMuon, 'j', j, entry.run, entry.luminosityBlock, entry.event, pairList[0].Pt(), pairList[1].Pt(), 'vs', entry.Muon_pt[j] 
	        if DRobj(entry.Muon_eta[j],entry.Muon_phi[j], pairList[0].Eta(), pairList[0].Phi())<0.1 or DRobj(entry.Muon_eta[j],entry.Muon_phi[j], pairList[1].Eta(), pairList[1].Phi())<0.1 :
                    #print 'found a match on findZ j now ', j, entry.Muon_pt[i]
                    continue

		if checkOverlapMuon(entry,i,j) : 
		    if printOn : print 'Muons ', i, j, 'overlaping ---->', entry.Muon_pt[i], entry.Muon_pt[j], DRobj(entry.Muon_eta[i],entry.Muon_phi[i], entry.Muon_eta[j],entry.Muon_phi[j]), entry.run, entry.luminosityBlock, entry.event
		    return []

    for i in range(entry.nMuon):

        # apply tau(mu) selections
        if mt['mu_type']:
            if not (entry.Muon_isGlobal[i] or entry.Muon_isTracker[i]) :
                if printOn : print("    fail mu_type Global or Tracker={0}".format(entry.Muon_isGlobal[i]))
                continue
        if mt['mu_ID']:
            #if not (entry.Muon_looseId[i] or entry.Muon_mediumId[i] or entry.Muon_tightId[i]) :
            if not (entry.Muon_tightId[i]):
                if printOn : print("    fail mu_ID mediumId={0}".format(entry.Muon_tightId[i]))
                continue
        if abs(entry.Muon_dxy[i]) > mt['mu_dxy']:
            if printOn : print("    fail mu_dxy={0:f}".format(entry.Muon_dxy[i]))
            continue
        if abs(entry.Muon_dz[i]) > mt['mu_dz']:
            if printOn : print("    fail mu_dz={0:f}".format(entry.Muon_dz[i]))
            continue
        mu_eta, mu_phi = entry.Muon_eta[i], entry.Muon_phi[i] 
        if entry.Muon_pt[i] < mt['mu_pt']:
            if printOn : print("    fail mu_pt={0:f}".format(entry.Muon_pt[i]))
            continue
        if abs(mu_eta) > mt['mu_eta']:
            if printOn : print("    fail mu_eta={0:f}".format(entry.Muon_eta[i]))
            continue 
        if  mt['mu_iso_f'] and entry.Muon_pfRelIso04_all[i] > mt['mu_iso']:
            if printOn : print("    fail mu_iso={0:f}".format(entry.Muon_pfRelIso04_all[i]))
            continue
        
        if isDCH2:
            if signC!=0 and entry.Muon_charge[i] == signC : continue
            DR0 = lTauDR(mu_eta,mu_phi,pairList[0]) # l1 vs. tau(mu)
            DR1 = lTauDR(mu_eta,mu_phi,pairList[1]) # l2 vs. tau(mu)
            if printOn : print 'some info----------------', i, entry.Muon_pt[i], entry.Muon_eta[i], entry.Muon_dxy[i], entry.Muon_dz[i], entry.Muon_looseId[i], entry.Muon_isGlobal[i], entry.Muon_isTracker[i], entry.Muon_pfRelIso04_all[i], DR0, DR1
            if DR0 < mt['ll_DR'] or DR1 < mt['ll_DR']:
                if printOn : print("    fail muon DR  DR0={0:f} DR1={1:f} for muon={2:s} vs 0.pT={3:f} vs 1.pT={4:f}".format(DR0,DR1, str(i), pairList[0].Pt(), pairList[1].Pt()))
                continue

        if printOn : print("    Good muon i={0:d} lumi={1:s} run={2:s} event={3:s}".format(i, str(entry.luminosityBlock), str(entry.run), str(entry.event)), entry.Muon_pt[i], entry.Muon_looseId[i], entry.Muon_pfRelIso04_all[i])

        for j in range(entry.nTau):
            if entry.Muon_charge[i] != entry.Tau_charge[j]: continue
            # apply tau(h) selections
            if printOn : print("        tau j={0:d}".format(j))
            if abs(entry.Tau_eta[j]) > mt['tau_eta']:
                if printOn : print("        fail tau eta={0:f}".format(entry.Tau_eta[j]))
                continue
            if entry.Tau_pt[j] < mt['tau_pt']:
                if printOn : print("        fail tau  pt={0:f}".format(entry.Tau_pt[j]))
                continue
            if abs(entry.Tau_dz[j]) > mt['tau_dz']:
                if printOn : print("        fail tau  dz={0:f}".format(entry.Tau_dz[j]))
                continue
            #if mt['tau_ID']:
            #    if not entry.Tau_idDecayModeNewDMs[j]:
            #        if printOn : print("        fail tau idDecayModeNewDMs={0}".format(entry.Tau_idDecayModeNewDMs[j]))
            #        continue
	    if  mt['tau_decayMode'] and (entry.Tau_decayMode[j] == 5 or entry.Tau_decayMode[j] == 6) :
                if printOn : print("        fail tau decayMode={0:d}".format(entry.Tau_decayMode[j]))
                continue

            if signC!=0 and entry.Tau_charge[j] == signC : continue

            ''' # this is the old (pre-DeepTau) selection
	    if ord(entry.Tau_idAntiMu[j]) <= mt['tau_antiMu']: continue
            if ord(entry.Tau_idAntiEle[j]) <= mt['tau_antiEle']: continue
            if dch== 'eemt':
                if ord(entry.Tau_idAntiMu[j]) < mt['tau_eemt_antiMu']: continue
	    '''

	    if mt['tau_vJet'] > 0  and not ord(entry.Tau_idDeepTau2017v2p1VSjet[j]) & mt['tau_vJet'] > 0 :
                if printOn : print("        fail DeepTau vs. Jet={0:d}".format(ord(entry.Tau_idDeepTau2017v2p1VSjet[j])))
                continue
	    if mt['tau_vEle'] > 0 and not ord(entry.Tau_idDeepTau2017v2p1VSe[j]) & mt['tau_vEle'] > 0 :
                if printOn : print("        fail DeepTau vs. ele={0:d}".format(ord(entry.Tau_idDeepTau2017v2p1VSe[j])))
                continue
            if mt['tau_vMu'] > 0 and not ord(entry.Tau_idDeepTau2017v2p1VSmu[j]) & mt['tau_vMu'] > 0 :
                if printOn : print("        fail DeepTau vs.  mu={0:d}".format(ord(entry.Tau_idDeepTau2017v2p1VSmu[j])))
                continue
            
            tau_eta, tau_phi = entry.Tau_eta[j], entry.Tau_phi[j]
            dPhi = min(abs(tau_phi-mu_phi),2.*pi-abs(tau_phi-mu_phi))
            DR = sqrt(dPhi**2 + (tau_eta-mu_eta)**2) # tau(mu) vs. tau(h)
            if DR < mt['mt_DR']:
                if printOn : print("        fail mtDR DR={0:f} for tau={1:d}".format(DR, j))
                continue
            if isDCH2:
                DR0 = lTauDR(tau_eta, tau_phi, pairList[0]) #l1 vs. tau(h)
                DR1 = lTauDR(tau_eta, tau_phi, pairList[1]) #l2 vs. tau(h)
                if DR0 < mt['mt_DR'] or DR1 < mt['mt_DR']:
                    if printOn : print("        fail DR  DR0={0:f} DR1={1:f} for tau={2:i}, mu={3:i}".format(DR0,DR1, j ))
                    continue

            if printOn: print("        Tau j={0:d} passes all cuts.".format(j))
            muTauPairs.append([i,j])

    return muTauPairs


def compareMuTauPair(entry,pair1,pair2) :
    # a return value of True means that pair2 is "better" than pair 1 
    i1, i2, j1, j2 = pair1[0], pair2[0], pair1[1], pair2[1]
    if entry.Muon_pfRelIso04_all[i2]  <  entry.Muon_pfRelIso04_all[i1] : return True
    if entry.Muon_pfRelIso04_all[i2] ==  entry.Muon_pfRelIso04_all[i1] :
        if entry.Muon_pt[i2] >  entry.Muon_pt[i1] : return True
        if entry.Muon_pt[i2] == entry.Muon_pt[i1] :
            if entry.Tau_rawDeepTau2017v2p1VSjet[j2] > entry.Tau_rawDeepTau2017v2p1VSjet[j1] : return True   
            if entry.Tau_rawDeepTau2017v2p1VSjet[j2] == entry.Tau_rawDeepTau2017v2p1VSjet[j1] :
                if entry.Tau_pt[j2] > entry.Tau_pt[j1] : return True
    return False

def compareMuTauPairvspT(entry,tauPairList, printOn=False) :

    SumList=[]
    for i in range(0,len(tauPairList)) :
         
        j1, j2 = tauPairList[i][0], tauPairList[i][1] # look at leading pt tau in each pair
        if printOn : print 'appending now', entry.Muon_pt[j1] + entry.Tau_pt[j2], j1, j2, tauPairList[i]
        SumList.append(entry.Muon_pt[j1] + entry.Tau_pt[j2])

    maxI=SumList.index(max(SumList))

    return maxI


def getBestMuTauPair(entry,dch='mt',pairList=[],printOn=False, isDCH2=False, signC=0) :

    # form all possible pairs that satisfy DR requirement
    if printOn : print("Entering getBestMuTauPair()") 
    tauPairList = getMuTauPairs(entry,dch=dch,pairList=pairList,printOn=printOn,isDCH2=isDCH2, signC=signC) 


    if len(tauPairList) == 0 : 
        if printOn : print 'failed to find good MuTau Pair', dch
        return []

    if len(tauPairList) >= 1 :
        '''for i in range(len(tauPairList)-1, 0, -1) :
            if compareMuTauPair(entry, tauPairList[i], tauPairList[i-1]) :
               tauPairList[i-1], tauPairList[i] = tauPairList[i], tauPairList[i-1]
        bestPair = tauPairList[0]#best pair is the first element
        '''
        bestPair=tauPairList[-1]
        for i in range(0,len(tauPairList)-1):
            if not compareMuTauPair(entry, tauPairList[i], bestPair) :
               bestPair=tauPairList[i]

    #maxI= compareMuTauPairvspT(entry, tauPairList, printOn)
    #bestPair=tauPairList[maxI]

    return bestPair


def getEMuTauPairs(entry,dch='em',pairList=[],printOn=False, isDCH2=False, signC=0) :
    """ tauFun.getEMuTauPairs(): returns a list of suitable
                                 H-> tau(mu) + tau(ele) cands 
    """

    if printOn : print("Entering getEMuTauPairs() nMuon={0:d} nElectron={1:d}".format(entry.nMuon,entry.nElectron)) 
    if entry.nElectron < dch.count('e') or entry.nMuon < dch.count('m'): return []    

    elmuTauPairs = []
    em = selections['em'] # selections for H->tau(ele)+tau(mu)
    for i in range(entry.nMuon):

        # selections for tau(mu)
        if printOn : print("Muon i={0:d}".format(i))
        if em['mu_ID']:
            #if not (entry.Muon_looseId[i] or entry.Muon_mediumId[i] or entry.Muon_tightId[i]) :
            if not (entry.Muon_tightId[i]):
                if printOn : print("    failed muID={0}".format(entry.Muon_tightId[i]))
                continue

        if em['mu_type']:
            if not (entry.Muon_isGlobal[i] or entry.Muon_isTracker[i]) :
                if printOn : print("    fail mu_type Global or Tracker={0}".format(entry.Muon_isGlobal[i]))
                continue
        if abs(entry.Muon_dxy[i]) > em['mu_dxy']:
            if printOn : print("    failed dxy={0:f}".format(entry.Muon_dxy[i]))
            continue
        if abs(entry.Muon_dz[i]) > em['mu_dz']:
            if printOn : print("    failed dz={0:f}".format(entry.Muon_dz[i]))
            continue
        mu_eta, mu_phi = entry.Muon_eta[i], entry.Muon_phi[i] 
        if entry.Muon_pt[i] < em['mu_pt']:
            if printOn : print("    failed pt={0:f}".format(entry.Muon_pt[i]))
            continue
        if abs(mu_eta) > em['mu_eta']:
            if printOn : print("    failed eta={0:f}".format(entry.Muon_eta[i]))
            continue   
        if em['mu_iso_f'] and entry.Muon_pfRelIso04_all[i] > em['mu_iso']:
            if printOn : print("    failed iso={0:f}".format(entry.Muon_pfRelIso04_all[i]))
            continue 
        
        if isDCH2:     
            if signC!=0 and entry.Muon_charge[i] == signC : continue   
            DR0 = lTauDR(mu_eta,mu_phi,pairList[0]) #l1 vs. tau(mu)
            DR1 = lTauDR(mu_eta,mu_phi,pairList[1]) #l2 vs. tau(mu)
            if DR0 < em['ll_DR'] or DR1 < em['ll_DR']:
                if printOn : print("    failed DR  DR0={0:f} DR1={1:f}".format(DR0,DR1))
                continue
                        
        for j in range(entry.nElectron):
            if entry.Muon_charge[i] != entry.Electron_charge[j] : continue
            # selections for tau(ele)
            if printOn: print("    electron={0:d}".format(j))
            if abs(entry.Electron_dxy[j]) > em['ele_dxy']:
                if printOn : print("        failed dxy={0:f}".format(entry.Electron_dxy[j]))
                continue
            if abs(entry.Electron_dz[j]) > em['ele_dz']:
                if printOn : print("        failed dz={0:f}".format(entry.Electron_dz[j]))
                continue
            ele_eta, ele_phi = entry.Electron_eta[j], entry.Electron_phi[j] 
            if entry.Electron_pt[j] < em['ele_pt']:
                if printOn : print("        failed pt={0:f}".format(entry.Electron_pt[j]))
                continue
            if abs(ele_eta) > em['ele_eta']:
                if printOn : print("        failed eta={0:f}".format(entry.Electron_eta[j]))
                continue
            if ord(entry.Electron_lostHits[j]) > em['ele_lostHits']:
                if printOn : print("        failed lost hits={0:d}".format(entry.Electron_lostHits[j]))
                continue 
            if em['ele_convVeto']:
                if not entry.Electron_convVeto[j]:
                    if printOn : print("        failed conv. veto={0}".format(entry.Electron_convVeto[j]))
                    continue
            if em['ele_ID']:    
                if not entry.Electron_mvaFall17V2noIso_WP90[j]:
                #if entry.Electron_cutBased[j] < 2 :
                    if printOn : print("        failed mvaWP90={0}".format(entry.entry.Electron_mvaFall17V2noIso_WP90[j]))
                    continue
            if em['ele_iso_f'] and entry.Electron_pfRelIso03_all[j] > em['ele_iso']:
                if printOn : print("        failed iso={0:f}".format(entry.Electron_pfRelIso03_all[j]))
                continue

            if signC!=0 and entry.Electron_charge[j] == signC : continue

            dPhi = min(abs(mu_phi-ele_phi),2.*pi-abs(mu_phi-ele_phi))
            DR = sqrt(dPhi**2 + (mu_eta-ele_eta)**2) # tau(mu) vs. tau(ele)
            if DR < em['em_DR']:
                if printOn : print("        failed emDR={0:f}".format(DR))
                continue
            if isDCH2:
                DR0 = lTauDR(ele_eta,ele_phi,pairList[0]) # l1 vs. tau(ele) 
                DR1 = lTauDR(ele_eta,ele_phi,pairList[1]) # l2 vs. tau(ele)
                if DR0 < em['ll_DR'] or DR1 < em['ll_DR']:
                    if printOn : print("        failed ltDR DR0={0:f} DR1={1:f}".format(DR0,DR1))
                    continue
            if printOn : print("        found a good pair i={0:i}, j={1:i}".format(i,j))
            elmuTauPairs.append([j,i])

    return elmuTauPairs

'''
def compareEMuTauPair(entry,pair1,pair2) :
    # a return value of True means that pair2 is "better" than pair 1 
    i1, i2, j1, j2 = pair1[0], pair2[0], pair1[1], pair2[1]
    #if entry.Electron_mvaFall17Iso[i2]  < entry.Electron_mvaFall17Iso[i2] : return True
    if entry.Electron_pfRelIso03_all[i2]  < entry.Electron_pfRelIso03_all[i1] : return True
    #if entry.Electron_mvaFall17Iso[i2] == entry.Electron_mvaFall17Iso[i2] :
    if entry.Electron_pfRelIso03_all[i1] == entry.Electron_mvaFall17V2noIso_WP90[i2] : 
        if entry.Electron_pt[i2]  > entry.Electron_pt[i1] : return True 
        if entry.Electron_pt[i2] == entry.Electron_pt[i1] : 
            if entry.Muon_pt[j2] < entry.Muon_pt[j1] : return True   
    return False 
'''
def compareEMuPairvspT(entry,tauPairList,printOn=False) :

    SumList=[]
    for i in range(0,len(tauPairList)) :
         
        j1, j2 = tauPairList[i][0], tauPairList[i][1] # look at leading pt tau in each pair
        if printOn : print 'appending now', entry.Electron_pt[j1] + entry.Muon_pt[j2], j1, j2, tauPairList[i]
        SumList.append(entry.Electron_pt[j1] + entry.Muon_pt[j2])

    
    maxI=SumList.index(max(SumList))

    return maxI




def compareEMuTauPair(entry,pair1,pair2) :
    # a return value of True means that pair2 is "better" than pair 1 
    i1, i2, j1, j2 = pair1[0], pair2[0], pair1[1], pair2[1]
    if entry.Muon_pfRelIso04_all[j2]  < entry.Muon_pfRelIso04_all[j1] : return True
    if entry.Muon_pfRelIso04_all[j1] == entry.Muon_pfRelIso04_all[j2] : 
	if entry.Muon_pt[j2]  > entry.Muon_pt[j1] : return True 
	if entry.Muon_pt[j2] == entry.Muon_pt[j1] :
            if entry.Electron_pfRelIso03_all[i2] < entry.Electron_pfRelIso03_all[i1] : return True
            if entry.Electron_pfRelIso03_all[i1] == entry.Electron_pfRelIso03_all[i2] :
	        if entry.Electron_pt[i2] > entry.Electron_pt[i1] : return True   
    return False 




def getBestEMuTauPair(entry,dch,pairList=[],printOn=False, isDCH2=False,signC=0) :

    if printOn : print("Entering getBestEMuTauPair")
    # form all possible pairs that satisfy DR requirement
    tauPairList = getEMuTauPairs(entry,dch=dch,pairList=pairList,printOn=printOn, isDCH2=isDCH2, signC=signC) 

    # Sort the pair list using a bubble sort
    # The list is not fully sorted, since only the top pairing is needed

    if len(tauPairList) == 0 : 
        if printOn : print 'failed to find good EMu Pair', dch
        return []
    if len(tauPairList) >= 1 :
        '''for i in range(len(tauPairList)-1, 0, -1) :
            if compareEMuTauPair(entry, tauPairList[i], tauPairList[i-1]) :
               tauPairList[i-1], tauPairList[i] = tauPairList[i], tauPairList[i-1]
        bestPair = tauPairList[0]#best pair is the first element
        '''
        bestPair=tauPairList[-1]
        for i in range(0,len(tauPairList)-1):
            if not compareEMuTauPair(entry, tauPairList[i], bestPair) :
               bestPair=tauPairList[i]
    #maxI= compareEMuPairvspT(entry, tauPairList, printOn)
    #bestPair=tauPairList[maxI]

    return bestPair

def checkOverlapElectron(entry,i,j, checkDR=False, printOn=False) :

    overlapEl = False

    et = selections['eeoverlap'] # selections for H->tau(ele)+tau(h)
  
    if printOn : print("Electron i={0:d}".format(i))

    if abs(entry.Electron_dxy[i]) < et['ele_dxy']:
        if abs(entry.Electron_dz[i]) < et['ele_dz']:
            if et['ele_ID'] and entry.Electron_mvaFall17V2noIso_WP90[i]:
            #if et['ele_ID'] and not entry.Electron_cutBased[i] < 2 :
                if ord(entry.Electron_lostHits[i]) < et['ele_lostHits']:
                    if et['ele_convVeto'] and  entry.Electron_convVeto[i]:
                        if et['ele_iso_f'] and  entry.Electron_pfRelIso03_all[i] < et['ele_iso']:
                            if entry.Electron_pt[i] > et['ele_pt']:
                                if abs(entry.Electron_eta[i]) < et['ele_eta']:

				    if abs(entry.Electron_dxy[j]) < et['ele_dxy']:
					if abs(entry.Electron_dz[j]) < et['ele_dz']:
					    if et['ele_ID'] and entry.Electron_mvaFall17V2noIso_WP90[j]:
                                            #if et['ele_ID'] and not entry.Electron_cutBased[j] < 2 :
						if ord(entry.Electron_lostHits[j]) < et['ele_lostHits']:
						    if et['ele_convVeto'] and  entry.Electron_convVeto[j]:
							if et['ele_iso_f'] and  entry.Electron_pfRelIso03_all[j] < et['ele_iso']:
							    if entry.Electron_pt[j] > et['ele_pt']:
								if abs(entry.Electron_eta[j]) < et['ele_eta']:

								    if printOn: print 'checking overlap....', i, j, DRobj(entry.Electron_eta[i],entry.Electron_phi[i], entry.Electron_eta[j],entry.Electron_phi[j])
								    if checkDR : 
                                                                        if DRobj(entry.Electron_eta[i],entry.Electron_phi[i], entry.Electron_eta[j],entry.Electron_phi[j])  < et['ll_DR'] : overlapEl = True
                                                                    else : overlapEl = True
    return overlapEl


def getETauPairs(entry,dch='et',pairList=[],printOn=False, isDCH2=False, signC=0) :
    """ tauFun.getETauPairs(): get suitable pairs of  
                               H -> tau(ele) + tau(h) 
    """

    if printOn : print("Entering getETauPairs() nElectron={0:d} nTau={1:d}".format(entry.nElectron,entry.nTau)) 
    if entry.nElectron < 1 or entry.nTau < 1: return []
    
    eTauPairs = []
    leptOverlap = False
    et = selections['et'] # selections for H->tau(ele)+tau(h)
    if dch ==0:
    #if entry.nElectron>1 and len(pairList)>0:#dropping too many events
	for i in range(entry.nElectron) :
	    #if entry.Electron_pt[i]==pairList[0].Pt() or entry.Electron_pt[i]==pairList[1].Pt() : continue # make sure that you don't consider the findZ leptons
	    if DRobj(entry.Electron_eta[i],entry.Electron_phi[i], pairList[0].Eta(), pairList[0].Phi())<0.1 or DRobj(entry.Electron_eta[i],entry.Electron_phi[i], pairList[1].Eta(), pairList[1].Phi())<0.1 : continue
	    for j in range(i+1,entry.nElectron) :
		#if entry.Electron_pt[j]==pairList[0].Pt() or entry.Electron_pt[j]==pairList[1].Pt() : continue# make sure that you don't consider the findZ leptons
	        if DRobj(entry.Electron_eta[j],entry.Electron_phi[j], pairList[0].Eta(), pairList[0].Phi())<0.1 or DRobj(entry.Electron_eta[j],entry.Electron_phi[j], pairList[1].Eta(), pairList[1].Phi())<0.1 : continue
		if checkOverlapElectron(entry,i,j): return []


    for i in range(entry.nElectron) :

        # selections for tau(ele)
        if printOn : print("Electron i={0:d}".format(i))
        if abs(entry.Electron_dxy[i]) > et['ele_dxy']:
            if printOn : print("    failed dxy={0:f}".format(entry.Electron_dxy[i]))
            continue
        if abs(entry.Electron_dz[i]) > et['ele_dz']:
            if printOn : print("    failed dz={0:f}".format(entry.Electron_dz[i]))
            continue
        if et['ele_ID']:
            if not entry.Electron_mvaFall17V2noIso_WP90[i]:
            #if entry.Electron_cutBased[i] < 2 :
                if printOn : print("    failed mva={0}".format(entry.Electron_mvaFall17V2noIso_WP90[i]))
                continue
        if ord(entry.Electron_lostHits[i]) > et['ele_lostHits']:
            if printOn : print("    failed losthits={0:s}".format(entry.Electron_lostHits[i]))
            continue 
        if et['ele_convVeto']:
            if not entry.Electron_convVeto[i]:
                if printOn : print("    failed convVeto={0}".format(entry.Electron_convVeto[i]))
                continue
        if et['ele_iso_f'] :
            if entry.Electron_pfRelIso03_all[i] > et['ele_iso']:
                if printOn : print("    failed convVeto={0:f}".format(entry.Electron_pfRelIso03_all[i]))
                continue

        if entry.Electron_pt[i] < et['ele_pt']:
            if printOn : print("    failed pt={0:f}".format(entry.Electron_pt[i]))
            continue
        
        ele_eta, ele_phi = entry.Electron_eta[i], entry.Electron_phi[i]
        if abs(ele_eta) > et['ele_eta']:
            if printOn : print("    failed eta={0:f}".format(entry.Electron_eta[i]))
            continue

        if isDCH2:
            if signC!=0 and entry.Electron_charge[i] == signC : continue
            DR0 = lTauDR(ele_eta,ele_phi,pairList[0]) # l1 vs. tau(ele)
            DR1 = lTauDR(ele_eta,ele_phi,pairList[1]) # l2 vs. tau(ele)
            if DR0 < et['ll_DR'] or DR1 < et['ll_DR']:
                if printOn : print("    failed ltDR DR0={0:f} DR1={1:f}".format(DR0,DR1))
                continue

        for j in range(entry.nTau) :
            if entry.Electron_charge[i] != entry.Tau_charge[j]: continue
            # selections for tau(h)
            if printOn : print("    tau={0:d}".format(j))
            if entry.Tau_pt[j] < et['tau_pt']:
                if printOn : print("        failed pt={0:f}".format(entry.Tau_pt[j]))
                continue
            if abs(entry.Tau_eta[j]) > et['tau_eta']:
                if printOn : print("        failed eta={0:f}".format(entry.Tau_eta[j]))
                continue

            #if et['tau_ID'] and not entry.Tau_idDecayModeNewDMs[j]:
            #    if printOn : print("        failed idDecayMode={0}".format(entry.idDecayModeNewDMs[j]))
            #    continue
	    if  et['tau_decayMode'] and (entry.Tau_decayMode[j] == 5 or entry.Tau_decayMode[j] == 6) :
                if printOn : print("        failed DecayMode={0}".format(entry.Tau_decayMode[j]))
                continue
            if abs(entry.Tau_dz[j]) > et['tau_dz']:
                if printOn : print("        failed dz={0:f}".format(entry.Tau_dz[j]))
                continue
            if abs(entry.Tau_charge[j]) != 1:
                if printOn : print("        failed tauCharge={0:d}".format(entry.Tau_charge[j]))
                continue

            if signC!=0 and entry.Tau_charge[j] == signC : continue

            '''
	    if ord(entry.Tau_idAntiMu[j]) <= et['tau_antiMu']: continue
            if ord(entry.Tau_idAntiEle[j]) <= et['tau_antiEle']: continue
            if dch == 'eeet':
                if ord(entry.Tau_idAntiEle[j]) < et['tau_eeet_antiEle']: continue
	    '''

            '''
            if not ord(entry.Tau_idDeepTau2017v2p1VSjet[j]) & 16 > 0 : continue
	    if not ord(entry.Tau_idDeepTau2017v2p1VSe[j]) & 32 > 0 : continue
            if not ord(entry.Tau_idDeepTau2017v2p1VSmu[j]) & 1 > 0 : continue
            '''
            
            if et['tau_vJet'] > 0  and not ord(entry.Tau_idDeepTau2017v2p1VSjet[j]) & et['tau_vJet'] > 0 :
                if printOn : print("        fail DeepTau vs. Jet={0:d}".format(ord(entry.Tau_idDeepTau2017v2p1VSjet[j])))
                continue
	    if et['tau_vEle'] > 0 and not ord(entry.Tau_idDeepTau2017v2p1VSe[j]) & et['tau_vEle'] > 0 :
                if printOn : print("        fail DeepTau vs. ele={0:d}".format(ord(entry.Tau_idDeepTau2017v2p1VSe[j])))
                continue
            if et['tau_vMu'] > 0 and not ord(entry.Tau_idDeepTau2017v2p1VSmu[j]) & et['tau_vMu'] > 0 :
                if printOn : print("        fail DeepTau vs.  mu={0:d}".format(ord(entry.Tau_idDeepTau2017v2p1VSmu[j])))
                continue


            tau_eta, tau_phi = entry.Tau_eta[j], entry.Tau_phi[j]
            dPhi = min(abs(tau_phi-ele_phi),2.*pi-abs(tau_phi-ele_phi))
            DR = sqrt(dPhi**2 + (tau_eta-ele_eta)**2)
            if DR < et['tt_DR']:
                if printOn : print("        failed ttDR={0:f}".format(DR))
                continue # tau(ele) vs. tau(h)

            if isDCH2:
                DR0 = lTauDR(tau_eta,tau_phi,pairList[0]) # l1 vs. tau(h)
                DR1 = lTauDR(tau_eta,tau_phi,pairList[1]) # l2 vs. tau(h)
                if DR0 < et['tt_DR'] or DR1 < et['tt_DR']:
                    if printOn : print("        failed ltDR DR0={0:f} DR1={1:f}".format(DR0,DR1))
                    continue

            eTauPairs.append([i,j])
    return eTauPairs

def compareETauPairvspT(entry,tauPairList, printOn=False) :

    SumList=[]
    for i in range(0,len(tauPairList)) :
         
        j1, j2 = tauPairList[i][0], tauPairList[i][1] # look at leading pt tau in each pair
        if printOn : print 'appending now', entry.Electron_pt[j1] + entry.Tau_pt[j2], j1, j2, tauPairList[i]
        SumList.append(entry.Electron_pt[j1] + entry.Tau_pt[j2])

    
    maxI=SumList.index(max(SumList))

    return maxI


def compareETauPair(entry,pair1,pair2) :
    # a return value of True means that pair2 is "better" than pair 1 
    i1, i2, j1, j2 = pair1[0], pair2[0], pair1[1], pair2[1]
    if entry.Electron_pfRelIso03_all[i2]  < entry.Electron_pfRelIso03_all[i1] : return True
    if entry.Electron_pfRelIso03_all[i1] == entry.Electron_pfRelIso03_all[i2] : 
        if entry.Electron_pt[i2]  > entry.Electron_pt[i1] : return True 
        if entry.Electron_pt[i2] == entry.Electron_pt[i1] : 
            if entry.Tau_rawDeepTau2017v2p1VSjet[j2] > entry.Tau_rawDeepTau2017v2p1VSjet[j1] : return True   
            if entry.Tau_rawDeepTau2017v2p1VSjet[j2] == entry.Tau_rawDeepTau2017v2p1VSjet[j1] :
                if entry.Tau_pt[j2] > entry.Tau_pt[j1] : return True
    return False 

def getBestETauPair(entry,dch,pairList=[],printOn=False, isDCH2=False, signC=0) :

    if printOn : print("Entering getBestETauPair")
    # form all possible pairs that satisfy DR requirement
    tauPairList = getETauPairs(entry,dch=dch,pairList=pairList,printOn=printOn, isDCH2=isDCH2, signC=signC) 


    if len(tauPairList) == 0 : 
        if printOn : print 'failed to find good ETau Pair', dch
        return []
    if len(tauPairList) >= 1 :
        '''
        for i in range(len(tauPairList)-1, 0, -1) :
            if compareETauPair(entry, tauPairList[i], tauPairList[i-1]) :
               tauPairList[i-1], tauPairList[i] = tauPairList[i], tauPairList[i-1]
        bestPair = tauPairList[0]#best pair is the first element
        '''
        bestPair=tauPairList[-1]
        for i in range(0,len(tauPairList)-1):
            if not compareETauPair(entry, tauPairList[i], bestPair) :
               bestPair=tauPairList[i]
        #print bestPair, entry.Electron_pt[bestPair[0]]+entry.Tau_pt[bestPair[1]]

    #maxI= compareETauPairvspT(entry, tauPairList, printOn)
    #bestPair=tauPairList[maxI]

    return bestPair



def getEEPairs(entry, dch='ee', pairList=[], printOn=False, isDCH2=False, signC=0):
    
    if printOn: print ("Entering getEEPairs(): nElectron={0:d}".format(entry.nElectron))

    # need a sufficient number of leptons
    if entry.nElectron < 2: return []
    
    selected_elecs = []
    ee = selections['et'] # impose selections for tau(ele) on each electron
    
    # get a list of suitable electrons
    for i in range(entry.nElectron):

        if printOn: print("Electron i={0:d}".format(i),entry.luminosityBlock,  entry.run, entry.event)
        
        if abs(entry.Electron_dxy[i]) > ee['ele_dxy']:
            if printOn: print("\t failed dxy={0:f}".format(entry.Electron_dxy[i]))
            continue
        if abs(entry.Electron_dz[i]) > ee['ele_dz']:
            if printOn: print("\t failed dz={0:f}".format(entry.Electron_dz[i]))
            continue
        if ee['ele_ID']:
            if not entry.Electron_mvaFall17V2noIso_WP90[i]:
            #if entry.Electron_cutBased[i] < 2 :
                if printOn: print("\t failed mva={0}".format(entry.Electron_mvaFall17V2noIso_WP90[i]))
                continue
        if ord(entry.Electron_lostHits[i]) > ee['ele_lostHits']:
            if printOn: print("\t failed losthits={0}".format(entry.Electron_lostHits[i]))
            continue
        if ee['ele_convVeto']:
            if not entry.Electron_convVeto[i]:
                if printOn: print("\t failed convVeto={0}".format(entry.Electron_convVeto[i]))
                continue
        if ee['ele_iso_f']:
            if entry.Electron_pfRelIso03_all[i] > ee['ele_iso']:
                if printOn: print("\t failed convVeto={0:f}".format(entry.Electron_pfRelIso03_all[i]))
                continue
        if entry.Electron_pt[i] < ee['ele_pt']:
            if printOn: print("\t failed pt={0:f}".format(entry.Electron_pt[i]))
            continue

        ele_eta, ele_phi = entry.Electron_eta[i], entry.Electron_phi[i]
        if abs(ele_eta) > ee['ele_eta']:
            if printOn: print("\t failed eta={0:f}".format(entry.Electron_eta[i]))
            continue


        if isDCH2:
            if signC!=0 and entry.Electron_charge[i] == signC : continue
            DR0 = lTauDR(ele_eta,ele_phi,pairList[0]) # l1 vs. e1
            DR1 = lTauDR(ele_eta,ele_phi,pairList[1]) # l2 vs. e2
            if DR0 < ee['ll_DR'] or DR1 < ee['ll_DR']:
                if printOn: print("\t failed DR(l_i,e) check: DR0={0:f} DR1={1:f}".format(DR0,DR1))
                continue

        selected_elecs.append(i)

    # pair up suitable electrons
    ee_pairs = []
    for i in selected_elecs:
        for j in selected_elecs:
            if (j <= i): continue
            if entry.Electron_charge[i] != entry.Electron_charge[j] : continue
            if printOn: print("\t considering (i,j)=({0:d}, {1:d})".format(i, j))
            
            e1_eta, e1_phi = entry.Electron_eta[i], entry.Electron_phi[i]
            e2_eta, e2_phi = entry.Electron_eta[j], entry.Electron_eta[j]
            dPhi = min(abs(e1_phi-e2_phi), 2.*pi-abs(e1_phi-e2_phi))
            DR = sqrt(dPhi**2 + (e1_eta-e2_eta)**2)

            if DR < ee['tt_DR']:
                if printOn: print("\t failed eeDR={0:f}".format(DR))
                continue
            
            # store in leading, sub-leading order
            if (entry.Electron_pt[i] > entry.Electron_pt[j]):
                ee_pairs.append([i,j])
            else: 
                ee_pairs.append([j,i])

    return ee_pairs


def compareEEPairs(entry, pair1, pair2):
    i1, i2, j1, j2 = pair1[0], pair2[0], pair1[1], pair2[1]
    
    if entry.Electron_pfRelIso03_all[i2]  < entry.Electron_pfRelIso03_all[i1]:
        if entry.Electron_pfRelIso03_all[j2] <= entry.Electron_pfRelIso03_all[j1]:
            return True

    if entry.Electron_pfRelIso03_all[j2] < entry.Electron_pfRelIso03_all[j1]:
         if entry.Electron_pfRelIso03_all[i2]  <= entry.Electron_pfRelIso03_all[i1]:
             return True

    if entry.Electron_pfRelIso03_all[i1] == entry.Electron_pfRelIso03_all[i2]:
        if entry.Electron_pfRelIso03_all[j1] == entry.Electron_pfRelIso03_all[j2]:
            if (entry.Electron_pt[i2]+entry.Electron_pt[j2]) > (entry.Electron_pt[i1]+entry.Electron_pt[j1]):
                return True
                
    return False


def getBestEEPair(entry, dch, pairList=[], printOn=False, isDCH2=False,signC=0):
    
    if printOn: print("Entering getBestEEPair")
    
    ee_pairs = getEEPairs(entry, dch=dch, pairList=pairList, printOn=printOn,isDCH2=isDCH2,signC=signC)
    if len(ee_pairs) == 0: return []
    if len(ee_pairs) >=1:
        '''
        for i in range(len(ee_pairs)-1, 0, -1):
            if compareEEPairs(entry, ee_pairs[i], ee_pairs[i-1]):
                ee_pairs[i-1], ee_pairs[i] = ee_pairs[i], ee_pairs[i-1]
        bestPair = ee_pairs[0]
        '''
        bestPair=ee_pairs[-1]
        for i in range(0,len(ee_pairs)-1):
            if not compareEEPairs(entry, ee_pairs[i], bestPair) :
                bestPair=ee_pairs[i]
        #print bestPair, entry.Electron_pt[bestPair[0]]+entry.Electron_pt[bestPair[1]]
           
    return bestPair


def getMuMuPairs(entry, dch='mm', pairList=[], printOn=False, isDCH2=False,signC=0):
    
    if entry.nMuon < 2:
        if printOn: print ("Entering getMuMuPairs, failing nMuon={0:d}".format(entry.nMuon))
        return []

    if printOn: print("Entering tauFun.getMuMuPairs() nMuon={0:d}".format(entry.nMuon))

    mm = selections['mt'] # inherit selections for tau(mu)
    selected_muons = []
    for i in range(entry.nMuon):
        
        if mm['mu_type']:
            if not (entry.Muon_isGlobal[i] or entry.Muon_isTracker[i]) :
                if printOn: print("\t fail mu_type Global or Tracker={0}".format(entry.Muon_isGlobal[i]))
                continue
        if mm['mu_ID']:
            #if not (entry.Muon_looseId[i] or entry.Muon_mediumId[i] or entry.Muon_tightId[i]) :
            if not (entry.Muon_tightId[i]):
                if printOn: print("\t fail mu_ID mediumId={0}".format(entry.Muon_tightId[i]))
                continue
        if abs(entry.Muon_dxy[i]) > mm['mu_dxy']:
            if printOn: print("\t fail mu_dxy={0:f}".format(entry.Muon_dxy[i]))
            continue
        if abs(entry.Muon_dz[i]) > mm['mu_dz']:
            if printOn: print("\t fail mu_dz={0:f}".format(entry.Muon_dz[i]))
            continue

        mu_eta, mu_phi = entry.Muon_eta[i], entry.Muon_phi[i]
        if entry.Muon_pt[i] < mm['mu_pt']:
            if printOn: print("\t fail mu_pt={0:f}".format(entry.Muon_pt[i]))
            continue
        if abs(mu_eta) > mm['mu_eta']:
            if printOn: print("\t fail mu_eta={0:f}".format(entry.Muon_eta[i]))
            continue
        if  mm['mu_iso_f'] and entry.Muon_pfRelIso04_all[i] > mm['mu_iso']:
            if printOn: print("\t fail mu_iso={0:f}".format(entry.Muon_pfRelIso04_all[i]))
            continue

        if isDCH2:
            if signC!=0 and entry.Muon_charge[i] == signC : continue
            DR0 = lTauDR(mu_eta,mu_phi,pairList[0]) # l1 vs. m
            DR1 = lTauDR(mu_eta,mu_phi,pairList[1]) # l2 vs. m
            if DR0 < mm['ll_DR'] or DR1 < mm['ll_DR']:
                if printOn : print("\t fail muon DR  DR0={0:f} DR1={1:f}".format(DR0,DR1))
                continue

        if printOn : print("\t Good muon i={0:d}".format(i))
        selected_muons.append(i)
        
    # pair up the selected muons
    mm_pairs = []
    for i in selected_muons:
        for j in selected_muons:    
            if j <= i: continue
            if entry.Muon_charge[i] != entry.Muon_charge[j] : continue
            if printOn: print("\t considering (i,j)=({0:d}, {1:d})".format(i, j))
            
            m1_eta, m1_phi = entry.Muon_eta[i], entry.Muon_phi[i]
            m2_eta, m2_phi = entry.Muon_eta[j], entry.Muon_phi[j]
            
            dPhi = min(abs(m1_phi-m2_phi), 2.*pi-abs(m1_phi-m2_phi))
            DR = sqrt(dPhi**2 + (m1_eta-m2_eta)**2)
            if DR < mm['mt_DR']:
                if printOn: print("\t fail mmDR DR={0:f}".format(DR))
                continue
            
            # store in [leading, sub-leading] order
            if entry.Muon_pt[i] > entry.Muon_pt[j]: 
                mm_pairs.append([i,j])
            else: 
                mm_pairs.append([j,i])
            
    return mm_pairs


def compareMuMuPairs(entry, pair1, pair2):
    i1, i2, j1, j2 = pair1[0], pair2[0], pair1[1], pair2[1]

    if entry.Muon_pfRelIso04_all[i2]  < entry.Muon_pfRelIso04_all[i1]:
        if entry.Muon_pfRelIso04_all[j2] <= entry.Muon_pfRelIso04_all[j1]:
            return True

    if entry.Muon_pfRelIso04_all[j2] < entry.Muon_pfRelIso04_all[j1]:
         if entry.Muon_pfRelIso04_all[i2]  <= entry.Muon_pfRelIso04_all[i1]:
             return True

    if entry.Muon_pfRelIso04_all[i1] == entry.Muon_pfRelIso04_all[i2]:
        if entry.Muon_pfRelIso04_all[j1] == entry.Muon_pfRelIso04_all[j2]:
            if (entry.Muon_pt[i2]+entry.Muon_pt[j2]) > (entry.Muon_pt[i1]+entry.Muon_pt[j1]):
                return True
    
    return False


def getBestMuMuPair(entry, dch='mm', pairList=[], printOn=False, isDCH2=False,signC=0):
    
    # form all possible pairs that satisfy DR requirement
    if printOn: print("Entering getBestMuMuPair()")
    mm_pairs = getMuMuPairs(entry,dch=dch, pairList=pairList,printOn=printOn,isDCH2=isDCH2,signC=signC)
    if len(mm_pairs) == 0: return []
    if len(mm_pairs) >=1:
        '''
        for i in range(len(mm_pairs)-1, 0, -1):
            if compareMuMuPairs(entry, mm_pairs[i], mm_pairs[i-1]):
                mm_pairs[i-1], mm_pairs[i] = mm_pairs[i], mm_pairs[i-1]
        bestPair = mm_pairs[0]
        '''
        bestPair=mm_pairs[-1]
        for i in range(0,len(mm_pairs)-1):
            if not compareMuMuPairs(entry, mm_pairs[i], bestPair) :
                bestPair=mm_pairs[i]
        #print bestPair, entry.Muon_pt[bestPair[0]]+entry.Muon_pt[bestPair[1]]

    return bestPair


def goodPhoton(entry, j ):
    
    gj = selections['gjets'] # selections for Z->mumu
    if entry.Photon_pt[j] < gj['photon_pt']: return False
    if abs(entry.Photon_eta[j]) > gj['photon_eta']: return False
    if entry.Photon_cutBased[j] !=  gj['photon_cutbased']: return False
             
    return True 

# select a muon for the Z candidate
def goodMuon(entry, j ):
    """ tauFun.goodMuon(): select good muons
                           for Z -> mu + mu
    """
    
    mm = selections['mm'] # selections for Z->mumu
    if entry.Muon_pt[j] < mm['mu_pt']: return False
    if abs(entry.Muon_eta[j]) > mm['mu_eta']: return False
    if mm['mu_iso_f'] and entry.Muon_pfRelIso04_all[j] >  mm['mu_iso']: return False
    if mm['mu_ID'] :
        #if not (entry.Muon_looseId[j] or entry.Muon_mediumId[j] or entry.Muon_tightId[j]): return False
        if not (entry.Muon_tightId[j]): return False
    if abs(entry.Muon_dxy[j]) > mm['mu_dxy']: return False 
    if abs(entry.Muon_dz[j]) > mm['mu_dz']: return False
    if mm['mu_type'] :
        if not (entry.Muon_isGlobal[j] or entry.Muon_isTracker[j]) : return False
             
    return True 

def goodMuonWjets(entry, j ):
    """ tauFun.goodMuon(): select good muons
                           for W -> mu + nu
    """
    
    mm = selections['mm'] # selections for Z->mumu
    if entry.Muon_pt[j] < mm['mu_pt']: return False
    if abs(entry.Muon_eta[j]) > mm['mu_eta']: return False
    if mm['mu_iso_f'] and entry.Muon_pfRelIso04_all[j] >  mm['mu_iso']: return False
    if mm['mu_ID'] :
        if not (entry.Muon_tightId[j]) or ord(entry.Muon_pfIsoId[j]) < 4: return False
        #if not (entry.Muon_looseId[j]): return False
    if abs(entry.Muon_dxy[j]) > mm['mu_dxy']: return False 
    if abs(entry.Muon_dz[j]) > mm['mu_dz']: return False
    if mm['mu_type'] :
        if not (entry.Muon_isGlobal[j] or entry.Muon_isTracker[j]) : return False
    return True 

def vetoMuon(entry, j ):
    #if mm['mu_type'] :
    #    if not (entry.Muon_isGlobal[j] or entry.Muon_isTracker[j]) : return False
    
    mm = selections['muveto'] # selections for Z->mumu
    if entry.Muon_pt[j]  > mm['mu_pt'] and \
    abs(entry.Muon_eta[j]) < mm['mu_eta'] and \
    mm['mu_iso_f'] and entry.Muon_pfRelIso04_all[j] <  mm['mu_iso'] and \
    mm['mu_ID'] and  (entry.Muon_mediumId[j] ) and \
    abs(entry.Muon_dxy[j]) < mm['mu_dxy'] and \
    abs(entry.Muon_dz[j]) < mm['mu_dz'] :  return True 

             

def makevetoMuonList(entry) :
    vetoMuonList = []
    for i in range(entry.nMuon) :
        if vetoMuon(entry, i) : vetoMuonList.append(i)
    #print("In tauFun.makeGoodMuonList = {0:s}".format(str(goodMuonList)))
    return vetoMuonList

def makeGoodPhotonList(entry) :
    goodPhotonList = []
    for i in range(entry.nPhoton) :
        if goodPhoton(entry, i) : goodPhotonList.append(i)
    return goodPhotonList


def makeGoodMuonList(entry) :
    goodMuonList = []
    for i in range(entry.nMuon) :
        if goodMuon(entry, i) : goodMuonList.append(i)
    #print("In tauFun.makeGoodMuonList = {0:s}".format(str(goodMuonList)))
    return goodMuonList

def makeGoodMuonListWjets(entry) :
    goodMuonList = []
    for i in range(entry.nMuon) :
        if goodMuonWjets(entry, i) : goodMuonList.append(i)
    #print("In tauFun.makeGoodMuonList = {0:s}".format(str(goodMuonList)))
    return goodMuonList

# select an electron for the Z candidate
def goodElectron(entry, j) :
    """ tauFun.goodElectron(): select good electrons 
                               for Z -> ele + ele
    """
    ee = selections['ee'] # selections for Z->ee
    if entry.Electron_pt[j] < ee['ele_pt']: return False
    if abs(entry.Electron_eta[j]) > ee['ele_eta']: return False
    if abs(entry.Electron_dxy[j]) > ee['ele_dxy']: return False
    if abs(entry.Electron_dz[j]) > ee['ele_dz']: return False

    #if entry.luminosityBlock==159 and entry.run==320853 and entry.event==197659973: print entry.nElectron, j, 'pt: ', entry.Electron_pt[j], 'eta :', abs(entry.Electron_eta[j]), 'hits: ', ord(entry.Electron_lostHits[j]), 'iso :', entry.Electron_pfRelIso03_all[j], 'veto: ',entry.Electron_convVeto[j], 'iD: ', entry.Electron_mvaFall17V2noIso_WP90[j], entry.event

    if ord(entry.Electron_lostHits[j]) > ee['ele_lostHits']: return False
    if ee['ele_iso_f'] and entry.Electron_pfRelIso03_all[j] >  ee['ele_iso']: return False
    if ee['ele_convVeto']:
        if not entry.Electron_convVeto[j]: return False
    if ee['ele_ID']:
        if not entry.Electron_mvaFall17V2noIso_WP90[j] : return False
    return True 

def goodElectronWjets(entry, j) :
    """ tauFun.goodElectron(): select good electrons 
                               for Z -> ele + ele
    """
    ee = selections['ee'] # selections for Z->ee
    if entry.Electron_pt[j] < ee['ele_pt']: return False
    if abs(entry.Electron_eta[j]) > ee['ele_eta']: return False
    if abs(entry.Electron_eta[j]) >  1.44 and abs(entry.Electron_eta[j]) < 1.57 : return False

    if abs(entry.Electron_dxy[j]) > ee['ele_dxy']: return False
    if abs(entry.Electron_dz[j]) > ee['ele_dz']: return False

    if ord(entry.Electron_lostHits[j]) > ee['ele_lostHits']: return False
    if ee['ele_iso_f'] and entry.Electron_pfRelIso03_all[j] >  ee['ele_iso']: return False
    if ee['ele_convVeto']:
        if not entry.Electron_convVeto[j]: return False
    if ee['ele_ID']:
        if not entry.Electron_mvaFall17V2Iso_WP90[j] : return False

    return True 


def vetoElectron(entry, j) :
    ee = selections['elveto'] # selections for Z->ee
    if entry.Electron_pt[j] > ee['ele_pt'] and \
    abs(entry.Electron_eta[j]) < ee['ele_eta'] and \
    abs(entry.Electron_dxy[j]) < ee['ele_dxy'] and \
    abs(entry.Electron_dz[j]) < ee['ele_dz'] and \
    ord(entry.Electron_lostHits[j]) < ee['ele_lostHits'] and \
    entry.Electron_convVeto[j] and \
    entry.Electron_pfRelIso03_all[j] <  ee['ele_iso'] and\
    entry.Electron_cutBased[j] == 1: return True
    #entry.Electron_mvaFall17V2Iso_WP90[j] : return True

    #entry.Electron_cutBased[j] == ee['ele_cut'] : return True

    #entry.Electron_mvaFall17V2noIso_WP90[j] and\
    #if entry.luminosityBlock==159 and entry.run==320853 and entry.event==197659973: print entry.nElectron, j, 'pt: ', entry.Electron_pt[j], 'eta :', abs(entry.Electron_eta[j]), 'hits: ', ord(entry.Electron_lostHits[j]), 'iso :', entry.Electron_pfRelIso03_all[j], 'veto: ',entry.Electron_convVeto[j], 'iD: ', entry.Electron_mvaFall17V2noIso_WP90[j], entry.event
    
    #if ord(entry.Electron_lostHits[j]) > ee['ele_lostHits']: return False
    #if ee['ele_iso_f'] and entry.Electron_pfRelIso03_all[j] >  ee['ele_iso']: return False
    #if ee['ele_convVeto']:
    #    if not entry.Electron_convVeto[j]: return False
    #if ee['ele_ID']:
    #    if not entry.Electron_mvaFall17V2noIso_WP90[j] : return False

    return False

def makevetoElectronList(entry) :
    vetoElectronList = []
    for i in range(entry.nElectron) :
        if vetoElectron(entry, i) : vetoElectronList.append(i)
    return vetoElectronList


def makeGoodElectronList(entry) :
    goodElectronList = []
    for i in range(entry.nElectron) :
        if goodElectron(entry, i) : goodElectronList.append(i)
    return goodElectronList

def makeGoodElectronListWjets(entry) :
    goodElectronList = []
    for i in range(entry.nElectron) :
        if goodElectronWjets(entry, i) : goodElectronList.append(i)
    return goodElectronList

#Double Charged Higgs cuts
def goodElectronDCH(entry, j,cutflow) :
    """ tauFun.goodElectron(): select good electrons 
                               for Z -> ele + ele
    """
    ee = selections['ee'] # selections for Z->ee
    cutflow.count('cut0')

    if entry.Electron_pt[j] < ee['ele_pt']: 
        cutflow.count('cut1')
        return False
    if abs(entry.Electron_eta[j]) > ee['ele_eta']:
        cutflow.count('cut2')
        return False
    if abs(entry.Electron_dxy[j]) > ee['ele_dxy']:
        cutflow.count('cut3')
        return False
    if abs(entry.Electron_dz[j]) > ee['ele_dz']:
        cutflow.count('cut4')
        return False

    #if ord(entry.Electron_lostHits[j]) > ee['ele_lostHits']: return False
    if ee['ele_iso_f'] and entry.Electron_pfRelIso03_all[j] >  ee['ele_iso']:
        cutflow.count('cut5')
        return False
    if ee['ele_convVeto']:
        if not entry.Electron_convVeto[j]:
            cutflow.count('cut6')
            return False
    if ee['ele_ID']:
        #if not entry.Electron_cutBased_HEEP[j] : return False
        if not entry.Electron_mvaFall17V2noIso_WP90[j] :
            cutflow.count('cut7')
            return False
        #if entry.Electron_cutBased[j] < 2 : return False
    return True 

def makeGoodElectronListDCH(entry,cutflow) :
    goodElectronList = []
    for i in range(entry.nElectron) :
        if goodElectronDCH(entry, i,cutflow) : goodElectronList.append(i)
    return goodElectronList

def goodMuonDCH(entry, j ,cutflow):
    """ tauFun.goodMuon(): select good muons
                           for Z -> mu + mu
    """
    
    mm = selections['mm'] # selections for Z->mumu
    cutflow.count('cut0')

    if entry.Muon_pt[j] < mm['mu_pt']:
        cutflow.count('cut1')
        return False
    if abs(entry.Muon_eta[j]) > mm['mu_eta']:
        cutflow.count('cut2')
        return False
    if mm['mu_iso_f'] and entry.Muon_pfRelIso04_all[j] >  mm['mu_iso']:
        cutflow.count('cut3')
        return False
    #if mm['mu_iso_f'] and entry.Muon_tkRelIso[j] >  mm['mu_iso']: return False
    #if mm['mu_ID'] and not entry.Muon_looseId[j] : return False
    if abs(entry.Muon_dxy[j]) > mm['mu_dxy']:
        cutflow.count('cut4')
        return False
    if abs(entry.Muon_dz[j]) > mm['mu_dz']:
        cutflow.count('cut5')
        return False
    if mm['mu_type'] :
        if not (entry.Muon_isGlobal[j] or entry.Muon_isTracker[j]) :
            cutflow.count('cut6')
            return False
    if mm['mu_ID'] :
        if not entry.Muon_tightId[j] :
            cutflow.count('cut7')
            return False 
    return True 

def makeGoodMuonListDCH(entry,cutflow) :
    goodMuonList = []
    for i in range(entry.nMuon) :
        if goodMuonDCH(entry, i,cutflow) : goodMuonList.append(i)
    #print("In tauFun.makeGoodMuonList = {0:s}".format(str(goodMuonList)))
    return goodMuonList

#Double Charged Higgs cuts


##############
def goodMuonExtraLepton(entry, j):
    """ tauFun.goodMuon(): select good muons
                           for Z -> mu + mu
    """
    
    mm = selections['mextra'] # selections for Z->mumu
    if entry.Muon_pt[j] < mm['mu_pt']: return False
    if abs(entry.Muon_eta[j]) > mm['mu_eta']: return False
    if mm['mu_iso_f'] and entry.Muon_pfRelIso04_all[j] >  mm['mu_iso']: return False
    if mm['mu_ID'] and not (entry.Muon_mediumId[j] or entry.Muon_tightId[j]): return False
    if abs(entry.Muon_dxy[j]) > mm['mu_dxy']: return False 
    if abs(entry.Muon_dz[j]) > mm['mu_dz']: return False
    if mm['mu_type'] :
        if not (entry.Muon_isGlobal[j] or entry.Muon_isTracker[j]) : return False
             
    return True 

def makeGoodMuonListExtraLepton(entry, listMu) :
    goodExtraMuonList = []
    for i in range(entry.nMuon) :
        if i not in listMu and goodMuonExtraLepton(entry, i) : goodExtraMuonList.append(i)
    #print("In tauFun.makeGoodMuonList = {0:s}".format(str(goodMuonList)))
    return goodExtraMuonList

# select an electron for the Z candidate
def goodElectronExtraLepton(entry, j) :
    """ tauFun.goodElectron(): select good electrons 
                               for Z -> ele + ele
    """
    ee = selections['eextra'] # selections for Z->ee
    if entry.Electron_pt[j] < ee['ele_pt']: return False
    if abs(entry.Electron_eta[j]) > ee['ele_eta']: return False
    if abs(entry.Electron_dxy[j]) > ee['ele_dxy']: return False
    if abs(entry.Electron_dz[j]) > ee['ele_dz']: return False
    if ord(entry.Electron_lostHits[j]) > ee['ele_lostHits']: return False
    if ee['ele_iso_f'] and entry.Electron_pfRelIso03_all[j] >  ee['ele_iso']: return False
    if ee['ele_convVeto']:
        if not entry.Electron_convVeto[j]: return False
    if ee['ele_ID']:
        if not entry.Electron_mvaFall17V2noIso_WP90[j] : return False
        #if entry.Electron_cutBased[j] < 2 : return False
    return True 

def makeGoodElectronListExtraLepton(entry, listEl) :
    goodExtraElectronList = []
    for i in range(entry.nElectron) :
        if i not in listEl and goodElectronExtraLepton(entry, i) : goodExtraElectronList.append(i)
    return goodExtraElectronList

def eliminateCloseTauAndLepton(entry, goodElectronList, goodMuonList, goodTauList) :

    badMuon, badElectron, badTau = [], [], []
    # check tau vs tau
    for tau1 in goodTauList :
        for tau2 in goodTauList :
            if tau1 == tau2 : continue
            dEta = abs(entry.Tau_eta[tau1] - entry.Tau_eta[tau2])
            if dEta > 0.3 : continue
            dPhi = abs(entry.Tau_phi[tau1] - entry.Tau_phi[tau2])
            if dPhi > 0.3 : continue
            if sqrt(dEta*dEta + dPhi*dPhi) > 0.3 : continue
            if not (tau1 in badTau) and entry.Tau_pt[tau1] < entry.Tau_pt[tau2] : badTau.append(tau1)
            if not (tau2 in badTau) and entry.Tau_pt[tau1] > entry.Tau_pt[tau2] : badTau.append(tau2)

        #check tau vs mu
        for mu2 in goodMuonList :
            dEta = abs(entry.Tau_eta[tau1] - entry.Muon_eta[mu2])
            if dEta > 0.3 : continue
            dPhi = abs(entry.Tau_phi[tau1] - entry.Muon_phi[mu2])
            if dPhi > 0.3 : continue
            if sqrt(dEta*dEta + dPhi*dPhi) > 0.3 : continue
            if not (tau1 in badTau) and entry.Tau_pt[tau1] < entry.Muon_pt[mu2] : badTau.append(tau1)
            if not (mu2 in badMuon) and entry.Tau_pt[tau1] > entry.Muon_pt[mu2] : badMuon.append(mu2)

        #check tau vs el
        for e2 in goodElectronList :
            dEta = abs(entry.Tau_eta[tau1] - entry.Electron_eta[e2])
            if dEta > 0.3 : continue
            dPhi = abs(entry.Tau_phi[tau1] - entry.Electron_phi[e2])
            if dPhi > 0.3 : continue
            if sqrt(dEta*dEta + dPhi*dPhi) > 0.3 : continue
            if not (tau1 in badTau) : badTau.append(tau1)
            if not (e2 in badElectron) : badElectron.append(e2)
            if not (tau1 in badTau) and entry.Tau_pt[tau1] < entry.Electron_pt[e2] : badTau.append(tau1)
            if not (e2 in badElectron) and entry.Tau_pt[tau1] > entry.Electron_pt[e2] : badElectron.append(e2)

    # check el vs el        
    for e1 in goodElectronList :
        for e2 in goodElectronList :
            if e1 == e2 : continue 
            dEta = abs(entry.Electron_eta[e1] - entry.Electron_eta[e2])
            if dEta > 0.3 : continue
            dPhi = abs(entry.Electron_phi[e1] - entry.Electron_phi[e2])
            if dPhi > 0.3 : continue
            if sqrt(dEta*dEta + dPhi*dPhi) > 0.3 : continue
            if not (e1 in badElectron) and entry.Electron_pt[e1] < entry.Electron_pt[e2] : badElectron.append(e1)
            if not (e2 in badElectron) and entry.Electron_pt[e1] > entry.Electron_pt[e2] : badElectron.append(e2)

    #check mu vs mu and vs el
    for mu1 in goodMuonList :
        for mu2 in goodMuonList :
            if mu1 == mu2 : continue
            dEta = abs(entry.Muon_eta[mu1] - entry.Muon_eta[mu2])
            if dEta > 0.3 : continue
            dPhi = abs(entry.Muon_phi[mu1] - entry.Muon_phi[mu2])
            if dPhi > 0.3 : continue
            if sqrt(dEta*dEta + dPhi*dPhi) > 0.3 : continue
            if not (mu1 in badMuon) and entry.Muon_pt[mu1] < entry.Muon_pt[mu2] : badMuon.append(mu1)
            if not (mu2 in badMuon) and entry.Muon_pt[mu1] > entry.Muon_pt[mu2] : badMuon.append(mu2)


        for e2 in goodElectronList :
            dEta = abs(entry.Muon_eta[mu1] - entry.Electron_eta[e2])
            if dEta > 0.3 : continue
            dPhi = abs(entry.Muon_phi[mu1] - entry.Electron_phi[e2])
            if dPhi > 0.3 : continue
            if sqrt(dEta*dEta + dPhi*dPhi) > 0.3 : continue
            if not (mu1 in badMuon) and entry.Muon_pt[mu1] < entry.Electron_pt[e2] : badMuon.append(mu1)
            if not (e2 in badElectron) and entry.Muon_pt[mu1] > entry.Electron_pt[e2] : badElectron.append(e2)



    for bade in badElectron : goodElectronList.remove(bade)
    for badmu in badMuon : goodMuonList.remove(badmu)
    for badtau in badTau : goodTauList.remove(badtau)

    return goodElectronList, goodMuonList, goodTauList


def eliminateCloseLeptons(entry, goodElectronList, goodMuonList) :
    badMuon, badElectron = [], []
    if len(goodMuonList)> 1 :
	for mu1 in goodMuonList :
	    for mu2 in goodMuonList :
		if mu1 == mu2 : continue
		dEta = abs(entry.Muon_eta[mu1] - entry.Muon_eta[mu2])
		if dEta > 0.3 : continue
		dPhi = abs(entry.Muon_phi[mu1] - entry.Muon_phi[mu2])
		if dPhi > 0.3 : continue
		if sqrt(dEta*dEta + dPhi*dPhi) > 0.3 : continue
		if not (mu1 in badMuon) and entry.Muon_pt[mu1] < entry.Muon_pt[mu2] : badMuon.append(mu1)
		if not (mu2 in badMuon) and entry.Muon_pt[mu1] > entry.Muon_pt[mu2] : badMuon.append(mu2)

	    for e2 in goodElectronList :
		dEta = abs(entry.Muon_eta[mu1] - entry.Electron_eta[e2])
		if dEta > 0.3 : continue
		dPhi = abs(entry.Muon_phi[mu1] - entry.Electron_phi[e2])
		if dPhi > 0.3 : continue
		if sqrt(dEta*dEta + dPhi*dPhi) > 0.3 : continue
		if not (mu1 in badMuon) and entry.Muon_pt[mu1] < entry.Electron_pt[e2] : badMuon.append(mu1)
		if not (e2 in badElectron) and entry.Muon_pt[mu1] > entry.Electron_pt[e2] : badElectron.append(e2)
    if len(goodElectronList)> 1 :        
	for e1 in goodElectronList :
	    for e2 in goodElectronList :
		if e1 == e2 : continue 
		dEta = abs(entry.Electron_eta[e1] - entry.Electron_eta[e2])
		if dEta > 0.3 : continue
		dPhi = abs(entry.Electron_phi[e1] - entry.Electron_phi[e2])
		if dPhi > 0.3 : continue
		if sqrt(dEta*dEta + dPhi*dPhi) > 0.3 : continue
		if not (e1 in badElectron) and entry.Electron_pt[e1] < entry.Electron_pt[e2] : badElectron.append(e1)
		if not (e2 in badElectron) and entry.Electron_pt[e1] > entry.Electron_pt[e2] : badElectron.append(e2)

    for bade in badElectron : goodElectronList.remove(bade)
    for badmu in badMuon : goodMuonList.remove(badmu)

    return goodElectronList, goodMuonList

def findETrigger(goodElectronList,entry,era):
    EltrigList =[]
    nElectron = len(goodElectronList)
    
    if nElectron > 1 :
	if era == '2016' and not entry.HLT_Ele27_WPTight_Gsf : return EltrigList
	if era == '2017' and not entry.HLT_Ele35_WPTight_Gsf : return EltrigList
        for i in range(nElectron) :
	    
            ii = goodElectronList[i]
            if era == '2016' and entry.Electron_pt[ii] < 29 : continue
            if era == '2017' and entry.Electron_pt[ii] < 37 : continue
            #print("Electron: pt={0:.1f} eta={1:.2f} phi={2:.2f}".format(entry.Electron_pt[ii], entry.Electron_eta[ii], entry.Electron_phi[ii]))
            #e1 = TLorentzVector()
            #e1.SetPtEtaPhiM(entry.Electron_pt[ii],entry.Electron_eta[ii],entry.Electron_phi[ii],0.0005)

            for iobj in range(0,entry.nTrigObj) :
	        dR = DRobj(entry.Electron_eta[ii],entry.Electron_phi[ii], entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj])
                #print("    Trg Obj: eta={0:.2f} phi={1:.2f} dR={2:.2f} bits={3:x}".format(
                    #entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj], dR, entry.TrigObj_filterBits[iobj]))
		if entry.TrigObj_filterBits[iobj] & 2  and dR < 0.5: ##that corresponds 0 WPTight
		    EltrigList.append(ii)
                    #print "======================= iobj", iobj, "entry.Trig",entry.TrigObj_id[iobj], "Bits", entry.TrigObj_filterBits[iobj]," dR", dR, "electron",i,"ii",ii,entry.TrigObj_id[iobj]

    return EltrigList


def findMuTrigger(goodMuonList,entry,era):
    MutrigList =[]
    nMuon = len(goodMuonList)
    
    if nMuon > 1 :
	if era == '2016' and not entry.HLT_IsoMu24 : return MutrigList
	if era == '2017' and not entry.HLT_IsoMu27 : return MutrigList
        for i in range(nMuon) :
	    

            ii = goodMuonList[i] 
	    if era == '2016' and entry.Muon_pt[ii] < 26 : continue
	    if era == '2017' and entry.Muon_pt[ii] < 29 : continue
            #print("Muon: pt={0:.1f} eta={1:.4f} phi={2:.4f}".format(entry.Muon_pt[ii], entry.Muon_eta[ii], entry.Muon_phi[ii]))
            for iobj in range(0,entry.nTrigObj) :
	        dR = DRobj(entry.Muon_eta[ii],entry.Muon_phi[ii], entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj])
                #print("    Trg Obj: eta={0:.4f} phi={1:.4f} dR={2:.4f} bits={3:x}".format(
                #    entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj], dR, entry.TrigObj_filterBits[iobj]))
		if entry.TrigObj_filterBits[iobj] & 8 or entry.TrigObj_filterBits[iobj] & 2 and dR < 0.5: ##that corresponds to Muon Trigger
		    MutrigList.append(ii)
                #print "======================= and === iobj", iobj, entry.TrigObj_id[iobj], "Bits", entry.TrigObj_filterBits[iobj]," dR", dR, "electron",i

    return MutrigList

def ComparepT(El, Mu, entry) :
    if entry.Electron_pt[El] > entry.Muon_pt[Mu] : return True
    else : return False


def make4Vec(lepList, dch, entry) :
    l1, l2, pairList = TLorentzVector(), TLorentzVector(), []
    if len(lepList) < 2  : 
       print 'FAILED lol'
       return []
    if dch=='ee':
       l1.SetPtEtaPhiM(entry.Electron_pt[lepList[0]],entry.Electron_eta[lepList[0]],entry.Electron_phi[lepList[0]],0.0005)
       l2.SetPtEtaPhiM(entry.Electron_pt[lepList[1]],entry.Electron_eta[lepList[1]],entry.Electron_phi[lepList[1]],0.0005)
    if dch=='em':
       l1.SetPtEtaPhiM(entry.Electron_pt[lepList[0]],entry.Electron_eta[lepList[0]],entry.Electron_phi[lepList[0]],0.0005)
       l2.SetPtEtaPhiM(entry.Muon_pt[lepList[1]],entry.Muon_eta[lepList[1]],entry.Muon_phi[lepList[1]],0.105)
    if dch=='et':
       l1.SetPtEtaPhiM(entry.Electron_pt[lepList[0]],entry.Electron_eta[lepList[0]],entry.Electron_phi[lepList[0]],0.0005)
       l2.SetPtEtaPhiM(entry.Tau_pt[lepList[1]],entry.Tau_eta[lepList[1]],entry.Tau_phi[lepList[1]],1.776)
    if dch=='mm':
       l1.SetPtEtaPhiM(entry.Muon_pt[lepList[0]],entry.Muon_eta[lepList[0]],entry.Muon_phi[lepList[0]],0.105)
       l2.SetPtEtaPhiM(entry.Muon_pt[lepList[1]],entry.Muon_eta[lepList[1]],entry.Muon_phi[lepList[1]],0.105)
    if dch=='mt':
       l1.SetPtEtaPhiM(entry.Muon_pt[lepList[0]],entry.Muon_eta[lepList[0]],entry.Muon_phi[lepList[0]],0.105)
       l2.SetPtEtaPhiM(entry.Tau_pt[lepList[1]],entry.Tau_eta[lepList[1]],entry.Tau_phi[lepList[1]],1.776)
    if dch=='tt':
       l2.SetPtEtaPhiM(entry.Tau_pt[lepList[0]],entry.Tau_eta[lepList[0]],entry.Tau_phi[lepList[0]],1.776)
       l2.SetPtEtaPhiM(entry.Tau_pt[lepList[1]],entry.Tau_eta[lepList[1]],entry.Tau_phi[lepList[1]],1.776)
    pairList=[l1,l2]
    return pairList

def findHpair(goodElectronList, goodMuonList, entry) :
    #mm = selections['mm'] # H->tau(mu)+tau(h) selections
    selpair,pairList, mH, bestDiff = [],[], 1000, 99999. 
    selpairP,pairListP, selpairM, pairListM = [],[], [], []
    #lepV = []
    nElectron = len(goodElectronList)
    #print 'going in tauFun', goodElectronList, nElectron, entry.event, entry.luminosityBlock, entry.run, bestDiff
    #if nElectron > 4 : continue
    if nElectron >= 2 :
        for i in range(nElectron) :
            ii = goodElectronList[i] 
            e1 = TLorentzVector()
            e1.SetPtEtaPhiM(entry.Electron_pt[ii],entry.Electron_eta[ii],entry.Electron_phi[ii],0.0005)
            if entry.Electron_charge[ii] >0 :
		for j in range(i+1,nElectron) :
		    jj = goodElectronList[j]
		    #print 'going in tauFun masses', goodElectronList, nElectron, entry.event, entry.luminosityBlock, entry.run, 'for', jj, ii, entry.Electron_charge[ii],  entry.Electron_charge[jj]
		    if entry.Electron_charge[ii] == entry.Electron_charge[jj] :
			e2 = TLorentzVector()
			e2.SetPtEtaPhiM(entry.Electron_pt[jj],entry.Electron_eta[jj],entry.Electron_phi[jj],0.0005)
			cand = e1 + e2
			mass = cand.M()
			#print 'going in tauFun masses', goodElectronList, nElectron, entry.event, entry.luminosityBlock, entry.run, bestDiff, 'is abs(mass-mZ > bestDiff', abs(mass-mZ), bestDiff, 'for', jj, ii
			#if abs(mass-mH) < bestDiff :
                        #    bestDiff = abs(mass-mH)
			#print 'masses', bestDiff, 'mass', mass, 'q_1', entry.Electron_charge[jj], 'q_2', entry.Electron_charge[ii], entry.event, entry.luminosityBlock, entry.run, 'elect', jj, ii, goodElectronList, entry.Electron_pt[ii], entry.Electron_pt[jj]
			if entry.Electron_pt[ii] > entry.Electron_pt[jj] :
			   pairListP= [e1,e2]
			   selpairP = [ii,jj]
			else : 
			   pairListP = [e2,e1]
			   selpairP = [jj,ii]

            if entry.Electron_charge[ii] <0 :
		for j in range(i+1,nElectron) :
		    jj = goodElectronList[j]
		    if entry.Electron_charge[ii] == entry.Electron_charge[jj] :
			e2 = TLorentzVector()
			e2.SetPtEtaPhiM(entry.Electron_pt[jj],entry.Electron_eta[jj],entry.Electron_phi[jj],0.0005)
			cand = e1 + e2
			mass = cand.M()
                        #print 'going in tauFun masses', goodElectronList, nElectron, entry.event, entry.luminosityBlock, entry.run, bestDiff, 'is abs(mass-mZ > bestDiff', abs(mass-mZ), bestDiff, 'for', jj, ii

                        #if abs(mass-mH) < bestDiff :
                        #    bestDiff = abs(mass-mH)
                        #print 'masses', bestDiff, 'mass', mass, 'q_1', entry.Electron_charge[jj], 'q_2', entry.Electron_charge[ii], entry.event, entry.luminosityBlock, entry.run, 'elect', jj, ii, goodElectronList, entry.Electron_pt[ii], entry.Electron_pt[jj]
                        if entry.Electron_pt[ii] > entry.Electron_pt[jj] :
                           pairListP= [e1,e2]
                           selpairP = [ii,jj]
                        else :
                           pairListP = [e2,e1]
                           selpairP = [jj,ii]
        #print "do this here Electrons=============>", selpairP, 'M', selpairM      
    nMuon = len(goodMuonList)
    if nMuon >= 2 :
        for i in range(nMuon) :
            ii = goodMuonList[i] 
            e1 = TLorentzVector()
            e1.SetPtEtaPhiM(entry.Muon_pt[ii],entry.Muon_eta[ii],entry.Muon_phi[ii],0.105)
            if entry.Muon_charge[ii] >0 :
		for j in range(i+1,nMuon) :
		    jj = goodMuonList[j]
		    #print 'going in tauFun masses', goodMuonList, nMuon, entry.event, entry.luminosityBlock, entry.run, 'for', jj, ii, entry.Muon_charge[ii],  entry.Muon_charge[jj]
		    if entry.Muon_charge[ii] == entry.Muon_charge[jj] :
			e2 = TLorentzVector()
			e2.SetPtEtaPhiM(entry.Muon_pt[jj],entry.Muon_eta[jj],entry.Muon_phi[jj],0.105)
			cand = e1 + e2
			mass = cand.M()
			#print 'going in tauFun masses', goodMuonList, nMuon, entry.event, entry.luminosityBlock, entry.run, bestDiff, 'is abs(mass-mZ > bestDiff', abs(mass-mZ), bestDiff, 'for', jj, ii
			#bestDiff = abs(mass)
			#print 'masses', bestDiff, 'mass', mass, 'q_1', entry.Muon_charge[jj], 'q_2', entry.Muon_charge[ii], entry.event, entry.luminosityBlock, entry.run, 'muons', jj, ii, goodMuonList, entry.Muon_pt[ii], entry.Muon_pt[jj]
			if entry.Muon_pt[ii] > entry.Muon_pt[jj] :
			    pairListP= [e1,e2]
			    selpairP = [ii,jj]

			else : 
			    pairListP = [e2,e1]
			    selpairP = [jj,ii]

            if entry.Muon_charge[ii] <0 :
		for j in range(i+1,nMuon) :
		    jj = goodMuonList[j]
		    if entry.Muon_charge[ii] == entry.Muon_charge[jj] :
			e2 = TLorentzVector()
			e2.SetPtEtaPhiM(entry.Muon_pt[jj],entry.Muon_eta[jj],entry.Muon_phi[jj],0.105)
			#cand = e1 + e2
			#mass = cand.M()
			#print 'going in tauFun masses', goodMuonList, nMuon, entry.event, entry.luminosityBlock, entry.run, bestDiff, 'is abs(mass-mZ > bestDiff', abs(mass-mZ), bestDiff, 'for', jj, ii
			#bestDiff = abs(mass)
			#print 'masses', bestDiff, 'mass', mass, 'q_1', entry.Muon_charge[jj], 'q_2', entry.Muon_charge[ii], entry.event, entry.luminosityBlock, entry.run, 'muons', jj, ii, goodMuonList, entry.Muon_pt[ii], entry.Muon_pt[jj]
			if entry.Muon_pt[ii] > entry.Muon_pt[jj] :
			    pairListM= [e1,e2]
			    selpairM = [ii,jj]

			else : 
			    pairListM = [e2,e1]
			    selpairM = [jj,ii]


        #print "do this here Muons=============>", selpairP, 'M', selpairM      


    # first particle of pair is positive
    #print 'returning', selpair,  'is muon', nMuon, goodMuonList, 'isEl', nElectron, goodElectronList, entry.event, entry.luminosityBlock, entry.run
    #lepV.append(e1)
    #lepV.append(e2)
    return pairListP, selpairP, pairListM, selpairM
                    


def findW(goodElectronList, goodMuonList, entry) :
    mm = selections['mm'] # H->tau(mu)+tau(h) selections
    selpair,pairList, mW, bestDiff = [],[], 80.38, 99999. 
    nElectron = len(goodElectronList)
    #print 'going in tauFun', goodElectronList, nElectron, entry.event, entry.luminosityBlock, entry.run, bestDiff
    if nElectron > 1 :
        for i in range(nElectron) :
            ii = goodElectronList[i] 
            e1 = TLorentzVector()
            e1.SetPtEtaPhiM(entry.Electron_pt[ii],entry.Electron_eta[ii],entry.Electron_phi[ii],0.0005)
            for j in range(i+1,nElectron) :
                jj = goodElectronList[j]
                #print 'going in tauFun masses', goodElectronList, nElectron, entry.event, entry.luminosityBlock, entry.run, 'for', jj, ii, entry.Electron_charge[ii],  entry.Electron_charge[jj]
                if entry.Electron_charge[ii] != entry.Electron_charge[jj] :
                    e2 = TLorentzVector()
                    e2.SetPtEtaPhiM(entry.Electron_pt[jj],entry.Electron_eta[jj],entry.Electron_phi[jj],0.0005)
                    cand = e1 + e2
                    mass = cand.M()
		    #if mass < 60 or mass > 120 : continue
                    #print 'going in tauFun masses', goodElectronList, nElectron, entry.event, entry.luminosityBlock, entry.run, bestDiff, 'is abs(mass-mZ > bestDiff', abs(mass-mZ), bestDiff, 'for', jj, ii
                    if abs(mass-mZ) < bestDiff :
                        bestDiff = abs(mass-mZ)
                        #print 'masses', bestDiff, mass, entry.Electron_charge[jj], entry.Electron_charge[ii], entry.event, entry.luminosityBlock, entry.run, 'elect', jj, ii, goodElectronList
                        if entry.Electron_charge[ii] > 0. :
                            pairList = [e1,e2]
                            selpair = [ii,jj]
                        else : 
                            pairList = [e2,e1]
                            selpair = [jj,ii]
                           
    nMuon = len(goodMuonList)
    if nMuon > 1 : 
        # find mass pairings
        for i in range(nMuon) :
            ii = goodMuonList[i]
            #if entry.Muon_pfRelIso04_all[ii] >  mm['mu_iso']: continue
            mu1 = TLorentzVector()
            mu1.SetPtEtaPhiM(entry.Muon_pt[ii],entry.Muon_eta[ii],entry.Muon_phi[ii],0.105)
            for j in range(i+1,nMuon) :
                jj = goodMuonList[j]
                if entry.Muon_charge[ii] != entry.Muon_charge[jj] :
                    mu2 = TLorentzVector()
                    mu2.SetPtEtaPhiM(entry.Muon_pt[jj],entry.Muon_eta[jj],entry.Muon_phi[jj],0.105)
                    cand = mu1 + mu2
                    mass = cand.M()
		    #if mass < 60 or mass > 120 : continue
                    if abs(mass-mZ) < bestDiff :
                        bestDiff = abs(mass-mZ)
                        if entry.Muon_charge[ii] > 0. :
                            pairList = [mu1,mu2]
                            selpair = [ii,jj]
                        else :
                            pairList = [mu2,mu1]
                            selpair = [jj,ii]

    # first particle of pair is positive
    #print 'returning', selpair,  'is muon', nMuon, goodMuonList, 'isEl', nElectron, goodElectronList, entry.event, entry.luminosityBlock, entry.run
    return pairList, selpair
                    







def findZ(goodElectronList, goodMuonList, entry) :
    mm = selections['mm'] # H->tau(mu)+tau(h) selections
    selpair,pairList, mZ, bestDiff = [],[], 91.19, 99999. 
    nElectron = len(goodElectronList)
    #print 'going in tauFun', goodElectronList, nElectron, entry.event, entry.luminosityBlock, entry.run, bestDiff
    if nElectron > 1 :
        for i in range(nElectron) :
            ii = goodElectronList[i] 
            e1 = TLorentzVector()
            e1.SetPtEtaPhiM(entry.Electron_pt[ii],entry.Electron_eta[ii],entry.Electron_phi[ii],0.0005)
            for j in range(i+1,nElectron) :
                jj = goodElectronList[j]
                #print 'going in tauFun masses', goodElectronList, nElectron, entry.event, entry.luminosityBlock, entry.run, 'for', jj, ii, entry.Electron_charge[ii],  entry.Electron_charge[jj]
                if entry.Electron_charge[ii] != entry.Electron_charge[jj] :
                    e2 = TLorentzVector()
                    e2.SetPtEtaPhiM(entry.Electron_pt[jj],entry.Electron_eta[jj],entry.Electron_phi[jj],0.0005)
                    cand = e1 + e2
                    mass = cand.M()
		    #if mass < 60 or mass > 120 : continue
                    #print 'going in tauFun masses', goodElectronList, nElectron, entry.event, entry.luminosityBlock, entry.run, bestDiff, 'is abs(mass-mZ > bestDiff', abs(mass-mZ), bestDiff, 'for', jj, ii
                    if abs(mass-mZ) < bestDiff :
                        bestDiff = abs(mass-mZ)
                        #print 'masses', bestDiff, mass, entry.Electron_charge[jj], entry.Electron_charge[ii], entry.event, entry.luminosityBlock, entry.run, 'elect', jj, ii, goodElectronList
                        if entry.Electron_charge[ii] > 0. :
                            pairList = [e1,e2]
                            selpair = [ii,jj]
                        else : 
                            pairList = [e2,e1]
                            selpair = [jj,ii]
                           
    nMuon = len(goodMuonList)
    if nMuon > 1 : 
        # find mass pairings
        for i in range(nMuon) :
            ii = goodMuonList[i]
            #if entry.Muon_pfRelIso04_all[ii] >  mm['mu_iso']: continue
            mu1 = TLorentzVector()
            mu1.SetPtEtaPhiM(entry.Muon_pt[ii],entry.Muon_eta[ii],entry.Muon_phi[ii],0.105)
            for j in range(i+1,nMuon) :
                jj = goodMuonList[j]
                if entry.Muon_charge[ii] != entry.Muon_charge[jj] :
                    mu2 = TLorentzVector()
                    mu2.SetPtEtaPhiM(entry.Muon_pt[jj],entry.Muon_eta[jj],entry.Muon_phi[jj],0.105)
                    cand = mu1 + mu2
                    mass = cand.M()
		    #if mass < 60 or mass > 120 : continue
                    if abs(mass-mZ) < bestDiff :
                        bestDiff = abs(mass-mZ)
                        if entry.Muon_charge[ii] > 0. :
                            pairList = [mu1,mu2]
                            selpair = [ii,jj]
                        else :
                            pairList = [mu2,mu1]
                            selpair = [jj,ii]

    # first particle of pair is positive
    #print 'returning', selpair,  'is muon', nMuon, goodMuonList, 'isEl', nElectron, goodElectronList, entry.event, entry.luminosityBlock, entry.run
    return pairList, selpair
                    
                    
def findZmumu(goodMuonList, entry) :
    pairList, mZ, bestDiff = [], 91.19, 99999.     
    nMuon = len(goodMuonList)
    if nMuon < 2 : return pairList 
    # find mass pairings
    for i in range(nMuon) :
        mu1 = TLorentzVector()
        mu1.SetPtEtaPhiM(entry.Muon_pt[i],entry.Muon_eta[i],entry.Muon_phi[i],0.105)
        for j in range(i+1,nMuon) :
            if entry.Muon_charge[i] != entry.Muon_charge[j] :
                mu2 = TLorentzVector()
                mu2.SetPtEtaPhiM(entry.Muon_pt[j],entry.Muon_eta[j],entry.Muon_phi[j],0.105)
                cand = mu1 + mu2
                mass = cand.M()
		#if mass < 60 or mass > 120 : continue
                if abs(mass-mZ) < bestDiff :
                    bestDiff = abs(mass-mZ)
                    pairList.append([mu1,mu2]) 

    return pairList

def findZee(goodElectronList, entry) :
    pairList, mZ, bestDiff = [], 91.19, 99999. 
    nElectron = len(goodElectronList)
    if nElectron < 2 : return pairList 
    # find mass pairings
    for i in range(nElectron) :
        e1 = TLorentzVector()
        e1.SetPtEtaPhiM(entry.Electron_pt[i],entry.Electron_eta[i],entry.Electron_phi[i],0.0005)
        for j in range(i+1,nElectron) :
            if entry.Electron_charge[i] != entry.Electron_charge[j] :
                e2 = TLorentzVector()
                e2.SetPtEtaPhiM(entry.Electron_pt[j],entry.Electron_eta[j],entry.Electron_phi[j],0.0005)
                cand = e1 + e2
                mass = cand.M()
		#if mass < 60 or mass > 120 : continue
                if abs(mass-mZ) < bestDiff :
                    bestDiff = abs(mass-mZ)
                    pairList.append([e1,e2]) 

    return pairList


def catToNumber(cat) :
    number = {'eeee':1,'eeem':2,'eeet':3,'eemm':4,'eemt':5,'eett':6,
               'emem':7,'emet':8,'emmm':9,'emmt':10,'emtt':11,
                      'etet':12,'etmm':13,'etmt':14,'ettt':15,
                             'mmmm':16,'mmmt':17,'mmtt':18,
                                    'mtmt':19,'mttt':20,
                                           'tttt':21}
    return number[cat]

def numberToCat(number) :
    cat = {1:'eeee',2:'eeem',3:'eeet',4:'eemm',5:'eemt',6:'eett',
               7:'emem',8:'emet',9:'emmm',10:'emmt',11:'emtt',
                      12:'etet',13:'etmm',14:'etmt',15:'ettt',
                             16:'mmmm',17:'mmmt',18:'mmtt',
                                    19:'mtmt',20:'mttt',
                                           21:'tttt'}
    return cat[number]


def catToNumber3L(cat) :
    #number = { 'eee':1, 'eem':2, 'eet':3, 'mme':4, 'mmm':5, 'mmt':6}
    #number = { 'ee':1, 'mm':2}
    number = { 'eee':22, 'eem':23, 'eet':24, 
               'eme':25, 'emm':26, 'emt':27, 
               'ete':28, 'etm':29, 'ett':30, 
               'mme':31, 'mmm':32, 'mmt':33, 
               'mte':34, 'mtm':35, 'mtt':36, 
               'tte':37, 'ttm':38, 'ttt':39}
    return number[cat]

def catToNumberW(cat) :
    number = { 'mnu':1, 'enu':2}
    return number[cat]

def numberToCat3L(number) :
    #cat = { 1:'eee', 2:'eem', 3:'eet', 4:'mme', 5:'mmm', 6:'mmt' }
    #cat = { 1:'ee', 2:'mm' }
    cat  = {   22:'eee', 23:'eem', 24:'eet',  
               25:'eme', 26:'emm', 27:'emt', 
               28:'ete', 29:'etm', 30:'ett', 
               31:'mme', 32:'mmm', 33:'mmt', 
               34:'mte', 35:'mtm', 36:'mtt',
               37:'tte', 38:'ttm', 39:'ttt'}
    return cat[number]



def catToNumber2Lep(cat) :
    number = { 'ee':1, 'mm':2}
    return number[cat]

def numberToCat2Lep(number) :
    cat = { 1:'ee', 2:'mm'}
    return cat[number]


