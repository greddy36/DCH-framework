# functions for H->tautau analysis CMSSW_10_2_X

from ROOT import TLorentzVector
from math import sqrt, sin, cos, pi
import subprocess

def getTauList(channel, entry) :
    # make a list of taus that pass the basic selection cuts
    if not channel in ['mmtt','eett','tt','mt','et'] :
        print("Invalid channel={0:s} in tauFun.getTauList()".format(channel))
        exit()
    tauList = []
    if entry.nTau == 0 : return tauList
    for j in range(entry.nTau) :
        if channel == 'tt' :
            if entry.Tau_pt[j] < 40. : continue
            if abs(entry.Tau_eta[j]) > 2.1 : continue
            if ord(entry.Tau_idMVAoldDM2017v2[j]) == 0 : continue
        elif channel == 'mmtt' or channel == 'eett' :
            if entry.Tau_pt[j] < 19. : continue
            if abs(entry.Tau_eta[j]) > 2.3 : continue
            if ord(entry.Tau_idAntiMu[j]) == 0 : continue
            if ord(entry.Tau_idAntiEle[j]) == 0 : continue
        if not entry.Tau_idDecayMode[j] : continue
        if abs(entry.Tau_dz[j]) > 0.2 : continue
        chg = abs(entry.Tau_charge[j])
        if chg < 0.5 or chg > 1.5 : continue 
        tauList.append(j)

    return tauList

def tauDR(entry, j1,j2) :
    if j1 == j2 : return 0. 
    phi1, eta1, phi2, eta2 = entry.Tau_phi[j1], entry.Tau_eta[j1], entry.Tau_phi[j2], entry.Tau_eta[j2]
    return sqrt( (phi2-phi1)**2 + (eta2-eta1)**2 )

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
        
def comparePair(entry, pair1, pair2) :
    # a return value of True means that pair2 is "better" than pair 1 
    j1, j2 = pair1[0], pair2[0]
    if entry.Tau_rawMVAoldDM2017v2[j2] < entry.Tau_rawMVAoldDM2017v2[j1] :
        return True
    elif  entry.Tau_rawMVAoldDM2017v2[j2] > entry.Tau_rawMVAoldDM2017v2[j1] :
        return False
    elif entry.Tau_pt[j2] > entry.Tau_pt[j1] :
        return True
    elif entry.Tau_pt[j2] < entry.Tau_pt[j1] :
        return False

    j1, j2 = pair1[1], pair2[1]
    if entry.Tau_rawMVAoldDM2017v2[j2] < entry.Tau_rawMVAoldDM2017v2[j1] :
        return True
    elif  entry.Tau_rawMVAoldDM2017v2[j2] > entry.Tau_rawMVAoldDM2017v2[j1] :
        return False
    elif entry.Tau_pt[j2] > entry.Tau_pt[j1] :
        return True
    else : 
        return False
        
def getBestTauPair(channel, entry, tauList) :
    if not channel in ['mmtt','eett','tt'] : 
        print("Invalid channel={0:s} in tauFun.getBestTauPair()".format(channel))
        exit()

    if len(tauList) < 2 : return [] 
    
    # form all possible pairs that satisfy DR requirement
    tauPairList = []
    for i1 in range(len(tauList)) :
        j1 = tauList[i1]
        for i2 in range(len(tauList)) :
            if i1 == i2 : continue
            j2 = tauList[i2]
            if tauDR(entry, j1, j2) < 0.5 : continue
            tauPairList.append([j1,j2])

    # Sort the pair list using a bubble sort
    # The list is not fully sorted, since only the top pairing is needed
    for i in range(len(tauPairList)-1,0,-1) :
        if comparePair(entry, tauPairList[i],tauPairList[i-1]) : 
            tauPairList[i-1], tauPairList[i] = tauPairList[i], tauPairList[i-1] 

    if len(tauPairList) == 0 : return []
    j1, j2 = tauPairList[0][0], tauPairList[0][1]
    if entry.Tau_pt[j2] > entry.Tau_pt[j1] : 
        temp = tauPairList[0][0]
        tauPairList[0][0] = tauPairList[0][1]
        tauPairList[0][1] = temp
        
    return tauPairList[0]

def getMuTauPairs(entry,printOn=False) :
    muTauPairs = [] 
    if entry.nMuon < 1 or entry.nTau < 1 :
        if printOn : print("Fail in getMuTauPairs nMuon={0:d} nTau={1:d}".format(entry.nMuon,entry.nTau))
        return muTauPairs
    
    for i in range(entry.nMuon) :
        if printOn: print("  Muon i={0:d}".format(i))
        if entry.Muon_pt[i] < 21. : continue
        if not entry.Muon_mediumId[i] : continue
        if printOn: print("  Muon ID OK")
        if abs(entry.Muon_dxy[i]) > 0.045 : continue
        if abs(entry.Muon_dz[i]) > 0.2 : continue
        eta1, phi1 = entry.Muon_eta[i], entry.Muon_phi[i] 
        if abs(eta1) > 2.1 : continue
        if printOn : print("  Good muon")
        for j in range(entry.nTau) :
            if printOn: print("     Tau  j={0:d}".format(j))
            if entry.Tau_pt[j] < 23. : continue
            if abs(entry.Tau_eta[j]) > 2.3 : continue
            if printOn: print("     Tau pt and eta OK") 
            if ord(entry.Tau_idMVAoldDM2017v2[j]) == 0 : continue
            if abs(entry.Tau_dz[j]) > 0.2 : continue
            if not entry.Tau_idDecayMode[j] : continue 
            if printOn: print("     Tau MVA OK") 
            eta2, phi2 = entry.Tau_eta[j], entry.Tau_phi[j]
            dPhi = min(abs(phi2-phi1),2.*pi-abs(phi2-phi1))
            DR = sqrt(dPhi**2 + (eta2-eta1)**2)
            if DR < 0.5 : continue
            if printOn: print("     Good mu-Tau pair") 
            muTauPairs.append([i,j])
    return muTauPairs

def compareMuTauPair(entry,pair1,pair2) :
    # a return value of True means that pair2 is "better" than pair 1 
    i1, i2, j1, j2 = pair1[0], pair2[0], pair1[1], pair2[1]
    if entry.Muon_pfRelIso04_all[i2]  <  entry.Muon_pfRelIso04_all[i1] : return True
    if entry.Muon_pfRelIso04_all[i2] ==  entry.Muon_pfRelIso04_all[i1] :
        if entry.Muon_pt[i2] >  entry.Muon_pt[i1] : return True
        if entry.Muon_pt[i2] == entry.Muon_pt[i1] :
            if entry.Tau_rawMVAoldDM2017v2[j2] < entry.Tau_rawMVAoldDM2017v2[j1] : return True   
    return False


def getBestMuTauPair(entry,printOn=False) :

    # form all possible pairs that satisfy DR requirement
    tauPairList = getMuTauPairs(entry,printOn=printOn) 

    # Sort the pair list using a bubble sort
    # The list is not fully sorted, since only the top pairing is needed
    for i in range(len(tauPairList)-1,0,-1) :
        if compareMuTauPair(entry, tauPairList[i],tauPairList[i-1]) : 
            tauPairList[i-1], tauPairList[i] = tauPairList[i], tauPairList[i-1] 

    if len(tauPairList) == 0 : return []
    return tauPairList[0]

        
def getETauPairs(entry,printOn=False) :
    if printOn : print("Enter getETauPairs nElectron={0:d} nTau={1:d}".format(entry.nElectron, entry.nTau))
    eTauPairs = [] 
    if entry.nElectron < 1 or entry.nTau < 1 : return eTauPairs 
    for i in range(entry.nElectron) :
        if printOn : print("   Try electron i={0:d}".format(i)) 
        if entry.Electron_pt[i] < 25. : continue
        if abs(entry.Electron_dxy[i]) > 0.045 : continue
        if abs(entry.Electron_dz[i]) > 0.2 : continue
        if printOn : print("   survived dz MVAnoIso={0:f}".format(entry.Electron_mvaFall17noIso[i]))
        #if not entry.Electron_mvaFall17noIso_WP90[i] : continue
        if not entry.Electron_mvaFall17V2noIso_WP90[i] : continue
        if printOn : print("   survived noIso") 
        if ord(entry.Electron_lostHits[i]) > 1 : continue 
        if not entry.Electron_convVeto[i] : continue
        eta1, phi1 = entry.Electron_eta[i], entry.Electron_phi[i] 
        if abs(eta1) > 2.1 : continue
        if printOn : print("   Good electron") 
        for j in range(entry.nTau) :
            if printOn : print("      Try tau j={0:d}".format(j)) 
            if entry.Tau_pt[j] < 23. : continue
            if abs(entry.Tau_eta[j]) > 2.3 : continue
            if ord(entry.Tau_idMVAoldDM2017v2[j]) == 0 : continue
            if not entry.Tau_idDecayMode[j] : continue
            if printOn : print("      Survived ID cuts")  
            if abs(entry.Tau_dz[j]) > 0.2 : continue
            #if abs(entry.Tau_dxy[j]) > 0.045 : continue
            eta2, phi2 = entry.Tau_eta[j], entry.Tau_phi[j]
            dPhi = min(abs(phi2-phi1),2.*pi-abs(phi2-phi1))
            DR = sqrt(dPhi**2 + (eta2-eta1)**2)
            if DR < 0.5 : continue
            chg = abs(entry.Tau_charge[j])
            if chg < 0.5 or chg > 1.5 : continue
            if printOn : print("      Good tau")  
            eTauPairs.append([i,j])
    return eTauPairs

def compareETauPair(entry,pair1,pair2) :
    # a return value of True means that pair2 is "better" than pair 1 
    i1, i2, j1, j2 = pair1[0], pair2[0], pair1[1], pair2[1]
    #if entry.Electron_mvaFall17Iso[i2]  < entry.Electron_mvaFall17Iso[i2] : return True
    if entry.Electron_mvaFall17V2Iso[i2]  < entry.Electron_mvaFall17V2Iso[i2] : return True
    #if entry.Electron_mvaFall17Iso[i2] == entry.Electron_mvaFall17Iso[i2] :
    if entry.Electron_mvaFall17V2Iso[i2] == entry.Electron_mvaFall17V2Iso[i2] : 
        if entry.Electron_pt[i2]  > entry.Electron_pt[i1] : return True 
        if entry.Electron_pt[i2] == entry.Electron_pt[i1] : 
            if entry.Tau_rawMVAoldDM2017v2[j2] < entry.Tau_rawMVAoldDM2017v2[j1] : return True   
    return False 

def getBestETauPair(entry,printOn=False) :

    if printOn : print("Entering getBestETauPair")
    # form all possible pairs that satisfy DR requirement
    tauPairList = getETauPairs(entry,printOn=printOn) 

    # Sort the pair list using a bubble sort
    # The list is not fully sorted, since only the top pairing is needed
    for i in range(len(tauPairList)-1,0,-1) :
        if compareETauPair(entry, tauPairList[i],tauPairList[i-1]) : 
            tauPairList[i-1], tauPairList[i] = tauPairList[i], tauPairList[i-1] 

    if len(tauPairList) == 0 : return []
    return tauPairList[0]


# select a muon for the Z candidate
def goodMuon(entry, j) :
    if entry.Muon_pt[j] < 10. : return False
    if abs(entry.Muon_eta[j]) > 2.4 : return False
    # drop this requirement for ZH
    #if not (entry.Muon_mediumId[j] or entry.Muon_tightId[j]): return False
    if entry.Muon_pfRelIso04_all[j] > 0.25 : return False
    if abs(entry.Muon_dxy[j]) > 0.045 : return False 
    if abs(entry.Muon_dz[j]) > 0.2 : return False
    return True 

def makeGoodMuonList(entry) :
    goodMuonList = []
    for i in range(entry.nMuon) :
        if goodMuon(entry, i) : goodMuonList.append(i)
    #print("In tauFun.makeGoodMuonList = {0:s}".format(str(goodMuonList)))
    return goodMuonList

# select an electron for the Z candidate
def goodElectron(entry, j) :
    if entry.Electron_pt[j] < 10. : return False
    if abs(entry.Electron_eta[j]) > 2.5 : return False
    if abs(entry.Electron_dxy[j]) > 0.045 : return False
    if abs(entry.Electron_dz[j]) > 0.2 : return False
    if ord(entry.Electron_lostHits[j]) > 1 : return False
    if not entry.Electron_convVeto[j]  : return False
    #if not entry.Electron_mvaFall17Iso_WP90[j] : return False
    if not entry.Electron_mvaFall17V2noIso_WP90[j] : return False
    return True 

def makeGoodElectronList(entry) :
    goodElectronList = []
    for i in range(entry.nElectron) :
        if goodElectron(entry, i) : goodElectronList.append(i)
    return goodElectronList

def eliminateCloseLeptons(entry, goodElectronList, goodMuonList) :
    badMuon, badElectron = [], []
    for mu1 in goodMuonList :
        for mu2 in goodMuonList :
            if mu1 == mu2 : continue
            dEta = abs(entry.Muon_eta[mu1] - entry.Muon_eta[mu2])
            if dEta > 0.3 : continue
            dPhi = abs(entry.Muon_phi[mu1] - entry.Muon_phi[mu2])
            if dPhi > 0.3 : continue
            if sqrt(dEta*dEta + dPhi*dPhi) > 0.3 : continue
            if not (mu1 in badMuon) : badMuon.append(mu1)
            if not (mu2 in badMuon) : badMuon.append(mu2)
        for e2 in goodElectronList :
            dEta = abs(entry.Muon_eta[mu1] - entry.Electron_eta[e2])
            if dEta > 0.3 : continue
            dPhi = abs(entry.Muon_phi[mu1] - entry.Electron_phi[e2])
            if dPhi > 0.3 : continue
            if sqrt(dEta*dEta + dPhi*dPhi) > 0.3 : continue
            if not (mu1 in badMuon) : badMuon.append(mu1)
            if not (e2 in badElectron) : badElectron.append(e2)
            
    for e1 in goodElectronList :
        for e2 in goodElectronList :
            if e1 == e2 : continue 
            dEta = abs(entry.Electron_eta[e1] - entry.Electron_eta[e2])
            if dEta > 0.3 : continue
            dPhi = abs(entry.Electron_phi[e1] - entry.Electron_phi[e2])
            if dPhi > 0.3 : continue
            if sqrt(dEta*dEta + dPhi*dPhi) > 0.3 : continue
            if not (e1 in badElectron) : badElectron.append(e1)
            if not (e2 in badElectron) : badElectron.append(e2)

    for bade in badElectron : goodElectronList.remove(bade)
    for badmu in badMuon : goodMuonList.remove(badmu)

    return goodElectronList, goodMuonList

def findZ(goodElectronList, goodMuonList, entry) :
    pairList, mZ, bestDiff = [], 91.19, 99999. 
    nElectron = len(goodElectronList)
    if nElectron > 1 :
        for i in range(nElectron) :
            ii = goodElectronList[i] 
            e1 = TLorentzVector()
            e1.SetPtEtaPhiM(entry.Electron_pt[ii],entry.Electron_eta[ii],entry.Electron_phi[ii],0.0005)
            for j in range(i+1,nElectron) :
                jj = goodElectronList[j]
                if entry.Electron_charge[ii] != entry.Electron_charge[jj] :
                    e2 = TLorentzVector()
                    e2.SetPtEtaPhiM(entry.Electron_pt[jj],entry.Electron_eta[jj],entry.Electron_phi[jj],0.0005)
                    cand = e1 + e2
                    mass = cand.M()
                    if abs(mass-mZ) < bestDiff :
                        bestDiff = abs(mass-mZ)
                        if entry.Electron_charge[ii] > 0. :
                            pairList = [e1,e2]
                        else : 
                            pairList = [e2,e1]
                            
    nMuon = len(goodMuonList)
    if nMuon > 1 : 
        # find mass pairings
        for i in range(nMuon) :
            ii = goodMuonList[i]
            mu1 = TLorentzVector()
            mu1.SetPtEtaPhiM(entry.Muon_pt[ii],entry.Muon_eta[ii],entry.Muon_phi[ii],0.105)
            for j in range(i+1,nMuon) :
                jj = goodMuonList[j]
                if entry.Muon_charge[ii] != entry.Muon_charge[jj] :
                    mu2 = TLorentzVector()
                    mu2.SetPtEtaPhiM(entry.Muon_pt[jj],entry.Muon_eta[jj],entry.Muon_phi[jj],0.105)
                    cand = mu1 + mu2
                    mass = cand.M()
                    if abs(mass-mZ) < bestDiff :
                        bestDiff = abs(mass-mZ)
                        if entry.Muon_charge[ii] > 0. :
                            pairList = [mu1,mu2]
                        else :
                            pairList = [mu2,mu1]

    # first particle of pair is positive 
    return pairList
                    
                    
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
                if abs(mass-mZ) < bestDiff :
                    bestDiff = abs(mass-mZ)
                    pairList.append([e1,e2]) 

    return pairList
                        
def mediumElectron(entry, j) :
    if entry.Electron_cutBased[j] < 4 : return False
    #if entry.Electron_cutBased_HLTPreSel[j] < 3 : return False
    if entry.Electron_miniPFRelIso_all[j] > 0.0588 : return False
    if entry.Electron_miniPFRelIso_all[j] > 0.0571 and entry.Electron_eta[j] > 1.479 : return False 
    if entry.Electron_tightCharge[j] < 2 : return False
    try :
        if ord(entry.Electron_lostHits[j]) > 0 : return False
    except TypeError :
        if entry.Electron_lostHits[j] > 0 : return False

    if abs(entry.Electron_eta[j]) < 1.479 :
        if abs(entry.Electron_dxy[j]) > 0.05 : return False
        if abs(entry.Electron_dz[j]) > 0.10 : return False
    else :
        if abs(entry.Electron_dxy[j]) > 0.10 : return False
        if abs(entry.Electron_dz[j]) > 0.20 : return False
    return True

def trigger(fName,entry) :
    try : 
        if entry.HLT_IsoMu20 : return True
    except AttributeError :
        pass

    try :
        if entry.HLT_IsoMu22 : return True
    except AttributeError :
        pass
    
    if entry.HLT_IsoMu24 : return True
    if entry.HLT_IsoMu27 : return True
    try:
        if entry.HLT_IsoTkMu18 : return True
    except AttributeError :
        pass 
    try:
        if entry.HLT_IsoTkMu20 : return True
    except AttributeError :
        pass 
    try:
        if entry.HLT_IsoTkMu24 : return True
    except AttributeError :
        pass
    try:
        if entry.HLT_IsoTkMu27 : return True
    except AttributeError :
        pass 
    try:
        if entry.HLT_Mu45_eta2p1 : return True
    except AttributeError :
        pass
    
    if entry.HLT_Mu50 : return True
        
    if entry.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL : return True 
    
    try :
        if entry.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL : return True 
    except AttributeError :
        pass

    if entry.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ : return True 

    try :
        if entry.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ : return True 
    except AttributeError :
        pass

    try:
        if entry.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ : return True
    except AttributeError :
        pass

    try :
        if entry.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL : return True
    except AttributeError :
        pass

    try :
        if entry.HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL : return True
    except AttributeError :
        pass

    try :
        if entry.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ : return True
    except AttributeError :
        pass

    try :
        if entry.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL : return True
    except AttributeError :
        pass

    try :
        if entry.HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ : return True
    except AttributeError :
        pass
        
    try:
        if entry.HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL : return True
    except AttributeError :
        pass 

    try :
        if entry.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ : return True
    except AttributeError :
        pass

    try :
        if entry.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL : return True        
    except AttributeError :
        pass

    try :
        if entry.HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL : return True    
    except AttributeError :
        pass
    
    try :
        if entry.HLT_Ele25_eta2p1_WPTight_Gsf : return True
    except AttributeError :
        pass

    try :
        if entry.HLT_Ele27_eta2p1_WPLoose_Gsf : return True        
    except AttributeError :
        pass

    try :
        if entry.HLT_Ele27_WPTight_Gsf : return True        
    except AttributeError :
        pass

    try :
        if entry.HLT_Ele30_WPTight_Gsf : return True        
    except AttributeError :
        pass

    try :
        if entry.HLT_Ele30_WPTight_Gsf : return True
    except AttributeError :
        pass
    
    if entry.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ : return True

    try :
        if entry.HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf : return True
    except AttributeError :
        pass
    
    return False





    
    
    

    
