# functions for H->tautau analysis 
from ROOT import TLorentzVector
from ROOT import TFile, TH1D, TCanvas, TGraph, kRed, kBlue, TLegend
from math import sqrt, sin, cos, pi
from termcolor import colored
import numpy as np
import json

def getLepIdxFrom4Vec(entry, lep_4vec, lep_type):
    """ Classification: genMatching
           - this function returns the index of a 4-vector (lep_4vec) belonging to
             the muon (lep_type = 'm') or electron (lep_type = 'e') collection
    """
    
    idx_lep = -1
    if lep_type == 'm':
        for i in range(entry.nMuon):
            if abs((lep_4vec.Pt() - entry.Muon_pt[i])/lep_4vec.Pt()) < 0.08:
                if abs((lep_4vec.Eta() - entry.Muon_eta[i])/lep_4vec.Eta()) < 0.08:
                    if abs((lep_4vec.Phi() - entry.Muon_phi[i])/lep_4vec.Phi()) < 0.08:
                        idx_lep = i

    if lep_type == 'e':
        for j in range(entry.nElectron):
            if abs((lep_4vec.Pt() - entry.Electron_pt[j])/lep_4vec.Pt()) < 0.08:
                    if abs((lep_4vec.Eta() - entry.Electron_eta[j])/lep_4vec.Eta()) < 0.08:
                        if abs((lep_4vec.Phi() - entry.Electron_phi[j])/lep_4vec.Phi()) < 0.08:
                            idx_lep = j

    return idx_lep

def printTriggerObjects(e) :
    print(' i    pt     phi   eta   ID   filter')
    for i in range(e.nTrigObj) :
        print('{0:3d}\t{1:5.1f}\t{2:5.2f}\t{3:5.2f}\t{4:5d}0x\t{5:04x}'.format(
            i,e.TrigObj_pt[i],  e.TrigObj_phi[i], e.TrigObj_eta[i], e.TrigObj_id[i], e.TrigObj_filterBits[i]))
    return

def genMatch(entry, jl,lep_type):#matches e/m/t to a gen paricle regardless of pID
    idx_match, dR_min, dPt_min = -99,99,99999
    if lep_type == 'e':
        for i in range(entry.nGenPart):
            if entry.GenPart_status[i] != 1: continue#make sure to look at visible gen particle
            dPt = abs(entry.GenPart_pt[i] - entry.Electron_pt[jl])
            dPhi = min(abs(entry.GenPart_phi[i] - entry.Electron_phi[jl]),
                       2.0*pi-abs(entry.GenPart_phi[i] - entry.Electron_phi[jl]))
            dEta = abs(entry.GenPart_eta[i] - entry.Electron_eta[jl])
            dR = sqrt(dPhi**2 + dEta**2)
            if dR <= dR_min and dPt <= dPt_min:
                idx_match, dR_min, dPt_min = i, dR, dPt
    if lep_type == 'm':
        for i in range(entry.nGenPart):
            if entry.GenPart_status[i] != 1: continue
            dPt = abs(entry.GenPart_pt[i] - entry.Muon_pt[jl])
            dPhi = min(abs(entry.GenPart_phi[i] - entry.Muon_phi[jl]),
                       2.0*pi-abs(entry.GenPart_phi[i] - entry.Muon_phi[jl]))
            dEta = abs(entry.GenPart_eta[i] - entry.Muon_eta[jl])
            dR = sqrt(dPhi**2 + dEta**2)
            if dR <= dR_min and dPt <= dPt_min:
                idx_match, dR_min, dPt_min = i, dR, dPt
    if lep_type == 't':
        for i in range(entry.nGenPart):
            #if entry.GenPart_status[i] != 2 : continue
            dPt = abs(entry.GenPart_pt[i] - entry.Tau_pt[jl])
            dPhi = min(abs(entry.GenPart_phi[i] - entry.Tau_phi[jl]),
                       2.0*pi-abs(entry.GenPart_phi[i] - entry.Tau_phi[jl]))
            dEta = abs(entry.GenPart_eta[i] - entry.Tau_eta[jl])
            dR = sqrt(dPhi**2 + dEta**2)
            if dR <= dR_min and dPt <= dPt_min:
                idx_match, dR_min, dPt_min = i, dR, dPt
    return idx_match

def genMatchTau(entry, jt, decayMode=''):
    """ Classification: genMatching
           - this function matches a hadronically-decaying tau (decayMode = 'had') to 
             a GenVisPart with  minimum dR(tau,GenVisPart[i])
           - it matches a leptonically-decaying tau (decayMode = 'lep') to a GenPart with
             abs(PDGID) = 15 with the smallest dR(tau,GenPart[i])
    """
    idx_match, dR_min = -99,99
    if decayMode == 'had':
        for i in range(entry.nGenVisTau):
            dPhi = min(abs(entry.GenVisTau_phi[i] - entry.Tau_phi[jt]),
                       2.0*pi-abs(entry.GenVisTau_phi[i] - entry.Tau_phi[jt]))
            dEta = abs(entry.GenVisTau_eta[i] - entry.Tau_eta[jt])
            dR = sqrt(dPhi**2 + dEta**2)
            if dR < dR_min:
                idx_match, dR_min = i, dR
    
    if decayMode == 'lep':
        for i in range(entry.nGenPart):
            if abs(entry.GenPart_pdgId[i]) != 15: continue
            dPhi = min(abs(entry.GenPart_phi[i] - entry.Tau_phi[jt]),
                       2.0*pi-abs(entry.GenPart_phi[i] - entry.Tau_phi[jt]))
            dEta = abs(entry.GenPart_eta[i] - entry.Tau_eta[jt])
            dR = sqrt(dPhi**2 + dEta**2)
            if dR < dR_min:
                idx_match, dR_min = i, dR

    return idx_match

def checkMETFlags(entry, year, isMC=False, proc="EOY") :
    METfilter = False
    '''
    if entry.Flag_BadChargedCandidateFilter or entry.Flag_BadChargedCandidateSummer16Filter or entry.Flag_BadPFMuonFilter
    or entry.Flag_BadPFMuonSummer16Filter or  entry.Flag_CSCTightHalo2015Filter or entry.Flag_CSCTightHaloFilter
    or entry.Flag_CSCTightHaloTrkMuUnvetoFilter or entry.Flag_EcalDeadCellBoundaryEnergyFilter	
    or entry.Flag_EcalDeadCellTriggerPrimitiveFilter or entry.Flag_HBHENoiseFilter or entry.Flag_HBHENoiseIsoFilter
    or entry.Flag_HcalStripHaloFilter	or entry.Flag_METFilters or entry.Flag_chargedHadronTrackResolutionFilter
    or entry.Flag_ecalBadCalibFilter	or entry.Flag_ecalBadCalibFilterV2 Flag_ecalLaserCorrFilter or entry.Flag_eeBadScFilter	
    or entry.Flag_globalSuperTightHalo2016Filter or entry.Flag_globalTightHalo2016Filter or entry.Flag_goodVertices 
    or entry.Flag_hcalLaserEventFilter	or entry.Flag_muonBadTrackFilter or entry.Flag_trkPOGFilters
    or entry.Flag_trkPOG_logErrorTooManyClusters or entry.Flag_trkPOG_manystripclus53X	or entry.Flag_trkPOG_toomanystripclus53X : 
    '''


    #if year== 2016 and (entry.Flag_goodVertices  == False or entry.Flag_HBHENoiseFilter == False or entry.Flag_HBHENoiseIsoFilter  == False or entry.Flag_EcalDeadCellTriggerPrimitiveFilter  == False or  entry.Flag_BadPFMuonFilter  == False or entry.Flag_ecalBadCalibFilter == False) : METfilter = True

    if proc=='EOY':
	if year== 2016 and (entry.Flag_goodVertices  == False or entry.Flag_globalSuperTightHalo2016Filter == False or entry.Flag_HBHENoiseFilter == False or entry.Flag_HBHENoiseIsoFilter  == False or entry.Flag_EcalDeadCellTriggerPrimitiveFilter  == False or  entry.Flag_BadPFMuonFilter  == False) : METfilter = True

	#if year== 2017 and (entry.Flag_goodVertices  == False or entry.Flag_globalSuperTightHalo2016Filter == False or entry.Flag_HBHENoiseFilter  == False or entry.Flag_HBHENoiseIsoFilter  == False or entry.Flag_EcalDeadCellTriggerPrimitiveFilter  == False or  entry.Flag_BadPFMuonFilter   == False or entry.Flag_ecalBadCalibFilterV2 == False) : METfilter = True

	#if year== 2018 and (entry.Flag_goodVertices  == False or entry.Flag_globalSuperTightHalo2016Filter == False or entry.Flag_HBHENoiseFilter  == False or entry.Flag_HBHENoiseIsoFilter  == False or entry.Flag_EcalDeadCellTriggerPrimitiveFilter  == False or  entry.Flag_BadPFMuonFilter == False  or  entry.Flag_ecalBadCalibFilterV2 == False) : METfilter = True
	#if year== 2018 and (Flag_METFilters == False ) : METfilter = True
	if year== 2018 and (entry.Flag_goodVertices  == False or entry.Flag_globalSuperTightHalo2016Filter == False or entry.Flag_HBHENoiseFilter  == False or entry.Flag_HBHENoiseIsoFilter  == False or entry.Flag_EcalDeadCellTriggerPrimitiveFilter  == False or  entry.Flag_BadPFMuonFilter == False  or entry.Flag_eeBadScFilter == False or entry.Flag_ecalBadCalibFilter == False) : METfilter = True


	#if not isMC and entry.Flag_eeBadScFilter == False : METfilter = True
    if proc=='UL' :
        if year ==2016  and (entry.Flag_goodVertices  == False or entry.Flag_globalSuperTightHalo2016Filter == False or entry.Flag_HBHENoiseFilter  == False or entry.Flag_HBHENoiseIsoFilter  == False or entry.Flag_EcalDeadCellTriggerPrimitiveFilter  == False or  entry.Flag_BadPFMuonFilter == False  or entry.Flag_eeBadScFilter == False ) : METfilter = True

        if year !=2016  and (entry.Flag_goodVertices  == False or entry.Flag_globalSuperTightHalo2016Filter == False or entry.Flag_HBHENoiseFilter  == False or entry.Flag_HBHENoiseIsoFilter  == False or entry.Flag_EcalDeadCellTriggerPrimitiveFilter  == False or  entry.Flag_BadPFMuonFilter == False  or entry.Flag_eeBadScFilter == False or entry.Flag_ecalBadCalibFilter == False) : METfilter = True

    #metdz = True
    #try : metdz = entry.Flag_BadPFMuonDzFilter 
    #except AttributeError : metdz = True

    #if not metdz : METfilter = True

    #print metdz, METfilter, 'filters', entry.Flag_goodVertices, entry.Flag_globalSuperTightHalo2016Filter, entry.Flag_HBHENoiseFilter,  entry.Flag_HBHENoiseIsoFilter,entry.Flag_EcalDeadCellTriggerPrimitiveFilter, entry.Flag_BadPFMuonFilter,  entry.Flag_eeBadScFilter,entry.Flag_ecalBadCalibFilter, entry.Flag_BadPFMuonDzFilter
    return METfilter

def printEvent(entry) :

    print("** Run={0:d} LS={1:d} Event={2:d} MET={3:.1f}".format(entry.run,entry.luminosityBlock,entry.event,entry.MET_pt))
    if entry.nMuon > 0 :
        print("Muons\n # Q    Pt   Eta   Phi   Iso  tkRelIso Medium Tight Soft    dxy     dz   MC     dR     Pt   eta   phi, genMatch")
        for j in range(entry.nMuon) :
            muSign = '+'
            if entry.Muon_charge[j] < 0 : muSign = '-'
            try : 
		print("{0:2d} {1:2s}{2:5.1f}{3:6.2f}{4:6.2f}{12:6.2f}{5:7.3f} {6:5s} {7:5s} {8:5s}{9:7.3f}{10:7.3f}{11:s}".format(
		    #j,muSign,entry.Muon_pt[j],entry.Muon_eta[j],entry.Muon_phi[j],entry.Muon_pfRelIso04_all[j],str(entry.Muon_mediumId[j]),str(entry.Muon_tightId[j]),
		    j,muSign,entry.Muon_pt[j],entry.Muon_eta[j],entry.Muon_phi[j],entry.Muon_tkRelIso[j],str(entry.Muon_mediumId[j]),str(entry.Muon_tightId[j]),
		    str(entry.Muon_softId[j]),entry.Muon_dxy[j],entry.Muon_dz[j],
		    getMCmatchString(entry.Muon_eta[j],entry.Muon_phi[j],entry), entry.Muon_tkRelIso[j])), ord(entry.Muon_genPartFlav[j])
            except AttributeError :  
		print("{0:2d} {1:2s}{2:5.1f}{3:6.2f}{4:6.2f}{5:7.3f} {6:5s} {7:5s} {8:5s}".format(
		    j,muSign,entry.Muon_pt[j],entry.Muon_eta[j],entry.Muon_phi[j],entry.Muon_pfRelIso04_all[j],str(entry.Muon_mediumId[j]),str(entry.Muon_tightId[j]),
		    str(entry.Muon_softId[j]),entry.Muon_dxy[j],entry.Muon_dz[j]))

    if entry.nElectron > 0 :
        print("Electrons                           Lost  \n # Q    Pt   Eta   Phi   Iso   Qual Hits  MVA  WP90    dxy     dz   MC     dR     Pt   eta   phi genMatch cutBased ")
        
        for j in range(entry.nElectron) :
            eSign = '+'
            if entry.Electron_charge[j] < 0 : eSign = '-'
            try :
		print("{0:2d} {1:2s}{2:5.1f}{3:6.2f}{4:6.2f}{5:7.3f}{6:6d}{7:5d}{8:7.3f} {9} {10:7.3f}{11:7.3f}{12:s} {13:f}".format(j,eSign,
		  entry.Electron_pt[j],entry.Electron_eta[j],entry.Electron_phi[j],entry.Electron_miniPFRelIso_all[j],
		  entry.Electron_cutBased[j],ord(entry.Electron_lostHits[j]),entry.Electron_mvaFall17V2noIso[j],entry.Electron_mvaFall17V2noIso_WP90[j],
		  entry.Electron_dxy[j],entry.Electron_dz[j],                                             
		  getMCmatchString(entry.Electron_eta[j],entry.Electron_phi[j],entry))), ord(entry.Electron_genPartFlav[j]), entry.Electron_cutBased[j]
            except AttributeError :  
		#print("{0:2d} {1:2s} {2:5.1f} {3:6.2f} {4:6.2f} {5:7.3f}{6:7.3f} {7} {8:7.3f} {9:f} {10:f}  {11:f}".format(j,eSign, entry.Electron_pt[j], entry.Electron_eta[j], entry.Electron_phi[j], entry.Electron_pfRelIso03_all[j],
		#  entry.Electron_cutBased[j],  ord(entry.Electron_lostHits[j]), entry.Electron_dxy[j], entry.Electron_dz[j]), ord(entry.Electron_genPartFlav[j]), entry.Electron_cutBased[j])
		print("{0:2d} {1:2s} {2:5.1f} {3:6.2f} {4:6.2f} {5:7.3f} {6:7.3f} {7:f} {8:f} ".format(j,eSign, entry.Electron_pt[j], entry.Electron_pfRelIso03_all[j], ord(entry.Electron_lostHits[j]), entry.Electron_dxy[j], entry.Electron_dz[j], ord(entry.Electron_genPartFlav[j]), entry.Electron_cutBased[j]))


    #print("Lepton List\n    Pt    Eta    Phi ")
    #for lepton in makeLeptonList(entry) :
    #    print("{0:7.1f} {1:7.3f} {2:7.3f}".format(lepton[0].Pt(),lepton[0].Eta(),lepton[0].Phi()))
        
    if entry.nJet > 0 :
        print("Jets\n #   Pt   Eta   Phi  jetId btagCSVV2  MC")
        for j in range(entry.nJet) :
            print("{0:2d} {1:5.1f}{2:6.2f}{3:6.2f}{4:6d}{5:8.3f}".format(
                j,entry.Jet_pt[j],entry.Jet_eta[j],entry.Jet_phi[j],entry.Jet_jetId[j],entry.Jet_btagCSVV2[j]))
        #print("JetList\n     Pt     Eta    Phi ")
        #i = 0 
        #for jet in makeJetList(entry) :
        #    print("{0:d} {1:7.1f} {2:6.2f} {3:6.2f}".format(i,jet[0].Pt(),jet[0].Eta(),jet[0].Phi()))
        #    i += 1
        
    if False and entry.nPhoton > 0 :
        print("Photons\n # Pt   Eta   Phi ")
        for j in range(entry.nPhoton) :
            print("{0:2d} {1:5.1f}{2:6.2f}{3:6.2f}".format(
                j,entry.Photon_pt[j],entry.Photon_eta[j],entry.Photon_phi[j]))

    if False and entry.nTau > 0:
        print("Taus                                    |-Deep Tau-||-------Iso------|")
        print(" #    Pt   Eta   Phi   Mode ID   DMID    vJ  vM  vE  Raw   Chg   Neu  jetIdx antiEl antiMu  dxy     dz  idMVA   rawIso  MC, genMatch")
        for j in range(entry.nTau) :
            try:
		print("{0:2d} {1:5.1f}{2:6.2f}{3:6.2f}{4:5d}  {5:5s} {6:5s}{18:4d}{19:4d}{20:4d} {7:6.2f}{8:6.2f}{9:6.2f}{10:6d}{11:6d}{12:6d}  {13:7.3f}{14:7.3f} {15:5d} {16:8.4f} {17:6s}".format(
		    j,entry.Tau_pt[j],entry.Tau_eta[j],entry.Tau_phi[j],entry.Tau_decayMode[j],
		    str(entry.Tau_idDecayMode[j]),str(entry.Tau_idDecayModeNewDMs[j]),
		    entry.Tau_rawIso[j],entry.Tau_chargedIso[j],entry.Tau_neutralIso[j],
		    entry.Tau_jetIdx[j],ord(entry.Tau_idAntiEle[j]),ord(entry.Tau_idAntiMu[j]),
		    entry.Tau_dxy[j],entry.Tau_dz[j],ord(entry.Tau_idMVAoldDM2017v2[j]),entry.Tau_rawMVAoldDM2017v2[j],
		    getMCmatchString(entry.Tau_eta[j],entry.Tau_phi[j],entry)[0:6],
		    ord(entry.Tau_idDeepTau2017v2p1VSjet[j]),ord(entry.Tau_idDeepTau2017v2p1VSmu[j]),ord(entry.Tau_idDeepTau2017v2p1VSe[j]))), ord(entry.Tau_genPartFlav[j])
            except AttributeError :  
		print("{0:2d} {1:5.1f}{2:6.2f}{3:6.2f}{4:5d}  {5:5s} {6:5s}{18:4d}{19:4d}{20:4d} {7:6.2f}{8:6.2f}{9:6.2f}{10:6d}{11:6d}{12:6d}  {13:7.3f}{14:7.3f} {15:5d} {16:8.4f}".format(
		    j,entry.Tau_pt[j],entry.Tau_eta[j],entry.Tau_phi[j],entry.Tau_decayMode[j],
		    str(entry.Tau_idDecayMode[j]),str(entry.Tau_idDecayModeNewDMs[j]),
		    entry.Tau_rawIso[j],entry.Tau_chargedIso[j],entry.Tau_neutralIso[j],
		    entry.Tau_jetIdx[j],ord(entry.Tau_idAntiEle[j]),ord(entry.Tau_idAntiMu[j]),
		    entry.Tau_dxy[j],entry.Tau_dz[j],ord(entry.Tau_idMVAoldDM2017v2[j]),entry.Tau_rawMVAoldDM2017v2[j],
		    getMCmatchString(entry.Tau_eta[j],entry.Tau_phi[j],entry)[0:6],
		    ord(entry.Tau_idDeepTau2017v2p1VSjet[j]),ord(entry.Tau_idDeepTau2017v2p1VSmu[j]),ord(entry.Tau_idDeepTau2017v2p1VSe[j])))

    if True and entry.nTrigObj > 0 :
        trigID = { 11:"Electr", 22:"Photon", 13:"  Muon",15:"   Tau", 1:"   Jet", 6:"FatJet", 2:"   MET", 3:"    HT" , 4:"   MHT" }
        print("Trigger Objects        Trig  Filt")
        print(" #    Pt   Eta   Phi     ID  Bits")
        for j in range(entry.nTrigObj) :
            tID = trigID[entry.TrigObj_id[j]]
            if tID == "   Tau": 
                print("{0:2d} {1:5.1f}{2:6.2f}{3:6.2f} {4:8s}{5:4X}".format(j,entry.TrigObj_pt[j],entry.TrigObj_eta[j],entry.TrigObj_phi[j],
                                                                            tID,entry.TrigObj_filterBits[j]))
        
    return

def getPDG_ID() :
    return { -24:"W-", 24:'W+', 22:"gamma", 23:"Z0", 9900041:'H++', -9900041:'H--',
                 1:'d',     2:'u',     3:'s',     4:'c',      5:'b',      6:'t',
                -1:'d_bar',-2:'u_bar',-3:'s_bar',-4:'c_bar', -5:'b_bar', -6:'t_bar',
                11:'e-',   12:'nue',  13:'mu-',  14:'nu_mu', 15:'tau-',  16:'nu_tau',
               -11:'e+',  -12:'nue', -13:'mu+', -14:'numu', -15:'tau+', -16:'nu_tau',
                21:'g' ,   25:'H', 111:'pi0',  211:'pi+', -211:'pi-' }

def printGenDecayMode(entry,printOn=False) :#works only for signal MC
    PDG_ID = getPDG_ID()
    cat =''
    if printOn: print("\n Run={0:d} Event={1:d}".format(entry.run,entry.event))
    try :
        if entry.nGenPart > 0 :
            if printOn: print("    \n #  Stat  ID  Mother")
            tau_idx = []
            for j in range(entry.nGenPart) :
                pID = entry.GenPart_pdgId[j]
                if pID in PDG_ID.keys() : pID = PDG_ID[pID]
                mother = entry.GenPart_genPartIdxMother[j]
                if mother > 0 and abs(entry.GenPart_pdgId[mother]) == 9900041:#making sure the particles come from DCH
                   if abs(entry.GenPart_pdgId[j]) == 11 :
                      cat += 'e'
                      if printOn: print("{0:2d}{1:4d}  {2:6s}{3:6d}".format(j,entry.GenPart_status[j],str(pID),mother))
                   if abs(entry.GenPart_pdgId[j]) == 13 :
                      cat += 'm'
                      if printOn: print("{0:2d}{1:4d}  {2:6s}{3:6d}".format(j,entry.GenPart_status[j],str(pID),mother))
                   if abs(entry.GenPart_pdgId[j]) == 15 :#hadronic tau selection. 
                      for i in range(j+1, entry.nGenPart, 1):
                         mom = entry.GenPart_genPartIdxMother[i]
                         if mom <= 0 or abs(entry.GenPart_pdgId[mom]) != 15: continue
                         if mom in tau_idx : continue#make sure this is not previously recognised tau
                         if abs(entry.GenPart_pdgId[i])!=11 and abs(entry.GenPart_pdgId[i])!=12 and abs(entry.GenPart_pdgId[i])!=13 and abs(entry.GenPart_pdgId[i])!=14 and abs(entry.GenPart_pdgId[i])!=15 and abs(entry.GenPart_pdgId[i])!=16 and entry.GenPart_pdgId[i]!=22:
                            cat += 't'
                            tau_idx.append(mom)
                            if printOn: print("{0:2d}{1:4d}  {2:6s}{3:6d}".format(j,entry.GenPart_status[j],str(pID),mother))
                            if printOn: print entry.GenPart_pdgId[i]
                            break
    except AttributeError :
        pass
    return cat

def printGenDecayModeBkg(entry,bkg,printOn=False) :#works only for bkg MC
    PDG_ID = getPDG_ID()
    cat =''
    if printOn: print("\n Run={0:d} Event={1:d}".format(entry.run,entry.event))
    try :
        if entry.nGenPart > 0 :
            tau_idx = []
            if printOn: print("    \n #  Stat  ID  Mother")
            for j in range(entry.nGenPart) :
                pID = entry.GenPart_pdgId[j]
                if pID in PDG_ID.keys() : pID = PDG_ID[pID]
                mother = entry.GenPart_genPartIdxMother[j]
                if bkg=='TTTo2L' and not (mother > 0 and abs(entry.GenPart_pdgId[mother]) == 24): continue#making sure leptons come from W (top decay)
                if bkg=='TTW' and not (mother > 0 and abs(entry.GenPart_pdgId[mother]) == 24): continue
                if bkg=='WZ' and not (mother > 0 and (abs(entry.GenPart_pdgId[mother]) == 24 or abs(entry.GenPart_pdgId[mother]) == 23)): continue#making sure leptons come from W and Z   
                if bkg=='DY' and not (mother in [0,1] and abs(entry.GenPart_pdgId[mother]) in [23, 1,2,3,4,5,6,21]): continue#making sure leptons come Z or from virutal photons  
                if abs(entry.GenPart_pdgId[j]) == 11 :
                   cat += 'e'
                   if printOn: print("{0:2d}{1:4d}  {2:6s}{3:6d}".format(j,entry.GenPart_status[j],str(pID),mother))
                if abs(entry.GenPart_pdgId[j]) == 13 :
                   cat += 'm'
                   if printOn: print("{0:2d}{1:4d}  {2:6s}{3:6d}".format(j,entry.GenPart_status[j],str(pID),mother))
                if abs(entry.GenPart_pdgId[j]) == 15 :#hadronic tau selection. 
                   for i in range(j+1, entry.nGenPart, 1):
                      mom = entry.GenPart_genPartIdxMother[i]
                      if mom <= 0 or abs(entry.GenPart_pdgId[mom]) != 15: continue
                      if mom in tau_idx : continue#make sure this is not previously recognised tau
                      if abs(entry.GenPart_pdgId[i])!=11 and abs(entry.GenPart_pdgId[i])!=12 and abs(entry.GenPart_pdgId[i])!=13 and abs(entry.GenPart_pdgId[i])!=14 and abs(entry.GenPart_pdgId[i])!=15 and abs(entry.GenPart_pdgId[i])!=16 and entry.GenPart_pdgId[i]!=22:
                         cat += 't'
                         tau_idx.append(mom)
                         if printOn: print("{0:2d}{1:4d}  {2:6s}{3:6d}".format(j,entry.GenPart_status[j],str(pID),mother))
                         if printOn: print entry.GenPart_pdgId[i]
                         break
    except AttributeError :
        pass
    return cat


def printMC(entry) :
    PDG_ID = getPDG_ID() 
    print("\n** MC ** Run={0:d} LS={1:d} Event={2:d} MET={3:.1f}".format(entry.run,entry.luminosityBlock,entry.event,entry.MET_pt))
    try :
        if entry.nGenPart > 0 :
            print("    \n #  Stat  ID  Mass  Mother   Pt      Eta     Phi   ")
            for j in range(entry.nGenPart) :
                pID = entry.GenPart_pdgId[j]
                if pID in PDG_ID.keys() : pID = PDG_ID[pID]
                mother = entry.GenPart_genPartIdxMother[j]
                #if mother in PDG_ID.keys() : mother = PDG_ID[mother]
                print("{0:2d}{1:4d}  {2:6s}{3:6.1f} {4:6d}{5:7.1f}{6:9.2f}{7:6.2f}".format(
                    j,entry.GenPart_status[j],str(pID),entry.GenPart_mass[j],mother,entry.GenPart_pt[j],entry.GenPart_eta[j],entry.GenPart_phi[j]))
        #print("goodMC={0:s}".format(str(goodMC(entry))))    
    except AttributeError :
        pass 
    return

def hasZmumu(entry) :
    # check gen collection to see if it has a Z->mumu decay 
    if entry.nGenPart < 3 : return False
    hasMuP, hasMuM = False, False
    for j in range(entry.nGenPart) :
        k = entry.GenPart_genPartIdxMother[j]
        if k < 1 : continue
        if entry.GenPart_pdgId[k] == 23 : 
            if entry.GenPart_pdgId[j] ==  13 : hasMuM = True
            if entry.GenPart_pdgId[j] == -13 : hasMuP = True
            if hasMuM and hasMuP : return True
    return False 

def hasZee(entry) :
    # check gen collection to see if it has a Z->ee decay 
    if entry.nGenPart < 3 : return False
    hasEP, hasEM = False, False
    for j in range(entry.nGenPart) :
        k = entry.GenPart_genPartIdxMother[j]
        if k < 1 : continue
        if entry.GenPart_pdgId[k] == 23 : 
            if entry.GenPart_pdgId[j] ==  11 : hasEM = True
            if entry.GenPart_pdgId[j] == -11 : hasEP = True
            if hasEM and hasEP : return True
    return False 

def getMCmatchString(eta, phi, entry) :
    jBest, smallestDeltaR = 999, 999.
    try :
        nGenPart = entry.nGenPart
    except AttributeError :
        return ' Not MC'
    
    for j in range(entry.nGenPart) :
        if entry.GenPart_status[j] != 1 : continue
        deltaR = sqrt((phi-entry.GenPart_phi[j])**2 + (eta-entry.GenPart_eta[j])**2)
        if deltaR < smallestDeltaR :
            jBest = j
            smallestDeltaR = deltaR

    PDG_ID = getPDG_ID() 
    pID = entry.GenPart_pdgId[jBest]
    if pID in PDG_ID.keys() :
        pID = PDG_ID[pID]
    else :
        pID = str(pID)
        
    return " {0:6s}{1:6.3f}{2:6.1f}{3:6.2f}{4:6.2f}".format(
                pID, smallestDeltaR,entry.GenPart_pt[jBest], entry.GenPart_eta[jBest], entry.GenPart_phi[jBest])
    if jBest == 999 : return '*'
    return '**'


def findDoubleLeptTrigger(goodLeptonList,entry,flavour,era):
    LepttrigList =[]
    nLepton = len(goodLeptonList)
    hltList = []
    hltListSubL = []
    leadL = -1
    subleadL = -1

    doubleLep = False
    singleLep1 = False
    singleLep2 = False
    isLfired = False
    issubLfired = False
    


    if 'ee' in flavour and nLepton > 1 :
        try : HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = entry.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ 
        except AttributeError : HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = False 

        try : HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = entry.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
        except AttributeError : HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = False
        
        try : HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL = entry.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL  
        except AttributeError : HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL = False 
        
        if era == '2016' and  not HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ : return LepttrigList, hltList
        if era != '2016'  and not HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ and not HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL  :
            return LepttrigList, hltList
     
        if entry.Electron_pt[goodLeptonList[0]] > entry.Electron_pt[goodLeptonList[1]]: 
            leadL = goodLeptonList[0]
            subleadL = goodLeptonList[1]
        else : 
            leadL = goodLeptonList[1]
            subleadL = goodLeptonList[0]

	if entry.Electron_pt[leadL] < 25 or entry.Electron_pt[subleadL] < 14 : return LepttrigList, hltList

    #if flavour == 'ee' :print 'pT ', entry.Electron_pt[leadL], entry.Electron_pt[subleadL], leadL, subleadL

    if  'mm' in flavour and nLepton > 1 :
        try : HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = entry.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
        except AttributeError : HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = False 
        try : HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ = entry.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ
        except AttributeError : HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ = False
        try : HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_Mass8 = entry.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_Mass8
        except AttributeError : HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_Mass8 = False  
        
        if era == '2016' and not  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ and not HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ : return LepttrigList, hltList
        if era != '2016'  and not HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_Mass8:  return LepttrigList, hltList
                                  

        if entry.Muon_pt[goodLeptonList[0]] > entry.Muon_pt[goodLeptonList[1]]: 
            leadL = goodLeptonList[0]
            subleadL = goodLeptonList[1]

        else : 
            leadL = goodLeptonList[1]
            subleadL = goodLeptonList[0]

        if entry.Muon_pt[leadL] < 19 or entry.Muon_pt[subleadL] < 10 : return LepttrigList, hltList
    #for electron triggers 1 = CaloIdL_TrackIdL_IsoVL, 2 = 1e (WPTight), 4 = 1e (WPLoose), 8 = OverlapFilter PFTau, 16 = 2e, 32 = 1e-1mu, 64 = 1e-1tau, 128 = 3e, 256 = 2e-1mu, 512 = 1e-2mu, 1024 = 1e (32_L1DoubleEG_AND_L1SingleEGOr), 2048 = 1e (CaloIdVT_GsfTrkIdT), 4096 = 1e (PFJet), 8192 = 1e (Photon175_OR_Photon200) for Electron (PixelMatched e/gamma)

    #for 2017, 2018 according to nAOD documentation  qualityBitsDoc = cms.string("1 = TrkIsoVVL, 2 = Iso, 4 = OverlapFilter PFTau, 8 = 1mu, 16 = 2mu, 32 = 1mu-1e, 64 = 1mu-1tau, 128 = 3mu, 256 = 2mu-1e, 512 =1mu-2e"),
    ## for 2016 in particular  "1 = TrkIsoVVL, 2 = Iso, 4 = OverlapFilter PFTau, 8 = IsoTkMu"
    dR=100.
    dRr=100.
    i_lead = -1
    i_trail = -1

    for iobj in range(0,entry.nTrigObj) :
        if 'ee' in flavour and abs(entry.TrigObj_id[iobj]) == 11 : 
	    dR = DRobj(entry.Electron_eta[leadL],entry.Electron_phi[leadL], entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj])
	    if dR  < 0.5 and entry.TrigObj_filterBits[iobj] & 1 : 
		hltList.append("LeadDEle")
		i_lead = iobj


            for iobjj in range(iobj,entry.nTrigObj) :
	        dRr = DRobj(entry.Electron_eta[subleadL],entry.Electron_phi[subleadL], entry.TrigObj_eta[iobjj], entry.TrigObj_phi[iobjj])

		if dRr  < 0.5 and entry.TrigObj_filterBits[iobjj] & 1 : 
		    hltList.append("TrailDEle")
		    i_trail = iobjj


	    if i_lead != i_trail and i_lead != -1 and i_trail != -1 : break
	


	if 'mm' in flavour and abs(entry.TrigObj_id[iobj]) == 13 : 

	    dR = DRobj(entry.Muon_eta[leadL],entry.Muon_phi[leadL], entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj])

            if era == '2016' :
		if dR  < 0.5 and entry.TrigObj_filterBits[iobj] & 1 : 
		    hltList.append("LeadDMu")
		    i_lead = iobj
            if era != '2016' :
		if dR  < 0.5 and entry.TrigObj_filterBits[iobj] & 16 : 
		    hltList.append("LeadDMu")
		    i_lead = iobj


            for iobjj in range(iobj,entry.nTrigObj) :
	        dRr = DRobj(entry.Muon_eta[subleadL],entry.Muon_phi[subleadL], entry.TrigObj_eta[iobjj], entry.TrigObj_phi[iobjj])

                if era == '2016' : 
		    if dRr  < 0.5 and entry.TrigObj_filterBits[iobjj] & 1 : 
			hltList.append("TrailDMu")
			i_trail = iobjj

                if era != '2016' : 
		    if dRr  < 0.5 and entry.TrigObj_filterBits[iobjj] & 16 : 
			hltList.append("TrailDMu")
			i_trail = iobjj

	    if i_lead != i_trail and i_lead != -1 and i_trail != -1 : break

			
    #if 'ee' in flavour  : print '=============', dR, dRr, i_trail, i_lead, hltList
    if i_lead != i_trail and i_lead != -1 and i_trail != -1  : 
         LepttrigList.append(leadL)
	 LepttrigList.append(subleadL)
	 hltList.append('BothLept')
    
    return LepttrigList, hltList



def findSingleLeptTrigger(goodLeptonList,entry,flavour,era, printOn=False):
    LepttrigList =[]
    nLepton = len(goodLeptonList)
    hltList = []
    hltListSubL = []
    objList=[]
    leadL = -1
    subleadL = -1

    doubleLep = False
    singleLep1 = False
    singleLep2 = False
    isLfired = False
    issubLfired = False
    
    HLT_IsoMu24= False
    HLT_IsoMu27= False
    HLT_IsoMu22= False
    HLT_IsoMu22_eta2p1= False
    HLT_IsoTkMu22= False
    HLT_IsoTkMu24= False
    HLT_IsoTkMu22_eta2p1= False
    HLT_Ele25_eta2p1_WPTight_Gsf = False
    HLT_Ele27_eta2p1_WPTight_Gsf = False
    HLT_Ele27_WPTight_Gsf = False
    HLT_Ele32_WPTight_Gsf = False
    HLT_Ele35_WPTight_Gsf = False


    if ('ee' in flavour or flavour=='enu') and nLepton > 0 :


        try : HLT_Ele25_eta2p1_WPTight_Gsf = entry.HLT_Ele25_eta2p1_WPTight_Gsf
        except AttributeError : HLT_Ele25_eta2p1_WPTight_Gsf = False
        try : HLT_Ele27_eta2p1_WPTight_Gsf = entry.HLT_Ele27_eta2p1_WPTight_Gsf
        except AttributeError : HLT_Ele27_eta2p1_WPTight_Gsf = False 
        try : HLT_Ele32_WPTight_Gsf = entry.HLT_Ele32_WPTight_Gsf
        except AttributeError : HLT_Ele32_WPTight_Gsf = False  
        try : HLT_Ele35_WPTight_Gsf = entry.HLT_Ele35_WPTight_Gsf
        except AttributeError : HLT_Ele35_WPTight_Gsf = False
        
        try : HLT_Ele27_WPTight_Gsf = entry.HLT_Ele27_WPTight_Gsf
        except AttributeError : HLT_Ele27_WPTight_Gsf = False

        #if era == '2016' and not HLT_Ele25_eta2p1_WPTight_Gsf and not HLT_Ele27_eta2p1_WPTight_Gsf : return LepttrigList, hltList
        #if era != '2016' and not HLT_Ele32_WPTight_Gsf and not HLT_Ele35_WPTight_Gsf :  return LepttrigList, hltList

        if '2016' in era and not HLT_Ele25_eta2p1_WPTight_Gsf : return LepttrigList, hltList, hltListSubL
        if '2017' in era and not HLT_Ele35_WPTight_Gsf :  return LepttrigList, hltList, hltListSubL
        if '2018' in era and not HLT_Ele35_WPTight_Gsf :  return LepttrigList, hltList, hltListSubL
        if 'ee' in flavour and nLepton<2 : return LepttrigList, hltList, hltListSubL
        if 'enu' in flavour and nLepton<1 : return LepttrigList, hltList, hltListSubL
	    
	#if era == '2016' and entry.Electron_pt[goodLeptonList[0]] < 29 and entry.Electron_pt[goodLeptonList[1]] < 29 : return LepttrigList, hltList
	#if era != '2016' and entry.Electron_pt[goodLeptonList[0]] < 37 and entry.Electron_pt[goodLeptonList[1]] < 37 : return LepttrigList, hltList

        


    #if flavour == 'ee' :print 'pT ', entry.Electron_pt[leadL], entry.Electron_pt[subleadL], leadL, subleadL

    if  ( flavour=='mm'  or flavour=='mnu') and nLepton > 0 :


        try : HLT_IsoMu24 = entry.HLT_IsoMu24
        except AttributeError : HLT_IsoMu24 = False
        try : HLT_IsoMu27 = entry.HLT_IsoMu27
        except AttributeError : HLT_IsoMu27 = False 

        try : HLT_IsoMu22 = entry.HLT_IsoMu22
        except AttributeError : HLT_IsoMu22 = False

        try : HLT_IsoMu22_eta2p1= entry.HLT_IsoMu22_eta2p1
        except AttributeError : HLT_IsoMu22_eta2p1 = False

        try : HLT_IsoTkMu22 = entry.HLT_IsoTkMu22
        except AttributeError : HLT_IsoTkMu22 = False

        try : HLT_IsoTkMu24 = entry.HLT_IsoTkMu24
        except AttributeError : HLT_IsoTkMu24 = False

        try : HLT_IsoTkMu22_eta2p1 = entry.HLT_IsoTkMu22_eta2p1
        except AttributeError : HLT_IsoTkMu22_eta2p1 = False



        #print ord(entry.Tau_idDeepTau2017v2p1VSmu[j]), ord(entry.Tau_idDeepTau2017v2p1VSe[j]), ord(entry.Tau_idDeepTau2017v2p1VSjet[j])
        #if era == '2016' and not HLT_IsoMu22 and not HLT_IsoMu22_eta2p1 and not HLT_IsoTkMu22 and not HLT_IsoTkMu22_eta2p1 and not HLT_IsoMu24 and not HLT_IsoTkMu24:  return LepttrigList, hltList, hltListSubL
        if '2016' in era and not HLT_IsoMu24 and not HLT_IsoTkMu24:  return LepttrigList, hltList, hltListSubL
        #if era != '2016' and not HLT_IsoMu24 and not HLT_IsoMu27 :  return LepttrigList, hltList, hltListSubL
        if '2016' not in era and not HLT_IsoMu27 and not HLT_IsoMu24:  return LepttrigList, hltList, hltListSubL
        if 'mm' in flavour and nLepton<2 : return LepttrigList, hltList, hltListSubL
        if 'mnu' in flavour and nLepton<1 : return LepttrigList, hltList, hltListSubL

    #for 2017, 2018 according to nAOD documentation  qualityBitsDoc = cms.string("1 = TrkIsoVVL, 2 = Iso, 4 = OverlapFilter PFTau, 8 = 1mu, 16 = 2mu, 32 = 1mu-1e, 64 = 1mu-1tau, 128 = 3mu, 256 = 2mu-1e, 512 =1mu-2e"),
    ## for 2016 in particular  "1 = TrkIsoVVL, 2 = Iso, 4 = OverlapFilter PFTau, 8 = IsoTkMu"

    #single mu 2016: HLT IsoMu22 v, HLT IsoMu22 eta2p1 v, HLT IsoTkMu22 v, HLT IsoTkMu22 eta2p1 v and cut pt(mu)>23, eta(mu)<2.1
    #single mu 2017: HLT IsoMu24 v, HLT IsoMu27 v and cut pt(mu)>25, eta(mu)<2.4
    #single mu 2018: HLT IsoMu24 v, HLT IsoMu27 v and cut pt(mu)>25, eta(mu)<2.4
    dR=100.
    dRr=100.
    hltList=[]
    hltListSubL=[]

    #lumiss=['509','1315','779','248','63','69']
    #if str(entry.luminosityBlock) in lumiss :  printOn=True

    for iobj in range(0,entry.nTrigObj) :
	dR=100.
	dRr=100.
	isbit2 = False
	isbit8 = False
        

        if ('ee' in flavour ) and abs(entry.TrigObj_id[iobj]) == 11 : 
	    if entry.Electron_pt[goodLeptonList[0]] > entry.Electron_pt[goodLeptonList[1]]: 
		leadL = goodLeptonList[0]
		subleadL = goodLeptonList[1]
	    else : 
		leadL = goodLeptonList[1]
		subleadL = goodLeptonList[0]

	    dR = DRobj(entry.Electron_eta[leadL],entry.Electron_phi[leadL], entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj])
	    dRr = DRobj(entry.Electron_eta[subleadL],entry.Electron_phi[subleadL], entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj])

	    if entry.TrigObj_filterBits[iobj] &2 > 0:  isbit2 = True
	    if entry.TrigObj_filterBits[iobj] &8 > 0 :  isbit8 = True

	if ('mm' in flavour )and abs(entry.TrigObj_id[iobj]) == 13 : 
	    if entry.Muon_pt[goodLeptonList[0]] > entry.Muon_pt[goodLeptonList[1]]: 
		leadL = goodLeptonList[0]
		subleadL = goodLeptonList[1]
	    else : 
		leadL = goodLeptonList[1]
		subleadL = goodLeptonList[0]

	    dR = DRobj(entry.Muon_eta[leadL],entry.Muon_phi[leadL], entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj])
	    dRr = DRobj(entry.Muon_eta[subleadL],entry.Muon_phi[subleadL], entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj])


	    if entry.TrigObj_filterBits[iobj] &2 > 0:  isbit2 = True
	    if entry.TrigObj_filterBits[iobj] &8 > 0 :  isbit8 = True

        
	    if printOn : 
		print ''
		print entry.luminosityBlock, entry.run, entry.event
	    #print("mm, iobj={7:d}, nTrigObj_id={0:d}, filter_bit={1:x}, dR_leading={2:f}, dR_subleading={3:f}, Muon_pT={4:f}, Muon_eta={5:f},  Muon_phi={6:f}, isbit2={8:b} isbit8={9:b}".format(entry.TrigObj_id[iobj], entry.TrigObj_filterBits[iobj], dR, dRr, entry.Muon_pt[leadL], abs(entry.Muon_eta[leadL]), entry.Muon_phi[leadL], iobj, isbit2, isbit8))
	    #print 'HLT_? ', HLT_IsoMu24, HLT_IsoMu27, entry.Muon_pt[leadL], abs(entry.Muon_eta[leadL]), isbit2


        if ('enu' in flavour ) and abs(entry.TrigObj_id[iobj]) == 11 : 
	    leadL = goodLeptonList[0]

	    dR = DRobj(entry.Electron_eta[leadL],entry.Electron_phi[leadL], entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj])
	    dRr = dR

	    if entry.TrigObj_filterBits[iobj] &2 > 0:  isbit2 = True
	    if entry.TrigObj_filterBits[iobj] &8 > 0 :  isbit8 = True

	if 'mnu' in flavour and abs(entry.TrigObj_id[iobj]) == 13 : 
            leadL = goodLeptonList[0]

	    dR = DRobj(entry.Muon_eta[leadL],entry.Muon_phi[leadL], entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj])
	    dRr = dR 

 	    if entry.TrigObj_filterBits[iobj] &2 > 0:  isbit2 = True
	    if entry.TrigObj_filterBits[iobj] &8 > 0 :  isbit8 = True

        
	    if printOn : 
		print ''
		print entry.luminosityBlock, entry.run, entry.event
	    #print("mm, iobj={7:d}, nTrigObj_id={0:d}, filter_bit={1:x}, dR_leading={2:f}, dR_subleading={3:f}, Muon_pT={4:f}, Muon_eta={5:f},  Muon_phi={6:f}, isbit2={8:b} isbit8={9:b}".format(entry.TrigObj_id[iobj], entry.TrigObj_filterBits[iobj], dR, dRr, entry.Muon_pt[leadL], abs(entry.Muon_eta[leadL]), entry.Muon_phi[leadL], iobj, isbit2, isbit8))
	    #print 'HLT_? ', HLT_IsoMu24, HLT_IsoMu27, entry.Muon_pt[leadL], abs(entry.Muon_eta[leadL]), isbit2

	if dR < 0.5 : 


	    if '2016' in era and (flavour == 'mm' or flavour == 'mnu'): 

                #        sel.qualityBitsDoc = cms.string("1 = TrkIsoVVL, 2 = Iso, 4 = OverlapFilter PFTau, 8 = IsoTkMu, 1024 = 1mu (Mu50)")
                if entry.Muon_pt[leadL] > 26 and abs(entry.Muon_eta[leadL]) < 2.4:
                    if printOn : 
			print ''
			print entry.luminosityBlock, entry.run, entry.event
			print("mm, iobj={7:d}, nTrigObj_id={0:d}, filter_bit={1:x}, dR_leading={2:f}, dR_subleading={3:f}, Muon_pT={4:f}, Muon_eta={5:f},  Muon_phi={6:f}, isbit2={8:b} isbit8={9:b}".format(entry.TrigObj_id[iobj], entry.TrigObj_filterBits[iobj], dR, dRr, entry.Muon_pt[leadL], abs(entry.Muon_eta[leadL]), entry.Muon_phi[leadL], iobj, isbit2, isbit8))
			#print 'HLT_IsoMu22:', HLT_IsoMu22, 'HLT_IsoMu22_eta2p1:', HLT_IsoMu22_eta2p1, 'HLT_IsoTkMu22:', HLT_IsoTkMu22, 'HLT_IsoTkMu22_eta2p1:', HLT_IsoTkMu22_eta2p1   
			#print 'HLT_IsoMu24:', HLT_IsoMu24, 'HLT_IsoTkMu24:', HLT_IsoTkMu24
                    #printTriggerObjects(entry)
                    '''
		    if (HLT_IsoMu22  and isbit2 )  :
			hltList.append(True)
                        if printOn: print 'HLT_IsoMu22:', HLT_IsoMu22

		    if (HLT_IsoMu22_eta2p1 and isbit2 and abs(entry.Muon_eta[leadL]) < 2.1)  :
			hltList.append(True)
                        if printOn: print 'HLT_IsoMu22_eta2p1:', HLT_IsoMu22_eta2p1
                    '''
		    if  HLT_IsoMu24 and isbit2 :
			hltList.append(True)
                        if printOn: print 'HLT_IsoMu24:', HLT_IsoMu24

		    if  HLT_IsoTkMu24 and isbit8 :
			hltListSubL.append(True)
                        if printOn: print 'HLT_IsoTkMu24:', HLT_IsoTkMu24

                    
		    #if  (HLT_IsoTkMu22_eta2p1 and isbit8 and abs(entry.Muon_eta[leadL]) < 2.1) :
	            #  	 hltList.append(True)
                    #    if printOn: print 'HLT_IsoTkMu22_eta2p1:', HLT_IsoTkMu22_eta2p1
                    #print 'HLT_IsoMu24:', HLT_IsoMu24, 'HLT_IsoTkMu24:', HLT_IsoTkMu24, hltList, hltListSubL

	    if  '2016' not in era and (flavour == 'mm' or flavour == 'mnu'): 

		if printOn : 
		    print ''
		    print entry.luminosityBlock, entry.run, entry.event
		    print("mm, iobj={7:d}, nTrigObj_id={0:d}, filter_bit={1:x}, dR_leading={2:f}, dR_subleading={3:f}, Muon_pT={4:f}, Muon_eta={5:f},  Muon_phi={6:f}, isbit2={8:b} isbit8={9:b}".format(entry.TrigObj_id[iobj], entry.TrigObj_filterBits[iobj], dR, dRr, entry.Muon_pt[leadL], abs(entry.Muon_eta[leadL]), entry.Muon_phi[leadL], iobj, isbit2, isbit8))
                    print 'HLT_? ', HLT_IsoMu27, entry.Muon_pt[leadL], abs(entry.Muon_eta[leadL]), isbit2, isbit8

		#if  HLT_IsoMu24  and entry.Muon_pt[leadL] > 25 and abs(entry.Muon_eta[leadL]) < 2.4 and isbit2 :
		#    hltList.append(True)

		if  HLT_IsoMu24 and entry.Muon_pt[leadL] > 26 and abs(entry.Muon_eta[leadL]) < 2.4 and isbit2 and isbit8:
		    hltList.append(True)

		if  HLT_IsoMu27 and entry.Muon_pt[leadL] > 29 and abs(entry.Muon_eta[leadL]) < 2.4 and isbit2 and isbit8:
		    hltListSubL.append(True)


	    if '2016' in era and (flavour == 'ee' or flavour=='enu'): 
                if printOn:
		    print ''
		    print entry.luminosityBlock, entry.run, entry.event

		    #print("ee, nTrigObj_id={0:d}, filter_bit={1:x}, dR_leading={2:f}, dR_subleading={3:f}, Electron_pT={4:f}, Electron_eta={5:f},  Electron_phi={9:f}, iobj={6:d}, obj_pt={10:f}, obj_eta={7:f}, obj_phi={8:f}".format(entry.TrigObj_id[iobj], entry.TrigObj_filterBits[iobj], dR, dRr, entry.Electron_pt[leadL], abs(entry.Electron_eta[leadL]), iobj, entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj], entry.Electron_phi[leadL], entry.TrigObj_pt[iobj])) 
		    print("ee, iobj={7:d}, nTrigObj_id={0:d}, filter_bit={1:x}, dR_leading={2:f}, dR_subleading={3:f}, Electron_pT={4:f}, Electron_eta={5:f},  Electron_phi={6:f}, , isbit2={8:b} isbit8={9:b}".format(entry.TrigObj_id[iobj], entry.TrigObj_filterBits[iobj], dR, dRr, entry.Electron_pt[subleadL], abs(entry.Electron_eta[subleadL]), entry.Electron_phi[subleadL], iobj, isbit2, isbit8)) 
                if HLT_Ele25_eta2p1_WPTight_Gsf and abs(entry.Electron_eta[leadL]) < 2.1 and entry.Electron_pt[leadL] > 27 :
		    if isbit2 : 
                        hltList.append(True)
                        if printOn:print 'HLT_Ele25_eta2p1_WPTight_Gsf: ',HLT_Ele25_eta2p1_WPTight_Gsf

	    if '2016' not in era and (flavour == 'ee' or flavour == 'enu'): 

                #if (HLT_Ele27_WPTight_Gsf) and abs(entry.Electron_eta[leadL]) < 2.1 and entry.Electron_pt[leadL] > 28 and isbit2 :
		#    hltList.append(True)
                #if (HLT_Ele32_WPTight_Gsf) and abs(entry.Electron_eta[leadL]) < 2.5 and entry.Electron_pt[leadL] > 34 and isbit2 :
		#    hltList.append(True)
                if (HLT_Ele35_WPTight_Gsf) and abs(entry.Electron_eta[leadL]) < 2.5 and entry.Electron_pt[leadL] > 37 and isbit2 :
		    hltList.append(True)


	if dRr < 0.5 : 


	    if '2016' in era and flavour == 'mm': 
                if entry.Muon_pt[subleadL] > 23 and abs(entry.Muon_eta[subleadL])<2.4:
                    #print("mm, subL nTrigObj_id={0:d}, filter_bit={1:x}, dR_leading={2:f}, dR_subleading={3:f}, Muon_pT={4:f}, Muon_eta={5:f},  Muon_phi={9:f}, iobj={6:d}, obj_pt={10:f}, obj_eta={7:f}, obj_phi={8:f}".format(entry.TrigObj_id[iobj], entry.TrigObj_filterBits[iobj], dR, dRr, entry.Muon_pt[subleadL], abs(entry.Muon_eta[subleadL]), iobj, entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj], entry.Muon_phi[subleadL], entry.TrigObj_pt[iobj])) 
                    if printOn: 
                        print("mm, subL iobj={7:d}, nTrigObj_id={0:d}, filter_bit={1:x}, dR_leading={2:f}, dR_subleading={3:f}, Muon_pT={4:f}, Muon_eta={5:f},  Muon_phi={6:f}, isbit2={8:b} isbit8={9:b}".format(entry.TrigObj_id[iobj], entry.TrigObj_filterBits[iobj], dR, dRr, entry.Muon_pt[subleadL], abs(entry.Muon_eta[subleadL]), entry.Muon_phi[subleadL], iobj, isbit2, isbit8)) 
                        print 'HLT_IsoMu22:', HLT_IsoMu22, 'HLT_IsoMu22_eta2p1:', HLT_IsoMu22_eta2p1, 'HLT_IsoTkMu22:', HLT_IsoTkMu22, 'HLT_IsoTkMu22_eta2p1:', HLT_IsoTkMu22_eta2p1   
		    if (HLT_IsoMu22  and isbit2)  :
			hltListSubL.append(True)
                        if printOn: print 'subL HLT_IsoMu22:', HLT_IsoMu22, 'dR', dRr

		    if (HLT_IsoMu22_eta2p1 and isbit2 and abs(entry.Muon_eta[subleadL])<2.1)  :
			hltListSubL.append(True)
                        if printOn: print 'subL HLT_IsoMu22_eta2p1:', HLT_IsoMu22_eta2p1, 'dR', dRr

		    if  (HLT_IsoTkMu22 and isbit8) :
			hltListSubL.append(True)
                        if printOn: print 'subL HLT_IsoTkMu22:', HLT_IsoTkMu22, 'dR', dRr

		    if  (HLT_IsoTkMu22_eta2p1 and isbit8 and abs(entry.Muon_eta[subleadL])<2.1) :
			hltListSubL.append(True)
                        if printOn: print 'subL HLT_IsoTkMu22_eta2p1:', HLT_IsoTkMu22_eta2p1, 'dR', dRr


	    if '2016' not in era and flavour == 'mm':
		if printOn : 
		    print ''
		    print entry.luminosityBlock, entry.run, entry.event
		    print("mm sub, iobj={7:d}, nTrigObj_id={0:d}, filter_bit={1:x}, dR_leading={2:f}, dR_subleading={3:f}, Muon_pT={4:f}, Muon_eta={5:f},  Muon_phi={6:f}, isbit2={8:b} isbit8={9:b}".format(entry.TrigObj_id[iobj], entry.TrigObj_filterBits[iobj], dR, dRr, entry.Muon_pt[subleadL], abs(entry.Muon_eta[subleadL]), entry.Muon_phi[subleadL], iobj, isbit2, isbit8))
                    print 'HLT_? sub', HLT_IsoMu24, HLT_IsoMu27, entry.Muon_pt[subleadL], abs(entry.Muon_eta[subleadL]), isbit2


		if  HLT_IsoMu24  and entry.Muon_pt[subleadL] > 25 and abs(entry.Muon_eta[subleadL]) < 2.4 and isbit2 and isbit8:
		    hltListSubL.append(True)

		if  HLT_IsoMu27 and entry.Muon_pt[subleadL] > 25 and abs(entry.Muon_eta[subleadL]) < 2.4 and isbit2 and isbit8:
		    hltListSubL.append(True)


	    if '2016' in era and flavour == 'ee': 
                if HLT_Ele25_eta2p1_WPTight_Gsf and abs(entry.Electron_eta[subleadL]) < 2.1 and entry.Electron_pt[subleadL] > 26 and isbit2:
		    hltListSubL.append(True)
                    #print("ee, subL nTrigObj_id={0:d}, filter_bit={1:x}, dR_leading={2:f}, dR_subleading={3:f}, Electron_pT={4:f}, Electron_eta={5:f},  Electron_phi={9:f}, iobj={6:d},  obj_pt={10:f}, obj_eta={7:f}, obj_phi={8:f}".format(entry.TrigObj_id[iobj], entry.TrigObj_filterBits[iobj], dR, dRr, entry.Electron_pt[subleadL], abs(entry.Electron_eta[subleadL]), iobj, entry.TrigObj_eta[iobj], entry.TrigObj_phi[iobj], entry.Electron_phi[subleadL], entry.TrigObj_eta[iobj])) 
                    if printOn: 
                        print("ee, subL iobj={7:d}, nTrigObj_id={0:d}, filter_bit={1:x}, dR_leading={2:f}, dR_subleading={3:f}, Electron_pT={4:f}, Electron_eta={5:f},  Electron_phi={6:f}".format(entry.TrigObj_id[iobj], entry.TrigObj_filterBits[iobj], dR, dRr, entry.Electron_pt[subleadL], abs(entry.Electron_eta[subleadL]), entry.Electron_phi[subleadL], iobj)) 
                        print 'subL HLT_Ele25_eta2p1_WPTight_Gsf:', HLT_Ele25_eta2p1_WPTight_Gsf, 'dR', dRr

	    if '2017' in era and flavour == 'ee': 

                if (HLT_Ele27_WPTight_Gsf) and abs(entry.Electron_eta[subleadL]) < 2.1 and entry.Electron_pt[subleadL] > 28 and isbit2 :
		    hltListSubL.append(True)
                if (HLT_Ele32_WPTight_Gsf) and abs(entry.Electron_eta[subleadL]) < 2.1 and entry.Electron_pt[subleadL] > 28 and isbit2 :
		    hltListSubL.append(True)
                if (HLT_Ele35_WPTight_Gsf) and abs(entry.Electron_eta[subleadL]) < 2.1 and entry.Electron_pt[subleadL] > 28 and isbit2 :
		    hltListSubL.append(True)

	    if '2018' in era  and flavour == 'ee': 

                if (HLT_Ele32_WPTight_Gsf) and abs(entry.Electron_eta[subleadL]) < 2.1 and entry.Electron_pt[subleadL] > 33 and isbit2 :
		    hltListSubL.append(True)
                if (HLT_Ele35_WPTight_Gsf) and abs(entry.Electron_eta[subleadL]) < 2.1 and entry.Electron_pt[subleadL] > 33 and isbit2 :
		    hltListSubL.append(True)



    #print 'check 2 HLT_IsoMu24 finished:', HLT_IsoMu24, 'HLT_IsoTkMu24:', HLT_IsoTkMu24

    if len(hltList)>0 : LepttrigList.append(leadL)
    if len(hltListSubL)>0 : LepttrigList.append(subleadL)

    #if flavour=='ee'  :  print '--------------->',  era, flavour, LepttrigList, 'hltL', hltList, 'hltSubL', hltListSubL, entry.luminosityBlock, entry.event
    return LepttrigList, hltList, hltListSubL




def DRobj(eta1,phi1,eta2,phi2) :
    dPhi = min(abs(phi2-phi1),2.*pi-abs(phi2-phi1))
    return sqrt(dPhi**2 + (eta2-eta1)**2)


class cutCounter() :

    def __init__(self):
        self.counter = {}
        self.counterGenWeight = {}
        self.nickNames = []
        self.yields = []
        self.labels = []

    def count(self,nickName) :
        try :
            self.counter[nickName] += 1
        except KeyError :
            self.nickNames.append(nickName) 
            self.counter[nickName] = 1

    def countGenWeight(self,nickName,w) :
        try :
            self.counterGenWeight[nickName] += float(w)
        except KeyError :
            self.nickNames.append(nickName) 
            self.counterGenWeight[nickName] = 1
            
    def printSummary(self) :
        #print("Cut summary:\n    Name      Events Fraction")
        nLast = 0.
        for nn in self.nickNames :
            fraction = 1.0
            if nLast > 0. : fraction = self.counter[nn]/nLast
            nLast = float(self.counter[nn])
            sFrac = '   N/A'
            if fraction < 1.0 : sFrac = "{0:6.1f}%".format(100.*fraction)
            if fraction < 0.01 : sFrac = "{0:6.3f}%".format(100.*fraction)
            print("{0:16s}{1:6d} {2:s}".format(nn,self.counter[nn],sFrac))

        return

    def getYield(self) :
        for nn in self.nickNames :
            self.yields.append(self.counter[nn])

        return self.yields

    def getYieldWeighted(self) :
        for nn in self.nickNames :
            self.yields.append(self.counterGenWeight[nn])

        return self.yields

    def getLabels(self) :
        for nn in self.nickNames :
            self.labels.append(nn)

        return self.labels

    def writeCSV(self,args) :
        outLines = [] 
        for nn in self.nickNames :
            outLines.append("{0:s},{1:d}\n".format(nn,self.counter[nn]))

        if len(args.csvFileName) > 1 :
            CSVfile = args.csvFileName 
        else :
            CSVfile = "./{0:s}.csv".format(args.inFileName.strip('.root'))
            
        if 'cmseos.fnal' in args.inFileName :
            CSVfile = CSVfile.split('/')[-1].replace('.root','.csv')

        print("Writing CSV to {0:s}".format(CSVfile))
        open(CSVfile,'w').writelines(outLines)
        return

def getOutFileName(args) :
    print("In generalFunctions.getOutFileName() args.outFileName={0:s}".format(args.outFileName))
    if len(args.outFileName) > 2 : return args.outFileName
    if ('cmseos.fnal' in args.inFileName) or ('cms-xrd-global' in args.inFileName) :
        outFileBase = args.inFileName.split('/')[-1]
        return "./outData/{0:s}_out.root".format(outFileBase[:-5])
        
    return 'temp_out.root'

class dupeDetector() :
    
    def __init__(self):
        self.nCalls = 0 
        self.runEventList = []

    def checkEvent(self,entry) :
        self.nCalls += 1 
        runEvent = "{0:d}:{1:d}".format(entry.run,entry.event)
        if runEvent in self.runEventList :
            #print("Event in list: runEventList={0:s}".format(str(self.runEventList)))
            return True
        else :
            self.runEventList.append(runEvent)
            #print("New event: runEventList={0:s}".format(str(self.runEventList)))
            return False

    def printSummary(self) :
        print("Duplicate Event Summary: Calls={0:d} Unique Events={1:d}".format(self.nCalls,len(self.runEventList)))
        return


class pileUpWeight() :
    
    def __init__(self):
        self.dummy = 0
        
    def calculateWeights(self,nickName,year) :
        # get data pileup histogram
        fData = TFile('pileup_{0:s}UL_data.root'.format(str(year)))
        hData = fData.Get('pileup')
        print("hData={0:s}".format(str(hData)))
        binWidth = hData.GetBinWidth(1)
        xMin = hData.GetBinLowEdge(1)
        nBins = hData.GetNbinsX()

        # lumi placeholder values 
        lumi = { 2016:35.9, 2017:41.5, 2018:59.7 }
        xSec = 1.
        # get MC cross section values
        #for line in open('MCsamples.csv','r').readlines() :
        #for line in open('MCsamples_{0:d}.csv'.format(year),'r').readlines() :
        #    if nickName == line.split(',')[0].strip() :
        #        xSec = 1000.*float(line.split(',')[2])
                 
        # get MC pileup histograms
        MCfile = "pileup_{0:s}UL_MC.root".format(str(year))
        print("Opening MC pileup file = {0:s}".format(MCfile))
        fMC = TFile(MCfile)
        print("fMC={0:s}".format(str(fMC)))

	#temp hack given that we dont have all histos
        #hMC = fData.Get('pileup')

        #hMC = fMC.Get('hMC_{0:s}'.format(nickName)) ##this is you would need to have one histo per process ## old naming , new is the Pileup_nTrueInt
        hMC = fMC.Get('hPileup_nTrueInt') ##this is you would need to have one histo per process

        #hMC = fMC.Get('h{0:s}'.format(nickName)) ##this is you would need to have one histo per process
        #hMC = fMC.Get('pileup') ##this is you would need to have one histo per process

        # check to be sure that data and MC histograms are commensurate
        if hData.GetBinWidth(1) != hMC.GetBinWidth(1) or hData.GetBinLowEdge(1) != hMC.GetBinLowEdge(1) or hData.GetNbinsX() != hMC.GetNbinsX() :
            print("Error in generalFunctions.pileUpWeight().calculateWeights()\nData and MC histograms not commensurate.") 

        nData = hData.GetSumOfWeights()
        pData = np.array(hData)[1:-1]/nData
        print("sum of pData={0:f}".format(np.sum(pData)))
        nMC = hMC.GetSumOfWeights()
        pMC = np.array(hMC)[1:-1]
        pMC /= nMC 
        pMC = np.maximum(1.e-5*np.ones_like(pMC),pMC)
        print("sum of pMC={0:f}".format(np.sum(pMC)))
        weights = np.divide(pData,pMC)
        self.PUweights = weights
	#print 'inside', weights, len(weights)
	#print '========', self.PU, hMC.FindBin(PU), hData.FindBin(PU), hData.GetBinContent(hData.FindBin(PU))/hMC.GetBinContent(hMC.FindBin(PU)), 'is it the same?', weights
        xMin = hData.GetBinLowEdge(1)
        xMax = xMin + hData.GetNbinsX()*hData.GetBinWidth(1) 
        bins = np.linspace(xMin+0.5*binWidth,xMax-0.5*binWidth,nBins)
        #self.sampleWeight = xSec*lumi[year]/nMC
        print("In generalFunctions.pileUpWeight.calculateWeights() :")
        #print(" nickName={0:s} year={1:d} lumi={2:.1f} /fb xSec={3:.3f} fb nMC={4:.1f} weight={5:f}".format(nickName,year,lumi[year],xSec,nMC,self.sampleWeight))
        
        if False :
            gData = TGraph(len(bins),bins,pData)
            gData.GetXaxis().SetTitle("PileUp") 
            gData.SetMarkerColor(kRed)
            gData.SetMarkerSize(1.0)
            gData.SetMarkerStyle(21)
            gMC = TGraph(len(bins),bins,pMC) 
            c1 = TCanvas("c1","c1",1000,750)
            gData.Draw('AP')
            gMC.SetMarkerColor(kBlue)
            gMC.SetMarkerSize(1.0)
            gMC.SetMarkerStyle(22)
            gMC.Draw('P')
            legend = TLegend(0.6,0.8,0.80,0.90);
            legend.AddEntry(gData,"Data") 
            legend.AddEntry(gMC,"MC") 
            legend.Draw()
            c1.Draw()
            #raw_input()
            
        return bins, weights

    def getWeight(self,PU) :
        iPU = min(99,int(PU))
        #weight = self.sampleWeight*self.PUweights[iPU]
        weight = self.PUweights[iPU] #the sampleWeight should be taken from the MC later on
        #print 'weights', iPU, weight, self.sampleWeight, self.PUweights[iPU], 'old', self.sampleWeight*self.PUweights[iPU]
        return weight 
        
    def displayWeights(self, bins, weights) :
        gWeights = TGraph(len(bins),bins,weights) 
        c1 = TCanvas("c1","c1",1000,750)
        gWeights.GetXaxis().SetTitle("PileUp")
        gWeights.GetYaxis().SetTitle("Weight")
        gWeights.SetLineWidth(2)
        gWeights.SetMarkerColor(kBlue)
        gWeights.SetMarkerSize(1.0)
        gWeights.SetMarkerStyle(22)
        gWeights.Draw("AP")
        c1.Draw()
        raw_input()
            
# find ID of first child of parent
def findFirst(e, vetoList, parent) :
    for i in range(parent+1,e.nGenPart,1) :
        if e.GenPart_genPartIdxMother[i] == parent :
            if not e.GenPart_pdgId[i] in vetoList :
                return i
    print("In generalFunction.findFirst() parent not found: parent={0:d} vetoList={1:s}".format(
        parent,str(vetoList)))
    return -1

# find last instance of a particle of a given parent  
def findLast(e, ID, parent) :
    last  = -1 
    for i in range(parent+1,e.nGenPart,1) :
        if e.GenPart_pdgId[i] == ID :
            last = i
    if last < 0 :
        print("In generalFunction.findLast() particle not found: ID={0:d} parent={1:d}".format(ID,parent))

    return last 
        
# look for ZH events
def eventID(e) :
    # find decay mode of Z0
    neutrinos = [12, 14, 16, -12, -14, -16] 
    elec, mu, tau, Z0, H = 11, 13, 15, 23, 25
    iZ0 = findLast(e,Z0,0)
    if iZ0 < 0 : return ''
    iLep = findFirst(e,[],iZ0)
    if iLep < 0 : return ''
    lepPDG = abs(e.GenPart_pdgId[iLep])
    if lepPDG == elec :
        cat = 'ee'
    elif lepPDG == mu :
        cat = 'mm'
    else :
        #print("In generalFunction.eventID() unrecognized Z0 decay mode lepPDG={0:d}".format(lepPDG))
        return ''
        
    # find decay mode of taus from Higgs
    iH = findLast(e,H,0)
    if iH < 0 : return ''
    for child in [-tau,tau] :
        iTauP = findLast(e,child,iH)
        tauChild = findFirst(e,neutrinos,iTauP)
        tauChildPDG = e.GenPart_pdgId[tauChild]
        if abs(tauChildPDG) == elec :
            cat += 'e'
        elif abs(tauChildPDG) == mu :
            cat += 'm'
        else :
            cat += 't'

    if cat[2:4] in ['ee','em','me','mm'] : return '' 
    cat = cat.replace('te','et')
    cat = cat.replace('tm','mt')
    return cat


class checkJSON() :
    
    def __init__(self,filein) :
        self.good, self.bad  = 0, 0
        print 'inside json function : will use the JSON', filein
        input_file = open (filein)
        self.json_array = json.load(input_file)  
      

    def checkJSON(self,LS,run) :
        try :
            LSlist = self.json_array[str(run)]
            for LSrange in LSlist :
		if int(LS) >= int(LSrange[0]) and int(LS) <= int(LSrange[1]) :
                    self.good += 1
                    return True
        except KeyError :
            pass
        
        self.bad += 1
        return False
        
    def printJSONsummary(self) :
        print("check JSON summary:  nCalls={0:d} nGood={1:d} nBad={2:d}".format(self.good+self.bad,self.good,self.bad))
        return
    

    
    

