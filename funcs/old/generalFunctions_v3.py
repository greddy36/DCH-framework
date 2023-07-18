# functions for H->tautau analysis 

from ROOT import TLorentzVector
from ROOT import TFile, TH1D, TCanvas, TGraph, kRed, kBlue, TLegend
from math import sqrt
import numpy as np

def printEvent(entry) :
    print("** Run={0:d} LS={1:d} Event={2:d} MET={3:.1f}".format(entry.run,entry.luminosityBlock,entry.event,entry.MET_pt))
    if entry.nMuon > 0 :
        print("Muons\n # Q    Pt   Eta   Phi   Iso  Medium Tight Soft    dxy     dz   MC     dR     Pt   eta   phi")
        for j in range(entry.nMuon) :
            muSign = '+'
            if entry.Muon_charge[j] < 0 : muSign = '-'
            print("{0:2d} {1:2s}{2:5.1f}{3:6.2f}{4:6.2f}{5:7.3f} {6:5s} {7:5s} {8:5s}{9:7.3f}{10:7.3f}{11:s}".format(
                j,muSign,entry.Muon_pt[j],entry.Muon_eta[j],entry.Muon_phi[j],entry.Muon_pfRelIso04_all[j],str(entry.Muon_mediumId[j]),str(entry.Muon_tightId[j]),
                str(entry.Muon_softId[j]),entry.Muon_dxy[j],entry.Muon_dz[j],
                getMCmatchString(entry.Muon_eta[j],entry.Muon_phi[j],entry)))
    if entry.nElectron > 0 :
        print("Electrons                           Lost  \n # Q    Pt   Eta   Phi   Iso   Qual Hits  MVA  WP90    dxy     dz   MC     dR     Pt   eta   phi")
        
        for j in range(entry.nElectron) :
            eSign = '+'
            if entry.Electron_charge[j] < 0 : eSign = '-'
            print("{0:2d} {1:2s}{2:5.1f}{3:6.2f}{4:6.2f}{5:7.3f}{6:6d}{7:5d}{8:7.3f} {9} {10:7.3f}{11:7.3f}{12:s}".format(j,eSign,
              entry.Electron_pt[j],entry.Electron_eta[j],entry.Electron_phi[j],entry.Electron_miniPFRelIso_all[j],
              entry.Electron_cutBased[j],ord(entry.Electron_lostHits[j]),entry.Electron_mvaFall17V2noIso[j],entry.Electron_mvaFall17V2noIso_WP90[j],
              entry.Electron_dxy[j],entry.Electron_dz[j],                                             
              getMCmatchString(entry.Electron_eta[j],entry.Electron_phi[j],entry)))


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

    if True and entry.nTau > 0:
        print("Taus                                    |-------Iso------|")
        print(" #    Pt   Eta   Phi   Mode ID   DMID    Raw   Chg   Neu  jetIdx antiEl antiMu  dxy     dz  idMVA   rawIso  MC")
        for j in range(entry.nTau) :
            print("{0:2d} {1:5.1f}{2:6.2f}{3:6.2f}{4:5d}  {5:5s} {6:5s} {7:6.2f}{8:6.2f}{9:6.2f}{10:6d}{11:6d}{12:6d}  {13:7.3f}{14:7.3f} {15:5d} {16:8.4f} {17:6s}".format(
                j,entry.Tau_pt[j],entry.Tau_eta[j],entry.Tau_phi[j],entry.Tau_decayMode[j],
                str(entry.Tau_idDecayMode[j]),str(entry.Tau_idDecayModeNewDMs[j]),
                entry.Tau_rawIso[j],entry.Tau_chargedIso[j],entry.Tau_neutralIso[j],
                entry.Tau_jetIdx[j],ord(entry.Tau_idAntiEle[j]),ord(entry.Tau_idAntiMu[j]),
                entry.Tau_dxy[j],entry.Tau_dz[j],ord(entry.Tau_idMVAoldDM2017v2[j]),entry.Tau_rawMVAoldDM2017v2[j],
                getMCmatchString(entry.Tau_eta[j],entry.Tau_phi[j],entry)[0:6]))
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
    return { -24:"W-", 24:'W+', 22:"gamma", 23:"Z0", 
                 1:'d',     2:'u',     3:'s',     4:'c',      5:'b',      6:'t',
                -1:'d_bar',-2:'u_bar',-3:'s_bar',-4:'c_bar', -5:'b_bar', -6:'t_bar',
                11:'e-',   12:'nue',  13:'mu-',  14:'nu_mu', 15:'tau-',  16:'nu_tau',
               -11:'e+',  -12:'nue', -13:'mu+', -14:'numu', -15:'tau+', -16:'nu_tau',
                21:'g' ,   25:'H',   211:'pi+', -211:'pi-' }

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

class cutCounter() :

    def __init__(self):
        self.counter = {}
        self.nickNames = []

    def count(self,nickName) :
        try :
            self.counter[nickName] += 1
        except KeyError :
            self.nickNames.append(nickName) 
            self.counter[nickName] = 1
            
    def printSummary(self) :
        print("Cut summary:\n    Name      Events Fraction")
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
        fData = TFile('data_pileup_{0:d}.root'.format(year))
        hData = fData.Get('pileup')
        print("hData={0:s}".format(str(hData)))
        binWidth = hData.GetBinWidth(1)
        xMin = hData.GetBinLowEdge(1)
        nBins = hData.GetNbinsX()

        # lumi placeholder values 
        lumi = { 2016:37.8, 2017:44.0, 2018:63.7 }
        # get MC cross section values
        for line in open('MCsamples.csv','r').readlines() :
            if nickName == line.split(',')[0].strip() :
                xSec = 1000.*float(line.split(',')[1])
                 
        # get MC pileup histograms
        fMC = TFile('MC_{0:d}.root'.format(year))
        hMC = fMC.Get('h{0:s}'.format(nickName))
        print("hMC={0:s}".format(str(hMC)))

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
        xMin = hData.GetBinLowEdge(1)
        xMax = xMin + hData.GetNbinsX()*hData.GetBinWidth(1) 
        bins = np.linspace(xMin+0.5*binWidth,xMax-0.5*binWidth,nBins)
        self.sampleWeight = xSec*lumi[year]/nMC
        print("In generalFunctions.pileUpWeight.calculateWeights() :")
        print(" nickName={0:s} year={1:d} lumi={2:.1f} /fb xSec={3:.3f} fb nMC={4:.1f} weight={5:f}".format(nickName,year,lumi[year],xSec,nMC,self.sampleWeight))
        
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
            raw_input()
            
        return bins, weights

    def getWeight(self,PU) :
        iPU = int(PU)
        weight = self.sampleWeight*self.PUweights[iPU]
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
            
