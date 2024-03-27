#!/usr/bin/env python


# import external modules 
import sys
import numpy as np
from ROOT import TFile, TTree, TH1, TH1D, TCanvas, TLorentzVector  
from math import sqrt, pi
#from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool

# import from ZH_Run2/funcs/
sys.path.insert(1,'../funcs/')
import tauFunDCH as TF
import generalFunctions as GF 
import Weights 
import outTuple
import time

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-v","--verbose",default=0,type=int,help="Print level.")
    parser.add_argument("-f","--inFileName",default='ZHtoTauTau.root',help="File to be analyzed.")
    parser.add_argument("-c","--category",default='none',help="Event category to analyze.")
    parser.add_argument("--nickName",default='',help="MC sample nickname") 
    parser.add_argument("-d","--dataType",default='MC',help="Data or MC") 
    parser.add_argument("-o","--outFileName",default='',help="File to be used for output.")
    parser.add_argument("-n","--nEvents",default=0,type=int,help="Number of events to process.")
    parser.add_argument("-m","--maxPrint",default=0,type=int,help="Maximum number of events to print.")
    parser.add_argument("-t","--testMode",default='',help="tau MVA selection")
    parser.add_argument("-y","--year",default=2017,type=int,help="Data taking period, 2016, 2017 or 2018")
    parser.add_argument("-s","--selection",default='ZH',help="is this for the ZH or the AZH analysis?")
    parser.add_argument("-u","--unique",default='none',help="CSV file containing list of unique events for sync studies.") 
    parser.add_argument("-w","--weights",default=False,type=int,help="to re-estimate Sum of Weights")
    parser.add_argument("-j","--doSystematics",type=str, default='false',help="do JME systematics")
    parser.add_argument("-e","--era",type=str, default='UL',help="EOY of UL")
    
    return parser.parse_args()

args = getArgs()
print("args={0:s}".format(str(args)))
maxPrint = args.maxPrint 
#file = open(args.nickName+'_'+args.category+'.txt', 'a')
#sys.stdout = file

cutCounter = {}
cutCounterGenWeight = {}
cutflow = {}
doJME  = args.doSystematics.lower() == 'true' or args.doSystematics.lower() == 'yes' or args.doSystematics == '1'

#doJME = True

#cats = ['mmtt', 'mmet', 'mmmt', 'mmem']
'''
cats_4l = ['eeee','eemm','mmmm']
cats_2l2t = ['eeee','eemm','mmmm',#tau(lep)+tau(lep)
             'eeet','eemt','mmet','mmmt',#tau(lep)+tau(had)
             'eett','mmtt'#tau(had)+tau(had)
            ]
cats_4t = ['tttt',#all hadr tau
           'ettt','mttt',#tau(lep)tau(had) + 2 tau(had)
           'eett','mmtt','emtt',#2 tau(lep) + 2 tau(had)
           'etet','etmt','mtmt',#2 tau(lep)tau(had) pairs
           'eeet','eemt','emmt','emet','mmmt','mmet',#2 tau(lep) + tau(lep)tau(had) 
           'eeee','eeem','emem','eemm','mmmm','mmem'#4 tau(lep)
           ]
cats = cats_4l+cats_2l2t+cats_4t
'''
cats = ['eeee','eeem','eeet','eemm','eemt','eett',
               'emem','emet','emmm','emmt','emtt',
                      'etet','etmm','etmt','ettt',
                             'mmmm','mmmt','mmtt',
                                    'mtmt','mttt',
                                           'tttt']#ggr

cats3L = ['eee', 'eem', 'eet',
          'eme', 'emm', 'emt',
          'ete', 'etm', 'ett',
          'mme', 'mmm', 'mmt',
          'mte', 'mtm', 'mtt',
          'tte', 'ttm', 'ttt']

for cat in cats : 
    cutCounter[cat] = GF.cutCounter()
    cutCounterGenWeight[cat] = GF.cutCounter()

cutflow['ele'] = GF.cutCounter()
cutflow['muon'] = GF.cutCounter()
cutflow['tau'] = GF.cutCounter()

inFileName = args.inFileName
print("Opening {0:s} as input.  Event category {1:s}".format(inFileName,cat))


inFile = TFile.Open(inFileName)
inFile.cd()

inTree = inFile.Get("Events")
nentries = inTree.GetEntries()
nMax = nentries
print("nentries={0:d} nMax={1:d}".format(nentries,nMax))
if args.nEvents > 0 : nMax = min(args.nEvents-1,nentries)


MC = len(args.nickName) > 0 
if args.dataType == 'Data' or args.dataType == 'data' : MC = False
if args.dataType == 'MC' or args.dataType == 'mc' : MC = True

if MC :
    #print "this is MC, will get PU etc", args.dataType
    PU = GF.pileUpWeight()
    PU.calculateWeights(args.nickName,args.year)
else :
    CJ = ''
    if args.year == 2016 : CJ = GF.checkJSON(filein='Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt')
    if args.year == 2017 : CJ = GF.checkJSON(filein='Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt')
    if args.year == 2018 : CJ = GF.checkJSON(filein='Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt')


#print 'systematics', doJME

era=str(args.year)

#outFileName = GF.getOutFileName(args).replace(".root",".ntup")
outFileName = GF.getOutFileName(args)

if MC : 
    if "WJetsToLNu" in outFileName and 'TWJets' not in outFileName:
	hWxGenweightsArr = []
	for i in range(5):
	    hWxGenweightsArr.append(TH1D("W"+str(i)+"genWeights",\
		    "W"+str(i)+"genWeights",1,-0.5,0.5))
    elif "DYJetsToLL" in outFileName:
	hDYxGenweightsArr = []
	for i in range(5):
	    hDYxGenweightsArr.append(TH1D("DY"+str(i)+"genWeights",\
		    "DY"+str(i)+"genWeights",1,-0.5,0.5))

#isZH = 'ZH' in outFileName or 'HZJ' in outFileName

if args.weights > 0 :
    hWeight = TH1D("hWeights","hWeights",1,-0.5,0.5)
    hWeight.Sumw2()
    hWeightScaleUp = TH1D("hWeightsScaleUp","hWeightsScaleUp",1,-0.5,0.5)
    hWeightScaleUp.Sumw2()
    hWeightScaleDown = TH1D("hWeightsScaleDown","hWeightsScaleDown",1,-0.5,0.5)
    hWeightScaleDown.Sumw2()


    hWeightScaleSTXS = TH1D("hWeightsScaleSTXS","hWeightsScaleSTXS",11,-0.5,10.5)
    hWeightScaleSTXS.Sumw2()
    hWeightScaleSTXSUp = TH1D("hWeightsScaleSTXSUp","hWeightsScaleSTXSUp",11,-0.5,10.5)
    hWeightScaleSTXSUp.Sumw2()
    hWeightScaleSTXSDown = TH1D("hWeightsScaleSTXSDown","hWeightsScaleSTXSDown",11,-0.5,10.5)
    hWeightScaleSTXSDown.Sumw2()



    for count, e in enumerate(inTree) :
        hWeight.Fill(0, e.genWeight)
    
        if "WJetsToLNu" in outFileName and 'TWJets' not in outFileName:

            npartons = ord(e.LHE_Njets)
	    if  npartons <= 4: 	hWxGenweightsArr[npartons].Fill(0, e.genWeight)
        if "DYJetsToLL" in outFileName :
            npartons = ord(e.LHE_Njets)
	    if  npartons <= 4 : hDYxGenweightsArr[npartons].Fill(0, e.genWeight)

    fName = GF.getOutFileName(args).replace(".root",".weights")
    fW = TFile( fName, 'recreate' )
    #print 'Will be saving the Weights in', fName
    fW.cd()

    if "WJetsToLNu" in outFileName and 'TWJets' not in outFileName:
        for i in range(len(hWxGenweightsArr)):
            hWxGenweightsArr[i].Write()
    elif "DYJetsToLL" in outFileName:
        for i in range(len(hDYxGenweightsArr)):
            hDYxGenweightsArr[i].Write()

    hWeight.Write()
    hWeightScaleUp.Write()
    hWeightScaleDown.Write()
    hWeightScaleSTXS.Write()
    hWeightScaleSTXSUp.Write()
    hWeightScaleSTXSDown.Write()
    if args.weights > 1 : 
        fW.Close()
        sys.exit()

#############end weights

# read a CSV file containing a list of unique events to be studied 
    
print("Opening {0:s} as output.".format(outFileName))

sysT = ["Central"]

sysall = ['scale_e', 'scale_m_etalt1p2', 'scale_m_eta1p2to2p1', 'scale_m_etagt2p1',
'scale_t_1prong', 'scale_t_1prong1pizero', 'scale_t_3prong', 'scale_t_3prong1pizero']


upS=sysall
downS=sysall

for i, sys in enumerate(sysall) : 
    sysT.append(sys+'Up')
    sysT.append(sys+'Down')


sysT=['Central']
#print sysT

isMC = True
if not MC : 
    sysT = ["Central"]
    isMC = False

#sysT = ["Central"]
doSyst= False
outTuple = outTuple.outTuple(outFileName, era, doSyst, sysT, isMC)



tStart = time.time()
countMod = 1000

#print outTuple.allsystMET

allMET=[]
for i,j in enumerate(outTuple.allsystMET):
    if 'MET' in j and 'T1_' in j and 'phi' not in j : allMET.append(j)



#Weights=Weights.Weights(args.year)

cat_yield = {} #empty dictoinary
n_lepton = {}
selected_evts, veto_evts , evts_3lep , evts_5lep = 0,0,0,0
pass_evts = 0
for cat in cats:
    cat_yield[cat] = 0
    n_lepton[cat]=[0,0,0]

for count, e in enumerate( inTree) :
    #if count != 738: continue #to run only over a single event
    if count % countMod == 0 :
        print("Count={0:d}".format(count))
        if count >= 10000 : countMod = 10000
    if count == nMax : break    
    printOn=False

    for cat in cats : 
        cutCounter[cat].count('All')
	if  MC :   cutCounterGenWeight[cat].countGenWeight('All', e.genWeight)
 
    isInJSON = False
    if not MC : isInJSON = CJ.checkJSON(e.luminosityBlock,e.run)
    if not isInJSON and not MC :
        #print("Event not in JSON: Run:{0:d} LS:{1:d}".format(e.run,e.luminosityBlock))
        continue

    for cat in cats: 
        cutCounter[cat].count('InJSON')
	if  MC :   cutCounterGenWeight[cat].countGenWeight('InJSON', e.genWeight)
    
    MetFilter = GF.checkMETFlags(e,args.year)
    if MetFilter : continue
    for cat in cats: 
        cutCounter[cat].count('METfilter') 
	if  MC :   cutCounterGenWeight[cat].countGenWeight('METfilter', e.genWeight)

    if not TF.goodTrigger(e,args.year) and printOn :   print cat, e.run, e.luminosityBlock,  e.event, 'Triggers not present...'
    if not TF.goodTrigger(e, args.year) : continue


    for cat in cats: 
        isTrig=False
	if cat[:2] =='ee' : isTrig = e.HLT_Ele27_WPTight_Gsf or e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ 
        if cat[:2] =='em' : isTrig = (e.HLT_Ele27_WPTight_Gsf or e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) and e.HLT_IsoMu24
        if cat[:2] =='et' : isTrig = e.HLT_Ele27_WPTight_Gsf or e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
        if cat[:2] =='mm' : isTrig = e.HLT_IsoMu24 or e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
        if cat[:2] =='mt' : isTrig = e.HLT_IsoMu24 or e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
	if isTrig : 
            cutCounter[cat].count('Lep_Trig1')
	    if  MC :   cutCounterGenWeight[cat].countGenWeight('Lep_Trig1', e.genWeight)

        isTrig=False
        if cat[2:] =='ee' : isTrig = e.HLT_Ele27_WPTight_Gsf or e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
        if cat[2:] =='em' : isTrig = (e.HLT_Ele27_WPTight_Gsf or e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) and e.HLT_IsoMu24
        if cat[2:] =='et' : isTrig = e.HLT_Ele27_WPTight_Gsf or e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
        if cat[2:] =='mm' : isTrig = e.HLT_IsoMu24 or e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
        if cat[2:] =='mt' : isTrig = e.HLT_IsoMu24 or e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
        if isTrig :
            cutCounter[cat].count('Lep_Trig2')
            if  MC :   cutCounterGenWeight[cat].countGenWeight('Lep_Trig2', e.genWeight)


    met_pt = float(e.MET_pt)
    met_phi = float(e.MET_phi)


    if doJME :  #default after JME systematics with Smear
        if era!='2017' :
	    try : 
		met_pt = float(e.MET_T1_pt)
		met_phi = float(e.MET_T1_phi)
	    except AttributeError : 
		met_pt = float(e.MET_pt)
		met_phi = float(e.MET_pt)

        if era=='2017' and args.era=='UL':
            try : 
		met_pt = float(e.MET_T1_pt)
		met_phi = float(e.MET_T1_phi)
	    except AttributeError : 
		met_pt = float(e.MET_pt)
		met_phi = float(e.MET_phi)

    ##print met_pt, 'smear', e.MET_T1Smear_pt, 'uncorrected?', e.MET_pt
    tauMass=[]
    tauPt=[]
    eleMass=[]
    elePt=[]
    muMass=[]
    muPt=[]
    metPtPhi=[]
    metPtPhi.append(float(met_pt))
    metPtPhi.append(float(met_phi))

    if MC : 
	if len(muMass) == 0 :
	    for j in range(e.nMuon):
		muMass.append(e.Muon_mass[j])
		muPt.append(e.Muon_pt[j])

	if len(eleMass) == 0 :
	    for j in range(e.nElectron):
		eleMass.append(e.Electron_mass[j])
		elePt.append(e.Electron_pt[j])

	if len(tauMass) == 0 :
	    for j in range(e.nTau):
		tauMass.append(e.Tau_mass[j])
		tauPt.append(e.Tau_pt[j])
    
    '''if 'Hpp' in  args.nickName :
       if 't' not in GF.printGenDecayMode(e): continue
       #if not GF.printGenDecayMode(e)=='eeee': continue
       GF.printMC(e)
       print 'Gen channel is', GF.printGenDecayMode(e)
    else:
       if not GF.printGenDecayModeBkg(e,bkg=args.nickName) == args.category : continue
       #GF.printMC(e)
       #print 'Gen channel is', GF.printGenDecayMode(e)
    '''
    for isyst, systematic in enumerate(sysT) : 
	if isyst>0 : #use the default pT/mass for Ele/Muon/Taus before doing any systematic
	#if 'Central' in systematic or 'prong' in systematic : #use the default pT/mass for Ele/Muon/Taus before doing the Central or the tau_scale systematics ; otherwise keep the correction

	    for j in range(e.nMuon): 
                e.Muon_pt[j] = muPt[j]
                e.Muon_mass[j] = muMass[j]
	    for j in range(e.nElectron): 
                e.Electron_pt[j] = elePt[j]
                e.Electron_mass[j] = eleMass[j]
	    for j in range(e.nTau): 
                e.Tau_pt[j] = tauPt[j]
                e.Tau_mass[j] = tauMass[j]

        
	'''if isMC: 

	    met_pt, met_phi, metlist, philist = Weights.applyES(e, args.year, systematic, metPtPhi, allMET)
	    
	    if systematic == 'Central' :
		for i, j in enumerate (metlist): 
		    outTuple.list_of_arrays[i][0] = metlist[i]
		for i, j in enumerate (philist): 
		    outTuple.list_of_arrays[i+len(metlist)][0] = philist[i]
        ''' 
	goodElectronList = TF.makeGoodElectronListDCH(e,cutflow['ele'])
	goodMuonList = TF.makeGoodMuonListDCH(e,cutflow['muon'])
	goodTauList = TF.makeGoodTauList(e,cutflow['tau'])
        dupl = 0 #to count duplicate channels
	selected_evts += 1
	#=============3-lep=====================================
	if len(goodElectronList)+len(goodMuonList)+len(goodTauList) == 3:
		evts_3lep += 1
		for cat3L in cats3L :
			if(e.nElectron < cat3L.count('e') or e.nMuon < cat3L.count('m') or e.nTau < cat3L.count('t')): continue
			if(len(goodElectronList) < cat3L.count('e') or len(goodMuonList) < cat3L.count('m') or len(goodTauList) < cat3L.count('t')): continue
			dch1 = cat3L[:2]
			dch2 = cat3L[2:]
			#if args.category != 'none' and not dch1 in args.category : continue
			#if args.category != 'none' and not dch2 in args.category : continue

			## pairListP -> containes the TLV of the ++ pair  
			## lepListP -> containes the indices of the ++ pair  
			#lepList, lepListP, lepListM = [], [], []
			signC=0
			netS=-99
			bestDCH1 = []#contains indices
			if dch1 == 'ee':
				bestDCH1 = TF.getBestEEPair(entry=e, dch=dch1, pairList=[])
				if len(bestDCH1) < 2: continue
				signC = e.Electron_charge[bestDCH1[0]]
				netS = e.Electron_charge[bestDCH1[0]] + e.Electron_charge[bestDCH1[1]]
				#print "DCH1 charge:",netS
			if dch1 == 'em':
				bestDCH1 = TF.getBestEMuTauPair(entry=e, dch=dch1, pairList=[])
				if len(bestDCH1) < 2: continue
				signC = e.Electron_charge[bestDCH1[0]]
				netS = e.Electron_charge[bestDCH1[0]] + e.Muon_charge[bestDCH1[1]]
				#print "DCH1 charge:",netS
			if dch1 == 'et':
				bestDCH1 = TF.getBestETauPair(entry=e, dch=dch1, pairList=[])
				if len(bestDCH1) < 2: continue
				signC = e.Electron_charge[bestDCH1[0]]
				netS = e.Electron_charge[bestDCH1[0]] + e.Tau_charge[bestDCH1[1]]
				#print "DCH1 charge:",netS
			if dch1 == 'mm':
				bestDCH1 = TF.getBestMuMuPair(entry=e, dch=dch1, pairList=[])
				if len(bestDCH1) < 2: continue 
				signC = e.Muon_charge[bestDCH1[0]]
				netS = e.Muon_charge[bestDCH1[0]] + e.Muon_charge[bestDCH1[1]]
				#print "DCH1 charge:",netS
			if dch1 == 'mt':
				bestDCH1 = TF.getBestMuTauPair(entry=e, dch=dch1, pairList=[])
				if len(bestDCH1) < 2: continue
				signC = e.Muon_charge[bestDCH1[0]]
				netS = e.Muon_charge[bestDCH1[0]] + e.Tau_charge[bestDCH1[1]]
				#print "DCH1 charge:",netS
			if dch1 == 'tt':
				tauList = TF.getTauList(dch1, entry=e, pairList=[])
				bestDCH1 = TF.getBestTauPair(dch1, entry=e, tauList=tauList)
				if len(bestDCH1) < 2: continue
				signC = e.Tau_charge[bestDCH1[0]]
				netS = e.Tau_charge[bestDCH1[0]] + e.Tau_charge[bestDCH1[1]]
				#print "DCH1 charge:",netS
			##print 'signC ',signC
			pairList1 = TF.make4Vec(bestDCH1,dch1,e)
                        list1 = goodElectronList+goodMuonList+goodTauList
			lep_3 = -99
			if dch2 =='e':
				for i in goodElectronList:
					if e.Electron_charge[i] == signC : continue 
					else: lep_3 = i 
                        elif dch2 =='m':
                                for i in goodMuonList:
                                        if e.Muon_charge[i] == signC : continue         
                                        else: lep_3 = i
                        elif dch2  =='t':
                                for i in goodTauList:
                                        if e.Tau_charge[i] == signC : continue
                                        else: lep_3 = i
                        if lep_3 == -99: continue

			SVFit = False
			if not MC : isMC = False 
			outTuple.Fill3L(e,SVFit,cat3L,bestDCH1,lep_3, isMC,era,doJME, met_pt, met_phi,  isyst, tauMass, tauPt, eleMass, elePt, muMass, muPt, args.era)
			#=========================================================	

        if len(goodElectronList)+len(goodMuonList)+len(goodTauList) > 4:
            evts_5lep += 1
            continue;#can be recovered for signal
        elif len(goodElectronList)+len(goodMuonList)+len(goodTauList) < 3:
            continue# remove to #print details of all the vetoed evts
        evt_charge = 0
        for i in goodElectronList:
	    evt_charge += e.Electron_charge[i]
        for i in goodMuonList:
            evt_charge += e.Muon_charge[i]
        for i in goodTauList:
            evt_charge += e.Tau_charge[i]
        #if evt_charge !=0: 
	#    continue
        #print 'Event: ', count, ' #e: ',len(goodElectronList), ' #mu: ',len(goodMuonList), ' #t: ', len(goodTauList), '--> Good candidates'
        #print 'Event: ', count, ' #e: ',e.nElectron, ' #mu: ',e.nMuon, ' #t: ', e.nTau, '--> All reco'
        #its faster to modify the pair functions to use these good lists rather than passing all the candidates through the cuts. 
	for cat in cats :
            if(e.nElectron < cat.count('e') or e.nMuon < cat.count('m') or e.nTau < cat.count('t')): continue
            if(len(goodElectronList) < cat.count('e') or len(goodMuonList) < cat.count('m') or len(goodTauList) < cat.count('t')): continue
            dch1 = cat[:2]
            dch2 = cat[2:]
	    #if args.category != 'none' and not dch1 in args.category : continue
            #if args.category != 'none' and not dch2 in args.category : continue

            '''if (dch1 == 'ee' or dch2 == 'ee') and len(goodElectronList) < 2 :
                if printOn:   #print cat, e.run, e.luminosityBlock,  e.event, 'failed as you have < 2 goodElectronList'
                continue
            if (dch1 == 'em' or dch2 == 'em') and (len(goodElectronList) < 1 or len(goodMuonList) < 1) :
                if printOn:   #print cat, e.run, e.luminosityBlock,  e.event, 'failed as you have < 1 goodElectron or < 1 goodMuon'
                continue
            if (dch1 == 'et' or dch2 == 'et') and (len(goodElectronList) < 1 or len(goodTauList) < 1) :
                if printOn:   #print cat, e.run, e.luminosityBlock,  e.event, 'failed as you have < 1 goodElectronList or < 1 goodTauList'
                continue
            if (dch1 == 'mm' or dch2 == 'mm') and len(goodMuonList) < 2 :
                if printOn:   #print cat, e.run, e.luminosityBlock,  e.event, 'failed as you have < 2 goodMuonList'
                continue
            if (dch1 == 'mt' or dch2 == 'mt') and (len(goodMuonList) < 1 or len(goodTauList) < 1) :
                if printOn:   #print cat, e.run, e.luminosityBlock,  e.event, 'failed as you have < 1 goodMuonList or < 1 goodTauList'
                continue
            if (dch1 == 'tt' or dch2 == 'tt') and  len(goodTauList) < 2 :
                if printOn:   #print cat, e.run, e.luminosityBlock,  e.event, 'failed as you have < 2 goodTauList'
                continue
            '''

            ## pairListP -> containes the TLV of the ++ pair  
            ## lepListP -> containes the indices of the ++ pair  
            #lepList, lepListP, lepListM = [], [], []
            signC=0
            netS=-99
            bestDCH1 = []#contains indices
            if dch1 == 'ee':
                bestDCH1 = TF.getBestEEPair(entry=e, dch=dch1, pairList=[])
                if len(bestDCH1) < 2: continue
                signC = e.Electron_charge[bestDCH1[0]]
                netS = e.Electron_charge[bestDCH1[0]] + e.Electron_charge[bestDCH1[1]]
                #print "DCH1 charge:",netS
            if dch1 == 'em':
                bestDCH1 = TF.getBestEMuTauPair(entry=e, dch=dch1, pairList=[])
                if len(bestDCH1) < 2: continue
                signC = e.Electron_charge[bestDCH1[0]]
                netS = e.Electron_charge[bestDCH1[0]] + e.Muon_charge[bestDCH1[1]]
                #print "DCH1 charge:",netS
            if dch1 == 'et':
                bestDCH1 = TF.getBestETauPair(entry=e, dch=dch1, pairList=[])
                if len(bestDCH1) < 2: continue
                signC = e.Electron_charge[bestDCH1[0]]
                netS = e.Electron_charge[bestDCH1[0]] + e.Tau_charge[bestDCH1[1]]
                #print "DCH1 charge:",netS
            if dch1 == 'mm':
                bestDCH1 = TF.getBestMuMuPair(entry=e, dch=dch1, pairList=[])
                if len(bestDCH1) < 2: continue 
                signC = e.Muon_charge[bestDCH1[0]]
                netS = e.Muon_charge[bestDCH1[0]] + e.Muon_charge[bestDCH1[1]]
                #print "DCH1 charge:",netS
            if dch1 == 'mt':
                bestDCH1 = TF.getBestMuTauPair(entry=e, dch=dch1, pairList=[])
                if len(bestDCH1) < 2: continue
                signC = e.Muon_charge[bestDCH1[0]]
                netS = e.Muon_charge[bestDCH1[0]] + e.Tau_charge[bestDCH1[1]]
                #print "DCH1 charge:",netS
            if dch1 == 'tt':
                tauList = TF.getTauList(dch1, entry=e, pairList=[])
                bestDCH1 = TF.getBestTauPair(dch1, entry=e, tauList=tauList)
                if len(bestDCH1) < 2: continue
                signC = e.Tau_charge[bestDCH1[0]]
                netS = e.Tau_charge[bestDCH1[0]] + e.Tau_charge[bestDCH1[1]]
                #print "DCH1 charge:",netS
            ##print 'signC ',signC

            pairList1 = TF.make4Vec(bestDCH1,dch1,e)
            bestDCH2 = []
            if dch2 == 'ee':
                bestDCH2 = TF.getBestEEPair(entry=e, dch=dch2, pairList=pairList1, isDCH2=True,signC=signC)
                if len(bestDCH2) < 2: continue
                signC = e.Electron_charge[bestDCH2[0]]
                netS = e.Electron_charge[bestDCH2[0]] + e.Electron_charge[bestDCH2[1]]
                #print "DCH2 charge:",netS
            if dch2 == 'em':
                bestDCH2 = TF.getBestEMuTauPair(entry=e, dch=dch2, pairList=pairList1, isDCH2=True,signC=signC)
                if len(bestDCH2) < 2: continue
                signC = e.Electron_charge[bestDCH2[0]]
                netS = e.Electron_charge[bestDCH2[0]] + e.Muon_charge[bestDCH2[1]]
                #print "DCH2 charge:",netS
            if dch2 == 'et':
                bestDCH2 = TF.getBestETauPair(entry=e, dch=dch2, pairList=pairList1, isDCH2=True,signC=signC)
                if len(bestDCH2) < 2: continue
                signC = e.Electron_charge[bestDCH2[0]]
                netS = e.Electron_charge[bestDCH2[0]] + e.Tau_charge[bestDCH2[1]]
                #print "DCH2 charge:",netS
            if dch2 == 'mm':
                bestDCH2 = TF.getBestMuMuPair(entry=e, dch=dch2, pairList=pairList1, isDCH2=True,signC=signC)
                if len(bestDCH2) < 2: continue
                signC = e.Muon_charge[bestDCH2[0]]
                netS = e.Muon_charge[bestDCH2[0]] + e.Muon_charge[bestDCH2[1]]
                #print "DCH2 charge:",netS
            if dch2 == 'mt':
                bestDCH2 = TF.getBestMuTauPair(entry=e, dch=dch2, pairList=pairList1, isDCH2=True,signC=signC)
                if len(bestDCH2) < 2: continue
                signC = e.Muon_charge[bestDCH2[0]]
                netS = e.Muon_charge[bestDCH2[0]] + e.Tau_charge[bestDCH2[1]]
                #print "DCH2 charge:",netS
            if dch2 == 'tt':
                tauList = TF.getTauList(dch2, entry=e, pairList=pairList1, isDCH2=True,signC=signC)
                bestDCH2 = TF.getBestTauPair(dch2, entry=e, tauList=tauList)
                if len(bestDCH2) < 2: continue
                signC = e.Tau_charge[bestDCH2[0]]
                netS = e.Tau_charge[bestDCH2[0]] + e.Tau_charge[bestDCH2[1]]
                #print "DCH2 charge:",netS

            pairList2 = TF.make4Vec(bestDCH2,dch2,e)
            '''if 1>0:# GF.printGenDecayModeBkg(e,bkg=args.nickName) == args.category :
                print count, cat, 'e:', len(goodElectronList), 'm:',len(goodMuonList), 't:',len(goodTauList)
                GF.printMC(e)
                for i in goodElectronList:
                    print 'Ele ',e.Electron_charge[i], i, GF.genMatch(e,i,'e'),e.Electron_pt[i],e.Electron_eta[i],e.Electron_phi[i]
                for i in goodMuonList:
                    print 'Mu ',e.Muon_charge[i], i, GF.genMatch(e,i,'m'),e.Muon_pt[i],e.Muon_eta[i],e.Muon_phi[i]
                for i in goodTauList:
                    print 'Tau ',e.Tau_charge[i], i, GF.genMatch(e,i,'t'),e.Tau_pt[i],e.Tau_eta[i],e.Tau_phi[i]
            '''
            cat_yield[cat] += 1
            n_lepton[cat] = np.array(n_lepton[cat]) + np.array([len(goodElectronList), len(goodMuonList), len(goodTauList)])

            cutCounter[cat].count(dch1+dch2)
            if  MC:   cutCounterGenWeight[cat].countGenWeight(dch1+dch2, e.genWeight)
            #continue
	    #GF.printMC(e)
            
            dupl+=1
	    if dupl>1:
               print 'AHA',count,'has a fake channel'
               continue
            pass_evts += 1
            #continue

	    '''if len(bestTauPair) < 1 : 
		if unique :
		    print("Tau Pair Fail: Event ID={0:d} cat={1:s}".format(e.event,cat))
		    bestTauPair = TF.getBestEMuTauPair(e,dch=dch1,pairList=LeptV,printOn=True)#change bestTauPair to bestDCH 
		    GF.printEvent(e)
		if False and maxPrint > 0 and (dch2 == GF.eventID(e)[2:4]) :
		    maxPrint -= 1
		    print("Failed tau-pair cut")
		    print("Event={0:d} cat={1:s}".format(e.event,cat))
		    print("goodMuonList={0:s} goodElectronList={1:s} Mll={3:.1f} bestTauPair={4:s}".format(
			str(goodMuonList),str(goodElectronList),str(pairList),M,str(bestTauPair)))
		    print("Lep1.pt() = {0:.1f} Lep2.pt={1:.1f}".format(pairList[0].Pt(),pairList[1].Pt()))
		    GF.printEvent(e)
		    GF.printMC(e)
		#continue
	    if len(bestTauPair) > 1 :
		jt1, jt2 = bestTauPair[0], bestTauPair[1]
	    else :
		continue
            print("ddddddddddddddddddddddddddddddddd")#ggr
	    cutCounter[cat].count("GoodTauPair")
	    if  MC:   cutCounterGenWeight[cat].countGenWeight('GoodTauPair', e.genWeight)
            '''

	    if MC :
		outTuple.setWeight(PU.getWeight(e.PV_npvs)) 
		outTuple.setWeightPU(PU.getWeight(e.Pileup_nPU)) 
		outTuple.setWeightPUtrue(PU.getWeight(e.Pileup_nTrueInt)) 
		##print 'nPU', e.Pileup_nPU, e.Pileup_nTrueInt, PU.getWeight(e.Pileup_nPU), PU.getWeight(e.Pileup_nTrueInt), PU.getWeight(e.PV_npvs), PU.getWeight(e.PV_npvsGood)
	    else : 
		outTuple.setWeight(1.) 
		outTuple.setWeightPU(1.) ##
		outTuple.setWeightPUtrue(1.)
            #print("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")#ggr

	    SVFit = False
            #continue	    
	    if not MC : isMC = False

            #outTuple.Fill(e,SVFit,cat,bestDCH2[0],bestDCH2[1],pairList1[0],pairList1[1],bestDCH1,isMC,era,doJME, met_pt, met_phi,  isyst, tauMass, tauPt, eleMass, elePt, muMass, muPt, args.era)
	    outTuple.Fill(e,SVFit,cat,bestDCH1,bestDCH2,isMC,era,doJME, met_pt, met_phi,  isyst, tauMass, tauPt, eleMass, elePt, muMass, muPt, args.era)
            '''
	    if maxPrint > 0 :
		maxPrint -= 1
		print("\n\nGood Event={0:d} cat={1:s}  MCcat={2:s}".format(e.event,cat,GF.eventID(e)))
		print("goodMuonList={0:s} goodElectronList={1:s} Mll={2:.1f} bestTauPair={3:s}".format(
		    str(goodMuonList),str(goodElectronList),M,str(bestTauPair)))
		print("Lep1.pt() = {0:.1f} Lep2.pt={1:.1f}".format(pairList[0].Pt(),pairList[1].Pt()))
		GF.printEvent(e)
		print("Event ID={0:s} cat={1:s}".format(GF.eventID(e),cat))
            '''
dT = time.time() - tStart
print("Run time={0:.2f} s  time/event={1:.1f} us".format(dT,1000000.*dT/count))

hLabels=[]
hLabels.append('All')
hLabels.append('inJSON')
hLabels.append('METfilter')
hLabels.append('Lep_Trig1')
hLabels.append('Lep_Trig2')
#hLabels.append('TwoLeptons')
#hLabels.append('GoodTauPair')
hLabels.append('channel')

hCutFlow=[]
hCutFlowW=[]

outTuple.writeTree()
fW = TFile( outFileName, 'update' )
fW.cd()
hNEvts = TH1D("hNEvts", "nEntries", 1, 0,1)
hNEvts.Fill(0,nentries) 
hNEvts.Write()
#print '------------------------->',fW, outFileName
for icat,cat in enumerate(cats) :
    print('\nSummary for {0:s}'.format(cat))
    cutCounter[cat].printSummary()
    hName="hCutFlow_"+str(cat)
    hNameW="hCutFlowWeighted_"+str(cat)
    lcount=len(hLabels)
    hCutFlow.append( TH1D(hName,hName,lcount,0.5,lcount+0.5))
    if MC  : hCutFlowW.append( TH1D(hNameW,hNameW,lcount,0.5,lcount+0.5))
    #print lcount, cat, icat
    for i in range(len(hLabels)) :
        hCutFlow[icat].GetXaxis().SetBinLabel(i+1,hLabels[i])
        if MC : hCutFlowW[icat].GetXaxis().SetBinLabel(i+1,hLabels[i])

    '''for i in range(lcount) :
        #hCutFlow[cat].Fill(1, float(cutCounter[cat].getYield()[i]))
        yields = cutCounter[cat].getYield()[i]
        hCutFlow[icat].Fill(i+1, float(yields))
    
        if MC : 
	    yieldsW = cutCounterGenWeight[cat].getYieldWeighted()[i]
            hCutFlowW[icat].Fill(i+1, float(yieldsW))
        ##print cutCounter[cat].getYield()[i], i, cutCounter[cat].getLabels()[i]
    '''
       
    hCutFlow[icat].Sumw2()
    hCutFlow[icat].Write()
    if MC : 
        hCutFlowW[icat].Sumw2()
        hCutFlowW[icat].Write()
    icat+=1
if not MC : CJ.printJSONsummary()

print '# Yields in each channel ', cat_yield, 'Total entries ',nentries
#print 'n_leptons in each channel [e, mu, tau]',n_lepton
'''cutflow_ele   = TH1D( 'cutflow_ele', 'Good Ele cutflow', 15, 0., 15 )
cutflow_mu   = TH1D( 'cutflow_mu', 'Good Muon cutflow', 15, 0., 15 )
cutflow_tau   = TH1D( 'cutflow_tau', 'Good Tau cutflow', 15, 0., 15 )
for icut in range(0,8): 
   cutflow_ele.Fill(icut,cutflow['ele'].getYield()[icut])
   #cutflow_mu.Fill(icut,cutflow['muon'].getYield()[icut])
   cutflow_tau.Fill(icut,cutflow['tau'].getYield()[icut])
cutflow_ele.Write()
cutflow_mu.Write()
cutflow_tau.Write()
'''
print '# of selected events', selected_evts,'\n# of 3 lep events', evts_3lep, '\n# of 5 lep evts',evts_5lep,'\n# of passed events',pass_evts 
#file.close()  # Close the file

