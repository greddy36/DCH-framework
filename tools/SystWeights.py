#!/usr/bin/env python

""" ZH.py: makes an nTuple for the ZH->tautau analysis """

__author__ = "Dan Marlow, Alexis Kalogeropoulos" 

# import external modules 
import sys
import numpy as np
from ROOT import TFile, TTree, TH1, TH1D, TCanvas, TLorentzVector  
from math import sqrt, pi
#from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool
#from TauPOG.TauIDSFs.TauIDSFTool import TauESTool
#from TauPOG.TauIDSFs.TauIDSFTool import TauFESTool

# import from ZH_Run2/funcs/
sys.path.insert(1,'../funcs/')
#import tauFun2 as TF
import generalFunctions as GF 
#import outTuple
import time

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-v","--verbose",default=0,type=int,help="Print level.")
    parser.add_argument("-f","--inFileName",default='ZHtoTauTau_test.root',help="File to be analyzed.")
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
    parser.add_argument("-e","--era",type=str, default='EOY',help="EOY of UL")
    
    return parser.parse_args()

args = getArgs()
print("args={0:s}".format(str(args)))
maxPrint = args.maxPrint 

cutCounter = {}
cutCounterGenWeight = {}

doJME  = args.doSystematics.lower() == 'true' or args.doSystematics.lower() == 'yes' or args.doSystematics == '1'

#doJME = True

cats = ['eeet','eemt','eett','eeem','mmet','mmmt','mmtt','mmem']

for cat in cats : 
    cutCounter[cat] = GF.cutCounter()
    cutCounterGenWeight[cat] = GF.cutCounter()

inFileName = args.inFileName
print("Opening {0:s} as input.  Event category {1:s}".format(inFileName,cat))

isAZH=False
if str(args.selection) == 'AZH' : isAZH = True
if isAZH : print 'You are running on the AZH mode !!!'

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
'''
if MC :
    print "this is MC, will get PU etc", args.dataType
    PU = GF.pileUpWeight()
    PU.calculateWeights(args.nickName,args.year)
    #else : PU.calculateWeights('TTJets_DiLept',args.year)
else :
    CJ = ''#GF.checkJSON(filein='Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt')
    if args.year == 2016 : CJ = GF.checkJSON(filein='Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt')
    if args.year == 2017 : CJ = GF.checkJSON(filein='Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt')
    if args.year == 2018 : CJ = GF.checkJSON(filein='Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt')
    if args.year == 2018 : CJ = GF.checkJSON(filein='Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt')

'''
print 'systematics', doJME

era=str(args.year)

outFileName = GF.getOutFileName(args).replace(".root",".ntup")

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

isZH = 'ZH' in outFileName or 'HZJ' in outFileName

if args.weights > 0 :
    hWeight = TH1D("hWeights","hWeights",1,-0.5,0.5)
    hWeightAll = TH1D("hWeightsAll","hWeightsAll",100,0,50)
    hWeight.Sumw2()
    hWeightAll.Sumw2()
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
        hWeightAll.Fill(e.genWeight)

        if isZH : 
	    hWeightScaleUp.Fill(0, e.genWeight*e.LHEScaleWeight[8])
	    hWeightScaleDown.Fill(0, e.genWeight*e.LHEScaleWeight[0])

	    hWeightScaleSTXS.Fill(0, e.genWeight)
	    hWeightScaleSTXSUp.Fill(0, e.genWeight*e.LHEScaleWeight[8])
	    hWeightScaleSTXSDown.Fill(0, e.genWeight*e.LHEScaleWeight[0])

	    if e.HTXS_stage1_1_cat_pTjet30GeV == 300 or e.HTXS_stage1_1_cat_pTjet30GeV == 400 or e.HTXS_stage1_1_cat_pTjet30GeV == 500 : 
	       hWeightScaleSTXS.Fill(1, e.genWeight)
	       hWeightScaleSTXSUp.Fill(1, e.genWeight*e.LHEScaleWeight[8])
	       hWeightScaleSTXSDown.Fill(1, e.genWeight*e.LHEScaleWeight[0])

	    if e.HTXS_stage1_1_cat_pTjet30GeV == 301 or e.HTXS_stage1_1_cat_pTjet30GeV == 401 or e.HTXS_stage1_1_cat_pTjet30GeV == 501 : 
	       hWeightScaleSTXS.Fill(2, e.genWeight)
	       hWeightScaleSTXSUp.Fill(2, e.genWeight*e.LHEScaleWeight[8])
	       hWeightScaleSTXSDown.Fill(2, e.genWeight*e.LHEScaleWeight[0])

	    if e.HTXS_stage1_1_cat_pTjet30GeV == 302 or e.HTXS_stage1_1_cat_pTjet30GeV == 402 or e.HTXS_stage1_1_cat_pTjet30GeV == 502 : 
	       hWeightScaleSTXS.Fill(3, e.genWeight)
	       hWeightScaleSTXSUp.Fill(3, e.genWeight*e.LHEScaleWeight[8])
	       hWeightScaleSTXSDown.Fill(3, e.genWeight*e.LHEScaleWeight[0])

	    if e.HTXS_stage1_1_cat_pTjet30GeV == 303 or e.HTXS_stage1_1_cat_pTjet30GeV == 403 or e.HTXS_stage1_1_cat_pTjet30GeV == 503 : 
	       hWeightScaleSTXS.Fill(4, e.genWeight)
	       hWeightScaleSTXSUp.Fill(4, e.genWeight*e.LHEScaleWeight[8])
	       hWeightScaleSTXSDown.Fill(4, e.genWeight*e.LHEScaleWeight[0])

	    if e.HTXS_stage1_1_cat_pTjet30GeV == 304 or e.HTXS_stage1_1_cat_pTjet30GeV == 404 or e.HTXS_stage1_1_cat_pTjet30GeV == 504 : 
	       hWeightScaleSTXS.Fill(5, e.genWeight)
	       hWeightScaleSTXSUp.Fill(5, e.genWeight*e.LHEScaleWeight[8])
	       hWeightScaleSTXSDown.Fill(5, e.genWeight*e.LHEScaleWeight[0])


	    if e.HTXS_stage1_1_cat_pTjet30GeV == 305 or e.HTXS_stage1_1_cat_pTjet30GeV == 405 or e.HTXS_stage1_1_cat_pTjet30GeV == 505 : 
	       hWeightScaleSTXS.Fill(6, e.genWeight)
	       hWeightScaleSTXSUp.Fill(6, e.genWeight*e.LHEScaleWeight[8])
	       hWeightScaleSTXSDown.Fill(6, e.genWeight*e.LHEScaleWeight[0])


    

        if "WJetsToLNu" in outFileName and 'TWJets' not in outFileName:

            npartons = ord(e.LHE_Njets)
	    if  npartons <= 4: 	hWxGenweightsArr[npartons].Fill(0, e.genWeight)
        if "DYJetsToLL" in outFileName :
            npartons = ord(e.LHE_Njets)
	    if  npartons <= 4 : hDYxGenweightsArr[npartons].Fill(0, e.genWeight)

    fName = GF.getOutFileName(args).replace(".root",".weights")
    fW = TFile( fName, 'recreate' )
    print 'Will be saving the Weights in', fName
    fW.cd()

    if "WJetsToLNu" in outFileName and 'TWJets' not in outFileName:
        for i in range(len(hWxGenweightsArr)):
            hWxGenweightsArr[i].Write()
    elif "DYJetsToLL" in outFileName:
        for i in range(len(hDYxGenweightsArr)):
            hDYxGenweightsArr[i].Write()

    hWeight.Write()
    hWeightAll.Write()
    hWeightScaleUp.Write()
    hWeightScaleDown.Write()
    hWeightScaleSTXS.Write()
    hWeightScaleSTXSUp.Write()
    hWeightScaleSTXSDown.Write()
    fW.Close()
    sys.exit()


