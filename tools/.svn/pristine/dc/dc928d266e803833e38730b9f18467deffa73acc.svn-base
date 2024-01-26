#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
##soon to be deprecated
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
##new way of using jme uncertainty
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *



#from  exampleModule import *

##Function parameters
##(isMC=True, dataYear=2016, runPeriod="B", jesUncert="Total", redojec=False, jetType = "AK4PFchs", noGroom=False)
##All other parameters will be set in the helper module
#sys.args are  sys.argv[1]= isdata, , sys.argv[2]= era, sys.argv[3]=period

#if sys.argv[1] =='True' :  jmeCorrections = createJMECorrector(isMC=True, dataYear=sys.argv[2], str(sys.argv[3]), "Total", "AK4PFchs", jesUncert="All", redojec=True)

#jetmetCorrector = createJMECorrector(isMC=True, dataYear=2016, jesUncert="All", redoJec=True)



metbranch = "MET"
if sys.argv[2]=='2017' : 
    metbranch="METFixEE2017"
    print 'this is 2017 MET v2', metbranch

if sys.argv[1] =='False' :  jmeCorrections = createJMECorrector(isMC=False, dataYear=str(sys.argv[2]), runPeriod=str(sys.argv[3]),  jesUncert="Total", jetType = "AK4PFchs",  applyHEMfix=True, splitJER=False, metBranchName = metbranch, applySmearing = False)
if sys.argv[1] =='True' :  jmeCorrections = createJMECorrector(isMC=True, dataYear=str(sys.argv[2]), runPeriod=str(sys.argv[3]),  jesUncert="Total", jetType = "AK4PFchs",  applyHEMfix=True, splitJER=False, metBranchName = metbranch, applySmearing = False)



#jmeCorrections = createJMECorrector(isMC=True, dataYear=2016, runPeriod="B", jesUncert="Total", jetType = "AK4PFchs",  applyHEMfix=True, splitJER=True, metBranchName = metbranch)
#jmeCorrections = createJMECorrector(isMC=True, dataYear=2018, runPeriod="B", jesUncert="Merged", jetType = "AK4PFchs",  applyHEMfix=True, splitJER=False, metBranchName = metbranch)

fnames=["./inFile.root"]
#p=PostProcessor(".",fnames,"","",[jmeCorrections()],provenance=True)
#p=PostProcessor("./",fnames,"Jet_pt>10 && abs(Jet_eta)<4.7001 &&    Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[jmeCorrections()],provenance=True)

p=''

# HLT_IsoMu22 and not HLT_IsoMu22_eta2p1 and not HLT_IsoTkMu22 and not HLT_IsoTkMu22_eta2p1
#if '2016' in sys.argv[2]  : p=PostProcessor("./",fnames, "(HLT_IsoMu22 || HLT_IsoMu22_eta2p1 || HLT_IsoTkMu22 || HLT_IsoTkMu22_eta2p1) && Muon_pt > 20. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.25 && Muon_looseId &&  Muon_pfRelIso04_all < 0.3 && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")
#else :p=PostProcessor("./",fnames, "(HLT_IsoMu24 || HLT_IsoMu27) && Muon_pt > 20. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.25 && Muon_looseId &&  Muon_pfRelIso04_all < 0.3 && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")
#if '2016' in sys.argv[2]  : p=PostProcessor("./",fnames, "(HLT_IsoMu22 || HLT_IsoMu22_eta2p1 || HLT_IsoTkMu22 || HLT_IsoTkMu22_eta2p1) && Muon_pt > 24. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.25 && Muon_tightId>0 && Muon_pfRelIso04_all< 0.5 && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")
#IsoTkMu24

#if '2016' in sys.argv[2]  : p=PostProcessor("./",fnames, "(((HLT_IsoMu24 || HLT_IsoTkMu24) && Muon_pt > 26. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.25 && Muon_tightId>0 && Muon_pfRelIso04_all< 0.5)  || ((HLT_Ele25_eta2p1_WPTight_Gsf) && Electron_pt > 27. && abs(Electron_eta) < 2.2 && abs(Electron_dxy) < 0.05 && abs    (Electron_dz)<0.25 && Electron_mvaFall17V2Iso_WP90>0 && Electron_pfRelIso03_all< 0.5 )) && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")

###not used anymore
#if '2016' in sys.argv[2]  : p=PostProcessor("./",fnames, "( (  (HLT_IsoMu24 || HLT_IsoTkMu24) && Muon_pt > 26. && abs(Muon_eta) < 2.5 &&  Muon_tightId>0 && Muon_pfRelIso04_all< 0.5  ) || (  HLT_Ele25_eta2p1_WPTight_Gsf && Electron_pt > 27. && abs(Electron_eta) < 2.2 && Electron_mvaFall17V2Iso_WP90>0 && Electron_pfRelIso03_all< 0.5 )) && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")
#else :p=PostProcessor("./",fnames, "( (HLT_Ele35_WPTight_Gsf && Electron_pt > 37. && abs(Electron_eta) < 2.2 && Electron_mvaFall17V2Iso_WP90>0 && Electron_pfRelIso03_all< 0.5)  ||   (HLT_IsoMu27 && Muon_pt > 29. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.25 && Muon_tightId>0  && Muon_pfRelIso04_all< 0.5) ) && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")
#else :p=PostProcessor("./",fnames, "(HLT_IsoMu27) && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")


#if is MC
if sys.argv[1] =='True':

    if '2016' in sys.argv[2]  : p=PostProcessor("./",fnames, "( (  (HLT_IsoMu24 || HLT_IsoTkMu24) && Muon_pt > 26. && abs(Muon_eta) < 2.5 &&  Muon_tightId>0 && Muon_pfRelIso04_all< 0.5  ) || (  HLT_Ele25_eta2p1_WPTight_Gsf && Electron_pt > 27. && abs(Electron_eta) < 2.2 && Electron_mvaFall17V2Iso_WP90>0 && Electron_pfRelIso03_all< 0.5 ))    && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")

    #else :p=PostProcessor("./",fnames, "( (HLT_Ele35_WPTight_Gsf && Electron_pt > 37. && abs(Electron_eta) < 2.6 && Electron_mvaFall17V2Iso_WP90>0 && Electron_pfRelIso03_all< 0.5)  ||   (  (HLT_IsoMu24 || HLT_IsoMu27) && Muon_pt > 29. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && Muon_tightId>0  && Muon_pfRelIso04_all< 0.5) ) && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")
    else :p=PostProcessor("./",fnames, "(   ( (HLT_IsoMu24 || HLT_IsoMu27) && Muon_pt > 29. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.25 && Muon_tightId>0  && Muon_pfRelIso04_all< 0.5  && (Muon_isGlobal || Muon_isTracker) && Muon_pfIsoId>2  ) ) && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")
    #else :p=PostProcessor("./",fnames, " (  (HLT_IsoMu24 || HLT_IsoMu27  || HLT_Ele35_WPTight_Gsf)  ) && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")


if sys.argv[1] =='False':

    if '2016' in sys.argv[2]  : p=PostProcessor("./",fnames, "( (  (HLT_IsoMu24 || HLT_IsoTkMu24) && Muon_pt > 26. && abs(Muon_eta) < 2.5 &&  Muon_tightId>0 && Muon_pfRelIso04_all< 0.5  ) &&  Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")

    #else :p=PostProcessor("./",fnames, "(   (HLT_IsoMu24 || HLT_IsoMu27 )  ) && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")
    else :p=PostProcessor("./",fnames, "(   ( (HLT_IsoMu24 || HLT_IsoMu27) && Muon_pt > 29. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.25 && Muon_tightId>0  && Muon_pfRelIso04_all< 0.5  && (Muon_isGlobal || Muon_isTracker) && Muon_pfIsoId>2  ) ) && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")
    #else :p=PostProcessor("./",fnames, "( (  (HLT_IsoMu24 || HLT_IsoMu27) && Muon_pt > 29. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && Muon_tightId>0  && Muon_pfRelIso04_all< 0.5) ) && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")

#else :p=PostProcessor("./",fnames, "(HLT_IsoMu24 || HLT_IsoMu27) && Muon_pt > 26. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.25 && Muon_tightId>0  && Muon_pfRelIso04_all< 0.5 && genWeight<15    && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")

p.run()
