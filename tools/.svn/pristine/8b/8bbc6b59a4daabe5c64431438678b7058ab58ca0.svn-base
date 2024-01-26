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
#from my_module import MyModule  # Import your custom module


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

json_file='Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'
if '2017' in sys.argv[2] :json_file='Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
if '2018' in sys.argv[2] :json_file='Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'

is_data = "False" in sys.argv[1] 

# If the input file is a data file, apply the JSON file
json_input = json_file if is_data else None

#json_input = None
p=''

#if '2016' in sys.argv[2]  : p=PostProcessor("./",fnames, "( (  (HLT_IsoMu24 || HLT_IsoTkMu24) && Muon_pt > 26. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.25 &&  Muon_tightId>0 && Muon_pfRelIso04_all< 0.5  && (Muon_isGlobal || Muon_isTracker) && Muon_pfIsoId>2  ) ) &&  Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")

#else :p=PostProcessor("./",fnames, "(   ( (HLT_IsoMu24 || HLT_IsoMu27) && Muon_pt > 29. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.25 && Muon_tightId>0  && Muon_pfRelIso04_all< 0.5  && (Muon_isGlobal || Muon_isTracker) && Muon_pfIsoId>2  ) ) && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")

#if '2016' in sys.argv[2]  : p=PostProcessor("./",fnames, "( (  (HLT_IsoMu24 || HLT_IsoTkMu24) && Muon_pt > 25. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.22 && Muon_pfRelIso04_all< 0.5  && (Muon_isGlobal || Muon_isTracker)   && Muon_pfIsoId>1 ) ) &&  Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt",jsonInput=json_input)

#else :p=PostProcessor("./",fnames, "(   ( (HLT_IsoMu27) && Muon_pt > 28. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.22 && Muon_pfRelIso04_all< 0.5  && (Muon_isGlobal || Muon_isTracker) && Muon_pfIsoId>1  ) ) && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt",jsonInput=json_input)
#else :p=PostProcessor("./",fnames, " Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",modules=[jmeCorrections(), MyModule()],provenance=True,  outputbranchsel="keep_and_drop.txt",jsonInput=json_input)
#AZH
if '2016' in sys.argv[2]  : p=PostProcessor("./",fnames, "( (  (HLT_IsoMu24 || HLT_IsoTkMu24) && Muon_pt > 26. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.25 && Muon_pfRelIso04_all< 0.5  && (Muon_isGlobal || Muon_isTracker) && Muon_pfIsoId>2  ) ) &&  Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")

else :p=PostProcessor("./",fnames, "(   ( (HLT_IsoMu24 || HLT_IsoMu27) && Muon_pt > 29. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.25 && Muon_pfRelIso04_all< 0.5  && (Muon_isGlobal || Muon_isTracker) && Muon_pfIsoId>2  ) ) && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt")


p.run()
