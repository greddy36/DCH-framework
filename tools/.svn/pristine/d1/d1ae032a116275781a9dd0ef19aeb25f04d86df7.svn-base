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

#jmeCorrections = createJMECorrector(True, 2016, "B", "Total", "AK4PFchs", False)


#def createJMECorrector(isMC=True, dataYear=2016, runPeriod="B", jesUncert="Total", jetType = "AK4PFchs", noGroom=False, metBranchName="MET", applySmearing=True, isFastSim=False, applyHEMfix=False, splitJER=False):


metbranch = "MET"
if sys.argv[2]=='2017' : 
    metbranch="METFixEE2017"
    print 'this is 2017 MET v2', metbranch

if sys.argv[1] =='False' :  jmeCorrections = createJMECorrector(isMC=False, dataYear=sys.argv[2], runPeriod=str(sys.argv[3]),  jesUncert="Merged", jetType = "AK4PFchs",  applyHEMfix=True, splitJER=False, metBranchName = metbranch)
if sys.argv[1] =='True' :  jmeCorrections = createJMECorrector(isMC=True, dataYear=sys.argv[2], runPeriod=str(sys.argv[3]),  jesUncert="Merged", jetType = "AK4PFchs",  applyHEMfix=True, splitJER=False, metBranchName = metbranch)


#jmeCorrections = createJMECorrector(isMC=True, dataYear=2016, runPeriod="B", jesUncert="Total", jetType = "AK4PFchs",  applyHEMfix=True, splitJER=True, metBranchName = metbranch)

#jmeCorrections = createJMECorrector(isMC=True, dataYear=2018, runPeriod="B", jesUncert="Merged", jetType = "AK4PFchs",  applyHEMfix=True, splitJER=False, metBranchName = metbranch)

fnames=["./inFile.root"]
#p=PostProcessor(".",fnames,"","",[jmeCorrections()],provenance=True)
p=PostProcessor("./",fnames,"Jet_pt>MIN && Jet_pt<MAX","",[jmeCorrections()],provenance=True)
p.run()
