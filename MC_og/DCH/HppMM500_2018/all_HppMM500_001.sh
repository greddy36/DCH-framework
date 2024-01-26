#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc820
eval `scramv1 project CMSSW CMSSW_10_6_5`
cd CMSSW_10_6_5/src
eval `scramv1 runtime -sh`
cmsenv
cd ${_CONDOR_SCRATCH_DIR}/CMSSW_10_6_5/src/
scram b -j 4
echo ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}/CMSSW_10_6_5/src/
cp ${_CONDOR_SCRATCH_DIR}/* .
tar -zxvf correctionlib.tar.gz
ls -altrh
echo 'this is the working dir' ${_CONDOR_SCRATCH_DIR}
cp cuts_DCH_2018.yaml cuts.yaml
xrdcp  root://cmseos.fnal.gov//store/mc/RunIISummer20UL18NanoAODv9/HPlusPlusHMinusMinusHTo2L_M-500_TuneCP5_13TeV_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2540000/8803D7C6-8501-AF49-A010-62BCD211144B.root inFile.root
if [ ! -f inFile.root ] ; 
 then 
 xrdcp root://cmsxrootdfnal.gov//store/mc/RunIISummer20UL18NanoAODv9/HPlusPlusHMinusMinusHTo2L_M-500_TuneCP5_13TeV_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2540000/8803D7C6-8501-AF49-A010-62BCD211144B.root inFile.root
 fi 
if [ ! -f inFile.root ] ; 
 then 
 xrdcp root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/HPlusPlusHMinusMinusHTo2L_M-500_TuneCP5_13TeV_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2540000/8803D7C6-8501-AF49-A010-62BCD211144B.root inFile.root
 fi 
python DCH.py -f inFile.root -o HppMM500_001.root --nickName HppMM500 -y 2018 -s DCH -w 0 -j no
rm inFile*.root
xrdcp HppMM500_001.ntup root://cmseos.fnal.gov//eos/uscms/store/user/greddy/DCH_files/DCH_out/HppMM500_2018/HppMM500_001.root

rm  *root *.ntup  *.so
rm *.pcm
rm *cc.d
cd ../.. ; rm -fr CMSSW_10_6_5
