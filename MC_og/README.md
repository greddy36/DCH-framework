This is how to create condor jobs to run the analyzer

# Step 1 - Make the csv listing the datasets you need - for example MCsamples_2018UL.csv. Please note the format

`WJetsToLNu_NLO, Other, 67350.7, 1, 1,,/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM`

in principle, the 1rst and the last arguments are used for the scripts, but the others can be used to define other parameters if needed be, like xsec, SFs etc

# Step 2- execute

`python makeMC.py -f MCsamples_2018UL.csv -y 2018 -s DCH -l no`

the `-s` flag : this is the executable that must be present in the src/DCH_KSU/DCH dir ie like DCH.py
the `-l` flag : if no, then it gets the dataset from DBS - if 'yes', then it can presume that you already have some ntuples stored locally etc


# Step 3 - take a look on the produced .sh from the previous step - it should have lines like
`mkdir -p DCH/HppM1000_2018 
cd DCH/HppM1000_2018
python ../../makeCondor.py --dataSet /HPlusPlusHMinusMinusHTo2L_M-1000_TuneCP5_13TeV_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM --nickName HppM1000 --mode anaXRD --year 2018 -c 5 -s DCH -j no -l no`

once executed, it will create inside the `./MC/DCH/HppM1000_2018` dir all the .sh and .jdl necessary for the work

inside the makeCondor.py there is the `-c `switch that groups input .root files into one job - default is 3 - **ATTENTION** - if you edit that and you re-execute the  runMC_2018_DCH.csh then you will end up in a situation where you have several .sh with overlapping input/output files. If you edit the concatenate, FIRST make sure to delete the corresponding MC/MYDATASET dir

a typical .sh set the CMSSW area and then executes the code and finally copies the .root to the final destinations
<pre>
cp cuts_DCH_2018.yaml cuts.yaml 
xrdcp  root://cmseos.fnal.gov//store/mc/RunIISummer20UL18NanoAODv9/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-
v2/80000/0B37487E-FC2E-D64C-89C1-F15ACD3F1904.root inFile.root 
python DCH.py -f inFile.root -o WJetsToLNu_NLO_001.root --nickName WJetsToLNu_NLO -y 2018 -s DCH -w 1 -j no 
rm inFile*.root 
xrdcp WJetsToLNu_NLO_001.ntup root://cmseos.fnal.gov//store/group/lpcsusyhiggs/ntuples/nAODv9/DCH_out/WJetsToLNu_NLO_2018/WJetsToLNu_NLO_001.root  
</pre>pre>



Please Note, that each file is copied, the DCH.py is executed (`-w 1` flag also saves the genWeights in a TH1F in the final file) and the .ntup is copied on a /eos area. This is an example, please edit the destination to your liking


# Step 4 - just cd inside the folder and submit .jdl 
