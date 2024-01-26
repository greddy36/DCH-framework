mkdir -p DCH/HppMM500_2018
cd DCH/HppMM500_2018
python ../../makeCondor.py --dataSet /HPlusPlusHMinusMinusHTo2L_M-500_TuneCP5_13TeV_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM --nickName HppMM500 --mode anaXRD --year 2018 -c 5 -s DCH -j no -l no -w 0
cd /eos/uscms/store/user/greddy/KSU/CMSSW_10_6_5/src/MC
