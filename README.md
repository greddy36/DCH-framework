## Set up working area
```
cmsrel CMSSW_13_0_10
cd CMSSW_13_0_10/src
cmsenv
git clone https://github.com/greddy36/DCH-framework.git
```
## Run the code locally
Once yyou've set up the working area, you can run the code as below.
eg:
```
cd DCH/
python DCH_noPairing.py -f  DY.root -o DY_out.root --nickName DY -y 2018 -s DCH -w 0 -j 0 -e UL
```
## For submitting condor jobs
```
cd MC/
```
Follow the instructions in the README in MC directory

## cloned SVfit from
# incomplete!!!
```
git clone https://github.com/SVfit/ClassicSVfit TauAnalysis/ClassicSVfit
export LIBRARY_PATH=$LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
make -f TauAnalysis/ClassicSVfit/Makefile -j4
```
## cloned TauIDSFs from 
```

```
