# ClassicSVfit4tau
Customized version of SVfit for hh->4tau analysis

# Installation instructions
The SVfitPerformanceStudies package has been tested with CMSSW 9_4_4.
It depends on the following other packages:
- TauAnalysis/ClassicSVfit
- TauAnalysis/SVfitTF
- VAMP (library for numeric integration, published in Comput.Phys.Commun. 120 (1999) 13)

In order to install the code, execute:

```
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src
cmsenv
echo "installing SVfit..."
git clone https://github.com/SVfit/ClassicSVfit4tau TauAnalysis/ClassicSVfit4tau
git clone https://github.com/SVfit/ClassicSVfit TauAnalysis/ClassicSVfit
git clone https://github.com/SVfit/SVfitTF TauAnalysis/SVfitTF

echo "installing VAMP library..."
mkdir $CMSSW_BASE/VAMP
wget http://whizard.hepforge.org/oldsrc/vamp-2.3.0.tar.gz -P $CMSSW_BASE/VAMP
tar zxvf $CMSSW_BASE/VAMP/vamp-2.3.0.tar.gz -C $CMSSW_BASE/VAMP
mkdir -p $CMSSW_BASE/VAMP/vamp-2.3.0/prefix/share/doc/vamp
cd $CMSSW_BASE/VAMP/vamp-2.3.0
./configure --prefix=$CMSSW_BASE/VAMP/vamp-2.3.0/prefix
make -j4
make install
cp $CMSSW_BASE/VAMP/vamp-2.3.0/prefix/lib/* $CMSSW_BASE/lib/$SCRAM_ARCH
cp $CMSSW_BASE/src/TauAnalysis/ClassicSVfit4tau/vamp.xml $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/vamp.xml
scram setup vamp
scram tool info vamp

echo "compiling SVfit..."
cd $CMSSW_BASE/src
scram b -j 4
```

In case of compilation problems, please sutmitt an issue on
https://github.com/SVfit/ClassicSVfit4tau/issues

# Without CMSSW software

It is possible to build the software without CMSSW framework if the following prerequisites are satisfied (oldest software version the instructions were tested with):
- ROOT (6.10/3 or newer)
- GCC (6.3 or newer)
- ClassicSVfit (instructions available [here](https://github.com/SVfit/ClassicSVfit#without-cmssw-software))
- VAMP 2.3.0

In order to build VAMP, make sure that you're in the same directory where you installed `ClassicSVfit`.
Then execute:
```bash
mkdir VAMP
cd $_
wget http://whizard.hepforge.org/oldsrc/vamp-2.3.0.tar.gz
tar zxvf vamp-2.3.0.tar.gz
rm vamp-2.3.0.tar.gz
cd vamp-2.3.0
./configure --prefix=$PWD/../install
make -j4
make install
cd ../..
```

In order to build the software, please execute the following lines in any directory with write access:
```bash
git clone https://github.com/SVfit/ClassicSVfit4tau TauAnalysis/ClassicSVfit4tau
export LIBRARY_PATH=$LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib:$PWD/VAMP/install/lib:$PWD/TauAnalysis/ClassicSVfit4tau/lib
make -f TauAnalysis/ClassicSVfit4tau/Makefile -j4
```

The test executables will be placed to `$PWD/TauAnalysis/ClassicSVfit4tau/exec`. In order to use them:
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib:$PWD/VAMP/install/lib:$PWD/TauAnalysis/ClassicSVfit4tau/lib

# either enter the full path to the executable, e.g.
./TauAnalysis/ClassicSVfit4tau/exec/testClassicSVfit4tau

# or make the executable available globally
export PATH=$PATH:$PWD/TauAnalysis/ClassicSVfit4tau/exec
testClassicSVfit4tau # run anywhere
```
You can add the export statements to your `$HOME/.bashrc` to make their effect permanent.

# Running instructions
- [Example(s)](https://github.com/SVfit/ClassicSVfit4tau/blob/master/bin/testClassicSVfit4tau.cc)

# Reference

If you use this code, please cite:                                                                                                    
```
   K. Ehataeht, L. Marzola, C. Veelken,
   "Reconstruction of the mass of Higgs boson pairs in events with Higgs boson pairs decaying into four tau leptons",
   arXiv:1809.06140
```
