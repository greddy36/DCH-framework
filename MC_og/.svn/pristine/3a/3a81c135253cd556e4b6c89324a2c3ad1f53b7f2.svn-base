
import os

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-v","--verbose",default=0,type=int,help="Print level.")
    defDS = '/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM '
    parser.add_argument("--dataSet",default=defDS,help="Data set name.") 
    parser.add_argument("--nickName",default='MCpileup',help="Data set nick name.") 
    parser.add_argument("-m","--mode",default='anaXRD',help="Mode (script to run).")
    parser.add_argument("-y","--year",default=2017,type=str,help="Data taking period, 2016, 2017 or 2018")
    parser.add_argument("-c","--concatenate",default=3,type=int,help="On how many files to run on each job")
    parser.add_argument("-s","--selection",default='ZH',type=str,help="select ZH or AZH")
    parser.add_argument("-j","--doSystematics",default='yes',type=str,help="do JME systematics")
    parser.add_argument("-l","--islocal",default='no',type=str,help="get list from /eos/ not DAS")
    parser.add_argument("-w","--weights",default='1',type=str,help="store genWeights in a TH1")
    return parser.parse_args()

def beginBatchScript(baseFileName, Systematics) :
    outLines = ['#!/bin/bash\n']
    outLines.append("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
    outLines.append("export SCRAM_ARCH=slc7_amd64_gcc820\n")
    outLines.append("eval `scramv1 project CMSSW CMSSW_10_6_5`\n")
    outLines.append("cd CMSSW_10_6_5/src\n")
    outLines.append("eval `scramv1 runtime -sh`\n")
    outLines.append("cmsenv\n")
    #outLines.append("git clone --recursive https://github.com/cms-nanoAOD/correctionlib.git\n")
    outLines.append("cd ${_CONDOR_SCRATCH_DIR}/CMSSW_10_6_5/src/\n")
    #outLines.append("git clone https://github.com/cms-tau-pog/TauIDSFs TauPOG/TauIDSFs\n")
    outLines.append("scram b -j 4\n")
    outLines.append("echo ${_CONDOR_SCRATCH_DIR}\n")
    outLines.append("cd ${_CONDOR_SCRATCH_DIR}/CMSSW_10_6_5/src/\n")
    outLines.append("cp ${_CONDOR_SCRATCH_DIR}/* .\n")
    outLines.append("tar -zxvf correctionlib.tar.gz\n")
    outLines.append("ls -altrh\n")
    outLines.append("echo 'this is the working dir' ${_CONDOR_SCRATCH_DIR}\n")
    return outLines

def getFileName(line) :
    tmp = line.split()[0].strip(',')
    fileName = tmp.strip()
    return fileName


args = getArgs()
era = str(args.year)
doJME  = args.doSystematics.lower() == 'true' or args.doSystematics.lower() == 'yes' or args.doSystematics == '1'

#weightSwitch = args.weights.lower() =='yes' or args.weights.lower()=='1' or  args.weights.lower() == 'true'

period="B"
if 'Run2016' in args.dataSet or 'Run2017' in args.dataSet or 'Run2018' in args.dataSet: 
    poss = args.dataSet.find("Run")
    period = args.dataSet[int(poss)+7:int(poss)+8]
    print 'will set up', poss, period

  

# sample query 
# dasgoclient --query="file dataset=/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8*/*/NANOAOD*" --limit=0   


query = '"file dataset={0:s}"'.format(args.dataSet)
if "USER" in str(args.dataSet) : query = '"file dataset={0:s}"'.format(args.dataSet+" instance=prod/phys03")
command = "dasgoclient --query={0:s} --limit=0  > fileList.txt".format(query)


#if str(args.islocal.lower())=='yes' : command ='ls  /eos/uscms/store/group/lpcsusyhiggs/ntuples/nAODv7/ZH_JECs_{0:s}/CRAB_PrivateMC/{1:s}_{0:s}/*/*/*root   > fileList.txt'.format(str(args.year), args.nickName) 

whate = str(args.year)
getNt = args.nickName
executable = str(args.selection)
eraD=era
if 'pre' in whate or 'preV' in era or 'pre' in getNt or 'preV' in str(args.dataSet):
    whate='2016'
    getNt=args.nickName+"_preVFP"
    eraD=eraD+'preVFP'
if str(args.islocal.lower())=='yes' : command ='ls  /eos/uscms/store/group/lpcsusyhiggs/ntuples/nAODv9/{0:s}/{1:s}/*root   > fileList.txt'.format(str(whate), getNt) 
if str(args.islocal.lower())=='yes' and executable == 'Gjets' : command ='ls  /eos/uscms/store/group/lpcsusyhiggs/ntuples/nAODv9/Gjets/{0:s}/{1:s}/*root   > fileList.txt'.format(str(whate), getNt) 
#if str(args.islocal.lower())=='yes' or str(args.islocal.lower())=='1': command ='ls  /eos/uscms/store/group/lpcsusyhiggs/ntuples/nAODv7/JEC_{0:s}/*/{1:s}_{0:s}/*/*/*root   > fileList.txt'.format(str(args.year), args.nickName)

print("Running in {0:s} mode.  Command={1:s}".format(args.mode,command))
os.system(command)
    
files = open('fileList.txt','r').readlines()
if len(files) < 1 :
    print("***In makeCondor.py: Empty fileList.txt")
    exit()

scriptList = [] 
file=[]
dataset=[]

mjobs=args.concatenate


for nFiles, file in enumerate(files) :
     
    fileName=getFileName(file)
    if str(args.islocal.lower())=='yes' : fileName = fileName.replace('/eos/uscms','')
    if '#' not in fileName :  dataset.append(fileName)


counter=0

dirout = str(args.selection)

runLocal = str(args.islocal.lower())=='yes'


for nFile in range(0, len(dataset),mjobs) :
    print("nFile={0:d} file[:80]={1:s}".format(nFile,file[:80]))
    fileName = getFileName(file)
    print dataset[nFile]
    scriptName = "{0:s}_{1:03d}_{2:s}.sh".format(args.nickName,nFile+1, args.year)


    print("scriptName is ={0:s}".format(scriptName))
    outLines = beginBatchScript(scriptName,doJME)

    outLines.append("cp cuts_{0:s}_{1:s}.yaml cuts.yaml\n".format(args.selection, args.year))


    maxx = mjobs
    if counter+mjobs > len(dataset) : 
	print 'should include', nFile, -nFile-mjobs + len(dataset)+1, 'from ', len(dataset), counter
	maxx = len(dataset)-counter
	#for j in range(0,mjobs) :
    for j in range(0,maxx) :
	print 'shoud see', nFile+maxx, maxx, len(dataset)
	fileloop=dataset[nFile:nFile+maxx][j]
        print 'file------------------', fileloop
	#query = '"file={0:s} | grep file.nevents"'.format(fileloop)
	#command = "dasgoclient --query={0:s} ".format(query)
	#else : 
	infile = "inFile.root"

        outLines.append("xrdcp  root://cmseos.fnal.gov/{0:s} {1:s}\n".format(fileloop, infile)) 
        outLines.append("if [ ! -f inFile.root ] ; \n then \n xrdcp root://cmsxrootdfnal.gov/{0:s} {1:s}\n fi \n".format(fileloop, infile))
        outLines.append("if [ ! -f inFile.root ] ; \n then \n xrdcp root://cms-xrd-global.cern.ch/{0:s} {1:s}\n fi \n".format(fileloop, infile))
        


	outFileName = "{0:s}_{1:03d}.root".format(args.nickName,nFile+j+1)


	#if doJME : 

        fend = fileloop.rsplit('/', 1)[-1]

	#if weightSwitch and not runLocal : 
	    #outLines.append("python SystWeights.py -f {4:s}  -o {0:s} --nickName {1:s} -y {2:s} -s {3:s} -w 1 -j {5:s}\n".format(outFileName,args.nickName, args.year, args.selection,infile, args.doSystematics, executable, str(weightSwitch)))

	if not runLocal : 
	    outLines.append("python {6:s}.py -f {4:s} -o {0:s} --nickName {1:s} -y {2:s} -s {3:s} -w 1 -j {5:s}\n".format(outFileName,args.nickName, args.year, args.selection,infile, args.doSystematics, executable ))
    
	outLines.append("rm inFile*.root\n")

        ntupfile = outFileName.replace(".root", ".ntup")
        weightsfile = outFileName.replace(".root", ".weights")
 
	outLines.append("xrdcp {5:s} root://cmseos.fnal.gov//store/group/lpcsusyhiggs/ntuples/nAODv9/{6:s}_out/{0:s}/{7:s}\n".format(args.nickName+'_'+eraD, nFile+1+j, executable, args.year, fend, ntupfile, dirout, outFileName))
	outLines.append("\n")

    outLines.append("rm  *root *.ntup  *.so\nrm *.pcm\nrm *cc.d\n")
    outLines.append("cd ../.. ; rm -fr CMSSW_10_6_5\n")
    fend = fend.replace(".root","")
    outFileNameNoRoot = outFileName.replace(".root","")
    scriptName = "all_{0:s}.sh".format(outFileNameNoRoot)

    print("Writing out file = {0:s}".format(scriptName))
    open(scriptName,'w').writelines(outLines)
    scriptList.append(scriptName)
    counter += mjobs


wdir='/uscms_data/d3/alkaloge/MetStudies/nAOD/CMSSW_10_6_5/src/'
wdir='/uscms_data/d3/alkaloge/DCH/CMSSW_10_6_5/src/DCH_KSU/'
dirMC = wdir+"/MC/"
dirCode = wdir+"/DCH/"
dirData = wdir+"/data/"
funcsDir = wdir+"/funcs/"
toolsDir = wdir+"/tools/"
pileupDir = wdir+"/pileup/"

print("dir={0:s}".format(dir))

for file in scriptList :
    base = file[:-3] 
    outLines = ['universe = vanilla\n']
    outLines.append('Executable = {0:s}\n'.format(file))
    outLines.append('Output = {0:s}.out\n'.format(base))
    outLines.append('Error = {0:s}.err\n'.format(base))
    outLines.append('Log = {0:s}.log\n'.format(base))
    #outLines.append('transfer_input_files = {0:s}ZH.py, {0:s}MC_{1:s}.root, {0:s}data_pileup_{1:s}.root, {0:s}MCsamples_{1:s}.csv, {0:s}ScaleFactor.py, {0:s}SFs.tar.gz, {0:s}cuts_{2:s}.yaml, '.format(dir,args.year, args.selection))
    outLines.append('transfer_input_files = {0:s}{1:s}.py, {0:s}SystWeights.py, {0:s}ScaleFactor.py,'.format(dirCode,executable))
    outLines.append('{0:s}pileup_{1:s}UL_MC.root, {0:s}pileup_{1:s}UL_data.root, {0:s}cuts_{2:s}_{1:s}.yaml, '.format(dirMC,args.year, args.selection, executable))
    #outLines.append('{0:s}*txt, '.format(dirData))
    outLines.append('{0:s}/METCorrections.py, {0:s}tauFunDCH.py, {0:s}generalFunctions.py, {0:s}outTuple.py, {0:s}Weights.py, '.format(funcsDir))
    outLines.append('{0:s}Electron_RunUL2016postVFP_Ele25_EtaLt2p1.root, {0:s}Electron_RunUL2016preVFP_Ele25_EtaLt2p1.root,  {0:s}Electron_RunUL2017_Ele35.root, {0:s}Electron_RunUL2018_Ele35.root, {0:s}muon_Z_{1:s}.json.gz, {0:s}electron_{1:s}.json.gz, {0:s}photon_{1:s}.json.gz, {0:s}correctionlib.tar.gz\n'.format(toolsDir, args.year))
    #outLines.append('{0:s}make_jme.py, {0:s}branchselection.py, {0:s}keep_and_drop.txt, {0:s}taupog.tar.gz\n'.format(toolsDir))
    outLines.append('priority = 2\n')
    outLines.append('should_transfer_files = YES\n')
    outLines.append('when_to_transfer_output = ON_EXIT\n')
    outLines.append('x509userproxy = $ENV(X509_USER_PROXY)\n')
    outLines.append('Queue 1\n')
    open('{0:s}.jdl'.format(base),'w').writelines(outLines)


    
    
