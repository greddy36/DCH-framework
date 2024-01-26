
# generate a runMC.csh script that creates the .csh and .jdl files
# to process MC data 

import os

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--inFile",default='MCsamples_2016.csv',help="Input file name.") 
    parser.add_argument("-y","--year",default=2017,type=str,help="Data taking period, 2016, 2017 or 2018")
    parser.add_argument("-s","--selection",default='ZH',type=str,help="Select analyser, like ZH, DCH etc")
    parser.add_argument("-j","--doSystematics",default='yes',type=str,help="do JME systematics")
    parser.add_argument("-l","--islocal",default='no',type=str,help="local /eos files or DBS")
    return parser.parse_args()

args = getArgs() 
era=str(args.year)
outLines = []
cwd = os.getcwd()
conc=20
if str(args.islocal.lower())=='yes' or str(args.islocal.lower())=='1' or str(args.islocal.lower())=='true': conc = 1
for line in open(args.inFile,'r').readlines() :
     
    nickname = line.split(',')[0]
    if '#' in nickname : continue

    #print("\n\n\n line.split(',')={0:s}".format(str(line.split(','))))
    dataset = line.split(',')[6].replace(' ','_').strip()
    if 'NANO' in dataset : conc=5
    if len(dataset) < 2 : continue
    #print("\n***line.split()={0:s}".format(str(line.split(','))))
    print("nickname={0:s} \n dataset={1:s}".format(nickname,dataset))

    mode = 'anaXRD'
    eraD=era
    if 'pre' in str(nickname) or 'pre' in str(line) : eraD=eraD+'preVFP'
    outLines.append("mkdir -p {0:s}/{1:s}_{2:s}\ncd {0:s}/{1:s}_{2:s}\n".format(args.selection,nickname,eraD))

    outLines.append("python ../../makeCondor.py --dataSet {0:s} --nickName {1:s} --mode {2:s} --year {3:s} -c {5:s} -s {4:s} -j {6:s} -l {7:s} -w 1\n".format(dataset,nickname, mode,era, args.selection, str(conc), args.doSystematics, str(args.islocal)))
    outLines.append("cd {0:s}\n".format(cwd))

fOut='runMC_{0:s}_{1:s}.csh'.format(str(args.year),args.selection)
open(fOut,'w').writelines(outLines)



    
    
