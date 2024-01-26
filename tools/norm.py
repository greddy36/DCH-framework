import ROOT
from ROOT import gSystem, gStyle, gROOT, kTRUE, gDirectory
from ROOT import TFile, TH1F, TH1D, TTree
import os
import os.path
import sys


def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--inFileName",default='file.csv',type=str,help="Data taking period, 2016, 2017 or 2018")
    parser.add_argument("-y","--year",default=2017,type=int,help="Year for data.")
    parser.add_argument("-o","--overideSyst",default='',type=str,help="overide systematics list")
    return parser.parse_args()

args = getArgs()

era=str(args.year)

scaleSyst = ["Central"]

scale = ['scale_e', 'scale_m_etalt1p2', 'scale_m_eta1p2to2p1', 'scale_m_etagt2p1',
'scale_t_1prong', 'scale_t_1prong1pizero', 'scale_t_3prong', 'scale_t_3prong1pizero']

for i, sys in enumerate(scale) :
    scaleSyst.append(sys+'Up')
    scaleSyst.append(sys+'Down')


#listsyst=['njets', 'nbtag', 'jpt', 'jeta', 'jflavour','MET_T1_pt', 'MET_T1_phi', 'MET_pt', 'MET_phi', 'MET_T1Smear_pt', 'MET_T1Smear_phi']

jes=['jesAbsolute', 'jesAbsolute_{0:s}'.format(str(era)), 'jesBBEC1', 'jesBBEC1_{0:s}'.format(str(era)), 'jesEC2', 'jesEC2_{0:s}'.format(str(era)), 'jesFlavorQCD', 'jesHF', 'jesHF_{0:s}'.format(str(era)), 'jesRelativeBal', 'jesRelativeSample_{0:s}'.format(str(era)), 'jesHEMIssue', 'jesTotal', 'jer']

jesSyst=[]
for i, sys in enumerate(jes) :
    jesSyst.append(sys+'Up')
    jesSyst.append(sys+'Down')


otherS=['NLOEWK','PreFire','tauideff_pt20to25', 'tauideff_pt25to30', 'tauideff_pt30to35', 'tauideff_pt35to40', 'tauideff_ptgt40','scale_met_unclustered', 'scale_lowpt', 'scale_highpt']
otherS=['PreFire','tauideff_pt20to25', 'tauideff_pt25to30', 'tauideff_pt30to35', 'tauideff_pt35to40', 'tauideff_ptgt40','scale_met_unclustered' ]
OtherSyst=[]
for i, sys in enumerate(otherS) :
    OtherSyst.append(sys+'Up')
    OtherSyst.append(sys+'Down')

#jesSyst hold the jes systematic  - we need njets, nbtag, nflavor, jpt

sysall = scaleSyst+jesSyst+OtherSyst
sysall = scaleSyst

if str(args.overideSyst) != '' : sysall=[str(args.overideSyst)]

#sysall = ["Central"]

nickNames, xsec, totalWeight, sampleWeight = {}, {}, {}, {}

totalW=[]

for line in open(args.inFileName,'r').readlines() :
    vals = line.split(',')
    if '#' in vals[0] : continue

    nickName = vals[0]
    #group = vals[1]
    #nickNames[group].append(nickName)
    if '*' in vals[2] :
        value1, value2 = map(float, vals[2].split("*"))
        xsec[nickName] = float(value1*value2)
    elif '+' in vals[2] :
        value1, value2 = map(float, vals[2].split("+"))
        xsec[nickName] = float(value1+value2)
    else : xsec[nickName] = float(str(vals[2]))

    #if '+' in vals[4] :
    #    value1, value2 = map(float, vals[4].split("+"))
    #	totalWeight[nickName] = float(value1+value2)
    #else : totalWeight[nickName] = float(vals[4])
	#sampleWeight[nickName]= Pblumi*weights['lumi']*xsec[nickName]/totalWeight[nickName]
    filein = '/eos/uscms/store/user/alkaloge/ZH/nAODv7/{0:s}/{1:s}_{0:s}/{1:s}_{0:s}.root'.format(era,vals[0])
    print filein
    fIn = TFile.Open(filein,"READ")
    totalW.append(float(fIn.Get("hWeights").GetSumOfWeights()))
    totalWeight[nickName] = float(fIn.Get("hWeights").GetSumOfWeights())



print totalWeight, sum(totalWeight.values()), sum(totalW)


listF=[]
listO=[]
for line in open(args.inFileName,'r').readlines() :
    vals = line.split(',')
    nickName = vals[0]
    for sys in sysall : 
        fin = '{0:s}_{1:s}_sys{2:s}.root'.format(vals[0],era, sys)
        fout = '{0:s}_{1:s}_norm_sys{2:s}.root'.format(vals[0],era, sys)
        listF.append(fin)
        listO.append(fin)

	cf = os.path.isfile('../{0:s}'.format(fin))
	if cf : 
	    command="mv ../{0:s} . ".format(fin)
            os.system(command)
            #print command

        fIn = TFile.Open(fin,'read')
        fOut = TFile.Open(fout,'recreate')
        #print 'will work for', fin, fout
        try :
	    fIn.cd()

	    dirList2 = gDirectory.GetListOfKeys()
	    for k2 in dirList2:
		h2 = k2.ReadObj()
		htest=h2.Clone()

		#if 'mmtt_m_sv_new_FMjall' in htest.GetName() : print 'for ', fin, nickName, totalWeight[nickName], sum(totalWeight.values()) , 'ratio of weights', totalWeight[nickName]/sum(totalWeight.values()), htest.GetSumOfWeights(), htest.GetSumOfWeights()* (totalWeight[nickName]/sum( (totalWeight.values()))), totalWeight[nickName]/sum( (totalWeight.values())) , htest.GetBinContent(3), htest.GetBinContent(3)*(totalWeight[nickName]/sum( (totalWeight.values())))
		
		htest.Scale(float(totalWeight[nickName]))
		htest.Scale(1/sum( (totalWeight.values())))

		fOut.cd()
		htest.Write()

        except ReferenceError: print 'PROBLEM !!!!!!!!!!!! CHECK THIS ONE', fin
	fOut.Write()
	fOut.Close()
	command="cp {0:s} ../.".format(fout)
	print command
	os.system(command)
   


#print listF, listO
