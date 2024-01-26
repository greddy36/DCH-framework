__author__ = "Alexis Kalogeropoulos"
__description__ = "Simple script to make LaTex table from a given .root file - it needs arguments as described below. It will produce several txt a) per channel b) per process c) per group (ZZ, DY, etc). "

'''
example  : python Yields_new.py -f MCsamples_2018_ZH.csv -y 2018 -u yes -i ../plotting/allGroups_2018_OS_LT00_16noSV_16brute_nom.root -n no -s ZH
'''

import os
import sys
from ROOT import TFile, TTree, TH1, TH1D, TCanvas, TLorentzVector
from array import *
import numpy as np
from decimal import *
#getcontext().prec = 3
from math import sqrt
import csv
import pandas as pn
import glob
import string

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plotting
from ROOT import kBlack, kBlue, kMagenta, kOrange, kAzure, kRed

pn.set_option('display.float_format', lambda x: '%.2f' % x)


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'



def latex_with_lines(df, *args, **kwargs):
    kwargs['column_format'] = '|'.join([''] + ['l'] * df.index.nlevels
                                            + ['r'] * df.shape[1] + [''])
    res = df.to_latex(*args, **kwargs)
    res.replace('\\\\\n', '\\\\ \\hline\n')
    res.replace('bottomrule', 'hline')
    return res.replace('\\\\\n', '\\\\ \\hline\n')

def WriteLatexGroup (fileOut,array) :
    f = open(fileOut,'w')
    print >> f,'\\documentclass[10pt]{report}'
    print >> f,'\\usepackage{adjustbox}'
    print >> f,'\\begin{document}'
    print >> f,'\\begin{table}[htp]'
    print >> f,'\\caption{' + '{0:s} {1:s}'.format(args.selection, era) + '}'  
    print >> f,'\\begin{center}'
    #print >>f, cc.to_latex()
    f.write(array.to_latex(column_format = "l | " + " | ".join(["r"] * len(array.columns))))
    print >> f,'\\end{center}'
    print >> f,'\\end{table}'
    print >> f,'\\end{document}'
    f.write(fileOut)
    f.close()

    s = open(fileOut).read()
    s = s.replace('\\bottomrule','')
    s = s.replace('\\midrule','')
    s = s.replace('\\toprule','')
    s = s.replace('\\\\','\\\\ \\hline')
    fa = open(fileOut,'w')
    fa.write(s)
    fa.close()


def column(matrix, i):
    return [row[i] for row in matrix]


def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-v","--verbose",default=0,type=int,help="Print level.")
    parser.add_argument("-f","--inFile",default='MCsamples_2016.csv',help="Input file name.")
    parser.add_argument("-y","--year",default='2016',type=str,help="Data taking period, 2016, 2017 or 2018")
    parser.add_argument("-s","--selection",default='ZH',type=str,help="Select ZH or AZH")
    parser.add_argument("-c","--category",default='all',type=str,help="Categories")
    parser.add_argument("-n","--normalize",default='yes',type=str,help="scale to lumi, if no/0 then pure yields, otherwise scale to 35.9/41.5/59.7 for 2016/2017/2018")
    parser.add_argument("-u","--usefromplots",default='no',type=str,help="what .root files to use. if is it yes, then it will look in the plotting dir for the file defined as with the -i argument")
    parser.add_argument("-i","--inplotsfile",default='allGroups_2016_OS_LT00.root',type=str,help="define the .root files to use for making the tables. Must be used with -u switch")
    return parser.parse_args()


args = getArgs()

incat=args.category

cats = [ 'eeet', 'eemt', 'eett', 'eeem', 'mmet', 'mmmt', 'mmtt', 'mmem']
if incat !='all' : 
    cats =[]
    cats.append(incat)

cats.insert(0,'Cuts')

era=str(args.year)

command = "mkdir txt"
os.system(command)


if args.category not in cats and incat != 'all' : sys.exit("There is not such channel --> {0:s} <-- ...".format(incat))

fInplot = args.inplotsfile
usePlotFile= args.usefromplots



lumi = {'2016':35920, '2017':41530, '2018':59740}

nickNames, xsec, totalWeight, sampleWeight = {}, {}, {}, {}
#groups = ['Signal','DY','Other','ZZ','data']
groups = ['Signal','ZZ','DY','Other','Top','data']
groups = ['Top','ZZ','WZ','DY','Other','Signal','data']


fIn=str(args.inFile)

if fIn == '' : fIn ='./MCsamples_'+era+'.csv'


for group in groups :
    nickNames[group] = []


#for line in open('./MCsamples_'+era+'_one.csv','r').readlines() :
for line in open(fIn,'r').readlines() :
    vals = line.split(',')
    if '#' in vals[0] : continue
    nickName = vals[0]
    group = vals[1]
    nickNames[group].append(nickName)
    xsec[nickName] = float(vals[2])
    totalWeight[nickName] = float(vals[4])
    sampleWeight[nickName]= float(lumi[era])*xsec[nickName]/totalWeight[nickName]
    print("group={0:10s} nickName={1:20s} xSec={2:10.3f} totalWeight={3:11.1f} sampleWeight={4:10.6f}".format(
        group,nickName,xsec[nickName],totalWeight[nickName],sampleWeight[nickName]))

# Stitch the DYJets and WJets samples
'''
for i in range(1,5) :
    nn = 'DY{0:d}JetsToLL'.format(i)
    if 'DYJetsToLL' in totalWeight : sampleWeight[nn] = float(lumi[era])/(totalWeight['DYJetsToLL']/xsec['DYJetsToLL'] + totalWeight[nn]/xsec[nn])

for i in range(1,4) :
    nn = 'W{0:d}JetsToLNu'.format(i)
    if 'WJetsToLNu' in totalWeight : sampleWeight[nn] = float(lumi[era])/(totalWeight['WJetsToLNu']/xsec['WJetsToLNu'] + totalWeight[nn]/xsec[nn])
'''

# now add the data - presumably it should the merged data
#for eras in ['2017B','2017C','2017D','2017E','2017F'] :
    #for dataset in ['SingleElectron','SingleMuon','DoubleEG','DoubleMuon'] :
#for dataset in ['SingleMuon', 'SingleElectron'] :
for dataset in ['data'] :
    #nickName = '{0:s}_Run{1:s}'.format(dataset,era)
    #nickName = '{0:s}_{1:s}'.format(dataset,era)
    nickName = 'data'
    totalWeight[nickName] = 1.
    sampleWeight[nickName] = 1.
    nickNames['data'].append(nickName)
    xsec[nickName] = 1.
#####################################3



channel=args.category
fIn=''
fIn_data = ''
fIn_mc = ''

b_usePlotFile = False

ccols=19 ##this must be equal to the nBins from the CutFlow histo
rrows=9


if usePlotFile.lower() == 'yes' : b_usePlotFile = True

for group in groups :
    rows, cols = (rrows, ccols) 
     
    
    for nickName in nickNames[group]:
	if b_usePlotFile:
	    fIn_data = '../plotting/'+fInplot
	    fIn_mc = '../plotting/'+fInplot
	else : 
	    fIn_data = '../data/condor/{0:s}/{1:s}_{2:s}/{1:s}_{2:s}.root'.format(args.selection, nickName, str(era))
	    fIn_mc = '../MC/condor/{0:s}/{1:s}_{2:s}/{1:s}_{2:s}.root'.format(args.selection, nickName, era)

	if group != 'data' : fIn = fIn_mc  #'{0:s}/{1:s}_{2:s}/{1:s}_{2:s}.root'.format(args.selection, nickName, era)
	else: fIn = fIn_data # '../../data/condor/{0:s}/{1:s}/{1:s}.root'.format(args.selection, nickName)
	inFile = TFile.Open(fIn)
	inFile.cd()
	#hW=''
	#if group != 'data' : hW = inFile.Get("hWeights")
	#inTree = inFile.Get("Events")
	
    header=nickName
    print ' opening ', fIn, nickName, group

    scale=str(args.normalize)
    ScaleToLumi= False

    if scale=="1" or scale.lower()=="true" or scale.lower() =="yes":  
	ScaleToLumi = True
	header=header+'_'+'{0:.3f}'.format(float(lumi[era]/1000))+'invfb'

    if ScaleToLumi : print "Events scaled to {0:.1f}/pb for xsec = {1:.3f} pb".format(float(lumi[era]),float(xsec[nickName]))
    if not ScaleToLumi : print "Events WILL NOT BE scaled to {0:.1f}/pb for xsec = {1:.3f} pb".format(float(lumi[era]),float(xsec[nickName]))



    arr = [[0 for i in range(cols)] for j in range(rows)] 
    cuts=[]
    cuts=['All', 'inJSON', 'METFilter',  'Trigger', 'LeptonCount', 'GoodLeptons', 'LeptonPair', 'FoundZ', 'GoodTauPair']
    if  b_usePlotFile : cuts=[]

    result=[]
    product=1.
    if ScaleToLumi : product = float(sampleWeight[nickName])
		    

    count=0
    yieldsdata=[]
    yieldpergroup = [[0 for i in range(cols)] for j in range(rows)] 
    for icat,cat in enumerate(cats[1:]):
	hist=''
	#yieldpergroup = [[0 for i in range(cols)] for j in range(rows)] 
	totalyield = [[0 for i in range(cols)] for j in range(rows)]
	if not b_usePlotFile :
	    if 'data' not in group : hist = "hCutFlowWeighted_"+cat
	    else : hist = "hCutFlow_"+cat

	else : 
	    hist = "hCutFlowPerGroup_"+group+'_'+cat


		#else : hist="hCutFlow_"+cat

	htest = inFile.Get(hist)
	#print 'will work on ',htest, hist
	for i in range(1,htest.GetNbinsX()) : 
	    #print ' test---------------', htest.GetXaxis().GetBinLabel(i), htest.GetBinContent(i), cat, nickName, htest.GetName()
	    arr[icat][i-1] = '{:.2f}'.format(htest.GetBinContent(i)*product)
	    #print '===================',arr[icat][i-1]
	    yieldpergroup[icat][i-1] +=  float('{0:.2f}'.format(htest.GetBinContent(i)*product))
	    #arr[icat][i-1] = 1
	    if len(cuts) == 0 :
		for i in range(1,htest.GetNbinsX()) : 
		    if '>' in htest.GetXaxis().GetBinLabel(i) : cuts.append('H_LT_gt_0')
		    elif htest.GetXaxis().GetBinLabel(i) != '' : cuts.append(htest.GetXaxis().GetBinLabel(i))
		    else : cuts.append('Cut')
	count+=1

	### this is to create a .txt per group

    with open('txt/All_'+group+'_'+args.selection+'_'+era+'_yields.txt', 'w') as f:
	#lines = [' & \t'.join([str(x[i]) if len(x) > i else ' ' for x in yieldpergroup]) for i in range(len(max(yieldpergroup)))]
	#lines = ['  \t'.join([str(x[i]) if len(x) > i else ' ' for x in yieldpergroup]) for i in range(len(max(yieldpergroup)))]
	lines = ['  \t'.join([str(x[i]) for x in yieldpergroup]) for i in range(0,len(max(yieldpergroup)))]
	#f.write("\n".join(["\t".join([str(groupp[index]) for groupp in yieldpergroup]) for index, entry in enumerate(range(0,len(max(yieldpergroup))))])) 
	#lines=["\n".join(["\t".join([str(groupp[index]) for groupp in yieldpergroup]) for index, entry in enumerate(range(0,len(max(yieldpergroup))))])]
	#lines = ["\n".join(["\t".join([str(groupp[index]) for groupp in yieldpergroup]) for index, entry in enumerate(range(0,len(max(yieldpergroup))))])]


	for i in range(len(cuts)):
	    #if cuts[i] != '' : 
	    #    print >> f,'{} & {} \\\\ \\hline'.format(cuts[i], lines[i])
	    #if len(cuts[i]) != 0 : 
	    #print cuts[i], lines[i]
	    #if cuts[i] == '' and group == 'data' : cuts[i] = 'NextCut'
	    print >> f, '{0:s}  {1:s} '.format(cuts[i], lines[i])

    #for i in range(rows):
    #    for j in range(columns):
top=[]
zz=[]
other=[]
signal=[]
Data=[]
cutlist=[]
allbkg=[]
e_allbkg=[]
wz=[] 
dy =[]

 
cutlist.append( np.genfromtxt('txt/All_Signal_{0:s}_{1:s}_yields.txt'.format(args.selection, args.year), dtype=None,usecols=(0)))


for counter in range( len(cats[1:])) :
#for counter in range(1,len(cats)) :
    #print counter, cats[counter]
    top.append( np.genfromtxt('txt/All_Top_{0:s}_{1:s}_yields.txt'.format(args.selection, args.year), dtype=None,usecols=(counter+1)))
    dy.append( np.genfromtxt('txt/All_DY_{0:s}_{1:s}_yields.txt'.format(args.selection, args.year), dtype=None,usecols=(counter+1)))
    wz.append( np.genfromtxt('txt/All_WZ_{0:s}_{1:s}_yields.txt'.format(args.selection, args.year), dtype=None,usecols=(counter+1)))
    zz.append( np.genfromtxt('txt/All_ZZ_{0:s}_{1:s}_yields.txt'.format(args.selection, args.year), dtype=None,usecols=(counter+1)))
    other.append( np.genfromtxt('txt/All_Other_{0:s}_{1:s}_yields.txt'.format(args.selection, args.year), dtype=None,usecols=(counter+1)))
    Data.append( np.genfromtxt('txt/All_data_{0:s}_{1:s}_yields.txt'.format(args.selection, args.year), dtype=None,usecols=(counter+1)))
    signal.append( np.genfromtxt('txt/All_Signal_{0:s}_{1:s}_yields.txt'.format(args.selection, args.year), dtype=None,usecols=(counter+1)))
    allbkg.append( top[counter] + dy[counter] +  wz[counter] +   zz[counter] + other[counter]) 
    #counter+=1

#data = {'names': cutlist, 'values': zz}

dftop = pn.DataFrame(data=top)
df2top_t = dftop.T
df2top_t.index=[i for i in cutlist]
df2top_t.columns=[i for i in cats[1:]]


dfdy = pn.DataFrame(data=dy)
df2dy_t = dfdy.T
df2dy_t.index=[i for i in cutlist]
df2dy_t.columns=[i for i in cats[1:]]

dfzz = pn.DataFrame(data=zz)
df2zz_t = dfzz.T
df2zz_t.index=[i for i in cutlist]
df2zz_t.columns=[i for i in cats[1:]]


dfwz = pn.DataFrame(data=wz)
df2wz_t = dfwz.T
df2wz_t.index=[i for i in cutlist]
df2wz_t.columns=[i for i in cats[1:]]


dfother = pn.DataFrame(data=other)
df2other_t = dfother.T
df2other_t.index=[i for i in cutlist]
df2other_t.columns=[i for i in cats[1:]]
#df2other_t.head()


dfallbkg = pn.DataFrame(data=allbkg)
df2allbkg_t = dfallbkg.T
df2allbkg_t.index=[i for i in cutlist]
df2allbkg_t.columns=[i for i in cats[1:]]
#df2allbkg_t.head()

dfdata = pn.DataFrame(data=Data)
df2data_t = dfdata.T
df2data_t.index=[i for i in cutlist]
df2data_t.columns=[i for i in cats[1:]]
#df2data_t.head()

dfsignal = pn.DataFrame(data=signal)
df2signal_t = dfsignal.T
df2signal_t.index=[i for i in cutlist]
df2signal_t.columns=[i for i in cats[1:]]

df2zz_t.index=[i for i in cutlist]
df2zz_t.columns=[i for i in cats[1:]]

df2wz_t.index=[i for i in cutlist]
df2wz_t.columns=[i for i in cats[1:]]
df2dy_t.index=[i for i in cutlist]
df2dy_t.columns=[i for i in cats[1:]]
#df2zz_t.head()


e_allbkg = pn.DataFrame(data=allbkg).pow(1./2)

e_dfallbkg = pn.DataFrame(data=e_allbkg)
e_df2allbkg_t = e_dfallbkg.T
e_df2allbkg_t.index=[i for i in cutlist]
e_df2allbkg_t.columns=[i for i in cats[1:]]



c=pn.concat([df2top_t, df2zz_t, df2wz_t, df2dy_t, df2other_t, df2allbkg_t, e_df2allbkg_t, df2data_t], axis=1)

sizes=[]


def func(pct, allvals):
    absolute = float(pct/100.*np.sum(allvals))
    return "{:.1f}\n({:.1f}%)".format(pct, absolute)

#colors = {kAzure-9, kMagenta-10,kBlue-8, 'lightcoral'}
plt.ioff()
for icat, cat in enumerate(cats[1:]) :
    #for i in range(len(cuts)-5) :
        #fig1, ax1 = plt.subplots()
        fig1, ax1 = plt.subplots(figsize=(5, 5), subplot_kw=dict(aspect="equal"))
        
        labels=['Other','DY','WZ','ZZ']
        #labels=['Obs %.1f'%Data[icat][12],   'allbkg %2.f'%allbkg[icat][12]]
                
        #    'DY','WZ','ZZ']
	sizes = [other[icat][12]+top[icat][12] , dy[icat][12],wz[icat][12] ,zz[icat][12] ]
        print '-->', sizes, cat
	#colors = ['gold','Orange', 'ForestGreen', 'LightCoral','Cyan']
	colors = ['gold','limegreen', 'LightCoral','Cyan']
        lab=["","","",""]
        #patches, texts, autotexts  = ax1.pie(sizes, colors = colors, labels=labels, autopct=lambda pct: func(pct, sizes),textprops=dict(color="black"), startangle=90, radius=1)
        patches, texts, autotexts  = ax1.pie(sizes, colors = colors, labels=labels, autopct=lambda pct: func(pct, sizes),textprops=dict(color="black"), startangle=90, radius=1)
        #patches, texts, autotexts  = ax1.pie(sizes, colors = colors,  autopct=lambda pct: func(pct, sizes),textprops=dict(color="black"), startangle=90, radius=1)
        
        plt.setp(autotexts, size=10, weight="bold")
        plt.text(0.6,1, 'Obs %1.f'%Data[icat][12]+'\nallbkg %.1f'%allbkg[icat][12], bbox=dict(facecolor='lightblue', alpha=0.5), fontsize=15)  
	plt.legend(patches, labels, loc="best")
        ax1.set_title(cat)
        ax1.get_legend().remove()
	plt.tight_layout()
        plt.savefig("pie_{0:s}_{1:s}.png".format(cat,era))

command="montage -auto-orient -title ZH_{0:s} -tile 2x4 -geometry +5+5 -page A4 pie_*_{0:s}.png Multi_pie_{0:s}.pdf;".format(str(era))

os.system(command)



groupss = ['Top','ZZ','WZ','DY','Other','Allbkg','data']
ngroup = groupss
#ngroup.insert(5,'All bkg')
ngroup.insert(6,'uncert bkg')
for cat in cats[1:]:
    cc=c[cat]
    #print cc, cat
    cc.to_latex()
    cc.columns=[i for i in ngroup]
    cc.index=[i for i in cutlist]

    fIn = 'txt/All_'+cat+'_'+args.selection+'_'+era+'yields.txt'

    with open(fIn, 'w') as f:
	print >> f,'\\documentclass[10pt]{report}'
	print >> f,'\\usepackage{adjustbox}'
	print >> f,'\\begin{document}'
	print >> f,'\\begin{table}[htp]'
	print >> f,'\\caption{' + '{0:s} {1:s} {2:s}'.format(cat, args.selection, era) +'}'  
	print >> f,'\\begin{center}'
        '''
	print >> f,'\\begin{adjustbox}{width=1\\textwidth}'
	#print >> f,'\\begin{tabular}{l r r r r r r r r }  \hline'
        f.write("\\begin{tabular}{l | " + " | ".join(["r"] * len(cc.columns)) + "} \\\hline \n")
        #f.write("\\begin{tabular}{l | " + " | ".join(["r"] * len(cc.head())) + "}\n")
        count=0
        for i, row in cc.iterrows():
            i_ = cc.index[count]
            f.write(" & ".join([str(x) for x in row.values]) + " \\\\\n")
            count+=1
	print >> f, '{} &'.format(i for i in cc.columns)
	print >> f, '\hline'
        '''
        f.write(cc.to_latex(column_format = "l | " + " | ".join(["r"] * len(cc.columns))))
        print >> f,'\\end{center}'
        print >> f,'\\end{table}'
        print >> f,'\\end{document}'

        '''
	print >> f,'\\end{tabular}'
	print >> f,'\\end{adjustbox}'
        '''


    s = open(fIn).read()
    s = s.replace('\\bottomrule','')
    s = s.replace('\\midrule','')
    s = s.replace('\\toprule','')
    s = s.replace('\\\\','\\\\ \\hline')
    f = open(fIn,'w')
    f.write(s)
    f.close()

WriteLatexGroup('txt/All_Top_'+args.selection+'_'+era+'_yields.txt',df2top_t)
WriteLatexGroup('txt/All_ZZ_'+args.selection+'_'+era+'_yields.txt',df2zz_t)
WriteLatexGroup('txt/All_WZ_'+args.selection+'_'+era+'_yields.txt',df2wz_t)
WriteLatexGroup('txt/All_DY_'+args.selection+'_'+era+'_yields.txt',df2dy_t)
WriteLatexGroup('txt/All_Other_'+args.selection+'_'+era+'_yields.txt',df2other_t)
WriteLatexGroup('txt/All_data_'+args.selection+'_'+era+'_yields.txt',df2data_t)
WriteLatexGroup('txt/All_Allbkg_'+args.selection+'_'+era+'_yields.txt',df2allbkg_t)
WriteLatexGroup('txt/All_Signal_'+args.selection+'_'+era+'_yields.txt',df2signal_t)






command="sed -i '/0.0 &  /d' txt/*txt"
command1="sed -i 's/0.00 \\\\/0 \\\\/g' txt/*txt"
command2="sed -i '/Cut & /d' txt/*txt"
command3="sed -i '11,20d' txt/*txt"
if  b_usePlotFile : 
    os.system(command)
os.system(command1)
os.system(command2)
os.system(command3)

print 'All txt files can be found in the ./txt dir....exiting'
