__author__ = "Alexis Kalogeropoulos"
__description__ = "Simple script to make LaTex table from a given .root file - it needs 3 arguments 1:input file 2:channel (all is preferred) 3:scale to lumi (0, 1, True, False)"


import sys
from ROOT import TFile, TTree, TH1, TH1D, TCanvas, TLorentzVector
from array import *
import numpy as np
from decimal import *
getcontext().prec = 2

def column(matrix, i):
    return [row[i] for row in matrix]

cats = [ 'eeet', 'eemt', 'eett', 'eeem', 'mmet', 'mmmt', 'mmtt', 'mmem']
cat=str(sys.argv[2])

xsec=0.01
lumi=37.8

inFile = TFile.Open(sys.argv[1])
inFile.cd()


hW = inFile.Get("hWeights")
inTree = inFile.Get("Events")

nentries = inTree.GetEntries()

channel=str(sys.argv[2])

header=sys.argv[1]

if '/' in header : header=header.split("/",1)[1]
if '.root' in header : header=header.split(".",1)[0]

scale=str(sys.argv[3])

ScaleToLumi=False

if scale=="1" or scale=="true" or scale=="True" : 
    ScaleToLumi = True
    header=header+'_'+str(lumi)+'invfb'

if ScaleToLumi : print "Events scaled to {0:.2f}/fb for xsec = {1:.2f} pb".format(float(lumi),float(xsec))


rows, cols = (8, 9) 
arr = [[0 for i in range(cols)] for j in range(rows)] 
cuts=[]
result=[]

product=1.
if ScaleToLumi : product = xsec*lumi*1000/hW.GetSumOfWeights()
        
#print 'The weight will be', product, xsec, lumi, hW.GetSumOfWeights()

if 'all' not in channel :
    hist="hCutFlow_"+cat
    h1 = inFile.Get(hist)
    
    for i in range(1,h1.GetNbinsX()) :
        print h1.GetXaxis().GetBinLabel(i), h1.GetBinContent(i)*product

else:
    count=0
    for cat in cats:
    
        hist="hCutFlow_"+cat
       
        htest = inFile.Get(hist)
        for i in range(1,htest.GetNbinsX()) :
            try :
                arr[count][i-1] = '{0:.2f}'.format(htest.GetBinContent(i)*product)
                if cat == 'mmmt' : cuts.append(htest.GetXaxis().GetBinLabel(i))
            except IndexError :
                print("Index error in histogram loop:  i={0:d}".format(i))
        count+=1


with open(header+'_yields.txt', 'w') as f:

        print >> f,'\\documentclass[10pt]{report}'
        print >> f,'\\begin{document}'
	print >> f,'\\begin{table}[htp]'
        hh = header.replace('_','\\_')
	print >> f,'\\caption{' + '{0:s}'.format(hh) + '}'  
	print >> f,'\\begin{center}'
	print >> f,'\\begin{tabular}{l r r r r r r r r }  \hline'
        
        arrDict = {}
        topLine = 'Cut'
        for i,cat in enumerate(cats) :
            print("**cat={0:s}".format(cat))
            topLine += '& {0:s}'.format(cat)
            row = arr[i]
            print("cat={0:s} row={1:s}".format(cat,str(row)))
            for j,cut in enumerate(cuts) :
                arrDict[cut+cat] = row[j]
                
        print >> f, topLine + '\\\\ \\hline'
        for cut in cuts :
            line = '{0:s}'.format(cut)
            
            for cat in cats :
                val = float(arrDict[cut+cat])
                if val > 100. :
                    line += " & {0:9.0f}".format(val)
                else :
                    line += " & {0:9.2f}".format(val)
                    
            line += '\\\\ \\hline '
	    print >> f,line

	print >> f,'\\end{tabular}'
        print >> f,'\\end{center}'
	print >> f,'\\end{table}'
        print >> f,'\\end{document}'


	#np.savetxt('text.txt',arr, delimiter=' & ', fmt='%.2f', newline=' \n')
	#y = np.loadtxt('text.txt')
	#print y

