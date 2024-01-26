from ROOT import gSystem, gStyle, gROOT, kTRUE, gDirectory
from ROOT import TFile, TBranch, TTree
import os

def catToNumber(cat) :
    number = { 'eeet':1, 'eemt':2, 'eett':3, 'eeem':4, 'mmet':5, 'mmmt':6, 'mmtt':7, 'mmem':8, 'et':9, 'mt':10, 'tt':11 }
    return number[cat]

def numberToCat(number) :
    cat = { 1:'eeet', 2:'eemt', 3:'eett', 4:'eeem', 5:'mmet', 6:'mmmt', 7:'mmtt', 8:'mmem', 9:'et', 10:'mt', 11:'tt' }
    return cat[number]


def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--inFileName",default='in.root',help="Input file")
    parser.add_argument("-o","--outFileName",default='out.root',help="File to be used for output.")
    parser.add_argument("-c","--channel",default='mmmt',type=str,help="Channel.")

    return parser.parse_args()



args = getArgs()
Channel = args.channel
inFile = args.inFileName
outFile = args.outFileName 

fIn = TFile.Open(inFile)
inFile_name, inFile_extension = os.path.splitext(inFile)

fIn.cd()
treeOld=fIn.Get("Events")
treeOld.SetBranchStatus("*",1)


if Channel.lower() == 'all' :
    channels = ['eeem','eeet','eemt','eett','mmem','mmet','mmmt','mmtt']
else :
    channels = [Channel]

for channel in channels :
    outFile = "{0:s}_{1:s}{2:s}".format(inFile_name,channel,inFile_extension)
    fOut = TFile.Open(outFile,"recreate")
    newTree = treeOld.CloneTree(0)
    fOut.cd()

    nentries = treeOld.GetEntries()

    for iev, e in enumerate(treeOld) :
        #print iev, treeOld.cat, numberToCat(treeOld.cat)
        if numberToCat(treeOld.cat) == channel : newTree.Fill()

    newTree.AutoSave()

    fIn.cd()
    dirList2 = gDirectory.GetListOfKeys()
    for k2 in dirList2:
        h2 = k2.ReadObj()
        if h2.IsA ().InheritsFrom ("TH1") :
            htest=h2.Clone()
            newname=htest.GetName()
            if channel in newname :
                fOut.cd()
                htest.Write()

    print 'Succesfully made', outFile, ' for', channel, 'channel'
    fOut.Close()
    


