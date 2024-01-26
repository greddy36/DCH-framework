# read an nutuple and remove duplicate events

from ROOT import TFile, TTree
import sys

class dupeDetector() :
    
    def __init__(self):
        self.nCalls = 0
        self.runEventList = {}

    def checkEvent(self,entry) :
        self.nCalls += 1 
        runEvent = "{0:d}:{1:d}:{2:d}".format(entry.run,entry.evt,entry.cat)
        try :
            return self.runEventList[runEvent] 
        except KeyError :
            self.runEventList[runEvent] = True
            return False 

    def printSummary(self) :
        print("Duplicate Event Summary: Calls={0:d} Unique Events={1:d}".format(self.nCalls,len(self.runEventList.keys())))
        return

DD = dupeDetector()

inFileName = sys.argv[1]
inFile = TFile.Open(inFileName)
inFile.cd()
inTree = inFile.Get("Events")
nentries = inTree.GetEntries()

outFileName = sys.argv[2]
outFile = TFile(outFileName,'recreate')
outTree = inTree.CloneTree(0)
print("Number of entries in input tree = {0:d}".format(inTree.GetEntries()))

for i, e in enumerate(inTree) :
    if i % 1000 == 0 : print("i={0:d}".format(i))
    if not DD.checkEvent(e) : outTree.Fill() 

outFile.cd()
outFile.Write()
outFile.Close()

DD.printSummary()


