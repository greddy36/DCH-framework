#Get old file, old tree and set top branch address
from ROOT import TFile, TTree
oldfile = TFile.Open("AZH.root");
oldfile.cd()
inTree = oldfile.Get("Events")
nentries = inTree.GetEntries()

print("Read {0:d} events on input file".format(nentries))

#Create a new file + a clone of old tree in new file
newfile = TFile.Open("small.root","recreate");
newtree = inTree.CloneTree(0);

for e in inTree :
    if e.cat == 2 or e.cat == 5 : newtree.Fill()

newfile.Write()
newfile.Close()

exit()


#newtree->Print();
#newtree->AutoSave();
#delete oldfile;
#delete newfile;




 
