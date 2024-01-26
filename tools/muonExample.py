## example how to read the muon format v2
## (Adapted from JMAR and EGM examples)
import sys
sys.path.insert(1,'../correctionlib/')
from correctionlib import _core

# Load CorrectionSet
fname = "../POG/MUO/2017_UL/muon_Z.json.gz"
fname = "./muon_Z_2017.json.gz"
if fname.endswith(".json.gz"):
    import gzip
    with gzip.open(fname,'rt') as file:
        data = file.read().strip()
        evaluator = _core.CorrectionSet.from_string(data)
else:
    evaluator = _core.CorrectionSet.from_file(fname)

# TrackerMuon Reconstruction UL scale factor
valsf = evaluator["NUM_MediumID_DEN_TrackerMuons"].evaluate("2017_UL", 1.1, 30.0, "sf")
print("sf 1 is: " + str(valsf))

valsf = evaluator["NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight"].evaluate("2017_UL", 1.1, 30.0, "sf")
print("sf IsoMu27 is: " + str(valsf))

# Medium ID UL scale factor, down variation
valsf = evaluator["NUM_MediumID_DEN_TrackerMuons"].evaluate("2017_UL", 1.1, 35.0, "systdown")
print("systdown is: " + str(valsf))

# Medium ID UL scale factor, up variation
valsf = evaluator["NUM_MediumID_DEN_TrackerMuons"].evaluate("2017_UL", 1.1, 35.0, "systup")
print("systup is: " + str(valsf))

# Trigger UL systematic uncertainty only 
#valsyst = evaluator["NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"].evaluate("2017_UL", 1.8, 54.0, "syst")
#print("syst is: " + str(valsyst))


