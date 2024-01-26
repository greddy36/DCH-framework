import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

class MyModule(Module):
    def __init__(self):
        self.counter = 0
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def analyze(self, event):
        # Accessing trigger branches
        #HLT_IsoMu24 = event.HLT_IsoMu24
        #HLT_IsoTkMu24 = event.HLT_IsoTkMu24
        HLT_IsoMu27 = event.HLT_IsoMu27
        muons = Collection(event, "Muon")

        #else :p=PostProcessor("./",fnames, "(   ( (HLT_IsoMu27) && Muon_pt > 28. && abs(Muon_eta) < 2.5 && abs(Muon_dxy) < 0.05 && abs(Muon_dz)<0.25 && Muon_pfRelIso04_all< 0.5  && (Muon_isGlobal || Muon_isTracker) && Muon_pfIsoId>1  ) ) && Entry$ >= Entries$*STARTEVENT &&  Entry$ < Entries$*FINISHEVENT","",[ jmeCorrections()],provenance=True,  outputbranchsel="keep_and_drop.txt",jsonInput=json_input)
        # Loop over muons in the event
        for muon in muons:

            # Check if the muon passes the selection criteria
            if (
                (HLT_IsoMu27 )
                and muon.pt > 28.
                and abs(muon.eta) < 2.5
                and abs(muon.dxy) < 0.05
                and abs(muon.dz) < 0.22
                and muon.pfRelIso04_all < 0.5
                and (muon.isGlobal or muon.isTracker)
                and muon.pfIsoId > 1
            ):
                self.counter += 1
                return True

        # If no muons pass the selection, discard the event
        return False

