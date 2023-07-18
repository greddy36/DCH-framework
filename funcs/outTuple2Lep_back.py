# output ntuple for H->tautau analysis for CMSSW_10_2_X

from ROOT import TLorentzVector, TH1
from math import sqrt, sin, cos, pi, fabs
import tauFun2
import ROOT, array
import os
import sys
import generalFunctions as GF


electronMass = 0.0005
muonMass  = 0.105
class outTuple2Lep() :
    
    def __init__(self,fileName, era, doSyst=False,shift=[], isMC=True, onlyNom=False, isW=False):
        from array import array
        from ROOT import TFile, TTree

        # Tau Decay types
       
        ########### JetMet systematics
	#self.listsyst=['njets', 'nbtag', 'jpt', 'jeta', 'jflavour','MET_T1_pt', 'MET_T1_phi', 'MET_pt', 'MET_phi', 'MET_T1Smear_pt', 'MET_T1Smear_phi']
        self.jessyst=['_nom']
	self.listsyst=['MET_T1_pt', 'MET_T1_phi', 'MET_pt', 'MET_phi', 'PuppiMET_pt', 'PuppiMET_phi', 'MET_T1Smear_pt', 'MET_T1Smear_phi']
        if doSyst :
	    self.jessyst=['_nom', '_jesTotal', '_jer','_jesHEMIssue']  

        if onlyNom :
	    self.jessyst=['_nom']


	varss=['Up','Down']
        self.n = array('f', [ 0 ])

        self.allsystMET = []
        self.allsystJets = []
        self.jetsVariations=[]
        self.list_of_arrays = []           
        self.list_of_arrays_noES = []           
        self.list_of_arraysJetsPt = []           
        self.list_of_arraysJetsEta = []           
        self.list_of_arraysJetsFlavour = []           
        self.list_of_arraysJetsNbtagDeep = []           
        self.list_of_arraysJetsNbtagL = []           
        self.list_of_arraysJetsNbtagM = []           
        self.list_of_arraysJetsNbtagT = []           
        self.list_of_arraysJetsNjets = []           
        self.list_of_arraysJetsFlavour = []           
	self.tauMass = 1.7768 

        if not isMC  :
        
	    self.listsyst=[ 'MET_pt', 'MET_phi']
	    self.jessyst=['_nom']
	    varss=[]

        if doSyst : 

	    #self.jetsVariations.append('_nom')
	    self.allsystMET = []
	    self.allsystJets = []
	    #create a list with Up/Down from the above combinations
	    
	    for i_ in self.listsyst :
		for jes in self.jessyst :
		    if 'nom' not in jes :
			for var in varss :
			    if 'MET' in i_ and 'T1' in i_: 
				self.allsystMET.append(i_+jes+var)
				self.list_of_arrays.append(array('f', [ 0 ]))
				self.list_of_arrays_noES.append(array('f', [ 0 ]))
 
                    '''
		    if 'nom' in jes :
			if 'MET' in i_ : continue
			    #self.allsystMET.append(i_+jes)
			    #self.list_of_arrays.append(array('f', [ 0 ]))

                    ''' 
            for jes in self.jessyst :
		    if 'nom' in jes :   
			self.allsystJets.append(jes)
			self.list_of_arraysJetsNjets.append( array('f',[0]))
			self.list_of_arraysJetsNbtagL.append( array('f',[0]))
			self.list_of_arraysJetsNbtagM.append( array('f',[0]))
			self.list_of_arraysJetsNbtagT.append( array('f',[0]))
			self.list_of_arraysJetsFlavour.append( array('f',[-9.99]*12))
			self.list_of_arraysJetsEta.append( array('f',[-9.99]*12))
			self.list_of_arraysJetsPt.append( array('f',[-9.99]*12))
			self.list_of_arraysJetsNbtagDeep.append( array('f',[-9.99]*12))
		    else :   
		        for var in varss :
			    self.allsystJets.append(jes+var)
			    self.list_of_arraysJetsNjets.append( array('f',[0]))
			    self.list_of_arraysJetsNbtagL.append( array('f',[0]))
			    self.list_of_arraysJetsNbtagM.append( array('f',[0]))
			    self.list_of_arraysJetsNbtagT.append( array('f',[0]))
			    self.list_of_arraysJetsFlavour.append( array('f',[-9.99]*12))
			    self.list_of_arraysJetsEta.append( array('f',[-9.99]*12))
			    self.list_of_arraysJetsPt.append( array('f',[-9.99]*12))
			    self.list_of_arraysJetsNbtagDeep.append( array('f',[-9.99]*12))
                     
                
	    #for i_ in self.allsystMET :  self.list_of_arrays.append(array('f', [ 0 ]))

	    #for i_ in self.allsystJets :  
		
             

        print '------>systematics list', self.allsystMET
        print '------>jetssystematics list', self.allsystJets


        self.f = TFile( fileName, 'recreate' )
        self.t = TTree( 'Events', 'Output tree' )

        self.entries          = 0 
        self.run              = array('l',[0])
        self.nElectron        = array('l',[0])
        self.nMuon            = array('l',[0])
        self.nTau            = array('l',[0])
        self.lumi             = array('l',[0])
        self.evnt              = array('l',[0])
        self.nPU              = array('l',[0])
        self.nPUEOOT              = array('l',[0])
        self.nPULOOT              = array('l',[0])
        self.nPUtrue              = array('f',[0])
        self.nPV              = array('l',[0])
        self.nPVGood              = array('l',[0])
        self.nPVchi2              = array('f',[0])
        self.nPVndof              = array('f',[0])
        self.nPVscore              = array('f',[0])
        self.cat              = array('l',[0])
        self.weight           = array('f',[0])
        self.weightPU           = array('f',[0])
        self.weightPUtrue           = array('f',[0])
        self.LHEweight        = array('f',[0])
        self.Generator_weight = array('f',[0])
        self.LHE_Njets        = array('l',[0])
        self.electronTriggerWord  = array('l',[0])
        self.muonTriggerWord  = array('l',[0])         
        self.whichTriggerWord  = array('l',[0])         
        self.whichTriggerWordSubL  = array('l',[0])         
        self.LHEScaleWeights        = array('f',[1]*9)
        
        self.nGoodElectron    = array('l',[0])
        self.nGoodMuon        = array('l',[0])
        self.Flag_hfNoisyHitsFilter        = array('I',[0])
        self.Flag_BadPFMuonDzFilter        = array('I',[0])
        self.Electron_convVeto        = array('f',[0])
        self.Electron_lostHits        = array('I',[0])

        self.L1PreFiringWeight_Nom        = array('f',[0])
        self.L1PreFiringWeight_Up        = array('f',[0])
        self.L1PreFiringWeight_Down        = array('f',[0])

        # di-lepton variables.   1 and 2 refer to plus and minus charge
        self.mll       = array('f',[0])
        self.zll_pt       = array('f',[0])
        self.zll_phi       = array('f',[0])

        self.d0_1        = array('f',[0])
        self.dZ_1        = array('f',[0])
        
        self.pt_uncor_1        = array('f',[0])
        self.m_uncor_1        = array('f',[0])

        self.Electron_mvaFall17V2noIso_WP90_1 = array('f',[0])
        self.Electron_cutBased_1 = array('f',[0])
        self.gen_match_1 = array('l',[0])

        # di-lepton variables.   1 and 2 refer to plus and minus charge
        self.H_LT       = array('f',[0])
        self.dRlH       = array('f',[0])
        self.mll       = array('f',[0])
        self.W_Pt       = array('f',[0])
        self.pt_1      = array('f',[0])
        self.IDSF      = array('f',[0])
        self.IsoSF      = array('f',[0])
        self.TrigSF      = array('f',[0])
        #self.muonTightiDsf_1      = array('f',[0])


        self.m_1_tr   = array('f',[0])
        self.pt_1_tr   = array('f',[0])
        self.GenPart_statusFlags_1   = array('l',[0])
        self.GenPart_status_1     = array('l',[0])
        self.phi_1     = array('f',[0])
        self.phi_1_tr  = array('f',[0])
        self.eta_1     = array('f',[0])
        self.eta_1_tr  = array('f',[0])
        self.iso_1       = array('f',[0])
        self.PFiso_1       = array('f',[0])
        self.q_1       = array('f',[0])
        self.Muon_Id_1       = array('f',[0])
        self.isGlobal_1       = array('f',[0])
        self.isTracker_1       = array('f',[0])

        self.tightId_1       = array('f',[0])
        self.mediumId_1       = array('f',[0])
        self.mediumPromptId_1       = array('f',[0])
        self.looseId_1       = array('f',[0])

        # second candidate
        self.d0_2        = array('f',[0])
        self.dZ_2        = array('f',[0])
        
        self.pt_uncor_2        = array('f',[0])
        self.m_uncor_2        = array('f',[0])

        self.Electron_mvaFall17V2noIso_WP90_2 = array('f',[0])
        self.Electron_cutBased_2 = array('f',[0])
        self.gen_match_2 = array('l',[0])

        self.pt_2      = array('f',[0])


        self.m_2_tr   = array('f',[0])
        self.pt_2_tr   = array('f',[0])
        self.GenPart_statusFlags_2   = array('l',[0])
        self.GenPart_status_2     = array('l',[0])
        self.phi_2     = array('f',[0])
        self.phi_2_tr  = array('f',[0])
        self.eta_2     = array('f',[0])
        self.eta_2_tr  = array('f',[0])
        self.iso_2       = array('f',[0])
        self.q_2       = array('f',[0])
        self.Muon_Id_2       = array('f',[0])
        self.isGlobal_2       = array('f',[0])
        self.isTracker_2       = array('f',[0])
        self.tightId_2       = array('f',[0])
        self.mediumId_2       = array('f',[0])
        self.mediumPromptId_2       = array('f',[0])
        self.looseId_2       = array('f',[0])


        
        # MET variables
        self.metcov00    = array('f',[0])
        self.metcov01    = array('f',[0])
        self.metcov10    = array('f',[0])
        self.metcov11    = array('f',[0])

	self.MET_pt = array('f',[0])
	self.MET_phi = array('f',[0])
	self.MET_significance = array('f',[0])

	self.RawMET_pt = array('f',[0])
	self.RawMET_phi = array('f',[0])
	self.RawPuppiMET_pt = array('f',[0])
	self.RawPuppiMET_phi = array('f',[0])

	self.PuppiMET_pt = array('f',[0])
	self.PuppiMET_phi = array('f',[0])
	self.PuppiMET_ptJESUp = array('f',[0])
	self.PuppiMET_ptJESDown = array('f',[0])
	self.PuppiMET_ptJERUp = array('f',[0])
	self.PuppiMET_ptJERDown = array('f',[0])
	self.PuppiMET_pt_UnclusteredUp = array('f',[0])
	self.PuppiMET_pt_UnclusteredDown = array('f',[0])

	self.PuppiMET_phiJESUp = array('f',[0])
	self.PuppiMET_phiJESDown = array('f',[0])
	self.PuppiMET_phiJERUp = array('f',[0])
	self.PuppiMET_phiJERDown = array('f',[0])
	self.PuppiMET_phi_UnclusteredUp = array('f',[0])
	self.PuppiMET_phi_UnclusteredDown = array('f',[0])

        #does not exist in nAOV9, has to make it from DeltaX,Y
	self.MET_pt_UnclusteredUp = array('f',[0])
	self.MET_pt_UnclusteredDown = array('f',[0])
	self.MET_phi_UnclusteredUp = array('f',[0])
	self.MET_phi_UnclusteredDown = array('f',[0])

	self.MET_ptJESUp = array('f',[0])
	self.MET_ptJESDown = array('f',[0])
	self.MET_ptJERUp = array('f',[0])
	self.MET_ptJERDown = array('f',[0])

	self.MET_phiJESUp = array('f',[0])
	self.MET_phiJESDown = array('f',[0])
	self.MET_phiJERUp = array('f',[0])
	self.MET_phiJERDown = array('f',[0])

	self.MET_T1_pt = array('f',[0])
	self.MET_T1_phi = array('f',[0])

	self.MET_T1_pt_UnclusteredUp = array('f',[0])
	self.MET_T1_pt_UnclusteredDown = array('f',[0])
	self.MET_T1_phi_UnclusteredUp = array('f',[0])
	self.MET_T1_phi_UnclusteredDown = array('f',[0])


        self.isTrig_1   = array('f',[0])
        self.isTrigObj   = array('l',[0])

        # jet variables
        #self.njetsold = array('f',[-1]*12)
        self.njets     = array('f',[0])
        self.nbtagL     = array('f',[0])
        self.nbtagM     = array('f',[0])
        self.nbtagT     = array('f',[0])

        self.jflavour     = array('f',[-9.99]*12)
        self.jeta     = array('f',[-9.99]*12)
        self.jpt     = array('f',[-9.99]*12)
        self.btagDeep     = array('f',[-9.99]*12)

        self.bpt_1     = array('f',[0]*12)
        self.bpt_1_tr  = array('f',[0]*12)
        self.beta_1    = array('f',[0]*12)
        self.beta_1_tr = array('f',[0]*12)
        self.bphi_1    = array('f',[0]*12)
        self.bphi_1_tr = array('f',[0]*12)
        self.bcsv_1    = array('f',[0]*12)
        self.bcsvfv_1    = array('f',[0]*12)
      
        self.t.Branch('run',              self.run,               'run/l' )
        self.t.Branch('nElectron',              self.nElectron,               'nElectron/l' )
        self.t.Branch('nMuon',              self.nMuon,               'nMuon/l' )
        self.t.Branch('nTau',              self.nTau,               'nTau/l' )
        self.t.Branch('lumi',             self.lumi,              'lumi/I' )
        self.t.Branch('evnt',              self.evnt,               'evnt/l' )
        self.t.Branch('nPU',              self.nPU,               'nPU/I' )
        self.t.Branch('nPUEOOT',              self.nPUEOOT,               'nPUEOOT/I' )
        self.t.Branch('nPULOOT',              self.nPULOOT,               'nPULOOT/I' )
        self.t.Branch('nPUtrue',              self.nPUtrue,               'nPUtrue/F' )
        self.t.Branch('nPV',              self.nPV,               'nPV/I' )
        self.t.Branch('nPVGood',              self.nPVGood,               'nPVGood/I' )
        self.t.Branch('cat',              self.cat,               'cat/I' )
        self.t.Branch('weight',           self.weight,            'weight/F' )
        self.t.Branch('weightPU',           self.weightPU,            'weightPU/F' )
        self.t.Branch('weightPUtrue',           self.weightPUtrue,            'weightPUtrue/F' )
        self.t.Branch('LHEweight',        self.LHEweight,         'LHEweight/F' )
        self.t.Branch('LHE_Njets',        self.LHE_Njets,         'LHE_Njets/I' )
        self.t.Branch('LHEScaleWeights',        self.LHEScaleWeights,         'LHEScaleWeights[9]/F' )
        self.t.Branch('Generator_weight', self.Generator_weight,  'Generator_weight/F' )
        self.t.Branch('electronTriggerWord',  self.electronTriggerWord, 'electronTriggerWord/I' )
        self.t.Branch('muonTriggerWord',      self.muonTriggerWord,  'muonTriggerWord/I' )
        self.t.Branch('whichTriggerWord',      self.whichTriggerWord,  'whichTriggerWord/I' )
        self.t.Branch('whichTriggerWordSubL',      self.whichTriggerWordSubL,  'whichTriggerWordSubL/I' )
        
        self.t.Branch('nGoodElectron',    self.nGoodElectron,     'nGoodElectron/I' )
        self.t.Branch('nGoodMuon',        self.nGoodMuon,         'nGoodMuon/I' )
        self.t.Branch('mll',        self.mll,         'mll/F' )
        self.t.Branch('zll_pt',        self.zll_pt,         'zll_pt/F' )
        self.t.Branch('zll_phi',        self.zll_phi,         'zll_phi/F' )
        
        self.t.Branch('GenPart_statusFlags_1',     self.GenPart_statusFlags_1,     'GenPart_statusFlags_1/I')
        self.t.Branch('GenPart_status_1',     self.GenPart_status_1,     'GenPart_status_1/I')
        self.t.Branch('pt_uncor_1',        self.pt_uncor_1,        'pt_uncor_1/F')
        self.t.Branch('m_uncor_1',        self.m_uncor_1,        'm_uncor_1/F')
        self.t.Branch('gen_match_1', self.gen_match_1, 'gen_match_1/l')


        self.t.Branch('pt_1',        self.pt_1,        'pt_1/F')
        self.t.Branch('m_1_tr',     self.m_1_tr,     'm_1_tr/F')
        self.t.Branch('pt_1_tr',     self.pt_1_tr,     'pt_1_tr/F')
        self.t.Branch('phi_1',       self.phi_1,       'phi_1/F')  
        self.t.Branch('phi_1_tr',    self.phi_1_tr,    'phi_1_tr/F')
        self.t.Branch('eta_1',       self.eta_1,       'eta_1/F')    
        self.t.Branch('iso_1',       self.iso_1,       'iso_1/F')
        self.t.Branch('q_1',       self.q_1,       'q_1/F')
        self.t.Branch('L1PreFiringWeight_Nom',        self.L1PreFiringWeight_Nom,        'L1PreFiringWeight_Nom/F')
        self.t.Branch('L1PreFiringWeight_Up',        self.L1PreFiringWeight_Up,        'L1PreFiringWeight_Up/F')
        self.t.Branch('L1PreFiringWeight_Down',        self.L1PreFiringWeight_Down,        'L1PreFiringWeight_Down/F')
        self.t.Branch('d0_1',        self.d0_1,        'd0_1/F')
        self.t.Branch('dZ_1',        self.dZ_1,        'dZ_1/F')
        self.t.Branch('Muon_Id_1',       self.Muon_Id_1,       'Muon_Id_1/F')
        self.t.Branch('isGlobal_1',       self.isGlobal_1,       'isGlobal_1/F')
        self.t.Branch('isTracker_1',       self.isTracker_1,       'isTracker_1/F')
        self.t.Branch('tightId_1', self.tightId_1, 'tightId_1/F')
        self.t.Branch('mediumId_1', self.mediumId_1, 'mediumId_1/F')
        self.t.Branch('mediumPromptId_1', self.mediumPromptId_1, 'mediumPromptId_1/F')
        self.t.Branch('looseId_1', self.looseId_1, 'looseId_1/F')


        self.t.Branch('GenPart_statusFlags_2',     self.GenPart_statusFlags_2,     'GenPart_statusFlags_2/I')
        self.t.Branch('GenPart_status_2',     self.GenPart_status_2,     'GenPart_status_2/I')
        self.t.Branch('pt_uncor_2',        self.pt_uncor_2,        'pt_uncor_2/F')
        self.t.Branch('m_uncor_2',        self.m_uncor_2,        'm_uncor_2/F')
        self.t.Branch('gen_match_2', self.gen_match_2, 'gen_match_2/l')


        self.t.Branch('pt_2',        self.pt_2,        'pt_2/F')
        self.t.Branch('m_2_tr',     self.m_2_tr,     'm_2_tr/F')
        self.t.Branch('pt_2_tr',     self.pt_2_tr,     'pt_2_tr/F')
        self.t.Branch('phi_2',       self.phi_2,       'phi_2/F')  
        self.t.Branch('phi_2_tr',    self.phi_2_tr,    'phi_2_tr/F')
        self.t.Branch('eta_2',       self.eta_2,       'eta_2/F')    
        self.t.Branch('iso_2',       self.iso_2,       'iso_2/F')
        self.t.Branch('q_2',       self.q_2,       'q_2/F')
        self.t.Branch('L1PreFiringWeight_Nom',        self.L1PreFiringWeight_Nom,        'L1PreFiringWeight_Nom/F')
        self.t.Branch('L1PreFiringWeight_Up',        self.L1PreFiringWeight_Up,        'L1PreFiringWeight_Up/F')
        self.t.Branch('L1PreFiringWeight_Down',        self.L1PreFiringWeight_Down,        'L1PreFiringWeight_Down/F')
        self.t.Branch('d0_2',        self.d0_2,        'd0_2/F')
        self.t.Branch('dZ_2',        self.dZ_2,        'dZ_2/F')
        self.t.Branch('Muon_Id_2',       self.Muon_Id_2,       'Muon_Id_2/F')
        self.t.Branch('isGlobal_2',       self.isGlobal_2,       'isGlobal_2/F')
        self.t.Branch('isTracker_2',       self.isTracker_2,       'isTracker_2/F')
        self.t.Branch('tightId_2', self.tightId_2, 'tightId_2/F')
        self.t.Branch('mediumId_2', self.mediumId_2, 'mediumId_2/F')
        self.t.Branch('mediumPromptId_2', self.mediumPromptId_2, 'mediumPromptId_2/F')
        self.t.Branch('looseId_2', self.looseId_2, 'looseId_2/F')


        
        # MET variables
        self.t.Branch('metcov00', self.metcov00, 'metcov00/F')
        self.t.Branch('metcov01', self.metcov01, 'metcov01/F')
        self.t.Branch('metcov10', self.metcov10, 'metcov10/F')
        self.t.Branch('metcov11', self.metcov11, 'metcov11/F')




        self.t.Branch('RawMET_pt', self.RawMET_pt, 'RawMET_pt /F')


        self.t.Branch('RawMET_pt', self.RawMET_pt, 'RawMET_pt /F')
        self.t.Branch('RawMET_phi', self.RawMET_phi, 'RawMET_phi /F')
        self.t.Branch('RawPuppiMET_pt', self.RawPuppiMET_pt, 'RawPuppiMET_pt /F')
        self.t.Branch('RawPuppiMET_phi', self.RawPuppiMET_phi, 'RawPuppiMET_phi /F')

        self.t.Branch('PuppiMET_pt', self.PuppiMET_pt, 'PuppiMET_pt /F')
        self.t.Branch('PuppiMET_phi', self.PuppiMET_phi, 'PuppiMET_phi /F')

        self.t.Branch('PuppiMET_ptJESUp', self.PuppiMET_ptJESUp, 'PuppiMET_ptJESUp /F')
        self.t.Branch('PuppiMET_ptJESDown', self.PuppiMET_ptJESDown, 'PuppiMET_ptJESDown /F')
        self.t.Branch('PuppiMET_ptJERUp', self.PuppiMET_ptJERUp, 'PuppiMET_ptJERUp /F')
        self.t.Branch('PuppiMET_ptJERDown', self.PuppiMET_ptJERDown, 'PuppiMET_ptJERDown /F')
        self.t.Branch('PuppiMET_pt_UnclusteredUp', self.PuppiMET_pt_UnclusteredUp, 'PuppiMET_pt_UnclusteredUp /F')
        self.t.Branch('PuppiMET_pt_UnclusteredDown', self.PuppiMET_pt_UnclusteredDown, 'PuppiMET_pt_UnclusteredDown /F')

        self.t.Branch('PuppiMET_phiJESUp', self.PuppiMET_phiJESUp, 'PuppiMET_phiJESUp /F')
        self.t.Branch('PuppiMET_phiJESDown', self.PuppiMET_phiJESDown, 'PuppiMET_phiJESDown /F')
        self.t.Branch('PuppiMET_phiJERUp', self.PuppiMET_phiJERUp, 'PuppiMET_phiJERUp /F')
        self.t.Branch('PuppiMET_phiJERDown', self.PuppiMET_phiJERDown, 'PuppiMET_phiJERDown /F')
        self.t.Branch('PuppiMET_phi_UnclusteredUp', self.PuppiMET_phi_UnclusteredUp, 'PuppiMET_phi_UnclusteredUp /F')
        self.t.Branch('PuppiMET_phi_UnclusteredDown', self.PuppiMET_phi_UnclusteredDown, 'PuppiMET_phi_UnclusteredDown /F')

        self.t.Branch('MET_ptJESUp', self.MET_ptJESUp, 'MET_ptJESUp /F')
        self.t.Branch('MET_ptJESDown', self.MET_ptJESDown, 'MET_ptJESDown /F')
        self.t.Branch('MET_ptJERUp', self.MET_ptJERUp, 'MET_ptJERUp /F')
        self.t.Branch('MET_ptJERDown', self.MET_ptJERDown, 'MET_ptJERDown /F')
        self.t.Branch('MET_pt_UnclusteredUp', self.MET_pt_UnclusteredUp, 'MET_pt_UnclusteredUp /F')
        self.t.Branch('MET_pt_UnclusteredDown', self.MET_pt_UnclusteredDown, 'MET_pt_UnclusteredDown /F')

        self.t.Branch('MET_phiJESUp', self.MET_phiJESUp, 'MET_phiJESUp /F')
        self.t.Branch('MET_phiJESDown', self.MET_phiJESDown, 'MET_phiJESDown /F')
        self.t.Branch('MET_phiJERUp', self.MET_phiJERUp, 'MET_phiJERUp /F')
        self.t.Branch('MET_phiJERDown', self.MET_phiJERDown, 'MET_phiJERDown /F')
        self.t.Branch('MET_phi_UnclusteredUp', self.MET_phi_UnclusteredUp, 'MET_phi_UnclusteredUp /F')
        self.t.Branch('MET_phi_UnclusteredDown', self.MET_phi_UnclusteredDown, 'MET_phi_UnclusteredDown /F')


        self.t.Branch('MET_phi', self.MET_phi, 'MET_phi /F')
        self.t.Branch('MET_pt', self.MET_pt, 'MET_pt /F')

        self.t.Branch('MET_significance', self.MET_significance, 'MET_significance /F')
        self.t.Branch('MET_T1_pt', self.MET_T1_pt, 'MET_T1_pt /F')
        self.t.Branch('MET_T1_phi', self.MET_T1_phi, 'MET_T1_phi /F')

        self.t.Branch('MET_T1_pt_UnclusteredUp', self.MET_T1_pt_UnclusteredUp, 'MET_T1_pt_UnclusteredUp /F')
        self.t.Branch('MET_T1_pt_UnclusteredDown', self.MET_T1_pt_UnclusteredDown, 'MET_T1_pt_UnclusteredDown /F')

        self.t.Branch('MET_T1_phi_UnclusteredUp', self.MET_T1_phi_UnclusteredUp, 'MET_T1_phi_UnclusteredUp /F')
        self.t.Branch('MET_T1_phi_UnclusteredDown', self.MET_T1_phi_UnclusteredDown, 'MET_T1_phi_UnclusteredDown /F')




        # trigger sf
        self.t.Branch('isTrig_1',  self.isTrig_1, 'isTrig_1/F' )
        self.t.Branch('isTrigObj',  self.isTrigObj, 'isTrigObj/l' )


        # jet variables
        #self.t.Branch('njetsold', self.njetsold, 'njetsold[8]/F') 
        #self.t.Branch('nbtagold', self.nbtagold, 'nbtagold[8]/F')
        self.t.Branch('njets', self.njets, 'njets/F')
        self.t.Branch('nbtagL', self.nbtagL, 'nbtagL/F')
        self.t.Branch('nbtagM', self.nbtagM, 'nbtagM/F')
        self.t.Branch('nbtagT', self.nbtagT, 'nbtagT/F')

        self.t.Branch('jflavour',     self.jflavour,     'jflavour[12]/F' )
        self.t.Branch('jeta',     self.jeta,     'jeta[12]/F' )
        self.t.Branch('jpt',     self.jpt,     'jpt[12]/F' )
        self.t.Branch('btagDeep', self.btagDeep, 'btagDeep[12]/F')

        if doSyst : 
                #Book the branches and the arrays needed to store variables
		for i, v in enumerate(self.allsystMET):
                 
                    if str(era)=='2017' : 
                        v = v.replace('MET','METFixEE2017')
                    iMET= v.replace('METFixEE2017','MET')
                    iiMET=iMET+'_noES'
	            self.t.Branch(iMET, self.list_of_arrays[i], '{0:s}/F'.format(iMET))
	            self.t.Branch(iiMET, self.list_of_arrays_noES[i], '{0:s}/F'.format(iiMET))

		for i, v in enumerate(self.allsystJets):
		    self.t.Branch('njets{0:s}'.format(v), self.list_of_arraysJetsNjets[i], 'njets{0:s}/F'.format(v))
		    self.t.Branch('nbtagL{0:s}'.format(v), self.list_of_arraysJetsNbtagL[i], 'nbtagL{0:s}/F'.format(v))
		    self.t.Branch('nbtagM{0:s}'.format(v), self.list_of_arraysJetsNbtagM[i], 'nbtagM{0:s}/F'.format(v))
		    self.t.Branch('nbtagT{0:s}'.format(v), self.list_of_arraysJetsNbtagT[i], 'nbtagT{0:s}/F'.format(v))
		    self.t.Branch('jflavour{0:s}'.format(v), self.list_of_arraysJetsFlavour[i], 'jflavour{0:s}[12]/F'.format(v))
		    self.t.Branch('jpt{0:s}'.format(v), self.list_of_arraysJetsPt[i], 'jpt{0:s}[12]/F'.format(v))
		    self.t.Branch('jeta{0:s}'.format(v), self.list_of_arraysJetsEta[i], 'jeta{0:s}[12]/F'.format(v))
		    self.t.Branch('btagDeep{0:s}'.format(v), self.list_of_arraysJetsNbtagDeep[i], 'bagDeep{0:s}[12]/F'.format(v))



        #self.MET_pt_jesEC2Up  = array('f',[0])
        #self.t.Branch('MET_pt_jesEC2Up', self.MET_pt_jesEC2Up, 'MET_pt_jesEC2Up/F' )
        self.tN=[]

	#self.t.SetBranchStatus("*Up",0)
	#self.t.SetBranchStatus("*Down",0)
	self.t.SetBranchStatus("GenPart*",0)
	self.t.SetBranchStatus("*_tr*",0)
	self.t.SetBranchStatus("*LHE*",1)
	#self.t.SetBranchStatus("*LHEScaleWeight",1)
	self.t.SetBranchStatus("*Up*",0)
	self.t.SetBranchStatus("*Down*",0)
	#self.t.SetBranchStatus("Smear",0)
        for i, isyst in enumerate(shift) : 
	    self.tN.append(isyst)

            #if isyst == "Events" : continue
            #else  : 
            if i > 0 : 
                self.tN[i-1]  = self.t.CloneTree()
                #self.t.SetBranchStatus("Smear",1)
                self.tN[i-1].SetName(isyst)

                print '====================>',self.tN[i-1], self.tN[i-1].GetName()

	self.t.SetBranchStatus("GenPart*",1)
	self.t.SetBranchStatus("*_tr*",1)
	#self.t.SetBranchStatus("*LHE*",1)
	self.t.SetBranchStatus("*LHEScaleWeight*",1)
	self.t.SetBranchStatus("*Up*",1)
	self.t.SetBranchStatus("*Down*",1)

    def get_mt(self,METtype,entry,tau) :
        if METtype == 'MVAMet' :
            # temporary choice 
            dphi = tau.Phi() - entry.MET_phi
            return sqrt(2.*tau.Pt()*entry.MET_pt*(1. - cos(dphi)))
        elif METtype == 'PFMet' :
            dphi = tau.Phi() - entry.MET_phi
            return sqrt(2.*tau.Pt()*entry.MET_pt*(1. - cos(dphi)))
        elif METtype == 'PUPPIMet' :
            dphi = tau.Phi() - entry.PuppiMET_phi
            return sqrt(2.*tau.Pt()*entry.PuppiMET_pt*(1. - cos(dphi)))
        else :
            print("Invalid METtype={0:s} in outTuple.get_mt().   Exiting".format(METtype))

    def getPt_tt(self,entry,tau1,tau2) :
        ptMiss = TLorentzVector() 
        ptMiss.SetPtEtaPhiM(entry.MET_pt,0.,entry.MET_phi,0.)
        return (tau1+tau2+ptMiss).Pt()

    def getMt_tot(self,entry,tau1,tau2) :
        pt1, pt2, met = tau1.Pt(), tau2.Pt(), entry.MET_pt
        phi1, phi2, metphi = tau1.Phi(), tau2.Phi(), entry.MET_phi
        arg = 2.*(pt1*met*(1. - cos(phi1-metphi)) + pt2*met*(1. - cos(phi2-metphi)) + pt1*pt2*(1. - cos(phi2-phi1)))
        return sqrt(arg)

    def getDR(self,entry, v1,v2) :

        dPhi = min(abs(v2.Phi()-v1.Phi()),2.*pi-abs(v2.Phi()-v1.Phi()))
        DR = sqrt(dPhi**2 + (v2.Eta()-v1.Eta())**2)
	return DR

    def getDRnV(self,entry, eta1,phi1, eta2,phi2) :

        dPhi = min(abs(phi2-phi1),2.*pi-abs(phi2-phi1))
        DR = sqrt(dPhi**2 + (eta2-eta1)**2)
	return DR

    def getdPhi(self, entry, v1,v2) :
        dPhi = min(abs(v2.Phi()-v1.Phi()),2.*pi-abs(v2.Phi()-v1.Phi()))
        return dPhi

    def getM_vis(self,entry,tau1,tau2) :
        return (tau1+tau2).M()

    def getJets(self,entry,tau1,tau2,era) :
	nJet30, jetList, bJetList, bJetListFlav = 0, [], [], []
        phi2_1, eta2_1 = tau1.Phi(), tau1.Eta() 
        phi2_2, eta2_2 = tau2.Phi(), tau2.Eta() 
	bjet_discr = 0.6321
	bjet_discrFlav = 0.0614
	if str(era) == '2017' : bjet_discr = 0.4941
	if str(era) == '2018' : bjet_discr = 0.4184

        for j in range(entry.nJet) :
            if entry.Jet_jetId[j]  < 2  : continue  #require tight jets
            if entry.Jet_pt[j]>20 and entry.Jet_pt[j] < 50 and entry.Jet_puId[j]  < 4  : continue #loose jetPU_iD
            if str(era) == '2017'  and entry.Jet_pt[j] > 20 and entry.Jet_pt[j] < 50 and abs(entry.Jet_eta[j]) > 2.65 and abs(entry.Jet_eta[j]) < 3.139 : continue  #remove noisy jets
            if entry.Jet_pt[j] < 20. : continue
            if abs(entry.Jet_eta[j]) > 4.7 : continue
            phi1, eta1 = entry.Jet_phi[j], entry.Jet_eta[j]
            dPhi = min(abs(phi2_1-phi1),2.*pi-abs(phi2_1-phi1))
            DR = sqrt(dPhi**2 + (eta2_1-eta1)**2)
            dPhi = min(abs(phi2_2-phi1),2.*pi-abs(phi2_2-phi1))
            DR = min(DR,sqrt(dPhi**2 + (eta2_2-eta1)**2))
            if DR < 0.5 : continue
            if entry.Jet_pt[j] > 30 :
		if abs(entry.Jet_eta[j]) < 2.4 and entry.Jet_btagDeepB[j] > bjet_discr : bJetList.append(j)
		if abs(entry.Jet_eta[j]) < 2.4 and entry.Jet_btagDeepFlavB[j] > bjet_discrFlav : bJetListFlav.append(j)
                jetList.append(j) 

        return jetList, bJetList,bJetListFlav



    def getJetsJMEMV(self,entry,LepList,era, syst) :
	jetList, jetListFlav, jetListEta, jetListPt, bTagListDeep, bJetListL, bJetListM, bJetListT, bJetListFlav = [], [], [], [], [], [], [], [], []
	#print 'will try', len(LepList)
	bjet_discrL = 0.2217
	bjet_discrM = 0.6321
	bjet_discrT = 0.8953
	bjet_discrFlav = 0.0614

	if str(era) == '2017' : 
	    bjet_discrL = 0.1522
	    bjet_discrM = 0.4941
	    bjet_discrT = 0.8001
	if str(era) == '2018' : 
	    bjet_discrL = 0.1241
	    bjet_discrM = 0.4184
	    bjet_discrT = 0.7527

	failJets=[]
        goodJets=[]
        bJetListL=[]
        bJetListM=[]
        bJetListT=[]
        bTagListDeep=[]
        #if syst !='' : syst="_"+syst
     
        if 'nom' in syst : syst='_nom'
          
        #if entry.event==18093 and syst=='_jesEC2Up' : print "Jet_pt{0:s}".format(str(syst)), entry.event, entry.nJet
        for j in range(entry.nJet) :

            try : 
		jpt = getattr(entry, "Jet_pt{0:s}".format(str(syst)), None)
                #if syst=='_nom' : print jpt[j],  entry.Jet_pt[j],  syst
                #if entry.event==18093 and syst=='_jesEC2Up' : print 'inside jets', jpt[j], syst, entry.event, "Jet_pt{0:s}".format(str(syst))

		if entry.Jet_jetId[j]  < 2  : continue  #require tight jets
		if jpt[j] > 30 and jpt[j] < 50 and entry.Jet_puId[j]  < 4  : continue #loose jetPU_iD
		if str(era) == '2017'  and jpt[j] > 20 and jpt[j] < 50 and abs(entry.Jet_eta[j]) > 2.65 and abs(entry.Jet_eta[j]) < 3.139 : continue  #remove noisy jets
		if jpt[j] < 25. : continue
		if abs(entry.Jet_eta[j]) > 4.7 : continue

		#for iv, lepv in enumerate(LepList) : 
		for iv, lv  in  enumerate(LepList) :
		    dr = self.getDRnV(entry, entry.Jet_eta[j], entry.Jet_phi[j], LepList[iv].Eta(), LepList[iv].Phi())
		    if float(dr) > 0.5 : 
			#print 'seems goodfor iv--->', iv, 'jet', j, entry.nJet, 'dr--', dr , LepList[iv].Eta(), LepList[iv].Phi(), LepList[iv].Pt()
			if j not in goodJets : goodJets.append(j)
		    if float(dr) < 0.5 : 
			#print ' failed for lepton--->', iv, 'jet', j, 'njets', entry.nJet, 'dr--', dr , LepList[iv].Eta(), LepList[iv].Phi(), LepList[iv].Pt()
			if j not in failJets : failJets.append(j)
			#continue
            except : continue

        #print 'will check failed jets',  entry.luminosityBlock, entry.event, entry.run, failJets, goodJets, 'from nJet i', j, entry.nJet
        for j in failJets : 
            if j in goodJets : goodJets.remove(j)


        for jj in goodJets : 
            #if isMC : 
            try : 
                jetListFlav.append(entry.Jet_partonFlavour[jj])
            except AttributeError  : jetListFlav.append(0)
            jetListEta.append(entry.Jet_eta[jj])
            jpt = getattr(entry, "Jet_pt{0:s}".format(str(syst)), None)
            jetListPt.append(jpt[jj])
            bTagListDeep.append(entry.Jet_btagDeepB[jj])

            #print 'will check',  entry.luminosityBlock, entry.event, entry.run, goodJets, jj, jpt[jj], 'flav', entry.Jet_partonFlavour[jj]
            if jpt[jj] > 25 : 
                
		if abs(entry.Jet_eta[jj]) < 2.4 : 
		    if entry.Jet_btagDeepB[jj] > bjet_discrL : bJetListL.append(jj)
		    if entry.Jet_btagDeepB[jj] > bjet_discrM : bJetListM.append(jj)
		    if entry.Jet_btagDeepB[jj] > bjet_discrT : bJetListT.append(jj)
		    if entry.Jet_btagDeepFlavB[jj] > bjet_discrFlav : bJetListFlav.append(jj)
            if jpt[jj] > 30 : 
                jetList.append(jj) 
                #print '--added ', jj, 'in good list', jpt[jj], abs(entry.Jet_eta[jj])

        #if entry.event==18093 and syst=='_jesEC2Up': print 'going out....', jetList, jetListPt, syst
        #if len(jetList)!=len(jetListPt) : print 'going out....', jetList, jetListPt, syst, len(jetList), len(jetListPt), entry.luminosityBlock, entry.event, entry.run
        #print 'going out....', jetList, jetListPt, syst, len(jetList), len(jetListPt), entry.luminosityBlock, entry.event, entry.run, btagWeightDeepCSVB
        #print ''
        return jetList, jetListFlav, jetListEta,  jetListPt, bTagListDeep, bJetListL,bJetListM,bJetListT,bJetListFlav



    def runSVFit(self, entry, channel, jt1, jt2, tau1, tau2, metpt, metphi) :
                      
        measuredMETx = metpt*cos(metphi)
        measuredMETy = metpt*sin(metphi)

        #define MET covariance
        covMET = ROOT.TMatrixD(2,2)
        covMET[0][0] = entry.MET_covXX
        covMET[1][0] = entry.MET_covXY
        covMET[0][1] = entry.MET_covXY
        covMET[1][1] = entry.MET_covYY
        #covMET[0][0] = 787.352
        #covMET[1][0] = -178.63
        #covMET[0][1] = -178.63
        #covMET[1][1] = 179.545

        #self.kUndefinedDecayType, self.kTauToHadDecay,  self.kTauToElecDecay, self.kTauToMuDecay = 0, 1, 2, 3

        if channel == 'et' :
            measTau1 = ROOT.MeasuredTauLepton(self.kTauToElecDecay, tau1.Pt(), tau1.Eta(), tau1.Phi(), 0.000511) 
        elif channel == 'mt' :
            measTau1 = ROOT.MeasuredTauLepton(self.kTauToMuDecay, tau1.Pt(), tau1.Eta(), tau1.Phi(), 0.106) 
        elif channel == 'tt' :
            measTau1 = ROOT.MeasuredTauLepton(self.kTauToHadDecay, tau1.Pt(), tau1.Eta(), tau1.Phi(), entry.Tau_mass[jt1])
                        
	if channel != 'em' :
            measTau2 = ROOT.MeasuredTauLepton(self.kTauToHadDecay, tau2.Pt(), tau2.Eta(), tau2.Phi(), entry.Tau_mass[jt2])

	if channel == 'em' :
            measTau1 = ROOT.MeasuredTauLepton(self.kTauToElecDecay, tau1.Pt(), tau1.Eta(), tau1.Phi(), 0.000511)
            measTau2 = ROOT.MeasuredTauLepton(self.kTauToMuDecay, tau2.Pt(), tau2.Eta(), tau2.Phi(), 0.106)

        VectorOfTaus = ROOT.std.vector('MeasuredTauLepton')
        instance = VectorOfTaus()
        instance.push_back(measTau1)
        instance.push_back(measTau2)

        FMTT = ROOT.FastMTT()
        FMTT.run(instance, measuredMETx, measuredMETy, covMET)
        ttP4 = FMTT.getBestP4()
        return ttP4.M(), ttP4.Mt() 
    

    def Fill2Lep(self, entry, cat, LepP, LepM, lepList, isMC, era, doUncertainties=False , proc="EOY") : 
    #def Fill(self, entry, SVFit, cat, jt1, jt2, LepP, LepM, lepList, isMC, era, doUncertainties=False ,  met_pt=-99, met_phi=-99, systIndex=0) : 
        SystIndex = 0

        #if SystIndex >0 : doUncertainties=False

        #channel_ll = 'mm' or 'ee'
        channel_ll = cat
	channel = cat

        if SystIndex ==0 : 

	    is_trig_1, is_trig_2, is_Dtrig_1 = 0., 0., 0.
	    TrigListLep = []
	    TrigListTau = []
	    hltListLep  = []
	    hltListLepSubL  = []

	    TrigListLep, hltListLep, hltListLepSubL  = GF.findSingleLeptTrigger(lepList, entry, channel_ll, era)

	    TrigListLep = list(dict.fromkeys(TrigListLep))
	    #if len(hltListLep) > 0 or len(hltListLepSubL)>0 :     print GF.printEvent(entry), SystIndex

	    if len(hltListLep) > 0 and  len(hltListLepSubL) == 0 :
		is_trig_1 = 1
	    if len(hltListLep) == 0 and len(hltListLepSubL) > 0 :
		is_trig_1 = -1
	    if len(hltListLep) > 0 and len(hltListLepSubL)>0 :
		is_trig_1 = 2

	    self.whichTriggerWord[0]=0
	    self.whichTriggerWordSubL[0]=0

	    #if len(TrigListLep) >0 : print 'TrigerList ===========>', TrigListLep, lepList, hltListLep, channel_ll, 'istrig_1', is_trig_1, 'istrig_2', is_trig_2, 'lenTrigList', len(TrigListLep),  'lenLept', len(lepList), 'lepList_0', lepList[0], 'TrigList_0', TrigListLep[0], hltListLep
	    
	    for i,bit in enumerate(hltListLep):
		    
		if bit : 
		    self.whichTriggerWord[0] += 2**i

	    for j,bitt in enumerate(hltListLepSubL):
		if bitt : self.whichTriggerWordSubL[0] += 2**j

	    #if channel_ll=='ee' and entry.luminosityBlock==90 and entry.event==8904: print self.whichTriggerWord[0], 'hlt', hltListLep, 'hltsub', hltListLepSubL
	    #print cat, self.whichTriggerWord
	    # channel = 'mt', 'et', 'tt', or 'em'
	    self.entries += 1

	    self.run[0]  = entry.run
	    self.nElectron[0]  = entry.nElectron
	    self.nMuon[0]  = entry.nMuon
	    self.nTau[0]  = entry.nTau
	    self.lumi[0] = entry.luminosityBlock 
	    self.evnt[0]  = entry.event
	    self.mll[0]  = (LepP+LepM).M()
	    self.zll_pt[0]  = (LepP+LepM).Pt()
	    self.zll_phi[0]  = (LepP+LepM).Phi()

	    self.iso_1[0]  = -99
	    self.q_1[0]  = -99
	    self.isGlobal_1[0]  = -99
	    self.iso_2[0]  = -99
	    self.q_2[0]  = -99
	    self.isGlobal_2[0]  = -99
	    try:
		self.L1PreFiringWeight_Nom[0] = entry.L1PreFiringWeight_Nom
		self.L1PreFiringWeight_Up[0] = entry.L1PreFiringWeight_Up
		self.L1PreFiringWeight_Down[0] = entry.L1PreFiringWeight_Dn
	    except AttributeError : 
		self.L1PreFiringWeight_Nom[0] = 1
		self.L1PreFiringWeight_Up[0] = 1
		self.L1PreFiringWeight_Down[0] = 1
		
	    self.tightId_1[0]       = -1 
	    self.mediumId_1[0]       = -1 
	    self.mediumPromptId_1[0]   = -1
	    self.looseId_1[0]       = -1
	    self.isGlobal_1[0]      = -1
	    self.isTracker_1[0]     = -1

	    self.GenPart_statusFlags_1[0]    = -1
	    self.GenPart_status_1[0]    = -1
	    self.gen_match_1[0] = -1

	    self.tightId_2[0]       = -1 
	    self.mediumId_2[0]       = -1 
	    self.mediumPromptId_2[0]   = -1
	    self.looseId_2[0]       = -1
	    self.isGlobal_2[0]      = -1
	    self.isTracker_2[0]     = -1

	    self.GenPart_statusFlags_2[0]    = -1
	    self.GenPart_status_2[0]    = -1
	    self.gen_match_2[0] = -1

	    try :
		self.weight[0]           = entry.genWeight
		self.LHEweight[0]        = entry.LHEWeight_originalXWGTUP
		self.Generator_weight[0] = entry.Generator_weight
		self.LHE_Njets[0]        = ord(entry.LHE_Njets)
                if SystIndex == 0 : 
		    for i in range(0, int(entry.nLHEScaleWeight)) : 
			self.LHEScaleWeights[i] = entry.LHEScaleWeight[i]

		self.nPU[0]  = entry.Pileup_nPU
		self.nPUEOOT[0]  = entry.Pileup_sumEOOT
		self.nPULOOT[0]  = entry.Pileup_sumLOOT
		self.nPUtrue[0]  = entry.Pileup_nTrueInt
		self.nPV[0]  = entry.PV_npvs
		self.nPVGood[0]  = entry.PV_npvsGood
			    
	    except AttributeError :
		self.weight[0]           = 1. 
		self.weightPU[0]         = -1
		self.weightPUtrue[0]     = -1
		self.LHEweight[0]        = 1. 
		self.Generator_weight[0] = 1.
		self.LHE_Njets[0] = -1
		self.nPU[0]  = -1
		self.nPUEOOT[0]  = -1
		self.nPULOOT[0]  = -1
		self.nPUtrue[0]  = -1
		self.nPV[0]  = -1
		self.nPVGood[0]  = -1

        e = entry

        '''
        List from Cecile 
        single ele 2016: HLT Ele25 eta2p1 WPTight Gsf v and cut pt(ele)>26, eta(ele)<2.1
        single ele 2017: HLT Ele27 WPTight Gsf v, HLT Ele32 WPTight Gsf v, HLT Ele35 WPTight Gsf v and cut pt(ele)>28, eta(ele)<2.1
        single ele 2018: HLT Ele32 WPTight Gsf v, HLT Ele35 WPTight Gsf v and cut pt(ele)>33, eta(ele)<2.1
        '''
        
        if int(SystIndex) ==0 : 
	    bits=[]
	    try : bits.append(e.HLT_Ele25_eta2p1_WPTight_Gsf)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_Ele27_WPTight_Gsf)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_Ele32_WPTight_Gsf)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_Ele35_WPTight_Gsf)
	    except AttributeError : bits.append(False)
	    # pad upper bits in this byte with zeros (False) 
	    #for i in range(4) :
	    #    bits.append(False)
		
	    try : bits.append(e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ)
	    except AttributeError : bits.append(False) 

	    self.electronTriggerWord[0] = 0
	    for i, bit in enumerate(bits) :
		if bit : self.electronTriggerWord[0] += 2**i

	    '''
	    List from Cecile 
	    single mu 2016: HLT IsoMu22 v, HLT IsoMu22 eta2p1 v, HLT IsoTkMu22 v, HLT IsoTkMu22 eta2p1 v and cut pt(mu)>23, eta(mu)<2.1
	    single mu 2017: HLT IsoMu24 v, HLT IsoMu27 v and cut pt(mu)>25, eta(mu)<2.4
	    single mu 2018: HLT IsoMu24 v, HLT IsoMu27 v and cut pt(mu)>25, eta(mu)<2.4
	    '''
	    bits=[]
	    try : bits.append(e.HLT_IsoMu22)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_IsoMu22_eta2p1)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_IsoTkMu22)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_IsoTkMu22_eta2p1)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_IsoMu24)
	    except AttributeError : bits.append(False) 
	    try : bits.append(e.HLT_IsoMu27)
	    except AttributeError : bits.append(False) 

	    #for i in range(2) :
	    #    bits.append(False)                             # pad remaining bit in this bit 
	   
	    try : bits.append(e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ)
	    except AttributeError : bits.append(False) 
	    try : bits.append(e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass12)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p12)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ)
	    except AttributeError : bits.append(False) 
	    try : bits.append(e.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_Mass12)
	    except AttributeError : bits.append(False) 

	    self.muonTriggerWord[0] = 0
	    for i, bit in enumerate(bits) :
		if bit : self.muonTriggerWord[0] += 2**i

            #neede for all systematics as jt1/jt2 may change per systematic
	self.cat[0]  = tauFun2.catToNumber2Lep(cat)


        self.pt_1[0]   = LepP.Pt()
        self.phi_1[0]  = LepP.Phi()
        self.eta_1[0]  = LepP.Eta()
        self.pt_2[0]   = LepM.Pt()
        self.phi_2[0]  = LepM.Phi()
        self.eta_2[0]  = LepM.Eta()

	lep_index_1 = lepList[0]
	lep_index_2 = lepList[1]

	if channel_ll == 'ee' : 
      
            self.iso_1[0]  = entry.Electron_pfRelIso03_all[lep_index_1]
            self.q_1[0]  = entry.Electron_charge[lep_index_1]
            self.d0_1[0]   = entry.Electron_dxy[lep_index_1]
            self.dZ_1[0]   = entry.Electron_dz[lep_index_1]
            self.Electron_mvaFall17V2noIso_WP90_1[0]  = entry.Electron_mvaFall17V2noIso_WP90[lep_index_1]
            self.Electron_cutBased_1[0]  = entry.Electron_cutBased[lep_index_1]
	    #if SystIndex ==0 and  isMC : 
	    #	self.pt_uncor_1[0] = ePt[lep_index_1]
	    #	self.m_uncor_1[0] = eMass[lep_index_1]

            if isMC :
		self.gen_match_1[0] = ord(entry.Electron_genPartFlav[lep_index_1])

            self.iso_2[0]  = entry.Electron_pfRelIso03_all[lep_index_2]
            self.q_2[0]  = entry.Electron_charge[lep_index_2]
            self.d0_2[0]   = entry.Electron_dxy[lep_index_2]
            self.dZ_2[0]   = entry.Electron_dz[lep_index_2]
            self.Electron_mvaFall17V2noIso_WP90_2[0]  = entry.Electron_mvaFall17V2noIso_WP90[lep_index_2]
            self.Electron_cutBased_2[0]  = entry.Electron_cutBased[lep_index_2]
	    #if SystIndex ==0 and  isMC : 
	    #	self.pt_uncor_2[0] = ePt[lep_index_2]
            #	self.m_uncor_2[0] = eMass[lep_index_2]

            if isMC :
		self.gen_match_2[0] = ord(entry.Electron_genPartFlav[lep_index_2])

	if channel_ll == 'mm' : 
            self.iso_1[0]  = entry.Muon_pfRelIso04_all[lep_index_1]
	    self.q_1[0]  = entry.Muon_charge[lep_index_1]
	    self.d0_1[0]   = entry.Muon_dxy[lep_index_1]
	    self.dZ_1[0]   = entry.Muon_dz[lep_index_1]
	    self.looseId_1[0]   = entry.Muon_looseId[lep_index_1] 
            self.tightId_1[0]      = entry.Muon_tightId[lep_index_1]
	    self.mediumId_1[0]   = entry.Muon_mediumId[lep_index_1] 
	    self.mediumPromptId_1[0]   = entry.Muon_mediumPromptId[lep_index_1] 
	    self.isGlobal_1[0]   = entry.Muon_isGlobal[lep_index_1] 
	    self.isTracker_1[0]   = entry.Muon_isTracker[lep_index_1] 
            if isMC :
		self.gen_match_1[0] = ord(entry.Muon_genPartFlav[lep_index_1])

            self.iso_2[0]  = entry.Muon_pfRelIso04_all[lep_index_2]
	    self.q_2[0]  = entry.Muon_charge[lep_index_2]
	    self.d0_2[0]   = entry.Muon_dxy[lep_index_2]
	    self.dZ_2[0]   = entry.Muon_dz[lep_index_2]
	    self.looseId_2[0]   = entry.Muon_looseId[lep_index_2] 
            self.tightId_2[0]      = entry.Muon_tightId[lep_index_2]
	    self.mediumId_2[0]   = entry.Muon_mediumId[lep_index_2] 
	    self.mediumPromptId_2[0]   = entry.Muon_mediumPromptId[lep_index_2] 
	    self.isGlobal_2[0]   = entry.Muon_isGlobal[lep_index_2] 
	    self.isTracker_2[0]   = entry.Muon_isTracker[lep_index_2] 
            if isMC :
		self.gen_match_2[0] = ord(entry.Muon_genPartFlav[lep_index_2])


        
        # genMatch the di-lepton variables
	if isMC :
	    idx_Lep1 = -1
	    idx_Lep1_tr = -1
	    if (LepP.M() > 0.05): # muon mass 
		idx_Lep1 = GF.getLepIdxFrom4Vec(entry, LepP, 'm')
		try :
		    idx_Lep1_tr = entry.Muon_genPartIdx[idx_Lep1]
		except IndexError : pass 
		    
	    elif (LepP.M() < 0.05 < 0.05): # electron mass
		idx_Lep1 = GF.getLepIdxFrom4Vec(entry, LepP, 'e')
		try :
		    idx_Lep1_tr = entry.Electron_genPartIdx[idx_Lep1]
		except IndexError : pass 
		    
	    if idx_Lep1_tr >= 0 :
		self.m_1_tr[0]  = entry.GenPart_mass[idx_Lep1_tr]
		self.pt_1_tr[0]  = entry.GenPart_pt[idx_Lep1_tr]
		self.eta_1_tr[0] = entry.GenPart_eta[idx_Lep1_tr]
		self.phi_1_tr[0] = entry.GenPart_phi[idx_Lep1_tr]
		self.GenPart_statusFlags_1[0]    = entry.GenPart_statusFlags[idx_Lep1_tr]
		self.GenPart_status_1[0]    = entry.GenPart_status[idx_Lep1_tr]

	    idx_Lep2 = -1
	    idx_Lep2_tr = -1
	    if (LepM.M() > 0.05): # muon mass 
		idx_Lep2 = GF.getLepIdxFrom4Vec(entry, LepM, 'm')
		try :
		    idx_Lep2_tr = entry.Muon_genPartIdx[idx_Lep2]
		except IndexError : pass 
		    
	    elif (LepM.M() < 0.05 < 0.05): # electron mass
		idx_Lep2 = GF.getLepIdxFrom4Vec(entry, LepM, 'e')
		try :
		    idx_Lep2_tr = entry.Electron_genPartIdx[idx_Lep2]
		except IndexError : pass 
		    
	    if idx_Lep2_tr >= 0 :
		self.m_2_tr[0]  = entry.GenPart_mass[idx_Lep2_tr]
		self.pt_2_tr[0]  = entry.GenPart_pt[idx_Lep2_tr]
		self.eta_2_tr[0] = entry.GenPart_eta[idx_Lep2_tr]
		self.phi_2_tr[0] = entry.GenPart_phi[idx_Lep2_tr]
		self.GenPart_statusFlags_2[0]    = entry.GenPart_statusFlags[idx_Lep2_tr]
		self.GenPart_status_2[0]    = entry.GenPart_status[idx_Lep2_tr]
        
        

        if True:  
            self.MET_significance[0]= entry.MET_significance

	    if not doUncertainties : 
		if proc =='UL' or str(era) != '2017': 
		    self.MET_pt[0]= entry.MET_pt
		    self.MET_phi[0]= entry.MET_phi

		if str(era) == '2017' and proc !='UL': 
		    try :
			self.MET_pt[0]= entry.METFixEE2017_pt
			self.MET_phi[0]= entry.METFixEE2017_phi

		    except AttributeError:
			self.MET_pt[0]= entry.MET_pt
			self.MET_phi[0]= entry.MET_phi

	    if  doUncertainties : 

		if proc =='UL' or str(era) != '2017': 
                    try : 
			self.MET_T1_pt[0]= entry.MET_T1_pt
			self.MET_T1_phi[0]= entry.MET_T1_phi
                    except AttributeError : 
			self.MET_T1_pt[0]= entry.MET_pt
			self.MET_T1_phi[0]= entry.MET_phi

		if str(era) == '2017' and proc !='UL': 
		    try : 
			self.MET_T1_pt[0]= entry.METFixEE2017_T1_pt
			self.MET_T1_phi[0]= entry.METFixEE2017_T1_phi
		    except AttributeError : 
			self.MET_T1_pt[0]= entry.METFixEE2017_pt
			self.MET_T1_phi[0]= entry.METFixEE2017_phi

        #metNoTauES holds the uncorrected TauES MET - if not doUncerta -> holds the default ucorrected MET, if doUncert the T1_corrected


        if str(era) != '2017' or proc=='UL': 

	    self.metcov00[0] = entry.MET_covXX
	    self.metcov01[0] = entry.MET_covXY
	    self.metcov10[0] = entry.MET_covXY
	    self.metcov11[0] = entry.MET_covYY

	    self.PuppiMET_pt[0]= entry.PuppiMET_pt
	    self.PuppiMET_phi[0]= entry.PuppiMET_phi
	    self.RawPuppiMET_pt[0] = entry.RawPuppiMET_pt
	    self.RawPuppiMET_phi[0] = entry.RawPuppiMET_phi
	    self.RawMET_pt[0]= entry.RawMET_pt
	    self.RawMET_phi[0]= entry.RawMET_phi

            metV, metUn =  TLorentzVector(), TLorentzVector()
            metUn.SetXYZT(entry.MET_MetUnclustEnUpDeltaX,entry.MET_MetUnclustEnUpDeltaY,0,0)
            metV.SetPtEtaPhiM(entry.MET_pt,0, entry.MET_phi,0)

	    self.MET_pt_UnclusteredUp[0] =  (metV+metUn).Pt()
	    self.MET_pt_UnclusteredDown[0] =  (metV-metUn).Pt()
	    self.MET_phi_UnclusteredUp[0] =  (metV+metUn).Phi()
	    self.MET_phi_UnclusteredDown[0] =  (metV-metUn).Phi()
 
	    self.PuppiMET_ptJESUp[0] = entry.PuppiMET_ptJESUp
	    self.PuppiMET_ptJESDown[0] = entry.PuppiMET_ptJESDown
	    self.PuppiMET_ptJERUp[0] = entry.PuppiMET_ptJERUp
	    self.PuppiMET_ptJERDown[0] = entry.PuppiMET_ptJERDown

	    self.PuppiMET_phiJESUp[0] = entry.PuppiMET_phiJESUp
	    self.PuppiMET_phiJESDown[0] = entry.PuppiMET_phiJESDown
	    self.PuppiMET_phiJERUp[0] = entry.PuppiMET_phiJERUp
	    self.PuppiMET_phiJERDown[0] = entry.PuppiMET_phiJERDown

	    self.PuppiMET_pt_UnclusteredUp[0] = entry.PuppiMET_ptUnclusteredUp
	    self.PuppiMET_pt_UnclusteredDown[0] = entry.PuppiMET_ptUnclusteredDown
	    self.PuppiMET_phi_UnclusteredUp[0] = entry.PuppiMET_phiUnclusteredUp
	    self.PuppiMET_phi_UnclusteredDown[0] = entry.PuppiMET_phiUnclusteredDown
 
            #do the default ones, without corrections

	    if doUncertainties : 

		if True : 
                    #print 'met entries', entry.MET_T1_pt, entry.MET_T1_pt_jesTotalUp, entry.MET_T1_pt_jesTotalDown, entry.MET_T1_pt_jerUp, entry.MET_T1_pt_jerDown, entry.event
                    '''
		    self.MET_T1_pt_jesTotalUp[0] = entry.MET_T1_pt_jesTotalUp
		    self.MET_T1_pt_jesTotalDown[0] = entry.MET_T1_pt_jesTotalDown
		    self.MET_T1_phi_jesTotalUp[0] = entry.MET_T1_phi_jesTotalUp
		    self.MET_T1_phi_jesTotalDown[0] = entry.MET_T1_phi_jesTotalDown
                    '''
                    '''
		    self.MET_T1_pt_jerUp[0] = entry.MET_T1_pt_jerUp
		    self.MET_T1_pt_jerDown[0] = entry.MET_T1_pt_jerDown
		    self.MET_T1_phi_jerUp[0] = entry.MET_T1_phi_jerUp
		    self.MET_T1_phi_jerDown[0] = entry.MET_T1_phi_jerDown
                    '''

		    if isMC: 
			self.MET_T1_pt_UnclusteredUp[0] = entry.MET_T1_pt_unclustEnUp
			self.MET_T1_pt_UnclusteredDown[0] = entry.MET_T1_pt_unclustEnDown
			self.MET_T1_phi_UnclusteredUp[0] = entry.MET_T1_phi_unclustEnUp
			self.MET_T1_phi_UnclusteredDown[0] = entry.MET_T1_phi_unclustEnDown



        if  str(era) == '2017' and proc=='EOY' :  

	    self.metcov00[0] = entry.METFixEE2017_covXX

	    self.metcov00[0] = entry.METFixEE2017_covXX
	    self.metcov01[0] = entry.METFixEE2017_covXY
	    self.metcov10[0] = entry.METFixEE2017_covXY
	    self.metcov11[0] = entry.METFixEE2017_covYY

	    self.PuppiMET_pt[0]= entry.PuppiMETFixEE2017_pt
	    self.PuppiMET_phi[0]= entry.PuppiMETFixEE2017_phi
	    self.RawPuppiMETFixEE2017_pt[0] = entry.RawPuppiMETFixEE2017_pt
	    self.RawPuppiMETFixEE2017_phi[0] = entry.RawPuppiMETFixEE2017_phi
	    self.RawMETFixEE2017_pt[0]= entry.RawMETFixEE2017_pt
	    self.RawMETFixEE2017_phi[0]= entry.RawMETFixEE2017_phi

            metV, metUn =  TLorentzVector(), TLorentzVector()
            metUn.SetXYZT(entry.METFixEE2017_MetUnclustEnUpDeltaX,entry.METFixEE2017_MetUnclustEnUpDeltaY,0,0)
            metV.SetPtEtaPhiM(entry.METFixEE2017_pt,0, entry.METFixEE2017_phi,0)

	    self.MET_pt_UnclusteredUp[0] =  (metV+metUn).Pt()
	    self.MET_pt_UnclusteredDown[0] =  (metV-metUn).Pt()
	    self.MET_phi_UnclusteredUp[0] =  (metV+metUn).Phi()
	    self.MET_phi_UnclusteredDown[0] =  (metV-metUn).Phi()
 
	    self.PuppiMET_ptJESUp[0] = entry.PuppiMETFixEE2017_ptJESUp
	    self.PuppiMET_ptJESDown[0] = entry.PuppiMETFixEE2017_ptJESDown
	    self.PuppiMET_ptJERUp[0] = entry.PuppiMETFixEE2017_ptJERUp
	    self.PuppiMET_ptJERDown[0] = entry.PuppiMETFixEE2017_ptJERDown

	    self.PuppiMET_phiJESUp[0] = entry.PuppiMETFixEE2017_phiJESUp
	    self.PuppiMET_phiJESDown[0] = entry.PuppiMETFixEE2017_phiJESDown
	    self.PuppiMET_phiJERUp[0] = entry.PuppiMETFixEE2017_phiJERUp
	    self.PuppiMET_phiJERDown[0] = entry.PuppiMETFixEE2017_phiJERDown

	    self.PuppiMET_pt_UnclusteredUp[0] = entry.PuppiMETFixEE2017_ptUnclusteredUp
	    self.PuppiMET_pt_UnclusteredDown[0] = entry.PuppiMETFixEE2017_ptUnclusteredDown
	    self.PuppiMET_phi_UnclusteredUp[0] = entry.PuppiMETFixEE2017_phiUnclusteredUp
	    self.PuppiMET_phi_UnclusteredDown[0] = entry.PuppiMETFixEE2017_phiUnclusteredDown


	    if doUncertainties : 
		if isMC : 
                    '''
		    self.MET_T1_pt_jesTotalUp[0] = entry.METFixEE2017_T1_pt_jesTotalUp
		    self.MET_T1_pt_jesTotalDown[0] = entry.METFixEE2017_T1_pt_jesTotalDown
		    self.MET_T1_phi_jesTotalUp[0] = entry.METFixEE2017_T1_phi_jesTotalUp
		    self.MET_T1_phi_jesTotalDown[0] = entry.METFixEE2017_T1_phi_jesTotalDown
                    '''

		    self.MET_T1_pt_UnclusteredUp[0] = entry.METFixEE2017_T1_pt_unclustEnUp
		    self.MET_T1_pt_UnclusteredDown[0] = entry.METFixEE2017_T1_pt_unclustEnDown
		    self.MET_T1_phi_UnclusteredUp[0] = entry.METFixEE2017_T1_phi_unclustEnUp
		    self.MET_T1_phi_UnclusteredDown[0] = entry.METFixEE2017_T1_phi_unclustEnDown

                    '''
		    self.MET_T1_pt_jerUp[0] = entry.METFixEE2017_T1_pt_jerUp
		    self.MET_T1_pt_jerDown[0] = entry.METFixEE2017_T1_pt_jerDown
		    self.MET_T1_phi_jerUp[0] = entry.METFixEE2017_T1_phi_jerUp
		    self.MET_T1_phi_jerDown[0] = entry.METFixEE2017_T1_phi_jerDown
                    '''



        # trig
        if SystIndex ==0 : 
	    self.isTrig_1[0]   = is_trig_1

        leplist=[]
        leplist.append(LepP)
        leplist.append(LepM)

        if doUncertainties: 
                ## this is not done from within ZH and the correctallMET function
                for i, v in enumerate(self.allsystMET) : 

                    if str(era)=='2017' and proc !='UL':
                        #i_ should be the righ-hand of the branch and should retain the METFixEE2017 if y=2017 
                        #iMET should appear always at the branch name...
                        v = v.replace('MET','METFixEE2017')
                    iMET= v.replace('METFixEE2017','MET')

                    try : j = getattr(entry, "{0:s}".format(str(v)))
                    except AttributeError : j = -9.99
		    self.list_of_arrays_noES[i][0] = j
                    #if '_pt_jerUp' in v  : print '=====================================while filling-----------------',j, self.list_of_arrays[i][0], i, v, entry.event 

                for i, v in enumerate(self.allsystJets) : 
                #njets_sys, nbtag_sys
		    jetList, jetListFlav, jetListEta, jetListPt, bTagListDeep, bJetListL,bJetListM, bJetListT, bJetListFlav = self.getJetsJMEMV(entry,leplist,era,v) 
                    #print 'jessyst', systematic, len(jetList), cat

	            self.list_of_arraysJetsNjets[i][0] = len(jetList)
	            self.list_of_arraysJetsNbtagL[i][0] = len(bJetListL)
	            self.list_of_arraysJetsNbtagM[i][0] = len(bJetListM)
	            self.list_of_arraysJetsNbtagT[i][0] = len(bJetListT)
		    for ifl in range(len(jetList)) :
			self.list_of_arraysJetsPt[i][ifl] = jetListPt[ifl]
			self.list_of_arraysJetsEta[i][ifl] = jetListEta[ifl]
			self.list_of_arraysJetsFlavour[i][ifl] = jetListFlav[ifl]
	                self.list_of_arraysJetsNbtagDeep[i][ifl] = bTagListDeep[ifl]


        #fill the un-corrected or just in the case you dont care to doUncertainties       
        nom_=''
	jetList, jetListFlav, jetListEta, jetListPt, bTagListDeep, bJetListL, bJetListM, bJetListT, bJetListFlav = self.getJetsJMEMV(entry,leplist,era,'') 
	self.njets[0] = len(jetList)
	self.nbtagL[0] = len(bJetListL)
	self.nbtagM[0] = len(bJetListM)
	self.nbtagT[0] = len(bJetListT)
	for ifl in range(len(jetListPt)) :
	    self.jflavour[ifl]  = jetListFlav[ifl]
	    self.jeta[ifl]  = jetListEta[ifl]
	    self.jpt[ifl]  = jetListPt[ifl]
	    self.btagDeep[ifl] = bTagListDeep[ifl]
        '''
        if self.nMuon[0] !=1 : print 'failed nMuon', entry.event, self.nMuon[0]
        if self.nTau[0]!=0 :  print 'failed nTau', entry.event, self.nTau[0]
        if (self.isGlobal_1[0]<1 and self.isTracker_1[0]<1) :  print 'failed muon Id', entry.event, self.isGlobal_1[0], self.isTracker_1[0]
        if self.pt_1[0]<29 :  print 'failed pt', entry.event, self.pt_1[0]
        if self.mediumId_1[0]<1 :  print 'failed mediumId', entry.event, self.mediumId_1[0]
        if fabs(self.eta_1[0])>2.4 or  fabs(self.dZ_1[0])>0.2 or fabs(self.d0_1[0])>0.045:  print 'failed eta, dZ, d0', entry.event, fabs(self.eta_1[0]), fabs(self.dZ_1[0]),  fabs(self.d0_1[0])
        if self.iso_1[0]>0.15 :  print 'failed iSo', entry.event, self.iso_1[0]
        if self.isTrig_1[0]<1 :  print 'failed Trig', entry.event,  self.isTrig_1[0]
        if self.nbtagT[0]>0 :  print 'failed bTag', entry.event, self.nbtagT[0]
        if fabs(self.q_1[0]) !=1 :  print 'failed charge', entry.event, fabs(self.q_1[0])
        '''
	SaveEventMu = False
	SaveEventEl = False
	if '2016' in era:  
            SaveEventMu =  cat=='mm' and (self.isGlobal_1[0]>0 or self.isTracker_1[0]>0) and self.pt_1[0]>24 and self.mediumId_1[0]>0 and fabs(self.eta_1[0])<2.4 and  fabs(self.dZ_1[0])<0.2 and  fabs(self.d0_1[0])<0.045 and self.iso_1[0]<0.2   and (self.isGlobal_2[0]>0 or self.isTracker_2[0]>0) and self.pt_2[0]>24 and self.mediumId_2[0]>0 and fabs(self.eta_2[0])<2.4 and  fabs(self.dZ_2[0])<0.2 and  fabs(self.d0_2[0])<0.045 and self.iso_2[0]<0.2 and self.isTrig_1[0]!=0 and self.nbtagT[0]==0

            #SaveEventEl = cat=='ee' and  self.pt_1[0]>27 and self.Electron_mvaFall17V2Iso_WP90_1[0]>0 and fabs(self.eta_1[0])<2.1 and  fabs(self.dZ_1[0])<0.2 and  fabs(self.d0_1[0])<0.045 and self.isTrig_1[0]>0 and self.iso_1[0] < 0.5 #and fabs(self.PVz[0])<26  and (self.PVy[0]*self.PVy[0] + self.PVx[0]*self.PVx[0])<3 and self.nPV[0]>2 

            SaveEventEl = cat=='ee' and  self.pt_1[0]>27 and self.Electron_mvaFall17V2Iso_WP90_1[0]>0 and fabs(self.eta_1[0])<2.1 and  fabs(self.dZ_1[0])<0.2 and  fabs(self.d0_1[0])<0.045 and self.iso_1[0] < 0.15 and  self.pt_2[0]>27 and self.Electron_mvaFall17V2Iso_WP90_2[0]>0 and fabs(self.eta_2[0])<2.1 and  fabs(self.dZ_2[0])<0.2 and  fabs(self.d0_2[0])<0.045 and self.iso_2[0] < 0.15        and self.isTrig_1[0]!=0

        #SaveEvent = self.nMuon[0]==1 and self.nTau[0]==0 and (self.isGlobal_1[0]>0 or self.isTracker_1[0]>0) and self.pt_1[0]>29 and self.mediumId_1[0]>0 and fabs(self.eta_1[0])<2.4 and  fabs(self.dZ_1[0])<0.2 and  fabs(self.d0_1[0])<0.045 and self.iso_1[0]<0.2 and self.isTrig_1[0]>0 and self.nbtagT[0]==0
        #SaveEvent = self.nTau[0]==0 and (self.isGlobal_1[0]>0 or self.isTracker_1[0]>0) and self.pt_1[0]>24 and self.mediumId_1[0]>0 and fabs(self.eta_1[0])<2.4 and  fabs(self.dZ_1[0])<0.2 and  fabs(self.d0_1[0])<0.045 and self.iso_1[0]<0.2   and (self.isGlobal_2[0]>0 or self.isTracker_2[0]>0) and self.pt_2[0]>24 and self.mediumId_2[0]>0 and fabs(self.eta_2[0])<2.4 and  fabs(self.dZ_2[0])<0.2 and  fabs(self.d0_2[0])<0.045 and self.iso_2[0]<0.2 and self.isTrig_1[0]>0 and self.nbtagT[0]==0
        SaveEvent =  (self.isGlobal_1[0]>0 or self.isTracker_1[0]>0) and self.pt_1[0]>24 and self.mediumId_1[0]>0 and fabs(self.eta_1[0])<2.4 and  fabs(self.dZ_1[0])<0.2 and  fabs(self.d0_1[0])<0.045 and self.iso_1[0]<0.2   and (self.isGlobal_2[0]>0 or self.isTracker_2[0]>0) and self.pt_2[0]>24 and self.mediumId_2[0]>0 and fabs(self.eta_2[0])<2.4 and  fabs(self.dZ_2[0])<0.2 and  fabs(self.d0_2[0])<0.045 and self.iso_2[0]<0.2 and self.isTrig_1[0]!=0 and self.nbtagT[0]==0
        #SaveEvent = True
        #print 'filling',  self.mll, SaveEvent, SystIndex
        #print self.nTau[0], self.isGlobal_1[0], self.isTracker_1[0], self.pt_1[0], self.mediumId_1[0], fabs(self.eta_1[0]), fabs(self.dZ_1[0]), fabs(self.d0_1[0]), self.iso_1[0], 'second',  (self.isGlobal_2[0]>0 or self.isTracker_2[0]>0), self.pt_2[0], self.mediumId_2[0], fabs(self.eta_2[0]), fabs(self.dZ_2[0]), fabs(self.d0_2[0]), self.iso_2[0], self.isTrig_1[0], self.nbtagT[0],  self.gen_match_1[0],  self.gen_match_2[0]

        if SaveEvent:
            #print 'I will save this......', entry.event
            #if  self.nbtag[0] == 0 : 
	    if SystIndex == 0 : 
		self.t.Fill()
	    else : 
		self.tN[SystIndex-1].Fill()

	return


    def setWeight(self,weight) :
        self.weight[0] = weight
        #print("outTuple.setWeight() weight={0:f}".format(weight))
        return
    def setWeightPU(self,weight) :
        self.weightPU[0] = weight
        #print("outTuple.setWeight() weight={0:f}".format(weight))
        return
    def setWeightPUtrue(self,weight) :
        self.weightPUtrue[0] = weight
        #print("outTuple.setWeight() weight={0:f}".format(weight))
        return

    def FillTree(self) :
        self.t.Fill()

    def writeTree(self) :
        print("In outTuple.writeTree() entries={0:d}".format(self.entries))
        self.f.Write()
        self.f.Close()
        return



