# output ntuple for H->tautau analysis for CMSSW_10_2_X

from ROOT import TLorentzVector, TH1
from math import sqrt, sin, cos, pi
import tauFunDCH
import ROOT, array
import os
import sys
import generalFunctions as GF


electronMass = 0.0005
muonMass  = 0.105
class outTuple() :
    
    def __init__(self,fileName, era, doSyst=False,shift=[], isMC=True, onlyNom=False):
        from array import array
        from ROOT import TFile, TTree

        # Tau Decay types
        self.kUndefinedDecayType, self.kTauToHadDecay,  self.kTauToElecDecay, self.kTauToMuDecay = 0, 1, 2, 3    
        ROOT.gInterpreter.ProcessLine(".include .")
        for baseName in ['MeasuredTauLepton','svFitAuxFunctions','FastMTT'] : 
            if os.path.isfile("{0:s}_cc.so".format(baseName)) :
                ROOT.gInterpreter.ProcessLine(".L {0:s}_cc.so".format(baseName))
            else :
                ROOT.gInterpreter.ProcessLine(".L {0:s}.cc++".format(baseName))   
                # .L is not just for .so files, also .cc
       
        ########### JetMet systematics
	#self.listsyst=['njets', 'nbtag', 'jpt', 'jeta', 'jflavour','MET_T1_pt', 'MET_T1_phi', 'MET_pt', 'MET_phi', 'MET_T1Smear_pt', 'MET_T1Smear_phi']
        self.jessyst=['_nom']
	self.listsyst=['njets', 'nbtagL', ',nbtagM', 'btagDeep','nbtagT','jpt', 'jeta', 'jflavour','MET_T1_pt', 'MET_T1_phi', 'MET_pt', 'MET_phi']
        if doSyst :
	    self.jessyst=['_nom','_jesAbsolute', '_jesAbsolute_{0:s}'.format(str(era)), '_jesBBEC1', '_jesBBEC1_{0:s}'.format(str(era)), '_jesEC2', '_jesEC2_{0:s}'.format(str(era)), '_jesFlavorQCD', '_jesHF', '_jesHF_{0:s}'.format(str(era)), '_jesRelativeBal', '_jesRelativeSample_{0:s}'.format(str(era)), '_jesHEMIssue', '_jesTotal', '_jer']  

        if onlyNom :
	    self.jessyst=['_nom']
        #shift are the ES basd systematics


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

        #if not isMC or 'ZHTo' in str(fileName):
        if not isMC  :
        
	    self.listsyst=['njets', 'nbtagL', ',nbtagM', 'nbtagT','btagDeep','jpt', 'jeta', 'jflavour', 'MET_pt', 'MET_phi']
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
        self.evt              = array('l',[0])
        self.nPU              = array('l',[0])
        self.nPUEOOT              = array('l',[0])
        self.nPULOOT              = array('l',[0])
        self.nPUtrue              = array('f',[0])
        self.nPV              = array('l',[0])
        self.nPVGood              = array('l',[0])
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

        self.L1PreFiringWeight_Nom        = array('f',[0])
        self.L1PreFiringWeight_Up        = array('f',[0])
        self.L1PreFiringWeight_Down        = array('f',[0])

        self.d0_1        = array('f',[0])
        self.dZ_1        = array('f',[0])
        self.d0_2        = array('f',[0])
        self.dZ_2        = array('f',[0])
        
        self.pt_3        = array('f',[0])
        self.pt_3_tr     = array('f',[0])
        self.pt_uncor_1        = array('f',[0])
        self.pt_uncor_2        = array('f',[0])
        self.pt_uncor_3        = array('f',[0])
        self.pt_uncor_4        = array('f',[0])
        self.m_uncor_1        = array('f',[0])
        self.m_uncor_2        = array('f',[0])
        self.m_uncor_3        = array('f',[0])
        self.m_uncor_4        = array('f',[0])

        self.GenPart_statusFlags_3     = array('l',[0])
        self.GenPart_status_3     = array('l',[0])
        self.phi_3       = array('f',[0])
        self.phi_3_tr    = array('f',[0])
        self.eta_3       = array('f',[0])
        self.eta_3_tr    = array('f',[0])
        self.m_3         = array('f',[0])
        self.q_3         = array('f',[0])
        self.d0_3        = array('f',[0])
        self.dZ_3        = array('f',[0])
        self.mt_3        = array('f',[0])
        self.pfmt_3      = array('f',[0])
        #self.puppimt_3   = array('f',[0])
        self.iso_3       = array('f',[0])
        self.Electron_mvaFall17V2noIso_WP90_1 = array('f',[0])
        self.Electron_mvaFall17V2noIso_WP90_2 = array('f',[0])
        self.Electron_mvaFall17V2noIso_WP90_3 = array('f',[0])
        self.Electron_mvaFall17V2noIso_WP90_4 = array('f',[0])
        self.gen_match_1 = array('l',[0])
        self.gen_match_2 = array('l',[0])
        self.gen_match_3 = array('l',[0])
        self.tightId_3       = array('f',[0])
        self.mediumId_3       = array('f',[0])
        self.mediumPromptId_3       = array('f',[0])
        self.looseId_3       = array('f',[0])
        self.isGlobal_3       = array('f',[0])
        self.isTracker_3       = array('f',[0])
        self.ip3d_3       = array('f',[0])

        self.idDecayModeNewDMs_3 = array('f',[0])
        self.idDeepTau2017v2p1VSe_3 = array('f',[0])
        self.idDeepTau2017v2p1VSjet_3 = array('f',[0])
        self.idDeepTau2017v2p1VSmu_3 = array('f',[0])
        self.idMVAnewDM2017v2_3 = array('f',[0])
        self.rawMVAnewDM2017v2_3 = array('f',[0])


        self.decayMode_3   = array('l',[0])

        self.pt_4        = array('f',[0])
        self.pt_4_tr     = array('f',[0])
        self.GenPart_statusFlags_4     = array('l',[0])
        self.GenPart_status_4     = array('l',[0])
        self.phi_4       = array('f',[0])
        self.phi_4_tr    = array('f',[0])
        self.eta_4       = array('f',[0])
        self.eta_4_tr    = array('f',[0])
        self.m_4         = array('f',[0])
        self.q_4         = array('f',[0])
        self.d0_4        = array('f',[0])
        self.dZ_4        = array('f',[0])
        self.mt_4        = array('f',[0])
        self.pfmt_4      = array('f',[0])
        #self.puppimt_4   = array('f',[0])
        self.iso_4       = array('f',[0])
        self.gen_match_4 = array('l',[0])
        self.tightId_4       = array('f',[0])
        self.mediumId_4       = array('f',[0])
        self.mediumPromptId_4       = array('f',[0])
        self.looseId_4       = array('f',[0])
        self.isGlobal_4       = array('f',[0])
        self.isTracker_4       = array('f',[0])
        self.ip3d_4       = array('f',[0])


        self.idDecayModeNewDMs_4 = array('f',[0])
        self.idDeepTau2017v2p1VSe_4 = array('f',[0])
        self.idDeepTau2017v2p1VSjet_4 = array('f',[0])
        self.idDeepTau2017v2p1VSmu_4 = array('f',[0])
        self.idMVAnewDM2017v2_4 = array('f',[0])
        self.rawMVAnewDM2017v2_4 = array('f',[0])


        '''
        self.pt_5        = array('f',[0])
        self.phi_5       = array('f',[0])
        self.eta_5       = array('f',[0])
        self.m_5         = array('f',[0])
        self.q_5         = array('f',[0])
        self.d0_5        = array('f',[0])
        self.dZ_5      = array('f',[0])
        self.gen_match_5 = array('l',[0])
        self.decayMode_5 = array('l',[0])

        self.idDecayModeNewDMs_5 = array('f',[0])
        self.idDeepTau2017v2p1VSe_5 = array('f',[0])
        self.idDeepTau2017v2p1VSjet_5 = array('f',[0])
        self.idDeepTau2017v2p1VSmu_5 = array('f',[0])
        self.idMVAnewDM2017v2_5 = array('f',[0])
        self.rawMVAnewDM2017v2_5 = array('f',[0])
        '''
        self.decayMode_4   = array('l',[0])

        # di-tau variables
        self.pt_tt  = array('f',[0])
        self.mt_tot = array('f',[0])
        self.m_vis  = array('f',[0])
        self.m_sv   = array('f',[0])
        self.mt_sv  = array('f',[0])
        self.H_DR  = array('f',[0])
        self.AMass   = array('f',[0])


        # di-lepton variables.   1 and 2 refer to plus and minus charge
        # ll_lmass is mass of decay lepton 
        self.H_LT       = array('f',[0])
        self.dRl1H       = array('f',[0])
        self.dRl2H       = array('f',[0])
        self.dRlH       = array('f',[0])
        self.dPhil1H       = array('f',[0])
        self.dPhil2H       = array('f',[0])
        self.dPhilH       = array('f',[0])
        self.mll       = array('f',[0])
        self.mll2       = array('f',[0])
        self.Z_Pt       = array('f',[0])
        self.Z_DR       = array('f',[0])
        self.Z_SS       = array('f',[0])
        self.pt_1      = array('f',[0])
        self.m_1_tr   = array('f',[0])
        self.pt_1_tr   = array('f',[0])
        self.GenPart_statusFlags_1   = array('l',[0])
        self.GenPart_status_1     = array('l',[0])
        self.phi_1     = array('f',[0])
        self.phi_1_tr  = array('f',[0])
        self.eta_1     = array('f',[0])
        self.eta_1_tr  = array('f',[0])
        self.pt_2      = array('f',[0])
        self.GenPart_statusFlags_2   = array('l',[0])
        self.GenPart_status_2     = array('l',[0])
        self.m_2_tr   = array('f',[0])
        self.pt_2_tr   = array('f',[0])
        self.phi_2     = array('f',[0])
        self.phi_2_tr  = array('f',[0])
        self.eta_2     = array('f',[0])
        self.eta_2_tr  = array('f',[0])
        self.iso_1       = array('f',[0])
        self.q_1       = array('f',[0])
        self.Muon_Id_1       = array('f',[0])
        self.Muon_Id_2       = array('f',[0])
        self.Muon_Id_3       = array('f',[0])
        self.isGlobal_1       = array('f',[0])
        self.isTracker_1       = array('f',[0])
        self.isTracker_2       = array('f',[0])
        self.isGlobal_2       = array('f',[0])
        self.tightId_1       = array('f',[0])
        self.mediumId_1       = array('f',[0])
        self.mediumPromptId_1       = array('f',[0])
        self.looseId_1       = array('f',[0])
        
        # MET variables
        self.met         = array('f',[0])
        self.metphi      = array('f',[0])
        self.metNoTauES      = array('f',[0])
        self.metphiNoTauES      = array('f',[0])
        self.metNoCor         = array('f',[0])
        self.metphiNoCor      = array('f',[0])
        #self.puppimet    = array('f',[0])
        #self.puppimetphi = array('f',[0])
        self.metcov00    = array('f',[0])
        self.metcov01    = array('f',[0])
        self.metcov10    = array('f',[0])
        self.metcov11    = array('f',[0])


        #systematics

        self.MET_pt_UnclUp = array('f',[0])
        self.MET_phi_UnclUp = array('f',[0])
        self.MET_pt_UnclDown = array('f',[0])
        self.MET_phi_UnclDown = array('f',[0])
        self.met_UnclX = array('f',[0])
        self.met_UnclY = array('f',[0])
        self.MET_T1Smear_pt= array('f',[0])
        self.MET_T1Smear_phi= array('f',[0])
        self.MET_pt_nom= array('f',[0])
        self.MET_pt_nom= array('f',[0])

        # trigger info
        self.isTrig_2   = array('f',[0])
        self.isTrig_1   = array('f',[0])
        self.isDoubleTrig   = array('f',[0])


        # jet variables
        #self.njetsold = array('f',[-1]*8)
        self.njets     = array('f',[0])
        self.nbtagL     = array('f',[0])
        self.nbtagM     = array('f',[0])
        self.nbtagT     = array('f',[0])

        self.HTXS_Higgs_cat     = array('l',[0])
        self.HTXS_Higgs_pt     = array('f',[0])

        '''
        self.jpt_1     = array('f',[0])
        self.jpt_1_tr  = array('f',[0])
        self.jeta_1    = array('f',[0])
        self.jeta_1_tr = array('f',[0])
        self.jphi_1    = array('f',[0])
        self.jphi_1_tr = array('f',[0])
        self.jcsv_1    = array('f',[0])
        self.jcsvfv_1    = array('f',[0])
        self.jpt_2     = array('f',[0])
        self.jpt_2_tr  = array('f',[0])
        self.jeta_2    = array('f',[0])
        self.jeta_2_tr = array('f',[0])
        self.jphi_2    = array('f',[0])
        self.jphi_2_tr = array('f',[0])
        self.jcsv_2    = array('f',[0])
        self.jcsvfv_2    = array('f',[0])
        '''
        self.iso_2       = array('f',[0])
        self.q_2       = array('f',[0])
        self.tightId_2       = array('f',[0])
        self.mediumId_2       = array('f',[0])
        self.mediumPromptId_2       = array('f',[0])
        self.looseId_2       = array('f',[0])
        self.jflavour     = array('f',[-9.99]*12)
        self.jeta     = array('f',[-9.99]*12)
        self.jpt     = array('f',[-9.99]*12)
        self.btagDeep     = array('f',[-9.99]*12)

        self.bpt_1     = array('f',[0]*8)
        self.bpt_1_tr  = array('f',[0]*8)
        self.beta_1    = array('f',[0]*8)
        self.beta_1_tr = array('f',[0]*8)
        self.bphi_1    = array('f',[0]*8)
        self.bphi_1_tr = array('f',[0]*8)
        self.bcsv_1    = array('f',[0]*8)
        self.bcsvfv_1    = array('f',[0]*8)
        self.bpt_2     = array('f',[0]*8)
        self.bpt_2_tr  = array('f',[0]*8)
        self.beta_2    = array('f',[0]*8)
        self.beta_2_tr = array('f',[0]*8)
        self.bphi_2    = array('f',[0]*8)
        self.bphi_2_tr = array('f',[0]*8)
        self.bcsv_2    = array('f',[0]*8)
        self.bcsvfv_2    = array('f',[0]*8)


      
        self.t.Branch('run',              self.run,               'run/l' )
        self.t.Branch('nElectron',              self.nElectron,               'nElectron/l' )
        self.t.Branch('nMuon',              self.nMuon,               'nMuon/l' )
        self.t.Branch('nTau',              self.nTau,               'nTau/l' )
        self.t.Branch('lumi',             self.lumi,              'lumi/I' )
        self.t.Branch('evt',              self.evt,               'evt/l' )
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
        
        self.t.Branch('GenPart_statusFlags_1',     self.GenPart_statusFlags_1,     'GenPart_statusFlags_1/I')
        self.t.Branch('GenPart_statusFlags_2',     self.GenPart_statusFlags_2,     'GenPart_statusFlags_2/I')
        self.t.Branch('GenPart_statusFlags_3',     self.GenPart_statusFlags_3,     'GenPart_statusFlags_3/I')
        self.t.Branch('GenPart_statusFlags_4',     self.GenPart_statusFlags_4,     'GenPart_statusFlags_4/I')
        self.t.Branch('GenPart_status_1',     self.GenPart_status_1,     'GenPart_status_1/I')
        self.t.Branch('GenPart_status_2',     self.GenPart_status_2,     'GenPart_status_2/I')
        self.t.Branch('GenPart_status_3',     self.GenPart_status_3,     'GenPart_status_3/I')
        self.t.Branch('GenPart_status_4',     self.GenPart_status_4,     'GenPart_status_4/I')
        self.t.Branch('pt_uncor_1',        self.pt_uncor_1,        'pt_uncor_1/F')
        self.t.Branch('pt_uncor_2',        self.pt_uncor_2,        'pt_uncor_2/F')
        self.t.Branch('pt_uncor_3',        self.pt_uncor_3,        'pt_uncor_3/F')
        self.t.Branch('pt_uncor_4',        self.pt_uncor_4,        'pt_uncor_4/F')
        self.t.Branch('m_uncor_1',        self.m_uncor_1,        'm_uncor_1/F')
        self.t.Branch('m_uncor_2',        self.m_uncor_2,        'm_uncor_2/F')
        self.t.Branch('m_uncor_3',        self.m_uncor_3,        'm_uncor_3/F')
        self.t.Branch('m_uncor_4',        self.m_uncor_4,        'm_uncor_4/F')
        self.t.Branch('pt_3',        self.pt_3,        'pt_3/F')
        self.t.Branch('pt_3_tr',     self.pt_3_tr,     'pt_3_tr/F')
        self.t.Branch('phi_3',       self.phi_3,       'phi_3/F')
        self.t.Branch('phi_3_tr',    self.phi_3_tr,    'phi_3_tr/F')
        self.t.Branch('eta_3',       self.eta_3,       'eta_3/F')
        self.t.Branch('eta_3_tr',    self.eta_3_tr,    'eta_3_tr/F')
        self.t.Branch('m_3',         self.m_3,         'm_3/F')
        self.t.Branch('q_3',         self.q_3,         'q_3/F')
        self.t.Branch('d0_3',        self.d0_3,        'd0_3/F')
        self.t.Branch('dZ_3',        self.dZ_3,        'dZ_3/F')
        self.t.Branch('mt_3',        self.mt_3,        'mt_3/F')
        self.t.Branch('pfmt_3',      self.pfmt_3,      'pfmt_3/F')
        #self.t.Branch('puppimt_3',   self.puppimt_3,   'puppimt_3/F')
        self.t.Branch('iso_3',       self.iso_3,       'iso_3/F')
        self.t.Branch('Electron_mvaFall17V2noIso_WP90_1', self.Electron_mvaFall17V2noIso_WP90_1, 'Electron_mvaFall17V2noIso_WP90_1/F')
        self.t.Branch('Electron_mvaFall17V2noIso_WP90_2', self.Electron_mvaFall17V2noIso_WP90_2, 'Electron_mvaFall17V2noIso_WP90_2/F')
        self.t.Branch('Electron_mvaFall17V2noIso_WP90_3', self.Electron_mvaFall17V2noIso_WP90_3, 'Electron_mvaFall17V2noIso_WP90_3/F')
        self.t.Branch('Electron_mvaFall17V2noIso_WP90_4', self.Electron_mvaFall17V2noIso_WP90_4, 'Electron_mvaFall17V2noIso_WP90_4/F')
        self.t.Branch('gen_match_1', self.gen_match_1, 'gen_match_1/l')
        self.t.Branch('gen_match_2', self.gen_match_2, 'gen_match_2/l')
        self.t.Branch('gen_match_3', self.gen_match_3, 'gen_match_3/l')
        self.t.Branch('tightId_3', self.tightId_3, 'tightId_3/F')
        self.t.Branch('mediumId_3', self.mediumId_3, 'mediumId_3/F')
        self.t.Branch('mediumPromptId_3', self.mediumPromptId_3, 'mediumPromptId_3/F')
        self.t.Branch('looseId_3', self.looseId_3, 'looseId_3/F')
        self.t.Branch('isGlobal_3', self.isGlobal_3, 'isGlobal_3/F')
        self.t.Branch('isTracker_3', self.isTracker_3, 'isTracker_3/F')
        self.t.Branch('ip3d_3', self.ip3d_3, 'ip3d_3/F')


        self.t.Branch('idDecayModeNewDMs_3', self.idDecayModeNewDMs_3, 'idDecayModeNewDMs_3/F')
        self.t.Branch('idDeepTau2017v2p1VSe_3', self.idDeepTau2017v2p1VSe_3, 'idDeepTau2017v2p1VSe_3/F')
        self.t.Branch('idDeepTau2017v2p1VSjet_3', self.idDeepTau2017v2p1VSjet_3, 'idDeepTau2017v2p1VSjet_3/F')
        self.t.Branch('idDeepTau2017v2p1VSmu_3', self.idDeepTau2017v2p1VSmu_3, 'idDeepTau2017v2p1VSmu_3/F')
        self.t.Branch('idMVAnewDM2017v2_3', self.idMVAnewDM2017v2_3, 'idMVAnewDM2017v2_3/F')
        self.t.Branch('rawMVAnewDM2017v2_3', self.rawMVAnewDM2017v2_3, 'rawMVAnewDM2017v2_3/F')

        self.t.Branch('decayMode_3',   self.decayMode_3,   'decayMode_3/I')

        self.t.Branch('pt_4',        self.pt_4,        'pt_4/F')
        self.t.Branch('pt_4_tr',     self.pt_4_tr,        'pt_4_tr/F')
        self.t.Branch('phi_4',       self.phi_4,       'phi_4/F')
        self.t.Branch('phi_4_tr',    self.phi_4_tr,    'phi_4_tr/F')
        self.t.Branch('eta_4',       self.eta_4,       'eta_4/F')
        self.t.Branch('eta_4_tr',    self.eta_4_tr,    'eta_4_tr/F')
        self.t.Branch('m_4',         self.m_4,         'm_4/F')
        self.t.Branch('q_4',         self.q_4,         'q_4/F')
        self.t.Branch('d0_4',        self.d0_4,        'd0_4/F')
        self.t.Branch('dZ_4',        self.dZ_4,        'dZ_4/F')
        self.t.Branch('mt_4',        self.mt_4,        'mt_4/F')
        self.t.Branch('pfmt_4',      self.pfmt_4,      'pfmt_4/F')
        #self.t.Branch('puppimt_4',   self.puppimt_4,   'puppimt_4/F')
        self.t.Branch('iso_4',       self.iso_4,       'iso_4/F')
        self.t.Branch('gen_match_4', self.gen_match_4, 'gen_match_4/l')
        self.t.Branch('tightId_4', self.tightId_4, 'tightId_4/F')
        self.t.Branch('mediumId_4', self.mediumId_4, 'mediumId_4/F')
        self.t.Branch('mediumPromptId_4', self.mediumPromptId_4, 'mediumPromptId_4/F')
        self.t.Branch('looseId_4', self.looseId_4, 'looseId_4/F')
        self.t.Branch('isGlobal_4', self.isGlobal_4, 'isGlobal_4/F')
        self.t.Branch('isTracker_4', self.isTracker_4, 'isTracker_4/F')
        self.t.Branch('ip3d_4', self.ip3d_4, 'ip3d_4/F')


        self.t.Branch('idDecayModeNewDMs_4', self.idDecayModeNewDMs_4, 'idDecayModeNewDMs_4/F')
        self.t.Branch('idDeepTau2017v2p1VSe_4', self.idDeepTau2017v2p1VSe_4, 'idDeepTau2017v2p1VSe_4/F')
        self.t.Branch('idDeepTau2017v2p1VSjet_4', self.idDeepTau2017v2p1VSjet_4, 'idDeepTau2017v2p1VSjet_4/F')
        self.t.Branch('idDeepTau2017v2p1VSmu_4', self.idDeepTau2017v2p1VSmu_4, 'idDeepTau2017v2p1VSmu_4/F')
        self.t.Branch('idMVAnewDM2017v2_4', self.idMVAnewDM2017v2_4, 'idMVAnewDM2017v2_4/F')
        self.t.Branch('rawMVAnewDM2017v2_4', self.rawMVAnewDM2017v2_4, 'rawMVAnewDM2017v2_4/F')

        self.t.Branch('decayMode_4',   self.decayMode_4,   'decayMode_4/I')

        '''
        self.t.Branch('pt_5',        self.pt_5,        'pt_5/F')
        self.t.Branch('phi_5',       self.phi_5,       'phi_5/F')
        self.t.Branch('eta_5',       self.eta_5,       'eta_5/F')
        self.t.Branch('m_5',         self.m_5,         'm_5/F')
        self.t.Branch('q_5',         self.q_5,         'q_5/F')
        self.t.Branch('dZ_5',        self.dZ_5,        'dZ_5/F')
        self.t.Branch('d0_5',        self.d0_5,        'd0_5/F')
        self.t.Branch('gen_match_5', self.gen_match_5, 'gen_match_5/l')
        self.t.Branch('decayMode_5',   self.decayMode_5,   'decayMode_5/I')


        self.t.Branch('idDecayModeNewDMs_5', self.idDecayModeNewDMs_5, 'idDecayModeNewDMs_5/F')
        self.t.Branch('idDeepTau2017v2p1VSe_5', self.idDeepTau2017v2p1VSe_5, 'idDeepTau2017v2p1VSe_5/F')
        self.t.Branch('idDeepTau2017v2p1VSjet_5', self.idDeepTau2017v2p1VSjet_5, 'idDeepTau2017v2p1VSjet_5/F')
        self.t.Branch('idDeepTau2017v2p1VSmu_5', self.idDeepTau2017v2p1VSmu_5, 'idDeepTau2017v2p1VSmu_5/F')
        self.t.Branch('idMVAnewDM2017v2_5', self.idMVAnewDM2017v2_5, 'idMVAnewDM2017v2_5/F')
        self.t.Branch('rawMVAnewDM2017v2_5', self.rawMVAnewDM2017v2_5, 'rawMVAnewDM2017v2_5/F')
        '''



        # di-tau variables
        self.t.Branch('pt_tt', self.pt_tt, 'pt_tt/F')
        self.t.Branch('mt_tot', self.mt_tot, 'mt_tot/F')
        self.t.Branch('m_vis', self.m_vis, 'm_vis/F')
        self.t.Branch('m_sv', self.m_sv, 'm_sv/F')
        self.t.Branch('mt_sv', self.mt_sv, 'mt_sv/F') 
        self.t.Branch('H_DR', self.H_DR, 'H_DR/F')
        self.t.Branch('AMass', self.AMass, 'AMass/F')

        # di-lepton variables. 
        self.t.Branch('H_LT',         self.H_LT,         'H_LT/F')   
        self.t.Branch('dRl1H',         self.dRl1H,         'dRl1H/F')   
        self.t.Branch('dRl2H',         self.dRl2H,         'dRl2H/F')   
        self.t.Branch('dRlH',         self.dRlH,         'dRlH/F')   
        self.t.Branch('dPhil1H',         self.dPhil1H,         'dPhil1H/F')   
        self.t.Branch('dPhil2H',         self.dPhil2H,         'dPhil2H/F')   
        self.t.Branch('dPhilH',         self.dPhilH,         'dPhilH/F')   

        self.t.Branch('mll',         self.mll,         'mll/F')   
        self.t.Branch('mll2',         self.mll2,         'mll2/F')   
        self.t.Branch('Z_Pt',       self.Z_Pt,       'Z_Pt/F')   
        self.t.Branch('Z_DR',       self.Z_DR,       'Z_DR/F')   
        self.t.Branch('Z_SS',       self.Z_SS,       'Z_SS/F')   
        self.t.Branch('pt_1',        self.pt_1,        'pt_1/F')
        self.t.Branch('m_1_tr',     self.m_1_tr,     'm_1_tr/F')
        self.t.Branch('pt_1_tr',     self.pt_1_tr,     'pt_1_tr/F')
        self.t.Branch('phi_1',       self.phi_1,       'phi_1/F')  
        self.t.Branch('phi_1_tr',    self.phi_1_tr,    'phi_1_tr/F')
        self.t.Branch('eta_1',       self.eta_1,       'eta_1/F')    
        self.t.Branch('eta_1_tr',    self.eta_1_tr,    'eta_1_tr/F')
        self.t.Branch('pt_2',        self.pt_2,        'pt_2/F')      
        self.t.Branch('m_2_tr',     self.m_2_tr,     'm_2_tr/F')
        self.t.Branch('pt_2_tr',     self.pt_2_tr,     'pt_2_tr/F')
        self.t.Branch('phi_2',       self.phi_2,       'phi_2/F')    
        self.t.Branch('phi_2_tr',    self.phi_2_tr,    'phi_2_tr/F')
        self.t.Branch('eta_2',       self.eta_2,       'eta_2/F')      
        self.t.Branch('eta_2_tr',    self.eta_2_tr,    'eta_2_tr/F')
        self.t.Branch('iso_1',       self.iso_1,       'iso_1/F')
        self.t.Branch('iso_2',       self.iso_2,       'iso_2/F')
        self.t.Branch('q_1',       self.q_1,       'q_1/F')
        self.t.Branch('q_2',       self.q_2,       'q_2/F')
        self.t.Branch('L1PreFiringWeight_Nom',        self.L1PreFiringWeight_Nom,        'L1PreFiringWeight_Nom/F')
        self.t.Branch('L1PreFiringWeight_Up',        self.L1PreFiringWeight_Up,        'L1PreFiringWeight_Up/F')
        self.t.Branch('L1PreFiringWeight_Down',        self.L1PreFiringWeight_Down,        'L1PreFiringWeight_Down/F')
        self.t.Branch('d0_1',        self.d0_1,        'd0_1/F')
        self.t.Branch('dZ_1',        self.dZ_1,        'dZ_1/F')
        self.t.Branch('d0_2',        self.d0_2,        'd0_2/F')
        self.t.Branch('dZ_2',        self.dZ_2,        'dZ_2/F')
        self.t.Branch('Muon_Id_1',       self.Muon_Id_1,       'Muon_Id_1/F')
        self.t.Branch('Muon_Id_2',       self.Muon_Id_2,       'Muon_Id_2/F')
        self.t.Branch('isGlobal_1',       self.isGlobal_1,       'isGlobal_1/F')
        self.t.Branch('isGlobal_2',       self.isGlobal_2,       'isGlobal_2/F')
        self.t.Branch('isTracker_1',       self.isTracker_1,       'isTracker_1/F')
        self.t.Branch('isTracker_2',       self.isTracker_2,       'isTracker_2/F')
        self.t.Branch('tightId_1', self.tightId_1, 'tightId_1/F')
        self.t.Branch('mediumId_1', self.mediumId_1, 'mediumId_1/F')
        self.t.Branch('mediumPromptId_1', self.mediumPromptId_1, 'mediumPromptId_1/F')
        self.t.Branch('looseId_1', self.looseId_1, 'looseId_1/F')
        self.t.Branch('tightId_2', self.tightId_2, 'tightId_2/F')
        self.t.Branch('mediumId_2', self.mediumId_2, 'mediumId_2/F')
        self.t.Branch('mediumPromptId_2', self.mediumPromptId_2, 'mediumPromptId_2/F')
        self.t.Branch('looseId_2', self.looseId_2, 'looseId_2/F')

        #systematics
        self.t.Branch('MET_pt_UnclUp', self.MET_pt_UnclUp, 'MET_pt_UnclUp/F')
        self.t.Branch('MET_phi_UnclUp', self.MET_phi_UnclUp, 'MET_phi_UnclUp/F')
        self.t.Branch('MET_pt_UnclDown', self.MET_pt_UnclDown, 'MET_pt_UnclDown/F')
        self.t.Branch('MET_phi_UnclDown', self.MET_phi_UnclDown, 'MET_phi_UnclDown/F')
        self.t.Branch('met_UnclX', self.met_UnclX, 'met_UnclX/F')
        self.t.Branch('met_UnclY', self.met_UnclY, 'met_UnclY/F')
        self.t.Branch('MET_T1Smear_pt', self.MET_T1Smear_pt, 'MET_T1Smear_pt/F')
        self.t.Branch('MET_T1Smear_phi', self.MET_T1Smear_phi, 'MET_T1Smear_phi/F')
        
        # MET variables
        self.t.Branch('met', self.met, 'met/F')
        self.t.Branch('metphi', self.metphi, 'metphi/F')
        self.t.Branch('metNoCor', self.metNoCor, 'metNoCor/F')
        self.t.Branch('metphiNoCor', self.metphiNoCor, 'metphiNoCor/F')
        self.t.Branch('metNoTauES', self.metNoTauES, 'metNoTauES/F')
        self.t.Branch('metphiNoTauES', self.metphiNoTauES, 'metphiNoTauES/F')
        #self.t.Branch('puppimet', self.puppimet, 'puppimet/F')
        #self.t.Branch('puppimetphi', self.puppimetphi, 'puppimetphi/F')
        self.t.Branch('metcov00', self.metcov00, 'metcov00/F')
        self.t.Branch('metcov01', self.metcov01, 'metcov01/F')
        self.t.Branch('metcov10', self.metcov10, 'metcov10/F')
        self.t.Branch('metcov11', self.metcov11, 'metcov11/F')

        # trigger sf
        self.t.Branch('isTrig_2',  self.isTrig_2, 'isTrig_2/F' )
        self.t.Branch('isTrig_1',  self.isTrig_1, 'isTrig_1/F' )
        self.t.Branch('isDoubleTrig',  self.isDoubleTrig, 'isDoubleTrig/F' )


        # jet variables
        #self.t.Branch('njetsold', self.njetsold, 'njetsold[8]/F') 
        #self.t.Branch('nbtagold', self.nbtagold, 'nbtagold[8]/F')
        self.t.Branch('njets', self.njets, 'njets/F')
        self.t.Branch('nbtagL', self.nbtagL, 'nbtagL/F')
        self.t.Branch('nbtagM', self.nbtagM, 'nbtagM/F')
        self.t.Branch('nbtagT', self.nbtagT, 'nbtagT/F')
        self.t.Branch('HTXS_Higgs_cat', self.HTXS_Higgs_cat, 'HTXS_Higgs_cat/l')
        self.t.Branch('HTXS_Higgs_pt', self.HTXS_Higgs_pt, 'HTXS_Higgs_pt/F')


        self.t.Branch('jflavour',     self.jflavour,     'jflavour[12]/F' )
        self.t.Branch('jeta',     self.jeta,     'jeta[12]/F' )
        self.t.Branch('jpt',     self.jpt,     'jpt[12]/F' )
        self.t.Branch('btagDeep', self.btagDeep, 'btagDeep[12]/F')

        '''
        self.t.Branch('jpt_1',     self.jpt_1,     'jpt_1/F' )
        self.t.Branch('jpt_2',     self.jpt_2,     'jpt_2/F' )

        self.t.Branch('jpt_1_tr',  self.jpt_1_tr,  'jpt_1_tr/F' )
        self.t.Branch('jeta_1',    self.jeta_1,    'jeta_1/F' ) 
        self.t.Branch('jeta_1_tr', self.jeta_1_tr, 'jeta_1_tr/F' )
        self.t.Branch('jphi_1',    self.jphi_1,    'jphi_1/F' )
        self.t.Branch('jphi_1_tr', self.jphi_1_tr, 'jphi_1_tr/F' )
        self.t.Branch('jcsv_1',    self.jcsv_1,    'jcsv_1/F' )
        self.t.Branch('jcsvfv_1', self.jcsvfv_1, 'jcsvfv_1/F' )
        self.t.Branch('jpt_2',     self.jpt_2,     'jpt_2/F' )
        self.t.Branch('jpt_2_tr',  self.jpt_2_tr,  'jpt_2_tr/F' )
        self.t.Branch('jeta_2',    self.jeta_2,    'jeta_2/F' ) 
        self.t.Branch('jeta_2_tr', self.jeta_2_tr, 'jeta_2_tr/F' )
        self.t.Branch('jphi_2',    self.jphi_2,    'jphi_2/F' )
        self.t.Branch('jphi_2_tr', self.jphi_2_tr, 'jphi_2_tr/F' )
        self.t.Branch('jcsv_2',    self.jcsv_2,    'jcsv_2/F' )
        self.t.Branch('jcsvfv_2', self.jcsvfv_2, 'jcsvfv_2/F' )

        self.t.Branch('bpt_1',     self.bpt_1,     'bpt_1/F' )
        self.t.Branch('bpt_1_tr',  self.bpt_1_tr,  'bpt_1_tr/F' )
        self.t.Branch('beta_1',    self.beta_1,    'beta_1/F' ) 
        self.t.Branch('beta_1_tr', self.beta_1_tr, 'beta_1_tr/F' )
        self.t.Branch('bphi_1',    self.bphi_1,    'bphi_1/F' )
        self.t.Branch('bphi_1_tr', self.bphi_1_tr, 'bphi_1_tr/F' )
        self.t.Branch('bcsv_1',    self.bcsv_1,    'bcsv_1/F' )
        self.t.Branch('bcsvfv_1', self.bcsvfv_1, 'bcsvfv_1/F' )
        self.t.Branch('bpt_2',     self.bpt_2,     'bpt_2[F' )
        self.t.Branch('bpt_2_tr',  self.bpt_2_tr,  'bpt_2_tr/F' )
        self.t.Branch('beta_2',    self.beta_2,    'beta_2/F' )
        self.t.Branch('beta_2_tr', self.beta_2_tr, 'beta_2_tr/F' )
        self.t.Branch('bphi_2',    self.bphi_2,    'bphi_2/F' )
        self.t.Branch('bphi_2_tr', self.bphi_2_tr, 'bphi_2_tr/F' )
        self.t.Branch('bcsv_2',    self.bcsv_2,    'bcsv_2/F' )
        self.t.Branch('bcsvfv_2', self.bcsvfv_2, 'bcsvfv_2/F' )
        '''
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
		    self.t.Branch('btagDeep{0:s}'.format(v), self.list_of_arraysJetsNbtagDeep[i], 'btagDeep{0:s}[12]/F'.format(v))



        #self.MET_pt_jesEC2Up  = array('f',[0])
        #self.t.Branch('MET_pt_jesEC2Up', self.MET_pt_jesEC2Up, 'MET_pt_jesEC2Up/F' )
        self.tN=[]

	#self.t.SetBranchStatus("*Up",0)
	#self.t.SetBranchStatus("*Down",0)
	self.t.SetBranchStatus("GenPart*",0)
	self.t.SetBranchStatus("*_tr*",0)
	self.t.SetBranchStatus("*LHE*",0)
	#self.t.SetBranchStatus("*LHEScaleWeight",1)
	self.t.SetBranchStatus("dR*",0)
	self.t.SetBranchStatus("dPhi*",0)
	self.t.SetBranchStatus("Z_*",0)
	self.t.SetBranchStatus("*ip3d*",0)
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

	#self.t.SetBranchStatus("*Up",1)
	#self.t.SetBranchStatus("*Down",1)
	self.t.SetBranchStatus("GenPart*",1)
	self.t.SetBranchStatus("*_tr*",1)
	#self.t.SetBranchStatus("*LHE*",1)
	self.t.SetBranchStatus("*LHEScaleWeight*",1)
	self.t.SetBranchStatus("dR*",1)
	self.t.SetBranchStatus("dPhi*",1)
	self.t.SetBranchStatus("Z_*",1)
	self.t.SetBranchStatus("*ip3d*",1)
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
    
    def Fill(self, entry, SVFit, cat, jt1, jt2, LepP, LepM, lepList, isMC, era, doUncertainties=False ,  met_pt=-99, met_phi=-99, systIndex=0, tMass=[], tPt=[], eMass=[], ePt=[], mMass=[], mPt=[], proc="EOY") : 
    #def Fill(self, entry, SVFit, cat, jt1, jt2, LepP, LepM, lepList, isMC, era, doUncertainties=False ,  met_pt=-99, met_phi=-99, systIndex=0) : 
    #def Fill(self, entry, SVFit, cat, jt1, jt2, LepP, LepM, lepList, isMC, era,  doUncertainties=False , sysVariations=[]) :
        ''' - jt1 and jt2 point to the selected tau candidates according to the table below.
            - if e.g., channel = 'et', the jt1 points to the electron list and jt2 points to the tau list.
            - LepP and LepM are TLorentz vectors for the positive and negative members of the dilepton pair
        '''
        SystIndex = int(systIndex)

        
        #if SystIndex >0 : doUncertainties=False

        #channel_ll = 'mm' or 'ee'
        channel_ll = cat[:-2]
	channel = cat[-2:]

        if SystIndex ==0 : 

	    is_trig_1, is_trig_2, is_Dtrig_1 = 0., 0., 0.
	    TrigListLep = []
	    TrigListTau = []
	    hltListLep  = []
	    hltListLepSubL  = []

	    TrigListLep, hltListLep, hltListLepSubL  = GF.findSingleLeptTrigger(lepList, entry, channel_ll, era)

	    TrigListLep = list(dict.fromkeys(TrigListLep))
	    #if len(hltListLep) > 0 or len(hltListLepSubL)>0 :     print GF.printEvent(entry), SystIndex

	    #TrigListLepD, hltListLepD  = GF.findDoubleLeptTrigger(lepList, entry, channel_ll, era)

	    #TrigListLepD = list(dict.fromkeys(TrigListLepD))

	    #if len(TrigListLepD) > 0 : print TrigListLepD, hltListLepD, TrigListLep, hltListLep
	    #if len(TrigListLepD) == 2 : 
	    #    if lepList[0] == TrigListLepD[0] :
	    #        is_Dtrig_1 = 1 #that means that the leading lepton 
	    #    else : 
	    #        is_Dtrig_1 = -1


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
	    self.evt[0]  = entry.event
	    self.iso_1[0]  = -99
	    self.iso_2[0]  = -99
	    self.q_1[0]  = -99
	    self.q_2[0]  = -99
	    self.isGlobal_1[0]  = -99
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

	    self.tightId_2[0]       = -1 
	    self.mediumId_2[0]       = -1 
	    self.mediumPromptId_2[0]   = -1
	    self.looseId_2[0]       = -1
	    self.isGlobal_2[0]      = -1
	    self.isTracker_2[0]     = -1

	    self.decayMode_3[0]        = -1
	    self.idDecayModeNewDMs_3[0]= -1
	    self.idDeepTau2017v2p1VSe_3[0] = -1
	    self.idDeepTau2017v2p1VSjet_3[0] = -1
	    self.idDeepTau2017v2p1VSmu_3[0] = -1
	    self.idMVAnewDM2017v2_3[0] = -1
	    self.rawMVAnewDM2017v2_3[0] = -1
	    self.mediumId_3[0]       = -1 
	    self.mediumPromptId_3[0]   = -1
	    self.looseId_3[0]       = -1
	    self.isGlobal_3[0]      = -1
	    self.isTracker_3[0]     = -1
	    self.ip3d_3[0]          = -1

	    self.decayMode_4[0]      = -1
	    self.idDecayModeNewDMs_4[0] = -1
	    self.idDeepTau2017v2p1VSe_4[0] = -1
	    self.idDeepTau2017v2p1VSjet_4[0] = -1
	    self.idDeepTau2017v2p1VSmu_4[0] = -1
	    self.idMVAnewDM2017v2_4[0] = -1
	    self.rawMVAnewDM2017v2_4[0] = -1
	    self.mediumId_4[0]       = -1 
	    self.mediumPromptId_4[0]   = -1
	    self.looseId_4[0]       = -1
	    self.isGlobal_4[0]      = -1
	    self.isTracker_4[0]     = -1
	    self.ip3d_4[0]          = -1
	    self.GenPart_statusFlags_1[0]    = -1
	    self.GenPart_status_1[0]    = -1
	    self.GenPart_statusFlags_2[0]    = -1
	    self.GenPart_status_2[0]    = -1
	    self.GenPart_statusFlags_3[0]    = -1
	    self.GenPart_status_3[0]    = -1
	    self.GenPart_statusFlags_4[0]    = -1
	    self.GenPart_status_4[0]    = -1
	    self.gen_match_1[0] = -1
	    self.gen_match_2[0] = -1
	    self.gen_match_3[0] = -1
	    self.gen_match_4[0] = -1
	    #self.gen_match_5[0] = -1


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
        '''
        goodElectronList = tauFunDCH.makeGoodElectronList(entry)
        goodMuonList = tauFunDCH.makeGoodMuonList(entry)
        
        self.nGoodElectron[0] = len(goodElectronList)
        self.nGoodMuon[0]     = len(goodMuonList)
        # pack trigger bits into integer word
        '''

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
	    try : bits.append(e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ)
	    except AttributeError : bits.append(False) 
	    try : bits.append(e.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_Mass8)
	    except AttributeError : bits.append(False) 

	    self.muonTriggerWord[0] = 0
	    for i, bit in enumerate(bits) :
		if bit : self.muonTriggerWord[0] += 2**i

        #neede for all systematics as jt1/jt2 may change per systematic
	if jt1 > -1 and jt2 > -1 : self.cat[0]  = tauFunDCH.catToNumber(cat)
        if jt1>-1 or jt2 >-1 :

	    tau1, tau2 = TLorentzVector(), TLorentzVector()

	    # Fill variables for Leg3, where 3->tau(ele) and 4->tau(had)
	    if channel == 'et' :
		self.pt_3[0]   = entry.Electron_pt[jt1]
		self.phi_3[0]  = entry.Electron_phi[jt1]
		self.eta_3[0]  = entry.Electron_eta[jt1]
		self.m_3[0]    = entry.Electron_mass[jt1]
		self.q_3[0]    = entry.Electron_charge[jt1]
		self.d0_3[0]   = entry.Electron_dxy[jt1]
		self.dZ_3[0]   = entry.Electron_dz[jt1]
		self.iso_3[0]  = entry.Electron_pfRelIso03_all[jt1]
		self.Electron_mvaFall17V2noIso_WP90_3[0]  = entry.Electron_mvaFall17V2noIso_WP90[jt1]

		if SystIndex ==0 and  isMC: 
		    self.pt_uncor_3[0] = ePt[jt1]
		    self.m_uncor_3[0] = eMass[jt1]
		    self.pt_uncor_4[0] = tPt[jt2]
		    self.m_uncor_4[0] = tMass[jt2]
		
		# Fill genMatch variables for tau(ele)
		if isMC:
		    idx_genEle = entry.Electron_genPartIdx[jt1]

		    # if idx_genMu = -1, no match was found
		    if idx_genEle >= 0:
			idx_genEle_mom      = entry.GenPart_genPartIdxMother[idx_genEle]
			self.pt_3_tr[0]     = entry.GenPart_pt[idx_genEle]
			self.phi_3_tr[0]    = entry.GenPart_phi[idx_genEle]
			self.eta_3_tr[0]    = entry.GenPart_eta[idx_genEle]
			self.GenPart_statusFlags_3[0]    = entry.GenPart_statusFlags[idx_genEle]
			self.GenPart_status_3[0]    = entry.GenPart_status[idx_genEle]

		    try: self.gen_match_3[0] = ord(entry.Electron_genPartFlav[jt1])
		    except AttributeError: self.gen_match_3[0] = -1
		
                #print '---------------------------->', self.pt_3[0], ePt[jt1], entry.Tau_pt[jt2], tPt[jt2] , jt1, jt2, cat, entry.event, SystIndex
		tau1.SetPtEtaPhiM(entry.Electron_pt[jt1],entry.Electron_eta[jt1], entry.Electron_phi[jt1], self.tauMass)
                tmass= self.tauMass
                if entry.Tau_decayMode[jt2] == 0 : tmass= 0.13960
		tau2.SetPtEtaPhiM(entry.Tau_pt[jt2],entry.Tau_eta[jt2],entry.Tau_phi[jt2],tmass)
		
		tauListE=[jt1]
	       

	    # Fill variables for Leg3 and Leg4, where 3->tau(ele) and 4->tau(mu)
	    elif channel == 'em' :
		self.pt_3[0]   = entry.Electron_pt[jt1]
		self.phi_3[0]  = entry.Electron_phi[jt1]
		self.eta_3[0]  = entry.Electron_eta[jt1]
		self.m_3[0]    = entry.Electron_mass[jt1]
		self.q_3[0]    = entry.Electron_charge[jt1]
		self.d0_3[0]   = entry.Electron_dxy[jt1]
		self.dZ_3[0]   = entry.Electron_dz[jt1]
		self.iso_3[0]  = entry.Electron_pfRelIso03_all[jt1]
		self.Electron_mvaFall17V2noIso_WP90_3[0]  = entry.Electron_mvaFall17V2noIso_WP90[jt1]

		if SystIndex ==0 and  isMC: 
		    self.pt_uncor_3[0] = ePt[jt1]
		    self.m_uncor_3[0] = eMass[jt1]
		    self.pt_uncor_4[0] = mPt[jt2]
		    self.m_uncor_4[0] = mMass[jt2]
		
		if isMC:
		    try : self.gen_match_3[0] = ord(entry.Electron_genPartFlav[jt1])
		    except AttributeError : self.gen_match_3[0] = -1
		
		tau1.SetPtEtaPhiM(entry.Electron_pt[jt1], entry.Electron_eta[jt1], entry.Electron_phi[jt1], self.tauMass)
													    #???
		# fill genMatch for tau(ele)
		if isMC:
		    idx_genEle = entry.Electron_genPartIdx[jt1]

		    # if idx_genEle = -1, no match was found
		    if idx_genEle >= 0:
			idx_genEle_mom      = entry.GenPart_genPartIdxMother[idx_genEle]
			self.pt_3_tr[0]     = entry.GenPart_pt[idx_genEle]
			self.phi_3_tr[0]    = entry.GenPart_phi[idx_genEle]
			self.eta_3_tr[0]    = entry.GenPart_eta[idx_genEle]
			self.GenPart_statusFlags_3[0]    = entry.GenPart_statusFlags[idx_genEle]
			self.GenPart_status_3[0]    = entry.GenPart_status[idx_genEle]

		self.pt_4[0]     = entry.Muon_pt[jt2]
		self.phi_4[0]    = entry.Muon_phi[jt2]
		self.eta_4[0]    = entry.Muon_eta[jt2]
		self.m_4[0]      = entry.Muon_mass[jt2]
		self.q_4[0]      = entry.Muon_charge[jt2]
		self.d0_4[0]     = entry.Muon_dxy[jt2]
		self.dZ_4[0]     = entry.Muon_dz[jt2]
		self.iso_4[0]    = entry.Muon_pfRelIso04_all[jt2]
		self.tightId_4[0]      = entry.Muon_tightId[jt2]
		self.mediumId_4[0]      = entry.Muon_mediumId[jt2]
		self.mediumPromptId_4[0]   = entry.Muon_mediumPromptId[jt2]
		self.looseId_4[0]       = entry.Muon_looseId[jt2]
		self.isGlobal_4[0]      = entry.Muon_isGlobal[jt2]
		self.isTracker_4[0]     = entry.Muon_isTracker[jt2]
		self.ip3d_4[0]       = entry.Muon_ip3d[jt2]
		if isMC:
		    try : self.gen_match_4[0] = ord(entry.Muon_genPartFlav[jt2]) 
		    except AttributeError : self.gen_match_4[0] = -1
		
		tau2.SetPtEtaPhiM(entry.Muon_pt[jt2], entry.Muon_eta[jt2], entry.Muon_phi[jt2], self.tauMass) 

		# fill genMatch for tau(mu)
		if isMC:
		    idx_genMu = entry.Muon_genPartIdx[jt2]
		    
		    # if idx_genMu = -1, no match was found
		    if idx_genMu >= 0:
			idx_genMu_mom       = entry.GenPart_genPartIdxMother[idx_genMu]
			self.pt_4_tr[0]     = entry.GenPart_pt[idx_genMu]
			self.phi_4_tr[0]    = entry.GenPart_phi[idx_genMu]
			self.eta_4_tr[0]    = entry.GenPart_eta[idx_genMu]
			self.GenPart_statusFlags_4[0]    = entry.GenPart_statusFlags[idx_genMu]
			self.GenPart_status_4[0]    = entry.GenPart_status[idx_genMu]


	    # Fill variables for Leg3, where 3->tau(mu) and 4->tau(had)
	    elif channel == 'mt' :
		self.pt_3[0]     = entry.Muon_pt[jt1]
		self.phi_3[0]    = entry.Muon_phi[jt1]
		self.eta_3[0]    = entry.Muon_eta[jt1]
		self.m_3[0]      = entry.Muon_mass[jt1]
		self.q_3[0]      = entry.Muon_charge[jt1]
		self.d0_3[0]     = entry.Muon_dxy[jt1]
		self.dZ_3[0]     = entry.Muon_dz[jt1]
		self.iso_3[0]    = entry.Muon_pfRelIso04_all[jt1]
		self.tightId_3[0]      = entry.Muon_tightId[jt1]
		self.mediumId_3[0]       = entry.Muon_mediumId[jt1]
		self.mediumPromptId_3[0]   = entry.Muon_mediumPromptId[jt1]
		self.looseId_3[0]       = entry.Muon_looseId[jt1]
		self.isGlobal_3[0]      = entry.Muon_isGlobal[jt1]
		self.isTracker_3[0]     = entry.Muon_isTracker[jt1]
		self.ip3d_3[0]       = entry.Muon_ip3d[jt1]
		if SystIndex ==0 and isMC : 
		    self.pt_uncor_3[0] = mPt[jt1]
		    self.m_uncor_3[0] = mMass[jt1]
		    self.pt_uncor_4[0] = tPt[jt2]
		    self.m_uncor_4[0] = tMass[jt2]
		
		if isMC:
		    try : self.gen_match_3[0] = ord(entry.Muon_genPartFlav[jt1])
		    except AttributeError : self.gen_match_3[0] = -1
		
		tau1.SetPtEtaPhiM(entry.Muon_pt[jt1], entry.Muon_eta[jt1], entry.Muon_phi[jt1], self.tauMass)
                tmass= self.tauMass
                if entry.Tau_decayMode[jt2] == 0 : tmass= 0.13960
		tau2.SetPtEtaPhiM(entry.Tau_pt[jt2],  entry.Tau_eta[jt2],  entry.Tau_phi[jt2],  tmass) 
		
		# fill genMatch for tau(mu)
		if isMC:
		    idx_genMu = entry.Muon_genPartIdx[jt1]
		    
		    # if idx_genMu = -1, no match was found
		    if idx_genMu >= 0:
			idx_genMu_mom       = entry.GenPart_genPartIdxMother[idx_genMu]
			self.pt_3_tr[0]     = entry.GenPart_pt[idx_genMu]
			self.phi_3_tr[0]    = entry.GenPart_phi[idx_genMu]
			self.eta_3_tr[0]    = entry.GenPart_eta[idx_genMu]
			self.GenPart_statusFlags_3[0]    = entry.GenPart_statusFlags[idx_genMu]
			self.GenPart_status_3[0]    = entry.GenPart_status[idx_genMu]
	    
	    # Fill variables for Leg3 and Leg4, where 3->tau(had) and 4->tau(had)
	    elif channel == 'tt' :
		self.pt_3[0]     = entry.Tau_pt[jt1]
		self.phi_3[0]    = entry.Tau_phi[jt1]
		self.eta_3[0]    = entry.Tau_eta[jt1]
		self.m_3[0]      = entry.Tau_mass[jt1]
		self.q_3[0]      = entry.Tau_charge[jt1]
		self.d0_3[0]     = entry.Tau_dxy[jt1]
		self.dZ_3[0]     = entry.Tau_dz[jt1]
		if SystIndex ==0 and isMC: 
		    self.pt_uncor_3[0] = tPt[jt1]
		    self.m_uncor_3[0] = tMass[jt1]
		    self.pt_uncor_4[0] = tPt[jt2]
		    self.m_uncor_4[0] = tMass[jt2]
                #print '=========================================--------------------------------> inside', entry.Tau_mass[jt1] , entry.Tau_pt[jt1], jt1, int(entry.Tau_decayMode[jt1])

		self.idDecayModeNewDMs_3[0] = entry.Tau_idDecayModeNewDMs[jt1]
		self.idDeepTau2017v2p1VSe_3[0] = ord(entry.Tau_idDeepTau2017v2p1VSe[jt1])
		self.idDeepTau2017v2p1VSjet_3[0] = ord(entry.Tau_idDeepTau2017v2p1VSjet[jt1])
		self.idDeepTau2017v2p1VSmu_3[0] = ord(entry.Tau_idDeepTau2017v2p1VSmu[jt1])
		self.idMVAnewDM2017v2_3[0] = ord(entry.Tau_idMVAnewDM2017v2[jt1])
		self.rawMVAnewDM2017v2_3[0] = entry.Tau_rawMVAnewDM2017v2[jt1]

	
		# genMatch the hadronic tau candidate
		if isMC:
		    idx_t1_gen = GF.genMatchTau(entry, jt1, 'had')
		    if idx_t1_gen >= 0:
			self.pt_3_tr[0]  = entry.GenVisTau_pt[idx_t1_gen]
			self.phi_3_tr[0] = entry.GenVisTau_phi[idx_t1_gen]
			self.eta_3_tr[0] = entry.GenVisTau_eta[idx_t1_gen]
			self.GenPart_statusFlags_3[0]    = entry.GenPart_statusFlags[idx_t1_gen]
			self.GenPart_status_3[0]    = entry.GenPart_status[idx_t1_gen]
		    else:
			self.pt_3_tr[0]  = 1.2*entry.Tau_pt[jt1]
			self.phi_3_tr[0] = 1.2*entry.Tau_phi[jt1]
			self.eta_3_tr[0] = 1.2*entry.Tau_eta[jt1]

		    try : self.gen_match_3[0] = ord(entry.Tau_genPartFlav[jt1])
		    except AttributeError : self.gen_match_3[0] = -1

		try : self.decayMode_3[0] = int(entry.Tau_decayMode[jt1])
		except AttributeError : self.decayMode_3[0] = -1

                tmass= self.tauMass
                if entry.Tau_decayMode[jt1] == 0 : tmass= 0.13960
		tau1.SetPtEtaPhiM(entry.Tau_pt[jt1], entry.Tau_eta[jt1], entry.Tau_phi[jt1], tmass)
                tmass= self.tauMass
                if entry.Tau_decayMode[jt2] == 0 : tmass= 0.13960
		tau2.SetPtEtaPhiM(entry.Tau_pt[jt2], entry.Tau_eta[jt2], entry.Tau_phi[jt2], tmass)
		
	    else :
		print("Invalid channel={0:s} in outTuple(). Exiting.".format(channel))
		exit()
		
	    self.mt_3[0]      = self.get_mt('MVAMet',   entry,tau1)
	    self.pfmt_3[0]    = self.get_mt('PFMet',    entry,tau1)
	    #self.puppimt_3[0] = self.get_mt('PUPPIMet', entry,tau1)

	    
	    # Fill variables for Leg4, where 4->tau(had)
	    if channel != 'em':
		self.pt_4[0]  = entry.Tau_pt[jt2]
		self.phi_4[0] = entry.Tau_phi[jt2]
		self.eta_4[0] = entry.Tau_eta[jt2]
		self.m_4[0]   = entry.Tau_mass[jt2]
		self.q_4[0]   = entry.Tau_charge[jt2]
		self.d0_4[0]  = entry.Tau_dxy[jt2]
		self.dZ_4[0]  = entry.Tau_dz[jt2]

		self.idDecayModeNewDMs_4[0] = entry.Tau_idDecayModeNewDMs[jt2]
		self.idDeepTau2017v2p1VSe_4[0] = ord(entry.Tau_idDeepTau2017v2p1VSe[jt2])
		self.idDeepTau2017v2p1VSjet_4[0] = ord(entry.Tau_idDeepTau2017v2p1VSjet[jt2])
		self.idDeepTau2017v2p1VSmu_4[0] = ord(entry.Tau_idDeepTau2017v2p1VSmu[jt2])
		self.idMVAnewDM2017v2_4[0] = ord(entry.Tau_idMVAnewDM2017v2[jt2])
		self.rawMVAnewDM2017v2_4[0] = entry.Tau_rawMVAnewDM2017v2[jt2]
		
		phi, pt = entry.Tau_phi[jt2], entry.Tau_pt[jt2]
		
		self.mt_4[0]      = self.get_mt('MVAMet',   entry, tau2) 
		self.pfmt_4[0]    = self.get_mt('PFMet',    entry, tau2)
		#self.puppimt_4[0] = self.get_mt('PUPPIMet', entry, tau2) 


		# genMatch the hadronic tau candidate
		if isMC:
		    idx_t2_gen = GF.genMatchTau(entry, jt2, 'had')
		    if idx_t2_gen >= 0:
			self.pt_4_tr[0]  = entry.GenVisTau_pt[idx_t2_gen]
			self.phi_4_tr[0] = entry.GenVisTau_phi[idx_t2_gen]
			self.eta_4_tr[0] = entry.GenVisTau_eta[idx_t2_gen]
			self.GenPart_statusFlags_4[0]    = entry.GenPart_statusFlags[idx_t2_gen]
			self.GenPart_status_4[0]    = entry.GenPart_status[idx_t2_gen]
		    else:
			self.pt_4_tr[0]  = 1.2*entry.Tau_pt[jt2]
			self.phi_4_tr[0] = 1.2*entry.Tau_phi[jt2]
			self.eta_4_tr[0] = 1.2*entry.Tau_eta[jt2]

		    try : self.gen_match_4[0] = ord(entry.Tau_genPartFlav[jt2])
		    except AttributeError: self.gen_match_4[0] = -1

		try : self.decayMode_4[0] = int(entry.Tau_decayMode[jt2])
		except AttributeError: self.decayMode_4[0] = -1


	    # di-tau variables
	    self.pt_tt[0]  = self.getPt_tt( entry, tau1, tau2)
	    self.H_DR[0] = self.getDR(entry,tau1,tau2)
	    self.mt_tot[0] = self.getMt_tot(entry, tau1, tau2)
	    self.m_vis[0]  = self.getM_vis( entry, tau1, tau2)
		
	    if SVFit :
		fastMTTmass, fastMTTtransverseMass = self.runSVFit(entry, channel, jt1, jt2, tau1, tau2,met_pt,met_phi) 
	    else :
		fastMTTmass, fastMTTtransverseMass = -999., -999.
		
	    self.m_sv[0] = fastMTTmass 
	    self.mt_sv[0] = fastMTTtransverseMass  


        # Sort the di-lepton system by Pt
        Lep1, Lep2 = TLorentzVector(), TLorentzVector()
        if (LepP.Pt() > LepM.Pt()): 
            Lep1 = LepP
            Lep2 = LepM
        else:
            Lep1 = LepM
            Lep2 = LepP


        # di-lepton variables.   _p and _m refer to plus and minus charge
        if jt1>-1 and jt2>-1 : self.AMass[0]       = (Lep1 + Lep2 + tau1 + tau2).M() 
        self.mll[0]       = (Lep1 + Lep2).M()
        '''
        self.Z_DR[0]       = self.getDR(entry,Lep1,Lep2)
       
        self.H_LT[0]       = Lep1.Pt() + Lep2.Pt()
        self.dRl1H[0]  = self.getDR(entry,Lep1,tau1+tau2)
        self.dRl2H[0]  = self.getDR(entry,Lep2,tau1+tau2)
        self.dRlH[0]  = self.getDR(entry,Lep1+Lep2,tau1+tau2)

        self.dPhil1H[0]  = self.getdPhi(entry,Lep1,tau1+tau2)
        self.dPhil2H[0]  = self.getdPhi(entry,Lep2,tau1+tau2)
        self.dPhilH[0]  = self.getdPhi(entry,Lep1+Lep2,tau1+tau2)
        '''
        self.pt_1[0]   = Lep1.Pt()
        self.phi_1[0]  = Lep1.Phi()
        self.eta_1[0]  = Lep1.Eta()
        self.pt_2[0]   = Lep2.Pt()
        self.phi_2[0]  = Lep2.Phi()
        self.eta_2[0]  = Lep2.Eta()

	lep_index_1 = lepList[0]
	lep_index_2 = lepList[1]

	if (LepP.Pt() < LepM.Pt()):
	    lep_index_1 = lepList[1]
	    lep_index_2 = lepList[0]
	#relIso 
	if channel_ll == 'ee' : 
      
            self.iso_1[0]  = entry.Electron_pfRelIso03_all[lep_index_1]
            self.iso_2[0]  = entry.Electron_pfRelIso03_all[lep_index_2]
            self.q_1[0]  = entry.Electron_charge[lep_index_1]
            self.q_2[0]  = entry.Electron_charge[lep_index_2]
            self.d0_1[0]   = entry.Electron_dxy[lep_index_1]
            self.dZ_1[0]   = entry.Electron_dz[lep_index_1]
            self.d0_2[0]   = entry.Electron_dxy[lep_index_2]
            self.dZ_2[0]   = entry.Electron_dz[lep_index_2]
            self.Electron_mvaFall17V2noIso_WP90_1[0]  = entry.Electron_mvaFall17V2noIso_WP90[lep_index_1]
            self.Electron_mvaFall17V2noIso_WP90_2[0]  = entry.Electron_mvaFall17V2noIso_WP90[lep_index_2]
	    if SystIndex ==0 and  isMC : 
		self.pt_uncor_1[0] = ePt[lep_index_1]
		self.m_uncor_1[0] = eMass[lep_index_1]
		self.pt_uncor_2[0] = ePt[lep_index_2]
		self.m_uncor_2[0] = eMass[lep_index_2]

            if isMC :
		self.gen_match_1[0] = ord(entry.Electron_genPartFlav[lep_index_1])
		self.gen_match_2[0] = ord(entry.Electron_genPartFlav[lep_index_2])


	if channel_ll == 'mm' : 
            self.iso_1[0]  = entry.Muon_pfRelIso04_all[lep_index_1]
	    self.iso_2[0]  = entry.Muon_pfRelIso04_all[lep_index_2]
	    self.q_1[0]  = entry.Muon_charge[lep_index_1]
	    self.q_2[0]  = entry.Muon_charge[lep_index_2]
	    self.d0_1[0]   = entry.Muon_dxy[lep_index_1]
	    self.dZ_1[0]   = entry.Muon_dz[lep_index_1]
	    self.d0_2[0]   = entry.Muon_dxy[lep_index_2]
	    self.dZ_2[0]   = entry.Muon_dz[lep_index_2]
	    self.looseId_1[0]   = entry.Muon_looseId[lep_index_1] 
	    self.looseId_2[0]   = entry.Muon_looseId[lep_index_2] 
            self.tightId_1[0]      = entry.Muon_tightId[lep_index_1]
            self.tightId_2[0]      = entry.Muon_tightId[lep_index_2]
	    self.mediumId_1[0]   = entry.Muon_mediumId[lep_index_1] 
	    self.mediumId_2[0]   = entry.Muon_mediumId[lep_index_2] 
	    self.mediumPromptId_1[0]   = entry.Muon_mediumPromptId[lep_index_1] 
	    self.mediumPromptId_2[0]   = entry.Muon_mediumPromptId[lep_index_2] 
	    self.isGlobal_1[0]   = entry.Muon_isGlobal[lep_index_1] 
	    self.isGlobal_2[0]   = entry.Muon_isGlobal[lep_index_2] 
	    self.isTracker_1[0]   = entry.Muon_isTracker[lep_index_1] 
	    self.isTracker_2[0]   = entry.Muon_isTracker[lep_index_2] 
	    if SystIndex ==0 and isMC: 
		self.pt_uncor_1[0] = mPt[lep_index_1]
		self.m_uncor_1[0] = mMass[lep_index_1]
		self.pt_uncor_2[0] = mPt[lep_index_2]
		self.m_uncor_2[0] = mMass[lep_index_2]
            if isMC :
		self.gen_match_1[0] = ord(entry.Muon_genPartFlav[lep_index_1])
		self.gen_match_2[0] = ord(entry.Muon_genPartFlav[lep_index_2])

        
        # genMatch the di-lepton variables
	if isMC :
	    idx_Lep1, idx_Lep2 = -1, -1
	    idx_Lep1_tr, idx_Lep2_tr = -1, -1
	    if (Lep1.M() > 0.05 and Lep2.M() > 0.05): # muon mass 
		idx_Lep1 = GF.getLepIdxFrom4Vec(entry, Lep1, 'm')
		idx_Lep2 = GF.getLepIdxFrom4Vec(entry, Lep2, 'm')
		try :
		    idx_Lep1_tr = entry.Muon_genPartIdx[idx_Lep1]
		    idx_Lep2_tr = entry.Muon_genPartIdx[idx_Lep2]
		except IndexError : pass 
		    
	    elif (Lep1.M() < 0.05 and Lep2.M() < 0.05): # electron mass
		idx_Lep1 = GF.getLepIdxFrom4Vec(entry, Lep1, 'e')
		idx_Lep2 = GF.getLepIdxFrom4Vec(entry, Lep2, 'e')
		try :
		    idx_Lep1_tr = entry.Electron_genPartIdx[idx_Lep1]
		    idx_Lep2_tr = entry.Electron_genPartIdx[idx_Lep2]
		except IndexError : pass 
		    
	    if idx_Lep1_tr >= 0 and idx_Lep2_tr >= 0:
		self.m_1_tr[0]  = entry.GenPart_mass[idx_Lep1_tr]
		self.pt_1_tr[0]  = entry.GenPart_pt[idx_Lep1_tr]
		self.m_2_tr[0]  = entry.GenPart_mass[idx_Lep2_tr]
		self.pt_2_tr[0]  = entry.GenPart_pt[idx_Lep2_tr]
		self.eta_1_tr[0] = entry.GenPart_eta[idx_Lep1_tr]
		self.eta_2_tr[0] = entry.GenPart_eta[idx_Lep2_tr]
		self.phi_1_tr[0] = entry.GenPart_phi[idx_Lep1_tr]
		self.phi_2_tr[0] = entry.GenPart_phi[idx_Lep2_tr]
		self.GenPart_statusFlags_1[0]    = entry.GenPart_statusFlags[idx_Lep1_tr]
		self.GenPart_status_1[0]    = entry.GenPart_status[idx_Lep1_tr]
		self.GenPart_statusFlags_2[0]    = entry.GenPart_statusFlags[idx_Lep2_tr]
		self.GenPart_status_2[0]    = entry.GenPart_status[idx_Lep2_tr]
        
        
        #self.btagWeightDeepCSVB[0]  = entry.btagWeight_DeepCSVB
        #print 'inside after filling----------------------->', entry.MET_pt,  self.met[0], met_pt
        #self.puppimet[0]    = entry.PuppiMET_pt
        #self.puppimetphi[0] = entry.PuppiMET_phi

        
	if isMC :
	    self.HTXS_Higgs_cat[0]         = entry.HTXS_stage1_1_cat_pTjet30GeV
	    self.HTXS_Higgs_pt[0]         = entry.HTXS_Higgs_pt
        
        
        # MET variables  at this point this is the TauES corrected MET

        #print 'let see', channel, self.pt_uncor_3[0], self.pt_3[0], self.pt_uncor_4[0], self.pt_4[0], entry.event, SystIndex

	if str(era) != '2017' : 
	    self.metNoCor[0]= entry.MET_pt
	    self.metphiNoCor[0]= entry.MET_phi
	if str(era) == '2017' :
            if proc=='EOY':
		try : 
		    self.metNoCor[0]= entry.METFixEE2017_pt
		    self.metphiNoCor[0]= entry.METFixEE2017_phi
		except AttributeError:
		    self.metNoCor[0]= entry.MET_pt
		    self.metphiNoCor[0]= entry.MET_phi

            if proc=='UL':
		try : 
		    self.metNoCor[0]= entry.MET_pt
		    self.metphiNoCor[0]= entry.MET_phi
		except AttributeError:
		    self.metNoCor[0]= -1
		    self.metphiNoCor[0]= -1

        #print 'inside', met_pt, entry.MET_pt, entry.MET_T1_pt, entry.event, entry.luminosityBlock, entry.run

        if met_pt != -99 : 
	    self.met[0]         = met_pt 
	    self.metphi[0]      = met_phi
            #if SystIndex==0 : print 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAa------------------------------>met corrected TES inside', met_pt, 'noTES corr', entry.MET_T1_pt, entry.event, entry.luminosityBlock, entry.run
          
        else : 
	    if not doUncertainties : 
		if str(era) != '2017' : 
		    self.met[0]= entry.MET_pt
		    self.metphi[0]= entry.MET_phi
		if str(era) == '2017' : 
                    if proc=='EOY':
			try :
			    self.met[0]= entry.METFixEE2017_pt
			    self.metphi[0]= entry.METFixEE2017_phi
			except AttributeError:
			    self.met[0]= entry.MET_pt
			    self.metphi[0]= entry.MET_phi
                    if proc=='UL':
			try :
			    self.met[0]= entry.MET_pt
			    self.metphi[0]= entry.MET_phi
			except AttributeError:
			    self.met[0]= -1
			    self.metphi[0]= -1

	    if  doUncertainties : 

		if str(era) != '2017' : 
                    try : 
			self.met[0]= entry.MET_T1_pt
			self.metphi[0]= entry.MET_T1_phi
                    except AttributeError : 
			self.met[0]= entry.MET_pt
			self.metphi[0]= entry.MET_phi

		if str(era) == '2017' : 
                    if proc=='EOY':
			try : 
			    self.met[0]= entry.METFixEE2017_T1_pt
			    self.metphi[0]= entry.METFixEE2017_T1_phi
			except AttributeError : 
			    self.met[0]= entry.METFixEE2017_pt
			    self.metphi[0]= entry.METFixEE2017_phi
                    if proc=='UL':
			try : 
			    self.met[0]= entry.MET_T1_pt
			    self.metphi[0]= entry.MET_T1_phi
			except AttributeError : 
			    self.met[0]= entry.MET_pt
			    self.metphi[0]= entry.MET_phi

        #metNoTauES holds the uncorrected TauES MET - if not doUncerta -> holds the default ucorrected MET, if doUncert the T1_corrected

        if str(era) != '2017' : 
	    self.metNoTauES[0]         = entry.MET_pt
	    self.metphiNoTauES[0]         = entry.MET_phi

	    if doUncertainties : 
                try : 
		    self.metNoTauES[0]         = entry.MET_T1_pt
		    self.metphiNoTauES[0]         = entry.MET_T1_phi
                except AttributeError : 
		    self.metNoTauES[0]         = entry.MET_pt
		    self.metphiNoTauES[0]         = entry.MET_phi

                if isMC :
		    try : 
			self.MET_T1Smear_pt[0]         = entry.MET_T1Smear_pt
			self.MET_T1Smear_phi[0]         = entry.MET_T1Smear_phi
		    except AttributeError : 
			self.MET_T1Smear_pt[0]         = -99
			self.MET_T1Smear_phi[0]         = -99

        if str(era) == '2017' : 
	    if proc=="EOY" :
		self.metNoTauES[0]         = entry.METFixEE2017_pt
		self.metphiNoTauES[0]         = entry.METFixEE2017_phi

		if doUncertainties : 
		    try :
			self.metNoTauES[0]         = entry.METFixEE2017_T1_pt
			self.metphiNoTauES[0]         = entry.METFixEE2017_T1_phi
		    except AttributeError : 
			self.metNoTauES[0]         = entry.METFixEE2017_pt_nom
			self.metphiNoTauES[0]         = entry.METFixEE2017_phi_nom
		    if isMC :
			try :
			    self.MET_T1Smear_pt[0]         = entry.METFixEE2017_T1Smear_pt
			    self.MET_T1Smear_phi[0]         = entry.METFixEE2017_T1Smear_phi
			except AttributeError : 
			    self.MET_T1Smear_pt[0]         = -1
			    self.MET_T1Smear_phi[0]         = -1

	    if proc=="UL" :
		self.metNoTauES[0]         = entry.MET_pt
		self.metphiNoTauES[0]         = entry.MET_phi

		if doUncertainties : 
		    try :
			self.metNoTauES[0]         = entry.MET_T1_pt
			self.metphiNoTauES[0]         = entry.MET_T1_phi
		    except AttributeError : 
			self.metNoTauES[0]         = entry.MET_pt_nom
			self.metphiNoTauES[0]         = entry.MET_phi_nom
		    if isMC :
			try :
			    self.MET_T1Smear_pt[0]         = entry.MET_T1Smear_pt
			    self.MET_T1Smear_phi[0]         = entry.MET_T1Smear_phi
			except AttributeError : 
			    self.MET_T1Smear_pt[0]         = -1
			    self.MET_T1Smear_phi[0]         = -1

        #print 'in NTUPLE ============================== met_pt', met_pt, 'met', self.met[0], 'metnoTauES', self.metNoTauES[0], 'met_T1', entry.MET_T1_pt, 'met_T1Smear', entry.MET_T1Smear_pt, 'doUncert ?', doUncertainties

        if str(era) != '2017' : 

	    self.metcov00[0] = entry.MET_covXX
	    self.metcov01[0] = entry.MET_covXY
	    self.metcov10[0] = entry.MET_covXY
	    self.metcov11[0] = entry.MET_covYY
	    self.met_UnclX = entry.MET_MetUnclustEnUpDeltaX
	    self.met_UnclY = entry.MET_MetUnclustEnUpDeltaY

	    if doUncertainties : 
		if isMC : 
		    self.MET_pt_UnclUp[0] = entry.MET_pt_unclustEnUp
		    self.MET_phi_UnclUp[0] = entry.MET_phi_unclustEnUp
		    self.MET_pt_UnclDown[0] = entry.MET_pt_unclustEnDown
		    self.MET_phi_UnclDown[0] = entry.MET_phi_unclustEnDown



        else :
            if proc=='EOY' :  
		self.metcov00[0] = entry.METFixEE2017_covXX
		self.metcov01[0] = entry.METFixEE2017_covXY
		self.metcov10[0] = entry.METFixEE2017_covXY
		self.metcov11[0] = entry.METFixEE2017_covYY
		self.met_UnclX = entry.METFixEE2017_MetUnclustEnUpDeltaX
		self.met_UnclY = entry.METFixEE2017_MetUnclustEnUpDeltaY

		if doUncertainties : 
		    if isMC : 
			self.MET_pt_UnclUp[0] = entry.METFixEE2017_pt_unclustEnUp
			self.MET_phi_UnclUp[0] = entry.METFixEE2017_phi_unclustEnUp
			self.MET_pt_UnclDown[0] = entry.METFixEE2017_pt_unclustEnDown
			self.MET_phi_UnclDown[0] = entry.METFixEE2017_phi_unclustEnDown
            if proc=='UL' :  
		self.metcov00[0] = entry.MET_covXX
		self.metcov01[0] = entry.MET_covXY
		self.metcov10[0] = entry.MET_covXY
		self.metcov11[0] = entry.MET_covYY
		self.met_UnclX = entry.MET_MetUnclustEnUpDeltaX
		self.met_UnclY = entry.MET_MetUnclustEnUpDeltaY

		if doUncertainties : 
		    if isMC : 
			self.MET_pt_UnclUp[0] = entry.MET_pt_unclustEnUp
			self.MET_phi_UnclUp[0] = entry.MET_phi_unclustEnUp
			self.MET_pt_UnclDown[0] = entry.MET_pt_unclustEnDown
			self.MET_phi_UnclDown[0] = entry.MET_phi_unclustEnDown

        # trig
        if SystIndex ==0 : 
	    self.isTrig_1[0]   = is_trig_1
	    self.isTrig_2[0]   = is_trig_2
	    self.isDoubleTrig[0]   = is_Dtrig_1

        leplist=[]
        leplist.append(LepP)
        leplist.append(LepM)
	if jt1>-1 and jt2>-1 :  
	    leplist.append(tau1)
	    leplist.append(tau2)

        if doUncertainties: 
                ## this is not done from within ZH and the correctallMET function
                for i, v in enumerate(self.allsystMET) : 

                    if str(era)=='2017' :
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
	    if len(jetList) > 0 :
		jpt1 = getattr(entry, "Jet_pt{0:s}".format(str(isys)), None)
		jj1 = jetList[0]
		self.jpt_1[ic]  = jpt1[jj1]
		self.jeta_1[ic] = entry.Jet_eta[jj1]
		self.jphi_1[ic] = entry.Jet_phi[jj1]
		self.jcsv_1[ic] = entry.Jet_btagDeepB[jj1]
		self.jcsvfv_1[ic] = entry.Jet_btagDeepFlavB[jj1]
                #print 'will use', ic, len(jetList), jetList, self.jpt_1[ic], self.njets[ic]
		
		# genMatch jet1
		if isMC:
		    idx_genJet = entry.Jet_genJetIdx[jj1]
		    if idx_genJet >= 0:
			try :
			    self.jpt_1_tr[ic]  = entry.GenJet_pt[idx_genJet]
			    self.jeta_1_tr[ic] = entry.GenJet_eta[idx_genJet]
			    self.jphi_1_tr[ic] = entry.GenJet_phi[idx_genJet]
			except IndexError : pass

	    self.jpt_2[ic], self.jeta_2[ic], self.jphi_2[ic], self.jcsv_2[ic],self.jcsvfv_2[ic] = -9.99, -9.99, -9.99, -9.99, -9.99
	    if len(jetList) > 1 :
		jpt2 = getattr(entry, "Jet_pt{0:s}".format(str(isys)), None)
		jj2 = jetList[1] 
		self.jpt_2[ic]  = jpt2[jj2]
		self.jeta_2[ic] = entry.Jet_eta[jj2]
		self.jphi_2[ic] = entry.Jet_phi[jj2]
		self.jcsv_2[ic] = entry.Jet_btagDeepB[jj2]
		self.jcsvfv_2[ic] = entry.Jet_btagDeepFlavB[jj2]
		
		# genMatch jet2
		if isMC:
		    idx_genJet = entry.Jet_genJetIdx[jj2]
		    if idx_genJet >= 0:
			try: 
			   self.jpt_2_tr[ic]  = entry.GenJet_pt[idx_genJet]
			   self.jeta_2_tr[ic] = entry.GenJet_eta[idx_genJet]
			   self.jphi_2_tr[ic] = entry.GenJet_phi[idx_genJet]
			except IndexError : pass 

	    self.bpt_1[ic], self.beta_1[ic], self.bphi_1[ic], self.bcsv_1[ic], self.bcsvfv_1[ic] = -9.99, -9.99, -9.99, -9.99, -9.99
	    if len(bJetList) > 0 :
		jpt1 = getattr(entry, "Jet_pt{0:s}".format(str(isys)), None)
		jbj1 = bJetList[0]
		self.bpt_1[ic] = jpt1[jbj1]
		self.beta_1[ic] = entry.Jet_eta[jbj1]
		self.bphi_1[ic] = entry.Jet_phi[jbj1]
		self.bcsv_1[ic] = entry.Jet_btagDeepB[jbj1] 
		self.bcsvfv_1[ic] = entry.Jet_btagDeepFlavB[jbj1]
		
	    self.bpt_2[ic], self.beta_2[ic], self.bphi_2[ic], self.bcsv_2[ic], self.bcsvfv_2[ic] = -9.99, -9.99, -9.99, -9.99, -9.99
	    if len(bJetList) > 1 :
		jpt2 = getattr(entry, "Jet_pt{0:s}".format(str(isys)), None)
		jbj2 = bJetList[1] 
		self.bpt_2[ic] = jpt2[jbj2]
		self.beta_2[ic] = entry.Jet_eta[jbj2]
		self.bphi_2[ic] = entry.Jet_phi[jbj2]
		self.bcsv_2[ic] = entry.Jet_btagDeepB[jbj2]
		self.bcsvfv_2[ic] = entry.Jet_btagDeepFlavB[jbj2]

		# genMatch bjet1
		if isMC:
		    idx_genJet = entry.Jet_genJetIdx[jbj2]
		    if idx_genJet >= 0:
			try :
			    self.bpt_2_tr[ic]  = entry.GenJet_pt[idx_genJet]
			    self.beta_2_tr[ic] = entry.GenJet_eta[idx_genJet]
			    self.bphi_2_tr[ic] = entry.GenJet_phi[idx_genJet]
			except IndexError : pass

        '''

        #if  self.nbtag[0] == 0 : 
	if SystIndex == 0 : 
            self.t.Fill()
	else : 
            self.tN[SystIndex-1].Fill()

	return

    def Fill3L(self, entry, cat, LepP, LepM, lepList, lepList_2, ElList, MuList, TauList, isMC, era, doUncertainties=False,met_pt=-99, met_phi=-99,systIndex=0) :
    #def Fill(self, entry, SVFit, cat, jt1, jt2, LepP, LepM, lepList, isMC, era, doUncertainties=False ,  met_pt=-99, met_phi=-99, systIndex=0, tMass=[], tPt=[], eMass=[], ePt=[], mMass=[], mPt=[]) : 
        ''' - jt1 point to the selected tau candidates according to the table below.
            - LepP and LepM are TLorentz vectors for the positive and negative members of the dilepton pair
        ''' 
        nelectrons = len(ElList)
        nmuons = len(MuList)
        ntaus = len(TauList)
        channel_ll = cat[:2]
        channel = cat[-2:]

        SystIndex = int(systIndex)

        if SystIndex ==0 : 
	    is_trig_1, is_trig_2, is_Dtrig_1 = 0,0,0
	    TrigListLep = []
	    TrigListTau = []
	    hltListLep  = []

	    TrigListLep, hltListLep, hltListLepSubL  = GF.findSingleLeptTrigger(lepList, entry, channel_ll, era, True)

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

        
	    self.entries += 1

	    self.run[0]  = entry.run
	    self.nElectron[0]  = nelectrons
	    self.nMuon[0]  = nmuons
	    self.nTau[0]  = ntaus
	    self.lumi[0] = entry.luminosityBlock 
	    self.evt[0]  = entry.event
	    self.iso_1[0]  = -99
	    self.iso_2[0]  = -99
	    self.q_1[0]  = -99
	    self.q_2[0]  = -99
	    self.isGlobal_1[0]  = -99
	    self.isGlobal_2[0]  = -99
	    try:
		self.L1PreFiringWeight_Nom[0] = entry.L1PreFiringWeight_Nom
		self.L1PreFiringWeight_Up[0] = entry.L1PreFiringWeight_Up
		self.L1PreFiringWeight_Down[0] = entry.L1PreFiringWeight_Down
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

	    self.tightId_2[0]       = -1 
	    self.mediumId_2[0]       = -1 
	    self.mediumPromptId_2[0]   = -1
	    self.looseId_2[0]       = -1
	    self.isGlobal_2[0]      = -1
	    self.isTracker_2[0]     = -1

	   
	    self.decayMode_3[0]        = -1
	    self.idDecayModeNewDMs_3[0]= -1
	    self.idDeepTau2017v2p1VSe_3[0] = -1
	    self.idDeepTau2017v2p1VSjet_3[0] = -1
	    self.idDeepTau2017v2p1VSmu_3[0] = -1
	    self.idMVAnewDM2017v2_3[0] = -1
	    self.rawMVAnewDM2017v2_3[0] = -1
	    self.mediumId_3[0]       = -1 
	    self.mediumPromptId_3[0]   = -1
	    self.looseId_3[0]       = -1
	    self.isGlobal_3[0]      = -1
	    self.isTracker_3[0]     = -1
	    self.ip3d_3[0]          = -1

	    self.decayMode_4[0]      = -1
	    self.idDecayModeNewDMs_4[0] = -1
	    self.idDeepTau2017v2p1VSe_4[0] = -1
	    self.idDeepTau2017v2p1VSjet_4[0] = -1
	    self.idDeepTau2017v2p1VSmu_4[0] = -1
	    self.idMVAnewDM2017v2_4[0] = -1
	    self.rawMVAnewDM2017v2_4[0] = -1
	    self.mediumId_4[0]       = -1 
	    self.mediumPromptId_4[0]   = -1
	    self.looseId_4[0]       = -1
	    self.isGlobal_4[0]      = -1
	    self.isTracker_4[0]     = -1
	    self.ip3d_4[0]          = -1
	    self.gen_match_1[0] = -1
	    self.gen_match_2[0] = -1
	    self.gen_match_3[0] = -1
	    self.gen_match_4[0] = -1

	    goodElectronList = tauFunDCH.makeGoodElectronList(entry)
	    goodMuonList = tauFunDCH.makeGoodMuonList(entry)
	    
	    self.nGoodElectron[0] =  nelectrons
	    self.nGoodMuon[0]     = nmuons
	    self.cat[0]  = tauFunDCH.catToNumber3L(cat)
	    

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

	    # pack trigger bits into integer word
	    year = int(era)
	    e = entry

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
	    for i in range(4) :
		bits.append(False)
		
	    try : bits.append(e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ)
	    except AttributeError : bits.append(False) 

	    self.electronTriggerWord[0] = 0
	    for i, bit in enumerate(bits) :
		if bit : self.electronTriggerWord[0] += 2**i


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

	    bits.append(False)                             # pad remaining bit in this bit 
	   
	    try : bits.append(e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ)
	    except AttributeError : bits.append(False) 
	    try : bits.append(e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)
	    except AttributeError : bits.append(False)
	    try : bits.append(e.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ)
	    except AttributeError : bits.append(False) 
	    try : bits.append(e.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_Mass8)
	    except AttributeError : bits.append(False) 


	    self.muonTriggerWord[0] = 0
	    for i, bit in enumerate(bits) :
		if bit : self.muonTriggerWord[0] += 2**i


        # Sort the di-lepton system by Pt
        Lep1, Lep2 = TLorentzVector(), TLorentzVector()
	lep_index_1 = lepList[0]
	lep_index_2 = lepList[1]

        if (LepP.Pt() > LepM.Pt()): 
            Lep1 = LepP
            Lep2 = LepM
        else:
            Lep1 = LepM
            Lep2 = LepP
	    lep_index_1 = lepList[1]
	    lep_index_2 = lepList[0]

        # di-lepton variables.   _p and _m refer to plus and minus charge
        self.mll[0]       = (Lep1 + Lep2).M()
        self.Z_DR[0]       = self.getDR(entry,Lep1,Lep2)

        self.H_LT[0]       = Lep1.Pt() + Lep2.Pt()
           
        self.pt_1[0]   = Lep1.Pt()
        self.phi_1[0]  = Lep1.Phi()
        self.eta_1[0]  = Lep1.Eta()
        self.pt_2[0]   = Lep2.Pt()
        self.phi_2[0]  = Lep2.Phi()
        self.eta_2[0]  = Lep2.Eta()

	if channel_ll == 'ee' : 
      
            self.iso_1[0]  = entry.Electron_pfRelIso03_all[lep_index_1]
            self.iso_2[0]  = entry.Electron_pfRelIso03_all[lep_index_2]
            self.q_1[0]  = entry.Electron_charge[lep_index_1]
            self.q_2[0]  = entry.Electron_charge[lep_index_2]
            self.d0_1[0]   = entry.Electron_dxy[lep_index_1]
            self.dZ_1[0]   = entry.Electron_dz[lep_index_1]
            self.d0_2[0]   = entry.Electron_dxy[lep_index_2]
            self.dZ_2[0]   = entry.Electron_dz[lep_index_2]
            self.Electron_mvaFall17V2noIso_WP90_1[0]  = entry.Electron_mvaFall17V2noIso_WP90[lep_index_1]
            self.Electron_mvaFall17V2noIso_WP90_2[0]  = entry.Electron_mvaFall17V2noIso_WP90[lep_index_2]
            if isMC :
		self.gen_match_1[0] = ord(entry.Electron_genPartFlav[lep_index_1])
		self.gen_match_2[0] = ord(entry.Electron_genPartFlav[lep_index_2])

	if channel_ll == 'mm' : 

            self.iso_1[0]  = entry.Muon_pfRelIso04_all[lep_index_1]
	    self.iso_2[0]  = entry.Muon_pfRelIso04_all[lep_index_2]
	    self.q_1[0]  = entry.Muon_charge[lep_index_1]
	    self.q_2[0]  = entry.Muon_charge[lep_index_2]
	    self.d0_1[0]   = entry.Muon_dxy[lep_index_1]
	    self.dZ_1[0]   = entry.Muon_dz[lep_index_1]
	    self.d0_2[0]   = entry.Muon_dxy[lep_index_2]
	    self.dZ_2[0]   = entry.Muon_dz[lep_index_2]
	    self.looseId_1[0]   = entry.Muon_looseId[lep_index_1] 
	    self.looseId_2[0]   = entry.Muon_looseId[lep_index_2] 
            self.tightId_1[0]      = entry.Muon_tightId[lep_index_1]
            self.tightId_2[0]      = entry.Muon_tightId[lep_index_2]
	    self.mediumId_1[0]   = entry.Muon_mediumId[lep_index_1] 
	    self.mediumId_2[0]   = entry.Muon_mediumId[lep_index_2] 
	    self.mediumPromptId_1[0]   = entry.Muon_mediumPromptId[lep_index_1] 
	    self.mediumPromptId_2[0]   = entry.Muon_mediumPromptId[lep_index_2] 
	    self.isGlobal_1[0]   = entry.Muon_isGlobal[lep_index_1] 
	    self.isGlobal_2[0]   = entry.Muon_isGlobal[lep_index_2] 
	    self.isTracker_1[0]   = entry.Muon_isTracker[lep_index_1] 
	    self.isTracker_2[0]   = entry.Muon_isTracker[lep_index_2] 
            if isMC :
		self.gen_match_1[0] = ord(entry.Muon_genPartFlav[lep_index_1])
		self.gen_match_2[0] = ord(entry.Muon_genPartFlav[lep_index_2])


        #print ElList, MuList, TauList
        eL1, eL2= TLorentzVector(), TLorentzVector()  
        if nelectrons > 0 :
	    ie = ElList[0]

            if len(ElList)>1 :

                if entry.Electron_pt[ElList[0]] > entry.Electron_pt[ElList[1]] : 
		    ie = ElList[0] 
		    iee = ElList[1] 
		else : 
		    ie = ElList[1] 
		    iee = ElList[0] 
	    self.pt_3[0]   = entry.Electron_pt[ie]
	    self.phi_3[0]  = entry.Electron_phi[ie]
	    self.eta_3[0]  = entry.Electron_eta[ie]
	    self.m_3[0]    = entry.Electron_mass[ie]
	    self.q_3[0]    = entry.Electron_charge[ie]
	    self.d0_3[0]   = entry.Electron_dxy[ie]
	    self.dZ_3[0]   = entry.Electron_dz[ie]
	    self.iso_3[0]  = entry.Electron_pfRelIso03_all[ie]
	    self.Electron_mvaFall17V2noIso_WP90_3[0]  = entry.Electron_mvaFall17V2noIso_WP90[ie]
	    try : self.gen_match_3[0] = ord(entry.Electron_genPartFlav[ie])
	    except AttributeError : self.gen_match_3[0] = -1


            if len(ElList)>1 :
		self.pt_4[0]   = entry.Electron_pt[iee]
		self.phi_4[0]  = entry.Electron_phi[iee]
		self.eta_4[0]  = entry.Electron_eta[iee]
		self.m_4[0]    = entry.Electron_mass[iee]
		self.q_4[0]    = entry.Electron_charge[iee]
		self.d0_4[0]   = entry.Electron_dxy[iee]
		self.dZ_4[0]   = entry.Electron_dz[iee]
		self.iso_4[0]  = entry.Electron_pfRelIso03_all[iee]
		self.Electron_mvaFall17V2noIso_WP90_4[0]  = entry.Electron_mvaFall17V2noIso_WP90[iee]
		try : self.gen_match_4[0] = ord(entry.Electron_genPartFlav[iee])
		except AttributeError : self.gen_match_4[0] = -1

		eL1.SetPtEtaPhiM(entry.Electron_pt[ElList[0]],entry.Electron_eta[ElList[0]], entry.Electron_phi[ElList[0]], electronMass)
		eL2.SetPtEtaPhiM(entry.Electron_pt[ElList[1]],entry.Electron_eta[ElList[1]], entry.Electron_phi[ElList[1]], electronMass)
		self.mll2[0] = (eL1 + eL2).M()


        if nmuons > 0:
            
	    im = MuList[0]

            if len(MuList)>1 :
                if  entry.Muon_pt[MuList[0]] > entry.Muon_pt[MuList[1]] : 
		    im = MuList[0] 
		    imm = MuList[1] 
		else : 
		    im = MuList[1] 
		    imm = MuList[0] 

	    self.pt_3[0]     = entry.Muon_pt[im]
	    self.phi_3[0]    = entry.Muon_phi[im]
	    self.eta_3[0]    = entry.Muon_eta[im]
	    self.m_3[0]      = entry.Muon_mass[im]
	    self.q_3[0]      = entry.Muon_charge[im]
	    self.d0_3[0]     = entry.Muon_dxy[im]
	    self.dZ_3[0]     = entry.Muon_dz[im]
	    self.iso_3[0]    = entry.Muon_pfRelIso04_all[im]
	    self.tightId_3[0]   = entry.Muon_tightId[im] 
	    self.mediumId_3[0]      = entry.Muon_mediumId[im]
	    self.tightId_3[0]      = entry.Muon_tightId[im]
	    self.mediumPromptId_3[0]   = entry.Muon_mediumPromptId[im]
	    self.looseId_3[0]       = entry.Muon_looseId[im]
	    self.isGlobal_3[0]      = entry.Muon_isGlobal[im]
	    self.isTracker_3[0]     = entry.Muon_isTracker[im]
	    self.ip3d_3[0]       = entry.Muon_ip3d[im]
	    try : self.gen_match_3[0] = ord(entry.Muon_genPartFlav[im])
	    except AttributeError : self.gen_match_3[0] = -1


            if len(MuList)>1 :

		self.pt_4[0]     = entry.Muon_pt[imm]
		self.phi_4[0]    = entry.Muon_phi[imm]
		self.eta_4[0]    = entry.Muon_eta[imm]
		self.m_4[0]      = entry.Muon_mass[imm]
		self.q_4[0]      = entry.Muon_charge[imm]
		self.d0_4[0]     = entry.Muon_dxy[imm]
		self.dZ_4[0]     = entry.Muon_dz[imm]
		self.iso_4[0]    = entry.Muon_pfRelIso04_all[imm]
		self.tightId_4[0]   = entry.Muon_tightId[imm] 
		self.mediumId_4[0]      = entry.Muon_mediumId[imm]
		self.mediumPromptId_4[0]   = entry.Muon_mediumPromptId[imm]
		self.looseId_4[0]       = entry.Muon_looseId[imm]
		self.isGlobal_4[0]      = entry.Muon_isGlobal[imm]
		self.isTracker_4[0]     = entry.Muon_isTracker[imm]
		self.ip3d_4[0]       = entry.Muon_ip3d[imm]
		try : self.gen_match_4[0] = ord(entry.Muon_genPartFlav[imm])
		except AttributeError : self.gen_match_4[0] = -1
		eL1.SetPtEtaPhiM(entry.Muon_pt[MuList[0]],entry.Muon_eta[MuList[0]], entry.Muon_phi[MuList[0]], muonMass)
		eL2.SetPtEtaPhiM(entry.Muon_pt[MuList[1]],entry.Muon_eta[MuList[1]], entry.Muon_phi[MuList[1]], muonMass)
		self.mll2[0] = (eL1 + eL2).M()

        # genMatch the di-lepton variables
	if isMC :
	    idx_Lep1, idx_Lep2 = -1, -1
	    idx_Lep1_tr, idx_Lep2_tr = -1, -1
	    if (Lep1.M() > 0.05 and Lep2.M() > 0.05): # muon mass 
		idx_Lep1 = GF.getLepIdxFrom4Vec(entry, Lep1, 'm')
		idx_Lep2 = GF.getLepIdxFrom4Vec(entry, Lep2, 'm')
		try :
		    idx_Lep1_tr = entry.Muon_genPartIdx[idx_Lep1]
		    idx_Lep2_tr = entry.Muon_genPartIdx[idx_Lep2]
		except IndexError : pass 
		    
	    elif (Lep1.M() < 0.05 and Lep2.M() < 0.05): # electron mass
		idx_Lep1 = GF.getLepIdxFrom4Vec(entry, Lep1, 'e')
		idx_Lep2 = GF.getLepIdxFrom4Vec(entry, Lep2, 'e')
		try :
		    idx_Lep1_tr = entry.Electron_genPartIdx[idx_Lep1]
		    idx_Lep2_tr = entry.Electron_genPartIdx[idx_Lep2]
		except IndexError : pass 
		    
	    if idx_Lep1_tr >= 0 and idx_Lep2_tr >= 0:
		self.m_1_tr[0]  = entry.GenPart_mass[idx_Lep1_tr]
		self.pt_1_tr[0]  = entry.GenPart_pt[idx_Lep1_tr]
		self.m_2_tr[0]  = entry.GenPart_mass[idx_Lep2_tr]
		self.pt_2_tr[0]  = entry.GenPart_pt[idx_Lep2_tr]
		self.eta_1_tr[0] = entry.GenPart_eta[idx_Lep1_tr]
		self.eta_2_tr[0] = entry.GenPart_eta[idx_Lep2_tr]
		self.phi_1_tr[0] = entry.GenPart_phi[idx_Lep1_tr]
		self.phi_2_tr[0] = entry.GenPart_phi[idx_Lep2_tr]
		self.GenPart_statusFlags_1[0]    = entry.GenPart_statusFlags[idx_Lep1_tr]
		self.GenPart_status_1[0]    = entry.GenPart_status[idx_Lep1_tr]
		self.GenPart_statusFlags_2[0]    = entry.GenPart_statusFlags[idx_Lep2_tr]
		self.GenPart_status_2[0]    = entry.GenPart_status[idx_Lep2_tr]
        
	if str(era) != '2017' : 
	    self.metNoCor[0]= entry.MET_pt
	    self.metphiNoCor[0]= entry.MET_phi
	if str(era) == '2017' : 
	    self.metNoCor[0]= entry.METFixEE2017_pt
	    self.metphiNoCor[0]= entry.METFixEE2017_phi

        #print 'inside', met_pt, entry.MET_pt, entry.MET_T1_pt, entry.event, entry.luminosityBlock, entry.run

        if met_pt != -99 : 
	    self.met[0]         = met_pt 
	    self.metphi[0]      = met_phi
          
        else : 
	    if not doUncertainties : 
		if str(era) != '2017' : 
		    self.met[0]= entry.MET_pt
		    self.metphi[0]= entry.MET_phi
		if str(era) == '2017' : 
		    self.met[0]= entry.METFixEE2017_pt
		    self.metphi[0]= entry.METFixEE2017_phi
	    if  doUncertainties : 

		if str(era) != '2017' : 
                    try : 
			self.met[0]= entry.MET_T1_pt
			self.metphi[0]= entry.MET_T1_phi
                    except AttributeError : 
			self.met[0]= entry.MET_pt
			self.metphi[0]= entry.MET_phi

		if str(era) == '2017' : 
                    try : 
			self.met[0]= entry.METFixEE2017_T1_pt
			self.metphi[0]= entry.METFixEE2017_T1_phi
                    except AttributeError : 
			self.met[0]= entry.METFixEE2017_pt
			self.metphi[0]= entry.METFixEE2017_phi

        #metNoTauES holds the uncorrected TauES MET - if not doUncerta -> holds the default ucorrected MET, if doUncert the T1_corrected

        if str(era) != '2017' : 
	    self.metNoTauES[0]         = entry.MET_pt
	    self.metphiNoTauES[0]         = entry.MET_phi

	    if doUncertainties : 
                try : 
		    self.metNoTauES[0]         = entry.MET_T1_pt
		    self.metphiNoTauES[0]         = entry.MET_T1_phi
                except AttributeError : 
		    self.metNoTauES[0]         = entry.MET_pt
		    self.metphiNoTauES[0]         = entry.MET_phi

                if isMC :
		    try : 
			self.MET_T1Smear_pt[0]         = entry.MET_T1Smear_pt
			self.MET_T1Smear_phi[0]         = entry.MET_T1Smear_phi
		    except AttributeError : 
			self.MET_T1Smear_pt[0]         = -99
			self.MET_T1Smear_phi[0]         = -99

        if str(era) == '2017' : 
	    self.metNoTauES[0]         = entry.METFixEE2017_pt
	    self.metphiNoTauES[0]         = entry.METFixEE2017_phi

	    if doUncertainties : 
                try :
		    self.metNoTauES[0]         = entry.METFixEE2017_T1_pt
		    self.metphiNoTauES[0]         = entry.METFixEE2017_T1_phi
                except AttributeError : 
		    self.metNoTauES[0]         = entry.METFixEE2017_pt_nom
		    self.metphiNoTauES[0]         = entry.METFixEE2017_phi_nom
                if isMC :
		    try :
			self.MET_T1Smear_pt[0]         = entry.METFixEE2017_T1Smear_pt
			self.MET_T1Smear_phi[0]         = entry.METFixEE2017_T1Smear_phi
		    except AttributeError : 
			self.MET_T1Smear_pt[0]         = -1
			self.MET_T1Smear_phi[0]         = -1


        #print 'in NTUPLE ============================== met_pt', met_pt, 'met', self.met[0], 'metnoTauES', self.metNoTauES[0], 'met_T1', entry.MET_T1_pt, 'met_T1Smear', entry.MET_T1Smear_pt, 'doUncert ?', doUncertainties

        if str(era) != '2017' : 

	    self.metcov00[0] = entry.MET_covXX
	    self.metcov01[0] = entry.MET_covXY
	    self.metcov10[0] = entry.MET_covXY
	    self.metcov11[0] = entry.MET_covYY
	    self.met_UnclX = entry.MET_MetUnclustEnUpDeltaX
	    self.met_UnclY = entry.MET_MetUnclustEnUpDeltaY

	    if doUncertainties : 
		if isMC : 
		    self.MET_pt_UnclUp[0] = entry.MET_pt_unclustEnUp
		    self.MET_phi_UnclUp[0] = entry.MET_phi_unclustEnUp
		    self.MET_pt_UnclDown[0] = entry.MET_pt_unclustEnDown
		    self.MET_phi_UnclDown[0] = entry.MET_phi_unclustEnDown



        else : 
	    self.metcov00[0] = entry.METFixEE2017_covXX
	    self.metcov01[0] = entry.METFixEE2017_covXY
	    self.metcov10[0] = entry.METFixEE2017_covXY
	    self.metcov11[0] = entry.METFixEE2017_covYY
	    self.met_UnclX = entry.METFixEE2017_MetUnclustEnUpDeltaX
	    self.met_UnclY = entry.METFixEE2017_MetUnclustEnUpDeltaY

	    if doUncertainties : 
		if isMC : 
		    self.MET_pt_UnclUp[0] = entry.METFixEE2017_pt_unclustEnUp
		    self.MET_phi_UnclUp[0] = entry.METFixEE2017_phi_unclustEnUp
		    self.MET_pt_UnclDown[0] = entry.METFixEE2017_pt_unclustEnDown
		    self.MET_phi_UnclDown[0] = entry.METFixEE2017_phi_unclustEnDown

        # trig
        if SystIndex ==0 : 
	    self.isTrig_1[0]   = is_trig_1
	    self.isTrig_2[0]   = is_trig_2
	    self.isDoubleTrig[0]   = is_Dtrig_1



        # jet variables
        leplist=[]
        leplist.append(LepP)
        leplist.append(LepM)
        if eL1.Pt()> 0 : leplist.append(eL1)
        if eL2.Pt()> 0 : leplist.append(eL2)


        if doUncertainties: 
                ## this is not done from within ZH and the correctallMET function
                for i, v in enumerate(self.allsystMET) : 

                    if str(era)=='2017' :
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
		    jetList, jetListFlav, jetListEta, jetListPt,bJetList, bJetListT, bJetListFlav = self.getJetsJMEMV(entry,leplist,era,v) 
                    #print 'jessyst', systematic, len(jetList), cat

	            self.list_of_arraysJetsNjets[i][0] = len(jetList)
	            self.list_of_arraysJetsNbtag[i][0] = len(bJetList)
		    for ifl in range(len(jetList)) :
			self.list_of_arraysJetsPt[i][ifl] = jetListPt[ifl]
			self.list_of_arraysJetsEta[i][ifl] = jetListEta[ifl]
			self.list_of_arraysJetsFlavour[i][ifl] = jetListFlav[ifl]


        #fill the un-corrected or just in the case you dont care to doUncertainties       
        nom_=''
	jetList, jetListFlav, jetListEta, jetListPt,bJetList, bJetListT, bJetListFlav = self.getJetsJMEMV(entry,leplist,era,'') 
	self.njets[0] = len(jetList)
	self.nbtag[0] = len(bJetList)
	self.nbtagT[0] = len(bJetList)
	for ifl in range(len(jetListPt)) :
	    self.jflavour[ifl]  = jetListFlav[ifl]
	    self.jeta[ifl]  = jetListEta[ifl]
	    self.jpt[ifl]  = jetListPt[ifl]


       
        #if self.isTrig_1 !=0 : 
	if SystIndex == 0 and ( self.isTrig_1 !=0 or  self.isTrig_2 !=0) : 
            self.t.Fill()
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

