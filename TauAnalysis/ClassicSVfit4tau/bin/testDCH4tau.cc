
/**
   \class testClassicSVfit4tau testClassicSVfit4tau.cc "TauAnalysis/ClassicSVfit/bin/testClassicSVfit4tau.cc"
   \brief Basic example of the use of the standalone version of the "classic" SVfit algorithm, customized for the hh->4tau analysis
*/

#include "TauAnalysis/ClassicSVfit4tau/interface/ClassicSVfit4tau.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // LorentzVector
#include "TauAnalysis/ClassicSVfit4tau/interface/svFitHistogramAdapter4tau.h"
#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"
#include "../../../MET_stuff/MyBranc.C"
#include "../../../MET_stuff/cat.h"
//#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBall2.h"

#include "Math/GenVector/Boost.h"

#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include "TH1F.h"
#include "TMath.h"
#include <string>

using namespace classic_svFit;


auto runSVfit4tau(std::string cat, double METx, double METy, TMatrixD covMET, ROOT::Math::PtEtaPhiMVector tau1, ROOT::Math::PtEtaPhiMVector tau2, ROOT::Math::PtEtaPhiMVector tau3, ROOT::Math::PtEtaPhiMVector tau4, int dm1, int dm2, int dm3, int dm4){
	
	float tmass1 = 1.21648, tmass2 = 1.21648, tmass3 = 1.21648, tmass4 = 1.21648;
	if (dm1 == 0) tmass1 = 0.13957;
	if (dm2 == 0) tmass2 = 0.13957;
	if (dm3 == 0) tmass3 = 0.13957;
	if (dm4 == 0) tmass4 = 0.13957;

	int eleType1 = MeasuredTauLepton::kTauToElecDecay, eleType2 = MeasuredTauLepton::kTauToElecDecay,
		muType1 = MeasuredTauLepton::kTauToMuDecay, muType2 = MeasuredTauLepton::kTauToMuDecay;
	
	/*for (int seq = 1; seq < 16; seq++){
		if (seq == 1){
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			eleType1 = MeasuredTauLepton::kPrompt;
			muType1 = MeasuredTauLepton::kPrompt; 
			muType2 = MeasuredTauLepton::kPrompt;
		}
		if (seq == 2){
			eleType1 = MeasuredTauLepton::kPrompt;
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			muType1 = MeasuredTauLepton::kPrompt; 
			muType2 = MeasuredTauLepton::kPrompt;
		}
		if (seq == 3){
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			muType1 = MeasuredTauLepton::kPrompt; 
			muType2 = MeasuredTauLepton::kPrompt;
		}
		if (seq == 4){
			eleType1 = MeasuredTauLepton::kPrompt;
			eleType1 = MeasuredTauLepton::kPrompt;
			muType1 = MeasuredTauLepton::kTauToMuDecay;
			muType2 = MeasuredTauLepton::kPrompt;
		}
		if (seq == 5){
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			eleType1 = MeasuredTauLepton::kPrompt;
			muType1 = MeasuredTauLepton::kTauToMuDecay;
			muType2 = MeasuredTauLepton::kPrompt;
		}
		if (seq == 6){
			eleType1 = MeasuredTauLepton::kPrompt;
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			muType1 = MeasuredTauLepton::kTauToMuDecay;
			muType2 = MeasuredTauLepton::kPrompt;
		}
		if (seq == 7){
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			muType1 = MeasuredTauLepton::kTauToMuDecay;
			muType2 = MeasuredTauLepton::kPrompt;
		}		
		if (seq == 8){
			eleType1 = MeasuredTauLepton::kPrompt;
			eleType1 = MeasuredTauLepton::kPrompt;
			muType1 = MeasuredTauLepton::kPrompt;
			muType2 = MeasuredTauLepton::kTauToMuDecay;
		}
		if (seq == 9){
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			eleType1 = MeasuredTauLepton::kPrompt;
			muType1 = MeasuredTauLepton::kPrompt; 
			muType2 = MeasuredTauLepton::kTauToMuDecay;
		}
		if (seq == 10){
			eleType1 = MeasuredTauLepton::kPrompt;
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			muType1 = MeasuredTauLepton::kPrompt; 
			muType2 = MeasuredTauLepton::kTauToMuDecay;
		}
		if (seq == 11){
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			muType1 = MeasuredTauLepton::kPrompt; 
			muType2 = MeasuredTauLepton::kTauToMuDecay;
		}
		if (seq == 12){
			eleType1 = MeasuredTauLepton::kPrompt;
			eleType1 = MeasuredTauLepton::kPrompt;
			muType1 = MeasuredTauLepton::kTauToMuDecay;
			muType2 = MeasuredTauLepton::kTauToMuDecay;
		}
		if (seq == 13){
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			eleType1 = MeasuredTauLepton::kPrompt;
			muType1 = MeasuredTauLepton::kTauToMuDecay;
			muType2 = MeasuredTauLepton::kTauToMuDecay;
		}
		if (seq == 14){
			eleType1 = MeasuredTauLepton::kPrompt;
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			muType1 = MeasuredTauLepton::kTauToMuDecay;
			muType2 = MeasuredTauLepton::kTauToMuDecay;
		}
		if (seq == 15){
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			eleType1 = MeasuredTauLepton::kTauToElecDecay;
			muType1 = MeasuredTauLepton::kTauToMuDecay;
			muType2 = MeasuredTauLepton::kTauToMuDecay;
		}*/

		std::vector<MeasuredTauLepton> measuredTauLeptons;
		if (cat.substr(0,2) == "ee"){
			measuredTauLeptons.push_back(MeasuredTauLepton(eleType1, tau1.Pt(), tau1.Eta(), tau1.Phi(), 0.000511));
			measuredTauLeptons.push_back(MeasuredTauLepton(eleType2, tau2.Pt(), tau2.Eta(), tau2.Phi(), 0.000511));
		}
		else if (cat.substr(0,2) == "em"){
			measuredTauLeptons.push_back(MeasuredTauLepton(eleType1, tau1.Pt(), tau1.Eta(), tau1.Phi(), 0.000511));
			measuredTauLeptons.push_back(MeasuredTauLepton(muType2, tau2.Pt(), tau2.Eta(), tau2.Phi(), 0.106));
		}
		else if (cat.substr(0,2) == "et"){
			measuredTauLeptons.push_back(MeasuredTauLepton(eleType1, tau1.Pt(), tau1.Eta(), tau1.Phi(), 0.000511));
			measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, tau2.Pt(), tau2.Eta(), tau2.Phi(), tmass2, dm2));
		}
		else if (cat.substr(0,2) == "mm"){
			measuredTauLeptons.push_back(MeasuredTauLepton(muType1, tau1.Pt(), tau1.Eta(), tau1.Phi(), 0.106));
			measuredTauLeptons.push_back(MeasuredTauLepton(muType2, tau2.Pt(), tau2.Eta(), tau2.Phi(), 0.106));
		}
		else if (cat.substr(0,2) == "mt"){
			measuredTauLeptons.push_back(MeasuredTauLepton(muType1, tau1.Pt(), tau1.Eta(), tau1.Phi(), 0.106));
			measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, tau2.Pt(), tau2.Eta(), tau2.Phi(), tmass2, dm2));
		}
		else if (cat.substr(0,2) == "tt"){
			measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, tau1.Pt(), tau1.Eta(), tau1.Phi(), tmass1, dm1));
			measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, tau2.Pt(), tau2.Eta(), tau2.Phi(), tmass2, dm2));
		}
		//######################################################
		
		if (cat.substr(2,2) == "ee"){
			measuredTauLeptons.push_back(MeasuredTauLepton(eleType1, tau3.Pt(), tau3.Eta(), tau3.Phi(), 0.000511));
			measuredTauLeptons.push_back(MeasuredTauLepton(eleType2, tau4.Pt(), tau4.Eta(), tau4.Phi(), 0.000511));
		}
		else if (cat.substr(2,2) == "em"){
			measuredTauLeptons.push_back(MeasuredTauLepton(eleType1, tau3.Pt(), tau3.Eta(), tau3.Phi(), 0.000511));
			measuredTauLeptons.push_back(MeasuredTauLepton(muType2, tau4.Pt(), tau4.Eta(), tau4.Phi(), 0.106));
		}
		else if (cat.substr(2,2) == "et"){
			measuredTauLeptons.push_back(MeasuredTauLepton(eleType1, tau3.Pt(), tau3.Eta(), tau3.Phi(), 0.000511));
			measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, tau4.Pt(), tau4.Eta(), tau4.Phi(), tmass4, dm4));
		}
		else if (cat.substr(2,2) == "mm"){
			measuredTauLeptons.push_back(MeasuredTauLepton(muType1, tau3.Pt(), tau3.Eta(), tau3.Phi(), 0.106));
			measuredTauLeptons.push_back(MeasuredTauLepton(muType2, tau4.Pt(), tau4.Eta(), tau4.Phi(), 0.106));
		}
		else if (cat.substr(2,2) == "mt"){
			measuredTauLeptons.push_back(MeasuredTauLepton(muType1, tau3.Pt(), tau3.Eta(), tau3.Phi(), 0.106));
			measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, tau4.Pt(), tau4.Eta(), tau4.Phi(), tmass4, dm4));
		}
		else if (cat.substr(2,2) == "tt"){
			measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, tau3.Pt(), tau3.Eta(), tau3.Phi(), tmass3, dm3));
			measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, tau4.Pt(), tau4.Eta(), tau4.Phi(), tmass4, dm4));
		}
		
		int verbosity = 0;
		ClassicSVfit4tau svFitAlgo_MarkovChain(ClassicSVfit4tau::kAlgoMarkovChain, verbosity);
		#ifdef USE_SVFITTF
		  //HadTauTFCrystalBall2* hadTauTF = new HadTauTFCrystalBall2();
		  //svFitAlgo_MarkovChain.setHadTauTF(hadTauTF);
		  //svFitAlgo_MarkovChain.enableHadTauTF();
		#endif
	  // run with mass constraint for each tau pair to match measured mass (125.06 GeV) of SM-like Higgs boson   
	  double massContraint = ((tau1+tau2).M()+(tau3+tau4).M())/2+ TMath::Sqrt(METx*METx+METy*METy);
	  //double kappa = 6.; // kappa parameter for log(M) term, cf. Eq. (41) in Nucl.Instrum.Meth. A862 (2017) 54-84
	  double kappa = 0.;
	  if ( kappa > 0. ) {
		svFitAlgo_MarkovChain.addLogM_fixed(true, kappa);
	  } else {
		svFitAlgo_MarkovChain.addLogM_fixed(false);
	  }
	  std::cout << "\n\nTesting Markov-Chain integration with ditau mass constraint set to " << massContraint << std::endl;
	  svFitAlgo_MarkovChain.setLikelihoodFileName("testClassicSVfit4tau_withMassContraint.root");
	  svFitAlgo_MarkovChain.setDiTau1MassConstraint(massContraint);
	  svFitAlgo_MarkovChain.setDiTau2MassConstraint(massContraint);
	  svFitAlgo_MarkovChain.integrate(measuredTauLeptons, METx, METy, covMET);
	  bool isValidSolution_1stRun = svFitAlgo_MarkovChain.isValidSolution();
	  //double dihiggs_mass_1stRun = svFitAlgo_MarkovChain.getMass();
	  //double dihiggs_massErr_1stRun = svFitAlgo_MarkovChain.getMassErr();
	  //double dihiggs_transverseMass_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->getTransverseMass();
	  //double dihiggs_transverseMassErr_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->getTransverseMassErr();
	  double ditau1_mass_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau1()->getMass();
	  double ditau1_massErr_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau1()->getMassErr();
	  double ditau2_mass_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau2()->getMass();
	  double ditau2_massErr_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau2()->getMassErr();
	  if ( isValidSolution_1stRun ) {
		std::cout << "(ditau1: mass = " << ditau1_mass_1stRun << " +/- " << ditau1_massErr_1stRun << ","
		          << " ditau2: mass = " << ditau2_mass_1stRun << " +/- " << ditau2_massErr_1stRun << ")" << std::endl;
	  } else {
		std::cout << "sorry, failed to find valid solution !!" << std::endl;
	  }
	
		// re-run without mass constraint for each tau pair;
		// this mode will set the mass of the second tau pair to match the mass of the first tau pair,
		// while the mass of the first tau pair is allowed to freely vary within the fit
	/*	svFitAlgo_MarkovChain.setLikelihoodFileName("testClassicSVfit4tau.root");
		svFitAlgo_MarkovChain.setDiTau1MassConstraint(-1.);
		svFitAlgo_MarkovChain.setDiTau2MassConstraint(-1.);
		svFitAlgo_MarkovChain.integrate(measuredTauLeptons, METx, METy, covMET);
		bool isValidSolution_2ndRun = svFitAlgo_MarkovChain.isValidSolution();
	
		printf("maxLikelihoood : %f \n", svFitAlgo_MarkovChain.getLmax());

		//double dihiggs_mass_2ndRun = svFitAlgo_MarkovChain.getMass();
		//double dihiggs_massErr_2ndRun = svFitAlgo_MarkovChain.getMassErr();
		//double dihiggs_transverseMass_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->getTransverseMass();
		//double dihiggs_transverseMassErr_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->getTransverseMassErr();
		ditau1_mass_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau1()->getMass();
		double ditau1_massErr_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau1()->getMassErr();
		ditau2_mass_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau2()->getMass();
		double ditau2_massErr_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau2()->getMassErr();
		if ( isValidSolution_2ndRun ) {
			std::cout << "(ditau1: mass = " << ditau1_mass_2ndRun << " +/- " << ditau1_massErr_2ndRun << ","
					<< " ditau2: mass = " << ditau2_mass_2ndRun << " +/- " << ditau2_massErr_2ndRun << ")" << std::endl;
		} else {
			std::cout << "sorry, failed to find valid solution !!" << std::endl;
		}*/
	//}//for
	struct mDCH_pair {double first; double second;};
	return mDCH_pair {ditau1_mass_1stRun, ditau2_mass_1stRun};
}


int main(int argc, char* argv[]){
	TFile *ifile = new TFile("MET_stuff/HppM1000_2018_1t.root","READ");
	TTree *tree = (TTree*)ifile->Get("Events");
	MyBranc(tree);
	TH1F* h_mll = new TH1F("h_mll", "mDCH1 Vis", 100, 0, 3000);
	TH1F* h_mll2 = new TH1F("h_mll2", "mDCH2 Vis", 100, 0, 3000);
	//TH1F* h_mDCH1_g = new TH1F("h_mDCH1_g", "mDCH1 MET split", 100, 0, 3000);
	//TH1F* h_mDCH2_g = new TH1F("h_mDCH2_g", "mDCH2 MET split", 100, 0, 3000);
	TH1F* h_mDCH1_sv = new TH1F("h_mDCH1_sv", "mDCH1 (MET split && FastMTT)", 100, 0, 3000);
	TH1F* h_mDCH2_sv = new TH1F("h_mDCH2_sv", "mDCH2 (MET split && FastMTT)", 100, 0, 3000);
	//TH1F* h_mDCH1_f = new TH1F("h_mDCH1_f", "mDCH1 FastMTT", 100, 0, 3000);
	//TH1F* h_mDCH2_f = new TH1F("h_mDCH2_f", "mDCH2 FastMTT", 100, 0, 3000);
	for (int i =0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		const char *cat_name = numberToCat(cat);	
		//int lep_count = cat_lepCount(cat_name,'e','m'); 
		if (strlen(cat_name) != 4) continue;
		//if (lep_count >2) continue;
		std::string	cat_str = cat_name;
		printf(" %s ; %f ; %f \n", cat_name, mll, mll2);
		
		h_mll->Fill(mll); h_mll2->Fill(mll2);

		// define MET covariance
		TMatrixD covMET(2, 2);
		covMET[0][0] = metcov00;
		covMET[1][0] = metcov10;
		covMET[0][1] = metcov01;
		covMET[1][1] = metcov11;
		//LVs of particles in CMS frame
		ROOT::Math::PtEtaPhiMVector MET(met, 0, metphi, 0);
		ROOT::Math::PtEtaPhiMVector Lep1(pt_1, eta_1, phi_1, m_1);//need to save m1
		ROOT::Math::PtEtaPhiMVector Lep2(pt_2, eta_2, phi_2, m_2);
		ROOT::Math::PtEtaPhiMVector Lep3(pt_3, eta_3, phi_3, m_3);
		ROOT::Math::PtEtaPhiMVector Lep4(pt_4, eta_4, phi_4, m_4);	
		
		auto [mDCH1_sv, mDCH2_sv] = runSVfit4tau(cat_str, MET.Px(), MET.Py(), covMET, Lep1, Lep2, Lep3, Lep4, decayMode_1, decayMode_2, decayMode_3, decayMode_4);
		
		//printf("dsdsd %f\t, %f\n",Lep4.M(), m_4);
		h_mDCH1_sv->Fill(mDCH1_sv); h_mDCH2_sv->Fill(mDCH2_sv);

	}//(event)for loop
	std::cout << std::endl;
	std::cout << "************************************************************************************************************************" << std::endl;
	std::cout << "* If you use this code, please cite:																																																																																		 *" << std::endl;
	std::cout << "*		 K. Ehataeht, L. Marzola, C. Veelken,																																																																														 *" << std::endl;
	std::cout << "*		 \"Reconstruction of the mass of Higgs boson pairs in events with Higgs boson pairs decaying into four tau leptons\", *" << std::endl;
	std::cout << "*		 arXiv:1809.06140																																																																																																		 *" << std::endl;
	std::cout << "************************************************************************************************************************" << std::endl;
	std::cout << std::endl;
	
	TFile *ofile = new TFile("MET_SVfit4tau.root", "RECREATE");
	if (!ofile || ofile->IsZombie()) {
		std::cerr << "Error: Could not open the output file " << std::endl;
		return 1; // Exit with an error code
	}
	h_mll->Write(); h_mll2->Write();
	//h_mDCH1_g->Write();h_mDCH2_g->Write();
	//h_mDCH1_f->Write(); h_mDCH2_f->Write();
	h_mDCH1_sv->Write(); h_mDCH2_sv->Write();
	ofile->Close(); 
	return 0;
}
