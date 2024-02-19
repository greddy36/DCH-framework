
/**
	 \class testClassicSVfitLFV testClassicSVfitLFV .cc "TauAnalysis/ClassicSVfit/bin/testClassicSVfitLFV .cc"
	 \brief Example for computing Higgs boson (H) mass in lepton-flavor-violating H decays, e.g H -> mu tau -> mu tau_h nu (cf. arXiv:1502.07400)
*/

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // LorentzVector
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"
#include "../../../MET_stuff/MyBranc.C"
#include "../../../MET_stuff/cat.h"

#include "Math/GenVector/Boost.h"

#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include "TH1F.h"
#include "TMath.h"
#include <string>
using namespace classic_svFit;

auto MET_split(ROOT::Math::PtEtaPhiMVector DCH1, ROOT::Math::PtEtaPhiMVector DCH2, ROOT::Math::PtEtaPhiMVector MET){

	ROOT::Math::PtEtaPhiMVector vis_total = DCH1+DCH2;
	// Define a boost vector along the z-axis with a beta such that Pz = 0
	double boostZ = -vis_total.Pz()/vis_total.E();
	ROOT::Math::Boost boostVec(0, 0, boostZ);
	boostVec(vis_total);
	boostVec(DCH1);
	boostVec(DCH2);
	boostVec(MET);

	//DCH1 and DCH2 form our new basis (_h). Basis transformation matrix is given by...
	TMatrixD basis_h(2,2);
	//double basis_arr[4]={DCH1.Pt()*TMath::Cos(DCH1.Phi()), DCH2.Pt()*TMath::Cos(DCH2.Phi()), DCH1.Pt()*TMath::Sin(DCH1.Phi()), DCH2.Pt()*TMath::Sin(DCH2.Phi())};
	double basis_arr[4]={TMath::Cos(DCH1.Phi()), TMath::Cos(DCH2.Phi()), TMath::Sin(DCH1.Phi()), TMath::Sin(DCH2.Phi())};
	basis_h.SetMatrixArray(basis_arr); basis_h.Invert();
	// Now we seperate the Azimuthal(transverse) component of boosted MET and hopefully this is our promptMET in boosted frame
	TMatrixD MET_xy(2,1), total_neutr(2,1);
	double MET_arr[2] = {MET.Pt()*TMath::Cos(MET.Phi()), MET.Pt()*TMath::Sin(MET.Phi())};
	MET_xy.SetMatrixArray(MET_arr);//transverse part of boosted MET
	total_neutr = basis_h * MET_xy; 
	//total_neutr.Print();

	float factor = 1, tmp = 1;
	if (abs(mll-mll2) > (mll+mll2)/2){
	//if (mll != mll2){
		factor = MET.Pt()*mll/(mll+mll2); //abs(mll-mll2)/abs(total_neutr(0,0)-total_neutr(1,0));
	}
	if (mll > mll2) tmp = 1-factor;
	else tmp = factor;
		/*neutr_leg1(MET.Pt()*tmp, DCH1.Eta(), DCH1.Phi(), abs(mll-mll2)*tmp);//in Transverse DCH basis
	neutr_leg2(MET.Pt()*(1-tmp), DCH2.Eta(), DCH2.Phi(),abs(mll-mll2)*(1-tmp));//in Transverse DCH basis*/
	
	double sh = TMath::SinH(DCH1.Eta())/TMath::SinH(DCH2.Eta());
	double nu_pt1 = MET.Pt()/(TMath::Sqrt(1+sh*sh+2*sh*TMath::Cos(DCH1.Phi()-DCH2.Phi())));
	double nu_pt2 = sh*nu_pt1;
	ROOT::Math::PtEtaPhiMVector	neutr_leg1(nu_pt1, DCH1.Eta(), DCH1.Phi(), abs(mll-mll2)*tmp);
	ROOT::Math::PtEtaPhiMVector neutr_leg2(nu_pt2, DCH2.Eta(), DCH2.Phi(), abs(mll-mll2)*(1-tmp));
	//now we need to un-boost these neutrino legs.
	ROOT::Math::Boost unboostVec(0, 0, -boostZ);
	unboostVec(neutr_leg1); unboostVec(neutr_leg2);
	
	struct neut_legs {ROOT::Math::PtEtaPhiMVector first; ROOT::Math::PtEtaPhiMVector second;};
	return neut_legs {neutr_leg1, neutr_leg2};
}


LorentzVector runFTT(std::string dch, double METx, double METy, TMatrixD covMET, ROOT::Math::PtEtaPhiMVector lep1_p4, ROOT::Math::PtEtaPhiMVector lep2_p4 ){
	// define lepton four vectors
	std::vector<MeasuredTauLepton> measuredTauLeptons;
	
	if (dch == "ee"){
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToElecDecay, lep1_p4.Pt(), lep1_p4.Eta(), lep1_p4.Phi(), 0.000511));
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToElecDecay, lep2_p4.Pt(), lep2_p4.Eta(), lep2_p4.Phi(), 0.000511));
	}
	else if (dch == "em"){
						measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToElecDecay, lep1_p4.Pt(), lep1_p4.Eta(), lep1_p4.Phi(), 0.000511));
										measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay, lep2_p4.Pt(), lep2_p4.Eta(), lep2_p4.Phi(), 0.106));
	}
	else if (dch == "et"){
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToElecDecay, lep1_p4.Pt(), lep1_p4.Eta(), lep1_p4.Phi(), 0.000511));
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, lep2_p4.Pt(), lep2_p4.Eta(), lep2_p4.Phi(), lep2_p4.M()));
	}
	else if (dch == "mm"){
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay, lep1_p4.Pt(), lep1_p4.Eta(), lep1_p4.Phi(), 0.106));
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay, lep2_p4.Pt(), lep2_p4.Eta(), lep2_p4.Phi(), 0.106));
	}
	else if (dch == "mt"){
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay, lep1_p4.Pt(), lep1_p4.Eta(), lep1_p4.Phi(), 0.106));
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, lep2_p4.Pt(), lep2_p4.Eta(), lep2_p4.Phi(), lep2_p4.M()));
	}
	else if (dch == "tt"){
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, lep1_p4.Pt(), lep1_p4.Eta(), lep1_p4.Phi(), lep1_p4.M()));
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, lep2_p4.Pt(), lep2_p4.Eta(), lep2_p4.Phi(), lep2_p4.M()));
	}
	FastMTT FMTT;
	FMTT.run(measuredTauLeptons, METx, METy, covMET);
	LorentzVector	DCH_4vec = FMTT.getBestP4();
	return DCH_4vec;
}


int main(int argc, char* argv[]){

	TFile *ifile = new TFile("MET_stuff/HppM1000_2018.root","READ");
	TTree *tree = (TTree*)ifile->Get("Events");
	MyBranc(tree);
	TH1F* h_mll = new TH1F("h_mll", "mDCH1 Vis", 100, 0, 3000);
	TH1F* h_mll2 = new TH1F("h_mll2", "mDCH2 Vis", 100, 0, 3000);
	TH1F* h_mDCH1_g = new TH1F("h_mDCH1_g", "mDCH1 MET split", 100, 0, 3000);
	TH1F* h_mDCH2_g = new TH1F("h_mDCH2_g", "mDCH2 MET split", 100, 0, 3000);
	TH1F* h_mDCH1_sv = new TH1F("h_mDCH1_sv", "mDCH1 (MET split && FastMTT)", 100, 0, 3000);
	TH1F* h_mDCH2_sv = new TH1F("h_mDCH2_sv", "mDCH2 (MET split && FastMTT)", 100, 0, 3000);
	TH1F* h_mDCH1_f = new TH1F("h_mDCH1_f", "mDCH1 FastMTT", 100, 0, 3000);
	TH1F* h_mDCH2_f = new TH1F("h_mDCH2_f", "mDCH2 FastMTT", 100, 0, 3000);
	for (int i =0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		const char *cat_name = numberToCat(cat);	
		int lep_count = cat_lepCount(cat_name,'e','m'); 
		if (strlen(cat_name) != 4) continue;
		if (lep_count !=3) continue;
		std::string	cat_str = cat_name;
		std::cout<< cat_name <<std::endl;
		h_mll->Fill(mll); h_mll2->Fill(mll2);

		// define MET covariance
		TMatrixD covMET(2, 2);
		covMET[0][0] = metcov00;
		covMET[1][0] = metcov10;
		covMET[0][1] = metcov01;
		covMET[1][1] = metcov11;

		//LVs of particles in CMS frame
		ROOT::Math::PtEtaPhiMVector MET(met, 0, metphi, 0);
		ROOT::Math::PtEtaPhiMVector Lep1(pt_1, eta_1, phi_1, m_1_tr);//need to save m1
		ROOT::Math::PtEtaPhiMVector Lep2(pt_2, eta_2, phi_2, m_2_tr);
		ROOT::Math::PtEtaPhiMVector Lep3(pt_3, eta_3, phi_3, m_3);
		ROOT::Math::PtEtaPhiMVector Lep4(pt_4, eta_4, phi_4, m_4);	
		ROOT::Math::PtEtaPhiMVector DCH1 = Lep1+Lep2;
		ROOT::Math::PtEtaPhiMVector DCH2 = Lep3+Lep4;
			
		LorentzVector DCH1_f = runFTT(cat_str.substr(0,2), MET.Px(), MET.Py(), covMET, Lep1, Lep2);
		LorentzVector DCH2_f = runFTT(cat_str.substr(2,2), MET.Px(), MET.Py(), covMET, Lep3, Lep4);
		h_mDCH1_f->Fill(DCH1_f.M()); h_mDCH2_f->Fill(DCH2_f.M());		
		
		
		auto [neut_leg1, neut_leg2] = MET_split(DCH1, DCH2, MET);
		h_mDCH1_g->Fill((DCH1+neut_leg1).M()); h_mDCH2_g->Fill((DCH2+neut_leg2).M());

		//boost each of the neutrino legs to transverse frame before giving to FTT
		std::cout<< neut_leg1.Pz() <<std::endl;		
		double boost1 = -neut_leg1.Pz()/neut_leg1.E();
		double boost2 = -neut_leg2.Pz()/neut_leg2.E();
		ROOT::Math::Boost boost_leg1(0, 0, boost1);
		ROOT::Math::Boost boost_leg2(0, 0, boost2);
		ROOT::Math::PtEtaPhiMVector met1 = boost_leg1(neut_leg1); 
		ROOT::Math::PtEtaPhiMVector met2 = boost_leg2(neut_leg2);
		std::cout<< met1.Pz() <<std::endl;		
		
		LorentzVector DCH1_sv = runFTT(cat_str.substr(0,2), met1.Px(), met1.Py(), covMET, boost_leg1(Lep1), boost_leg1(Lep2));
		LorentzVector DCH2_sv = runFTT(cat_str.substr(2,2), met2.Px(), met2.Py(), covMET, boost_leg2(Lep3), boost_leg2(Lep4));
		h_mDCH1_sv->Fill(DCH1_sv.M()); h_mDCH2_sv->Fill(DCH2_sv.M());
		

	}//evt loop	
	TFile *ofile = new TFile("MET_SVfit.root", "RECREATE");
	if (!ofile || ofile->IsZombie()) {
		std::cerr << "Error: Could not open the output file " << std::endl;
		return 1; // Exit with an error code
	}
	h_mll->Write(); h_mll2->Write();
	h_mDCH1_g->Write();h_mDCH2_g->Write();
	h_mDCH1_f->Write(); h_mDCH2_f->Write();
	h_mDCH1_sv->Write(); h_mDCH2_sv->Write();
	ofile->Close(); 
	return 0;
}
