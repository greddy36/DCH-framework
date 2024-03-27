
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

	//ROOT::Math::PtEtaPhiMVector vis_total = DCH1+DCH2+MET;
	// Define a boost vector along the z-axis with a beta such that Pz = 0
	double boostZ = 0;//-vis_total.Pz()/vis_total.E();
	ROOT::Math::Boost boostVec(0, 0, boostZ);
	//boostVec(vis_total);
	//ROOT::Math::PtEtaPhiMVector dch1_b = boostVec(DCH1);
	//ROOT::Math::PtEtaPhiMVector dch2_b = boostVec(DCH2);
	ROOT::Math::PtEtaPhiMVector met_b = boostVec(MET);

	//dch1_b and dch2_b form our new basis (_h). Basis transformation matrix is given by...
	/*TMatrixD basis_h(2,2);
	//double basis_arr[4]={dch1_b.Pt()*TMath::Cos(dch1_b.Phi()), dch2_b.Pt()*TMath::Cos(dch2_b.Phi()), dch1_b.Pt()*TMath::Sin(dch1_b.Phi()), dch2_b.Pt()*TMath::Sin(dch2_b.Phi())};
	double basis_arr[4]={TMath::Cos(dch1_b.Phi()), TMath::Cos(dch2_b.Phi()), TMath::Sin(dch1_b.Phi()), TMath::Sin(dch2_b.Phi())};
	basis_h.SetMatrixArray(basis_arr); basis_h.Invert();
	// Now we seperate the Azimuthal(transverse) component of boosted met_b and hopefully this is our promptmet_b in boosted frame
	TMatrixD met_b_xy(2,1), total_neutr(2,1);
	double met_b_arr[2] = {met_b.Pt()*TMath::Cos(met_b.Phi()), met_b.Pt()*TMath::Sin(met_b.Phi())};
	met_b_xy.SetMatrixArray(met_b_arr);//transverse part of boosted met_b
	total_neutr = basis_h * met_b_xy; 
	//total_neutr.Print();*/

	float factor = 1, tmp = 1;
	//if (abs(mll-mll2) > (mll+mll2)/2){
	if (mll != mll2){
		factor = mll/(mll+mll2); //abs(mll-mll2)/abs(total_neutr(0,0)-total_neutr(1,0));
	}
	if (mll > mll2) tmp = 1-factor;
	else tmp = factor;
		ROOT::Math::PtEtaPhiMVector	neutr_leg1(met_b.Pt()*tmp, 0, met_b.Phi(), 0);
		ROOT::Math::PtEtaPhiMVector	neutr_leg2(met_b.Pt()*(1-tmp), 0, met_b.Phi(), 0);
	
	/*double sh = -TMath::SinH(dch1_b.Eta())/TMath::SinH(dch2_b.Eta());
	double nu_pt1 = met_b.Pt()/(TMath::Sqrt(1+sh*sh+2*sh*TMath::Cos(dch1_b.Phi()-dch2_b.Phi())));// split in Transverse DCH basis
	double nu_pt2 = sh*nu_pt1;
	ROOT::Math::PtEtaPhiMVector	neutr_leg1(nu_pt1, dch1_b.Eta(), dch1_b.Phi(), abs(mll-mll2)*tmp);
	ROOT::Math::PtEtaPhiMVector neutr_leg2(nu_pt2, dch2_b.Eta(), dch2_b.Phi(), abs(mll-mll2)*(1-tmp));*/
	//now we need to un-boost these neutrino legs.
	ROOT::Math::Boost unboostVec(0, 0, -boostZ);
	
	struct neut_legs {ROOT::Math::PtEtaPhiMVector first; ROOT::Math::PtEtaPhiMVector second;};
	return neut_legs {unboostVec(neutr_leg1), unboostVec(neutr_leg2)};
}

double tmass(int DM){
	double mass = 1.0;
	if (DM == 0) mass = 0.1396;
	return mass;
}  

LorentzVector runFTT(std::string dch, double METx, double METy, TMatrixD covMET, ROOT::Math::PtEtaPhiMVector lep1_p4, ROOT::Math::PtEtaPhiMVector lep2_p4, double dm1, double dm2 ){
	// define lepton four vectors
	std::vector<MeasuredTauLepton> measuredTauLeptons;
	//std::cout<<dch<< "\t"<<0.1395<<std::endl;
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
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, lep2_p4.Pt(), lep2_p4.Eta(), lep2_p4.Phi(),  tmass(dm2),dm2));
	}
	else if (dch == "mm"){
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay, lep1_p4.Pt(), lep1_p4.Eta(), lep1_p4.Phi(), 0.106));
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay, lep2_p4.Pt(), lep2_p4.Eta(), lep2_p4.Phi(), 0.106));
	}
	else if (dch == "mt"){
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay, lep1_p4.Pt(), lep1_p4.Eta(), lep1_p4.Phi(), 0.106));
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, lep2_p4.Pt(), lep2_p4.Eta(), lep2_p4.Phi(),  tmass(dm2),dm2));
	}
	else if (dch == "tt"){
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, lep1_p4.Pt(), lep1_p4.Eta(), lep1_p4.Phi(),  tmass(dm1),dm1));
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, lep2_p4.Pt(), lep2_p4.Eta(), lep2_p4.Phi(), tmass(dm2),dm2));
	}
	FastMTT FMTT;
	FMTT.run(measuredTauLeptons, METx, METy, covMET);
	LorentzVector	DCH_4vec = FMTT.getBestP4();
	return DCH_4vec;
}

double runSVfit(std::string dch, double METx, double METy, TMatrixD covMET, ROOT::Math::PtEtaPhiMVector lep1_p4, ROOT::Math::PtEtaPhiMVector lep2_p4, double dm1, double dm2 ){
	// define lepton four vectors
	std::vector<MeasuredTauLepton> measuredTauLeptons;
	//std::cout<<dch<< "\t"<<0.1395<<std::endl;
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
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, lep2_p4.Pt(), lep2_p4.Eta(), lep2_p4.Phi(), tmass(dm2),dm2));
	}
	else if (dch == "mm"){
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay, lep1_p4.Pt(), lep1_p4.Eta(), lep1_p4.Phi(), 0.106));
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay, lep2_p4.Pt(), lep2_p4.Eta(), lep2_p4.Phi(), 0.106));
	}
	else if (dch == "mt"){
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay, lep1_p4.Pt(), lep1_p4.Eta(), lep1_p4.Phi(), 0.106));
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, lep2_p4.Pt(), lep2_p4.Eta(), lep2_p4.Phi(), tmass(dm2),dm2));
	}
	else if (dch == "tt"){
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, lep1_p4.Pt(), lep1_p4.Eta(), lep1_p4.Phi(), tmass(dm1),dm1));
					measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, lep2_p4.Pt(), lep2_p4.Eta(), lep2_p4.Phi(), tmass(dm2),dm2));
	}
	ClassicSVfit SVfit;
	SVfit.integrate(measuredTauLeptons, METx, METy, covMET);
	double mDCH = static_cast<HistogramAdapterDiTau*>(SVfit.getHistogramAdapter())->getMass();
	//double massErr = static_cast<HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getMassErr();
	return mDCH;
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
		if (lep_count !=4) continue;
		std::string	cat_str = cat_name;
		//std::cout<< cat_name <<std::endl;
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
		ROOT::Math::PtEtaPhiMVector DCH1 = Lep1+Lep2;
		ROOT::Math::PtEtaPhiMVector DCH2 = Lep3+Lep4;
			
		double mDCH1_f = runSVfit(cat_str.substr(0,2), MET.Px(), MET.Py(), covMET, Lep1, Lep2, decayMode_1, decayMode_2);
		double mDCH2_f = runSVfit(cat_str.substr(2,2), MET.Px(), MET.Py(), covMET, Lep3, Lep4,decayMode_3, decayMode_4);
		h_mDCH1_f->Fill(mDCH1_f); h_mDCH2_f->Fill(mDCH2_f);
		
		
		auto [neut_leg1, neut_leg2] = MET_split(DCH1, DCH2, MET);
		h_mDCH1_g->Fill((DCH1+neut_leg1).M()); h_mDCH2_g->Fill((DCH2+neut_leg2).M());

		//boost each of the neutrino legs to transverse frame before giving to FTT
		//std::cout<< (DCH2+neut_leg2).M() <<std::endl;		
		/*double boost1 = -neut_leg1.Pz()/neut_leg1.E();
		double boost2 = -neut_leg2.Pz()/neut_leg2.E();
		ROOT::Math::Boost boost_leg1(0, 0, boost1);
		ROOT::Math::Boost boost_leg2(0, 0, boost2);
		ROOT::Math::Boost boost_leg1(0, 0, 0);
		ROOT::Math::Boost boost_leg2(0, 0, 0);
		ROOT::Math::PtEtaPhiMVector met1 = boost_leg1(neut_leg1); 
		ROOT::Math::PtEtaPhiMVector met2 = boost_leg2(neut_leg2);*/
		//std::cout<<"boosts: "<<boost1<<"\t"<<boost2<<std::endl; 
		//std::cout<< met1.Pz() <<std::endl;		
		
		double mDCH1_sv = runSVfit(cat_str.substr(0,2), neut_leg1.Px(), neut_leg1.Py(), covMET, Lep1, Lep2,decayMode_1, decayMode_2);
		double mDCH2_sv = runSVfit(cat_str.substr(2,2), neut_leg2.Px(), neut_leg2.Py(), covMET, Lep3, Lep4,decayMode_3, decayMode_4);
		//printf("dsdsd %f\t, %f\n",Lep4.M(), m_4);
		h_mDCH1_sv->Fill(mDCH1_sv); h_mDCH2_sv->Fill(mDCH2_sv);
		


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
