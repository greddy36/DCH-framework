#include "TMath.h"
#include <cmath>
#include <vector>
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMatrix.h"
#include "MyBranch.C"
#include "cat.h"
void MET_incl() {
TCanvas *can= new TCanvas("can","can",700,500);gStyle->SetOptStat(0); 
TLegend *leg = new TLegend(0.7, 0.7, .9, .9);
	TFile *ifile = new TFile("HppM1000_2018.root","READ");
	TTree *tree = (TTree*)ifile->Get("Events");
	MyBranch(tree);
	TFile* ofile = new TFile("MET_exercise.root", "RECREATE");
	TH1F* h_mll = new TH1F("h_mll", "mll", 100, 0, 3000);
	TH1F* h_mll2 = new TH1F("h_mll2", "mll2", 100, 0, 3000);
	TH1F* h_mDCH1 = new TH1F("h_mDCH1", "mDCH1", 100, 0, 3000);
	TH1F* h_mDCH2 = new TH1F("h_mDCH2", "mDCH2", 100, 0, 3000);
	for (int i =0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		const char *cat_name = numberToCat(cat);
		int lep_count = cat_lepCount(cat_name,'e','m'); 
		if (strlen(cat_name) != 4) continue;
		if (lep_count !=3) continue;
		//if (cat !=1) continue; 
		TLorentzVector Lep1, Lep2, Lep3, Lep4, MET, DCH1, DCH2; //TLV for CMS frame particles
		MET.SetPtEtaPhiM(met, 0, metphi, 0);
		Lep1.SetPtEtaPhiM(pt_1, eta_1, phi_1, m_1_tr);//need to save m1
		Lep2.SetPtEtaPhiM(pt_2, eta_2, phi_2, m_2_tr);
		Lep3.SetPtEtaPhiM(pt_3, eta_3, phi_3, m_3);
		Lep4.SetPtEtaPhiM(pt_4, eta_4, phi_4, m_4);	
		DCH1 = Lep1+Lep2;
		DCH2 = Lep3+Lep4;
		TLorentzVector vis_total = DCH1+DCH2; 
		// Define a boost vector along the z-axis with a beta such that Pz = 0
		TVector3 boostVec(0, 0, -vis_total.Pz()/vis_total.E());
		//TVector3 boostVec(0, 0, 0);
		vis_total.Boost(boostVec); 
		DCH1.Boost(boostVec); 
		DCH2.Boost(boostVec); 
		MET.Boost(boostVec); 
		//cout<<MET.Pt()<<"\t"<<MET.Pz()<<endl;

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
		TLorentzVector neutr_leg1, neutr_leg2;
		float factor = 1, tmp = 1;
		if (abs(mll-mll2) > (mll+mll2)/2){
		//if (mll != mll2){
			factor = MET.Pt()*mll/(mll+mll2); //abs(mll-mll2)/abs(total_neutr(0,0)-total_neutr(1,0));
		}
		if (mll > mll2) tmp = 1-factor;
		else tmp = factor;

		/*neutr_leg1.SetPtEtaPhiM(MET.Pt()*tmp, DCH1.Eta(), DCH1.Phi(), abs(mll-mll2)*tmp);//in Transverse DCH basis
		neutr_leg2.SetPtEtaPhiM(MET.Pt()*(1-tmp), DCH2.Eta(), DCH2.Phi(),abs(mll-mll2)*(1-tmp));//in Transverse DCH basis*/
		
		double sh = TMath::SinH(DCH1.Eta())/TMath::SinH(DCH2.Eta());
		double nu_pt1 = MET.Pt()/(TMath::Sqrt(1+sh*sh+2*sh*TMath::Cos(DCH1.Phi()-DCH2.Phi())));
		double nu_pt2 = sh*nu_pt1;
		neutr_leg1.SetPtEtaPhiM(nu_pt1, DCH1.Eta(), DCH1.Phi(), abs(mll-mll2)*tmp);//in Transverse DCH basis
		neutr_leg2.SetPtEtaPhiM(nu_pt2, DCH2.Eta(), DCH2.Phi(), abs(mll-mll2)*(1-tmp));//in Transverse DCH basis
		
		
		
		/*neutr_leg1.SetPtEtaPhiM(total_neutr(0,0), TMath::ASinH(-DCH1.Pz()/total_neutr(0,0)), DCH1.Phi(),0);//in Transverse DCH basis
		neutr_leg2.SetPtEtaPhiM(total_neutr(1,0), TMath::ASinH(-DCH2.Pz()/total_neutr(1,0)), DCH2.Phi(),0);//in Transverse DCH basis*/
		/*##########################################
		Assumptions: 
		1) neutr_leg1 and neutr_leg2 are the full neutrino vectors of each DCH leg. 
		2) neutr_legs are massless. 
		3) neutrinos are in the same direction as DCH
	         ##########################################*/     
	        /*double C = (pow(mll,2)-pow(mll2,2))/2, A = DCH1.E()-DCH1.P(), B = DCH2.E()-DCH2.P() ;
	        TVector3 dummy(MET.Px(),MET.Py(),MET.Pz());
	        double cosx = TMath::Cos( DCH1.Angle(dummy) );
	        double var1 = 2*A*C-B*B*MET.P()*cosx, var2 = A*A-B*B, var3 = B*B*MET.P()*MET.P()-C*C;
	        double P1_met = (var1+TMath::Sqrt(4*var2*var3+var1*var1))/(2*var2);
	        double P2_met = TMath::Sqrt(MET.P()*MET.P()+P1_met*P1_met-2*MET.P()*P1_met*cosx);
	        neutr_leg1.SetPtEtaPhiM(P1_met*1/TMath::CosH(DCH1.Eta()), DCH1.Eta(), DCH1.Phi(),0);//in Transverse DCH basis
		neutr_leg2.SetPtEtaPhiM(P2_met*1/TMath::CosH(DCH2.Eta()), DCH2.Eta(), DCH2.Phi(),0);*/

		
		//##########################################
		cout<<neutr_leg1.Pz()<<"\t"<<neutr_leg2.Pz()<<endl;
		//cout<<neutr_leg1.Eta()<<"\t"<<neutr_leg2.Eta()<<endl;
		cout<< (neutr_leg1+neutr_leg2).Pz() <<"\t"<< (MET).Pz() <<endl;
		
		h_mll->Fill(mll); h_mll2->Fill(mll2);
		h_mDCH1->Fill((DCH1+neutr_leg1).M());
		h_mDCH2->Fill((DCH2+neutr_leg2).M());
	}
	h_mll->Write();
	h_mll2->Write();
	h_mDCH1->Write();
	h_mDCH2->Write();
	h_mll->SetLineColor(4);h_mll2->SetLineColor(4);h_mDCH1->SetLineColor(2);h_mDCH2->SetLineColor(2);
	h_mll->SetLineWidth(2);h_mll2->SetLineWidth(2);h_mDCH1->SetLineWidth(2);h_mDCH2->SetLineWidth(2);
	h_mDCH1->Draw();h_mll->Draw("same");leg->AddEntry(h_mll, "no MET");leg->AddEntry(h_mDCH1, "w/ MET split"); 
	h_mDCH2->Draw(); h_mll2->Draw("same"); leg->AddEntry(h_mll2, "no MET");leg->AddEntry(h_mDCH2, "w/ MET split"); 

	leg->Draw();
	can->SaveAs("0mll.png");	
}

