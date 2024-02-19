#include "Kinematics.h"
#include "cat.h" //cat sort funs

float* SortPt(string cat, string tau){
	float pt[] = {pt_1,pt_2,pt_3,pt_4};
	int n = cat_lepCount(cat, 'e','m');
	std::vector<float> sort_pt;
	if (tau=="t"){
		for (int j = 0; j < cat.length(); j++){
			if (cat[j] == 't')
				sort_pt.push_back(pt[j]);
		}
	}	
	else{
		for (int j = 0; j < cat.length(); j++){
			if (cat[j] != 't')
				sort_pt.push_back(pt[j]);
		}
	}	
	float *pt_arr = new float[0]; //any number to initialize it.
	std::sort(sort_pt.begin(), sort_pt.end(), greater());//descending order
	std::copy(sort_pt.begin(), sort_pt.begin()+n, pt_arr);//vector->array	
	return pt_arr;
}
	
float ST(string cat){//for only leptons
	float st = 0;
	if (cat.find("e")+1==1 or cat.find("m")+1==1)
		st += pt_1;
	if (cat.find("e")+1==2 or cat.find("m")+1==2)
		st += pt_2; 
	if (cat.find("e")+1==3 or cat.find("m")+1==3)
		st += pt_3;
	if (cat.find("e")+1==4 or cat.find("m")+1==4)
		st += pt_4; 
	return st;
}
float getDR(float eta1, float phi1, float eta2, float phi2) {
    float pi = TMath::Pi();
    float dPhi = fmin(fabs(phi2 - phi1), 2.0 * pi - fabs(phi2 - phi1));
    float DR = sqrt(pow(dPhi, 2) + pow(eta2 - eta1, 2));
    return DR;
}

TLorentzVector LepV(int n){
	TLorentzVector lepV;
	if (n==1)
		lepV.SetPtEtaPhiM(pt_1, eta_1, phi_1, m_uncor_1);
	if (n==2)
		lepV.SetPtEtaPhiM(pt_2, eta_2, phi_2, m_uncor_1);
	if (n==3)
		lepV.SetPtEtaPhiM(pt_3, eta_3, phi_3, m_3);
	if (n==4)
		lepV.SetPtEtaPhiM(pt_4, eta_4, phi_4, m_4);
	return lepV;
}			
