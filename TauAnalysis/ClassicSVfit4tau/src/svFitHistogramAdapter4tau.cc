#include "TauAnalysis/ClassicSVfit4tau/interface/svFitHistogramAdapter4tau.h"

#include <TMath.h>
#include <TLorentzVector.h>

#include <numeric>

#include <boost/algorithm/string/replace.hpp>

using namespace classic_svFit;

//-------------------------------------------------------------------------------------------------
// auxiliary classes to reconstruct pT, eta, phi, mass, and transverse mass of Higgs boson pairs (built from 4 taus)
SVfitQuantityDiHiggs::SVfitQuantityDiHiggs(const std::string& label)
  : SVfitQuantity(label)
{}

void SVfitQuantityDiHiggs::bookHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met)
{
  if ( histogram_ != nullptr ) delete histogram_;
  histogram_ = createHistogram(vis1P4, vis2P4, vis3P4, vis4P4, met);
  histogram_->SetName(std::string(histogram_->GetName() + uniqueName_).c_str());
}

SVfitQuantityDiHiggsPt::SVfitQuantityDiHiggsPt(const std::string& label)
  : SVfitQuantityDiHiggs(label)
{}

TH1* SVfitQuantityDiHiggsPt::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met) const
{
  return HistogramTools::makeHistogram_logBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramPt", 1., 1.e+3, 1.025);
}

SVfitQuantityDiHiggsEta::SVfitQuantityDiHiggsEta(const std::string& label)
  : SVfitQuantityDiHiggs(label)
{}

TH1* SVfitQuantityDiHiggsEta::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met) const
{
  return HistogramTools::makeHistogram_linBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramEta", 198, -9.9, +9.9);
}

SVfitQuantityDiHiggsPhi::SVfitQuantityDiHiggsPhi(const std::string& label)
  : SVfitQuantityDiHiggs(label)
{}

TH1* SVfitQuantityDiHiggsPhi::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met) const
{
  return HistogramTools::makeHistogram_linBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramPhi", 180, -TMath::Pi(), +TMath::Pi());
}

SVfitQuantityDiHiggsMass::SVfitQuantityDiHiggsMass(const std::string& label)
  : SVfitQuantityDiHiggs(label)
{}

TH1* SVfitQuantityDiHiggsMass::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met) const
{
  double visMass = (vis1P4 + vis2P4 + vis3P4 + vis4P4).mass();
  double minMass = visMass/1.0125;
  double maxMass = TMath::Max(1.e+4, 1.e+1*minMass);
  return HistogramTools::makeHistogram_logBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramMass", minMass, maxMass, 1.025);
}

SVfitQuantityDiHiggsTransverseMass::SVfitQuantityDiHiggsTransverseMass(const std::string& label)
  : SVfitQuantityDiHiggs(label)
{}

TH1* SVfitQuantityDiHiggsTransverseMass::createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met) const
{
  classic_svFit::LorentzVector measuredDiTauSystem1 = vis1P4 + vis2P4;
  classic_svFit::LorentzVector measuredDiTauSystem2 = vis3P4 + vis4P4;
  classic_svFit::LorentzVector measuredDiHiggsSystem = measuredDiTauSystem1 + measuredDiTauSystem2;
  double visTransverseMass2 = square(measuredDiTauSystem1.Et() + measuredDiTauSystem2.Et()) - (square(measuredDiHiggsSystem.px()) + square(measuredDiHiggsSystem.py()));
  double visTransverseMass = TMath::Sqrt(TMath::Max(1., visTransverseMass2));
  double minTransverseMass = visTransverseMass/1.0125;
  double maxTransverseMass = TMath::Max(1.e+4, 1.e+1*minTransverseMass);
  return HistogramTools::makeHistogram_logBinWidth("ClassicSVfitIntegrand_" + label_ + "_histogramTransverseMass", minTransverseMass, maxTransverseMass, 1.025);
}
    
HistogramAdapterDiHiggs::HistogramAdapterDiHiggs(const std::string& label)
  : HistogramAdapter(label)
  , quantity_pt_(nullptr)
  , quantity_eta_(nullptr)
  , quantity_phi_(nullptr)
  , quantity_mass_(nullptr)
  , quantity_transverseMass_(nullptr)
  , adapter_ditau1_(nullptr)
  , adapter_ditau2_(nullptr)
{
  quantity_pt_ = new SVfitQuantityDiHiggsPt(label_);
  quantities_.push_back(quantity_pt_);
  quantity_eta_ = new SVfitQuantityDiHiggsEta(label_);
  quantities_.push_back(quantity_eta_);
  quantity_phi_ = new SVfitQuantityDiHiggsPhi(label_);
  quantities_.push_back(quantity_phi_);
  quantity_mass_ = new SVfitQuantityDiHiggsMass(label_);
  quantities_.push_back(quantity_mass_);
  quantity_transverseMass_ = new SVfitQuantityDiHiggsTransverseMass(label_);
  quantities_.push_back(quantity_transverseMass_);
  adapter_ditau1_ = new HistogramAdapterDiTau(label_ + "_ditau1");
  adapter_ditau2_ = new HistogramAdapterDiTau(label_ + "_ditau2");
}

HistogramAdapterDiHiggs::~HistogramAdapterDiHiggs()
{
  delete adapter_ditau1_;
  delete adapter_ditau2_;
}

void HistogramAdapterDiHiggs::setMeasurement(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met)
{
  vis1P4_ = vis1P4;
  vis2P4_ = vis2P4;
  vis3P4_ = vis3P4;
  vis4P4_ = vis4P4;
  met_ = met;
  adapter_ditau1_->setMeasurement(vis1P4, vis2P4, met);
  adapter_ditau2_->setMeasurement(vis3P4, vis4P4, met);
}

void HistogramAdapterDiHiggs::setTau123And4P4(const LorentzVector& tau1P4, const LorentzVector& tau2P4, const LorentzVector& tau3P4, const LorentzVector& tau4P4) 
{
  tau1P4_ = tau1P4;
  tau2P4_ = tau2P4;
  ditau1P4_ = tau1P4_ + tau2P4_;
  tau3P4_ = tau3P4;
  tau4P4_ = tau4P4;
  ditau2P4_ = tau3P4_ + tau4P4_;
  dihiggsP4_ = ditau1P4_ + ditau2P4_;
  adapter_ditau1_->setTau1And2P4(tau1P4, tau2P4);
  adapter_ditau2_->setTau1And2P4(tau3P4, tau4P4);
}

void HistogramAdapterDiHiggs::bookHistograms(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met)
{
  quantity_pt_->bookHistogram(vis1P4, vis2P4, vis3P4, vis4P4, met);
  quantity_eta_->bookHistogram(vis1P4, vis2P4, vis3P4, vis4P4, met);
  quantity_phi_->bookHistogram(vis1P4, vis2P4, vis3P4, vis4P4, met);
  quantity_mass_->bookHistogram(vis1P4, vis2P4, vis3P4, vis4P4, met);
  quantity_transverseMass_->bookHistogram(vis1P4, vis2P4, vis3P4, vis4P4, met);
  adapter_ditau1_->bookHistograms(vis1P4, vis2P4, met);
  adapter_ditau2_->bookHistograms(vis3P4, vis4P4, met);
}

void HistogramAdapterDiHiggs::fillHistograms(const LorentzVector& tau1P4, const LorentzVector& tau2P4, const LorentzVector& ditau1P4,
					     const LorentzVector& tau3P4, const LorentzVector& tau4P4, const LorentzVector& ditau2P4,
					     const LorentzVector& dihiggsP4,
					     const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met) const
{
  quantity_pt_->fillHistogram(dihiggsP4.pt());
  quantity_eta_->fillHistogram(dihiggsP4.eta());
  quantity_phi_->fillHistogram(dihiggsP4.phi());
  quantity_mass_->fillHistogram(dihiggsP4.mass());
  double transverseMass2 = square(ditau1P4.Et() + ditau2P4.Et()) - (square(dihiggsP4.px()) + square(dihiggsP4.py()));
  quantity_transverseMass_->fillHistogram(TMath::Sqrt(TMath::Max(1., transverseMass2)));
  adapter_ditau1_->fillHistograms(tau1P4, tau2P4, ditau1P4, vis1P4, vis2P4, met);
  adapter_ditau2_->fillHistograms(tau3P4, tau4P4, ditau2P4, vis3P4, vis4P4, met);
}

HistogramAdapterDiTau* HistogramAdapterDiHiggs::ditau1() const 
{ 
  return adapter_ditau1_; 
}
 
HistogramAdapterDiTau* HistogramAdapterDiHiggs::ditau2() const 
{
  return adapter_ditau2_; 
}

double HistogramAdapterDiHiggs::getPt() const
{
  return extractValue(quantity_pt_);
}

double HistogramAdapterDiHiggs::getPtErr() const
{
  return extractUncertainty(quantity_pt_);
}

double HistogramAdapterDiHiggs::getPtLmax() const
{
  return extractLmax(quantity_pt_);
}

double HistogramAdapterDiHiggs::getEta() const
{
  return extractValue(quantity_eta_);
}

double HistogramAdapterDiHiggs::getEtaErr() const
{
  return extractUncertainty(quantity_eta_);
}

double HistogramAdapterDiHiggs::getEtaLmax() const
{
  return extractLmax(quantity_eta_);
}

double HistogramAdapterDiHiggs::getPhi() const
{
  return extractValue(quantity_phi_);
}

double HistogramAdapterDiHiggs::getPhiErr() const
{
  return extractUncertainty(quantity_phi_);
}

double HistogramAdapterDiHiggs::getPhiLmax() const
{
  return extractLmax(quantity_phi_);
}

double HistogramAdapterDiHiggs::getMass() const
{
  return extractValue(quantity_mass_);
}

double HistogramAdapterDiHiggs::getMassErr() const
{
  return extractUncertainty(quantity_mass_);
}

double HistogramAdapterDiHiggs::getMassLmax() const
{
  return extractLmax(quantity_mass_);
}

double HistogramAdapterDiHiggs::getTransverseMass() const
{
  return extractValue(quantity_transverseMass_);
}

double HistogramAdapterDiHiggs::getTransverseMassErr() const
{
  return extractUncertainty(quantity_transverseMass_);
}

double HistogramAdapterDiHiggs::getTransverseMassLmax() const
{
  return extractLmax(quantity_transverseMass_);
}

classic_svFit::LorentzVector HistogramAdapterDiHiggs::getP4() const
{
  TLorentzVector p4;
  p4.SetPtEtaPhiM(this->getPt(), this->getEta(), this->getPhi(), this->getMass());
  return classic_svFit::LorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E());
}

double HistogramAdapterDiHiggs::DoEval(const double* x) const
{
  fillHistograms(tau1P4_, tau2P4_, ditau1P4_, tau3P4_, tau4P4_, ditau2P4_, dihiggsP4_, vis1P4_, vis2P4_, vis3P4_, vis4P4_, met_);
  return 0.;
}
//-------------------------------------------------------------------------------------------------
