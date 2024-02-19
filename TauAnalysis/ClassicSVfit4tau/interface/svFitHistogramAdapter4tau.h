#ifndef TauAnalysis_ClassicSVfit4tau_svFitHistogramAdapter4tau_h
#define TauAnalysis_ClassicSVfit4tau_svFitHistogramAdapter4tau_h

#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h" // HistogramTools, HistogramAdapterDiTau

#include <Math/Functor.h>
#include <TH1.h>

namespace classic_svFit
{
  //-------------------------------------------------------------------------------------------------
  // auxiliary classes to reconstruct pT, eta, phi, mass, and transverse mass of Higgs boson pairs (built from 4 taus)
  class SVfitQuantityDiHiggs : public SVfitQuantity
  {
   public:
    SVfitQuantityDiHiggs(const std::string& label);

    virtual TH1* createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met) const = 0;

    void bookHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met);
  };

  class SVfitQuantityDiHiggsPt : public SVfitQuantityDiHiggs
  {
   public:
    SVfitQuantityDiHiggsPt(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met) const;
  };

  class SVfitQuantityDiHiggsEta : public SVfitQuantityDiHiggs
  {
   public:
    SVfitQuantityDiHiggsEta(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met) const;
  };

  class SVfitQuantityDiHiggsPhi : public SVfitQuantityDiHiggs
  {
   public:
    SVfitQuantityDiHiggsPhi(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met) const;
  };

  class SVfitQuantityDiHiggsMass : public SVfitQuantityDiHiggs
  {
   public:
    SVfitQuantityDiHiggsMass(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met) const;
  };

  class SVfitQuantityDiHiggsTransverseMass : public SVfitQuantityDiHiggs
  {
   public:
    SVfitQuantityDiHiggsTransverseMass(const std::string& label);
    virtual TH1* createHistogram(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met) const;
  };

  class HistogramAdapterDiHiggs : public HistogramAdapter
  {
   public:
    HistogramAdapterDiHiggs(const std::string& label = "dihiggs");
    ~HistogramAdapterDiHiggs();

    void bookHistograms(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met);

    void setMeasurement(const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met);
    void setTau123And4P4(const LorentzVector& tau1P4,  const LorentzVector& tau2P4, const LorentzVector& tau3P4, const LorentzVector& tau4P4);

    void fillHistograms(const LorentzVector& tau1P4, const LorentzVector& tau2P4, const LorentzVector& ditau1P4,
			const LorentzVector& tau3P4, const LorentzVector& tau4P4, const LorentzVector& ditau2P4,
			const LorentzVector& dihiggsP4,
			const LorentzVector& vis1P4, const LorentzVector& vis2P4, const LorentzVector& vis3P4, const LorentzVector& vis4P4, const Vector& met) const;

    HistogramAdapterDiTau* ditau1() const;
    HistogramAdapterDiTau* ditau2() const;

    /// get pT, eta, phi, mass and transverse mass of di-Higgs system (built from 4 taus)
    double getPt() const;
    double getPtErr() const;
    double getPtLmax() const;
    double getEta() const;
    double getEtaErr() const;
    double getEtaLmax() const;
    double getPhi() const;
    double getPhiErr() const;
    double getPhiLmax() const;
    double getMass() const;
    double getMassErr() const;
    double getMassLmax() const;
    double getTransverseMass() const;
    double getTransverseMassErr() const;
    double getTransverseMassLmax() const;

    /// convenient access to four-vector of di-Higgs system 
    LorentzVector getP4() const;

  private:
    double DoEval(const double* x) const;

   protected:
    LorentzVector vis1P4_;
    LorentzVector vis2P4_;
    LorentzVector vis3P4_;
    LorentzVector vis4P4_;
    Vector met_;
    LorentzVector tau1P4_;
    LorentzVector tau2P4_;    
    LorentzVector ditau1P4_;
    LorentzVector tau3P4_;
    LorentzVector tau4P4_;    
    LorentzVector ditau2P4_;
    LorentzVector dihiggsP4_;

    SVfitQuantityDiHiggsPt* quantity_pt_;
    SVfitQuantityDiHiggsEta* quantity_eta_;
    SVfitQuantityDiHiggsPhi* quantity_phi_;
    SVfitQuantityDiHiggsMass* quantity_mass_;
    SVfitQuantityDiHiggsTransverseMass* quantity_transverseMass_;

    HistogramAdapterDiTau* adapter_ditau1_;
    HistogramAdapterDiTau* adapter_ditau2_;
  };
  //-------------------------------------------------------------------------------------------------
}

#endif
