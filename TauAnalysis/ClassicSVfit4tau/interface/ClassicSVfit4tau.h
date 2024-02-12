#ifndef TauAnalysis_ClassicSVfit4tau_ClassicSVfit4tau_h
#define TauAnalysis_ClassicSVfit4tau_ClassicSVfit4tau_h

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitBase.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit4tau/interface/svFitHistogramAdapter4tau.h"

class ClassicSVfit4tau : public ClassicSVfitBase
{
 public:
  enum { kAlgoMarkovChain, kAlgoVAMP };
  ClassicSVfit4tau(int = kAlgoMarkovChain, int = 0);
  ~ClassicSVfit4tau();

  void setDiTau1MassConstraint(double diTauMass);
  void setDiTau2MassConstraint(double diTauMass);

  /// set and get histogram adapter
  void setHistogramAdapter(classic_svFit::HistogramAdapterDiHiggs* histogramAdapter);
  classic_svFit::HistogramAdapterDiHiggs* getHistogramAdapter() const;

  /// prepare the integrand
  void prepareIntegrand();

  /// prepare input measurements
  void prepareLeptonInput(const std::vector<classic_svFit::MeasuredTauLepton>& measuredTauLeptons);

  /// run integration with Markov Chain
  void integrate(const std::vector<classic_svFit::MeasuredTauLepton>&, double, double, const TMatrixD&);

  /// access reconstructed mass value and uncertainty
  double getMass() const;
  double getMassErr() const;

  /// access likelihood for reconstructed mass value
  double getLmax() const;

 protected:
  /// set integration algorithm
  /// (either VAMP or Markov Chain integration)
  int intAlgoFlag_; 

  void integrateMC();
  void integrateVAMP();

  /// initialize Markov Chain integrator class
  void initializeMCIntegrator();

  /// set integration indices and ranges for both legs
  /// when useMassConstraint is true reduce number of
  /// dimension by using the mass contraint
  void setIntegrationParams(bool useDiTau1MassConstraint=false, bool useDiTau2MassConstraint=true);

  double diTau1MassConstraint_;
  double diTau2MassConstraint_;
  
  /// histograms for evaluation of pT, eta, phi, mass and transverse mass of di-Higgs system (built from 4 taus)
  mutable classic_svFit::HistogramAdapterDiHiggs* histogramAdapter_;

  /// reconstructed mass value and uncertainty
  double dihiggsMass_;
  double dihiggsMassErr_;

  /// likelihood for reconstructed mass value
  double Lmax_;
};

#endif
