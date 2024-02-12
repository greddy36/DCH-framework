#include "TauAnalysis/ClassicSVfit4tau/interface/ClassicSVfitIntegrand4tau.h"

#include "TauAnalysis/ClassicSVfit/interface/SVfitIntegratorMarkovChain.h"

#include <TMath.h>
#include <TString.h> // Form
#include <Math/VectorUtil.h>

#include <math.h>

using namespace classic_svFit;

/// global function pointer, needed for Markov Chain integration
const ClassicSVfitIntegrand4tau* ClassicSVfitIntegrand4tau::gSVfitIntegrand = 0;

ClassicSVfitIntegrand4tau::ClassicSVfitIntegrand4tau(int verbosity)
  : ClassicSVfitIntegrandBase(verbosity)
  , fittedTauLepton1_(0, verbosity)
  , fittedTauLepton2_(1, verbosity)
  , fittedTauLepton3_(2, verbosity)
  , fittedTauLepton4_(3, verbosity)
  , diTau1MassConstraint_(-1.)
  , diTau2MassConstraint_(-1.)
  , histogramAdapter_(nullptr)
{
  if ( verbosity_ >= 2 ) {
    std::cout << "<ClassicSVfitIntegrand4tau::ClassicSVfitIntegrand4tau>:" << std::endl;
  }

  numTaus_ = 4;
  legIntegrationParams_.resize(numTaus_);

  maxNumberOfDimensions_ = 3*numTaus_;
  xMin_ = new double[maxNumberOfDimensions_];
  xMax_ = new double[maxNumberOfDimensions_];
  x_ = new double[maxNumberOfDimensions_];

  // CV: disable log(M) term unless explicitely requested by user
  addLogM_fixed_ = false;
  addLogM_fixed_power_ = 0.; 

  fittedTauLeptons_.resize(numTaus_);
  fittedTauLeptons_[0] = &fittedTauLepton1_;
  fittedTauLeptons_[1] = &fittedTauLepton2_;
  fittedTauLeptons_[2] = &fittedTauLepton3_;
  fittedTauLeptons_[3] = &fittedTauLepton4_;

  // set global function pointer to this
  gSVfitIntegrand = this;
}

ClassicSVfitIntegrand4tau::~ClassicSVfitIntegrand4tau()
{
  if ( verbosity_ >= 2 ) {
    std::cout << "<ClassicSVfitIntegrand4tau::~ClassicSVfitIntegrand4tau>:" << std::endl;
  }
}

void ClassicSVfitIntegrand4tau::setDiTau1MassConstraint(double diTauMass)
{
  diTau1MassConstraint_ = diTauMass;
  diTau1MassConstraint2_ = square(diTau1MassConstraint_);
}

void ClassicSVfitIntegrand4tau::setDiTau2MassConstraint(double diTauMass)
{
  diTau2MassConstraint_ = diTauMass;
  diTau2MassConstraint2_ = square(diTau2MassConstraint_);
}

void ClassicSVfitIntegrand4tau::setDiHiggsMassConstraint(double diHiggsMass)
{
  diHiggsMassConstraint_ = diHiggsMass;
  diHiggsMassConstraint2_ = square(diHiggsMassConstraint_);
}

void ClassicSVfitIntegrand4tau::setHistogramAdapter(HistogramAdapterDiHiggs* histogramAdapter)
{
  histogramAdapter_ = histogramAdapter;
}

void ClassicSVfitIntegrand4tau::setLeptonInputs(const std::vector<MeasuredTauLepton>& measuredTauLeptons)
{
  if ( verbosity_ >= 2 ) {
    std::cout << "<ClassicSVfitIntegrand4tau::setLeptonInputs>:" << std::endl;
  }

  ClassicSVfitIntegrandBase::setLeptonInputs(measuredTauLeptons);

  // set momenta of visible tau decay products, reset momenta of reconstructed tau leptons
  measuredTauLepton1_ = measuredTauLeptons[0];
  fittedTauLepton1_.setMeasuredTauLepton(measuredTauLepton1_);
  leg1isLeptonicTauDecay_ = measuredTauLepton1_.isLeptonicTauDecay();
  leg1isHadronicTauDecay_ = measuredTauLepton1_.isHadronicTauDecay();
  leg1isPrompt_ = measuredTauLepton1_.isPrompt();
  measuredTauLepton2_ = measuredTauLeptons[1];
  fittedTauLepton2_.setMeasuredTauLepton(measuredTauLepton2_);
  leg2isLeptonicTauDecay_ = measuredTauLepton2_.isLeptonicTauDecay();
  leg2isHadronicTauDecay_ = measuredTauLepton2_.isHadronicTauDecay();
  leg2isPrompt_ = measuredTauLepton2_.isPrompt();
  measuredTauLepton3_ = measuredTauLeptons[2];
  fittedTauLepton3_.setMeasuredTauLepton(measuredTauLepton3_);
  leg3isLeptonicTauDecay_ = measuredTauLepton3_.isLeptonicTauDecay();
  leg3isHadronicTauDecay_ = measuredTauLepton3_.isHadronicTauDecay();
  leg3isPrompt_ = measuredTauLepton3_.isPrompt();
  measuredTauLepton4_ = measuredTauLeptons[3];
  fittedTauLepton4_.setMeasuredTauLepton(measuredTauLepton4_);
  leg4isLeptonicTauDecay_ = measuredTauLepton4_.isLeptonicTauDecay();
  leg4isHadronicTauDecay_ = measuredTauLepton4_.isHadronicTauDecay();
  leg4isPrompt_ = measuredTauLepton4_.isPrompt();

  ditau1_mVis_measured_ = (measuredTauLepton1_.p4() + measuredTauLepton2_.p4()).mass();
  if ( verbosity_ >= 1 ) {
    std::cout << "mVis(1st ditau) = " << ditau1_mVis_measured_ << std::endl;
  }
  ditau1_mVis2_measured_ = square(ditau1_mVis_measured_);
  ditau2_mVis_measured_ = (measuredTauLepton3_.p4() + measuredTauLepton4_.p4()).mass();
  if ( verbosity_ >= 1 ) {
    std::cout << "mVis(2nd ditau) = " << ditau2_mVis_measured_ << std::endl;
  }
  ditau2_mVis2_measured_ = square(ditau2_mVis_measured_);
  if ( verbosity_ >= 2 ) {
    std::cout << "mVis(4tau) = " << (measuredTauLepton1_.p4() + measuredTauLepton2_.p4() + measuredTauLepton3_.p4() + measuredTauLepton4_.p4()).mass() << std::endl;
  }

  //-----------------------------------------------------------------------------
  // CV: only used in conjunction with VAMP integration
  mVis12_ = ditau1_mVis_measured_;
  //std::cout << "mVis12 = " << mVis12_ << std::endl;
  mVis13_ = (measuredTauLepton1_.p4() + measuredTauLepton3_.p4()).mass();
  //std::cout << "mVis13 = " << mVis13_ << std::endl;
  mVis14_ = (measuredTauLepton1_.p4() + measuredTauLepton4_.p4()).mass();
  //std::cout << "mVis14 = " << mVis14_ << std::endl;
  mVis23_ = (measuredTauLepton2_.p4() + measuredTauLepton3_.p4()).mass();
  //std::cout << "mVis23 = " << mVis23_ << std::endl;
  mVis24_ = (measuredTauLepton2_.p4() + measuredTauLepton4_.p4()).mass();
  //std::cout << "mVis24 = " << mVis24_ << std::endl;
  mVis34_ = ditau2_mVis_measured_;
  //std::cout << "mVis34 = " << mVis34_ << std::endl;
  //-----------------------------------------------------------------------------
}

double ClassicSVfitIntegrand4tau::EvalPS(const double* q) const
{
  rescaleX(q);

  if ( verbosity_ >= 2 ) {
    std::cout << "<ClassicSVfitIntegrand4tau::EvalPS(const double*)>:" << std::endl;
    std::cout << " x = { ";
    for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
      std::cout << x_[iDimension];
      if ( iDimension < (numDimensions_ - 1) ) std::cout << ", ";
    }
    std::cout << " }" << std::endl;
  }

  // in case of initialization errors don't start to do anything
  if ( errorCode_ & MatrixInversion ||
       errorCode_ & LeptonNumber    ||
       errorCode_ & TestMass        ) {
    return 0.; 
  }

  double visPtShift1 = 1.;
  double visPtShift2 = 1.;
  double visPtShift3 = 1.;
  double visPtShift4 = 1.;
#ifdef USE_SVFITTF
  int idx_visPtShift1 = legIntegrationParams_[0].idx_VisPtShift_;
  if( useHadTauTF_ && idx_visPtShift1 != -1 && !leg1isLeptonicTauDecay_ ) visPtShift1 = (1./x_[idx_visPtShift1]);
  int idx_visPtShift2 = legIntegrationParams_[1].idx_VisPtShift_;
  if( useHadTauTF_ && idx_visPtShift2 != -1 && !leg2isLeptonicTauDecay_ ) visPtShift2 = (1./x_[idx_visPtShift2]);
  int idx_visPtShift3 = legIntegrationParams_[2].idx_VisPtShift_;
  if( useHadTauTF_ && idx_visPtShift3 != -1 && !leg3isLeptonicTauDecay_ ) visPtShift3 = (1./x_[idx_visPtShift3]);
  int idx_visPtShift4 = legIntegrationParams_[3].idx_VisPtShift_;
  if( useHadTauTF_ && idx_visPtShift4 != -1 && !leg4isLeptonicTauDecay_ ) visPtShift4 = (1./x_[idx_visPtShift4]);
#endif 
  if ( visPtShift1 < 1.e-2 || visPtShift2 < 1.e-2 || visPtShift3 < 1.e-2 || visPtShift4 < 1.e-2 ) return 0.;

  // scale momenta of visible tau decays products
  fittedTauLepton1_.updateVisMomentum(visPtShift1);
  fittedTauLepton2_.updateVisMomentum(visPtShift2);
  fittedTauLepton3_.updateVisMomentum(visPtShift3);
  fittedTauLepton4_.updateVisMomentum(visPtShift4);

  // compute visible energy fractions for both taus of first tau pair
  double x1_dash = 1.;  
  bool x1isConstrained = true;
  if ( !leg1isPrompt_ ) {
    int idx_x1 = legIntegrationParams_[0].idx_X_;
    assert(idx_x1 != -1);
    x1_dash = x_[idx_x1];
    x1isConstrained = false;
  }
  double x1 = x1_dash/visPtShift1;
  if ( !(x1 >= 1.e-5 && x1 <= 1.) ) return 0.;

  double x2_dash = 1.;
  bool x2isConstrained = true;
  if ( !leg2isPrompt_ ) {
    int idx_x2 = legIntegrationParams_[1].idx_X_;
    if ( idx_x2 != -1 ) {
      x2_dash = x_[idx_x2];
      x2isConstrained = false;
    } else {
      x2_dash = (ditau1_mVis2_measured_/diTau1MassConstraint2_)/x1_dash;
    }
  }
  double x2 = x2_dash/visPtShift2;
  if ( !(x2 >= 1.e-5 && x2 <= 1.) ) return 0.;

  // compute visible energy fractions for both taus of second tau pair
  std::vector<double> x3_dash_solutions;
  bool x3isConstrained = true;
  double term_b = 0.;
  double term_c = 0.;
  double term_d = 0.;
  double term_e = 0.;
  if ( !leg3isPrompt_ ) {
    int idx_x3 = legIntegrationParams_[2].idx_X_;
    if ( idx_x3 != -1 ) {
      double x3_dash = x_[idx_x3];
      x3_dash_solutions.push_back(x3_dash);
      x3isConstrained = false;
    } else {
      //-------------------------------------------------------------------------
      // CV: only used in conjunction with VAMP integration
      //std::cout << "diHiggsMassConstraint2 = " << diHiggsMassConstraint2_ << std::endl;
      assert(diTau1MassConstraint2_ > 0. && diTau2MassConstraint2_ > 0.);
      double term_a1 = diTau1MassConstraint2_;
      double term_a2 = diTau2MassConstraint2_;
      term_b = square(mVis13_);
      //std::cout << "term_b = " << term_b << std::endl;
      assert(mVis34_ > 0.);
      term_c = square(mVis14_/mVis34_)*diTau2MassConstraint2_;
      //std::cout << "term_c = " << term_c << std::endl;
      assert(mVis12_ > 0.);
      term_d = square(mVis23_/mVis12_)*diTau1MassConstraint2_;
      //std::cout << "term_d = " << term_d << std::endl;
      term_e = square(mVis24_/(mVis12_*mVis34_))*diTau1MassConstraint2_*diTau2MassConstraint2_;
      //std::cout << "term_e = " << term_e << std::endl;
      double term1 = (diHiggsMassConstraint2_ - (term_a1 + term_a2))*x1_dash;
      //std::cout << "term1 = " << term1 << std::endl;
      double x1_dash2 = square(x1_dash);
      //std::cout << "x1_dash2 = " << x1_dash2 << std::endl;
      double term2_2 = square(term1) - 4.*(term_b + term_d*x1_dash2)*(term_c + term_e*x1_dash2);
      //std::cout << "term2^2 = " << term2_2 << std::endl;
      //assert(term2_2 >= 0.);
      //if ( term2_2 > 0. ) std::cout << "term2^2 = " << term2_2 << std::endl;
      if ( !(term2_2 >= 0.) ) return 0.;
      double term2 = TMath::Sqrt(term2_2);
      double term3 = 2.*(term_c + term_e*x1_dash2);      
      double x3_dash_plus = (term1 + term2)/term3;
      double x3_dash_minus = (term1 - term2)/term3;
      //std::cout << "x3: '+' solution = " << x3_dash_plus << ", '-' solution = " << x3_dash_minus << std::endl;
      x3_dash_solutions.push_back(x3_dash_plus);
      x3_dash_solutions.push_back(x3_dash_minus);
      //-------------------------------------------------------------------------
    }
  }

  double prob = 0.;
  for ( std::vector<double>::const_iterator x3_dash = x3_dash_solutions.begin();
	x3_dash != x3_dash_solutions.end(); ++x3_dash ) {
    double x3 = (*x3_dash)/visPtShift3;
    if ( !(x3 >= 1.e-5 && x3 <= 1.) ) continue;
    
    double x4_dash = 1.;
    bool x4isConstrained = true;
    if ( !leg4isPrompt_ ) {
      int idx_x4 = legIntegrationParams_[3].idx_X_;
      if ( idx_x4 != -1 ) {
	x4_dash = x_[idx_x4];
	x4isConstrained = false;
      } else {
	if ( diTau2MassConstraint_ > 0. ) {
	  x4_dash = (ditau2_mVis2_measured_/diTau2MassConstraint2_)/(*x3_dash);
	} else {
	  double ditau1_mTauTau2_approx = ditau1_mVis2_measured_/(x1_dash*x2_dash);
	  x4_dash = (ditau2_mVis2_measured_/ditau1_mTauTau2_approx)/(*x3_dash);
	}
      }
    }
    double x4 = x4_dash/visPtShift4;
    if ( !(x4 >= 1.e-5 && x4 <= 1.) ) continue;

    prob += ClassicSVfitIntegrand4tau::EvalPS(
      x1,  x1_dash, visPtShift1, x1isConstrained,
      x2,  x2_dash, visPtShift2, x2isConstrained,
      x3, *x3_dash, visPtShift3, x3isConstrained,
      x4,  x4_dash, visPtShift4, x4isConstrained,
      term_b, term_c, term_d, term_e);
  }

  return prob;
}

double ClassicSVfitIntegrand4tau::EvalPS(double x1, double x1_dash, double visPtShift1, bool x1isConstrained,
					 double x2, double x2_dash, double visPtShift2, bool x2isConstrained, 
					 double x3, double x3_dash, double visPtShift3, bool x3isConstrained, 
					 double x4, double x4_dash, double visPtShift4, bool x4isConstrained, 
					 double term_b, double term_c, double term_d, double term_e) const
{
  // compute neutrino and tau lepton momenta 
  if ( !leg1isPrompt_ ) {
    int idx_phiNu1 = legIntegrationParams_[0].idx_phi_;
    assert(idx_phiNu1 != -1);
    double phiNu1 = x_[idx_phiNu1];
    int idx_nu1Mass = legIntegrationParams_[0].idx_mNuNu_;
    double nu1Mass = ( idx_nu1Mass != -1 ) ? TMath::Sqrt(x_[idx_nu1Mass]) : 0.;
    fittedTauLepton1_.updateTauMomentum(x1, phiNu1, nu1Mass);
    //std::cout << "fittedTauLepton1: errorCode = " << fittedTauLepton1_.errorCode() << std::endl;
    if ( fittedTauLepton1_.errorCode() != FittedTauLepton::None ) {
      errorCode_ |= TauDecayParameters;
      return 0.;
    }
  }

  if ( !leg2isPrompt_ ) {
    int idx_phiNu2 = legIntegrationParams_[1].idx_phi_;
    assert(idx_phiNu2 != -1);
    double phiNu2 = x_[idx_phiNu2];
    int idx_nu2Mass = legIntegrationParams_[1].idx_mNuNu_;
    double nu2Mass = ( idx_nu2Mass != -1 ) ? TMath::Sqrt(x_[idx_nu2Mass]) : 0.;
    fittedTauLepton2_.updateTauMomentum(x2, phiNu2, nu2Mass);
    //std::cout << "fittedTauLepton2: errorCode = " << fittedTauLepton2_.errorCode() << std::endl;
    if ( fittedTauLepton2_.errorCode() != FittedTauLepton::None ) {
      errorCode_ |= TauDecayParameters;
      return 0.;
    }
  }

  if ( !leg3isPrompt_ ) {
    int idx_phiNu3 = legIntegrationParams_[2].idx_phi_;
    assert(idx_phiNu3 != -1);
    double phiNu3 = x_[idx_phiNu3];
    int idx_nu3Mass = legIntegrationParams_[2].idx_mNuNu_;
    double nu3Mass = ( idx_nu3Mass != -1 ) ? TMath::Sqrt(x_[idx_nu3Mass]) : 0.;
    fittedTauLepton3_.updateTauMomentum(x3, phiNu3, nu3Mass);
    //std::cout << "fittedTauLepton3: errorCode = " << fittedTauLepton3_.errorCode() << std::endl;
    if ( fittedTauLepton3_.errorCode() != FittedTauLepton::None ) {
      errorCode_ |= TauDecayParameters;
      return 0.;
    }
  }

  if ( !leg4isPrompt_ ) {
    int idx_phiNu4 = legIntegrationParams_[3].idx_phi_;
    assert(idx_phiNu4 != -1);
    double phiNu4 = x_[idx_phiNu4];
    int idx_nu4Mass = legIntegrationParams_[3].idx_mNuNu_;
    double nu4Mass = ( idx_nu4Mass != -1 ) ? TMath::Sqrt(x_[idx_nu4Mass]) : 0.;
    fittedTauLepton4_.updateTauMomentum(x4, phiNu4, nu4Mass);
    //std::cout << "fittedTauLepton4: errorCode = " << fittedTauLepton4_.errorCode() << std::endl;
    if ( fittedTauLepton4_.errorCode() != FittedTauLepton::None ) {
      errorCode_ |= TauDecayParameters;
      return 0.;
    }
  }

  if ( verbosity_ >= 2 ) {
    for ( unsigned iTau = 0; iTau < numTaus_; ++iTau ) {
      const FittedTauLepton* fittedTauLepton = fittedTauLeptons_[iTau];
      const LorentzVector& visP4 = fittedTauLepton->visP4();
      const LorentzVector& nuP4 = fittedTauLepton->nuP4();
      const LorentzVector& tauP4 = fittedTauLepton->tauP4();
      std::cout << "leg" << (iTau + 1) << ": En = " << visP4.E() << ", Px = " << visP4.px()
		<< ", Py = " << visP4.py() << ", Pz = " << visP4.pz() << ";"
		<< " Pt = " << visP4.pt() << ", eta = " << visP4.eta()
		<< ", phi = " << visP4.phi() << ", mass = " << visP4.mass()
		<< " (x = " << fittedTauLepton->x() << ")" << std::endl;
      std::cout << "tau" << (iTau + 1) << ": En = " << tauP4.E() << ", Px = " << tauP4.px() << ", Py = " << tauP4.py() << ", Pz = " << tauP4.pz() << ";"
		<< " Pt = " << tauP4.pt() << ", eta = " << tauP4.eta() << ", phi = " << tauP4.phi() << std::endl;
      std::cout << "nu" << (iTau + 1) << ": En = " << nuP4.E() << ", Px = " << nuP4.px() << ", Py = " << nuP4.py() << ", Pz = " << nuP4.pz() << ";"
		<< " Pt = " << nuP4.pt() << ", eta = " << nuP4.eta() << ", phi = " << nuP4.phi() << ", mass = " << nuP4.mass() << std::endl;
    }
  }

  double prob_PS_and_tauDecay = classic_svFit::constFactor;
  double prob_tauDecay = 1.;
  double prob_TF = 1.;
  for ( unsigned iTau = 0; iTau < numTaus_; ++iTau ) {
    const FittedTauLepton* fittedTauLepton = fittedTauLeptons_[iTau];
    const MeasuredTauLepton& measuredTauLepton = fittedTauLepton->getMeasuredTauLepton();
    double x = fittedTauLepton->x();
    double nuMass = fittedTauLepton->nuMass();
    const LorentzVector& visP4 = fittedTauLepton->visP4();
    const LorentzVector& nuP4 = fittedTauLepton->nuP4();

    // evaluate tau decay matrix elements
    double prob = 1.;
    if      ( measuredTauLepton.isLeptonicTauDecay() ) prob = compPSfactor_tauToLepDecay(x, visP4.E(), visP4.P(), measuredTauLepton.mass(), nuP4.E(), nuP4.P(), nuMass);
    else if ( measuredTauLepton.isHadronicTauDecay() ) prob = compPSfactor_tauToHadDecay(x, visP4.E(), visP4.P(), measuredTauLepton.mass(), nuP4.E(), nuP4.P());
    prob_tauDecay *= prob;

    // evaluate transfer functions for tau energy reconstruction
#ifdef USE_SVFITTF
    if ( useHadTauTF_ && legIntegrationParams_[iTau].idx_VisPtShift_ != -1 && measuredTauLepton.isHadronicTauDecay() ) {
      double prob = (*hadTauTFs_[iTau])(measuredTauLepton.pt(), visP4.pt(), visP4.eta());
      if ( verbosity_ >= 2 ) {
	std::cout << "TF(leg" << iTau << "): recPt = " << measuredTauLepton.pt() << ", genPt = " << visP4.pt()
		  << ", genEta = " << visP4.eta() << " --> prob = " << prob << std::endl;
      }
      prob_TF *= prob;
    }
#endif
  }
  prob_PS_and_tauDecay *= prob_tauDecay;
  prob_PS_and_tauDecay *= classic_svFit::matrixElementNorm;

  double mhh = (fittedTauLepton1_.tauP4() + fittedTauLepton2_.tauP4() + fittedTauLepton3_.tauP4() + fittedTauLepton4_.tauP4()).mass();
  double prob_logM = 1.;
  if ( addLogM_fixed_ ) {
    prob_logM = 1./TMath::Power(TMath::Max(1., mhh), addLogM_fixed_power_);
  }
  if ( addLogM_dynamic_ ) {    
    double addLogM_power = addLogM_dynamic_formula_->Eval(mhh);
    prob_logM = 1./TMath::Power(TMath::Max(1., mhh), TMath::Max(0., addLogM_power));
  }

  double jacobiFactor = 1.;
  jacobiFactor *= 1./(visPtShift1*visPtShift2); // product of derrivatives dx1/dx1' and dx2/dx2' for parametrization of x1, x2 by x1', x2'  
  if ( x2isConstrained ) {
    if ( diTau1MassConstraint_ > 0. ) {
      jacobiFactor *= (x2/square(diTau1MassConstraint_));
    } else {
      double diTau1Mass = (fittedTauLepton1_.tauP4() + fittedTauLepton2_.tauP4()).mass();
      jacobiFactor *= (x2/square(diTau1Mass));
    }
  }
  jacobiFactor *= 1./(visPtShift3*visPtShift4); // product of derrivatives dx3/dx3' and dx4/dx4' for parametrization of x3, x4 by x3', x4'  
  //-----------------------------------------------------------------------------
  // CV: only used in conjunction with VAMP integration
  if ( x3isConstrained ) {
    double x3_dash2 = square(x3_dash);
    jacobiFactor *= 1./TMath::Abs(term_c/x1_dash + term_e*x1_dash - term_b/(x1_dash*x3_dash2) - term_d*x1_dash/x3_dash2);
  }
  //-----------------------------------------------------------------------------
  if ( x4isConstrained ) {
    if ( diTau2MassConstraint_ > 0. ) {
      jacobiFactor *= (x4/square(diTau2MassConstraint_));
    } else {
      double diTau2Mass = (fittedTauLepton3_.tauP4() + fittedTauLepton4_.tauP4()).mass();
      jacobiFactor *= (x4/square(diTau2Mass));
    }
  }

  //static int numCalls = 0;
  //++numCalls;
  //std::cout << "call #" << numCalls << ":" << std::endl;
  //if ( numCalls > 1000 ) assert(0);

  double prob = prob_PS_and_tauDecay*prob_TF*prob_logM*jacobiFactor;
  if ( verbosity_ >= 2 ) {
    std::cout << "mTauTau(1) = " << (fittedTauLepton1_.tauP4() + fittedTauLepton2_.tauP4()).mass() << ","
	      << " mTauTau(2) = " << (fittedTauLepton3_.tauP4() + fittedTauLepton4_.tauP4()).mass() << ","
	      << " mhh = " << mhh << std::endl;
    std::cout << "prob: PS+decay = " << prob_PS_and_tauDecay << ","
              << " TF = " << prob_TF << ", log(M) = " << prob_logM << ", Jacobi = " << jacobiFactor 
	      << " --> returning " << prob << std::endl;
  }
  if ( TMath::IsNaN(prob) ) {
    prob = 0.;
  }

  return prob;
}

double ClassicSVfitIntegrand4tau::Eval(const double* x, unsigned int iComponent) const
{
  if ( iComponent == 0 ) {
    phaseSpaceComponentCache_ = EvalPS(x);
  }
  if ( phaseSpaceComponentCache_ < 1.e-300 ) return 0.;
  double prob_metTF = EvalMET_TF(iComponent);
  double prob = phaseSpaceComponentCache_*prob_metTF;
  if ( verbosity_ >= 2 ) {
    std::cout << " metTF: " << prob_metTF << ","
	      << " phaseSpaceComponentCache: " << phaseSpaceComponentCache_
	      << " --> returning " << prob << std::endl;
  }
  if ( histogramAdapter_ && prob > 1.e-300 ){
    histogramAdapter_->setTau123And4P4(fittedTauLepton1_.tauP4(), fittedTauLepton2_.tauP4(), fittedTauLepton3_.tauP4(), fittedTauLepton4_.tauP4());
  }
  return prob;
}
