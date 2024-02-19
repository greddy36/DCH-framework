#include "TauAnalysis/ClassicSVfit4tau/interface/ClassicSVfit4tau.h"

#include "TauAnalysis/ClassicSVfit4tau/interface/ClassicSVfitIntegrand4tau.h"
#include "TauAnalysis/ClassicSVfit/interface/SVfitIntegratorMarkovChain.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h" // HistogramTools

#include <TGraphErrors.h>
#include <TH1.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>

#include <algorithm>
#include <assert.h>

using namespace classic_svFit;

namespace
{
  double g_C(const double* x, size_t dim, void* param)
  {
    return ClassicSVfitIntegrand4tau::gSVfitIntegrand->Eval(x);
  }

  double g_Fortran(const double** x, size_t dim, void** param)
  {
    return ClassicSVfitIntegrand4tau::gSVfitIntegrand->Eval(*x);
  }
}

typedef double (*gPtr_Fortran)(const double**, size_t, void**);

extern "C" 
{
  void vamp_integrate_(gPtr_Fortran, const double*, const double*, int*, int*, int*, double*, double*);
}

ClassicSVfit4tau::ClassicSVfit4tau(int intAlgo, int verbosity)
  : ClassicSVfitBase(verbosity)
  , intAlgoFlag_(intAlgo)
  , diTau1MassConstraint_(-1.)
  , diTau2MassConstraint_(-1.)
  , histogramAdapter_(new HistogramAdapterDiHiggs("dihiggs"))
{
  if ( !(intAlgoFlag_ == kAlgoMarkovChain || intAlgoFlag_ == kAlgoVAMP) ) {
    std::cerr << "Error: Invalid integration algorithm " << intAlgoFlag_ << " declared !!" << std::endl;
    assert(0);
  }
  integrand_ = new ClassicSVfitIntegrand4tau(verbosity_);
  legIntegrationParams_.resize(4);
  xl_ = new double[12];
  xh_ = new double[12];
}

ClassicSVfit4tau::~ClassicSVfit4tau()
{
  delete histogramAdapter_;
}

void ClassicSVfit4tau::setDiTau1MassConstraint(double diTauMass)
{
  diTau1MassConstraint_ = diTauMass;
  (static_cast<ClassicSVfitIntegrand4tau*>(integrand_))->setDiTau1MassConstraint(diTau1MassConstraint_);
}

void ClassicSVfit4tau::setDiTau2MassConstraint(double diTauMass)
{
  diTau2MassConstraint_ = diTauMass;
  (static_cast<ClassicSVfitIntegrand4tau*>(integrand_))->setDiTau2MassConstraint(diTau2MassConstraint_);
}

void ClassicSVfit4tau::initializeMCIntegrator()
{
  ClassicSVfitBase::initializeMCIntegrator();
  intAlgo_->registerCallBackFunction(*histogramAdapter_);
}

void ClassicSVfit4tau::setIntegrationParams(bool useDiTau1MassConstraint, bool useDiTau2MassConstraint)
{
  numDimensions_ = 0;
  legIntegrationParams_[0].reset();
  legIntegrationParams_[1].reset();
  legIntegrationParams_[2].reset();
  legIntegrationParams_[3].reset();
  setLegIntegrationParams(0, false);
  setLegIntegrationParams(1, useDiTau1MassConstraint);
  if ( intAlgoFlag_ == kAlgoMarkovChain ) {
    setLegIntegrationParams(2, false);
  } else if ( intAlgoFlag_ == kAlgoVAMP ) {
    setLegIntegrationParams(2, true);
  } else assert(0); 
  setLegIntegrationParams(3, useDiTau2MassConstraint);
  if ( verbosity_ >= 1 ) printIntegrationRange();
}

void ClassicSVfit4tau::prepareIntegrand()
{
  integrand_->setLeptonInputs(measuredTauLeptons_);
  (static_cast<ClassicSVfitIntegrand4tau*>(integrand_))->setHistogramAdapter(histogramAdapter_);
#ifdef USE_SVFITTF
  if ( useHadTauTF_ ) integrand_->enableHadTauTF();
  else integrand_->disableHadTauTF();
#endif
  for ( unsigned iLeg = 0; iLeg < legIntegrationParams_.size(); ++iLeg ) {
    integrand_->setLegIntegrationParams(iLeg, legIntegrationParams_[iLeg]);
  }
  integrand_->setNumDimensions(numDimensions_);
  integrand_->setIntegrationRanges(xl_, xh_);
  ClassicSVfitIntegrand4tau::gSVfitIntegrand = static_cast<ClassicSVfitIntegrand4tau*>(integrand_);
}

void ClassicSVfit4tau::prepareLeptonInput(const std::vector<MeasuredTauLepton>& measuredTauLeptons)
{
  // CV: sort MeasuredTauLeptons such that:
  //      - association of first two MeasuredTauLeptons to first tau pair and
  //      - association of last two MeasuredTauLeptons to second tau pair
  //     remains intact
  assert(measuredTauLeptons.size() == 4);
  std::vector<MeasuredTauLepton> ditau1_measuredTauLeptons;
  ditau1_measuredTauLeptons.push_back(measuredTauLeptons[0]);
  ditau1_measuredTauLeptons.push_back(measuredTauLeptons[1]);
  for ( std::vector<MeasuredTauLepton>::iterator measuredTauLepton = ditau1_measuredTauLeptons.begin();
        measuredTauLepton != ditau1_measuredTauLeptons.end(); ++measuredTauLepton ) {
    measuredTauLepton->roundToNdigits();
  }
  std::sort(ditau1_measuredTauLeptons.begin(), ditau1_measuredTauLeptons.end(), sortMeasuredTauLeptons());
  std::vector<MeasuredTauLepton> ditau2_measuredTauLeptons;
  ditau2_measuredTauLeptons.push_back(measuredTauLeptons[2]);
  ditau2_measuredTauLeptons.push_back(measuredTauLeptons[3]);
  for ( std::vector<MeasuredTauLepton>::iterator measuredTauLepton = ditau2_measuredTauLeptons.begin();
        measuredTauLepton != ditau2_measuredTauLeptons.end(); ++measuredTauLepton ) {
    measuredTauLepton->roundToNdigits();
  }
  std::sort(ditau2_measuredTauLeptons.begin(), ditau2_measuredTauLeptons.end(), sortMeasuredTauLeptons());
  measuredTauLeptons_.clear();
  measuredTauLeptons_.insert(measuredTauLeptons_.end(), ditau1_measuredTauLeptons.begin(), ditau1_measuredTauLeptons.end());
  measuredTauLeptons_.insert(measuredTauLeptons_.end(), ditau2_measuredTauLeptons.begin(), ditau2_measuredTauLeptons.end());
  if ( verbosity_ >= 1 ) printLeptons();
}

void ClassicSVfit4tau::integrate(const std::vector<MeasuredTauLepton>& measuredTauLeptons,
				 double measuredMETx, double measuredMETy,
				 const TMatrixD& covMET)
{
  if ( verbosity_ >= 1 ) std::cout << "<ClassicSVfit4tau::integrate>:" << std::endl;

  clock_->Reset();
  clock_->Start("<ClassicSVfit4tau::integrate>");

  prepareLeptonInput(measuredTauLeptons);
  integrand_->clearMET();
  addMETEstimate(measuredMETx, measuredMETy, covMET);
  bool useDiTau1MassConstraint = (diTau1MassConstraint_ > 0);
  // CV: require mass of second tau pair to match mass of first tau pair, 
  //     in case mass of first and second tau pair is not explicitely constrained by user
  //    (NOTE: logic needs to match code for setting x4_dash variable in ClassicSVfitIntegrand4tau::EvalPS function)
  bool useDiTau2MassConstraint = (diTau2MassConstraint_ > 0 || !useDiTau1MassConstraint);
  setIntegrationParams(useDiTau1MassConstraint, useDiTau2MassConstraint);
  prepareIntegrand();

  // CV: book histograms for evaluation of pT, eta, phi, mass and transverse mass of di-tau system
  if ( measuredTauLeptons_.size() == 4 ) {
    met_.SetX(measuredMETx);
    met_.SetY(measuredMETy);
    histogramAdapter_->setMeasurement(measuredTauLeptons_[0].p4(), measuredTauLeptons_[1].p4(), measuredTauLeptons_[2].p4(), measuredTauLeptons_[3].p4(), met_);
    histogramAdapter_->bookHistograms(measuredTauLeptons_[0].p4(), measuredTauLeptons_[1].p4(), measuredTauLeptons_[2].p4(), measuredTauLeptons_[3].p4(), met_);
  } else assert(0);

  if ( intAlgoFlag_ == kAlgoMarkovChain ) {
    integrateMC();
  } else if ( intAlgoFlag_ == kAlgoVAMP ) {
    integrateVAMP();
  } else assert(0); 
  
  clock_->Stop("<ClassicSVfit4tau::integrate>");
  numSeconds_cpu_ = clock_->GetCpuTime("<ClassicSVfit4tau::integrate>");
  numSeconds_real_ = clock_->GetRealTime("<ClassicSVfit4tau::integrate>");
  
  if ( verbosity_ >= 1 ) {
    clock_->Show("<ClassicSVfit4tau::integrate>");
  }
}

void ClassicSVfit4tau::integrateMC()
{
  if ( !intAlgo_ ) {
    initializeMCIntegrator();
  }
  
  double theIntegral, theIntegralErr;
  intAlgo_->integrate(&g_C, xl_, xh_, numDimensions_, theIntegral, theIntegralErr);
  isValidSolution_ = histogramAdapter_->isValidSolution();
  if ( isValidSolution_ ) {
    dihiggsMass_ = histogramAdapter_->getMass();
    dihiggsMassErr_ = histogramAdapter_->getMassErr();
    Lmax_ = getProbMax();
  } else {
    dihiggsMass_ = -1.;
    dihiggsMassErr_ = -1.;
    Lmax_ = -1.;
  }

  if ( verbosity_ >= 1 ) {
    if ( isValidSolution_ ) {
      std::cout << "found valid solution: mass = " << dihiggsMass_ << " +/- " << dihiggsMassErr_ << std::endl;
      std::cout << "(ditau1: mass = " << histogramAdapter_->ditau1()->getMass() << " +/- " << histogramAdapter_->ditau1()->getMassErr() << ","
		<< " ditau2: mass = " << histogramAdapter_->ditau2()->getMass() << " +/- " << histogramAdapter_->ditau2()->getMassErr() << ","
		<< " Lmax = " << Lmax_ << ")" << std::endl;
    } else {
      std::cout << "sorry, failed to find valid solution !!" << std::endl;
    }
  }

  if ( likelihoodFileName_ != "" ) {
    histogramAdapter_->writeHistograms(likelihoodFileName_);
  }
}

void ClassicSVfit4tau::integrateVAMP()
{
  ClassicSVfitIntegrand4tau* integrand4tau = dynamic_cast<ClassicSVfitIntegrand4tau*>(integrand_);
  assert(integrand4tau);
  
  assert(measuredTauLeptons_.size() == 4);
  double diHiggsVisMass = (measuredTauLeptons_[0].p4() + measuredTauLeptons_[1].p4() + measuredTauLeptons_[2].p4() + measuredTauLeptons_[3].p4()).mass();
  std::cout << "diHiggsVisMass = " << diHiggsVisMass << std::endl;
  double minMass = TMath::Max(diTau1MassConstraint_ + diTau2MassConstraint_, diHiggsVisMass);
  double maxMass = TMath::Max(1.e+4, 1.e+1*minMass);
  TH1* histogramMass = HistogramTools::makeHistogram_logBinWidth("ClassicSVfit4tau_integrateVAMP_histogramMass", minMass, maxMass, 1.025);

  int count = 0;
  assert(diTau1MassConstraint_ > 0. && diTau2MassConstraint_ > 0.);
  double m_test = 1.0125*minMass;
  double m = 0.;
  double pMax = 0.;
  bool skip_high_mass_tail = false;
  for ( int i = 0; i < 100 && (!skip_high_mass_tail); ++i ) {

    integrand4tau->setDiHiggsMassConstraint(m_test);

    int numDimensions   = numDimensions_;
    int numCallsGridOpt =  4000;
    int numCallsIntEval = 16000;
    double theIntegral, theIntegralErr;
    vamp_integrate_(&g_Fortran, xl_, xh_, &numDimensions, &numCallsGridOpt, &numCallsIntEval, &theIntegral, &theIntegralErr);
    
    double p = theIntegral;
    double pErr = theIntegralErr;
    if ( verbosity_ >= 1 ) {
      std::cout << "--> scan idx = " << i << ": m_test = " << m_test << ", p = " << p << " +/- " << pErr << " (pMax = " << pMax << ")" << std::endl;
    }    
    if ( p > pMax ) {
      m = m_test;      
      pMax = p;
      count = 0;
    } else {
      if ( p < (1.e-3*pMax) ) {
	++count;
	if ( count>= 5 ) {
	  skip_high_mass_tail = true;
	}
      } else {
	count = 0;
      }
    }
    double m_test_step = 0.025*m_test;
    int idxBin = histogramMass->FindBin(m_test);
    double binWidth = histogramMass->GetBinWidth(idxBin);
    histogramMass->SetBinContent(idxBin, p*binWidth);  // CV: multiply by bin-width to counteract division by bin-width 
    histogramMass->SetBinError(idxBin, pErr*binWidth); //     in functions 'extractHistogramProperties' (called by function 'extractLmax') and 'extractUncertainty'
    m_test += m_test_step;
  }
  double Lmax = HistogramTools::extractLmax(histogramMass);
  isValidSolution_ = ( Lmax > 0. ) ? true : false;
  if ( isValidSolution_ ) {
    dihiggsMass_ = m;
    dihiggsMassErr_ = HistogramTools::extractUncertainty(histogramMass);
    Lmax_ = Lmax;
  } else {
    dihiggsMass_ = -1.;
    dihiggsMassErr_ = -1.;
    Lmax_ = -1.;
  }

  if ( verbosity_ >= 1 ) {
    if ( isValidSolution_ ) {
      std::cout << "found valid solution: mass = " << dihiggsMass_ << " +/- " << dihiggsMassErr_ << ","
		<< " Lmax = " << Lmax_ << ")" << std::endl;
    } else {
      std::cout << "sorry, failed to find valid solution !!" << std::endl;
    }
  }

  delete histogramMass;
}

void ClassicSVfit4tau::setHistogramAdapter(classic_svFit::HistogramAdapterDiHiggs* histogramAdapter)
{
  if ( histogramAdapter_ ) delete histogramAdapter_;
  histogramAdapter_ = histogramAdapter;
}

classic_svFit::HistogramAdapterDiHiggs* ClassicSVfit4tau::getHistogramAdapter() const
{
  return histogramAdapter_;
}

double ClassicSVfit4tau::getMass() const
{
  return dihiggsMass_;
}

double ClassicSVfit4tau::getMassErr() const
{
  return dihiggsMassErr_;
}

double ClassicSVfit4tau::getLmax() const
{
  return Lmax_;
}
