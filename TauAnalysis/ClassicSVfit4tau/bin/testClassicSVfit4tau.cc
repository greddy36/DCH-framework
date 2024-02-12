
/**
   \class testClassicSVfit4tau testClassicSVfit4tau.cc "TauAnalysis/ClassicSVfit/bin/testClassicSVfit4tau.cc"
   \brief Basic example of the use of the standalone version of the "classic" SVfit algorithm, customized for the hh->4tau analysis
*/

#include "TauAnalysis/ClassicSVfit4tau/interface/ClassicSVfit4tau.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit4tau/interface/svFitHistogramAdapter4tau.h"
//#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBall2.h"

#include "TH1F.h"

using namespace classic_svFit;

int main(int argc, char* argv[])
{
  /*
     This is a single event for testing purposes.
  */

  // define MET
  double measuredMETx = -0.727501;
  double measuredMETy =  2.22865;

  // define MET covariance
  TMatrixD covMET(2, 2);
  covMET[0][0] = 415.234;
  covMET[1][0] =  66.9038;
  covMET[0][1] =  66.9038;
  covMET[1][1] = 329.008;

  // define lepton four vectors
  std::vector<MeasuredTauLepton> measuredTauLeptons;
  // decay products of first Higgs bosons
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay,  14.9606,  0.465224,  2.54672,   0.10566));    // tau -> muon decay (Pt, eta, phi, mass)
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, 67.5723,  1.00655,  -2.49902,   1.21648, 2)); // tau -> 1prong2pi0 hadronic decay (Pt, eta, phi, mass)
  // decay products of second Higgs bosons
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay,  21.4138,  1.43247,   0.0136936, 0.10566));    // tau -> muon decay (Pt, eta, phi, mass)
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, 47.3219, -0.0773025, 0.310924,  0.13957, 0)); // tau -> 1prong1pi0 hadronic decay (Pt, eta, phi, mass)
  /*
     tauDecayModes:  0 one-prong without neutral pions
                     1 one-prong with neutral pions
                    10 three-prong without neutral pions
  */

  int verbosity = 1;
  ClassicSVfit4tau svFitAlgo_MarkovChain(ClassicSVfit4tau::kAlgoMarkovChain, verbosity);
#ifdef USE_SVFITTF
  //HadTauTFCrystalBall2* hadTauTF = new HadTauTFCrystalBall2();
  //svFitAlgo_MarkovChain.setHadTauTF(hadTauTF);
  //svFitAlgo_MarkovChain.enableHadTauTF();
#endif
  // run with mass constraint for each tau pair to match measured mass (125.06 GeV) of SM-like Higgs boson   
  double massContraint = 125.06;
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
  svFitAlgo_MarkovChain.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
  bool isValidSolution_1stRun = svFitAlgo_MarkovChain.isValidSolution();
  double dihiggs_mass_1stRun = svFitAlgo_MarkovChain.getMass();
  double dihiggs_massErr_1stRun = svFitAlgo_MarkovChain.getMassErr();
  double dihiggs_transverseMass_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->getTransverseMass();
  double dihiggs_transverseMassErr_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->getTransverseMassErr();
  double ditau1_mass_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau1()->getMass();
  double ditau1_massErr_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau1()->getMassErr();
  double ditau2_mass_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau2()->getMass();
  double ditau2_massErr_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau2()->getMassErr();
  if ( isValidSolution_1stRun ) {
    std::cout << "found valid solution: mass = " << dihiggs_mass_1stRun << " +/- " << dihiggs_massErr_1stRun << " (expected value = 405.091 +/- 10.6929),"
              << " transverse mass = " << dihiggs_transverseMass_1stRun << " +/- " << dihiggs_transverseMassErr_1stRun << " (expected value = 361.319 +/- 12.8636)" << std::endl;
    std::cout << "(ditau1: mass = " << ditau1_mass_1stRun << " +/- " << ditau1_massErr_1stRun << ","
	      << " ditau2: mass = " << ditau2_mass_1stRun << " +/- " << ditau2_massErr_1stRun << ")" << std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }
  if ( kappa == 0. ) {
    if (std::abs((dihiggs_mass_1stRun - 405.091) / 405.091) > 0.001) return 1;
    if (std::abs((dihiggs_massErr_1stRun - 10.6929) / 10.6929) > 0.001) return 1;
    if (std::abs((dihiggs_transverseMass_1stRun - 361.319) / 361.319) > 0.001) return 1;
    if (std::abs((dihiggs_transverseMassErr_1stRun - 12.8636) / 12.8636) > 0.001) return 1;
  }
  
  // re-run without mass constraint for each tau pair;
  // this mode will set the mass of the second tau pair to match the mass of the first tau pair,
  // while the mass of the first tau pair is allowed to freely vary within the fit
  std::cout << "\n\nTesting Markov-Chain integration without fixed ditau mass constraint" << std::endl;
  svFitAlgo_MarkovChain.setLikelihoodFileName("testClassicSVfit4tau.root");
  svFitAlgo_MarkovChain.setDiTau1MassConstraint(-1.);
  svFitAlgo_MarkovChain.setDiTau2MassConstraint(-1.);
  svFitAlgo_MarkovChain.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
  bool isValidSolution_2ndRun = svFitAlgo_MarkovChain.isValidSolution();
  double dihiggs_mass_2ndRun = svFitAlgo_MarkovChain.getMass();
  double dihiggs_massErr_2ndRun = svFitAlgo_MarkovChain.getMassErr();
  double dihiggs_transverseMass_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->getTransverseMass();
  double dihiggs_transverseMassErr_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->getTransverseMassErr();
  double ditau1_mass_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau1()->getMass();
  double ditau1_massErr_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau1()->getMassErr();
  double ditau2_mass_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau2()->getMass();
  double ditau2_massErr_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo_MarkovChain.getHistogramAdapter())->ditau2()->getMassErr();
  if ( isValidSolution_2ndRun ) {
    std::cout << "found valid solution: mass = " << dihiggs_mass_2ndRun << " +/- " << dihiggs_massErr_2ndRun << " (expected value = 366.992 +/- 88.1767),"
              << " transverse mass = " << dihiggs_transverseMass_2ndRun << " +/- " << dihiggs_transverseMassErr_2ndRun << " (expected value = 343.908 +/- 71.5488)" << std::endl;
    std::cout << "(ditau1: mass = " << ditau1_mass_2ndRun << " +/- " << ditau1_massErr_2ndRun << ","
	      << " ditau2: mass = " << ditau2_mass_2ndRun << " +/- " << ditau2_massErr_2ndRun << ")" << std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }
  if ( kappa == 0. ) {
    if (std::abs((dihiggs_mass_2ndRun - 366.992) / 366.992) > 0.001) return 1;
    if (std::abs((dihiggs_massErr_2ndRun - 88.1767) / 88.1767) > 0.001) return 1;
    if (std::abs((dihiggs_transverseMass_2ndRun - 343.908) / 343.908) > 0.001) return 1;
    if (std::abs((dihiggs_transverseMassErr_2ndRun - 71.5488) / 71.5488) > 0.001) return 1;
  }

  ClassicSVfit4tau svFitAlgo_VAMP(ClassicSVfit4tau::kAlgoVAMP, verbosity);
#ifdef USE_SVFITTF
  //HadTauTFCrystalBall2* hadTauTF = new HadTauTFCrystalBall2();
  //svFitAlgo_VAMP.setHadTauTF(hadTauTF);
  //svFitAlgo_VAMP.enableHadTauTF();
#endif
  if ( kappa > 0. ) {
    svFitAlgo_VAMP.addLogM_fixed(true, kappa);
  } else {
    svFitAlgo_VAMP.addLogM_fixed(false);
  }

  std::cout << "\n\nTesting VAMP integration with ditau mass constraint set to " << massContraint << std::endl;
  svFitAlgo_VAMP.setDiTau1MassConstraint(massContraint);
  svFitAlgo_VAMP.setDiTau2MassConstraint(massContraint);  
  svFitAlgo_VAMP.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
  bool isValidSolution_3rdRun = svFitAlgo_VAMP.isValidSolution();
  double dihiggs_mass_3rdRun = svFitAlgo_VAMP.getMass();
  double dihiggs_massErr_3rdRun  = svFitAlgo_VAMP.getMassErr();
  if ( isValidSolution_3rdRun ) {
    std::cout << "found valid solution: mass = " << dihiggs_mass_3rdRun << " +/- " << dihiggs_massErr_3rdRun << " (expected value = 404.853 +/- 10.1509)" << std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }
  if ( kappa == 0. ) {
    if (std::abs((dihiggs_mass_3rdRun - 404.853) / 404.853) > 0.001) return 1;
    if (std::abs((dihiggs_massErr_3rdRun - 10.1509) / 10.1509) > 0.001) return 1;
  }

  std::cout << std::endl;
  std::cout << "************************************************************************************************************************" << std::endl;
  std::cout << "* If you use this code, please cite:                                                                                   *" << std::endl;
  std::cout << "*   K. Ehataeht, L. Marzola, C. Veelken,                                                                               *" << std::endl;
  std::cout << "*   \"Reconstruction of the mass of Higgs boson pairs in events with Higgs boson pairs decaying into four tau leptons\", *" << std::endl;
  std::cout << "*   arXiv:1809.06140                                                                                                   *" << std::endl;
  std::cout << "************************************************************************************************************************" << std::endl;
  std::cout << std::endl;

  return 0;
}
