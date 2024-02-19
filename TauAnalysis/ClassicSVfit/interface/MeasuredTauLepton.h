#ifndef TauAnalysis_ClassicSVfit_MeasuredTauLepton_h
#define TauAnalysis_ClassicSVfit_MeasuredTauLepton_h

#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // Vector, LorentzVector

namespace classic_svFit
{
  class MeasuredTauLepton
  {
   public:
    /**
       \enum    MeasuredTauLepton::kDecayType
       \brief   enumeration of all tau decay types
    */
    enum kDecayType {
      kUndefinedDecayType,
      kTauToHadDecay,  /* < hadronic tau lepton decay                                                        */
      kTauToElecDecay, /* < tau lepton decay to electron                                                     */
      kTauToMuDecay,   /* < tau lepton decay to muon                                                         */
      kPrompt          /* < electron or muon directly originating from LFV Higgs boson decay                 */
    };

    MeasuredTauLepton();
    MeasuredTauLepton(int, double, double, double, double, int = -1);
    MeasuredTauLepton(const MeasuredTauLepton&);
    ~MeasuredTauLepton();

    /// return decay type of the tau lepton
    int type() const;

    /// return pt of the measured tau lepton in labframe
    double pt() const;
    /// return pseudo-rapidity of the measured tau lepton in labframe
    double eta() const;
    /// return azimuthal angle of the measured tau lepton in labframe
    double phi() const;
    /// return visible mass in labframe
    double mass() const;

    /// return visible energy in labframe
    double energy() const;
    /// return px of the measured tau lepton in labframe
    double px() const;
    /// return py of the measured tau lepton in labframe
    double py() const;
    /// return pz of the measured tau lepton in labframe
    double pz() const;

    /// return the measured tau lepton momentum in labframe
    double p() const;

    /// return decay mode of the reconstructed hadronic tau decay
    int decayMode() const;

    /// return the lorentz vector in the labframe
    LorentzVector p4() const;

    /// return the momentum vector in the labframe
    Vector p3() const;

    /// return auxiliary data-members to speed-up numerical computations
    double cosPhi_sinTheta() const;
    double sinPhi_sinTheta() const;
    double cosTheta() const;

    void roundToNdigits(unsigned int nDigis = 3);

    /// flag indicating if tau decayed leptonically or hadronically 
    bool isLeptonicTauDecay() const;
    bool isHadronicTauDecay() const;

    /// flag indicating electrons or muons directly originating from LFV Higgs boson decay
    bool isPrompt() const;

   protected:
    /// set visible momentum in all coordinates systems
    void initialize();

   private:
    /// decay type
    int type_;

    /// visible momentum in labframe (in polar coordinates)
    double pt_;
    double eta_;
    double phi_;
    double mass_;

    /// visible momentum in labframe (in cartesian coordinates)
    double energy_;
    double px_;
    double py_;
    double pz_;

    /// visible momentum in labframe (magnitude);
    double p_;

    /// decay mode (hadronic tau decays only)
    int decayMode_;

    /// visible momentum in labframe (four-vector)
    LorentzVector p4_;

    /// visible momentum in labframe
    Vector p3_;

    /// mass of visible tau decay products (recomputed to reduce rounding errors)
    double preciseVisMass_;

    /// auxiliary data-members to speed-up numerical computations
    double cosPhi_sinTheta_;
    double sinPhi_sinTheta_;
    double cosTheta_;
    bool isLeptonicTauDecay_;
    bool isHadronicTauDecay_;
    bool isPrompt_;
  };

  // auxiliary class for sorting MeasuredTauLeptons
  struct sortMeasuredTauLeptons
  {
    bool operator() (const MeasuredTauLepton& measuredTauLepton1, const MeasuredTauLepton& measuredTauLepton2);
  };
}

#endif
