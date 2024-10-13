// -*- C++ -*-
// Base off of code: ATLAS_2022_I2152933.cc

// Library from hepforge rivite
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"

#include "Rivet/Projections/ChargedLeptons.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"

#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/SmearedMET.hh"

// Setup Rivit anmespace
namespace Rivet {

  /// @brief bjet structure investigation
  class ATLAS_2023_BJETS : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2023_BJETS);

    // Declare constants ETA
    constexpr static double ETA = 3.5;
    // Declare constants beta and Rcut
    // Larger Rcut for captureing entire b-quark structure
    // Beta value of 1 for groomed jets
    constexpr static double beta = 1, Rcut = 1;

    /// Book histograms and initialise projections before the run
    void init() {
      // Test if Signal "SIG" (ttbar or diboson) or Background "BKG" (Z+jet) mode from cmnd file
      _signalmode = (getOption("MODE", "SIG") == "SIG");

      // Initialise and register projections
    	// Photons
    	PromptFinalState photons(Cuts::abspid == PID::PHOTON);

    	// Muons
    	Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;
    	PromptFinalState bare_mu(Cuts::abspid == PID::MUON, true);
    	DressedLeptons all_dressed_mu(photons, bare_mu, 0.1, lepton_cuts, true);
    	declare(all_dressed_mu, "muons");

    	// Electrons
    	PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, true);
    	DressedLeptons all_dressed_el(photons, bare_el, 0.1, lepton_cuts, true);
    	declare(all_dressed_el, "electrons");

    	// Jet forming
    	//const InvisibleFinalState neutrinos(true, true);
      // Ignore neutrinos (depends on purpose)

      // Probe finalstates in a larger angle than ETA
    	VetoedFinalState vfs(FinalState(Cuts::abseta < ETA + 0.4));
    	vfs.addVetoOnThisFinalState(all_dressed_el);
    	vfs.addVetoOnThisFinalState(all_dressed_mu);
    	//vfs.addVetoOnThisFinalState(neutrinos);

      FastJets jetfs(vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL);
    	declare(jetfs, "jets");

      // FinalState charged particles
      ChargedFinalState tracks(Cuts::pT > 0.5*GeV && Cuts::abseta < 2.5);
      declare(tracks, "tracks");

      // Create MET from undedected particles
      declare(MissingMomentum(), "ETmiss");

      // Book jet count and efficency scatter
      book(_c["JetsNum"], "JetsNum");
      book(_c["bJetsNum"], "bJetsNum");
      Scatter1D s("/ATLAS_2023_BJETS/Efficiency");
      _s["Eff"] = registerAO(s);

      // Book jet count and efficiency histograms and 2D scatter for barchart representation
      book(_hjcount["NJets"],"NJets", 2, 1.5, 3.5);
      book(_hjcount["NbJets"],"NbJets", 4, -0.5, 3.5);
      book(_hjcount["fbJets"],"fbJets", 5, 0.0, 1.0);
      Scatter2D nj("/ATLAS_2023_BJETS/NJetsEff");
      Scatter2D nb("/ATLAS_2023_BJETS/NbJetsEff");
      _s2["NJetsEff"] = registerAO(nj);
      _s2["NbJetsEff"] = registerAO(nb);

      // Book bjets substructure histograms
      book(_h["bJetsNSJ"],"bJetsNSJ", 10, -0.5, 9.5);
      book(_h["bJetsTAU21"],"bJetsTAU21", 40, 0.0, 2.0);
      book(_h["bJetsLHA"],"bJetsLHA", 35, 0.0, 0.7);
      book(_h["bJetsC2"],"bJetsC2", 30, 0.0, 0.3);
      book(_h["bJetsD2"],"bJetsD2", 60, 0.0, 4.0);
      book(_h["bJetsECF2"],"bJetsECF2", 30, 0.0, 0.2);
      book(_h["bJetsECF3"],"bJetsECF3", 30, 0.0, 0.01);

      // Book light jets substructure histograms
      book(_h["LJetsNSJ"],"LJetsNSJ", 10, -0.5, 9.5);
      book(_h["LJetsTAU21"],"LJetsTAU21", 40, 0.0, 2.0);
      book(_h["LJetsLHA"],"LJetsLHA", 35, 0.0, 0.7);
      book(_h["LJetsC2"],"LJetsC2", 30, 0.0, 0.3);
      book(_h["LJetsD2"],"LJetsD2", 60, 0.0, 4.0);
      book(_h["LJetsECF2"],"LJetsECF2", 30, 0.0, 0.2);
      book(_h["LJetsECF3"],"LJetsECF3", 30, 0.0, 0.01);
    };


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Retrieve dressed leptons, sorted by pT
      vector<DressedLepton> electrons = apply<DressedLeptons>(event, "electrons").dressedLeptons();
      vector<DressedLepton> muons = apply<DressedLeptons>(event, "muons").dressedLeptons();

      // Retrieve clustered jets, sorted by pT, with a minimum pT and maximum ETA cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < ETA);

      // Retrieve charge final state particles, sorted by pT
      const Particles& tracks = apply<ChargedFinalState>(event, "tracks").particlesByPt();

      // OVERLAP REMOVAL
      idiscardIfAnyDeltaRLess(muons, jets, 0.4);
      idiscardIfAnyDeltaRLess(electrons, jets, 0.4);

      // Select jets ghost-associated to B-hadrons with a certain fiducial selection
      Jets bjets = filter_select(jets, [](const Jet& jet) {
        return jet.bTagged(Cuts::pT > 5*GeV);
      });

      /// Dilepton constain (modify based on ttbar-dilep condition)
      /*
      // Veto event if there are not exactly one electron and one muon
      if (electrons.size() != 1 || muons.size() != 1)  vetoEvent;

      // Veto event if the selected electron and muon are not OS
      if (electrons[0].charge() == muons[0].charge())  vetoEvent;

      // Veto event with dilepton mass < 15 GeV
      FourMomentum ll = electrons[0].momentum() + muons[0].momentum();
      if (ll.mass() <= 15*GeV ) vetoEvent;
      */

      if (_signalmode) {
        // Selection for ttbar or diboson events
        // Veto event if there are no 2 jets or 3 jets
        if (jets.size() < 2 || jets.size() > 3)  vetoEvent;
      } else {
        // Selection for Z + jet events
        // Veto event if there is not at least 1 jet
        if (jets.size() < 1)  vetoEvent;
      };

      // Fill total jet counter before 2-bjet-selections
      _c["JetsNum"] -> fill(jets.size());
      _c["bJetsNum"] -> fill(bjets.size());

      // Jet multiplicity histograms
      _hjcount["NJets"] -> fill(jets.size());
      _hjcount["NbJets"] -> fill(bjets.size());
      _hjcount["fbJets"] -> fill(bjets.size()/double(jets.size()));

      // Filter jets within ETA 2.5 for analysis
      Jets FilterJets = filter_select(jets, Cuts::abseta < 2.5);
      Jets FilterbJets = filter_select(bjets, Cuts::abseta < 2.5);

      // Further selection for ttbar or diboson events
      if (_signalmode) {
        // Veto even if there is no b-jets
        if (FilterJets.empty()) vetoEvent;
        // Veto event if the number of b-jets not exactly 2
        if (FilterbJets.size() != 2)  vetoEvent;
      };

      // Jet structure base off of code: ATLAS_2019_I1724098

      /// Fill substructure histograms for b-jets and light (non-b) jets seperately
      for (const Jet& j : FilterJets) {
        // Test if a b-jet or not
        bool isBjet = j.bTagged(Cuts::pT > 5*GeV);

        const fastjet::PseudoJet &anaJet = j;

        // Declare jet structure variables
        double nsj, lha, ecf2, ecf3, c2, d2, tau21;
        lha = 0.0;

        // NSubjettiness
        JetDefinition subjet_def(fastjet::kt_algorithm, 0.2);
        ClusterSequence subjet_cs(anaJet.constituents(), subjet_def);
        PseudoJets subjets = sorted_by_pt(subjet_cs.inclusive_jets(10.0));
        nsj = subjets.size();
        
        fastjet::contrib::Nsubjettiness nSub1(1, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta,Rcut));
        fastjet::contrib::Nsubjettiness nSub2(2, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta,Rcut));
        double tau1 = nSub1.result(anaJet);
        double tau2 = nSub2.result(anaJet);
        if(tau1 != 0) tau21 = tau2/tau1;
        else tau21 = -99;

        // LHA
        for (const PseudoJet& p : anaJet.constituents()){
          double pt = p.pt();
          double theta = p.squared_distance(anaJet);
          lha += pow(pt, 1.0) * pow(theta, 0.25);
        }
        double lterm = pow(anaJet.pt(), 1.0) * pow(1.0, 0.5);
        if (lterm)  lha /= lterm;
        else        lha = -99;

        // C2
        fastjet::contrib::EnergyCorrelator ECF3(3,beta,fastjet::contrib::EnergyCorrelator::pt_R);
        fastjet::contrib::EnergyCorrelator ECF2(2,beta,fastjet::contrib::EnergyCorrelator::pt_R);
        fastjet::contrib::EnergyCorrelator ECF1(1,beta,fastjet::contrib::EnergyCorrelator::pt_R);

        double recf3 = ECF3(anaJet);
        double recf2 = ECF2(anaJet);
        double recf1 = ECF1(anaJet);

        c2 = (recf2 != 0 ? recf3 * recf1 / (recf2*recf2) : -1);
        d2 = (recf2 != 0 ? recf3 * (recf1*recf1*recf1) /(recf2*recf2*recf2) : -1);

        ecf2 = (recf1 !=0 ? recf2 /(recf1*recf1) : -1);
        ecf3 = (recf1 !=0 ? recf3 / (recf1*recf1*recf1) : -1);

        // Fill parallel histogram collections for b and light jets, using isBjet
        if (isBjet) {
          _h["bJetsNSJ"]->fill(nsj);
          _h["bJetsTAU21"]->fill(tau21);
          _h["bJetsC2"]->fill(c2);
          _h["bJetsD2"]->fill(d2);
          _h["bJetsLHA"]->fill(lha);
          _h["bJetsECF2"]->fill(ecf2);
          _h["bJetsECF3"]->fill(ecf3);
        } else {
          _h["LJetsNSJ"]->fill(nsj);
          _h["LJetsTAU21"]->fill(tau21);
          _h["LJetsC2"]->fill(c2);
          _h["LJetsD2"]->fill(d2);
          _h["LJetsLHA"]->fill(lha);
          _h["LJetsECF2"]->fill(ecf2);
          _h["LJetsECF3"]->fill(ecf3);
        };
      };
    };
    
    void finalize() {
      // Return the Efficiency dividing _c["bJetNum"] by _c["JetsNum"] and store into Scatter _s["Eff"]
      divide(_c["bJetsNum"], _c["JetsNum"], _s["Eff"]);
      scale(_c, crossSection()/picobarn/sumOfWeights());

      scale(_hjcount["NJets"], 1/sumOfWeights());
      scale(_hjcount["NbJets"], 1/sumOfWeights());
      normalize(_hjcount["fbJets"]);
      barchart(_hjcount["NJets"], _s2["NJetsEff"]);
      barchart(_hjcount["NbJets"], _s2["NbJetsEff"]);

      scale(_h, crossSection()/picobarn/sumOfWeights());
    };

    /// @name Histograms, Counter, Scatter
    /// @{
    map<string, Histo1DPtr> _h, _hjcount;
    map<string, CounterPtr> _c;
    map<string, Scatter1DPtr> _s;
    map<string, Scatter2DPtr> _s2;
    /// @}

    bool _signalmode = true;

  };

  DECLARE_RIVET_PLUGIN(ATLAS_2023_BJETS);

};
