// -*- C++ -*-
// Base off of code: ATLAS_2022_I2152933.cc
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

      declare(MissingMomentum(), "ETmiss");

      // Trimming function
      //_trimmer = fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05));

      // Book counter pointers and efficiency in scatter
      book(_c["JetsNum"], "JetsNumber");
      book(_c["bJetNum"], "BjetsNumber");
      Scatter1D s("/ATLAS_2023_BJETS/Efficiency");
      _s["Eff"] = registerAO(s);

      // Book jet count histograms
      book(_hjcount["NJets"],"Njets", 2, 1.5, 3.5);
      book(_hjcount["NBJets"],"NBjets", 4, -0.5, 3.5);
      book(_hjcount["fBJets"],"fBjets", 5, 0.0, 1.0);
      Scatter1D snj("/ATLAS_2023_BJETS/NJetsEff");
      Scatter1D snb("/ATLAS_2023_BJETS/NBJetsEff");
      _s["NJetsEff"] = registerAO(snj);
      _s["NBJetsEff"] = registerAO(snb);

      // Book jet structure histograms
      book(_h["nsj"],"Nsubjets", 10, -0.5, 9.5);
      book(_h["tau21"],"tau21", 40, 0.0, 2.0);
      book(_h["lha"],"LesHouchesAngularity", 35, 0.0, 0.7);
      book(_h["c2"],"C2correlation", 30, 0.0, 0.3);
      book(_h["d2"],"D2correlation", 60, 0.0, 4.0);
      book(_h["ecf2"],"EnergyCorrelationFunction2", 30, 0.0, 0.2);
      book(_h["ecf3"],"EnergyCorrelationFunction3", 30, 0.0, 0.01);
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
        return jet.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5);
      });

      /// Dilepton constain
      /*
      // Veto event if there are not exactly one electron and one muon
      if (electrons.size() != 1 || muons.size() != 1)  vetoEvent;

      // Veto event if the selected electron and muon are not OS
      if (electrons[0].charge() == muons[0].charge())  vetoEvent;

      // Veto event with dilepton mass < 15 GeV
      FourMomentum ll = electrons[0].momentum() + muons[0].momentum();
      if (ll.mass() <= 15*GeV ) vetoEvent;
      */

      // Veto event if there are not exactly 2 jets
      //if (jets.size() != 2) vetoEvent;

      // Veto event if there are no 2 jets or 3 jets
      if (jets.size() < 2 || jets.size() > 3)  vetoEvent;

      // Fill total jet counter before 2-bjet-selections
      _c["JetsNum"] -> fill(jets.size());
      _c["bJetNum"] -> fill(bjets.size());

      // Jet multiplicity histograms
      _hjcount["NJets"] -> fill(jets.size());
      _hjcount["NBJets"] -> fill(bjets.size());
      _hjcount["fBJets"] -> fill(bjets.size()/jets.size());

      // Veto even if there is no b-jets
      if (bjets.empty()) vetoEvent;
      
      // Veto event if the number of b-jets not exactly 2
      if (bjets.size() != 2)  vetoEvent;

      // Jet structure base off of code: ATLAS_2019_I1724098

      // Trim the retrieved clustered jets
      /*
      PseudoJets tr_jets;
      for (size_t i = 0; i < jets.size(); ++i) {
        tr_jets += _trimmer(jets[i]);
        tr_jets[tr_jets.size()-1].set_user_index(i);
      }
      
      size_t nBaseline = count(tr_jets, [](const Jet &j) { return j.pT() > 200*GeV && j.abseta() < 2.5; });
      if (nBaseline < 2)  return;

      ifilter_select(tr_jets, [](const PseudoJet &j) { return j.perp() > 450*GeV; });
      if (tr_jets.size() > 1)  tr_jets = sorted_by_pt(tr_jets);
      else if (tr_jets.empty())  return;

      if (abs(tr_jets[0].eta()) > 1.5)  return;
      */

      /// @todo Add a loop to compute and fill substructure for all b-jets
      /// @todo Then fill separate histograms for light (non-b) jets and b-jets
      for (unsigned i=0; i <= bjets.size(); i++) {

        // Extract the leading jet to perform calculation
        // Replace jets with tr_jets if trimming is applied
        const fastjet::PseudoJet &LJet = jets[i];

        // Declare jet structure variables
        double nsub, lha, ecf2, ecf3, c2, d2, tau21;
        lha = 0.0;

        // NSubjettiness
        JetDefinition subjet_def(fastjet::kt_algorithm, 0.2);
        ClusterSequence subjet_cs(LJet.constituents(), subjet_def);
        PseudoJets subjets = sorted_by_pt(subjet_cs.inclusive_jets(10.0));
        nsub = subjets.size();
        
        fastjet::contrib::Nsubjettiness nSub1(1, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta,Rcut));
        fastjet::contrib::Nsubjettiness nSub2(2, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta,Rcut));
        double tau1 = nSub1.result(LJet);
        double tau2 = nSub2.result(LJet);
        if(tau1 != 0) tau21 = tau2/tau1;
        else tau21 = -99;

        // LHA
        for (const PseudoJet& p : LJet.constituents()){
          double pt = p.pt();
          double theta = p.squared_distance(LJet);
          lha += pow(pt, 1.0) * pow(theta, 0.25);
        }
        double lterm = pow(LJet.pt(), 1.0) * pow(1.0, 0.5);
        if (lterm)  lha /= lterm;
        else        lha = -99;

        // C2
        fastjet::contrib::EnergyCorrelator ECF3(3,beta,fastjet::contrib::EnergyCorrelator::pt_R);
        fastjet::contrib::EnergyCorrelator ECF2(2,beta,fastjet::contrib::EnergyCorrelator::pt_R);
        fastjet::contrib::EnergyCorrelator ECF1(1,beta,fastjet::contrib::EnergyCorrelator::pt_R);

        double recf3 = ECF3(LJet);
        double recf2 = ECF2(LJet);
        double recf1 = ECF1(LJet);

        c2 = (recf2 != 0 ? recf3 * recf1 / (recf2*recf2) : -1);
        d2 = (recf2 != 0 ? recf3 * (recf1*recf1*recf1) /(recf2*recf2*recf2) : -1);

        ecf2 = (recf1 !=0 ? recf2 /(recf1*recf1) : -1);
        ecf3 = (recf1 !=0 ? recf3 / (recf1*recf1*recf1) : -1);

        // Fill bjet counter and structure histograms
        _h["nsj"]->fill(nsub);
        _h["tau21"]->fill(tau21);
        _h["c2"]->fill(c2);
        _h["d2"]->fill(d2);
        _h["lha"]->fill(lha);
        _h["ecf2"]->fill(ecf2);
        _h["ecf3"]->fill(ecf3);
      };
    };
    
    void finalize() {
      // Return the Efficiency dividing _c["bJetNum"] by _c["JetsNum"] and store into Scatter _s["Eff"]
      divide(_c["bJetNum"], _c["JetsNum"], _s["Eff"]);
      scale(_c, crossSection()/picobarn/sumOfWeights());

      scale(_hjcount["Njets"], 1/sumOfWeights());
      scale(_hjcount["NBjets"], 1/sumOfWeights());
      normalize(_hjcount["fBjets"]);
      //barchart(_hjcount["NJets"], _s["NJetsEff"]);
      //barchart(_hjcount["NBJets"], _s["NBJetsEff"]);

      scale(_h, crossSection()/picobarn/sumOfWeights());
    };

    /// @name Histograms, Counter, Scatter
    /// @{
    map<string, Histo1DPtr> _h, _hjcount;
    map<string, CounterPtr> _c;
    map<string, Scatter1DPtr> _s;
    /// @}

    //fastjet::Filter _trimmer;
    //size_t _mode;

  };

  DECLARE_RIVET_PLUGIN(ATLAS_2023_BJETS);

};
