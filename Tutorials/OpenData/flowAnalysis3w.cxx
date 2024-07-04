// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief Flow analysis.
///        Run as:
///        o2-analysis-timestamp --aod-file AO2D.root -b | o2-analysis-event-selection -b | o2-analysis-multiplicity-table -b | o2-analysis-centrality-table -b | o2-analysis-track-propagation -b | o2-analysis-trackextension -b | o2-analysis-trackselection -b | o2-analysis-pid-tpc-full -b | o2-analysis-pid-tof-full -b | o2-analysis-pid-tof-beta -b |  o2-analysis-pid-tof-base -b | o2-analysis-fdd-converter - b | o2-analysistutorial-flow-analysis3 -b
///        
///        o2-analysis-timestamp --aod-file AO2D.root --configuration json://config-file.json | o2-analysis-event-selection --configuration json://config-file.json | o2-analysis-multiplicity-table -b | o2-analysis-centrality-table -b | o2-analysis-track-propagation -b | o2-analysis-trackextension --configuration json://config-file.json | o2-analysis-trackselection -b | o2-analysis-pid-tpc-full -b | o2-analysis-pid-tof-full -b | o2-analysis-pid-tof-beta -b |  o2-analysis-pid-tof-base -b | o2-analysis-fdd-converter - b | o2-analysistutorial-flow-analysis3 -b
/// \author
/// \since

#include <vector>

#include <TF1.h>
#include <TH1.h>
#include <TH3.h>
#include <TFile.h>
#include <Framework/runDataProcessing.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/AnalysisDataModel.h>
#include <Common/DataModel/EventSelection.h>
#include <Common/CCDB/EventSelectionParams.h>
#include <Common/CCDB/TriggerAliases.h>
#include <Common/DataModel/Centrality.h>
#include <Common/DataModel/Multiplicity.h>
#include <Common/DataModel/TrackSelectionTables.h>
#include <Common/DataModel/PIDResponse.h>
#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPObject.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::evsel;

struct flow_base {

  using BCsInfos = soa::Join<aod::BCs, o2::aod::Timestamps>;
  using Colls = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>;

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  using FilteredTracks = soa::Filtered<TrackCandidates>;

  //Filter trackFilter = ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));
  Filter trackFilter = (requireGlobalTrackInFilter());

  Configurable<int> eventSelection{"eventSelection", 1, "event selection"};
  Configurable<bool> phiCut{"phiCut", false, "activate phi cut"};
  Configurable<bool> crsRowsFrcShCls{"crsRowsFrcShCls", false, "crsRowsFrcShCl"};
  Configurable<bool> hasQA{"hasQA", true, "Activate QA"};
  Configurable<float> vtxCut{"vtxCut", 10.0, "Z vertex cut"};
  Configurable<float> etaCut{"etaCut", 0.8, "Eta cut"};
  Configurable<float> dcazCut{"dcazCut", 3.2, "DCAz cut"};
  Configurable<float> dcaxyCut{"dcaxyCut", 2.4, "DCAxy cut"};
  Configurable<float> etaGap{"etaGap", 0.5, "Eta gap"};
  Configurable<int> noClus{"noClus", 70, "Number of clusters"};
  Configurable<int> nHarm{"nHarm", 2, "Number of harmonics"};
  Configurable<float> minPt{"minPt", 0.2, "Minimum pt"};
  Configurable<float> maxPt{"maxPt", 20.0, "Maximum pt"};
  Configurable<float> magField{"magField", 99999, "Configurable magnetic field; default CCDB will be queried"};
    Configurable<bool> rofCut{"rofCut", true, "activate rof cut"};

  Configurable<std::string> weightFile{"weightFile", "./histCenWght.root", "File with weights"};
  Configurable<std::string> weightTHname{"weightTHname", "hCenWg", "Hist with weights"};
  Configurable<std::string> weightPhiTHname{"weightPhiTHname", "hEtaPhiWeights", "Hist phi with weights"};
    Configurable<std::string> recQxAmTHname{"meanQxA", "fQxnAm", "Hist QxA mean"};
    Configurable<std::string> recQyAmTHname{"meanQyA", "fQynAm", "Hist QyA mean"};
    Configurable<std::string> recQxCmTHname{"meanQxC", "fQxnCm", "Hist QxC mean"};
    Configurable<std::string> recQyCmTHname{"meanQyC", "fQynCm", "Hist QyC mean"};
    Configurable<std::string> recQxAsTHname{"sigmaQxA", "fQxnAs", "Hist QxA sigma"};
    Configurable<std::string> recQyAsTHname{"sigmaQyA", "fQynAs", "Hist QyA sigma"};
    Configurable<std::string> recQxCsTHname{"sigmaQxC", "fQxnCs", "Hist QxC sigma"};
    Configurable<std::string> recQyCsTHname{"sigmaQyC", "fQynCs", "Hist QyC sigma"};
    

  TFile* wTF = nullptr;
  TH1D*  wH = nullptr;
  std::vector<TH3D*> wPhiH;
    std::vector<TH1D*> hQxAm;
    std::vector<TH1D*> hQyAm;
    std::vector<TH1D*> hQxCm;
    std::vector<TH1D*> hQyCm;
    std::vector<TH1D*> hQxAs;
    std::vector<TH1D*> hQyAs;
    std::vector<TH1D*> hQxCs;
    std::vector<TH1D*> hQyCs;


  TF1* fPhiCutLow = nullptr;
  TF1* fPhiCutHigh = nullptr;

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;

  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;

  //TF1* fMultMultPVCut = nullptr;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
    
    static constexpr int nVtxZ_bins = 10;
      static constexpr const char* QAEtaPhiAft[nVtxZ_bins] = {"QA/QAEtaPhiAft_0", "QA/QAEtaPhiAft_1", "QA/QAEtaPhiAft_2", "QA/QAEtaPhiAft_3", "QA/QAEtaPhiAft_4", "QA/QAEtaPhiAft_5", "QA/QAEtaPhiAft_6", "QA/QAEtaPhiAft_7", "QA/QAEtaPhiAft_8", "QA/QAEtaPhiAft_9"};
    
    
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void fillAPt(double trackpt, double cent, double vn, double sinHarm, double cosHarm, double wght) {
    histos.fill(HIST("PtA"), cent, trackpt);
    histos.fill(HIST("VnAPt"), trackpt, cent, vn, wght);
    histos.fill(HIST("SinnAPt"), trackpt, cent, sinHarm, wght);
    histos.fill(HIST("CosnAPt"), trackpt, cent, cosHarm, wght);
    }

  void fillCPt(double trackpt, double cent, double vn, double sinHarm, double cosHarm, double wght) {
    histos.fill(HIST("PtC"), cent, trackpt);
    histos.fill(HIST("VnCPt"), trackpt, cent, vn, wght);
    histos.fill(HIST("SinnCPt"), trackpt, cent, sinHarm, wght);
    histos.fill(HIST("CosnCPt"), trackpt, cent, cosHarm, wght);
    }


  void init(InitContext&) {
        AxisSpec axisVtxcounts{2, -0.5f, 1.5f, "Vtx info (0=no, 1=yes)"};
        AxisSpec axisZvert{120, -30.f, 30.f, "Vtx z (cm)"};
        AxisSpec axisCent{100, 0.f, 100.f, "centrality"};
        AxisSpec axisCentBins{{0, 5., 10., 20., 30., 40., 50., 60., 70., 80.}, "centrality percentile"};
        AxisSpec axisPtBins{{0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10., 13., 16., 20.}, "p_{T} (GeV/c)"};
      //AxisSpec axisQvec{220, -1.1, 1.1, "Q"};
      AxisSpec axisQvec{1000, -100., 100., "Q"};
      AxisSpec axisVtxZ{10, -10.f, 10.f, "Vtx z (cm)"};

      //AxisSpec axisEta{72, -0.9f, 0.9f, "#eta"};
      //AxisSpec axisPhi{144, 0.f, TMath::TwoPi(), "varphi"};
      AxisSpec axisEta{36, -0.9f, 0.9f, "#eta"};
        AxisSpec axisPhi{72, 0.f, TMath::TwoPi(), "varphi"};
        AxisSpec axisDCAz{100, -2.5f, 2.5f, "DCA_{z}"};
        AxisSpec axisDCAxy{100, -0.5f, 0.5f, "DCA_{xy}"};
        AxisSpec axisMultFw{1000, 0, 200000, "mult"};
        AxisSpec axisMult{1000, 0.f, 4500.f, "multiplicity"};
        AxisSpec axisPhiMod{100, 0.f, 0.4f, "phiMod"};
        
        histos.add("vtx", "Vtx info (0=no, 1=yes); Vtx; Counts", kTH1I, {axisVtxcounts});
        
        histos.add("vtxCutsBef", "Vtx distribution; Vtx z [cm]; Counts", kTH1F, {axisZvert});
        histos.add("multvsCentBef", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
        histos.add("multvsMultT0CBef", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
        histos.add("multvsmultV0ABef", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
        histos.add("multvsmultT0ABef", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
        histos.add("multvsmultTrkPVBef", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
        histos.add("multTrkPVvsCentBef", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
        histos.add("multTrkPVvsMultT0CBef", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
        histos.add("multTrkPVvsmultV0ABef", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
        histos.add("multTrkPVvsmultT0ABef", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
        histos.add("multV0AvsCentBef", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
        histos.add("multT0CvsmultT0ABef", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
        histos.add("multV0AvsmultT0ABef", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
        histos.add("multV0AvsmultT0CBef", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
        histos.add("multT0CvsCentBef", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});
        histos.add("centBef", "centrality T0C", kTH1F, {axisCent});
        
        histos.add("vtxCutsAft", "Vtx distribution; Vtx z [cm]; Counts", kTH1F, {axisZvert});
        histos.add("multvsCentAft", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
        histos.add("multvsMultT0CAft", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
        histos.add("multvsmultV0AAft", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
        histos.add("multvsmultT0AAft", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
        histos.add("multvsmultTrkPVAft", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
        histos.add("multTrkPVvsCentAft", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
        histos.add("multTrkPVvsMultT0CAft", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
        histos.add("multTrkPVvsmultV0AAft", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
        histos.add("multTrkPVvsmultT0AAft", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
        histos.add("multV0AvsCentAft", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
        histos.add("multT0CvsmultT0AAft", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
        histos.add("multV0AvsmultT0AAft", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
        histos.add("multV0AvsmultT0CAft", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
        histos.add("multT0CvsCentAft", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});
        histos.add("centAft", "centrality T0C", kTH1F, {axisCent});
        
        
        histos.add("res", "centrality percentile vs Resolution", kTProfile, {axisCentBins});
        histos.add("QxnA", "centrality percentile vs #LT Q_{x}^{nA} #GT", kTProfile, {axisCentBins});
        histos.add("QxnC", "centrality percentile vs #LT Q_{x}^{nC} #GT", kTProfile, {axisCentBins});
        histos.add("QynA", "centrality percentile vs #LT Q_{y}^{nA} #GT", kTProfile, {axisCentBins});
        histos.add("QynC", "centrality percentile vs #LT Q_{y}^{nC} #GT", kTProfile, {axisCentBins});
      
      histos.add("QxnAm", "centrality percentile vs #LT Q_{x}^{nA} #GT", kTProfile2D, {axisCent, axisVtxZ});
      histos.add("QxnCm", "centrality percentile vs #LT Q_{x}^{nC} #GT", kTProfile2D, {axisCent, axisVtxZ});
      histos.add("QynAm", "centrality percentile vs #LT Q_{y}^{nA} #GT", kTProfile2D, {axisCent, axisVtxZ});
      histos.add("QynCm", "centrality percentile vs #LT Q_{y}^{nC} #GT", kTProfile2D, {axisCent, axisVtxZ});
      
      histos.add("QxnAs", "centrality percentile vs Q_{x}^{nA}", kTH3F, {axisCent, axisQvec, axisVtxZ});
      histos.add("QynAs", "centrality percentile vs Q_{y}^{nA}", kTH3F, {axisCent, axisQvec, axisVtxZ});
      histos.add("QxnCs", "centrality percentile vs Q_{x}^{nC}", kTH3F, {axisCent, axisQvec, axisVtxZ});
      histos.add("QynCs", "centrality percentile vs Q_{y}^{nC}", kTH3F, {axisCent, axisQvec, axisVtxZ});
      
      
      histos.add("QxnAmCor", "centrality percentile vs #LT Q_{x}^{nA} #GT", kTProfile2D, {axisCent, axisVtxZ});
      histos.add("QxnCmCor", "centrality percentile vs #LT Q_{x}^{nC} #GT", kTProfile2D, {axisCent, axisVtxZ});
      histos.add("QynAmCor", "centrality percentile vs #LT Q_{y}^{nA} #GT", kTProfile2D, {axisCent, axisVtxZ});
      histos.add("QynCmCor", "centrality percentile vs #LT Q_{y}^{nC} #GT", kTProfile2D, {axisCent, axisVtxZ});
      
      histos.add("QxnAsCor", "centrality percentile vs Q_{x}^{nA}", kTH3F, {axisCent, axisQvec, axisVtxZ});
      histos.add("QynAsCor", "centrality percentile vs Q_{y}^{nA}", kTH3F, {axisCent, axisQvec, axisVtxZ});
      histos.add("QxnCsCor", "centrality percentile vs Q_{x}^{nC}", kTH3F, {axisCent, axisQvec, axisVtxZ});
      histos.add("QynCsCor", "centrality percentile vs Q_{y}^{nC}", kTH3F, {axisCent, axisQvec, axisVtxZ});
      
      
        
        histos.add("VnAPt", "v_{n} A", kTProfile2D, {{axisPtBins}, {axisCentBins}});
        histos.add("VnCPt", "v_{n} C", kTProfile2D, {{axisPtBins}, {axisCentBins}});
        histos.add("SinnAPt", "#LT sin(n*#phi) #GT A", kTProfile2D, {{axisPtBins}, {axisCentBins}});
        histos.add("SinnCPt", "#LT sin(n*#phi) #GT C", kTProfile2D, {{axisPtBins}, {axisCentBins}});
        histos.add("CosnAPt", "#LT cos(n*#phi) #GT A", kTProfile2D, {{axisPtBins}, {axisCentBins}});
        histos.add("CosnCPt", "#LT cos(n*#phi) #GT C", kTProfile2D, {{axisPtBins}, {axisCentBins}});
        histos.add("PtA", "p_{T} A", kTH2F, {{axisCentBins}, {axisPtBins}});
        histos.add("PtC", "p_{T} C", kTH2F, {{axisCentBins}, {axisPtBins}});
        

        histos.add("QA/QAEtaPhi", "#eta #varphi", kTH3F, {{axisEta}, {axisPhi}, {axisCentBins}});
        histos.add("QA/QAEtaPhiAft", "#eta #varphi (after cuts)", kTH3F, {{axisEta}, {axisPhi}, {axisCentBins}});
        histos.add("QA/QADCAz", "DCAz", kTH2F, {{axisDCAz}, {axisCentBins}});
        histos.add("QA/QADCAxy", "DCAxy", kTH2F, {{axisDCAxy}, {axisCentBins}});
        histos.add("QA/QADCAzAft", "DCAz (after cuts)", kTH2F, {{axisDCAz}, {axisCentBins}});
        histos.add("QA/QADCAxyAft", "DCAxy (after cuts)", kTH2F, {{axisDCAxy}, {axisCentBins}});
        histos.add("QA/QAPhiWAft", "#varphi (after cuts)", kTH2F, {{axisCentBins}, {axisPhi}});
        
        histos.add("QA/QAPhiModPtBef", "PhiMod (before cuts)", kTH2F, {{axisPtBins}, {axisPhiMod}});
        histos.add("QA/QAPhiModPtAft", "PhiMod (after cuts)", kTH2F, {{axisPtBins}, {axisPhiMod}});
        
      
      for (Int_t i = 0 ; i < 10; i++){
          histos.add(QAEtaPhiAft[i], "#eta #varphi p_{T} (after cuts)", kTH3F, {{axisEta}, {axisPhi}, {axisCentBins}});
      }

      
      //lhc23zzg_apass2
      fPhiCutLow = new TF1("fPhiCutLow", "0.06/x+pi/18.0-0.06", 0, 100);
      fPhiCutHigh = new TF1("fPhiCutHigh", "0.1/x+pi/18.0+0.06", 0, 100);
      
      
      /*
      //zzg pas3
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(3275.37, -122.097, 1.98076, -0.0171413, 6.4579e-05, 149.824, -1.64183, -0.029711, 0.000611661, -3.07289e-06);
      
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(3275.37, -122.097, 1.98076, -0.0171413, 6.4579e-05, 149.824, -1.64183, -0.029711, 0.000611661, -3.07289e-06);
      
        
      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutLow->SetParameters(1561.74, -43.3343, 0.392998, -0.00113134, 156.515, -3.67186, 0.0505422, -0.00060284, 3.2461e-06);
      
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutHigh->SetParameters(1561.74, -43.3343, 0.392998, -0.00113134, 156.515, -3.67186, 0.0505422, -0.00060284, 3.2461e-06);
      */
      
      
      //zzh pas3
      /*
      //old cuts Igor
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(3307.55, -122.1, 1.97653, -0.0172405, 6.57892e-05, 147.955, -2.40658, 0.00626944, 7.19406e-05, -3.92605e-07);
      
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(3307.55, -122.1, 1.97653, -0.0172405, 6.57892e-05, 147.955, -2.40658, 0.00626944, 7.19406e-05, -3.92605e-07);
      
        
      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutLow->SetParameters(1770.02, -49.5537, 0.460941, -0.00140622, 105.477, -1.58301, 0.011655, -0.000190804, 1.36003e-06);
      
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutHigh->SetParameters(1770.02, -49.5537, 0.460941, -0.00140622, 105.477, -1.58301, 0.011655, -0.000190804, 1.36003e-06);
      
      
      //new cuts Igor
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(3311.9, -121.502, 1.95845, -0.0171276, 6.59171e-05, 144.697, -3.1424, 0.0410312, -0.000433032, 2.00146e-06);
      
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(3311.9, -121.502, 1.95845, -0.0171276, 6.59171e-05, 144.697, -3.1424, 0.0410312, -0.000433032, 2.00146e-06);
      
        
      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutLow->SetParameters(1844.3, -51.4177, 0.476717, -0.00145503, 69.9419, -0.421047, -0.00841163, 5.49217e-05, 3.40521e-08);
      
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutHigh->SetParameters(1844.3, -51.4177, 0.476717, -0.00145503, 69.9419, -0.421047, -0.00841163, 5.49217e-05, 3.40521e-08);
       */
      
      
      //zzh apass4
      /*
      //test - new cuts Igor
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(3268.75, -119.296, 1.89914, -0.0163405, 6.2029e-05, 151.18, -3.22006, 0.0341881, -0.000265065, 9.47782e-07);
      
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(3268.75, -119.296, 1.89914, -0.0163405, 6.2029e-05, 151.18, -3.22006, 0.0341881, -0.000265065, 9.47782e-07);
      
        
      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutLow->SetParameters(1583.95, -42.1168, 0.354222, -0.000874759, 70.118, -0.629755, -0.00389294, 2.646e-05, 5.54674e-08);
      
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutHigh->SetParameters(1583.95, -42.1168, 0.354222, -0.000874759, 70.118, -0.629755, -0.00389294, 2.646e-05, 5.54674e-08);
      */
      
      
      //test3
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(3264.28, -118.681, 1.87824, -0.0160638, 6.07412e-05, 155.703, -3.73923, 0.0512274, -0.000490943, 2.02637e-06);
      
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(3264.28, -118.681, 1.87824, -0.0160638, 6.07412e-05, 155.703, -3.73923, 0.0512274, -0.000490943, 2.02637e-06);
      
        
      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutLow->SetParameters(2002.87, -60.3593, 0.631186, -0.00230718, 100.991, -2.25966, 0.0364445, -0.000469859, 2.41897e-06);
      
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutHigh->SetParameters(2002.87, -60.3593, 0.631186, -0.00230718, 100.991, -2.25966, 0.0364445, -0.000469859, 2.41897e-06);
      
        
        ccdb->setURL("http://alice-ccdb.cern.ch");
        ccdb->setCaching(true);
        ccdb->setLocalObjectValidityChecking();
        
        long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        ccdb->setCreatedNotAfter(now); // TODO must become global parameter from the train creation time
        
        std::string weight_file_name = weightFile.value;
        wTF = TFile::Open(weight_file_name.c_str(), "READ");

        std::string weight_histogram_name = weightTHname.value;
        std::string weightPhi_histogram_name = weightPhiTHname.value;
      std::string recMeanQxA_name = recQxAmTHname.value;
      std::string recMeanQyA_name = recQyAmTHname.value;
      std::string recMeanQxC_name = recQxCmTHname.value;
      std::string recMeanQyC_name = recQyCmTHname.value;
      std::string recSigmaQxA_name = recQxAsTHname.value;
      std::string recSigmaQyA_name = recQyAsTHname.value;
      std::string recSigmaQxC_name = recQxCsTHname.value;
      std::string recSigmaQyC_name = recQyCsTHname.value;

        if (wTF){
            wTF->GetObject(weight_histogram_name.c_str(), wH);

            for (int i = 0; i < 10; i++) {
                
                TH3D* tPPhi = nullptr;
                wTF->GetObject(Form("%s_%d", weightPhi_histogram_name.c_str(), i), tPPhi);
                wPhiH.push_back(tPPhi);
                
                TH1D* tPMeanQxA = nullptr;
                wTF->GetObject(Form("%s_%d", recMeanQxA_name.c_str(), i), tPMeanQxA);
                hQxAm.push_back(tPMeanQxA);
                
                TH1D* tPMeanQyA = nullptr;
                wTF->GetObject(Form("%s_%d", recMeanQyA_name.c_str(), i), tPMeanQyA);
                hQyAm.push_back(tPMeanQyA);
                
                TH1D* tPMeanQxC = nullptr;
                wTF->GetObject(Form("%s_%d", recMeanQxC_name.c_str(), i), tPMeanQxC);
                hQxCm.push_back(tPMeanQxC);
                
                TH1D* tPMeanQyC = nullptr;
                wTF->GetObject(Form("%s_%d", recMeanQyC_name.c_str(), i), tPMeanQyC);
                hQyCm.push_back(tPMeanQyC);
                
                
                TH1D* tPSigmaQxA = nullptr;
                wTF->GetObject(Form("%s_%d", recSigmaQxA_name.c_str(), i), tPSigmaQxA);
                hQxAs.push_back(tPSigmaQxA);
                
                TH1D* tPSigmaQyA = nullptr;
                wTF->GetObject(Form("%s_%d", recSigmaQyA_name.c_str(), i), tPSigmaQyA);
                hQyAs.push_back(tPSigmaQyA);
                
                TH1D* tPSigmaQxC = nullptr;
                wTF->GetObject(Form("%s_%d", recSigmaQxC_name.c_str(), i), tPSigmaQxC);
                hQxCs.push_back(tPSigmaQxC);
                
                TH1D* tPSigmaQyC = nullptr;
                wTF->GetObject(Form("%s_%d", recSigmaQyC_name.c_str(), i), tPSigmaQyC);
                hQyCs.push_back(tPSigmaQyC);
                
            }
        }
    }

  int getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPObject* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

    
  void process(Colls::iterator const& collision, BCsInfos const& bcs, FilteredTracks const& tracks)
  {

      if ((eventSelection == 1) && (!collision.sel8())) {
        // LOGF(info, "Collision index : %d skipped not selected", collision.index());
        return;
          }
      
      
    // new cut to remove collisions close to ROF
      if (rofCut){
          if (!collision.selection_bit(kNoITSROFrameBorder)){
              //close to rof
              return;
          }
          
          if (!collision.selection_bit(kNoSameBunchPileup)){
              //rejects collisions which are associated with the same "found-by-T0" bunch crossing
              return;
          }
          
          if (!collision.selection_bit(kIsGoodZvtxFT0vsPV)){
              //removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
              return;
          }
          
          if (!collision.selection_bit(kIsVertexITSTPC)){
              //selects collisions with at least one ITS-TPC track, and thus rejects vertices built from ITS-only tracks
              return;
          }
          
          if (!collision.selection_bit(kNoTimeFrameBorder)) {
              return;
          }
          
                    
          //new cuts
          int occupancy = collision.trackOccupancyInTimeRange();
          if (occupancy < 0 || occupancy >= 500)
              return;
          
          if (!collision.selection_bit(kNoCollInTimeRangeStandard))
              return;
          
      }
      
      
    //new cut: remove collisions with TRD trigger
    if (collision.alias_bit(kTVXinTRD)) {
        // TRD triggered
        return;
        }



    float zvtx = -999;
    if (collision.numContrib() > 1) {
        zvtx = collision.posZ();
        float zRes = TMath::Sqrt(collision.covZZ());
        if (zRes > 0.25 && collision.numContrib() < 20) { zvtx = -999; }
        }

    // if (zvtx < -990) then 0 else 1
    int is_vtx = (zvtx < -990) ? 0 : 1 ;
    histos.fill(HIST("vtx"), is_vtx);

    auto bc = collision.bc_as<BCsInfos>();
    auto field = (magField == 99999) ? getMagneticField(bc.timestamp()) : magField ;

    auto t0cCentr = collision.centFT0C();
      if (t0cCentr >= 80. || t0cCentr < 0) { return; }
      
    
      
    Short_t cen = -1;
    if (t0cCentr < 5) {
        cen = 0;
        }
    else if (t0cCentr >= 5 && t0cCentr < 10) {
        cen = 1;
        }
    else if (t0cCentr >= 10 && t0cCentr < 20) {
        cen = 2;
        }
    else if (t0cCentr >= 20 && t0cCentr < 30) {
        cen = 3;
        }
    else if (t0cCentr >= 30 && t0cCentr < 40) {
        cen = 4;
        }
    else if (t0cCentr >= 40 && t0cCentr < 50) {
        cen = 5;
        }
    else if (t0cCentr >= 50 && t0cCentr < 60) {
        cen = 6;
        } 
    else if (t0cCentr >= 60 && t0cCentr < 70){
        cen = 7;
        }
    else if (t0cCentr >= 70 && t0cCentr < 80){
        cen = 8;
        }

    if (cen < 0) { return; }
      
    // Int_t t0cCen = Int_t(t0cCentr); // warnings throw errors now

    // Use weigts TH1D
    auto t0cCentrW = 1.;
    if (wH){ t0cCentrW = wH->GetBinContent(wH->FindBin(t0cCentr)); }

    auto multV0A = collision.multFV0A();
    auto multT0A = collision.multFT0A();
    auto multT0C = collision.multFT0C();
    auto multNTracksPV = collision.multNTracksPV();

    Int_t multTrk = tracks.size();

    histos.fill(HIST("vtxCutsBef"), zvtx);
    histos.fill(HIST("multvsCentBef"), t0cCentr, multTrk);
    histos.fill(HIST("multvsMultT0CBef"), multT0C, multTrk);
    histos.fill(HIST("multvsmultV0ABef"), multV0A, multTrk);
    histos.fill(HIST("multvsmultT0ABef"), multT0A, multTrk);
    histos.fill(HIST("multvsmultTrkPVBef"), multNTracksPV, multTrk);
    histos.fill(HIST("multTrkPVvsCentBef"), t0cCentr, multNTracksPV);
    histos.fill(HIST("multTrkPVvsMultT0CBef"), multT0C, multNTracksPV);
    histos.fill(HIST("multTrkPVvsmultV0ABef"), multV0A, multNTracksPV);
    histos.fill(HIST("multTrkPVvsmultT0ABef"), multT0A, multNTracksPV);
    histos.fill(HIST("multV0AvsCentBef"), t0cCentr, multV0A);
    histos.fill(HIST("multT0CvsmultT0ABef"), multT0A, multT0C);
    histos.fill(HIST("multV0AvsmultT0ABef"), multT0A, multV0A);
    histos.fill(HIST("multV0AvsmultT0CBef"), multT0C, multV0A);
    histos.fill(HIST("multT0CvsCentBef"), t0cCentr, multT0C);
    histos.fill(HIST("centBef"), t0cCentr);

    if (TMath::Abs(zvtx) >= vtxCut) { return; }


    if (multNTracksPV < fMultPVCutLow->Eval(t0cCentr)) { return; }
    if (multNTracksPV > fMultPVCutHigh->Eval(t0cCentr)) { return; }

    if (multTrk < fMultCutLow->Eval(t0cCentr)) { return; }
    if (multTrk > fMultCutHigh->Eval(t0cCentr)) { return; }

      
    //new cut
    //if (multTrk < fMultMultPVCut->Eval(multNTracksPV))
    //if (multTrk > fMultMultPVCut->Eval(multNTracksPV)) { return; }
      
            
      Short_t zvt = -1;
      if (zvtx > -10. && zvtx < -8.)
          zvt = 0;
      else if (zvtx >= -8. && zvtx < -6.)
          zvt = 1;
      else if (zvtx >= -6. && zvtx < -4.)
          zvt = 2;
      else if (zvtx >= -4. && zvtx < -2.)
          zvt = 3;
      else if (zvtx >= -2. && zvtx < 0.)
          zvt = 4;
      else if (zvtx >= 0. && zvtx < 2.)
          zvt = 5;
      else if (zvtx >= 2. && zvtx < 4.)
          zvt = 6;
      else if (zvtx >= 4. && zvtx < 6.)
          zvt = 7;
      else if (zvtx >= 6. && zvtx < 8.)
          zvt = 8;
      else if (zvtx >= 8. && zvtx < 10.)
          zvt = 9;
      
      if (zvt < 0)
          return;

      
      histos.fill(HIST("vtxCutsAft"), zvtx);
      histos.fill(HIST("multvsCentAft"), t0cCentr, multTrk);
      histos.fill(HIST("multvsMultT0CAft"), multT0C, multTrk);
      histos.fill(HIST("multvsmultV0AAft"), multV0A, multTrk);
      histos.fill(HIST("multvsmultT0AAft"), multT0A, multTrk);
      histos.fill(HIST("multvsmultTrkPVAft"), multNTracksPV, multTrk);
      histos.fill(HIST("multTrkPVvsCentAft"), t0cCentr, multNTracksPV);
      histos.fill(HIST("multTrkPVvsMultT0CAft"), multT0C, multNTracksPV);
      histos.fill(HIST("multTrkPVvsmultV0AAft"), multV0A, multNTracksPV);
      histos.fill(HIST("multTrkPVvsmultT0AAft"), multT0A, multNTracksPV);
      histos.fill(HIST("multV0AvsCentAft"), t0cCentr, multV0A);
      histos.fill(HIST("multT0CvsmultT0AAft"), multT0A, multT0C);
      histos.fill(HIST("multV0AvsmultT0AAft"), multT0A, multV0A);
      histos.fill(HIST("multV0AvsmultT0CAft"), multT0C, multV0A);
      histos.fill(HIST("multT0CvsCentAft"), t0cCentr, multT0C);
      histos.fill(HIST("centAft"), t0cCentr, t0cCentrW);

    // process the tracks of a given collision
    Double_t QxnGapA = 0., QynGapA = 0.;
    Double_t QxnGapC = 0., QynGapC = 0.;

    Int_t multGapA = 0, multGapC = 0;

    for (auto& track : tracks) {
        
        //new cut to remove tracks without PV
        if (!track.isPVContributor()) {
            continue;
        }
        
        //if (!track.hasTRD())
        //    continue;
        
        //if (track.itsNCls() < 7)
        //    continue;

        Double_t trackpt = track.pt();
        Double_t tracketa = track.eta();
        Double_t trackdcaz = track.dcaZ();
        Double_t trackdcaxy = track.dcaXY();
        Double_t trackphi = track.phi();

      if (TMath::Abs(tracketa) >= etaCut || track.tpcNClsFound() < noClus || trackpt < minPt || trackpt >= maxPt || TMath::Abs(trackdcaz) >= dcazCut || TMath::Abs(trackdcaxy) >= dcaxyCut)
        { continue; }

        if (phiCut) {
          Double_t phimodn = trackphi;
          if (field < 0) // for negative polarity field
            phimodn = TMath::TwoPi() - phimodn;
          if (track.sign() < 0) // for negative charge
            phimodn = TMath::TwoPi() - phimodn;
          if (phimodn < 0)
            LOGF(warning, "phi < 0: %g", phimodn);

          phimodn += TMath::Pi() / 18.0; // to center gap in the middle
          phimodn = fmod(phimodn, TMath::Pi() / 9.0);

          if (phimodn < fPhiCutHigh->Eval(trackpt) && phimodn > fPhiCutLow->Eval(trackpt))
              continue; // reject track
            
            //if (phimodn > 0.15 && trackpt < 1.5)
            //    continue;
        }
        
        if (crsRowsFrcShCls) {
          Float_t nrowscr = track.tpcNClsCrossedRows();
          if (nrowscr < 120) { continue; }

          Float_t clsFind = track.tpcNClsFindable();
          if (clsFind <= 0) { continue; }

          if (track.tpcCrossedRowsOverFindableCls() < 0.9) { continue; }
          }
        
        //Double_t phiW = 1.0;
        Double_t phiW = wPhiH[zvt]->GetBinContent(wPhiH[zvt]->FindBin(tracketa, trackphi, t0cCentr));
        
        
        Double_t sinHarm = phiW*TMath::Sin(nHarm * trackphi);
        Double_t cosHarm = phiW*TMath::Cos(nHarm * trackphi);
        


        if (tracketa > etaGap) {
            QxnGapC += cosHarm;
            QynGapC += sinHarm;
            multGapC++;
        }

        if (tracketa < -etaGap) {
            QxnGapA += cosHarm;
            QynGapA += sinHarm;
            multGapA++;
        }
    }
     
      
      if (multGapA <= 0 || multGapC <= 0)
          return;
      
      
      Double_t multGapAf = multGapA;
      Double_t multGapCf = multGapC;

      /*
      Double_t QxnGapAcor = (QxnGapA - hQxAm[zvt]->GetBinContent(t0cCen+1))/hQxAs[zvt]->GetBinContent(t0cCen+1);
      Double_t QynGapAcor = (QynGapA - hQyAm[zvt]->GetBinContent(t0cCen+1))/hQyAs[zvt]->GetBinContent(t0cCen+1);
      Double_t QxnGapCcor = (QxnGapC - hQxCm[zvt]->GetBinContent(t0cCen+1))/hQxCs[zvt]->GetBinContent(t0cCen+1);
      Double_t QynGapCcor = (QynGapC - hQyCm[zvt]->GetBinContent(t0cCen+1))/hQyCs[zvt]->GetBinContent(t0cCen+1);
       */
      
      
      Double_t QxnGapAcor = QxnGapA;
      Double_t QynGapAcor = QynGapA;
      Double_t QxnGapCcor = QxnGapC;
      Double_t QynGapCcor = QynGapC;
      
       
      
      Double_t resGap = (QxnGapAcor * QxnGapCcor + QynGapAcor * QynGapCcor) / (multGapAf * multGapCf);
      histos.fill(HIST("res"), t0cCentr, resGap, t0cCentrW);

      histos.fill(HIST("QxnA"), t0cCentr, QxnGapAcor, t0cCentrW);
      histos.fill(HIST("QxnC"), t0cCentr, QxnGapCcor, t0cCentrW);

      histos.fill(HIST("QynA"), t0cCentr, QynGapAcor, t0cCentrW);
      histos.fill(HIST("QynC"), t0cCentr, QynGapCcor, t0cCentrW);
        
        
        histos.fill(HIST("QxnAm"), t0cCentr, zvtx, QxnGapA);
        histos.fill(HIST("QxnCm"), t0cCentr, zvtx, QxnGapC);

        histos.fill(HIST("QynAm"), t0cCentr, zvtx, QynGapA);
        histos.fill(HIST("QynCm"), t0cCentr, zvtx, QynGapC);
        
        histos.fill(HIST("QxnAs"), t0cCentr, QxnGapA, zvtx);
        histos.fill(HIST("QxnCs"), t0cCentr, QxnGapC, zvtx);

        histos.fill(HIST("QynAs"), t0cCentr, QynGapA, zvtx);
        histos.fill(HIST("QynCs"), t0cCentr, QynGapC, zvtx);
      
      
      histos.fill(HIST("QxnAmCor"), t0cCentr, zvtx, QxnGapAcor);
      histos.fill(HIST("QxnCmCor"), t0cCentr, zvtx, QxnGapCcor);

      histos.fill(HIST("QynAmCor"), t0cCentr, zvtx, QynGapAcor);
      histos.fill(HIST("QynCmCor"), t0cCentr, zvtx, QynGapCcor);
      
      histos.fill(HIST("QxnAsCor"), t0cCentr, QxnGapAcor, zvtx);
      histos.fill(HIST("QxnCsCor"), t0cCentr, QxnGapCcor, zvtx);

      histos.fill(HIST("QynAsCor"), t0cCentr, QynGapAcor, zvtx);
      histos.fill(HIST("QynCsCor"), t0cCentr, QynGapCcor, zvtx);


      
    for (auto& track : tracks) {
        
        //new cut to remove tracks without PV
        if (!track.isPVContributor()) {
            continue;
        }
        
        //if (!track.hasTRD())
        //    continue;
        
        //if (track.itsNCls() < 7)
        //    continue;

      Double_t trackpt = track.pt();
      Double_t tracketa = track.eta();
        Double_t trackdcaz = track.dcaZ();
        Double_t trackdcaxy = track.dcaXY();
        Double_t trackphi = track.phi();
        
        histos.fill(HIST("QA/QAEtaPhi"), tracketa, trackphi, t0cCentr);
        histos.fill(HIST("QA/QADCAz"), trackdcaz, t0cCentr);
        histos.fill(HIST("QA/QADCAxy"), trackdcaxy, t0cCentr);

      if (TMath::Abs(tracketa) >= etaCut || track.tpcNClsFound() < noClus || trackpt < minPt || trackpt >= maxPt || TMath::Abs(trackdcaz) >= dcazCut || TMath::Abs(trackdcaxy) >= dcaxyCut)
        continue;


      if (phiCut) {
        Double_t phimod = trackphi;

        if (field < 0) // for negative polarity field
          phimod = TMath::TwoPi() - phimod;

        if (track.sign() < 0) // for negative charge
          phimod = TMath::TwoPi() - phimod;

        if (phimod < 0) { LOGF(warning, "phi < 0: %g", phimod); }

        phimod += TMath::Pi() / 18.0; // to center gap in the middle
        phimod = fmod(phimod, TMath::Pi() / 9.0);

        histos.fill(HIST("QA/QAPhiModPtBef"), trackpt, phimod);

        if (phimod < fPhiCutHigh->Eval(trackpt) && phimod > fPhiCutLow->Eval(trackpt)) { continue; } // reject track
          
          //if (phimod > 0.15 && trackpt < 1.5)
          //    continue;

        histos.fill(HIST("QA/QAPhiModPtAft"), trackpt, phimod);
          
    }

      if (crsRowsFrcShCls) {
        Float_t nrowscr = track.tpcNClsCrossedRows();
        if (nrowscr < 120) { continue; }

        Float_t clsFind = track.tpcNClsFindable();
        if (clsFind <= 0) { continue; }

        if (track.tpcCrossedRowsOverFindableCls() < 0.9) { continue; }
        }

      histos.fill(HIST("QA/QAEtaPhiAft"), tracketa, trackphi, t0cCentr);
        
        switch (zvt) {
            case 0: { histos.fill(HIST(QAEtaPhiAft[0]), tracketa, trackphi, t0cCentr); break; }
            case 1: { histos.fill(HIST(QAEtaPhiAft[1]), tracketa, trackphi, t0cCentr); break; }
            case 2: { histos.fill(HIST(QAEtaPhiAft[2]), tracketa, trackphi, t0cCentr); break; }
            case 3: { histos.fill(HIST(QAEtaPhiAft[3]), tracketa, trackphi, t0cCentr); break; }
            case 4: { histos.fill(HIST(QAEtaPhiAft[4]), tracketa, trackphi, t0cCentr); break; }
            case 5: { histos.fill(HIST(QAEtaPhiAft[5]), tracketa, trackphi, t0cCentr); break; }
            case 6: { histos.fill(HIST(QAEtaPhiAft[6]), tracketa, trackphi, t0cCentr); break; }
            case 7: { histos.fill(HIST(QAEtaPhiAft[7]), tracketa, trackphi, t0cCentr); break; }
            case 8: { histos.fill(HIST(QAEtaPhiAft[8]), tracketa, trackphi, t0cCentr); break; }
            case 9: { histos.fill(HIST(QAEtaPhiAft[9]), tracketa, trackphi, t0cCentr); break; }
        }
        
      histos.fill(HIST("QA/QADCAzAft"), trackdcaz, t0cCentr);
      histos.fill(HIST("QA/QADCAxyAft"), trackdcaxy, t0cCentr);

        //Double_t phiWn = 1.;
      Double_t phiWn = wPhiH[zvt]->GetBinContent(wPhiH[zvt]->FindBin(tracketa, trackphi, t0cCentr));
      

      histos.fill(HIST("QA/QAPhiWAft"), t0cCentr, trackphi, phiWn);

      Double_t sinHarmn = phiWn*TMath::Sin(nHarm * trackphi);
      Double_t cosHarmn = phiWn*TMath::Cos(nHarm * trackphi);

      Double_t harmGapC = cosHarmn * QxnGapCcor + sinHarmn * QynGapCcor;
      Double_t harmGapA = cosHarmn * QxnGapAcor + sinHarmn * QynGapAcor;

      if (tracketa > etaGap && multGapA > 0) {
        Double_t vnC = harmGapA / multGapAf;
        fillCPt(trackpt, t0cCentr, vnC, sinHarmn, cosHarmn, t0cCentrW);
      }

      if (tracketa < -etaGap && multGapC > 0) {
        Double_t vnA = harmGapC / multGapCf;
        fillAPt(trackpt, t0cCentr, vnA, sinHarmn, cosHarmn, t0cCentrW);
      }

    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flow_base>(cfgc),
  };
}
