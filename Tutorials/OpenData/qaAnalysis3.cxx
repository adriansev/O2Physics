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

#include <TF1.h>
#include <TH3.h>
#include <Framework/runDataProcessing.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/AnalysisDataModel.h>
#include <Common/DataModel/EventSelection.h>
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

struct qa3 {

    using BCsInfos = soa::Join<aod::BCs, o2::aod::Timestamps>;
    using Colls = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
    
    using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
    using FilteredTracks = soa::Filtered<TrackCandidates>;
    Filter trackFilter = (requireGlobalTrackInFilter());

    
    Configurable<int> eventSelection{"eventSelection", 1, "event selection"};
    Configurable<float> vtxCut{"vtxCut", 10.0, "Z vertex cut"};
    Configurable<float> etaCut{"etaCut", 0.8, "Eta cut"};
    Configurable<float> dcazCut{"dcazCut", 3.2, "DCAz cut"};
    Configurable<float> dcaxyCut{"dcaxyCut", 2.4, "DCAxy cut"};
    Configurable<int> noClus{"noClus", 70, "Number of clusters"};
    Configurable<float> minPt{"minPt", 0.2, "Minimum pt"};
    Configurable<float> maxPt{"maxPt", 20.0, "Maximum pt"};

        
    Service<o2::ccdb::BasicCCDBManager> ccdb;

    HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};


void init(InitContext&)
{
    AxisSpec axisVtxcounts{2, -0.5f, 1.5f, "Vtx info (0=no, 1=yes)"};
    AxisSpec axisZvert{120, -30.f, 30.f, "Vtx z (cm)"};
    AxisSpec axisPtBins{{0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10., 13., 16., 20.}, "p_{T} (GeV/c)"};

    AxisSpec axisEta{72, -0.9f, 0.9f, "#eta"};
    AxisSpec axisPhi{144, 0.f, TMath::TwoPi(), "varphi"};
    AxisSpec axisDCAz{700, -3.5f, 3.5f, "DCA_{z}"};
    AxisSpec axisDCAxy{500, -2.5f, 2.5f, "DCA_{xy}"};
    AxisSpec axisMultFw{1000, 0, 200000, "mult"};
    AxisSpec axisMult{1000, 0.f, 10000.f, "multiplicity"};
    AxisSpec axisQ{4, -2.0f, 2.0f, "charge"};
    AxisSpec axisNclITS{9, 0.f, 9.0f, "ncls_{ITS}"};
    AxisSpec axisNclTPC{160, 0.f, 160.f, "ncls_{TPC}"};

      histos.add("vtx", "Vtx info (0=no, 1=yes); Vtx; Counts", kTH1I, {axisVtxcounts});
      
      histos.add("vtxCutsBef", "Vtx distribution; Vtx z [cm]; Counts", kTH1F, {axisZvert});
      histos.add("multvsMultT0CBef", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
      histos.add("multvsmultV0ABef", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
      histos.add("multvsmultT0ABef", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
      histos.add("multvsmultTrkPVBef", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
      histos.add("multTrkPVvsMultT0CBef", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
      histos.add("multTrkPVvsmultV0ABef", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
      histos.add("multTrkPVvsmultT0ABef", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
      histos.add("multT0CvsmultT0ABef", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
      histos.add("multV0AvsmultT0ABef", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
      histos.add("multV0AvsmultT0CBef", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});

        
      histos.add("vtxCutsAft", "Vtx distribution; Vtx z [cm]; Counts", kTH1F, {axisZvert});
      histos.add("multvsMultT0CAft", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
      histos.add("multvsmultV0AAft", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
      histos.add("multvsmultT0AAft", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
      histos.add("multvsmultTrkPVAft", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
      histos.add("multTrkPVvsMultT0CAft", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
      histos.add("multTrkPVvsmultV0AAft", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
      histos.add("multTrkPVvsmultT0AAft", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
      histos.add("multT0CvsmultT0AAft", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
      histos.add("multV0AvsmultT0AAft", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
      histos.add("multV0AvsmultT0CAft", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
      histos.add("multITSvsMultITSTPCAft", " multiplicity ITS vs multiplicity ITS+TPC", kTH2F, {axisMult, axisMult});
      histos.add("multITSvsMultITSTPCNAft", " multiplicity ITS vs multiplicity ITS+TPC", kTH2F, {axisMult, axisMult});
    histos.add("multTPCvsMultITSTPCAft", " multiplicity TPC vs multiplicity ITS+TPC", kTH2F, {axisMult, axisMult});
    histos.add("multTRDvsMultITSTPCAft", " multiplicity TRD vs multiplicity ITS+TPC", kTH2F, {axisMult, axisMult});
    histos.add("multTOFvsMultITSTPCAft", " multiplicity TOF vs multiplicity ITS+TPC", kTH2F, {axisMult, axisMult});

    
      histos.add("QAEtaPhi", "#eta vs #varphi", kTH3F, {{axisEta}, {axisPhi}, {axisQ}});
      histos.add("QAEtaPhiAft", "#eta vs #varphi (after cuts)", kTH3F, {{axisEta}, {axisPhi}, {axisQ}});
      histos.add("QADCAz", "DCAz", kTH3F, {{axisDCAz}, {axisPtBins}, {axisQ}});
      histos.add("QADCAxy", "DCAxy", kTH3F, {{axisDCAxy}, {axisPtBins}, {axisQ}});
      histos.add("QADCAzAft", "DCAz (after cuts)", kTH3F, {{axisDCAz}, {axisPtBins}, {axisQ}});
      histos.add("QADCAxyAft", "DCAxy (after cuts)", kTH3F, {{axisDCAxy}, {axisPtBins}, {axisQ}});
    
        histos.add("QANclsITS", "n_{cls}^{ITS} vs p_{T}", kTH3F, {{axisNclITS}, {axisPtBins}, {axisQ}});
        histos.add("QANclsITSAft", "n_{cls}^{ITS} vs p_{T}", kTH3F, {{axisNclITS}, {axisPtBins}, {axisQ}});
    
    histos.add("QANclsTPC", "n_{cls}^{TPC} vs p_{T}", kTH3F, {{axisNclTPC}, {axisPtBins}, {axisQ}});
    histos.add("QANclsTPCAft", "n_{cls}^{TPC} vs p_{T}", kTH3F, {{axisNclTPC}, {axisPtBins}, {axisQ}});
    
    histos.add("QANcrRowTPC", "n_{cls}^{TPC} vs p_{T}", kTH3F, {{axisNclTPC}, {axisPtBins}, {axisQ}});
    histos.add("QANcrRowTPCAft", "n_{cls}^{TPC} vs p_{T}", kTH3F, {{axisNclTPC}, {axisPtBins}, {axisQ}});
    
    histos.add("QAPtEta", "p_{T} vs #eta", kTH3F, {{axisEta}, {axisPtBins}, {axisQ}});
    histos.add("QAPtEtaAft", "p_{T} vs #eta", kTH3F, {{axisEta}, {axisPtBins}, {axisQ}});
    

      ccdb->setURL("http://alice-ccdb.cern.ch");
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();

      long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
      ccdb->setCreatedNotAfter(now); // TODO must become global parameter from the train creation time
  }


void process(Colls::iterator const& collision, BCsInfos const& bcs, FilteredTracks const& tracks)
{
    
    
    if ((eventSelection == 1) && (!collision.sel8())) {
        // LOGF(info, "Collision index : %d skipped not selected", collision.index());
        return;
    }
     
    
    float zvtx = -999;
    if (collision.numContrib() > 1) {
          
        zvtx = collision.posZ();
          
        float zRes = TMath::Sqrt(collision.covZZ());
        if (zRes > 0.25 && collision.numContrib() < 20)
            zvtx = -999.;
          
    }
      
    if (zvtx < -990)
        histos.fill(HIST("vtx"), 0);
    else
        histos.fill(HIST("vtx"), 1);
      
      
      auto multV0A = collision.multFV0A();
      auto multT0A = collision.multFT0A();
      auto multT0C = collision.multFT0C();
      auto multNTracksPV = collision.multNTracksPV();
      
      Int_t multTrk = tracks.size();
      
      
      histos.fill(HIST("vtxCutsBef"), zvtx);
      histos.fill(HIST("multvsMultT0CBef"), multT0C, multTrk);
      histos.fill(HIST("multvsmultV0ABef"), multV0A, multTrk);
      histos.fill(HIST("multvsmultT0ABef"), multT0A, multTrk);
      histos.fill(HIST("multvsmultTrkPVBef"), multNTracksPV, multTrk);
      histos.fill(HIST("multTrkPVvsMultT0CBef"), multT0C, multNTracksPV);
      histos.fill(HIST("multTrkPVvsmultV0ABef"), multV0A, multNTracksPV);
      histos.fill(HIST("multTrkPVvsmultT0ABef"), multT0A, multNTracksPV);
      histos.fill(HIST("multT0CvsmultT0ABef"), multT0A, multT0C);
      histos.fill(HIST("multV0AvsmultT0ABef"), multT0A, multV0A);
      histos.fill(HIST("multV0AvsmultT0CBef"), multT0C, multV0A);
      

      
      if (TMath::Abs(zvtx) > vtxCut)
          return;
      
      
      histos.fill(HIST("vtxCutsAft"), zvtx);
      histos.fill(HIST("multvsMultT0CAft"), multT0C, multTrk);
      histos.fill(HIST("multvsmultV0AAft"), multV0A, multTrk);
      histos.fill(HIST("multvsmultT0AAft"), multT0A, multTrk);
      histos.fill(HIST("multvsmultTrkPVAft"), multNTracksPV, multTrk);
      histos.fill(HIST("multTrkPVvsMultT0CAft"), multT0C, multNTracksPV);
      histos.fill(HIST("multTrkPVvsmultV0AAft"), multV0A, multNTracksPV);
      histos.fill(HIST("multTrkPVvsmultT0AAft"), multT0A, multNTracksPV);
      histos.fill(HIST("multT0CvsmultT0AAft"), multT0A, multT0C);
      histos.fill(HIST("multV0AvsmultT0AAft"), multT0A, multV0A);
      histos.fill(HIST("multV0AvsmultT0CAft"), multT0C, multV0A);
      

      
    Int_t multITS = 0, multITSn = 0, multITSTPC = 0, multTPC = 0, multTRD = 0, multTOF = 0;
    
      for (auto& track : tracks) {
          
          if (track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF())
              multITS++;
          
          if (track.hasITS())
              multITSn++;
          
          if (track.hasTPC())
              multTPC++;
          
          if (track.hasITS() && track.hasTPC() && !track.hasTRD() && !track.hasTOF())
              multITSTPC++;
          
          if (track.hasTRD())
              multTRD++;
          
          if (track.hasTOF())
              multTOF++;
          
          
          Double_t trackpt = track.pt();
          Double_t tracketa = track.eta();
          Double_t trackdcaz = track.dcaZ();
          Double_t trackdcaxy = track.dcaXY();
          Double_t trackphi = track.phi();
          Double_t chargeQ = track.sign();
          Double_t nclsITS =  track.itsNCls();
          Double_t nclsTPC = track.tpcNClsFound();
          Double_t ncrsRowTPC = track.tpcNClsCrossedRows();
          
          histos.fill(HIST("QAEtaPhi"), tracketa, trackphi, chargeQ);
          histos.fill(HIST("QADCAz"), trackdcaz, trackpt, chargeQ);
          histos.fill(HIST("QADCAxy"), trackdcaxy, trackpt, chargeQ);
          
          histos.fill(HIST("QANclsITS"), nclsITS, trackpt, chargeQ);
          histos.fill(HIST("QANclsTPC"), nclsTPC, trackpt, chargeQ);
          histos.fill(HIST("QANcrRowTPC"), ncrsRowTPC, trackpt, chargeQ);
          histos.fill(HIST("QAPtEta"), tracketa, trackpt, chargeQ);
          
          
          if (TMath::Abs(tracketa) >= etaCut || nclsTPC < noClus || trackpt < minPt || trackpt >= maxPt || TMath::Abs(trackdcaz) >= dcazCut || TMath::Abs(trackdcaxy) >= dcaxyCut)
              continue;
          
          histos.fill(HIST("QAEtaPhiAft"), tracketa, trackphi, chargeQ);
          histos.fill(HIST("QADCAzAft"), trackdcaz, trackpt, chargeQ);
          histos.fill(HIST("QADCAxyAft"), trackdcaxy, trackpt, chargeQ);
          
          histos.fill(HIST("QANclsITSAft"), nclsITS, trackpt, chargeQ);
          histos.fill(HIST("QANclsTPCAft"), nclsTPC, trackpt, chargeQ);
          histos.fill(HIST("QANcrRowTPCAft"), ncrsRowTPC, trackpt, chargeQ);
          histos.fill(HIST("QAPtEtaAft"), tracketa, trackpt, chargeQ);
                    
      }
    
    
    histos.fill(HIST("multITSvsMultITSTPCAft"), multITSTPC, multITS);
    histos.fill(HIST("multITSvsMultITSTPCNAft"), multITSTPC, multITSn);
    histos.fill(HIST("multTPCvsMultITSTPCAft"), multITSTPC, multTPC);
    histos.fill(HIST("multTRDvsMultITSTPCAft"), multITSTPC, multTRD);
    histos.fill(HIST("multTOFvsMultITSTPCAft"), multITSTPC, multTOF);
    
  }
    
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qa3>(cfgc),
  };
}
