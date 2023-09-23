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
///        o2-analysis-timestamp --aod-file AO2D.root -b | o2-analysis-event-selection -b | o2-analysis-multiplicity-table -b | o2-analysis-centrality-table -b | o2-analysis-trackextension -b | o2-analysis-trackselection -b | o2-analysistutorial-flow-analysis -b
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
#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPObject.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct flow_base {

  using BCsWithRun2Infos = soa::Join<aod::BCs, aod::Run2BCInfos, o2::aod::Timestamps>;
  using Colls_EvSels_Mults_Cents = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentRun2V0Ms, aod::CentRun2CL0s, aod::CentRun2CL1s>;
  using FilteredCollisions = soa::Filtered<Colls_EvSels_Mults_Cents>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;
  using FilteredTracks = soa::Filtered<TrackCandidates>;

  Configurable<int> eventSelection{"eventSelection", 1, "event selection"};
  Configurable<bool> phiCut{"phiCut", false, "activate phi cut"};
  Configurable<bool> crsRowsFrcShCls{"crsRowsFrcShCls", false, "crsRowsFrcShCl"};
  Configurable<float> vtxCut{"vtxCut", 10.0, "Z vertex cut"};
  Configurable<float> etaCut{"etaCut", 0.8, "Eta cut"};
  Configurable<float> etaGap{"etaGap", 0.5, "Eta gap"};
  Configurable<int> noClus{"noClus", 70, "Number of clusters"};
  Configurable<int> nHarm{"nHarm", 2, "Number of harmonics"};
  Configurable<float> minPt{"minPt", 0.2, "Minimum pt"};
  Configurable<float> maxPt{"maxPt", 20.0, "Maximum pt"};

  TF1* fPhiCutLow = nullptr;
  TF1* fPhiCutHigh = nullptr;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

    
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Filter collisionFilter = (aod::collision::flags & (uint16_t)aod::collision::CollisionFlagsRun2::Run2VertexerTracks) == (uint16_t)aod::collision::CollisionFlagsRun2::Run2VertexerTracks;
  Filter trackFilter = ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));

  void fillAPt(double trackpt, double cent, double vn, double sinHarm, double cosHarm)
  {
    histos.fill(HIST("VnAPt"), trackpt, cent, vn);
    histos.fill(HIST("SinnAPt"), trackpt, cent, sinHarm);
    histos.fill(HIST("CosnAPt"), trackpt, cent, cosHarm);
  }

  void fillCPt(double trackpt, double cent, double vn, double sinHarm, double cosHarm)
  {
    histos.fill(HIST("VnCPt"), trackpt, cent, vn);
    histos.fill(HIST("SinnCPt"), trackpt, cent, sinHarm);
    histos.fill(HIST("CosnCPt"), trackpt, cent, cosHarm);
  }



  void init(InitContext&)
  {
    AxisSpec axisVtxcounts{2, -0.5f, 1.5f, "Vtx info (0=no, 1=yes)"};
    AxisSpec axisZvert{120, -30.f, 30.f, "Vtx z (cm)"};
    AxisSpec axisCent{100, 0.f, 100.f, "centrality V0M"};
    AxisSpec axisCentCL0{100, 0.f, 100.f, "centrality CL0"};
    AxisSpec axisCentCL1{100, 0.f, 100.f, "centrality CL1"};
    AxisSpec axisMult{1000, -0.5f, 3999.5f, "multiplicity"};
    AxisSpec axisTracklets{1000, -0.5f, 6999.5f, "SPD N_{tracklets}"};
    AxisSpec axisClusters{1000, -0.5f, 24999.5f, "SPD N_{clusters}"};
    AxisSpec axismultV0on{1000, 0, 50000, "multV0on"};
    AxisSpec axismultV0of{1000, 0, 50000, "multV0of"};
    AxisSpec axisCentBins{{0, 5., 10., 20., 30., 40., 50., 60., 70., 80.}, "centrality percentile"};
    AxisSpec axisPtBins{{0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10., 13., 16., 20.}, "p_{T} (GeV/c)"};

    histos.add("vtx", "Vtx info (0=no, 1=yes); Vtx; Counts", kTH1I, {axisVtxcounts});
    histos.add("vtxCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", kTH1F, {axisZvert});
    histos.add("multvsCent", "centrality vs multiplicity", kTH2F, {axisCent, axisMult});
    histos.add("cenCL0vsV0M", "centrality V0M vs centrality CL0", kTH2F, {axisCent, axisCentCL0});
    histos.add("cenCL1vsV0M", "centrality V0M vs centrality CL1", kTH2F, {axisCent, axisCentCL1});
    histos.add("cenCL1vsCL0", "centrality CL1 vs centrality CL0", kTH2F, {axisCentCL1, axisCentCL0});
    histos.add("SPclsvsSPDtrks", "SPD N_{tracklets} vs SPD N_{clusters}", kTH2I, {axisTracklets, axisClusters});
    histos.add("multV0onvsMultV0of", "V0 offline vs V0 online", kTH2F, {axismultV0of, axismultV0on});
    histos.add("res", "centrality percentile vs Resolution", kTProfile, {axisCentBins});
    histos.add("QxnA", "centrality percentile vs #LT Q_{x}^{nA} #GT", kTProfile, {axisCentBins});
    histos.add("QxnC", "centrality percentile vs #LT Q_{x}^{nC} #GT", kTProfile, {axisCentBins});
    histos.add("QynA", "centrality percentile vs #LT Q_{y}^{nA} #GT", kTProfile, {axisCentBins});
    histos.add("QynC", "centrality percentile vs #LT Q_{x}^{nC} #GT", kTProfile, {axisCentBins});

    histos.add("VnAPt", "v_{n} A", kTProfile2D, {{axisPtBins}, {axisCentBins}});
    histos.add("VnCPt", "v_{n} C", kTProfile2D, {{axisPtBins}, {axisCentBins}});
    histos.add("SinnAPt", "#LT sin(n*#phi) #GT A", kTProfile2D, {{axisPtBins}, {axisCentBins}});
    histos.add("SinnCPt", "#LT sin(n*#phi) #GT C", kTProfile2D, {{axisPtBins}, {axisCentBins}});
    histos.add("CosnAPt", "#LT cos(n*#phi) #GT A", kTProfile2D, {{axisPtBins}, {axisCentBins}});
    histos.add("CosnCPt", "#LT cos(n*#phi) #GT C", kTProfile2D, {{axisPtBins}, {axisCentBins}});



    fPhiCutLow = new TF1("fPhiCutLow", "0.1/x/x+pi/18.0-0.025", 0, 100);
    fPhiCutHigh = new TF1("fPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 100);

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now); // TODO must become global parameter from the train creation time
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

  void process(FilteredCollisions::iterator const& collision, BCsWithRun2Infos const& bcs, FilteredTracks const& tracks)
  {

    if ((eventSelection == 1) && (!collision.alias()[kINT7] || !collision.sel7())) {
      // LOGF(info, "Collision index : %d skipped not kINT7", collision.index());
      return;
    }

    float zvtx = -999;
      if (collision.numContrib() > 1) {
          
          zvtx = collision.posZ();
          
          float zRes = TMath::Sqrt(collision.covZZ());
          bool vertexerZ = collision.flags() == aod::collision::Run2VertexerZ;
          if (vertexerZ && zRes > 0.25 && collision.numContrib() < 20)
              zvtx = -999;
                    
      }
      
      if (zvtx < -990)
          histos.fill(HIST("vtx"), 0);
      else
          histos.fill(HIST("vtx"), 1);
      
      if (TMath::Abs(zvtx) > vtxCut)
          return;
      
      auto v0Centr = collision.centRun2V0M();
      auto cl1Centr = collision.centRun2CL1();
      auto cl0Centr = collision.centRun2CL0();
      
      if (v0Centr >= 80. || v0Centr < 0)
          return;
      
      // cannot use vertex quality comparing SPD and Trk vertices as below
      // (Should we use ncontrib, chi2)?
      /*
       if (TMath::Abs(dz)>0.5 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)
       return; // bad vertexing
       */
      // float errTot = collision.covZZ();
      
      auto bc = collision.bc_as<BCsWithRun2Infos>();
      auto field = getMagneticField(bc.timestamp());
      
      auto nITSClsLy0 = bc.spdClustersL0();
      auto nITSClsLy1 = bc.spdClustersL1();
      auto nITSCls = nITSClsLy0 + nITSClsLy1;
      
      auto nITSTrkls = collision.multTracklets();
      
      auto multV0a = collision.multFV0A();
      auto multV0c = collision.multFV0C();
      auto multV0Tot = multV0a + multV0c;
      auto multV0aOn = bc.v0TriggerChargeA();
      auto multV0cOn = bc.v0TriggerChargeC();
      auto multV0On = multV0aOn + multV0cOn;
      
      
      histos.fill(HIST("vtxCuts"), zvtx);
      histos.fill(HIST("cenCL0vsV0M"), v0Centr, cl0Centr);
      histos.fill(HIST("cenCL1vsV0M"), v0Centr, cl1Centr);
      histos.fill(HIST("cenCL1vsCL0"), cl1Centr, cl0Centr);
      histos.fill(HIST("SPclsvsSPDtrks"), nITSTrkls, nITSCls);
      histos.fill(HIST("multV0onvsMultV0of"), multV0Tot, multV0On);
      
      
      // process the tracks of a given collision
      Double_t QxnGapA = 0., QynGapA = 0.;
      Double_t QxnGapC = 0., QynGapC = 0.;
      
      Int_t multGapA = 0, multGapC = 0;
      
      // Tracks are already filtered with GlobalTrack || GlobalTrackSDD
      Int_t multTrk = tracks.size();
      
      for (auto& track : tracks) {
          
          Double_t trackpt = track.pt();
          Double_t tracketa = track.eta();
          
          if (TMath::Abs(tracketa) >= etaCut ||
              track.tpcNClsFound() < noClus ||
              trackpt < minPt || trackpt >= maxPt)
              continue;
          
          Double_t sinHarm = TMath::Sin(nHarm * track.phi());
          Double_t cosHarm = TMath::Cos(nHarm * track.phi());
          
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
      
      histos.fill(HIST("multvsCent"), v0Centr, multTrk);
      
      if (multGapA > 0 && multGapC > 0) {
          Double_t resGap = (QxnGapA * QxnGapC + QynGapA * QynGapC) / (multGapA * multGapC);
          histos.fill(HIST("res"), v0Centr, resGap);
          
          histos.fill(HIST("QxnA"), v0Centr, QxnGapA / multGapA);
          histos.fill(HIST("QxnC"), v0Centr, QxnGapC / multGapC);
          
          histos.fill(HIST("QynA"), v0Centr, QynGapA / multGapA);
          histos.fill(HIST("QynC"), v0Centr, QynGapC / multGapC);
      }
      
      for (auto& track : tracks) {
          
          Double_t trackpt = track.pt();
          Double_t tracketa = track.eta();
          
          if (TMath::Abs(tracketa) >= etaCut ||
              track.tpcNClsFound() < noClus ||
              trackpt < minPt || trackpt >= maxPt)
              continue;
          
          if (phiCut) {
              Double_t phimod = track.phi();
              if (field < 0) // for negative polarity field
                  phimod = TMath::TwoPi() - phimod;
              if (track.sign() < 0) // for negative charge
                  phimod = TMath::TwoPi() - phimod;
              if (phimod < 0)
                  LOGF(warning, "phi < 0: %g", phimod);
              
              phimod += TMath::Pi() / 18.0; // to center gap in the middle
              phimod = fmod(phimod, TMath::Pi() / 9.0);
              if (phimod < fPhiCutHigh->Eval(trackpt) && phimod > fPhiCutLow->Eval(trackpt))
                  continue; // reject track
          }
          
          if (crsRowsFrcShCls) {
              Float_t nrowscr = track.tpcNClsCrossedRows();
              if (nrowscr < 120)
                  continue;
              
              Float_t clsFind = track.tpcNClsFindable();
              if (clsFind <= 0)
                  continue;
              
              if (track.tpcCrossedRowsOverFindableCls() < 0.9)
                  continue;
          }
          
          Double_t sinHarmn = TMath::Sin(nHarm * track.phi());
          Double_t cosHarmn = TMath::Cos(nHarm * track.phi());
          
          Double_t harmGapC = cosHarmn * QxnGapC + sinHarmn * QynGapC;
          Double_t harmGapA = cosHarmn * QxnGapA + sinHarmn * QynGapA;
          
          if (tracketa > etaGap && multGapA > 0) {
              Double_t vnC = harmGapA / multGapA;
              fillCPt(trackpt, v0Centr, vnC, sinHarmn, cosHarmn);
          }
          
          if (tracketa < -etaGap && multGapC > 0) {
              Double_t vnA = harmGapC / multGapC;
              fillAPt(trackpt, v0Centr, vnA, sinHarmn, cosHarmn);
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
