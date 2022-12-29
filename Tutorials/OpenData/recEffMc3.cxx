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
///        o2-analysis-timestamp --aod-file AO2D.root -b  | o2-analysis-track-propagation -b | o2-analysis-trackselection -b | o2-analysistutorial-flow-analysis3 -b
/// \author
/// \since

#include <TH1.h>
#include <TH3.h>
#include <TH2.h>
#include <TPDGCode.h>
#include <TMath.h>

#include <Framework/runDataProcessing.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/AnalysisDataModel.h>
#include <Common/DataModel/TrackSelectionTables.h>
#include <Common/DataModel/Multiplicity.h>
#include <Common/DataModel/EventSelection.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;


struct runEffMc {

    using Colls = soa::Join<aod::Collisions, aod::Mults, aod::McCollisionLabels>;
    using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>;
    using FilteredTracks = soa::Filtered<TrackCandidates>;
    Filter trackFilter = (requireGlobalTrackInFilter());

    
    Configurable<bool> crsRowsFrcShCls{"crsRowsFrcShCls", false, "crsRowsFrcShCl"};
    Configurable<float> vtxCut{"vtxCut", 10.0, "Z vertex cut"};
    Configurable<float> etaCut{"etaCut", 0.8, "Eta cut"};
    Configurable<int> noClus{"noClus", 70, "Number of clusters"};
    Configurable<float> minPt{"minPt", 0.2, "Minimum pt"};
    Configurable<float> maxPt{"maxPt", 20.0, "Maximum pt"};


    HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

    
    void init(InitContext&)
    {
      
        AxisSpec axisVtxcounts{2, -0.5f, 1.5f, "Vtx info (0=no, 1=yes)"};
        AxisSpec axisZvert{120, -30.f, 30.f, "Vtx z (cm)"};
        AxisSpec axisPtBins{{0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10., 13., 16., 20.}, "p_{T} (GeV/c)"};
        AxisSpec axisEta{16, -0.8f, 0.8f, "#eta"};
        AxisSpec axisPid{7, 0.f, 7.f, "pid (pdg code)"};
        AxisSpec axisEtaQA{72, -0.9f, 0.9f, "#eta"};
        AxisSpec axisPhiQA{144, 0.f, TMath::TwoPi(), "varphi"};
        AxisSpec axisDCAz{500, -10.f, 10.f, "DCA_{z}"};
        AxisSpec axisDCAxy{500, -10.f, 10.f, "DCA_{xy}"};
        AxisSpec axisMultFw{1000, 0.f, 200000.f, "multiplicity_fw"};
        AxisSpec axisMult{1000, 0.f, 4000.f, "multiplicity"};
        AxisSpec axisDif{250, -2.5f, 2.5f, "truth - rec"};
        AxisSpec axisPhiDif{100, -TMath::Pi(), TMath::Pi(), "truth - rec"};


        histos.add("rec/vtx", "Vtx info (0=no, 1=yes); Vtx; Counts", kTH1I, {axisVtxcounts});
      
        histos.add("rec/vtxCutsBef", "Vtx distribution; Vtx z [cm]; Counts", kTH1F, {axisZvert});
        histos.add("rec/multvsMultT0CBef", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
        histos.add("rec/multvsmultV0ABef", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
        histos.add("rec/multvsmultT0ABef", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
        histos.add("rec/multvsmultTrkPVBef", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
        histos.add("rec/multTrkPVvsMultT0CBef", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
        histos.add("rec/multTrkPVvsmultV0ABef", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
        histos.add("rec/multTrkPVvsmultT0ABef", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
        histos.add("rec/multT0CvsmultT0ABef", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
        histos.add("rec/multV0AvsmultT0ABef", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
        histos.add("rec/multV0AvsmultT0CBef", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});

        histos.add("rec/vtxCutsAft", "Vtx distribution; Vtx z [cm]; Counts", kTH1F, {axisZvert});
        histos.add("rec/multvsMultT0CAft", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
        histos.add("rec/multvsmultV0AAft", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
        histos.add("rec/multvsmultT0AAft", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
        histos.add("rec/multvsmultTrkPVAft", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
        histos.add("rec/multTrkPVvsMultT0CAft", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
        histos.add("rec/multTrkPVvsmultV0AAft", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
        histos.add("rec/multTrkPVvsmultT0AAft", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
        histos.add("rec/multT0CvsmultT0AAft", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
        histos.add("rec/multV0AvsmultT0AAft", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
        histos.add("rec/multV0AvsmultT0CAft", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
      
        histos.add("rec/ptEtaPhiQAPrim", "p_{T} vs #eta vs #varphi", kTH3F, {{axisPtBins}, {axisEtaQA}, {axisPhiQA}});
        histos.add("rec/DcazQAPrim", "DCAz", kTH1F, {axisDCAz});
        histos.add("rec/DcaxyQAPrim", "DCAxy", kTH1F, {axisDCAxy});
        histos.add("rec/ptEtaPidPrim", "p_{T} vs #eta vs pdg code", kTH3F, {{axisPtBins}, {axisEta}, {axisPid}});
        histos.add("rec/ptEtaPhiDifPrim", "p_{T} vs #eta vs #varphi (truth - rec)", kTH3F, {{axisDif}, {axisDif}, {axisPhiDif}});
      
        histos.add("rec/ptEtaPhiQASec", "p_{T} vs #eta vs #varphi", kTH3F, {{axisPtBins}, {axisEtaQA}, {axisPhiQA}});
        histos.add("rec/DcazQASec", "DCAz", kTH1F, {axisDCAz});
        histos.add("rec/DcaxyQASec", "DCAxy", kTH1F, {axisDCAxy});
        histos.add("rec/ptEtaPidSec", "p_{T} vs #eta vs pdg code", kTH3F, {{axisPtBins}, {axisEta}, {axisPid}});
  
        
        histos.add("truth/vtxCorTruthRecBef", "Vtx distribution;  Vtx z (cm); Vtx z (cm)", kTH2F, {{axisZvert}, {axisZvert}});
        histos.add("truth/vtxCorTruthRecAft", "Vtx distribution;  Vtx z (cm); Vtx z (cm)", kTH2F, {{axisZvert}, {axisZvert}});
        
        histos.add("truth/vtxCutsBefMC", "Vtx distribution; Vtx z [cm]; Counts", kTH1F, {axisZvert});
        histos.add("truth/vtxCutsAftMC", "Vtx distribution; Vtx z [cm]; Counts", kTH1F, {axisZvert});
        
        histos.add("truth/ptEtaPhiQAMC", "#eta vs #varphi", kTH3F, {{axisPtBins}, {axisEtaQA}, {axisPhiQA}});
        histos.add("truth/ptEtaPidMC", "p_{T} vs #eta vs pdg code", kTH3F, {{axisPtBins}, {axisEta}, {axisPid}});
      
    }

    int getPidCode(int pdgCode)
    {
        int pidCode = 0;
        switch (pdgCode) {
            case 211:
                pidCode = 1; // pion
                break;
            case 321:
                pidCode = 2; // kaon
                break;
            case 2212:
                pidCode = 3; // proton
                break;
            case 11:
                pidCode = 4; // electron
                break;
            case 13:
                pidCode = 5; // muon
                break;
            default:
                pidCode = 6;  // something else?
        };
        return pidCode;
    }
 
    
    
    void process(Colls::iterator const& collision, FilteredTracks const& tracks, aod::McCollisions const& mcCollisions, aod::McParticles& mcParticles)
    {

        float zvtx = collision.posZ();
        /*
        float zRes = TMath::Sqrt(collision.covZZ());
        if ((collision.numContrib() < 2) || (zRes > 0.25 && collision.numContrib() < 20))
            zvtx = -999;
         */
        if (collision.numContrib() < 2)
            zvtx = -999;
        
        if (zvtx < -990)
            histos.fill(HIST("rec/vtx"), 0);
        else
            histos.fill(HIST("rec/vtx"), 1);
      
      
        auto multV0A = collision.multFV0A();
        auto multT0A = collision.multFT0A();
        auto multT0C = collision.multFT0C();
        auto multNTracksPV = collision.multNTracksPV();
        auto multTrk = tracks.size();
      
        histos.fill(HIST("rec/vtxCutsBef"), zvtx);
        histos.fill(HIST("rec/multvsMultT0CBef"), multT0C, multTrk);
        histos.fill(HIST("rec/multvsmultV0ABef"), multV0A, multTrk);
        histos.fill(HIST("rec/multvsmultT0ABef"), multT0A, multTrk);
        histos.fill(HIST("rec/multvsmultTrkPVBef"), multNTracksPV, multTrk);
        histos.fill(HIST("rec/multTrkPVvsMultT0CBef"), multT0C, multNTracksPV);
        histos.fill(HIST("rec/multTrkPVvsmultV0ABef"), multV0A, multNTracksPV);
        histos.fill(HIST("rec/multTrkPVvsmultT0ABef"), multT0A, multNTracksPV);
        histos.fill(HIST("rec/multT0CvsmultT0ABef"), multT0A, multT0C);
        histos.fill(HIST("rec/multV0AvsmultT0ABef"), multT0A, multV0A);
        histos.fill(HIST("rec/multV0AvsmultT0CBef"), multT0C, multV0A);


        if (TMath::Abs(zvtx) <= vtxCut) {
            
            histos.fill(HIST("rec/vtxCutsAft"), zvtx);
            histos.fill(HIST("rec/multvsMultT0CAft"), multT0C, multTrk);
            histos.fill(HIST("rec/multvsmultV0AAft"), multV0A, multTrk);
            histos.fill(HIST("rec/multvsmultT0AAft"), multT0A, multTrk);
            histos.fill(HIST("rec/multvsmultTrkPVAft"), multNTracksPV, multTrk);
            histos.fill(HIST("rec/multTrkPVvsMultT0CAft"), multT0C, multNTracksPV);
            histos.fill(HIST("rec/multTrkPVvsmultV0AAft"), multV0A, multNTracksPV);
            histos.fill(HIST("rec/multTrkPVvsmultT0AAft"), multT0A, multNTracksPV);
            histos.fill(HIST("rec/multT0CvsmultT0AAft"), multT0A, multT0C);
            histos.fill(HIST("rec/multV0AvsmultT0AAft"), multT0A, multV0A);
            histos.fill(HIST("rec/multV0AvsmultT0CAft"), multT0C, multV0A);

            
            for (auto& track : tracks) {

                Double_t trackpt = track.pt();
                Double_t tracketa = track.eta();

                if (TMath::Abs(tracketa) >= etaCut || track.tpcNClsFound() < noClus || trackpt < minPt || trackpt >= maxPt)
                    continue;
                
                
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
                
                Double_t trackdcaz = track.dcaZ();
                Double_t trackdcaxy = track.dcaXY();
                Double_t trackphi = track.phi();
                
                int pdgCode = TMath::Abs(track.mcParticle().pdgCode());
                Int_t pidCode = getPidCode(pdgCode);
                                    
                if (track.mcParticle().isPhysicalPrimary()) {
                 
                    histos.fill(HIST("rec/ptEtaPhiQAPrim"), trackpt, tracketa, trackphi);
                    histos.fill(HIST("rec/DcazQAPrim"), trackdcaz);
                    histos.fill(HIST("rec/DcaxyQAPrim"), trackdcaxy);
                    histos.fill(HIST("rec/ptEtaPidPrim"), trackpt, tracketa, pidCode);
                      
                    Double_t parteta = track.mcParticle().eta();
                    Double_t partpt = track.mcParticle().pt();
                    Double_t partphi = track.mcParticle().phi();
                    
                    Double_t etaDif = parteta - tracketa;
                    Double_t ptDif = partpt - trackpt;
                    Double_t phiDif = partphi - trackphi;
                    if (phiDif > TMath::Pi())
                        phiDif -= 2.*TMath::Pi();
                    if (phiDif < -TMath::Pi())
                        phiDif += 2.*TMath::Pi();
                    histos.fill(HIST("rec/ptEtaPhiDifPrim"), ptDif, etaDif, phiDif);
                    
                } else {
                    
                    histos.fill(HIST("rec/ptEtaPhiQASec"), trackpt, tracketa, trackphi);
                    histos.fill(HIST("rec/DcazQASec"), trackdcaz);
                    histos.fill(HIST("rec/DcaxyQASec"), trackdcaxy);
                    histos.fill(HIST("rec/ptEtaPidSec"), trackpt, tracketa, pidCode);
                    
                }
                
            }
            
        }
        
        
        
        //truth
        float zvtxMC = collision.mcCollision().posZ();

        histos.fill(HIST("truth/vtxCorTruthRecBef"), zvtxMC, zvtx);
        histos.fill(HIST("truth/vtxCutsBefMC"), zvtxMC);
        
        if (TMath::Abs(zvtxMC) <= vtxCut) {
            
            histos.fill(HIST("truth/vtxCorTruthRecAft"), zvtxMC, zvtx);
            histos.fill(HIST("truth/vtxCutsAftMC"), zvtxMC);
            
            for (auto& mcPart : mcParticles) {
            
                Double_t ptPart = mcPart.pt();
                Double_t etaPart = mcPart.eta();

                if (TMath::Abs(etaPart) >= etaCut || ptPart < minPt || ptPart >= maxPt)
                    continue;
                
                if (!mcPart.isPhysicalPrimary())
                    continue;
                
                Double_t phiPart = mcPart.phi();
                int pdgCodeMC = TMath::Abs(mcPart.pdgCode());
                Int_t pidCodeMC = getPidCode(pdgCodeMC);
                                               
                histos.fill(HIST("truth/ptEtaPhiQAMC"), ptPart, etaPart, phiPart);
                histos.fill(HIST("truth/ptEtaPidMC"), ptPart, etaPart, pidCodeMC);

            }
                    
        }
        
    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
    return WorkflowSpec{
      adaptAnalysisTask<runEffMc>(cfgc),
    };
}
