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
/// \brief Flow analysis :: converter 2tree
///        Run as:
///        o2-analysis-timestamp --aod-file AO2D.root -b | o2-analysis-event-selection -b | o2-analysis-multiplicity-table -b | o2-analysis-centrality-table -b | o2-analysis-trackextension -b | o2-analysis-trackselection -b | o2-analysis-pid-tpc-full -b | o2-analysis-pid-tof-full -b | o2-analysis-pid-tof-beta -b | o2-analysistutorial-v0zdc-treecalib -b
/// \author
/// \since

#include <memory>
#include <TreeCalibV0ZDC.h>
#include <Framework/runDataProcessing.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/AnalysisDataModel.h>
#include <Common/DataModel/EventSelection.h>
#include <Common/CCDB/TriggerAliases.h>
#include <Common/DataModel/Centrality.h>
#include <Common/DataModel/Multiplicity.h>
#include <Common/DataModel/TrackSelectionTables.h>
#include <Common/Core/PID/PIDResponse.h>
#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPObject.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;


// define task content/structure
struct tree_calib {

using BCsWithRun2Infos = soa::Join<aod::BCs, aod::Run2BCInfos, o2::aod::Timestamps>;
using Colls_EvSels_Mults_Cents = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentRun2V0Ms, aod::CentRun2CL0s, aod::CentRun2CL1s>;
using FilteredCollisions = soa::Filtered<Colls_EvSels_Mults_Cents>;
using TracksPID = soa::Join<aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, TracksPID>;
using FilteredTracks = soa::Filtered<TrackCandidates>;

Filter collisionFilter = (aod::collision::flags & (uint16_t)aod::collision::CollisionFlagsRun2::Run2VertexerTracks) == (uint16_t)aod::collision::CollisionFlagsRun2::Run2VertexerTracks;
Filter trackFilter = ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));

Service<o2::ccdb::BasicCCDBManager> ccdb;
Configurable<std::string> outFile {"outFile", std::string("calib.root"), "Name of the output file"};
Configurable<std::string> treeName {"treeName", std::string("tree"), "Name of the output TTree"};
Configurable<float> vtxCut {"vtxCut", 10.0, "Z vertex cut"};

// Output settings
std::unique_ptr<TFile> output_tfile( TFile::Open(outFile, "UPDATE") );
auto fTree = std::make_unique<TTree>(treeName, "Event data");
Event* fEvent = nullptr;
TClonesArray* fVzeroArray = nullptr;
static const int fVzeroArray_size = 64;

void init(InitContext&) {
    if (!output_tfile) {
        LOGF(fatal, "Could not open the output TFile");
        }

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now); // TODO must become global parameter from the train creation time

    fEvent = new Event();
    fVzeroArray = new TClonesArray("VZero", fVzeroArray_size);
    fTree->Branch("event", &fEvent);
    fTree->Branch("vzero", "TClonesArray", &fVzeroArray);
    }

void process(FilteredCollisions::iterator const& collision, BCsWithRun2Infos const& bcs, FilteredTracks const& tracks) {

    float zvtx = -999;
    if (collision.numContrib() > 1) {
        float zRes = TMath::Sqrt(collision.covZZ());
        bool vertexerZ = collision.flags() == aod::collision::Run2VertexerZ;
        if (vertexerZ && zRes < 0.25 && collision.numContrib() > 20) { zvtx = collision.posZ(); }
        }

    if ( TMath::Abs(zvtx) > vtxCut ) { return; }

    auto v0Centr = collision.centRun2V0M();
    if (v0Centr > 90) { return; }

    auto cl1Centr = collision.centRun2CL1();
    auto cl0Centr = collision.centRun2CL0();

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

    // Number of total TPC clusters and total ITS for layers 2-6 not available.
    // but the task only fills some histos on these
    /*
    int tpcClsTot = aod->GetNumberOfTPCClusters();

    //clusters SDD+SSD
    AliVMultiplicity* mult = aod->GetMultiplicity();
    int nCluSDDSSD=0;
    for(int iLay = 2; iLay < 6; iLay++)
        nCluSDDSSD += mult->GetNumberOfITSClusters(iLay);
    */

    //new vertex selection
    /*
    const AliAODVertex* vtTrc = aodEvent->GetPrimaryVertex();
    const AliAODVertex* vtSPD = aodEvent->GetPrimaryVertexSPD();

    if (vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1) { return; } // one of vertices is missing

    double covTrc[6], covSPD[6];
    vtTrc->GetCovarianceMatrix(covTrc);
    vtSPD->GetCovarianceMatrix(covSPD);
    
    double dz = vtTrc->GetZ() - vtSPD->GetZ();
    
    double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
    double errTrc = TMath::Sqrt(covTrc[5]);
    double nsigTot = dz/errTot;
    double nsigTrc = dz/errTrc;
    
    if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20) { return; } // bad vertexing
    
    TString vtxTyp = vtSPD->GetTitle();
    Double_t zRes = TMath::Sqrt(covSPD[5]);
    if ((vtxTyp.Contains("vertexer:Z")) && (zRes>0.25) && (vtSPD->GetNContributors() < 20)) { return; } // bad vertexing

    // if (((AliAODHeader*)aodEvent->GetHeader())->GetRefMultiplicityComb08() < 0) { return; }
    
    //new function for 2015 to remove incomplete events
    // if (aodEvent->IsIncompleteDAQ()) { return; }
    
*/    


/*    
    const Int_t nTracks = aodEvent->GetNumberOfTracks();
    Int_t multEsd = ((AliAODHeader*)aodEvent->GetHeader())->GetNumberOfESDTracks();
    Int_t multEsdBefpu = aodEvent->GetNTPCTrackBeforeClean();
    

    Int_t multTrkBefC = 0;
    Int_t multTrkTOFBefC = 0;
    Int_t multTPC = 0;
    Int_t multTPCout = 0;
    Int_t multITS = 0;
    
    for (Int_t it = 0; it < nTracks; it++) {
        
        AliAODTrack* aodTrk = (AliAODTrack*)aodEvent->GetTrack(it);
        
        if (!aodTrk){
            delete aodTrk;
            continue;
        }
        
        if (aodTrk->GetFlags()&AliESDtrack::kTPCout)
            multTPCout++;
        
        
        if ((aodTrk->GetFlags()&(AliESDtrack::kITSin | AliESDtrack::kTPCin))==AliESDtrack::kITSin)
            multITS++;
        
        
        if (aodTrk->TestFilterBit(32)){
            
            multTrkBefC++;
            
            if ( TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10. && aodTrk->GetTOFsignal() >= 12000. && aodTrk->GetTOFsignal() <= 25000.)
                multTrkTOFBefC++;
      
            }
        
        if (aodTrk->TestFilterBit(128)) multTPC++;
        
    }
*/    

/*

    Int_t tpcClsTot = aodEvent->GetNumberOfTPCClusters();
    UInt_t timeSt = aodEvent->GetTimeStamp();
    Int_t run = aodEvent->GetRunNumber();
    ULong64_t eventId = -1;
    if(aodEvent->GetHeader()) eventId = GetEventIdAsLong((AliAODHeader*)aodEvent->GetHeader());
    

    AliAODTracklets* aodTrkl = (AliAODTracklets*)aodEvent->GetTracklets();
    Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();
    
    Int_t nITSClsLy0 = aodEvent->GetNumberOfITSClusters(0);
    Int_t nITSClsLy1 = aodEvent->GetNumberOfITSClusters(1);
    Int_t nITSCls = nITSClsLy0 + nITSClsLy1;
    
    //clusters SDD+SSD
    AliVMultiplicity* mult = aodEvent->GetMultiplicity();
    Int_t nCluSDDSSD=0;
    for(Int_t iLay = 2; iLay < 6; iLay++) nCluSDDSSD += mult->GetNumberOfITSClusters(iLay);
    
    //V0 info
    AliAODVZERO* aodV0 = (AliAODVZERO*)aodEvent->GetVZEROData();
    Float_t multV0a = aodV0->GetMTotV0A();
    Float_t multV0c = aodV0->GetMTotV0C();
    Float_t multV0 = multV0a + multV0c;
    UShort_t multV0aOn = aodV0->GetTriggerChargeA();
    UShort_t multV0cOn = aodV0->GetTriggerChargeC();
    UShort_t multV0On = multV0aOn + multV0cOn;
    
    Int_t nAddVzero = 0;
  
    for (Int_t iv0 = 0; iv0 < 64; iv0++) {
    
        Float_t multv0 = aodV0->GetMultiplicity(iv0);
        Float_t multv0eq = aodEvent->GetVZEROEqMultiplicity(iv0);
    
        Vzero* v0 = new((*fVzeroArray)[nAddVzero]) Vzero();
        nAddVzero++;
    
        v0->cIdV0     = iv0;
        v0->cMultV0   = multv0;
        v0->cMultV0Eq = multv0eq;
        
    }
*/    

/*
    AliAODZDC* aodZDC = aodEvent->GetZDCData();
    
    Double_t centroidZNC[2] = {999., 999.};
    Double_t centroidZNA[2] = {999., 999.};
    Double_t beamEn = 2511.;
    
    Double_t towZNC[5] = {999., 999., 999., 999., 999.};
    Double_t towZNA[5] = {999., 999., 999., 999., 999.};
    
    if (aodZDC->GetZNCentroidInPbPb(beamEn, centroidZNC, centroidZNA)){
        for (Int_t i = 0; i < 5; i++) {
            towZNC[i] = aodZDC->GetZNCTowerEnergy()[i];
            towZNA[i] = aodZDC->GetZNATowerEnergy()[i];
        }
    }
*/

    //Global variables
    fEvent->eventid      = eventId;
    fEvent->run          = run;
    fEvent->v0Cent       = v0Centr;
    fEvent->cl0Cent      = cl0Centr;
    fEvent->cl1Cent      = cl1Centr;
    fEvent->zvtx         = zVtx;
    fEvent->xvtx         = vtTrc->GetX();
    fEvent->yvtx         = vtTrc->GetY();
    fEvent->itsCls       = nITSCls;
    fEvent->multV0off    = multV0;
    fEvent->multV0on     = multV0On;
    fEvent->multTrklAll  = nITSTrkls;
    fEvent->centroidZNCX = centroidZNC[0];
    fEvent->centroidZNCY = centroidZNC[1];
    fEvent->centroidZNAX = centroidZNA[0];
    fEvent->centroidZNAY = centroidZNA[1];
    fEvent->towEnZNC0 = towZNC[0];
    fEvent->towEnZNC1 = towZNC[1];
    fEvent->towEnZNC2 = towZNC[2];
    fEvent->towEnZNC3 = towZNC[3];
    fEvent->towEnZNC4 = towZNC[4];
    fEvent->towEnZNA0 = towZNA[0];
    fEvent->towEnZNA1 = towZNA[1];
    fEvent->towEnZNA2 = towZNA[2];
    fEvent->towEnZNA3 = towZNA[3];
    fEvent->towEnZNA4 = towZNA[4];
    fEvent->timeStamp = timeSt;
    fEvent->totTPCCls = tpcClsTot;
    fEvent->multESD   = multEsd;
    fEvent->multTPCout = multTPCout;
    fEvent->multITSTPCTOF = multTrkTOFBefC;
    fEvent->multITSTPC = multTrkBefC;
    fEvent->multITS = multITS;
    fEvent->multTPC = multTPC;
    fEvent->multESDBefPU = multEsdBefpu;
    fEvent->itsClsSDDSSD = nCluSDDSSD;
    fEvent->multAOD = nTracks;


    for (auto& track : tracks) {
      Double_t trackpt = track.pt();
      Double_t tracketa = track.eta();
      }


    fTree->Fill();
    fVzeroArray->Clear();

    }

int getMagneticField(uint64_t timestamp) {
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

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
    return WorkflowSpec { adaptAnalysisTask<tree_calib>(cfgc), };
    }
