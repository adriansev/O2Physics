#include <TObject.h>

class Vzero: public TObject {
public:
    Float_t    cMultV0;
    Float_t    cMultV0Eq;
    Int_t      cIdV0;
    Vzero(): cMultV0(-999.99), cMultV0Eq(-999.99), cIdV0(-1) {};
    ClassDef(Vzero, 1);
    };

class Event: public TObject {
public:
    Float_t   v0Cent;        // centrality V0
    Float_t   cl0Cent;       // centrality CL0
    Float_t   cl1Cent;       // centrality CL1
    Float_t   zvtx;          // rec vertex z tracks
    Float_t   xvtx;          // rec vertex x tracks
    Float_t   yvtx;          // rec vertex y tracks
    Float_t   centroidZNCX;        // centroid ZNC x
    Float_t   centroidZNCY;        // centroid ZNC y
    Float_t   centroidZNAX;        // centroid ZNA x
    Float_t   centroidZNAY;        // centroid ZNA y
    Float_t   towEnZNC0;           // ZNC tower energy channel 0
    Float_t   towEnZNC1;           // ZNC tower energy channel 1
    Float_t   towEnZNC2;           // ZNC tower energy channel 2
    Float_t   towEnZNC3;           // ZNC tower energy channel 3
    Float_t   towEnZNC4;           // ZNC tower energy channel 4
    Float_t   towEnZNA0;           // ZNA tower energy channel 0
    Float_t   towEnZNA1;           // ZNA tower energy channel 1
    Float_t   towEnZNA2;           // ZNA tower energy channel 2
    Float_t   towEnZNA3;           // ZNA tower energy channel 3
    Float_t   towEnZNA4;           // ZNA tower energy channel 4
    Float_t   multV0off;          // V0 multiplicity offline
    UInt_t    multV0on;           // V0 multiplicity online
    Int_t     multESD;             // multiplicity ESD
    Int_t     multTPCout;          // multiplicity TPC out
    Int_t     multITSTPCTOF;       // multiplicity ITS+TPC+TOF
    Int_t     multITSTPC;          // multiplicity ITS+TPC
    Int_t     multITS;             // multiplicity ITS
    Int_t     multTPC;             // multiplicity TPC
    Int_t     multESDBefPU;        // multiplicity ESD before PU removal
    Int_t     itsCls;              // its clusters layer 0+layer 1
    Int_t     itsClsSDDSSD;        // number of clusters from SDD+SSD
    Int_t     totTPCCls;           // number of TPC clusters per event
    Int_t     multAOD;             // multiplicity AOD
    Int_t     multTrklAll;         // mult tracklets all
    Int_t     run;                 // run number
    UInt_t    timeStamp;           // event time stamp
    ULong64_t eventid;       // unique event id

    Event(): v0Cent(-999.999), cl0Cent(-999.999), cl1Cent(-999.999), zvtx(-999.999), xvtx(-999.999), yvtx(-999.999),
    centroidZNCX(-999.999), centroidZNCY(-999.999), centroidZNAX(-999.999), centroidZNAY(-999.999),
    towEnZNC0(-999.999), towEnZNC1(-999.999), towEnZNC2(-999.999), towEnZNC3(-999.999), towEnZNC4(-999.999),
    towEnZNA0(-999.999), towEnZNA1(-999.999), towEnZNA2(-999.999), towEnZNA3(-999.999), towEnZNA4(-999.999),
    multV0off(-999.999),
    multV0on(0),  // unsigned int
    multESD(-999), multTPCout(-999), multITSTPCTOF(-999), multITSTPC(-999), multITS(-999), multTPC(-999), multESDBefPU(-999),
    itsCls(-999), itsClsSDDSSD(-999), totTPCCls(-999), multAOD(-999), multTrklAll(-999),
    run(-1), timeStamp(-1), eventid(-1) {};
    ClassDef(Event, 1);
    };

