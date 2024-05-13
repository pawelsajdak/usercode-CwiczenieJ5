#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//My contribution from 15/12/2022
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/MuonDetId/interface/DTSuperLayerId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTChamber.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCChamber.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTSuperLayer.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TEfficiency.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>





using namespace std;

typedef struct{
    int  event;
    int  track;
    int  pid;
    int  superlayer;
    int  station;
    int  wheel;
    int  sector;
    int  layer;
    double localX;
    double localY;
    double localZ;
    double globalX;
    double globalY;
    double globalZ;
    double pt;
    double phi;
    double r;
    double eta;
    double tof;
    double beta;
    double distanceVTX;
    double distanceIP;
    double nextLocalX;
    double nextLocalY;
    double nextLocalZ;
    double nextGlobalX;
    double nextGlobalY;
    double nextGlobalZ;
    double globalHitDistDiff;
    double betaNextHit;
    double diffTOFGeant;
    double diffTOFHitB;
    double diffTOFNextB;
    bool samewheel;
    bool samesl;
    bool samesect;
    bool samest;
  } SimHitData;

typedef struct{
    int  pdgId;
    int  event;
    double phi;
    double theta;
    double eta;
    double pT;
    double mass;
    double vx;
    double vy;
    double vz;
    double beta;
    double invbeta;
}RecCandidateData;  

typedef struct{
    int  pdgId;
    int  event;
    double phi;
    double theta;
    double eta;
    double pT;
    double mass;
    double beta;
    double invbeta;
}GenCandidateData; 

typedef struct{
    double pgen;
    double preco;
    double invbetagen;
    double invbetamanual;
    double invbetareco;
    double invbetarecoerror;
    double event;
    double deltaR;
    double genPt;
    double innerrecoPt;
    double outerrecoPt;
    double recoPt;
    double genEta;
    double recoEta;
    double recoPhi;
    double genPhi;
    int  stationMatches;
}InverseBetaData;

typedef struct{
    double pt;
    double eta;
    double p;
    double invbeta;
}TotalMCData;

typedef struct{
    double pt;
    double ptGen;
    double eta;
    double etaGen;
    double p;
    double pGen;
    double invbeta;
    double invbetaGen;
    double deltaR;
}GMTMuonData;

//object definition
class Projekt : public edm::EDAnalyzer {
    public:

        //constructor, function is called when new object is created
        explicit Projekt(const edm::ParameterSet& conf);

        //destructor, function is called when object is destroyed
        ~Projekt();

        //edm filter plugin specific functions
        virtual void beginJob();
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob();

    private:

        edm::ParameterSet theConfig;
        unsigned int theEventCount;
        unsigned int hscpCount;
        unsigned int nlines;
        double hscpMass;
        unsigned int genMCCount;
        unsigned int genMCPassCount;
        unsigned int recCount;
        unsigned int recPassCount;
        unsigned int counter;
        unsigned int COUNTGEN;

        TFile *myRootFile;
        TH1F *histoGenDeltaR;
        TH1F *histoGenDeltaPhi;
        TH1F *histoGenDeltaEta;
        TH1F *histoDeltaR;
        TH1F *histoMinDeltaR;

        TH2F *histoMuonStationsTotal;
        TH2F *histoMuonStationsSelected;
        TH2F *histoEtaCheck;
        TH2F *histoPtCheck;
        TH2F *histoPCheck;
        TH2F *histoInvbetaCheck;

        TH2F *EfficiencyRecoInvbetaEtaAll;
        TH2F *EfficiencyRecoInvbetaEtaSelected;

        TH2F *EfficiencyRecoInvbetaPhiAllBarrel;
        TH2F *EfficiencyRecoInvbetaPhiSelectedBarrel;

        TH2F *EfficiencyRecoInvbetaPhiAllEndcap;
        TH2F *EfficiencyRecoInvbetaPhiSelectedEndcap;

        TH2F *EfficiencyRecoInvbetaPhiAllOverlap;
        TH2F *EfficiencyRecoInvbetaPhiSelectedOverlap;

        TH2F *EfficiencyInvbetaEtaAll;
        TH2F *EfficiencyInvbetaEtaSelected;

        TH2F *EfficiencyInvbetaPhiAllBarrel;
        TH2F *EfficiencyInvbetaPhiSelectedBarrel;

        TH2F *EfficiencyInvbetaPhiAllEndcap;
        TH2F *EfficiencyInvbetaPhiSelectedEndcap;

        TH2F *EfficiencyInvbetaPhiAllOverlap;
        TH2F *EfficiencyInvbetaPhiSelectedOverlap;

        TEfficiency *histoPtEfficiency;
        TEfficiency *histoPEfficiency;
        TEfficiency *histoEtaEfficiency;
        TEfficiency *histoInvbetaEfficiency;
        TEfficiency *histoPtEfficiencyL1T;
        TEfficiency *histoPEfficiencyL1T;
        TEfficiency *histoEtaEfficiencyL1T;
        TEfficiency *histoInvbetaEfficiencyL1T;

        TEfficiency *histoPtEfficiencyL1T_BMTF;
        TEfficiency *histoInvbetaEfficiencyL1T_BMTF;
        TEfficiency *histoPEfficiencyL1T_BMTF;
        
        TEfficiency *histoPtEfficiencyL1T_OMTF;
        TEfficiency *histoInvbetaEfficiencyL1T_OMTF;
        TEfficiency *histoPEfficiencyL1T_OMTF;
        
        TEfficiency *histoPtEfficiencyL1T_EMTF;
        TEfficiency *histoInvbetaEfficiencyL1T_EMTF;
        TEfficiency *histoEtaEfficiencyL1T_EMTF_A;
        TEfficiency *histoEtaEfficiencyL1T_EMTF_B;
        TEfficiency *histoPEfficiencyL1T_EMTF;

        TTree *hscpTree;
        TTree *recoTree;
        TTree *generatedTree;
        TTree *invbetaTree;
        TTree *totalMCTree;
        TTree *gmtMuonTree;

        SimHitData hitData;
        RecCandidateData recData;
        GenCandidateData generatedData;
        InverseBetaData invbetaData;
        TotalMCData totalMCData;
        GMTMuonData gmtMuonData;

  ////////////////////////////////
// Definitions of various inputs // 
  //////////////////////////////
        edm::EDGetTokenT<edm::SimTrackContainer> inputSim;
        edm::EDGetTokenT<edm::SimVertexContainer> inputVtx;
        edm::EDGetTokenT<TrackingParticleCollection> inputTP;
        edm::EDGetTokenT<TrackingVertexCollection> inputTV, inputTV0;
        edm::EDGetTokenT<vector<reco::GenParticle> > inputGP;

//My simHits contribution to the code:
       //edm::EDGetTokenT<vector<PSimHit>> inputHitsDT;
        edm::EDGetTokenT<vector<reco::Muon>> inputRecMuons;
        edm::EDGetTokenT<edm::ValueMap<reco::MuonTimeExtra>> inputRecMuonsExtra;
        const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> theGeometryToken;
        const edm::ESGetToken<DTGeometry, MuonGeometryRecord> theDTGeomToken;
        const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> theCSCGeomToken;
        const edm::ESGetToken<RPCGeometry, MuonGeometryRecord> theRPCGeomToken;

        edm::EDGetTokenT<BXVector<l1t::Muon>> theGmtToken;  

};

//Definitions of various functions

/*
bool match(const TrackingParticle & tp, const l1t::TrackerMuon & gmt) {
  return (   (fabs(tp.pt()-gmt.trkPtr()->momentum().perp()) )/tp.pt() < 0.1
            && fabs(tp.phi()-gmt.trkPtr()->momentum().phi()) < 0.1
            && fabs(tp.eta()-gmt.trkPtr()->momentum().eta()) < 0.1 );
}*/

string print(const TrackingParticle & tp) {
    stringstream ss;
    ss << tp.pdgId()
      <<" pt: "<<tp.pt()
      <<" eta: "<<tp.eta()
      <<" phi: "<<tp.phi()
      <<" vtx[r,z]:  ["<<tp.parentVertex()->position().Rho() <<", "<<tp.parentVertex()->position().z()<<"]"
      <<" time: "<<tp.parentVertex()->position().T(); 
    return ss.str();
}

string print(const reco::GenParticle & gp) {
    stringstream ss;
    ss << gp.pdgId()
       <<" pt: "<<gp.pt()
       <<" eta: "<<gp.eta()
       <<" phi: "<<gp.phi()
       <<" eta test: " << gp.momentum().eta(); 
    return ss.str();
}

double Invbeta(double momentum, double mass){
    double invbeta = sqrt(mass*mass+momentum*momentum)/momentum;
    return invbeta;
}

bool isTightMuon(const reco::Muon& muon){
    if(!muon.isPFMuon() || !muon.isGlobalMuon() ) return false;
    bool muID = muon::isGoodMuon(muon,muon::GlobalMuonPromptTight) &&
                                  (muon.numberOfMatchedStations() > 1);
    bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
                        muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
    bool ip = true;
    return muID && hits && ip;
}

const TrackingParticle & ancestor(const TrackingParticle & particle) {

    const TrackingVertexRef&  tpv = particle.parentVertex(); 
    if (tpv->nSourceTracks() == 0) return particle;
    const TrackingParticle & parent =  **(tpv->sourceTracks_begin());
    return ancestor(parent);
}

Projekt::Projekt(const edm::ParameterSet& conf) 
  : theConfig(conf), theEventCount(0), hscpCount(0), nlines(0), genMCCount(0), genMCPassCount(0), recCount(0), recPassCount(0), counter(0), COUNTGEN(0), theGeometryToken(esConsumes()), theDTGeomToken(esConsumes()), theCSCGeomToken(esConsumes()), theRPCGeomToken(esConsumes())
{
    cout <<" CTORXX" << endl;
//  inputOMTF = consumes<l1t::RegionalMuonCandBxCollection>(theConfig.getParameter<edm::InputTag>("inputOMTF") );
 // inputSim =  consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));
 // inputVtx =  consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits"));
//  inputGMT =  consumes< vector<l1t::TrackerMuon> >(edm::InputTag("gmtMuons"));
    inputTP  =   consumes<TrackingParticleCollection>(edm::InputTag("mix","MergedTrackTruth"));
    inputTV  =   consumes<TrackingVertexCollection>(edm::InputTag("mix","MergedTrackTruth"));
    inputTV0 =   consumes<TrackingVertexCollection>(edm::InputTag("mix","InitialVertices"));
    inputGP  =  consumes< vector<reco::GenParticle> >(edm::InputTag("genParticles"));

  //nputHitsDT = consumes<vector<PSimHit>>(edm::InputTag("g4SimHits","MuonDTHits"));
  inputRecMuons = consumes<vector<reco::Muon>>(edm::InputTag("muons",""));
  inputRecMuonsExtra = consumes<edm::ValueMap<reco::MuonTimeExtra>>(edm::InputTag("muons","combined"));
  
  theGmtToken = consumes<BXVector<l1t::Muon>>(edm::InputTag("gmtStage2Digis","Muon"));
 //Order of input tags needs to be the same as in edmDumpEventContent

}


Projekt::~Projekt() 
{ 
  cout <<" DTOR" << endl;
}

void Projekt::beginJob()
{

  myRootFile = new TFile("stau_M432_analysis.root","RECREATE"); //remember to change name when changing datafiles!
  histoDeltaR = new TH1F("histoDeltaR",";#DeltaR^{GEN/GEN};#Entries",300,0,6);
  histoMinDeltaR = new TH1F("histoMinDeltaR",";#DeltaR;#Events",300,0,6);
  histoGenDeltaR = new TH1F("histoGenDeltaR","histoGenDeltaR",300,0,6);
  histoGenDeltaPhi = new TH1F("histoGenDeltaPhi","histoGenDeltaPhi;#left|#Delta#phi^{GEN/GEN}#right| [rad];#Entries",100,0,3.5);
  histoGenDeltaEta = new TH1F("histoGenDeltaEta","histoGenDeltaEta;#left|#Delta#eta^{GEN/GEN}#right|;#Entries",100,0,6);

  histoMuonStationsTotal = new TH2F("histoMuonStationsTotal","histoMuonStationsTotal", 30, 0, 6, 6 , 0 , 6);
  histoMuonStationsSelected = new TH2F("histoMuonStationsSelected","histoMuonStationsSelected", 30, 0, 6 , 6, 0 , 6);
  
  histoEtaCheck = new TH2F("histoEtaCheck","histoEtaCheck",56,-2.8,2.8,56,-2.8,2.8);
  histoPtCheck = new TH2F("histoPtCheck","histoPtCheck",350,0,3500,350,0,3500);
  histoPCheck = new TH2F("histoPCheck","histoPCheck",500,0,5000,500,0,5000);
  histoInvbetaCheck = new TH2F("histoInvbetaCheck","histoInvbetaCheck",130,0,6.5,130,0,6.5);

//Gen refers to total gen candidates, reco refers to gen candidates matched with reco candidates

  const int  binsPT = 32;
  double binsPTEdges[binsPT+1] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,200.,250.,300.,350.,400.,450.,500.,600.,700.,800.,900.,1000.,1200.,1400.,1600.,1800.,2000.,3000.};

  const int  binsP = 26;
  double binsPEdges[binsP+1] = {0.,20.,40.,60.,80.,100.,150.,200.,250.,300.,350.,400.,450.,500.,600.,700.,800.,900.,1000.,1200.,1400.,1600.,1800.,2000.,3000.,4000.,5000.};

  const int  binsINVB = 19;
  double binsINVBEdges[binsINVB+1] = {1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,5.0,6.0};
  histoPtEfficiency = new TEfficiency("histoPtEfficiency",";p_{T}^{GEN} [GeV];Efficiency",binsPT,binsPTEdges);
  histoPEfficiency = new TEfficiency("histoPEfficiency", ";p^{GEN} [GeV];Efficiency", binsP, binsPEdges);
  histoEtaEfficiency = new TEfficiency("histoEtaEfficiency", ";#eta^{GEN};Efficiency", 100, -2.8, 2.8);
  histoInvbetaEfficiency = new TEfficiency("histoInvbetaEfficiency", ";1/#beta^{GEN};Efficiency", binsINVB, binsINVBEdges);

  histoPtEfficiencyL1T = new TEfficiency("histoPtEfficiencyL1T",";p_{T}^{GEN} [GeV];Efficiency", binsPT, binsPTEdges);
  histoPEfficiencyL1T = new TEfficiency("histoPEfficiencyL1T", ";p^{GEN} [GeV];Efficiency", binsP, binsPEdges);
  histoEtaEfficiencyL1T = new TEfficiency("histoEtaEfficiencyL1T", ";#eta^{GEN};Efficiency", 100, -2.8, 2.8);
  histoInvbetaEfficiencyL1T = new TEfficiency("histoInvbetaEfficiencyL1T", ";1/#beta^{GEN};Efficiency", binsINVB, binsINVBEdges);

  histoPtEfficiencyL1T_BMTF = new TEfficiency("histoPtEfficiencyL1T_BMTF",";p_{T}^{GEN} [GeV];Efficiency", binsPT, binsPTEdges);
  histoInvbetaEfficiencyL1T_BMTF = new TEfficiency("histoInvbetaEfficiencyL1T_BMTF", ";1/#beta^{GEN};Efficiency", binsINVB, binsINVBEdges);
  histoPEfficiencyL1T_BMTF = new TEfficiency("histoPEfficiencyL1T_BMTF", ";p^{GEN} [GeV];Efficiency", binsP, binsPEdges);

  histoPtEfficiencyL1T_OMTF = new TEfficiency("histoPtEfficiencyL1T_OMTF",";p_{T}^{GEN} [GeV];Efficiency", binsPT, binsPTEdges);
  histoInvbetaEfficiencyL1T_OMTF = new TEfficiency("histoInvbetaEfficiencyL1T_OMTF", ";1/#beta^{GEN};Efficiency", binsINVB, binsINVBEdges);
  histoPEfficiencyL1T_OMTF = new TEfficiency("histoPEfficiencyL1T_OMTF", ";p^{GEN} [GeV];Efficiency", binsP, binsPEdges);

  histoPtEfficiencyL1T_EMTF = new TEfficiency("histoPtEfficiencyL1T_EMTF",";p_{T}^{GEN} [GeV];Efficiency", binsPT, binsPTEdges);
  histoInvbetaEfficiencyL1T_EMTF = new TEfficiency("histoInvbetaEfficiencyL1T_EMTF", ";1/#beta^{GEN};Efficiency", binsINVB, binsINVBEdges);
  histoEtaEfficiencyL1T_EMTF_A = new TEfficiency("histoEtaEfficiencyL1T_EMTF_A", ";#eta^{GEN};Efficiency", 100, -2.8, 2.8);
  histoEtaEfficiencyL1T_EMTF_B = new TEfficiency("histoEtaEfficiencyL1T_EMTF_B", ";#eta^{GEN};Efficiency", 100, -2.8, 2.8);
  histoPEfficiencyL1T_EMTF = new TEfficiency("histoPEfficiencyL1T_EMTF", ";p^{GEN} [GeV];Efficiency", binsP, binsPEdges);

  EfficiencyRecoInvbetaEtaAll = new TH2F("EfficiencyRecoInvbetaEtaAll",";1/#beta^{GEN};#eta^{GEN};Efficiency",100,-2.8,2.8,binsINVB, binsINVBEdges);
  EfficiencyRecoInvbetaEtaSelected = new TH2F("EfficiencyRecoInvbetaEtaSelected",";1/#beta^{GEN};#eta^{GEN};Efficiency",100,-2.8,2.8,binsINVB, binsINVBEdges);
  
  EfficiencyRecoInvbetaPhiAllBarrel = new TH2F("EfficiencyRecoInvbetaPhiAllBarrel",";1/#beta^{GEN};#phi^{GEN};Efficiency",12,0,2*M_PI,binsINVB, binsINVBEdges);
  EfficiencyRecoInvbetaPhiSelectedBarrel = new TH2F("EfficiencyRecoInvbetaPhiSelectedBarrel",";1/#beta^{GEN};#phi^{GEN};Efficiency",12,0,2*M_PI,binsINVB, binsINVBEdges);

  EfficiencyRecoInvbetaPhiAllEndcap = new TH2F("EfficiencyRecoInvbetaPhiAllEndcap",";1/#beta^{GEN};#phi^{GEN};Efficiency",18,0,2*M_PI,binsINVB, binsINVBEdges);
  EfficiencyRecoInvbetaPhiSelectedEndcap = new TH2F("EfficiencyRecoInvbetaPhiSelectedEndcap",";1/#beta^{GEN};#phi^{GEN};Efficiency",18,0,2*M_PI,binsINVB, binsINVBEdges);
  
  EfficiencyRecoInvbetaPhiAllOverlap = new TH2F("EfficiencyRecoInvbetaPhiAllOverlap",";1/#beta^{GEN};#phi^{GEN};Efficiency",6,0,2*M_PI,binsINVB, binsINVBEdges);
  EfficiencyRecoInvbetaPhiSelectedOverlap = new TH2F("EfficiencyRecoInvbetaPhiSelectedOverlap",";1/#beta^{GEN};#phi^{GEN};Efficiency",6,0,2*M_PI,binsINVB, binsINVBEdges);


  EfficiencyInvbetaEtaAll = new TH2F("EfficiencyInvbetaEtaAll",";1/#beta^{GEN};#eta^{GEN};Efficiency",100,-2.8,2.8,binsINVB, binsINVBEdges);
  EfficiencyInvbetaEtaSelected = new TH2F("EfficiencyInvbetaEtaSelected",";1/#beta^{GEN};#eta^{GEN};Efficiency",100,-2.8,2.8,binsINVB, binsINVBEdges);
  
  EfficiencyInvbetaPhiAllBarrel = new TH2F("EfficiencyInvbetaPhiAllBarrel",";1/#beta^{GEN};#phi^{GEN};Efficiency",12,0,2*M_PI,binsINVB, binsINVBEdges);
  EfficiencyInvbetaPhiSelectedBarrel = new TH2F("EfficiencyInvbetaPhiSelectedBarrel",";1/#beta^{GEN};#phi^{GEN};Efficiency",12,0,2*M_PI,binsINVB, binsINVBEdges);

  EfficiencyInvbetaPhiAllEndcap = new TH2F("EfficiencyInvbetaPhiAllEndcap",";1/#beta^{GEN};#phi^{GEN};Efficiency",18,0,2*M_PI,binsINVB, binsINVBEdges);
  EfficiencyInvbetaPhiSelectedEndcap = new TH2F("EfficiencyInvbetaPhiSelectedEndcap",";1/#beta^{GEN};#phi^{GEN};Efficiency",18,0,2*M_PI,binsINVB, binsINVBEdges);

  EfficiencyInvbetaPhiAllOverlap = new TH2F("EfficiencyInvbetaPhiAllOverlap",";1/#beta^{GEN};#phi^{GEN};Efficiency",6,0,2*M_PI,binsINVB, binsINVBEdges);
  EfficiencyInvbetaPhiSelectedOverlap = new TH2F("EfficiencyInvbetaPhiSelectedOverlap",";1/#beta^{GEN};#phi^{GEN};Efficiency",6,0,2*M_PI,binsINVB, binsINVBEdges);

  ///////////////////
// TTree definitions //
  //////////////////

  recoTree = new TTree("recoTree", "recoTree");
  recoTree->Branch("PDGID", &(recData.pdgId),32000, 99);
  recoTree->Branch("EVENT", &(recData.event),32000, 99);
  recoTree->Branch("PHI", &(recData.phi),32000, 99);
  recoTree->Branch("THETA", &(recData.theta),32000, 99);
  recoTree->Branch("PT", &(recData.pT),32000, 99);
  recoTree->Branch("ETA", &(recData.eta),32000, 99);
  recoTree->Branch("MASS", &(recData.mass),32000, 99);
  recoTree->Branch("VX", &(recData.vx),32000, 99);
  recoTree->Branch("VY", &(recData.vy),32000, 99);
  recoTree->Branch("VZ", &(recData.vz),32000, 99);
  recoTree->Branch("BETA",&(recData.beta),32000,99);
  recoTree->Branch("INVBETA",&(recData.invbeta),32000,99);

  generatedTree = new TTree("generatedTree", "generatedTree");
  generatedTree->Branch("PDGID", &(generatedData.pdgId),32000, 99);
  generatedTree->Branch("EVENT", &(generatedData.event),32000, 99);
  generatedTree->Branch("PHI", &(generatedData.phi),32000, 99);
  generatedTree->Branch("THETA", &(generatedData.theta),32000, 99);
  generatedTree->Branch("PT", &(generatedData.pT),32000, 99);
  generatedTree->Branch("ETA", &(generatedData.eta),32000, 99);
  generatedTree->Branch("MASS", &(generatedData.mass),32000, 99);
  generatedTree->Branch("BETA", &(generatedData.beta),32000, 99);
  generatedTree->Branch("INVBETA", &(generatedData.invbeta),32000, 99);

  invbetaTree = new TTree("invbetaTree", "invbetaTree");
  invbetaTree->Branch("PGEN",&(invbetaData.pgen), 32000, 99);
  invbetaTree->Branch("PRECO",&(invbetaData.preco), 32000, 99);
  invbetaTree->Branch("INVBETAGEN", &(invbetaData.invbetagen), 32000, 99);
  invbetaTree->Branch("INVBETAMANUAL", &(invbetaData.invbetamanual), 32000, 99);
  invbetaTree->Branch("INVBETARECO", &(invbetaData.invbetareco), 32000,99);
  invbetaTree->Branch("INVBETARECOERROR",&(invbetaData.invbetarecoerror), 32000,99);
  invbetaTree->Branch("EVENT", &(invbetaData.event), 32000, 99);
  invbetaTree->Branch("DELTAR", &(invbetaData.deltaR), 32000, 99);
  invbetaTree->Branch("GENPT",&(invbetaData.genPt), 32000, 99);
  invbetaTree->Branch("INNERRECOPT", &(invbetaData.innerrecoPt), 32000,99);
  invbetaTree->Branch("OUTERRECOPT", &(invbetaData.outerrecoPt), 32000,99);
  invbetaTree->Branch("RECOPT",&(invbetaData.recoPt), 32000,99);
  invbetaTree->Branch("RECOETA",&(invbetaData.recoEta), 32000,99);
  invbetaTree->Branch("GENETA",&(invbetaData.genEta), 32000,99);
  invbetaTree->Branch("RECOPHI",&(invbetaData.recoPhi), 32000,99);
  invbetaTree->Branch("GENPHI",&(invbetaData.genPhi), 3200,99);
  invbetaTree->Branch("STMATCH",&(invbetaData.stationMatches), 32000,99);

  totalMCTree = new TTree("totalMCTree", "totalMCTree");
  totalMCTree->Branch("PT", &(totalMCData.pt), 32000, 99);
  totalMCTree->Branch("P", &(totalMCData.p), 32000, 99);
  totalMCTree->Branch("ETA", &(totalMCData.eta), 32000, 99);
  totalMCTree->Branch("INVBETA", &(totalMCData.invbeta), 32000, 99);

  gmtMuonTree = new TTree("gmtMuonTree", "gmtMuonTree");
  gmtMuonTree->Branch("PT", &(gmtMuonData.pt), 32000, 99);
  gmtMuonTree->Branch("PTGEN", &(gmtMuonData.ptGen), 32000, 99);
  gmtMuonTree->Branch("P", &(gmtMuonData.p), 32000, 99);
  gmtMuonTree->Branch("PGEN", &(gmtMuonData.pGen), 32000, 99);
  gmtMuonTree->Branch("ETA", &(gmtMuonData.eta), 32000, 99);
  gmtMuonTree->Branch("ETAGEN", &(gmtMuonData.etaGen), 32000, 99);
  gmtMuonTree->Branch("INVBETA", &(gmtMuonData.invbeta), 32000, 99);
  gmtMuonTree->Branch("INVBETAGEN", &(gmtMuonData.invbetaGen), 32000);
  gmtMuonTree->Branch("DELTAR",&(gmtMuonData.deltaR), 32000, 99);

  ////////////////////////////  //////////////////////////////
// PSEUDORAPIDITY HISTOGRAMS SPECIFYING MUON SYSTEM REGIONS //
  //////////////////////////////////////////////////////////
  //cout << "HERE Projekt::beginJob()" << endl;
}

void Projekt::endJob()
{
  cout << "MC HSCP: " << genMCCount << "PASS: " << genMCPassCount << endl;
  cout << "REC MUON: " << recCount << "PASS: " << recPassCount << endl;
  cout << "TEfficiency entries: " << counter << endl;
  cout << "Stau gen in detector: " << COUNTGEN << endl;
  cout << "Total nr of events: " << theEventCount << endl;

  recoTree->Write();
  generatedTree->Write();
  invbetaTree->Write();
  totalMCTree->Write();
  gmtMuonTree->Write();

  histoDeltaR->Write();
  histoMinDeltaR->Write();
  histoGenDeltaR->Write();
  histoGenDeltaPhi->Write();
  histoGenDeltaEta->Write();

  histoEtaCheck->Write();
  histoPtCheck->Write();
  histoPCheck->Write();
  histoInvbetaCheck->Write();

  EfficiencyRecoInvbetaEtaAll->Write();
  EfficiencyRecoInvbetaEtaSelected->Write();

  EfficiencyRecoInvbetaPhiAllBarrel->Write();
  EfficiencyRecoInvbetaPhiSelectedBarrel->Write();

  EfficiencyRecoInvbetaPhiAllEndcap->Write();
  EfficiencyRecoInvbetaPhiSelectedEndcap->Write();

  EfficiencyRecoInvbetaPhiAllOverlap->Write();
  EfficiencyRecoInvbetaPhiSelectedOverlap->Write(); 
  
  EfficiencyInvbetaEtaAll->Write();
  EfficiencyInvbetaEtaSelected->Write();

  EfficiencyInvbetaPhiAllBarrel->Write();
  EfficiencyInvbetaPhiSelectedBarrel->Write();

  EfficiencyInvbetaPhiAllEndcap->Write();
  EfficiencyInvbetaPhiSelectedEndcap->Write();

  EfficiencyInvbetaPhiAllOverlap->Write();
  EfficiencyInvbetaPhiSelectedOverlap->Write();   

  histoPtEfficiency->Write();
  histoPEfficiency->Write();
  histoEtaEfficiency->Write();
  histoInvbetaEfficiency->Write();

  histoPtEfficiencyL1T->Write();
  histoPEfficiencyL1T->Write();
  histoEtaEfficiencyL1T->Write();
  histoInvbetaEfficiencyL1T->Write();

  histoPtEfficiencyL1T_BMTF->Write();
  histoInvbetaEfficiencyL1T_BMTF->Write();
  histoPEfficiencyL1T_BMTF->Write();
  
  histoPtEfficiencyL1T_OMTF->Write();
  histoInvbetaEfficiencyL1T_OMTF->Write();
  histoPEfficiencyL1T_OMTF->Write();
    
  histoPtEfficiencyL1T_EMTF->Write();
  histoInvbetaEfficiencyL1T_EMTF->Write();
  histoPEfficiencyL1T_EMTF->Write();
  histoEtaEfficiencyL1T_EMTF_A->Write();
  histoEtaEfficiencyL1T_EMTF_B->Write();

  histoMuonStationsTotal->Write();
  histoMuonStationsSelected->Write();         

  myRootFile->Close();

  delete myRootFile;
  //cout << "HERE Cwiczenie::endJob()" << endl;
}

void Projekt::analyze(
    const edm::Event& ev, const edm::EventSetup& es){
    //cout << " -------------------------------- HERE Cwiczenie::analyze "<< endl;
    cout <<"*** Cwiczenie, analyze event: " << ev.id()<<" useful event count:"<<++theEventCount << endl;

  //  edm::Handle<vector<l1t::TrackerMuon> > gmtColl;
  //  ev.getByToken(inputGMT, gmtColl);
  //  const vector<l1t::TrackerMuon> & gmtMuons = *gmtColl.product();
  //  histo->Fill(gmtMuons.size());

    auto const & globalGeometry = es.getData(theGeometryToken);  //assign global geometry

      //////////////////////
    // INITIAL GEN ANALYSIS //
      /////////////////////

    const vector<reco::GenParticle> & genParticles = ev.get(inputGP);
    for (const auto & gp : genParticles) {
        if (abs(gp.pdgId())>1000000 && gp.status()==1){ //BSM+stable
            genMCCount++;
            double hscpID = gp.pdgId();
            double hscpP = gp.p();
            double hscpPt = gp.pt();
            double hscpPhi = gp.phi();
            double hscpTheta = gp.theta();
            double hscpEta = gp.eta();
            double hscpMass = gp.mass(); //predefined!
            double hscpInvbeta = Invbeta(hscpP,hscpMass);
            if(hscpPhi <= 0) hscpPhi += 2*M_PI;    
            if(hscpPt<50 || abs(hscpEta)>2.1) continue; //2.1 is regarded as the point up to which the trigger works well
            genMCPassCount++;  
            generatedData.pdgId = hscpID;
            generatedData.event = theEventCount;
            generatedData.phi = hscpPhi;
            generatedData.invbeta = hscpInvbeta;
            generatedData.theta = hscpTheta;
            generatedData.pT = hscpPt;
            generatedData.eta = hscpEta;
            generatedData.mass = hscpMass;
            generatedTree->Fill();
        }
    }

    //////////////////
  // RECHIT ANALYSIS //
    ////////////////

    edm::Handle<vector<reco::Muon>> muonHandle;
    ev.getByToken(inputRecMuons, muonHandle);
    const vector<reco::Muon> & recoMuons = *muonHandle.product();

    edm::Handle<edm::ValueMap<reco::MuonTimeExtra>> muonExtraHandle;
    ev.getByToken(inputRecMuonsExtra,muonExtraHandle);
    edm::ValueMap<reco::MuonTimeExtra> muonValueMap = *muonExtraHandle;

    for(vector<reco::Muon>::const_iterator iter = recoMuons.begin();iter<recoMuons.end();iter++){
        recCount++;
        const reco::Muon & recoCandidate = *iter;
        double recoPhi = recoCandidate.phi();
        if(recoPhi <= 0) recoPhi += 2*M_PI;
        double recoTheta = recoCandidate.theta();
        double recoPt = recoCandidate.pt();
        double recoP = recoCandidate.p();
        double recoEta = recoCandidate.eta();
        int  recoID = recoCandidate.pdgId();
        double recoInvbeta = Invbeta(recoP,hscpMass);
        if(recoPt<50 || abs(recoEta)>2.1 || isTightMuon(recoCandidate)==false) continue;
        recPassCount++;
        recData.pdgId = recoID;
        recData.event = theEventCount;
        recData.phi = recoPhi;
        recData.theta = recoTheta;
        recData.pT = recoPt;
        recData.eta = recoEta;
        recData.invbeta = recoInvbeta;
        recData.vx = recoCandidate.vx();
        recData.vy = recoCandidate.vy();
        recData.vz = recoCandidate.vz();
        recoTree->Fill();
    }

    ////////////////////////////////
  //  MATCHING (GEN/GEN + GEN/RECO) // 
    ///////////////////////////////
    bool gpMatchCount=false;  
    for (const auto & gp : genParticles){          
        double hscpID = gp.pdgId();
        double hscpP = gp.p();
        double hscpPt = gp.pt();
        double hscpPhi = gp.phi();
        double hscpTheta = gp.theta();
        double hscpEta = gp.eta();
        double hscpMass = gp.mass(); //predefined!
        double hscpInvbeta = Invbeta(hscpP,hscpMass);
        if(hscpPt<50 || abs(hscpEta)>2.1 || abs(hscpID)<1000000 || gp.status()!=1) continue;
        for (const auto & gp2 : genParticles){
            double hscpID_2 = gp2.pdgId();
            double hscpPt_2 = gp2.pt();
            double hscpPhi_2 = gp2.phi();
            double hscpEta_2 = gp2.eta();
            if(hscpPt_2<50 || abs(hscpEta_2)>2.1 || abs(hscpID_2)<1000000 || gp2.status()!=1 || hscpID==hscpID_2) continue;
            if(gpMatchCount==true) continue;
            gpMatchCount=true;
            double deltaR = reco::deltaR(gp,gp2);
            double deltaPhi = abs(hscpPhi - hscpPhi_2);
            if(deltaPhi>M_PI) deltaPhi -= 2*M_PI;
            deltaPhi = abs(deltaPhi);
            double deltaEta = abs(hscpEta - hscpEta_2);
            histoGenDeltaR->Fill(deltaR);
            histoGenDeltaPhi->Fill(deltaPhi);
            histoGenDeltaEta->Fill(deltaEta);
        }
    }
        
    for (const auto & gp : genParticles){ 
        double minDeltaR = 100;
        int  recoMuonCount = 0;
        int hscpID = gp.pdgId();
        double hscpP = gp.p();
        double hscpPt = gp.pt();
        double hscpEta = gp.eta();
        double hscpPhi = gp.phi();
        double hscpMass = gp.mass(); //predefined!
        double hscpInvbeta = Invbeta(hscpP,hscpMass);  
        if(hscpPt<50 || abs(hscpEta)>2.1 || abs(hscpID)<1000000 || gp.status()!=1) continue;
        for(vector<reco::Muon>::const_iterator iter = recoMuons.begin();iter<recoMuons.end();iter++){
            reco::MuonRef muref(muonHandle,recoMuonCount);
            recoMuonCount++;
            reco::MuonTimeExtra muonExtraData = (muonValueMap)[muref];
            const reco::Muon & recoCandidate = *iter;
            double recoPt = recoCandidate.pt();
            double recoEta = recoCandidate.eta();
            double recoP = recoCandidate.p();
            double recoPhi = recoCandidate.phi();
            double recInvBetaMan = Invbeta(recoP,hscpMass);
            double recInvBetaRec = muonExtraData.inverseBeta();
            double recInvBetaRecError = muonExtraData.inverseBetaErr();
            if(recoPt<50 || abs(recoEta)>2.1 || isTightMuon(recoCandidate)==false || abs(recoCandidate.pdgId())!=13) continue;
            double deltaR = reco::deltaR(gp,recoCandidate);
            histoDeltaR->Fill(deltaR);
            if(deltaR<minDeltaR) minDeltaR = deltaR;
            invbetaData.event = theEventCount;
            invbetaData.pgen = hscpP;
            invbetaData.preco = recoP;
            invbetaData.genPt = hscpPt;
            invbetaData.recoPt = recoPt;
            invbetaData.invbetagen = hscpInvbeta;
            invbetaData.invbetamanual = recInvBetaMan;
            invbetaData.invbetareco = recInvBetaRec;
            invbetaData.invbetarecoerror = recInvBetaRecError;
            invbetaData.deltaR = deltaR;
            invbetaData.genEta = hscpEta;
            invbetaData.recoEta = recoEta;
            invbetaData.recoPhi = recoPhi;
            invbetaData.genPhi = hscpPhi;
            int matchedStations = recoCandidate.numberOfMatchedStations();
            invbetaData.stationMatches = matchedStations;
            invbetaTree->Fill();       
        }
        histoMinDeltaR->Fill(minDeltaR);
    }

    //////////////////
  // EFFICIENCY PLOTS //
    //////////////////
    for (const auto & gp : genParticles) {
        if(abs(gp.pdgId())<1000000 || gp.status()!=1 || abs(gp.eta())>2.1 || gp.pt()<50) continue;
        COUNTGEN++;
        double hscpP = gp.p();
        double hscpPt = gp.pt();
        double hscpMass = gp.mass();
        double hscpInvbeta = Invbeta(hscpP,hscpMass);
        double hscpPhi = gp.phi();
        if (hscpPhi <= 0) hscpPhi += 2*M_PI;
        double hscpEta = gp.eta();
        histoMuonStationsTotal->Fill(hscpInvbeta,0);
        //histoMuonStationsTotal->Fill(invbeta,1);
        histoMuonStationsTotal->Fill(hscpInvbeta,2);
        histoMuonStationsTotal->Fill(hscpInvbeta,3);
        histoMuonStationsTotal->Fill(hscpInvbeta,4);
        histoMuonStationsTotal->Fill(hscpInvbeta,5);
        bool rec = false;
        for(vector<reco::Muon>::const_iterator iter = recoMuons.begin(); iter<recoMuons.end();iter++){
            const reco::Muon & recoCandidate = *iter;  
            double deltaR = reco::deltaR(gp,recoCandidate);
            double recoEta = recoCandidate.eta();
            double recoPt = recoCandidate.pt();
            if(recoPt<50 || abs(recoEta)>2.1 || isTightMuon(recoCandidate)==false || abs(recoCandidate.pdgId())!=13 || deltaR>0.02) continue;
            rec=true;      
        }
        histoPtEfficiency->Fill(rec,hscpPt);
        histoPEfficiency->Fill(rec,hscpP);
        histoEtaEfficiency->Fill(rec,hscpEta);
        histoInvbetaEfficiency->Fill(rec,hscpInvbeta);
        counter++;
        if(abs(hscpEta)<2.1){ //not going to consider 2.4 at all
            EfficiencyRecoInvbetaEtaAll->Fill(hscpEta,hscpInvbeta); //regardless of eta value here
            if(rec==true)EfficiencyRecoInvbetaEtaSelected->Fill(hscpEta,hscpInvbeta);
        }

        if(abs(hscpEta)<0.83){
            EfficiencyRecoInvbetaPhiAllBarrel->Fill(hscpPhi,hscpInvbeta);
            if(rec==true) EfficiencyRecoInvbetaPhiSelectedBarrel->Fill(hscpPhi,hscpInvbeta);
        }

        else if(abs(hscpEta)>0.83 && abs(hscpEta)<1.24){
            EfficiencyRecoInvbetaPhiAllOverlap->Fill(hscpPhi,hscpInvbeta);
            if(rec==true) EfficiencyRecoInvbetaPhiSelectedOverlap->Fill(hscpPhi,hscpInvbeta);
        }

        else if(abs(hscpEta)>1.24 && abs(hscpEta)<2.1){
            EfficiencyRecoInvbetaPhiAllEndcap->Fill(hscpPhi,hscpInvbeta);
            if(rec==true) EfficiencyRecoInvbetaPhiSelectedEndcap->Fill(hscpPhi,hscpInvbeta);
        }
    }

    //////////////////////////////
  // TRIGGER EFFICIENCY ANALYSIS //
    ////////////////////////////

    int bxNumber = 0;
    const l1t::MuonBxCollection & gmts = ev.get(theGmtToken); 
    edm::Handle<BXVector<l1t::Muon>> l1tMuonHandle;
    ev.getByToken(theGmtToken, l1tMuonHandle);
    const BXVector<l1t::Muon> & l1tMuons = *l1tMuonHandle.product();

    for (const auto & gp : genParticles){
        double hscpP = gp.p();
        double hscpMass = gp.mass();
        double hscpInvbeta = Invbeta(hscpP,hscpMass);
        double hscpPhi = gp.phi();
        double hscpPt = gp.pt();
        double hscpEta = gp.eta();
        if (hscpPhi <= 0) hscpPhi += 2*M_PI;
        bool l1t = false;
        if(abs(gp.pdgId())<1000000 || gp.status()!=1 || hscpPt<50) continue;
        for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it){
            const l1t::Muon & l1tMuon = *it;  
            double l1tP = l1tMuon.p();
            double l1tPt = l1tMuon.pt();
            double l1tEta = l1tMuon.eta();
            double l1tMuonInvbeta = Invbeta(l1tP,hscpMass);
            double deltaR = reco::deltaR(gp,l1tMuon);
            int l1tQuality = l1tMuon.hwQual();
            gmtMuonData.pt = l1tPt;
            gmtMuonData.ptGen = hscpPt;
            gmtMuonData.p = l1tP;
            gmtMuonData.pGen = hscpP;
            gmtMuonData.eta = l1tEta;
            gmtMuonData.etaGen = hscpEta;
            gmtMuonData.invbeta = l1tMuonInvbeta;
            gmtMuonData.invbetaGen = hscpInvbeta;
            gmtMuonData.deltaR = deltaR;
            gmtMuonTree->Fill();
            if(l1tQuality<12 || deltaR > 0.26) continue;
            l1t = true;

        }
        if(abs(hscpEta)<2.1){ //not going to consider 2.4 at all, keeping to 2.1 even on eta efficiency plot!
            histoEtaEfficiencyL1T->Fill(l1t,hscpEta);
            histoPtEfficiencyL1T->Fill(l1t,hscpPt);
            histoPEfficiencyL1T->Fill(l1t,hscpP);
            histoInvbetaEfficiencyL1T->Fill(l1t,hscpInvbeta);
        }
        
        if(abs(hscpEta)<2.1){ //not going to consider 2.4 at all
            EfficiencyInvbetaEtaAll->Fill(hscpEta,hscpInvbeta); //regardless of eta value here
            if(l1t==true)EfficiencyInvbetaEtaSelected->Fill(hscpEta,hscpInvbeta);
            if(hscpInvbeta>1.0 && hscpInvbeta < 1.1) histoEtaEfficiencyL1T_EMTF_A->Fill(l1t,hscpEta);
            if(hscpInvbeta>1.2 && hscpInvbeta < 1.3) histoEtaEfficiencyL1T_EMTF_B->Fill(l1t,hscpEta);

        }

        if(abs(hscpEta)<0.83){
            histoPtEfficiencyL1T_BMTF->Fill(l1t,hscpPt);
            histoPEfficiencyL1T_BMTF->Fill(l1t,hscpP);
            histoInvbetaEfficiencyL1T_BMTF->Fill(l1t,hscpInvbeta);
            EfficiencyInvbetaPhiAllBarrel->Fill(hscpPhi,hscpInvbeta);
            if(l1t==true) EfficiencyInvbetaPhiSelectedBarrel->Fill(hscpPhi,hscpInvbeta);
        }

        else if(abs(hscpEta)>0.83 && abs(hscpEta)<1.24){
            histoPtEfficiencyL1T_OMTF->Fill(l1t,hscpPt);
            histoPEfficiencyL1T_OMTF->Fill(l1t,hscpP);
            histoInvbetaEfficiencyL1T_OMTF->Fill(l1t,hscpInvbeta);
            EfficiencyInvbetaPhiAllOverlap->Fill(hscpPhi,hscpInvbeta);
            if(l1t==true) EfficiencyInvbetaPhiSelectedOverlap->Fill(hscpPhi,hscpInvbeta);
        }

        else if(abs(hscpEta)>1.24 && abs(hscpEta)<2.1){
            histoPtEfficiencyL1T_EMTF->Fill(l1t,hscpPt);
            histoPEfficiencyL1T_EMTF->Fill(l1t,hscpP);
            histoInvbetaEfficiencyL1T_EMTF->Fill(l1t,hscpInvbeta);
            EfficiencyInvbetaPhiAllEndcap->Fill(hscpPhi,hscpInvbeta);
            if(l1t==true) EfficiencyInvbetaPhiSelectedEndcap->Fill(hscpPhi,hscpInvbeta);
        }
    }   
}

DEFINE_FWK_MODULE(Projekt);

