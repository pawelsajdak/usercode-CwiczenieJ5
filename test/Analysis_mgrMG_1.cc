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
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
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
  double distanceVTX;
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
} SimHitCompData;

typedef struct{
  int event;
  int track;
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
  double distanceIP;
  double distanceVTX;
  double VTX_x;
  double VTX_y;
  double VTX_z;
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
  //added variables
  TFile *myRootFile;

  TTree *hscpTree;
  TTree *generatedTree;
  TTree *totalMCTree;
  TTree *hscpCompTree;

  SimHitData hitData;
  SimHitCompData hitCompData;
  GenCandidateData generatedData;
  TotalMCData totalMCData;

  ////////////////////////////////
// Definitions of various inputs // 
  //////////////////////////////
  edm::EDGetTokenT<edm::SimTrackContainer> inputSim;
  edm::EDGetTokenT<edm::SimVertexContainer> inputVtx;
  edm::EDGetTokenT<TrackingParticleCollection> inputTP;
  edm::EDGetTokenT<TrackingVertexCollection> inputTV, inputTV0;
  edm::EDGetTokenT<vector<reco::GenParticle> > inputGP;

//My simHits contribution to the code:
  edm::EDGetTokenT<vector<PSimHit>> inputHitsDT;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> theGeometryToken;
  const edm::ESGetToken<DTGeometry, MuonGeometryRecord> theDTGeomToken;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> theCSCGeomToken;
  const edm::ESGetToken<RPCGeometry, MuonGeometryRecord> theRPCGeomToken;

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
     <<" time: "<<tp.parentVertex()->position().T() 
     ; 
  return ss.str();
}

string print(const reco::GenParticle & gp) {
  stringstream ss;
  ss << gp.pdgId()
     <<" pt: "<<gp.pt()
     <<" eta: "<<gp.eta()
     <<" phi: "<<gp.phi()
     <<" eta test: " << gp.momentum().eta();
     ; 
  return ss.str();
}

const TrackingParticle & ancestor(const TrackingParticle & particle) {

  const TrackingVertexRef&  tpv = particle.parentVertex(); 
  if (tpv->nSourceTracks() == 0) return particle;
  const TrackingParticle & parent =  **(tpv->sourceTracks_begin());
  return ancestor(parent);
}

Projekt::Projekt(const edm::ParameterSet& conf) 
  : theConfig(conf), theEventCount(0), hscpCount(0), nlines(0), hscpMass(0), genMCCount(0), genMCPassCount(0), recCount(0), recPassCount(0), counter(0), COUNTGEN(0), theGeometryToken(esConsumes()), theDTGeomToken(esConsumes()), theCSCGeomToken(esConsumes()), theRPCGeomToken(esConsumes())
{
  cout <<" CTORXX" << endl;
//  inputOMTF = consumes<l1t::RegionalMuonCandBxCollection>(theConfig.getParameter<edm::InputTag>("inputOMTF") );
  inputSim =  consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));
  inputVtx =  consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits"));
  inputTP  =   consumes<TrackingParticleCollection>(edm::InputTag("mix","MergedTrackTruth"));
  inputTV  =   consumes<TrackingVertexCollection>(edm::InputTag("mix","MergedTrackTruth"));
  inputTV0 =   consumes<TrackingVertexCollection>(edm::InputTag("mix","InitialVertices"));
  inputGP  =  consumes< vector<reco::GenParticle> >(edm::InputTag("genParticles"));

  inputHitsDT = consumes<vector<PSimHit>>(edm::InputTag("g4SimHits","MuonDTHits"));
   //Order of input tags needs to be the same as in edmDumpEventContent

}


Projekt::~Projekt() 
{ 
  cout <<" DTOR" << endl;
}

void Projekt::beginJob()
{

  /////////////////////////////////////////////////////////
// FILES FOR PSEUDORAPIDITY ANALYSIS - files use simTracks //
  /////////////////////////////////////////////////////////

//FEVTSIM_stau_M200_full.root - FEVTSIM.root analysis
//stau_M432_analysis.root - stau_432.root analysis (GEN AND RECO ONLY)

  myRootFile = new TFile("stau_M432_analysis.root","RECREATE"); //remember to change name when changing datafiles!

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

  hscpCompTree = new TTree("simCompTree","simCompTree");
  hscpCompTree->Branch("EV", &(hitCompData.event), 32000, 99);
  hscpCompTree->Branch("TR", &(hitCompData.track), 32000, 99);
  hscpCompTree->Branch("NX_L", &(hitCompData.nextLocalX), 32000, 99);
  hscpCompTree->Branch("NY_L", &(hitCompData.nextLocalY), 32000, 99);
  hscpCompTree->Branch("NZ_L", &(hitCompData.nextLocalZ), 32000, 99);
  hscpCompTree->Branch("NX_G", &(hitCompData.nextGlobalX), 32000, 99);
  hscpCompTree->Branch("NY_G", &(hitCompData.nextGlobalY), 32000, 99);
  hscpCompTree->Branch("NZ_G", &(hitCompData.nextGlobalZ), 32000, 99);
  hscpCompTree->Branch("DISTDIFF_G", &(hitCompData.globalHitDistDiff), 32000, 99);
  hscpCompTree->Branch("NBETA", &(hitCompData.betaNextHit), 32000, 99);
  hscpCompTree->Branch("GEANTTOFDIFF", &(hitCompData.diffTOFGeant), 32000, 99);
  hscpCompTree->Branch("HITPTOFDIFF", &(hitCompData.diffTOFHitB), 32000, 99);
  hscpCompTree->Branch("NEXTPTOFFDIFF", &(hitCompData.diffTOFNextB), 32000, 99);
  hscpCompTree->Branch("SAMEWHEEL", &(hitCompData.samewheel), 32000, 99);
  hscpCompTree->Branch("SAMESL", &(hitCompData.samesl), 32000, 99);
  hscpCompTree->Branch("SAMESECT", &(hitCompData.samesect), 32000, 99);
  hscpCompTree->Branch("SAMEST", &(hitCompData.samest), 32000, 99);
  hscpCompTree->Branch("DISTVTX", &(hitCompData.distanceVTX), 32000, 99);


  hscpTree = new TTree("simTree", "simTree");
  hscpTree->Branch("EV", &(hitData.event), 32000, 99);
  hscpTree->Branch("TR", &(hitData.track), 32000, 99);
  hscpTree->Branch("PID", &(hitData.pid), 32000, 99);
  hscpTree->Branch("SL", &(hitData.superlayer), 32000, 99);
  hscpTree->Branch("ST", &(hitData.station), 32000, 99);
  hscpTree->Branch("WH", &(hitData.wheel), 32000, 99);
  hscpTree->Branch("SECT", &(hitData.sector), 32000, 99);
  hscpTree->Branch("LYR", &(hitData.layer), 32000, 99);
  hscpTree->Branch("X_L", &(hitData.localX), 32000, 99);
  hscpTree->Branch("Y_L", &(hitData.localY), 32000, 99);
  hscpTree->Branch("Z_L", &(hitData.localZ), 32000, 99);
  hscpTree->Branch("X_G", &(hitData.globalX), 32000, 99);
  hscpTree->Branch("Y_G", &(hitData.globalY), 32000, 99);
  hscpTree->Branch("Z_G", &(hitData.globalZ), 32000, 99);
  hscpTree->Branch("PT", &(hitData.pt), 32000, 99);
  hscpTree->Branch("DISTVTX", &(hitData.distanceVTX), 32000, 99);
  hscpTree->Branch("PHI", &(hitData.phi), 32000, 99);
  hscpTree->Branch("R", &(hitData.r), 32000, 99);
  hscpTree->Branch("ETA", &(hitData.eta), 32000, 99);
  hscpTree->Branch("TOF", &(hitData.tof), 32000, 99);
  hscpTree->Branch("BETA", &(hitData.beta), 32000, 99);
  hscpTree->Branch("DISTIP", &(hitData.distanceIP), 32000, 99);
  hscpTree->Branch("VTX_X", &(hitData.VTX_x), 32000, 99);
  hscpTree->Branch("VTX_Y", &(hitData.VTX_y), 32000, 99);
  hscpTree->Branch("VTX_Z", &(hitData.VTX_z), 32000, 99);

  totalMCTree = new TTree("totalMCTree", "totalMCTree");
  totalMCTree->Branch("PT", &(totalMCData.pt), 32000, 99);
  totalMCTree->Branch("P", &(totalMCData.p), 32000, 99);
  totalMCTree->Branch("ETA", &(totalMCData.eta), 32000, 99);
  totalMCTree->Branch("INVBETA", &(totalMCData.invbeta), 32000, 99);
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

  totalMCTree->Write();
  hscpTree->Write();  
  hscpCompTree->Write();
  generatedTree->Write();  

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
    
  ////////////////////////////////////////
// Assigning global geometry to analysis //
  //////////////////////////////////////

  auto const & globalGeometry = es.getData(theGeometryToken);  

  ////////////////////////////////////////////////////////////////////////////////////////////////////
// GENERATED PARTICLE ANALYSIS (MC LEVEL)                                                             //
// Definition of genParticle vector, which is needed at this stage to extract generated particle mass //
// (Value of mass for all HSCPs should be the same)                                                   //
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  const vector<reco::GenParticle> & genParticles = ev.get(inputGP);
  int  genStauCandidates = 0;

  for (const auto & gp : genParticles) {
  	if (abs(gp.pdgId()) == 13 && gp.status() == 1){ //generated particle is a muon (which must be stable but that's a given)
      double muon_p = gp.pt()*cosh(gp.eta());
			double muon_phi = gp.phi();
      if (muon_phi <= 0){
				muon_phi += 2*M_PI;
			}
			double muon_beta = muon_p/gp.energy();
    }
    else if (abs(gp.pdgId())>1000000) { //generated particle is BSM particle (e.g. stau->pdgID=1000015
			if(gp.status()==1){ // particle is stable
        genMCCount++;
        hscpMass = gp.mass();
				double hscp_p = gp.pt()*cosh(gp.eta());
				double hscp_phi = gp.phi();
        double hscp_theta = gp.theta();
				if(hscp_phi <= 0){
					hscp_phi += 2*M_PI;
				}
				double hscp_beta = hscp_p/gp.energy();
      
      if(gp.pt()>50 && abs(gp.eta())<2.1 ){
      genMCPassCount++;  
      genStauCandidates++;
      generatedData.pdgId = gp.pdgId();
      generatedData.event = theEventCount;
      generatedData.phi = hscp_phi;
      generatedData.beta = hscp_beta;
      generatedData.invbeta = 1/hscp_beta;
      generatedData.theta = hscp_theta;
      generatedData.pT = gp.pt();
      generatedData.eta = gp.eta();
      generatedData.mass = hscpMass;
      generatedTree->Fill();
      }
    }
  }
  if(gp.status()==1){
    totalMCData.pt = gp.pt();
    totalMCData.p = gp.p();
    totalMCData.eta = gp.eta();
    totalMCData.invbeta = gp.energy()/gp.p();
    totalMCTree->Fill();
  }
}
   

  ///////////////////////////////////
// SimTrack/SimVtx vector definitions //
  //////////////////////////////////

  edm::Handle<edm::SimTrackContainer> simTrk;
  ev.getByToken(inputSim, simTrk);
  const vector<SimTrack>  & mySimTracks = *(simTrk.product());
  //cout <<" SIMULATED TRACKS: "<<mySimTracks.size()<<endl;

  edm::Handle<edm::SimVertexContainer> simVtx;
  ev.getByToken(inputVtx, simVtx);
  const vector<SimVertex> & mySimVerts= *(simVtx.product());
  //cout <<" SIMULATED VERTICES: "<<mySimVerts.size()<<endl;

  //////////////////
// SIMHIT ANALYSIS //
  ////////////////

 const vector<PSimHit> & simDTHits = ev.get(inputHitsDT);
  //cout <<"Number of simulated DT hits in event: "<<simDTHits.size() << endl;

  for(const auto &track:mySimTracks){
    //cout << track.trackId();
    //cout << "Track: " <<track<< endl;
    //cout << track.trackId() << ", " << track.vertIndex() << ", " << endl;
    //cout << "Position rho = " << mySimVerts[track.vertIndex()].position().Rho() << ", Position z = " << mySimVerts[track.vertIndex()].position().z() << endl;
    //cout << "Parent index: " <<  mySimVerts[track.vertIndex()].parentIndex() << endl; 
  }
  //////////////
// DT CHAMBERS //
  ////////////

  int  hitCount = 0;

  for (vector<PSimHit>::const_iterator iter=simDTHits.begin();iter<simDTHits.end();iter++){ //Iterator is from 0
    
    const PSimHit & hit = *iter;

    if(abs(hit.particleType())<1000000) continue;
      DTLayerId dtDetLayerId(hit.detUnitId()); 
      //cout << "======================================HIT SUMMARY==========================================" << std::endl;
      
      cout << dtDetLayerId << endl;
      int  eventNr = theEventCount;
	    int  pid = hit.particleType();
      int  trackNr = hit.trackId(); 
      cout << "TRACK ID: " << trackNr << endl;
      Local3DPoint localPosition = hit.localPosition();
      GlobalPoint globalPosition = globalGeometry.idToDet(dtDetLayerId)->toGlobal(hit.localPosition());
      double r = sqrt((globalPosition.x()*globalPosition.x())+(globalPosition.y()*globalPosition.y()));
      GlobalVector p = globalGeometry.idToDet(dtDetLayerId)->toGlobal(hit.momentumAtEntry());
      double pt = sqrt((p.y()*p.y())+(p.z()*p.z()));
      double phi = globalPosition.phi();
      if(phi<0){phi = phi + 2*M_PI;}
      double theta = globalPosition.theta();
      double eta = -log(tan(abs(theta)/2));
      double tof = hit.timeOfFlight();
      double distanceIP = globalPosition.mag();

      cout << "DEBUG1" << endl;
      int vtxIndex = mySimTracks[trackNr].vertIndex();
      if(vtxIndex<0) continue;
      math::XYZTLorentzVector hscpVertexV = mySimVerts[vtxIndex].position();
      cout << "DEBUG2" << endl;

	   double betaHit = p.mag()/sqrt((p.mag()*p.mag())+(hscpMass*hscpMass)); //.pabs() calculates the length of the vector pointing to the hit from (0,0,0) 
    															//mass and momenta in GeV (no unit conversion required)
      hitData.event = eventNr;
      hitData.track = trackNr;
      hitData.pid = pid;
      hitData.superlayer = dtDetLayerId.superlayer();
      hitData.station = dtDetLayerId.station();
      hitData.wheel = dtDetLayerId.wheel();
      hitData.sector = dtDetLayerId.sector();
      hitData.layer = dtDetLayerId.layer();
      hitData.localX = localPosition.x();
      hitData.localY = localPosition.y();
      hitData.localZ = localPosition.z();
      hitData.globalX = globalPosition.x();
      hitData.globalY = globalPosition.y();
      hitData.globalZ = globalPosition.z();
      hitData.pt = pt;
      hitData.phi = phi;
      hitData.r = r;
      hitData.eta = eta;
      hitData.tof = tof;
      hitData.beta = betaHit;
      hitData.distanceVTX = sqrt((globalPosition.x() - hscpVertexV.Px())*(globalPosition.x() - hscpVertexV.Px())+(globalPosition.y() - hscpVertexV.Py())*(globalPosition.y() - hscpVertexV.Py())+(globalPosition.z() - hscpVertexV.Pz())*(globalPosition.z() - hscpVertexV.Pz()));
      hitData.distanceIP = distanceIP;
      hitData.VTX_x = hscpVertexV.Px();
      hitData.VTX_y = hscpVertexV.Py();
      hitData.VTX_z = hscpVertexV.Pz();
      hscpTree->Fill();

  //////////////////
// TOF BETWEEN HITS // 
  //////////////////

        vector<PSimHit>::const_iterator next = std::next(iter,1);
        if(next==simDTHits.end()) break;
        const PSimHit & nextHit = *next;
        DTLayerId nextDTDetLayerId(nextHit.detUnitId());

        Local3DPoint nextHitLocal = nextHit.localPosition();
        GlobalPoint nextHitGlobal = globalGeometry.idToDet(nextDTDetLayerId)->toGlobal(nextHitLocal);  
        cout << "this hit position" << globalPosition << "|nextHitGlobal" << nextHitGlobal << endl << endl;
        Vector3DBase lHitDistV = nextHitLocal - localPosition;
        Vector3DBase gHitDistV = nextHitGlobal - globalPosition;

        double gHitDist = sqrt(gHitDistV.x()*gHitDistV.x() + gHitDistV.y()*gHitDistV.y() + gHitDistV.z()*gHitDistV.z());
        double lHitDist = sqrt(lHitDistV.x()*lHitDistV.x() + lHitDistV.y()*lHitDistV.y() + lHitDistV.z()*lHitDistV.z());
        cout << lHitDist << endl;
      
        double betaNextHit = nextHit.pabs()/sqrt((nextHit.pabs()*nextHit.pabs())+(hscpMass*hscpMass));

        double tofGivenDifference = nextHit.timeOfFlight() - tof;
        double tofBetweenHitsCurrent = gHitDist*1e9*0.01/((betaHit)*TMath::C());
        double tofBetweenHitsNext = gHitDist*1e9*0.01/((betaNextHit)*TMath::C());

        if(hit.trackId() != nextHit.trackId()) continue;

        nlines++;
      
      hitCompData.event = eventNr;
      hitCompData.track = trackNr;
      hitCompData.distanceVTX = sqrt((globalPosition.x() - hscpVertexV.Px())*(globalPosition.x() - hscpVertexV.Px())+(globalPosition.y() - hscpVertexV.Py())*(globalPosition.y() - hscpVertexV.Py())+(globalPosition.z() - hscpVertexV.Pz())*(globalPosition.z() - hscpVertexV.Pz()));
      hitCompData.nextLocalX = nextHitLocal.x();
      hitCompData.nextLocalY = nextHitLocal.y();
      hitCompData.nextLocalZ = nextHitLocal.z();
      hitCompData.nextGlobalX = nextHitGlobal.x();
      hitCompData.nextGlobalY = nextHitGlobal.y();
      hitCompData.nextGlobalZ = nextHitGlobal.z();
      hitCompData.globalHitDistDiff = gHitDist;
      hitCompData.betaNextHit = betaNextHit;
      hitCompData.diffTOFGeant = tofGivenDifference;
      hitCompData.diffTOFHitB = tofBetweenHitsCurrent;
      hitCompData.diffTOFNextB = tofBetweenHitsNext;

      bool sect = false;
      bool supl = false;
      bool wh = false;
      bool st = false;

      if(nextDTDetLayerId.wheel() == dtDetLayerId.wheel()) wh = true;
      if(nextDTDetLayerId.superlayer() == dtDetLayerId.superlayer()) supl = true;
      if(nextDTDetLayerId.sector() == dtDetLayerId.sector()) sect = true;
      if(nextDTDetLayerId.station() == dtDetLayerId.station()) st = true;

      hitCompData.samesect = sect;
      hitCompData.samesl = supl;
      hitCompData.samewheel = wh;
      hitCompData.samest = st;

      hscpCompTree->Fill();
      hitCount++; //has to be here in order for hit count to match numbering in hit ntuples :)
    
  }
  //////////////////////////
// Simulated track analysis //
  //////////////////////////

  for (std::vector<SimTrack>::const_iterator it=mySimTracks.begin(); it<mySimTracks.end(); it++) {
    const SimTrack & track = *it;
    if ( track.type() == -99) continue;
    if ( track.vertIndex() != 0) continue;

    double phi_sim = track.momentum().phi(); //momentum azimutal angle
    double pt_sim = track.momentum().pt(); //transverse momentum
    double eta_sim = track.momentum().eta(); //pseudorapidity

    if(track.type()>1000000) hscpCount++;

  }


}
DEFINE_FWK_MODULE(Projekt);

