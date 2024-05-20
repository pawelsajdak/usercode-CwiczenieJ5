#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
// #include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include <Math/Vector4D.h>

#include <sstream>


using namespace std;


double muonMass = 0.105658;
double kaonMass = 0.493677;
double jpsiMass = 3.096916;
ROOT::Math::PxPyPzEVector lorentzVector(const math::XYZVector & mom, double mass) {
  return ROOT::Math::PxPyPzEVector( mom.x(), mom.y(), mom.z(), sqrt( mass*mass+mom.mag2()));
}

//object definition
class Analysis : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit Analysis(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~Analysis();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;
  TH1D *histo;
  TH1D *hMjpsi, *hMBpm;

  edm::EDGetTokenT< vector<pat::Muon> > theMuonToken;
  edm::EDGetTokenT< vector<pat::Muon> > theMuonDsplToken;
  edm::EDGetTokenT< vector<reco::Track> > theDsplTrkToken;
  edm::EDGetTokenT< vector<reco::Track> > theDsplMuToken;
  edm::EDGetTokenT< vector<reco::Vertex> > theVertexToken, theVertexWithBSToken;
  edm::EDGetTokenT< vector<reco::VertexCompositePtrCandidate> > theVertexCPCToken; 
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTrackBuilderToken;
  edm::EDGetTokenT< vector<pat::PackedCandidate> > theCandidateToken;

};


Analysis::Analysis(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;
  theMuonToken = consumes< vector<pat::Muon> >( theConfig.getParameter<edm::InputTag>("muonSrc"));
  theMuonDsplToken = consumes< vector<pat::Muon> >( edm::InputTag("slimmedDisplacedMuons"));
  theDsplMuToken = consumes< vector<reco::Track> >( edm::InputTag("displacedGlobalMuons"));
  theDsplTrkToken = consumes< vector<reco::Track> >( edm::InputTag("displacedTracks"));
  theVertexToken = consumes< vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"));
  theVertexWithBSToken = consumes< vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVerticesWithBS"));
  theVertexCPCToken = consumes< vector<reco::VertexCompositePtrCandidate> >(edm::InputTag("slimmedSecondaryVertices"));
  theCandidateToken     = consumes< vector<pat::PackedCandidate> > (edm::InputTag("packedPFCandidates"));
  theTrackBuilderToken = esConsumes(edm::ESInputTag("", "TransientTrackBuilder"));
}

Analysis::~Analysis()
{
  cout <<" DTOR" << endl;
}

void Analysis::beginJob()
{
  //create a histogram
  histo =new TH1D("histo","test; X; #events",10, 0., 10.);
  hMjpsi=new TH1D("hMjpsi","test; X; #events",1000, 2., 12.);
  hMBpm=new TH1D("hMBpm","test; X; #events",200, 4., 6.);
  cout << "HERE Analysis::beginJob()" << endl;
}

void Analysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");
  //write histogram data
  histo->Write();
  hMjpsi->Write();
  hMBpm->Write();
  myRootFile.Close();
  delete histo;
  cout << "HERE Cwiczenie::endJob()" << endl;
}

void Analysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  bool debug = false;
  if (debug) std::cout << " -------------------------------- HERE Cwiczenie::analyze "<< std::endl;
  const vector<pat::Muon> & muons = ev.get(theMuonToken);
  const vector<pat::PackedCandidate> & candidates = ev.get(theCandidateToken);

  if (debug) std::cout <<" number of      muons: " << muons.size() <<std::endl;
  const auto & trackBuilder = es.getData(theTrackBuilderToken);

  std::vector< std::pair<reco::TransientTrack, reco::TransientTrack> > jpsis;
  for (std::vector<pat::Muon>::const_iterator im1 = muons.begin(); im1 < muons.end(); im1++) {
    const pat::Muon & muon = *im1;
    if(!im1->isGlobalMuon() || !im1->isTrackerMuon()) continue; // skip if mu1 is not GLB or TRK (PASS IF LooseMuId)
    if (im1->pt() < 3.) continue;
    ROOT::Math::PxPyPzEVector lm1(lorentzVector(im1->momentum(),muonMass));
    histo->Fill(muon.eta());
    if (debug) {
      std::cout <<"pat muon pt: "<<muon.charge()*muon.pt()<<" mom: "<<muon.momentum()
                <<" vertex: "<<muon.vertex(); 
      const reco::Track * tr=muon.bestTrack(); if(tr!=nullptr) std::cout <<" dxy: "<<tr->dxy();
      std::cout <<std::endl;
    }
    //reco::TrackRef mu1Ref = im1->get<reco::TrackRef>();
    reco::TrackRef mu1Ref = im1->track();
    if (!mu1Ref)continue;
    for (std::vector<pat::Muon>::const_iterator im2 = im1+1; im2 < muons.end(); im2++) {
      if(!im2->isGlobalMuon() || !im2->isTrackerMuon()) continue; // skip if mu1 is not GLB or TRK (PASS IF LooseMuId)
      //reco::TrackRef mu2Ref = im1->get<reco::TrackRef>();
      reco::TrackRef mu2Ref = im2->track();
      if (!mu2Ref)continue;
      if (im1->charge()*im2->charge() != -1) continue;
      if(im2->pt() < 3.) continue;
      if (fabs(im1->vz()-im2->vz())> 0.3) continue;
      std::vector<reco::TransientTrack> trackTTs;
      trackTTs.push_back(trackBuilder.build(mu1Ref));
      trackTTs.push_back(trackBuilder.build(mu2Ref));
      KalmanVertexFitter kvf(true);
      reco::Vertex vjp(TransientVertex(kvf.vertex(trackTTs)));
      double prob = TMath::Prob(vjp.chi2(),vjp.ndof());
      ROOT::Math::PxPyPzEVector lm2(lorentzVector(im2->momentum(), muonMass));
      ROOT::Math::PxPyPzEVector ljp=lm1+lm2;
      if(debug) { 
        std::cout <<" fitted vertex: "<<vjp.position()<<" prob.: "<<prob<<std::endl;
        std::cout <<" lm1: "<<lm1.E()<<" "<<lm1.mass()<<std::endl;
        std::cout <<" lm2: "<<lm2.E()<<" "<<lm2.mass()<<std::endl;
        std::cout <<" ljp: "<<ljp.E()<<" "<<ljp.mass()<<std::endl;
      }
      if (prob>0.1) hMjpsi->Fill(ljp.mass());
      if (fabs(ljp.mass()-jpsiMass)>0.06) continue;
      for (std::vector<pat::PackedCandidate>::const_iterator ic1 = candidates.begin(); ic1 < candidates.end(); ic1++) {
        if( abs(ic1->pdgId()) != 211 || !ic1->hasTrackDetails() || ic1->pt() < 1. || ic1->charge()==0 ) continue;
        //reco::TrackRef mu1Ref = im1->get<reco::TrackRef>();
        const reco::Track & trk1 = ic1->pseudoTrack();
        if (fabs(vjp.position().z()- trk1.vz())>0.3)continue;
        trackTTs.push_back(trackBuilder.build(trk1));
        reco::Vertex vBpm(TransientVertex(kvf.vertex(trackTTs)));
        double probBpm = TMath::Prob(vBpm.chi2(),vBpm.ndof());  
        ROOT::Math::PxPyPzEVector lc1(lorentzVector(ic1->momentum(), kaonMass));
        ROOT::Math::PxPyPzEVector lBpm = ljp+lc1;
        if (debug) std::cout<<" candidate: prob:"<<probBpm<<" mass: "<<lBpm.mass()<<std::endl;
        if (probBpm>0.1) hMBpm->Fill(lBpm.mass());
      }
    }
  }
  /*
  const vector<pat::Muon> & muonsDspl = ev.get(theMuonDsplToken);
  std::cout <<" number of Dspl muons: " << muonsDspl.size() <<std::endl;
  for (const auto & muon : muonsDspl) {
    if (debug) {
      std::cout <<"pat muon pt: "<<muon.charge()*muon.pt()<<" mom: "<<muon.momentum()
                <<" vertex: "<<muon.vertex(); 
      const reco::Track * tr=muon.bestTrack(); if(tr!=nullptr) std::cout <<" dxy: "<<tr->dxy();
      std::cout <<std::endl;
    }
  }
  */

/*
  const vector<reco::Track> & dsplTrks = ev.get(theDsplTrkToken);
  std::cout <<": number of dspl trks: "<<dsplTrks.size() << std::endl;
  for (const auto & trk : dsplTrks) {
    std::cout <<"trk pt: "<<trk.charge()*trk.pt()<<" mom: "<<trk.momentum()
              <<" refPoint; "<<trk.referencePoint()<< " dxy: "<<trk.dxy()<< std::endl;
  }
  const vector<reco::Track> & dsplGlMu = ev.get(theDsplMuToken);
  std::cout <<": number of dspl GlMu: "<<dsplGlMu.size() << std::endl;
  for (const auto & trk : dsplGlMu) {
    std::cout <<"trk pt: "<<trk.charge()*trk.pt()<<" mom: "<<trk.momentum()
              <<" refPoint; "<<trk.referencePoint()<< " dxy: "<<trk.dxy()<< std::endl;
  }
  */
  
 /* 
  const vector<reco::Vertex> & vtcs = ev.get(theVertexToken);
  if (debug) std::cout<<" number of vertices:         " << vtcs.size() << std::endl;
  for (const auto & vtx : vtcs) { std::cout <<" vertex:"<< vtx.position()<<std::endl; }
  */

//  const vector<reco::Vertex> & vtcsBS = ev.get(theVertexWithBSToken);
// if (debug) std::cout<<" number of vertices with BS: " << vtcsBS.size() << std::endl;
//  for (const auto & vtx : vtcsBS) { std::cout <<" vertex:"<< vtx.position()<<std::endl; }

/*
  const vector<reco::VertexCompositePtrCandidate> & vCPC = ev.get(theVertexCPCToken);
  std::cout <<"number of vCPC : "<<vCPC.size() <<std::endl;
  for (const auto & vtx : vCPC) {
    std::cout << " position: "<< vtx.vertex()<< " id: "<<vtx.pdgId() << std::endl;
    std::cout << " number of mothers:  "<<vtx.numberOfMothers()<<std::endl;
    std::cout << " number of dauthers: "<<vtx.numberOfDaughters()<< std::endl;
    for (unsigned int i=0; i < vtx.numberOfDaughters(); ++i) {
      const reco::Candidate * d = vtx.daughter(i);
      std::cout <<"pdgid: "<<d->pdgId()<<" mom: "<< d->momentum()<<" vtx: "<<d->vertex()<<std::endl; 
    }
  }
  */

  //write std io
  if (debug) cout <<"*** Analyze event: " << ev.id()<<" analysed event count:"<<++theEventCount << endl;
}

DEFINE_FWK_MODULE(Analysis);

