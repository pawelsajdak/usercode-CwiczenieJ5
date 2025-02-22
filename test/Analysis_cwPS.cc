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
double jpsiMass = 3.096900;
double pionMass = 0.139570;
double protonMass = 0.938272;
double lambdaMass = 1.115683;
double phiMass = 1.019461;

template <typename T> T sqr(T v) { return v*v; }

ROOT::Math::PxPyPzEVector lorentzVector(const math::XYZVector & mom, double mass) {
  return ROOT::Math::PxPyPzEVector( mom.x(), mom.y(), mom.z(), sqrt( sqr(mass)+mom.mag2()));
}
ROOT::Math::PxPyPzEVector lorentzVector(const ROOT::Math::PxPyPzEVector & orig, double mass) {
  return ROOT::Math::PxPyPzEVector(orig).SetE(sqrt(mass*mass+orig.P2()));
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
  bool debug;
  unsigned int theEventCount;
  TH1D *histoK, *histoPi, *histoPr;

  edm::EDGetTokenT< vector<pat::Muon> > theMuonToken;
  edm::EDGetTokenT< vector<pat::PackedCandidate> > theCandidateToken;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTrackBuilderToken;
};


Analysis::Analysis(const edm::ParameterSet& conf)
  : theConfig(conf), debug(false),  theEventCount(0)
{
  cout <<" CTORXX" << endl;
  theMuonToken = consumes< vector<pat::Muon> >( theConfig.getParameter<edm::InputTag>("muonSrc"));
  theCandidateToken     = consumes< vector<pat::PackedCandidate> > (edm::InputTag("packedPFCandidates"));
  theTrackBuilderToken = esConsumes(edm::ESInputTag("", "TransientTrackBuilder"));
  if(theConfig.exists("debug")) debug = theConfig.getParameter<bool>("debug"); 
}

Analysis::~Analysis()
{
  cout <<" DTOR" << endl;
}

void Analysis::beginJob()
{
  //create a histogram
  histoK =new TH1D("histoK","kaon; Minv; #events",2200, 3.8,6.);
  histoPi =new TH1D("histoPi","pion; Minv; #events",2200, 3.8, 6.);
  histoPr =new TH1D("histoPr","proton; Minv; #events",2200, 3.8, 6.);
  cout << "HERE Analysis::beginJob()" << endl;
}

void Analysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");
  //write histogram data
  histoK->Write();
  cout << "Wrote histoK \n";
  histoPi->Write();
  cout << "Wrote histoPi \n";
  histoPr->Write();
  cout << "Wrote histoPr \n";
  myRootFile.Close();
  delete histoK;
  delete histoPi;
  delete histoPr;
  cout << "HERE Cwiczenie::endJob()" << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::analyze(const edm::Event& ev, const edm::EventSetup& es)
{
  if (debug) std::cout << " -------------------------------- HERE Cwiczenie::analyze "<< std::endl;
  const vector<pat::Muon> & muons = ev.get(theMuonToken);
  const vector<pat::PackedCandidate> & candidates = ev.get(theCandidateToken);
  const auto & trackBuilder = es.getData(theTrackBuilderToken);

  if (debug) std::cout <<" number of      muons: " << muons.size() <<std::endl;
 
  //std::vector< std::pair<reco::TransientTrack, reco::TransientTrack> > jpsis;
  for (std::vector<pat::Muon>::const_iterator im1 = muons.begin(); im1 < muons.end(); im1++)
  {
    const pat::Muon & muon = *im1;
    if(!im1->isGlobalMuon() || !im1->isTrackerMuon()) continue;
    if(muon.pt()<3) continue;
    reco::TrackRef mu1Ref = im1->track();
    if (!mu1Ref)continue;

    for (std::vector<pat::Muon>::const_iterator im2 = im1+1; im2 < muons.end(); im2++)
    {
      if(!im2->isGlobalMuon() || !im2->isTrackerMuon()) continue;
      const pat::Muon & muon2 = *im2;
      if(muon2.pt()<3 || muon.charge()*muon2.charge()!=-1) continue;
      reco::TrackRef mu2Ref = im2->track();
      if (!mu2Ref)continue;
      if(fabs(muon.vz()-muon2.vz())>0.3) continue;

      ROOT::Math::PxPyPzEVector lMuonsVector = muon.p4()+muon2.p4();
      //Minv of two muons close to the J/psi peak
      if(fabs(lMuonsVector.M()-jpsiMass)>0.1) continue;

      // Could the two muons have a common vertex - vjp?
      std::vector<reco::TransientTrack> trackTTs;
      trackTTs.push_back(trackBuilder.build(mu1Ref));
      trackTTs.push_back(trackBuilder.build(mu2Ref));
      KalmanVertexFitter kvf(true);
      reco::Vertex vjp(TransientVertex(kvf.vertex(trackTTs)));
      double prob = TMath::Prob(vjp.chi2(),vjp.ndof());
      if (prob<0.1) continue;



      for (std::vector<pat::PackedCandidate>::const_iterator ic1 = candidates.begin(); ic1 < candidates.end(); ic1++) 
      {
        if(abs(ic1->pdgId()) != 211 || !ic1->hasTrackDetails() || ic1->pt() < 2. || ic1->charge()==0) continue;
        
        // Could J/psi and the candidate (kaon,pion,proton) come from a common vertex - vBX?
        const reco::Track & trk1 = ic1->pseudoTrack();
        if (fabs(vjp.position().z()- trk1.vz())>0.3)continue;
        /////
        //cout<< "before adding the candidate "<<trackTTs.size()<<endl;
        trackTTs.push_back(trackBuilder.build(trk1));
        //cout<< "after adding the candidate "<<trackTTs.size()<<endl;
        reco::Vertex vBX(TransientVertex(kvf.vertex(trackTTs)));
        double probvBX = TMath::Prob(vBX.chi2(),vBX.ndof());
        trackTTs.pop_back();  
        if (probvBX<0.15) continue;
        
        // Kaon
        math::XYZVector candMom = ic1->momentum();
        ROOT::Math::PxPyPzEVector lFullVectorK = lMuonsVector+lorentzVector(candMom, kaonMass);
        histoK->Fill(lFullVectorK.M());

        // Pion
        ROOT::Math::PxPyPzEVector lFullVectorPi = lMuonsVector+lorentzVector(candMom, pionMass);
        histoPi->Fill(lFullVectorPi.M());

        // Proton
        ROOT::Math::PxPyPzEVector lFullVectorPr = lMuonsVector+lorentzVector(candMom, protonMass);
        histoPr->Fill(lFullVectorPr.M());

      }  
    }
  } 
    
   
  
  cout << "\n";


  if (debug) cout <<"*** Analyze event: " << ev.id()<<" analysed event count:"<<++theEventCount << endl;
}

DEFINE_FWK_MODULE(Analysis);
