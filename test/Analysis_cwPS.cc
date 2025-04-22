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
double psi2SMass = 3.686097;
double B0sMass = 5.3669;

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
  TH1D *hKaonKaon;

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
  hKaonKaon = new TH1D("hKaonKaon","K+K- from Jpsi vertex;Minv;Counts",10000,0.,15.);
  
  cout << "HERE Analysis::beginJob()" << endl;
}

void Analysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");
  //write histogram data
  hKaonKaon->Write();
  

  myRootFile.Close();
  delete hKaonKaon;
  
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

      // rescale muom momenta for exact jpsi mass
      double alpha=1.;
      math::XYZVector mom1 = im1->momentum();
      math::XYZVector mom2 = im2->momentum();
      {
        double a = mom1.mag2()*mom2.mag2()-sqr(mom1.Dot(mom2));
        double b = -sqr(jpsiMass)*mom1.Dot(mom2)+sqr(muonMass)*(mom1+mom2).mag2();
        double c = -sqr(jpsiMass)*(sqr(jpsiMass)/4.-sqr(muonMass));
        double delta= sqr(b)-4*a*c;
        alpha = sqrt((-b+sqrt(delta))/2./a);
      } 
      lMuonsVector = lorentzVector((mom1+mom2)*alpha, jpsiMass);

      


      /////////////FIRST PACKED CANDIDATE///////////
      for (std::vector<pat::PackedCandidate>::const_iterator ic1 = candidates.begin(); ic1 < candidates.end(); ic1++) 
      {
        if(abs(ic1->pdgId()) != 211 || !ic1->hasTrackDetails() || ic1->pt() < 2. || ic1->charge()==0) continue;
        
        // Could J/psi and the candidate (kaon,pion,proton) come from a common vertex - vBX?
        const reco::Track & trk1 = ic1->pseudoTrack();
        if (fabs(vjp.position().z()- trk1.vz())>0.3)continue;
        trackTTs.push_back(trackBuilder.build(trk1));
        reco::Vertex vBX(TransientVertex(kvf.vertex(trackTTs)));
        double probvBX = TMath::Prob(vBX.chi2(),vBX.ndof());
        trackTTs.pop_back();  
        if (probvBX<0.15) continue;

        // deltaR check
        if(std::min(deltaR(trk1,*mu1Ref),deltaR(trk1,*mu2Ref))<0.0003) continue;

        ///////////SECOND PACKED CANDIDATE///////////////
        trackTTs.push_back(trackBuilder.build(trk1)); // was removed, now added again
        math::XYZVector cand1Mom = ic1->momentum();
        for (std::vector<pat::PackedCandidate>::const_iterator ic2 = ic1+1; ic2 < candidates.end(); ic2++) 
        {
          if(abs(ic2->pdgId()) != 211 || !ic2->hasTrackDetails() || ic2->pt() < 2. || ic2->charge()*ic1->charge() !=-1) continue;

          // Could J/psi and both candidates come from a common vertex - vJXX?
          const reco::Track & trk2 = ic2->pseudoTrack();
          if (fabs(vBX.position().z()- trk2.vz())>0.3)continue;
          
          trackTTs.push_back(trackBuilder.build(trk2));
          reco::Vertex vJXX(TransientVertex(kvf.vertex(trackTTs)));
          double probvJXX = TMath::Prob(vJXX.chi2(),vJXX.ndof());
          trackTTs.pop_back();  
          if (probvJXX<0.15) continue;

          // deltaR check
          if(std::min(deltaR(trk2,*mu1Ref),deltaR(trk2,*mu2Ref))<0.0003) continue;

          // HISTOGRAMS
          math::XYZVector cand2Mom = ic2->momentum();

          // two Kaons
          ROOT::Math::PxPyPzEVector lVectorKK = lorentzVector(cand1Mom, kaonMass)+lorentzVector(cand2Mom,kaonMass);
          ROOT::Math::PxPyPzEVector lVectorJKK = lMuonsVector + lVectorKK;
          if(fabs(lVectorJKK.M()-B0sMass)>0.12) continue;
          hKaonKaon->Fill(lVectorKK.M());

        }
        trackTTs.pop_back();  // removes trk1
        ///////////////////////

      }  
    }
  } 
    
   
  
  cout << "\n";


  if (debug) cout <<"*** Analyze event: " << ev.id()<<" analysed event count:"<<++theEventCount << endl;
}

DEFINE_FWK_MODULE(Analysis);