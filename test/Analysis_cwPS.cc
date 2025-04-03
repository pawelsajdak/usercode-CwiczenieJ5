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
#include "DataFormats/TrackReco/interface/Track.h"

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
double BpmMass = 5.27941;

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
  TH1D* htrack_dR;
  //TH1D *hKaonKaon,*hPionPion,*hKaonPion;

  edm::EDGetTokenT< vector<pat::Muon> > theMuonToken;
  edm::EDGetTokenT< vector<pat::PackedCandidate> > theCandidateToken;
  edm::EDGetTokenT< vector<reco::Vertex> > thePrimaryVertexToken;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTrackBuilderToken;
};


Analysis::Analysis(const edm::ParameterSet& conf)
  : theConfig(conf), debug(false),  theEventCount(0)
{
  cout <<" CTORXX" << endl;
  theMuonToken = consumes< vector<pat::Muon> >( theConfig.getParameter<edm::InputTag>("muonSrc"));
  theCandidateToken     = consumes< vector<pat::PackedCandidate> > (edm::InputTag("packedPFCandidates"));
  thePrimaryVertexToken = consumes< vector<reco::Vertex> > (edm::InputTag("offlineSlimmedPrimaryVertices"));
  theTrackBuilderToken = esConsumes(edm::ESInputTag("", "TransientTrackBuilder"));
  if(theConfig.exists("debug")) debug = theConfig.getParameter<bool>("debug"); 
}

Analysis::~Analysis()
{
  cout <<" DTOR" << endl;
}

void Analysis::beginJob()
{
  htrack_dR = new TH1D("htrack_dR","track dR",1000,0.0,0.1);
  /*/create a histogram
  hKaonKaon = new TH1D("hKaonKaon","K+K- from Jpsi vertex;Minv;Counts",10000,0.,15.);
  hPionPion = new TH1D("hPionPion","#pi+#pi- from Jpsi vertex;Minv;Counts",10000,0.,15.);
  hKaonPion = new TH1D("hKaonPion","K#pm#pi#pm (opposite signs) from Jpsi vertex;Minv;Counts",10000,0.,15.);
  */
  cout << "HERE Analysis::beginJob()" << endl;
}

void Analysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");
  //write histogram data
  //hKaonKaon->Write();
  //hPionPion->Write();
  //hKaonPion->Write();

  htrack_dR->Write();
  myRootFile.Close();
  delete htrack_dR;
  /*delete hKaonKaon;
  delete hPionPion;
  delete hKaonPion;
  */
  cout << "HERE Cwiczenie::endJob()" << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////
void Analysis::analyze(const edm::Event& ev, const edm::EventSetup& es)
{
  if (debug) std::cout << " -------------------------------- HERE Cwiczenie::analyze "<< std::endl;
  const vector<pat::Muon> & muons = ev.get(theMuonToken);
  const vector<pat::PackedCandidate> & candidates = ev.get(theCandidateToken);
  const vector <reco::Vertex> & primVertices = ev.get(thePrimaryVertexToken);
  const auto & trackBuilder = es.getData(theTrackBuilderToken);

  if (debug) std::cout <<" number of      muons: " << muons.size() <<std::endl;
 
        for(std::vector<reco::Vertex>::const_iterator ipv1 = primVertices.begin();ipv1<primVertices.end();ipv1++)
        {
          cout << "New Primary Vertex: ("<<ipv1->x()<<", "<<ipv1->y()<<", "<<ipv1->z()<<")" <<endl;
          ROOT::Math::XYZPoint testPoint(1.0,1.0,5.0);
          ROOT::Math::XYZVector bpmDispl (testPoint-ipv1->position());
          cout << "Displacement to (1.0,1.0,5.0): ("<<bpmDispl.x()<<", "<<bpmDispl.y()<<", "<<bpmDispl.z()<<")" <<endl;

          //if(fabs(ipv1->z()-vBX.z())> 0.2) continue;  //ct~0.05 cm
          //cout << "I survived"<<endl;
          
          /*
          for(reco::Vertex::trackRef_iterator itr1 = ipv1->tracks_begin();itr1<ipv1->tracks_end();itr1++) //loop over tracks of a primary vertex
          {
            reco::TrackBaseRef track = *itr1;
            const math::XYZVector & trackMom = track->momentum();
            std::cout << trackMom.mag2()<<std::endl;
            double dR = deltaR(trackMom,BpmMom);
            cout << "dR: "<<dR<<endl;
            if(dR<0.1) htrack_dR->Fill(dR);
          }
          */
        }       
    
  
  cout << "\n";

  if (debug) cout <<"*** Analyze event: " << ev.id()<<" analysed event count:"<<++theEventCount << endl;
}

DEFINE_FWK_MODULE(Analysis);
