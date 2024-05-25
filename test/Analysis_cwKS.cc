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

ROOT::Math::PxPyPzEVector lorentzVector(const math::XYZVector & mom, double mass) {
  return ROOT::Math::PxPyPzEVector( mom.x(), mom.y(), mom.z(), sqrt( mass*mass+mom.mag2()));
}
ROOT::Math::PxPyPzEVector lorentzVector(const ROOT::Math::PxPyPzEVector & orig, double mass) {
  return ROOT::Math::PxPyPzEVector(orig).SetE(sqrt(mass*mass+orig.P2()));
}



template <typename T> T sqr(T v) { return v*v; }
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
  TH1D *histo;
  TH1D *hmJpsi, *hmBpm, *hmBpm_mm, *hmLam, *hTest1, *hTest2, *hTest3, *hTest4A, *hTest4B;
  TH1D *hmKK, *hmBs, *hmLamB;

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
  : theConfig(conf), debug(false),  theEventCount(0)
{
  cout <<" CTORXX" << endl;
  theMuonToken = consumes< vector<pat::Muon> >( theConfig.getParameter<edm::InputTag>("muonSrc"));
  theMuonDsplToken = consumes< vector<pat::Muon> >( edm::InputTag("slimmedDisplacedMuons"));
  theDsplMuToken = consumes< vector<reco::Track> >( edm::InputTag("displacedGlobalMuons"));
  theDsplTrkToken = consumes< vector<reco::Track> >( edm::InputTag("displacedTracks"));
  theVertexToken = consumes< vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"));
  theVertexWithBSToken = consumes< vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVerticesWithBS"));
  theVertexCPCToken = consumes< vector<reco::VertexCompositePtrCandidate> >(edm::InputTag("slimmedLambdaVertices"));
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
  histo =new TH1D("histo","test; X; #events",10, 0., 10.);
  hmJpsi=new TH1D("hmJpsi","test; X; #events",1000, 2., 12.);
  hmKK=new TH1D("hmKK","test; X; #events",1000, 0., 7.);
  hmBpm=new TH1D("hmBpm","test; X; #events",300, 4., 7.);
  hmBpm_mm=new TH1D("hmBpm_mm","test; X; #events",300, 4., 7.);
  hmBs=new TH1D("hmBs","test; X; #events",300, 4., 7.);
  hmLam=new TH1D("hmLam","test; X; #events",100, 0.9, 1.9);
  hmLamB=new TH1D("hmLamB","test; X; #events",100, 4., 7);
  hTest1=new TH1D("hTest1","test; X; #events",300, 4., 7.);
  hTest2=new TH1D("hTest2","test; X; #events",300, 4., 7.);
  hTest3=new TH1D("hTest3","test; X; #events",300, 4., 10.);
  hTest4A=new TH1D("hTest4A","test; X; #events",300, 4., 10.);
  hTest4B=new TH1D("hTest4B","test; X; #events",300, 4., 10.);
  cout << "HERE Analysis::beginJob()" << endl;
}

void Analysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");
  //write histogram data
  histo->Write();
  hmJpsi->Write();
  hmBpm->Write();
  hmBpm_mm->Write();
  hmLam->Write();
  hmLamB->Write();
  hTest1->Write();
  hTest2->Write();
  hTest3->Write();
  hTest4A->Write();
  hTest4B->Write();
  hmBs->Write();
  hmKK->Write();
  myRootFile.Close();
  delete histo;
  cout << "HERE Cwiczenie::endJob()" << endl;
}

void Analysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
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
      ROOT::Math::PxPyPzEVector lmm=lm1+lm2;
      if(debug) { 
        std::cout <<" fitted vertex: "<<vjp.position()<<" prob.: "<<prob<<std::endl;
        std::cout <<" lm1: "<<lm1.E()<<" "<<lm1.mass()<<std::endl;
        std::cout <<" lm2: "<<lm2.E()<<" "<<lm2.mass()<<std::endl;
        std::cout <<" lmm: "<<lmm.E()<<" "<<lmm.mass()<<std::endl;
      }
      if (prob<0.1) continue;
      hmJpsi->Fill(lmm.mass());
      if (fabs(lmm.mass()-jpsiMass)>0.15) continue;

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
      if (debug) {
        ROOT::Math::PxPyPzEVector lm1p(lorentzVector(mom1*alpha, muonMass));
        ROOT::Math::PxPyPzEVector lm2p(lorentzVector(mom2*alpha, muonMass));
        std::cout <<" mumu mass: "<<lmm.mass()<<" rescaled mumu wrt jpsi"<< (lm1p+lm2p).mass()-jpsiMass<<" alpha: "<<alpha<<std::endl;
      }  
      ROOT::Math::PxPyPzEVector ljp(lorentzVector((mom1+mom2)*alpha, jpsiMass));

      if (ljp.pt()<8) continue;

      bool hasLambda = false;
      ROOT::Math::PxPyPzEVector lLam_pi;
      ROOT::Math::PxPyPzEVector lLam_pr;
      ROOT::Math::PxPyPzEVector lLam;
      double diff = 1.;
      const vector<reco::VertexCompositePtrCandidate> & vCPC = ev.get(theVertexCPCToken);
      for (const auto & vtx : vCPC) {
//        if (fabs(vtx.vertex().z()-vjp.vz()) > 0.3) continue; 
        if(vtx.numberOfDaughters()!=2) continue; 
        const reco::Candidate * d1 = vtx.daughter(0);
        const reco::Candidate * d2 = vtx.daughter(1);
        ROOT::Math::PxPyPzEVector lpi_A(lorentzVector(d1->momentum(), pionMass));
        ROOT::Math::PxPyPzEVector lpi_B(lorentzVector(d2->momentum(), pionMass));
        ROOT::Math::PxPyPzEVector lpr_A(lorentzVector(d2->momentum(), protonMass));
        ROOT::Math::PxPyPzEVector lpr_B(lorentzVector(d1->momentum(), protonMass));
        std::vector<reco::TransientTrack> tracksLB = trackTTs;
        if (d1->bestTrack() != nullptr && d2->bestTrack() != nullptr) {
          tracksLB.push_back(trackBuilder.build(d1->bestTrack()));
          tracksLB.push_back(trackBuilder.build(d2->bestTrack()));
        } 
//        if (d2->bestTrack() != nullptr) tracksLB.push_back(trackBuilder(*(d2->bestTrack())));

        if (tracksLB.size()!=4) continue;
        double diffA = fabs((lpi_A+lpr_A).mass()- lambdaMass);
        if (diffA < 0.055 && diffA < diff) {
          reco::Vertex vLB(TransientVertex(kvf.vertex(tracksLB)));
          double prob = TMath::Prob(vLB.chi2(),vLB.ndof());
          if (prob >0.1) { 
            diff = diffA;
            lLam_pi = lpi_A;
            lLam_pr = lpr_A; 
            lLam = lpi_A+lpr_A;
            hasLambda = true;
          }
        }
        double diffB = fabs((lpi_B+lpr_B).mass()- lambdaMass);
        if (diffB < 0.055 && diffB < diff) {
          reco::Vertex vLB(TransientVertex(kvf.vertex(tracksLB)));
          double prob = TMath::Prob(vLB.chi2(),vLB.ndof());
          if (prob>0.1) {
            diff = diffB;
            lLam_pi = lpi_A;
            lLam_pr = lpr_A; 
            lLam = lpi_B+lpr_B;
            hasLambda = true;
          }
        }
        hmLam->Fill((lpi_A+lpr_A).mass());
        hmLam->Fill((lpi_B+lpr_B).mass());
        if (debug) std::cout <<" Lambda mass 1: "<<(lpi_A+lpr_A).mass()<<" mass2: "<<(lpi_B+lpr_B).mass()<<std::endl;
      }
      if (hasLambda) hmLamB->Fill((ljp+lLam).mass());

      for (std::vector<pat::PackedCandidate>::const_iterator ic1 = candidates.begin(); ic1 < candidates.end(); ic1++) {
        if( abs(ic1->pdgId()) != 211 || !ic1->hasTrackDetails() || ic1->pt() < 2. || ic1->charge()==0 ) continue;
        //reco::TrackRef mu1Ref = im1->get<reco::TrackRef>();
        const reco::Track & trk1 = ic1->pseudoTrack();
        if (fabs(vjp.position().z()- trk1.vz())>0.3)continue;
        trackTTs.push_back(trackBuilder.build(trk1));
        reco::Vertex vBpm(TransientVertex(kvf.vertex(trackTTs)));
        double probBpm = TMath::Prob(vBpm.chi2(),vBpm.ndof());  
        if (probBpm<0.15) continue;
        
        ROOT::Math::PxPyPzEVector lK1(lorentzVector(ic1->momentum(), kaonMass));
        ROOT::Math::PxPyPzEVector lBpm = ljp+lK1;
        ROOT::Math::PxPyPzEVector lpi1(lorentzVector(ic1->momentum(), pionMass));
        ROOT::Math::PxPyPzEVector lX=lpi1+ljp;

        if (lBpm.pt()>10) { 
          hmBpm->Fill(lBpm.mass());
          hmBpm_mm->Fill((lmm+lK1).mass());
          hTest1->Fill((ROOT::Math::PxPyPzEVector(lorentzVector(ic1->momentum(),protonMass))+ljp).mass());
          hTest2->Fill(lX.mass());
          if (hasLambda) hTest4A->Fill((lLam+lX).mass());
          if (hasLambda && fabs(lX.mass()-4.38)<0.04) hTest4B->Fill((lLam+lX).mass());
        }

        for (std::vector<pat::PackedCandidate>::const_iterator ic2 = ic1+1; ic2 < candidates.end(); ic2++) {
          if(abs(ic2->pdgId()) != 211 || !ic2->hasTrackDetails() || ic2->pt() < 1. || ic2->charge()==0 ) continue;
          if( ic1->charge()*ic2->charge() >= 0) continue;
          const reco::Track & trk2 = ic2->pseudoTrack();
          if (fabs(vjp.position().z()- trk2.vz())>0.3)continue;
          std::vector<reco::TransientTrack> tracks4=trackTTs;
          tracks4.push_back(trackBuilder.build(trk2));
          reco::Vertex vtx4(TransientVertex(kvf.vertex(tracks4)));
          if (TMath::Prob(vtx4.chi2(),vtx4.ndof())<0.1) continue;  
          // check for Phi and Bs
          ROOT::Math::PxPyPzEVector lK2(lorentzVector(trk2.momentum(), kaonMass));
          ROOT::Math::PxPyPzEVector lKK=lK1+lK2;
          hmKK->Fill((lK1+lK2).mass());
          if (fabs(lKK.mass()-phiMass)<0.02) hmBs->Fill((lKK+ljp).mass());
        }
//      double massBpm = sqrt( sqr(ljp2.E()+lc1.E())-sqr(ljp2.px()+lc1.px())-sqr(ljp2.py()+lc1.py())-sqr(ljp2.pz()+lc1.pz()));
//      if (debug) std::cout<<" candidate: prob:"<<probBpm<<" mass: "<<lBpm.mass()<<" vs: "<<massBpm<<std::endl;
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
    if(vtx.numberOfDaughters()==2) {
      const reco::Candidate * d1 = vtx.daughter(0);
      const reco::Candidate * d2 = vtx.daughter(1);
      ROOT::Math::PxPyPzEVector lpi_A(lorentzVector(d1->momentum(), pionMass));
      ROOT::Math::PxPyPzEVector lpi_B(lorentzVector(d2->momentum(), pionMass));
      ROOT::Math::PxPyPzEVector lpr_A(lorentzVector(d2->momentum(), protonMass));
      ROOT::Math::PxPyPzEVector lpr_B(lorentzVector(d1->momentum(), protonMass));
      hmLam->Fill((lpi_A+lpr_A).mass());
      hmLam->Fill((lpi_B+lpr_B).mass());
      if (debug) std::cout <<" Lambda mass 1: "<<(lpi_A+lpr_A).mass()<<" mass2: "<<(lpi_B+lpr_B).mass()<<std::endl;

    } 
  }
  */

  //write std io
  if (debug) cout <<"*** Analyze event: " << ev.id()<<" analysed event count:"<<++theEventCount << endl;
}

DEFINE_FWK_MODULE(Analysis);

