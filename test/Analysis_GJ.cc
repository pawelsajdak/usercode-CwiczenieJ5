// system include files
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <bitset>
#include <map>
#include <string>

// user include files
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

//#include "DataFormats/MuonDetId/interface/DtDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"


#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "TH2D.h"
#include "TFile.h"

const TrackingParticle & ancestor(const TrackingParticle & particle) {

  const TrackingVertexRef&  tpv = particle.parentVertex();
  if (tpv->nSourceTracks() == 0) return particle;
  const TrackingParticle & parent =  **(tpv->sourceTracks_begin());
  return ancestor(parent);
}

std::string print(const TrackingParticle & tp) {
  std::stringstream ss;
  ss <<" pid: "<<tp.pdgId()
     <<" pt: "<<tp.pt()
     <<" eta: "<<tp.eta()
     <<" phi: "<<tp.phi()
     <<" vtx[r,z]:  ["<<tp.parentVertex()->position().Rho() <<", "<<tp.parentVertex()->position().z()<<"]"
     <<" time: "<<tp.parentVertex()->position().T()
     ;
  return ss.str();
}



class Analysis : public edm::one::EDAnalyzer<> {
public:
  Analysis (const edm::ParameterSet & cfg);
  virtual ~Analysis(){ 
    TFile f("histos.root","RECREATE");
    f.cd();
    hLicExample->Write(); 
    f.Write();
    std::cout << "AT THE END: "; printStat(); 
  }

  virtual void analyze(const edm::Event &ev, const edm::EventSetup& es) {
    theEventCnt++;
    debug = 1;
    analyzeDT(ev,es);

    if (theEventCnt%100==1) printStat(); 
  }
  void analyzeDT(const edm::Event&, const edm::EventSetup& es);
  void printStat();

private:
  edm::EDGetTokenT<L1MuDTChambPhContainer> inputDTPh_leg;
  edm::EDGetTokenT<L1MuDTChambThContainer> inputDTTh_leg;
  edm::EDGetTokenT<TrackingParticleCollection> inputTP;
//  edm::EDGetTokenT<TrackingVertexCollection> inputTV, inputTV0;

 edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> theGeomteryToken;
 edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> theFieldToken;
 edm::ESGetToken<Propagator, TrackingComponentsRecord> thePropagatorToken;


  bool debug;
  unsigned int theEventCnt;
  unsigned int theAllDtPDigisCnt, theAllDtTDigisCnt;
  TH2D *hLicExample;
};

void Analysis::printStat()
{
   std::cout<<std::dec <<"=========> Analyzed #"<<theEventCnt
            <<" DtPh: "<<theAllDtPDigisCnt <<"  DtTh: "<<theAllDtTDigisCnt <<std::endl; 
}

Analysis::Analysis(const edm::ParameterSet & cfg) 
  : debug (false), theEventCnt(0),
    theAllDtPDigisCnt(0),
    theAllDtTDigisCnt(0)
{
  inputDTPh_leg = consumes<L1MuDTChambPhContainer>(cfg.getParameter<edm::InputTag>("srcDTPh_leg"));
  inputDTTh_leg = consumes<L1MuDTChambThContainer>(cfg.getParameter<edm::InputTag>("srcDTTh_leg"));
  inputTP  =   consumes<TrackingParticleCollection>(edm::InputTag("mix","MergedTrackTruth"));
//  inputTV  =   consumes<TrackingVertexCollection>(edm::InputTag("mix","MergedTrackTruth"));
//  inputTV0 =   consumes<TrackingVertexCollection>(edm::InputTag("mix","InitialVertices"));
//  : theGeomteryToken(cColl.esConsumes()), theFieldToken(cColl.esConsumes()),
  theGeomteryToken=esConsumes<GlobalTrackingGeometry, GlobalTrackingGeometryRecord>();
  theFieldToken=esConsumes<MagneticField, IdealMagneticFieldRecord>();   
  thePropagatorToken=esConsumes<Propagator, TrackingComponentsRecord>(edm::ESInputTag("","SteppingHelixPropagatorAny")); 


  hLicExample = new TH2D("hLicExample","hLicExample", 12,0.5,12.5, 8,-0.5,7.5);
}

void Analysis::analyzeDT( const edm::Event &ev, const edm::EventSetup& es) {
  if (debug) std::cout << "-------- Tracking Particles -----------" << std::endl;

  GlobalTrajectoryParameters gtp;
  const MagneticField * field = &es.getData(theFieldToken);
  const GlobalTrackingGeometry & globalGeometry = es.getData(theGeomteryToken);



  edm::Handle<TrackingParticleCollection> tpColl;
  ev.getByToken(inputTP, tpColl);
  const TrackingParticleCollection & myTP = *(tpColl.product());
  if (debug) std::cout<<" TRACKING PARTICLES: " << myTP.size() << std::endl;
  for (const auto & tp : myTP) {
//    if ( abs( tp.pdgId())!=13  || tp.pt() < 1. || tp.parentVertex()->position().Rho()>200. ||  fabs(tp.parentVertex()->position().T()) > 24.) continue;
    if ( abs( tp.pdgId())==13  && tp.pt() > 1. && gtp.charge()==0) {
      const auto & vtx = tp.parentVertex()->position();
      const auto & mom = tp.momentum();
      gtp=GlobalTrajectoryParameters(GlobalPoint(vtx.x(), vtx.y(), vtx.z()), GlobalVector(mom.x(), mom.y(), mom.z()), tp.charge(), field);
    } 
    if (debug) std::cout << print(tp)<<std::endl;
  }


  if (debug) std::cout << "-------- HERE DIGI COMPARE DT ---------" << std::endl;
  edm::Handle<L1MuDTChambPhContainer> digiCollectionDTPh_leg;
  ev.getByToken(inputDTPh_leg, digiCollectionDTPh_leg);
  const L1MuDTChambPhContainer& dtphDigisLeg= *digiCollectionDTPh_leg.product();
  if (debug) std::cout <<" DTPh digis from BMTF " << dtphDigisLeg.getContainer()->size()<< std::endl;
  for (const auto &  chDigi : *dtphDigisLeg.getContainer() ) {
    if (abs(chDigi.whNum()) != 2) continue;
    if (chDigi.stNum() ==4) continue;
    if (chDigi.bxNum() != 0) continue;
    if (chDigi.code()==7) continue;
    theAllDtPDigisCnt++;
    if (debug) std::cout <<"DtDataWord64 BMTF    " 
        <<" bxNum: "<<chDigi.bxNum()
        <<" whNum: "<<chDigi.whNum()
        <<" station: "<< chDigi.stNum()
        <<" sector: "<< chDigi.scNum()
        <<" phi:   "<<chDigi.phi()
        <<" phiB:   "<<chDigi.phiB()
        <<" code(q): "<< chDigi.code()
        << std::endl;
    hLicExample->Fill(chDigi.scNum()+1,chDigi.code());
    if (gtp.charge() !=0) {
      FreeTrajectoryState fts(gtp);
      DTChamberId dtChamberId(chDigi.whNum(),chDigi.stNum(),chDigi.scNum()+1);
      // GlobalPoint detPosition = globalGeometry.idToDet(dtChamberId)->position();
      //const DTChamber * chamber = geomDt.chamber(dtChamberId);
      const GeomDet * geomDet = globalGeometry.idToDet(dtChamberId);
      const Propagator & propagator = es.getData(thePropagatorToken);
      TrajectoryStateOnSurface stateAtDet =  propagator.propagate(fts, geomDet->surface());
      std::cout <<" Is PROPAGATION VALID: "<<stateAtDet.isValid() <<std::endl;
    }
  }

/*
  edm::Handle<L1MuDTChambThContainer> digiCollectionDTTh_BMTF;
  ev.getByToken(inputDTTh_BMTF, digiCollectionDTTh_BMTF);
  const L1MuDTChambThContainer& dtthDigisBMTF = *digiCollectionDTTh_BMTF.product();
  if (debug) std::cout <<" DTTh digis from BMTF " << dtthDigisBMTF.getContainer()->size()<< std::endl;
  for (const auto &  chDigi : *dtthDigisBMTF.getContainer() ) {
    unsigned int eta = 0;
    unsigned int etaQ = 0;
    for (unsigned int ipos=0; ipos <7; ipos++) {
     if (chDigi.position(ipos) >1 ) if (debug) std::cout <<" HERE PROBLEM !!!!" << std::endl;
     if (chDigi.position(ipos)==1) eta |= (1 <<ipos);
     if (chDigi.quality(ipos)==1) etaQ |= (1 <<ipos);
    }
    if (chDigi.bxNum() != 0) continue;
//    if (abs(chDigi.bxNum()) >2) continue;
    if (eta==0 || etaQ==0) continue;
    if (abs(chDigi.whNum()) != 2) continue;
    if (debug) std::cout <<"DtDataWord64 BMTF TH " 
        <<" bxNum: "<<chDigi.bxNum()
        <<" whNum: "<<chDigi.whNum()
        <<" station: "<< chDigi.stNum()
        <<" sector: "<< chDigi.scNum()
        <<" eta: " << eta
        <<" etaQ: " << etaQ
        << std::endl;
  }
*/
}

DEFINE_FWK_MODULE(Analysis);
