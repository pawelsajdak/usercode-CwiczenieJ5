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
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include "L1Trigger/RPCTrigger/interface/RPCConst.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <sstream>



using namespace std;

//object definition
class Cwiczenie : public edm::EDAnalyzer {
public:

  //constructor, function is called when new object is created
  explicit Cwiczenie(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~Cwiczenie();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;
  //added variables
  TFile *myRootFile;
  TH1D *histo;

  edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection> inputOMTF;
  edm::EDGetTokenT<edm::SimTrackContainer> inputSim;

};


Cwiczenie::Cwiczenie(const edm::ParameterSet& conf) 
  : theConfig(conf), theEventCount(0) 
{
  cout <<" CTORXX" << endl;
  inputOMTF = consumes<l1t::RegionalMuonCandBxCollection>(edm::InputTag("simOmtfDigis","OMTF"));
  inputSim =  consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));
}


Cwiczenie::~Cwiczenie() 
{ 
  cout <<" DTOR" << endl;
}

void Cwiczenie::beginJob()
{
  //make a new Root file
  myRootFile=new TFile("test.root","RECREATE");
  //create a histogram
  histo =new TH1D("histo","test; phi omtf; #events",360, 0., 2*M_PI);
  cout << "HERE Cwiczenie::beginJob()" << endl;
}

void Cwiczenie::endJob()
{
  //write histogram data
  histo->Write();
  myRootFile->Close();
  delete histo;
  delete myRootFile;
  cout << "HERE Cwiczenie::endJob()" << endl;
}

void Cwiczenie::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " HERE Cwiczenie::analyze "<< std::endl;

  float pt_sim, phi_sim, eta_sim;
  float pt_omtf, phi_omtf, eta_omtf;
  unsigned int simMuonCount=0; 
  unsigned int omtfMuonCount=0;

  std::cout <<" SIMULATED MUONS: "<<std::endl;
  edm::Handle<edm::SimTrackContainer> simTk;
  ev.getByToken(inputSim, simTk);
  std::vector<SimTrack> mySimTracks = *(simTk.product());
  for (std::vector<SimTrack>::const_iterator it=mySimTracks.begin(); it<mySimTracks.end(); it++) {
    const SimTrack & track = *it;
    if ( track.type() == -99) continue;
    if ( track.vertIndex() != 0) continue;

    //sucessful muon, add to count
    simMuonCount++;

    phi_sim = track.momentum().phi(); //momentum azimutal angle
    pt_sim = track.momentum().pt(); //transverse momentum
    eta_sim = track.momentum().eta(); //pseudorapidity
    std::cout <<" trackId: " <<track.trackId() 
          << " pt_sim: " << pt_sim <<" eta_sim: "<<eta_sim<<" phi_sim: "<<phi_sim
          <<" vtx: "<<track.vertIndex()<<" type: "<<track.type() // 13 or -13 is a muon
          << std::endl; 
  }
  if(simMuonCount!=1) {
     cout<<"    Simulated muon count != 1"<<endl;
     return;
  }

  std::cout <<" L1 MUONS: "<<std::endl;
  edm::Handle<l1t::RegionalMuonCandBxCollection> l1Omtf;
  ev.getByToken(inputOMTF, l1Omtf);
  int bxNumber = 0;
  for (l1t::RegionalMuonCandBxCollection::const_iterator it = l1Omtf.product()->begin(bxNumber);
       it != l1Omtf.product()->end(bxNumber); ++it) {
    omtfMuonCount++;
    pt_omtf  =  (it->hwPt()-1.)/2.;
    phi_omtf = ( (15.+it->processor()*60.)/360. + it->hwPhi()/576. ) *2*M_PI; 
    if (phi_omtf > 2*M_PI) phi_omtf -=  2*M_PI;
    eta_omtf = it->hwEta()/240.*2.26;
    std::cout<<" Processor : "<<it->processor()  <<" pT: "<<it->hwPt()<<" eta: "<<it->hwEta()<<" phi: "<<it->hwPhi()<< std::endl;
    std::cout<<" pT: "<< pt_omtf <<" phi: "<<phi_omtf<<" eta: "<<eta_omtf<<std::endl;
  }
  if(omtfMuonCount != 1) cout<<"    OMTF muon count != 1"<<endl;

  if(omtfMuonCount==1) histo->Fill(phi_omtf);
  
  //write std io
  cout <<"*** Cwiczenie, analyze event: " << ev.id()<<" useful event count:"<<++theEventCount << endl;
}




DEFINE_FWK_MODULE(Cwiczenie);

