// -*- C++ -*-
//
// Package:    UserData/nTupleMaker
// Class:      nTupleMaker
// 
/**\class nTupleMaker nTupleMaker.cc UserData/nTupleMaker/plugins/nTupleMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Halil Gamsizkan
//         Created:  Wed, 27 Jan 2016 13:18:02 GMT
//
//

//see also: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookWriteFrameworkModule?LOCALSHELL=bash

//#define ELECSCOLLTYPE reco::TrackCollection
#define ELECSCOLLTYPE pat::ElectronCollection
#define MUONSCOLLTYPE pat::MuonCollection
#define PHOTSCOLLTYPE pat::PhotonCollection
#define JETSCOLLTYPE pat::JetCollection
#define METSCOLLTYPE pat::METCollection

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//PAT Data types
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/PatCandidates/interface/Jet.h>

// Added for tracks
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Added for TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TTree.h"

// Added for accessing generator information
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

// Own source
#include "event.h"

//
// class declaration
//

class nTupleMaker : public edm::EDAnalyzer {
   public:
      explicit nTupleMaker(const edm::ParameterSet&);
      ~nTupleMaker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      int _debug;
      
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // added for token handling
      edm::EDGetTokenT< MUONSCOLLTYPE > tok_muons_;
      edm::EDGetTokenT< ELECSCOLLTYPE > tok_electrons_;
      edm::EDGetTokenT< PHOTSCOLLTYPE > tok_photons_;
      edm::EDGetTokenT< JETSCOLLTYPE > tok_jets_;
      edm::EDGetTokenT< METSCOLLTYPE > tok_METs_;
      edm::EDGetTokenT< GenEventInfoProduct > tok_gen_;
      
      // cfg communication
      //unsigned int minTracks_;
      bool addElec_;
      bool addPhot_;
      bool addMuons_;
      bool addJets_;
      bool addMET_;
      bool addBTag_;

      std::string labelElec_;
      std::string labelPhot_;
      std::string labelMuons_;
      std::string labelJets_;
      std::string labelMET_;
      std::string algoBTag_;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      //TH1D * histo;
      
      event myEvent;
      TTree * eventTree;
      //static const int maxReco = 30; // need tp check!!!!
      //#define maxReco 30
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
nTupleMaker::nTupleMaker(const edm::ParameterSet& iConfig):
   //minTracks_(iConfig.getUntrackedParameter<unsigned int>("minTracks",0)),
   addElec_(iConfig.getUntrackedParameter<bool>("addElec",0)),
   addPhot_(iConfig.getUntrackedParameter<bool>("addPhot",0)),
   addMuons_(iConfig.getUntrackedParameter<bool>("addMuons",0)),
   addJets_(iConfig.getUntrackedParameter<bool>("addJets",0)),
   addMET_(iConfig.getUntrackedParameter<bool>("addMET",0)),
   addBTag_(iConfig.getUntrackedParameter<bool>("addBTag",0)),
   labelElec_ ( iConfig.getUntrackedParameter<std::string>( "labelElec" ) ),
   labelPhot_ ( iConfig.getUntrackedParameter<std::string>( "labelPhot" ) ),
   labelMuons_ ( iConfig.getUntrackedParameter<std::string>( "labelMuons" ) ),
   labelJets_ ( iConfig.getUntrackedParameter<std::string>( "labelJets" ) ),
   labelMET_ ( iConfig.getUntrackedParameter<std::string>( "labelMET" ) ),
   algoBTag_ ( iConfig.getUntrackedParameter<std::string>( "algoBTag" ) )
{
   _debug=0;
   //now do what ever initialization is needed
   edm::Service<TFileService> fs;
   eventTree = fs->make<TTree>("Events","Events");
//   histo = fs->make<TH1D>("charge" , "Charges" , 200 , -2 , 2 );
   //charge, n_reco;
   //px, py, pz, phi, theta, eta;

   tok_gen_ = consumes< GenEventInfoProduct >(edm::InputTag("generator"));
   eventTree->Branch("gen_weight",&myEvent.gen_weight,"gen_weight/D");

   if (addElec_){
      tok_electrons_ = consumes< ELECSCOLLTYPE >(edm::InputTag(labelElec_));
      eventTree->Branch("ele_n",&myEvent.ele_n,"ele_n/I");
      eventTree->Branch("ele_charge",myEvent.ele_charge,"charge[ele_n]/I");
      eventTree->Branch("ele_px",myEvent.ele_px,"ele_px[ele_n]/D");
      eventTree->Branch("ele_py",myEvent.ele_py,"ele_py[ele_n]/D");
      eventTree->Branch("ele_pz",myEvent.ele_pz,"ele_pz[ele_n]/D");
      eventTree->Branch("ele_phi",myEvent.ele_phi,"ele_phi[ele_n]/D");
      eventTree->Branch("ele_theta",myEvent.ele_theta,"ele_theta[ele_n]/D");
      eventTree->Branch("ele_eta",myEvent.ele_eta,"ele_eta[ele_n]/D");
      eventTree->Branch("ele_id",myEvent.ele_id,"ele_id[ele_n]/D");
      
   }
   
   if (addPhot_){
      tok_photons_ = consumes< PHOTSCOLLTYPE >(edm::InputTag(labelPhot_));
      eventTree->Branch("phot_n",&myEvent.phot_n,"phot_n/I");
      eventTree->Branch("phot_px",myEvent.phot_px,"phot_px[phot_n]/D");
      eventTree->Branch("phot_py",myEvent.phot_py,"phot_py[phot_n]/D");
      eventTree->Branch("phot_pz",myEvent.phot_pz,"phot_pz[phot_n]/D");
      eventTree->Branch("phot_phi",myEvent.phot_phi,"phot_phi[phot_n]/D");
      eventTree->Branch("phot_theta",myEvent.phot_theta,"phot_theta[phot_n]/D");
      eventTree->Branch("phot_eta",myEvent.phot_eta,"phot_eta[phot_n]/D");
   }    
   
   if (addMuons_){
      tok_muons_ = consumes< MUONSCOLLTYPE >(edm::InputTag(labelMuons_));
      eventTree->Branch("mu_n",&myEvent.mu_n,"mu_n/I");
      eventTree->Branch("mu_charge",myEvent.mu_charge,"charge[mu_n]/I");
      eventTree->Branch("mu_px",myEvent.mu_px,"mu_px[mu_n]/D");
      eventTree->Branch("mu_py",myEvent.mu_py,"mu_py[mu_n]/D");
      eventTree->Branch("mu_pz",myEvent.mu_pz,"mu_pz[mu_n]/D");
      eventTree->Branch("mu_phi",myEvent.mu_phi,"mu_phi[mu_n]/D");
      eventTree->Branch("mu_theta",myEvent.mu_theta,"mu_theta[mu_n]/D");
      eventTree->Branch("mu_eta",myEvent.mu_eta,"mu_eta[mu_n]/D");
   }
   
   if (addJets_){
      tok_jets_ = consumes< JETSCOLLTYPE >(edm::InputTag(labelJets_));
      eventTree->Branch("jet_n",&myEvent.jet_n,"jet_n/I");
      eventTree->Branch("jet_px",myEvent.jet_px,"jet_px[jet_n]/D");
      eventTree->Branch("jet_py",myEvent.jet_py,"jet_py[jet_n]/D");
      eventTree->Branch("jet_pz",myEvent.jet_pz,"jet_pz[jet_n]/D");
      eventTree->Branch("jet_phi",myEvent.jet_phi,"jet_phi[jet_n]/D");
      eventTree->Branch("jet_theta",myEvent.jet_theta,"jet_theta[jet_n]/D");
      eventTree->Branch("jet_eta",myEvent.jet_eta,"jet_eta[jet_n]/D");
      
      
      //jet_nhf = Branch("jet_nhf",myEvent.jet_nhf,"jet_nhf[jet_n]/D");
      //jet_nef = Branch("jet_nef",myEvent.jet_nef,"jet_nef[jet_n]/D");
      //jet_chf = Branch("jet_chf",myEvent.jet_chf,"[jet_n]/D");
      //jet_cef = Branch("jet_cef",myEvent.jet_cef,"jet_cef[jet_n]/D");
      //jet_nconstituents = Branch("jet_nconstituents",myEvent.jet_nconstituents,"jet_nconstituents[jet_n]/I");
      //jet_nch = Branch("jet_nch",myEvent.jet_nch,"jet_nch[jet_n]/I");
   }
   
   if (addMET_){
      tok_METs_ = consumes< METSCOLLTYPE >(edm::InputTag(labelMET_));
      eventTree->Branch("MET_n",&myEvent.MET_n,"MET_n/I");
      eventTree->Branch("MET_px",myEvent.MET_px,"MET_px[MET_n]/D");
      eventTree->Branch("MET_py",myEvent.MET_py,"MET_py[MET_n]/D");
      eventTree->Branch("MET_pz",myEvent.MET_pz,"MET_pz[MET_n]/D");
      eventTree->Branch("MET_phi",myEvent.MET_phi,"MET_phi[MET_n]/D");
      eventTree->Branch("MET_theta",myEvent.MET_theta,"MET_theta[MET_n]/D");
      eventTree->Branch("MET_eta",myEvent.MET_eta,"MET_eta[MET_n]/D");
   }
   
}


nTupleMaker::~nTupleMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
nTupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // Access the gen information
   {
      std::string labelGen_("generator");
      Handle<GenEventInfoProduct> genEvtInfo;
      try{iEvent.getByToken( tok_gen_, genEvtInfo );}
      catch( cms::Exception& ex ) { LogError("nTupleMaker") << "Product: " << labelGen_ << " not found"; }
      if (!genEvtInfo.isValid())
      throw cms::Exception("ProductNotFound: ") << labelGen_ << " product not found";
   
      //double qScale = genEvtInfo->qScale();  // in case of Pythia6, this will be pypars/pari(23)
      //std::vector<double>& evtWeights = genEvtInfo->weights();
      myEvent.gen_weight = genEvtInfo->weight();
      if (_debug) std::cout << "Weight: " << myEvent.gen_weight << std::endl;
   }
   
   if (addElec_){
      Handle<ELECSCOLLTYPE> elec;
      
      try {iEvent.getByToken( tok_electrons_, elec );} //try {iEvent.getByLabel(labelElec_, elec); }
      catch( cms::Exception& ex ) { LogError("nTupleMaker") << "Product: " << labelElec_ << " not found"; }
      if (!elec.isValid())
      throw cms::Exception("ProductNotFound") << labelElec_ << " product not found";
   
      if (_debug) std::cout << elec->size() << std::endl;
      myEvent.ele_n=elec->size();
      size_t idx = 0;
      for( ELECSCOLLTYPE::const_iterator ele = elec->begin(); 
            ele != elec->end(); ++ ele, ++ idx ) {
            if (idx > maxReco) break;
            // dump electron IDs
            //std::cout << ele->electronIDs().size() << std::endl;
            //for (unsigned i=0; i < ele->electronIDs().size(); i++) {
            //   pat::Electron::IdPair idp = ele->electronIDs().at(i);
            //   std::cout << idp.first << ", " << idp.second << std::endl;
            //}
            //std::cout << ele->isElectronIDAvailable("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V2-standalone-veto") << std::endl;
            if (_debug) LogInfo("nTupleMaker") << "Electron #" << idx << ": " 
                 << "charge: " << ele->charge() 
                 << ", px: " << ele->px()
                 << ", py: " << ele->py()
                 << ", pz: " << ele->pz()
                 << ", phi: " << ele->phi()
                 << ", theta: " << ele->theta()
                 << ", eta: " << ele->eta();
            myEvent.ele_charge[idx]=ele->charge();
            myEvent.ele_px[idx]=ele->px();
            myEvent.ele_py[idx]=ele->py();
            myEvent.ele_pz[idx]=ele->pz();
            myEvent.ele_phi[idx]=ele->phi();
            myEvent.ele_theta[idx]=ele->theta();
            myEvent.ele_eta[idx]=ele->eta();
            myEvent.ele_id[idx]=ele->electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight");
      }
   } // end electron code block
   
   if (addPhot_){
      Handle<PHOTSCOLLTYPE> particles;
      try {iEvent.getByToken( tok_photons_, particles );} //try {iEvent.getByLabel(labelPhot_, particles); }
      catch( cms::Exception& ex ) { LogError("nTupleMaker") << "Product: " << labelPhot_ << " not found"; }
      if (!particles.isValid())
      throw cms::Exception("ProductNotFound") << labelPhot_ << " product not found";
   
      if (_debug) std::cout << particles->size() << std::endl;
      myEvent.phot_n=particles->size();
      size_t idx = 0;
      for( PHOTSCOLLTYPE::const_iterator particle = particles->begin(); 
            particle != particles->end(); ++ particle, ++ idx ) {
            if (idx > maxReco) break;
            if (_debug) LogInfo("nTupleMaker") << "Photon #" << idx << ": " 
                 << ", px: " << particle->px()
                 << ", py: " << particle->py()
                 << ", pz: " << particle->pz()
                 << ", phi: " << particle->phi()
                 << ", theta: " << particle->theta()
                 << ", eta: " << particle->eta();
            myEvent.phot_px[idx]=particle->px();
            myEvent.phot_py[idx]=particle->py();
            myEvent.phot_pz[idx]=particle->pz();
            myEvent.phot_phi[idx]=particle->phi();
            myEvent.phot_theta[idx]=particle->theta();
            myEvent.phot_eta[idx]=particle->eta();
      }
   } // end photon code block
   
   if (addMuons_){
      Handle<MUONSCOLLTYPE> muons;
      
      try {iEvent.getByToken( tok_muons_, muons );} //iEvent.getByLabel(labelMuons_, muons);  }
      catch( cms::Exception& ex ) { LogError("nTupleMaker") << "Product: " << labelMuons_ << " not found"; }
      if (!muons.isValid())
      throw cms::Exception("ProductNotFound") << labelMuons_ << " product not found";      
      
      if (_debug) std::cout << muons->size() << std::endl;
      //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMuonAnalysis
      //http://cmslxr.fnal.gov/lxr/source/DataFormats/Candidate/interface/LeafCandidate.h?v=CMSSW_7_6_1
      myEvent.mu_n=muons->size();
      size_t idx = 0;
      for( MUONSCOLLTYPE::const_iterator muon = muons->begin(); 
            muon != muons->end(); ++ muon, ++ idx ) {
            if (idx > maxReco) break;
            // do something, cand is a pointer to your object
            if (_debug) LogInfo("nTupleMaker") << "Muon #" << idx << ": " 
                 << "charge: " << muon->charge() 
                 << ", px: " << muon->px()
                 << ", py: " << muon->py()
                 << ", pz: " << muon->pz()
                 << ", phi: " << muon->phi()
                 << ", theta: " << muon->theta()
                 << ", eta: " << muon->eta();
            myEvent.mu_charge[idx]=muon->charge();
            myEvent.mu_px[idx]=muon->px();
            myEvent.mu_py[idx]=muon->py();
            myEvent.mu_pz[idx]=muon->pz();
            myEvent.mu_phi[idx]=muon->phi();
            myEvent.mu_theta[idx]=muon->theta();
            myEvent.mu_eta[idx]=muon->eta();
      }
   } // end muon code block
   
   if (addJets_){
      Handle<JETSCOLLTYPE> particles;
      try {iEvent.getByToken( tok_jets_, particles );} //try {iEvent.getByLabel(labelJets_, particles); }
      catch( cms::Exception& ex ) { LogError("nTupleMaker") << "Product: " << labelJets_ << " not found"; }
      if (!particles.isValid())
      throw cms::Exception("ProductNotFound") << labelJets_ << " product not found";
   
      if (_debug) std::cout << particles->size() << std::endl;
      myEvent.phot_n=particles->size();
      size_t idx = 0;
      for( JETSCOLLTYPE::const_iterator particle = particles->begin(); 
            particle != particles->end(); ++ particle, ++ idx ) {
            if (idx > maxReco) break;
            if (_debug) LogInfo("nTupleMaker") << "Jet #" << idx << ": " 
                 << ", px: " << particle->px()
                 << ", py: " << particle->py()
                 << ", pz: " << particle->pz()
                 << ", phi: " << particle->phi()
                 << ", theta: " << particle->theta()
                 << ", eta: " << particle->eta();
            myEvent.jet_px[idx]=particle->px();
            myEvent.jet_py[idx]=particle->py();
            myEvent.jet_pz[idx]=particle->pz();
            myEvent.jet_phi[idx]=particle->phi();
            myEvent.jet_theta[idx]=particle->theta();
            myEvent.jet_eta[idx]=particle->eta();
            //myEvent.jet_nhf = jet.neutralHadronEnergy() / uncorrJet.E();
            //myEvent.jet_nef = jet.neutralEmEnergy() / uncorrJet.E();
            //myEvent.jet_chf = jet.chargedHadronEnergy() / uncorrJet.E();
            //myEvent.jet_cef = jet.chargedEmEnergy() / uncorrJet.E();
            //myEvent.jet_nconstituents = particle->numberOfDaughters();
            //myEvent.jet_nch = particle->chargedMultiplicity();
      }
   } // end jet code block   

   if (addMET_){
      Handle<METSCOLLTYPE> particles;
      try {iEvent.getByToken( tok_METs_, particles );} //try {iEvent.getByLabel(labelMET_, particles); }
      catch( cms::Exception& ex ) { LogError("nTupleMaker") << "Product: " << labelMET_ << " not found"; }
      if (!particles.isValid())
      throw cms::Exception("ProductNotFound") << labelMET_ << " product not found";
   
      if (_debug) std::cout << particles->size() << std::endl;
      myEvent.phot_n=particles->size();
      size_t idx = 0;
      for( METSCOLLTYPE::const_iterator particle = particles->begin(); 
            particle != particles->end(); ++ particle, ++ idx ) {
            if (idx > maxReco) break;
            if (_debug) LogInfo("nTupleMaker") << "MET #" << idx << ": " 
                 << ", px: " << particle->px()
                 << ", py: " << particle->py()
                 << ", pz: " << particle->pz()
                 << ", phi: " << particle->phi()
                 << ", theta: " << particle->theta()
                 << ", eta: " << particle->eta();
            myEvent.MET_px[idx]=particle->px();
            myEvent.MET_py[idx]=particle->py();
            myEvent.MET_pz[idx]=particle->pz();
            myEvent.MET_phi[idx]=particle->phi();
            myEvent.MET_theta[idx]=particle->theta();
            myEvent.MET_eta[idx]=particle->eta();
      }
   } // end jet code block
   
   eventTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
nTupleMaker::beginJob()
{
   edm::LogInfo("nTupleMaker") << "Starting";
   std::cout << "Ntuplize: e " << addElec_ << " : gm " << addPhot_ << " : mu " << addMuons_ << " : j " << addJets_ << " : MET " << addMET_ << " : bT" << addBTag_ << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
nTupleMaker::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
nTupleMaker::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
nTupleMaker::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
nTupleMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
nTupleMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
nTupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(nTupleMaker);