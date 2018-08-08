#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <ValidationHHtobbtt/MyEventLoop.h>

#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEvent.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

#include "TLorentzVector.h"

// this is needed to distribute the algorithm to the workers
ClassImp(MyEventLoop)



MyEventLoop :: MyEventLoop ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



EL::StatusCode MyEventLoop :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  xAOD::Init( "ValidationHHtobbtt" ).ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyEventLoop :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  h_graviton_mass = new TH1F("h_graviton_mass","; G Mass (GeV); #countes", 100,0,4000);
  h_hh_mass = new TH1F("h_hh_mass","; hh Mass (GeV); #counts",100,0,4000);
  h_higgs_pT = new TH1F("h_higgs_pT","; Higgs p_{T} (GeV); #counts",100,0,2000);
  h_bquark_pT = new TH1F("h_bquark_pT","; b-quark p_{T} (GeV); #counts",100,0,2000);
  h_bquark_eta = new TH1F("h_bquark_eta","; b-quark #eta; #counts", 100,-4,4);
  h_bb_pT = new TH1F("h_bb_pT","; p_{T}^{bb} (GeV); #counts",100,0,2000);

  h_tauh_pT = new TH1F("h_tauh_pT","; Visible Tau p_{T} (GeV); #counts",100,0,2000);
  h_tauh_eta = new TH1F("h_tauh_eta","; Visible Tau #eta; #counts",100,-4,4);

  h_tau_pT = new TH1F("h_tau_pT","; Tau p_{T} (GeV); #counts",100,0,2000);
  h_tau_eta = new TH1F("h_tau_eta","; Tau #eta; #counts",100,-4,4);

  h_el_pT = new TH1F("h_el_pT","; electron p_{T} (GeV); #counts",100,0,2000);
  h_el_eta = new TH1F("h_el_eta","; electron #eta; #counts",100,-4,4);

  h_mu_pT = new TH1F("h_mu_pT","; muon pT (GeV); #counts",100,0,2000);
  h_mu_eta = new TH1F("h_mu_eta","; muon #eta; #counts",100,-4,4);

  h_tautau_pT = new TH1F("h_tautau_pT","; p_{T}^{tautau, full} (GeV); #counts",100,0,2000);
  h_tautau_vis_pT = new TH1F("h_tautau_vis_pT","; p_{vis}^{tautau, vis} (GeV); #counts",100,0,2000);
  h_tautau_hh_mass = new TH1F("h_tautau_hh_mass","; H->tautau mass (GeV); #counts",100,0,300);
  h_tautau_hh_vis_mass = new TH1F("h_tautau_hh_vis_mass","; H->tautau visible mass (GeV); #counts",100,0,300);

  h_ntaus_higgs = new TH1F("h_ntaus_higgs","; N taus; #counts",5,-0.5,4.5);
  h_npho_taurad = new TH1F("h_npho_tau","; N photons; #counts",5,-0.5,4.5);
  h_pho_taurad_pT = new TH1F("h_pho_taurad_pT","; p_{T}^{photon} (GeV); #counts",100,0,2000);
  h_tau_taurad_pT = new TH1F("h_tau_taurad_pT","; p_{T}^{tautau, full} (GeV); #counts",100,0,2000);
  h_tau_notaurad_pT = new TH1F("h_tau_notaurad_pT","; p_{T}^{tautau, full} (GeV); #counts",100,0,2000);

  wk()->addOutput(h_graviton_mass);
  wk()->addOutput(h_hh_mass);
  wk()->addOutput(h_higgs_pT);
  wk()->addOutput(h_bquark_pT);
  wk()->addOutput(h_bquark_eta);
  wk()->addOutput(h_bb_pT);
  wk()->addOutput(h_tautau_pT);
  wk()->addOutput(h_tautau_vis_pT);
  wk()->addOutput(h_tau_pT);
  wk()->addOutput(h_tau_eta);
  wk()->addOutput(h_tauh_pT);
  wk()->addOutput(h_tauh_eta);
  wk()->addOutput(h_tautau_hh_mass);
  wk()->addOutput(h_el_pT);
  wk()->addOutput(h_el_eta);
  wk()->addOutput(h_mu_pT);
  wk()->addOutput(h_mu_eta);
  wk()->addOutput(h_tautau_hh_vis_mass);
  wk()->addOutput(h_ntaus_higgs);
  wk()->addOutput(h_npho_taurad);
  wk()->addOutput(h_pho_taurad_pT);
  wk()->addOutput(h_tau_taurad_pT);
  wk()->addOutput(h_tau_notaurad_pT);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyEventLoop :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyEventLoop :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyEventLoop :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.
  m_event = wk()->xaodEvent();
  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

  // count number of events
  m_eventCounter = 0;

  nhtaus=0;
  neltaus=0;
  nmutaus=0;
  ntausall=0;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyEventLoop :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  if(m_eventCounter%1000 ==0) Info("execute()", "Event number = %i", m_eventCounter );
  m_eventCounter++;

  //----------------------------
  // Event information
  //---------------------------
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //----------------------------
  //    Truth information
  //----------------------------
  const xAOD::TruthParticleContainer* m_TruthParticleContainer = 0;
  if( m_event->contains<xAOD::TruthParticleContainer>("TruthParticles") ){
     if( !m_event->retrieve( m_TruthParticleContainer, "TruthParticles" ) ){
      Error("execute()", "Failed to retrieve TruthParticle. Exiting." );
      return EL::StatusCode::FAILURE;
     }
  } else {
     Warning("execute()", "Did not find TruthParticle.");
  }

  //----------------------------
  //  Truth Jet information
  //----------------------------
  /*const xAOD::JetContainer* jet_vec = 0;
  //const DataVector<xAOD::TauJet>* tau_vec = 0;
  if( ! m_event->retrieve( jet_vec, "Jets").isSuccess() ){
    Error("execute()", "Failed to retrieve JetContainer collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }*/


  //const xAOD::JetContainer* m_AntiKt4TruthContainer = 0;
  //if( m_event->contains<xAOD::JetContainer>("AntiKt4TruthJets") ){
  //   if( !m_event->retrieve( m_AntiKt4TruthContainer, "AntiKt4TruthJets" ) ){
  //    Error("execute()", "Failed to retrieve AntiKt4TruthJets. Exiting." );
  //    return EL::StatusCode::FAILURE;
  //   }
  //} else {
  //   Warning("execute()", "Did not find AntiKt4TruthJets.");
  //}

  xAOD::TruthParticleContainer::const_iterator tpItr = m_TruthParticleContainer->begin();
  xAOD::TruthParticleContainer::const_iterator tpItrE = m_TruthParticleContainer->end();

  //xAOD::JetContainer::const_iterator jetItr = m_AntiKt4TruthContainer->begin();
  //xAOD::JetContainer::const_iterator jetItrE = m_AntiKt4TruthContainer->end();

  //std::cout << "Seg fault is below this line " << std::endl;

  int ntaus=0,nbs=0;
  int nphotons=0;


  int HiggsDescendants[100]={0};
  //std::map<int,int> HiggsDescendants;
  int nDescendants = 0;

  TLorentzVector SMhiggs(0,0,0,0),
                 bquarks(0,0,0,0), 
                 taus(0,0,0,0), 
                 visTaus(0,0,0,0);

  std::vector<TLorentzVector> sortedTaus;

  int nhadtaus=0;
  //HiggsDescendants.clear();

  // Looping over the truth particles in the event
  for( ; tpItr != tpItrE; ++tpItr ) {

     double pT  = (*tpItr)->p4().Pt();

     if(!(*tpItr)->parent(0) || pT <1e-3) continue;

     // Set few varibles for being used below
     //int barcode = (*tpItr)->barcode();
     int pdgId  = abs((*tpItr)->pdgId());
     int status = (*tpItr)->status();
     int ishadTau = -1;
     double mass   = (*tpItr)->p4().M();
     double eta = (*tpItr)->p4().Eta();
     double phi = (*tpItr)->p4().Phi();

     TLorentzVector Part4V(0,0,0,0);
     Part4V.SetPtEtaPhiM(pT,eta,phi,mass);
     sortedTaus.clear();

     int parentpdgId = abs((*tpItr)->parent(0)->pdgId());
     //int parentStat = abs((*tpItr)->parent(0)->status());
     //int parentsize = (*tpItr)->nParents();
     int nchilds = (*tpItr)->nChildren();

     // Write all particles in the 1st event
     if(m_eventCounter == 1) std::cout << " pdgId : "    << pdgId 
                                       << " Status: "    << status
                                       << " Mass: "      << mass 
                                       << " PT: "        << pT
                                       << " Parent: "    << parentpdgId 
                                       //<< " N parents: " << parentsize 
                                       << " N childs: "  << nchilds
                                       << " Barcode: "   << (*tpItr)->barcode()
                                       << " Parent: "    << (*tpItr)->parent(0)->barcode()
                                       << std::endl;

     if(pdgId == 39 && status == 62) h_graviton_mass->Fill(mass*1e-3);

     // Save Higgs pT
     if(abs(pdgId)==25) {
         h_higgs_pT->Fill(pT*1e-3);

         //Loop over the Higgs Cildrens
         for(int nch=0;nch<nchilds;nch++){
            int chId = abs((*tpItr)->child(nch)->pdgId());
            //Skip if the child is not a tau
            if(chId!=15) continue;
            //std::cout << "Higgs child: " << (*tpItr)->child(nch)->pdgId() << " chId " <<chId<< " barcode " << (*tpItr)->child(nch)->barcode() << std::endl;
            HiggsDescendants[nDescendants]=(*tpItr)->child(nch)->barcode();
            //HiggsDescendants[nDescendants]=(*tpItr)->child(nch)->pdgId();
            //std::cout << "Recording descendants " << HiggsDescendants[nDescendants] <<std::endl;
            nDescendants++;
         }
     }

     double isHiggsDescendant=0;

     if(parentpdgId !=15 && pdgId==15 && parentpdgId!=25) {
        //HiggsDescendants[nDescendants]=(*tpItr)->child(0)->barcode();
        //std::cout << "b-child: " << (*tpItr)->pdgId() << " barcode " << (*tpItr)->barcode() << " parent: " << parentpdgId << " pT: " << pT << std::endl;
        ntaus++;
        //std::cout << "Adding further descendants tau->tau decay " << HiggsDescendants[nDescendants] << std::endl;
     }

     // Only Higgs/Tau/b childs from this point
     if(parentpdgId == 25 || parentpdgId == 5 || parentpdgId == 15){ 

        // This is needed for Herwig events
        if(parentpdgId==15 && nchilds==1) {
          HiggsDescendants[nDescendants]=(*tpItr)->child(0)->barcode();
          //std::cout << "Adding further descendants tau->tau decay " << HiggsDescendants[nDescendants] << std::endl;
        	 nDescendants++;
        }

        //std::cout << "nDescendants: " << nDescendants << std::endl;
        for(int i=0;i<nDescendants;i++){
           //Check if the particle is a Higgs descendant and save their childrens
           /*std::cout << " part barcode "     << (*tpItr)->barcode() 
                     << " pdgId "            << (*tpItr)->pdgId()
                     << " status "            << (*tpItr)->status()
                     << " nDesc "            << nDescendants 
                     << " isHiggs Desc ? "   << isHiggsDescendant
                     << " N childs "         << nchilds << std::endl;*/
           //if(HiggsDescendants[i]==(*tpItr)->pdgId()){
           if((*tpItr)->barcode() == HiggsDescendants[i]) {
             // If particle barcode is equal to what is stored on HiggsDescendants
             isHiggsDescendant = 1;
             //std::cout << " is Higgs Descendant ? " << HiggsDescendants[i] << " pdgId: " << (*tpItr)->pdgId() << " parent: " << abs((*tpItr)->parent(0)->pdgId()) << std::endl;
             if(pdgId == 22 && abs((*tpItr)->parent(0)->pdgId())==15) {nphotons++; h_pho_taurad_pT->Fill(pT*1e-3); /*std::cout << "Summing photons: " << nphotons << std::endl;*/} //Do not add further gg pi0 decays
             for(int nch=0;nch<nchilds;nch++){
                HiggsDescendants[nDescendants]=(*tpItr)->child(nch)->pdgId();
                //std::cout << "Childs pdgId: " << (*tpItr)->child(nch)->pdgId() << " " << (*tpItr)->child(nch)->barcode() << std::endl;
                HiggsDescendants[nDescendants]=(*tpItr)->child(nch)->barcode();
                //std::cout << "Adding further descendants " << HiggsDescendants[nDescendants] << std::endl;
                nDescendants++;
             }
           }
        }


        TLorentzVector visTau(0,0,0,0);
        int flavour = -1;

        //std::cout << "Is Higgs descendant: ? " << isHiggsDescendant << std::endl;
        for(Int_t j=0;j<nDescendants;j++){
          //std::cout << "Looping over HiggsDescendants: " << HiggsDescendants[j] << std::endl;
        }

        // There should be only two of these taus by event
        if(abs(pdgId)==15 && status==2 && (parentpdgId==15 || parentpdgId==25) && isHiggsDescendant){
          const xAOD::TruthParticle *tau = (*tpItr);    
          ishadTau = IsHadronicTau(tau,visTau,flavour);
          /*std::cout << "Is hadronic tau ? " << ishadTau 
                    << " vis tau pT " << visTau.Pt() 
                    << " full tau pT " << pT 
                    << " flavour " << flavour  << std::endl;
                    
          std::cout << "Tau from Higgs decay: " << pT 
                    << " with status " << status 
                    << " is a Higgs Descendant ? " << isHiggsDescendant
                    << std::endl; */
          if(parentpdgId==15) h_tau_taurad_pT->Fill(pT*1e-3);
          else if(parentpdgId==25) h_tau_notaurad_pT->Fill(pT*1e-3);
          taus += Part4V;
          visTaus += visTau;
          sortedTaus.push_back(visTau);
          h_tau_pT->Fill(pT*1e-3); 
          h_tau_eta->Fill(eta);
          ntausall++;
          if(ishadTau==1){
             h_tauh_pT->Fill(visTau.Pt()*1e-3); 
             h_tauh_eta->Fill(visTau.Eta());
             nhadtaus++;
             nhtaus++;
          }
          if(ishadTau==0){
             if(flavour==11) {
                 h_el_pT->Fill(visTau.Pt()*1e-3);
                 h_el_eta->Fill(visTau.Eta()); 
                 neltaus++;
             }
             if(flavour==13) {
                 h_mu_pT->Fill(visTau.Pt()*1e-3);
                 h_mu_eta->Fill(visTau.Eta());
                 nmutaus++;
             }
          }
          ntaus++;
        }

        if(abs(pdgId) == 5 && (status == 23 || status == 11) && parentpdgId==25) {
          std::cout << "b from Higgs decay: " << pT << " with status " << status << std::endl;
          bquarks += Part4V; 
          h_bquark_pT->Fill(pT*1e-3);
          h_bquark_eta->Fill(eta);
          nbs++;
        }
     }
     //if(abs(pdgId) == 25 || abs(pdgId)==3 || abs(pdgId)==15 || abs(pdgId)==93){

  }

  
  //for( ; jetItr != jetItrE; ++jetItr ){
  //      std::cout << "Jet pT:= " << (*jetItr)->pt() << std::endl;
  //}

  sortVectorByPt(sortedTaus);

  h_npho_taurad->Fill(nphotons);
  h_ntaus_higgs->Fill(ntaus);

  std::cout << "Ntaus: " << ntaus << "NBs " << nbs << " NHadTaus " << nhadtaus << std::endl;

  std::cout << "TauTau pT: " << taus.Pt() << std::endl;
  std::cout << "BB pT: " << bquarks.Pt() << std::endl;

  h_bb_pT->Fill(bquarks.Pt()*1e-3);

  h_tautau_hh_mass->Fill(taus.M()*1e-3);
  h_tautau_hh_vis_mass->Fill(visTaus.M()*1e-3);
  h_tautau_pT->Fill(taus.Pt()*1e-3);
  h_tautau_vis_pT->Fill(visTaus.Pt()*1e-3);

  SMhiggs += bquarks;
  SMhiggs += taus;
  
  //std::cout << "mhh: " << SMhiggs.M() << std::endl;
  h_hh_mass->Fill(SMhiggs.M()*1e-3);

  //Select offline taus which satisfy the recommended properties
  //std::vector<const xAOD::TruthParticle*> truthTaus;
  //std::vector<const xAOD::TruthParticle*> truthTaus;


  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyEventLoop :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyEventLoop :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  //std::cout << "seg fault is below this line " << std::endl;

  std::cout << " tau->el fraction " << (float)neltaus/ntausall <<  std::endl;
  std::cout << " tau->mu fraction " << (float)nmutaus/ntausall << std::endl;
  std::cout << " tau->hadrons fraction " << (float)nhtaus/ntausall << std::endl;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyEventLoop :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}


int MyEventLoop::IsHadronicTau(const xAOD::TruthParticle *part,TLorentzVector &visTau, int &flavour){

    std::vector<const xAOD::TruthParticle*> m_truetaus_had;
    std::vector<const xAOD::TruthParticle*> m_truetaus_lep;
    m_truetaus_had.clear();
    m_truetaus_lep.clear();
    flavour = -1;
    if( fabs(part->pdgId()) != 15 ) return -1;
	    if(part->status() ==3 ) return -1;
       bool ishad=true;
       TLorentzVector invisTau(0,0,0,0);

	    if( fabs(part->pdgId()) == 15){
	       // std::cout<<"----- "<<part->nChildren() <<std::endl;
	        for(unsigned int iChild=0; iChild< part->nChildren(); ++iChild){
	            xAOD::TruthParticle* child = ( xAOD::TruthParticle*) part->child(iChild);

	            //std::cout<<"Child status:pdgid"<<child->status()  <<child->pdgId()  <<std::endl;

	            if( child->status() == 3 ) continue;

               // Sum the neutrinos 4Vector in the tau decay
               if( fabs( child->pdgId() ) == 12 || fabs( child->pdgId() ) == 14 || fabs( child->pdgId() ) == 16 ){
                invisTau += child->p4();
               }

               // Check if the tau decay to leptons
	            if( fabs( child->pdgId() ) == 11 || fabs( child->pdgId() ) == 13 ){
                ishad=false;
                flavour = abs(child->pdgId());
                //break;
	            } 
               

	        }

           //std::cout << " Invisible pT " << invisTau.Pt() << std::endl; 
           visTau = part->p4() - invisTau;

    }//taus

    if (ishad) {
        flavour = 15;
        return 1;
    }
    else {
        return 0;
    }

}

bool TLorentzVectorLessThan(const TLorentzVector &a, const TLorentzVector &b) {
   return a.Perp() > b.Perp(); // inverted to invert sorting
}

void MyEventLoop::sortVectorByPt(std::vector<TLorentzVector> &vect)
{
   std::sort(vect.begin(), vect.end(), TLorentzVectorLessThan);
   return;
}


//void MyEventLoop::TraceBackTauParent()
//{  
//  parent_index
//}
