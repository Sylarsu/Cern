#ifndef ValidationHHtobbtt_MyEventLoop_H
#define ValidationHHtobbtt_MyEventLoop_H

#include <EventLoop/Algorithm.h>

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

#include "TH1.h"
#include "TLorentzVector.h"

#include <iostream>

namespace xAOD{
  class TruthParticle_v1;
  typedef TruthParticle_v1 TruthParticle;
  class Jets_v1;
  typedef Jets_v1 Jets;
}


class MyEventLoop : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!
  xAOD::TEvent *m_event;  //!
 
  int m_eventCounter; //!

  TH1* h_graviton_mass;  //!
  TH1* h_hh_mass;  //!
  TH1* h_higgs_pT; //!
  TH1* h_bquark_pT; //!
  TH1* h_bquark_eta; //!
  TH1* h_tauh_pT; //!
  TH1* h_tauh_eta; //!
  TH1* h_tau_pT; //!
  TH1* h_tau_eta; //!
  TH1* h_tautau_hh_mass; //!
  TH1* h_bb_pT; //!
  TH1* h_tautau_pT; //!
  TH1* h_tautau_vis_pT; //!
  TH1* h_el_pT; //!
  TH1* h_el_eta; //!
  TH1* h_mu_pT; //!
  TH1* h_mu_eta; //!
  TH1* h_tautau_hh_vis_mass; //! 
  TH1* h_ntaus_higgs; //!
  TH1* h_npho_taurad; //!
  TH1* h_pho_taurad_pT; //!
  TH1* h_tau_taurad_pT; //!
  TH1* h_tau_notaurad_pT; //!

  int nhtaus; //!
  int ntausall; //!
  int neltaus; //!
  int nmutaus; //!

  // this is a standard constructor
  MyEventLoop ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  int IsHadronicTau(const xAOD::TruthParticle *, TLorentzVector &visTau, int &flavour);
  void sortVectorByPt(std::vector<TLorentzVector> &vect);

  // this is needed to distribute the algorithm to the workers
  ClassDef(MyEventLoop, 1);
};

#endif
