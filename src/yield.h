//
// Created by mikhail on 11/30/20.
//

#ifndef HADES_RAPIDITY_SRC_RAPIDITY_H_
#define HADES_RAPIDITY_SRC_RAPIDITY_H_

#include <memory>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TProfile2D.h>
#include <TProfile3D.h>

#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/Matching.hpp>

#include <at_task/Task.h>

class Yield : public UserFillTask {

public:
  void UserInit(std::map<std::string, void *> &Map) override;
  void UserExec() override;
  void UserFinish() override;
  boost::program_options::options_description GetBoostOptions() override;
  void PreInit() override;
  void PostFinish() override {
    UserTask::PostFinish();
  }
private:
  double AngleDifference( double phi, double psi ){
    auto delta_phi = phi - psi;
    if( delta_phi < -M_PI )
      delta_phi+=2*M_PI;
    if( delta_phi > M_PI )
      delta_phi-=2*M_PI;
    return delta_phi;
  }
  void LoopRecTracks();
  void LoopTruParticles();

  ATI2::Branch* event_header_{nullptr};
  ATI2::Branch* sim_header_{nullptr};
  ATI2::Branch* tracks_{nullptr};
  ATI2::Branch* sim_particles_{nullptr};

  AnalysisTree::Matching* rec_sim_matching_{nullptr};

  TH1F* h1_centrality_;

  // True gen particles
  TH3F* h3_tru_delta_phi_theta_centrality_pid_{nullptr};
  TH3F* h3_tru_delta_phi_theta_centrality_all_{nullptr};
  TH2F* h2_rec_theta_centrality_all_{nullptr};
  TProfile2D* p2_tru_v1_pid_{nullptr};
  TProfile2D* p2_rec_v1_pid_{nullptr};

  // Reconstructed particles
  TH3F* h3_rec_delta_phi_theta_centrality_pid_{nullptr};
  TH3F* h3_rec_delta_phi_theta_centrality_all_{nullptr};
  TH2F* h2_tru_theta_centrality_all_{nullptr};
  TProfile2D* p2_tru_v1_all_{nullptr};
  TProfile2D* p2_rec_v1_all_{nullptr};

  TProfile3D* p3_dtheta_dphi_dpT_loss_{nullptr};
  TH2F* h2_theta_phi_sector_lost_population_{nullptr};

  TProfile2D*p2_dphi_dtheta_efficiency_{nullptr};
  TH2F* h2_dphi_dtheta_{nullptr};

  double beta_cm_;
  double ref_mass_;

  int reference_pdg_code_;

TASK_DEF(Yield, 0)
};

#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
