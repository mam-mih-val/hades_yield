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
#include <TLinearFitter.h>

class SuperEvent : public UserFillTask {

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
  void LoopRecTracks();
  void LoopTruParticles();

  ATI2::Branch* event_header_{nullptr};
  ATI2::Branch* sim_header_{nullptr};
  ATI2::Branch* tracks_{nullptr};
  ATI2::Branch* sim_particles_{nullptr};

  TH1F* h1_rec_all_theta_distribution_;

  TH1F* h1_rec_pid_theta_distribution_;
  TH1F* h1_tru_pid_theta_distribution_;

  std::vector<TLinearFitter*> linear_fitters_;


  int reference_pdg_code_;
  int ref_n_events_;
  int n_events_{0};

TASK_DEF(SuperEvent, 0)
};

#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
