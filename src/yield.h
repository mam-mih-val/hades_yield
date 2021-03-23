//
// Created by mikhail on 11/30/20.
//

#ifndef HADES_RAPIDITY_SRC_RAPIDITY_H_
#define HADES_RAPIDITY_SRC_RAPIDITY_H_

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>

#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/EventHeader.hpp>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <at_task/Task.h>
#include <memory>
#include <string>

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
  ATI2::Branch* event_header_;
  ATI2::Branch* tracks_{nullptr};
  ATI2::Branch* sim_particles_{nullptr};
  ATI2::Variable centrality_;
  ATI2::Variable pdg_code_;
  ATI2::Variable dca_xy_;
  ATI2::Variable dca_z_;
  ATI2::Variable chi2_;

  ATI2::Variable sim_pdg_code_;
  ATI2::Variable is_primary_;

  std::vector<TH2F*> yields_;
  TH1F* pt_distribution_reco_;
  TH1F* pt_distribution_sim_;
  TH2F* rapidity_true_mass_;
  TH1F* centrality_classes_;

  TH2F* ecm_pT_y_protons_;
  TH2F* ecm_pT_y_pions_;

TASK_DEF(Yield, 0)
};

#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
