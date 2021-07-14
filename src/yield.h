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
#include <AnalysisTree/Matching.hpp>
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
  bool is_mc_;

  ATI2::Branch* event_header_;
  ATI2::Branch* tracks_{nullptr};
  ATI2::Branch* meta_hits_{nullptr};
  ATI2::Branch* sim_particles_{nullptr};
  AnalysisTree::Matching* mdc2meta_matching_{nullptr};
  AnalysisTree::Matching* mdc2sim_matching_{nullptr};
  // event variables
  ATI2::Variable centrality_var_;
  // mdc variables
  ATI2::Variable charge_var_;
  ATI2::Variable pdg_code_var_;
  ATI2::Variable dca_xy_var_;
  ATI2::Variable dca_z_var_;
  ATI2::Variable chi2_var_;
  ATI2::Variable dedx_mdc_var_;
  // meta variables
  ATI2::Variable mass2_var_;
  ATI2::Variable beta_var_;
  ATI2::Variable dedx_meta_var_;
  ATI2::Variable is_rpc_hit_var_;
  // sim particles variables
  ATI2::Variable sim_pdg_code_var_;
  ATI2::Variable is_primary_var_;

  // All
  std::vector<TH2F*> m2_vs_pq_all_;
  std::vector<TH2F*> beta_vs_pq_all_;
  std::vector<TH2F*> dedx_mdc_vs_pq_all_;
  std::vector<TH2F*> dedx_meta_vs_pq_all_;
  std::vector<TH2F*> pt_eta_all_;
  std::vector<TH2F*> sim_pt_eta_all_;
  // PID Reco
  std::vector<TH2F*> m2_vs_pq_pid_reco_;
  std::vector<TH2F*> beta_vs_pq_pid_reco_;
  std::vector<TH2F*> dedx_mdc_vs_pq_pid_reco_;
  std::vector<TH2F*> dedx_meta_vs_pq_pid_reco_;
  // Mismatch
  std::vector<TH2F*> m2_vs_pq_mismatch_;
  std::vector<TH2F*> beta_vs_pq_mismatch_;
  std::vector<TH2F*> dedx_mdc_vs_pq_mismatch_;
  std::vector<TH2F*> dedx_meta_vs_pq_mismatch_;


  double beta_cm_;

  TH1F* centrality_classes_;
  TH1F* pt_distribution_reco_;
  TH1F* pt_distribution_sim_;
  TH2F* rapidity_true_mass_;
  TH3F* pT_rapidity_eta_;
  std::vector<TH3F*> pT_rapidity_eta_matricies_;

TASK_DEF(Yield, 0)
};

#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
