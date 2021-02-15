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
#include <at_task/Task.h>
#include <memory>
#include <string>

class Yield : public UserFillTask {

public:
  void Init(std::map<std::string, void *> &Map) override;
  void Exec() override;
  void Finish() override;
  boost::program_options::options_description GetBoostOptions() override;
  void PreInit() override;
  void PostFinish() override {
    UserTask::PostFinish();
  }

private:
  int pdg_code_;
  std::string tracks_branch_;
  AnalysisTree::Particles *tracks_{nullptr};
  AnalysisTree::BranchConfig tracks_config_;
  AnalysisTree::EventHeader *event_header_{nullptr};
  AnalysisTree::BranchConfig event_header_config_;

  std::vector<TH2F*> yields_;
  TH1F* centrality_classes_;
  TH1F* vtx_z_;

TASK_DEF(Yield, 0)
};

#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
