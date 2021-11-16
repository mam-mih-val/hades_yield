//
// Created by mikhail on 11/30/20.
//

#ifndef HADES_RAPIDITY_SRC_RAPIDITY_H_
#define HADES_RAPIDITY_SRC_RAPIDITY_H_

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TProfile3D.h>
#include <TF1.h>

#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/Matching.hpp>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
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
  std::vector<double> CalculateThetaRange( const std::vector<double>& y_range, const std::vector<double>& pT_range, double m ){
    std::vector<double> theta_vector;
    for( auto y : y_range ){
      for( auto pT : pT_range ){
        auto sinh_eta = sqrt( m*m + pT*pT )/pT * sinh( y );
        auto eta = asinh( sinh_eta );
        auto theta = 2 * atan( exp(-eta) );
        theta_vector.push_back(theta);
      }
    }
    auto theta_min = *std::min_element( theta_vector.begin(), theta_vector.end() );
    auto theta_max = *std::max_element( theta_vector.begin(), theta_vector.end() );
    return {theta_min, theta_max};
  }
  int CalculateNumberOfChargedTracks( std::vector<double> );
private:
  void InitEfficiency();

  ATI2::Branch* event_header_;
  ATI2::Branch* sim_header_;
  ATI2::Branch* tracks_{nullptr};
  ATI2::Branch* sim_particles_{nullptr};

  TH1F* h1_centrality_;
  TH2F* h1_phi_theta_all_;
  TH3F* h3_tru_theta_pT_phi_;
  TH3F* h3_rec_theta_pT_phi_;
  TProfile3D* p3_theta_pT_npart_sector_;

  double beta_cm_;
  double ref_mass_;

  int reference_pdg_code_;

TASK_DEF(Yield, 0)
};

#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
