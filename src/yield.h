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
//    for( auto theta : theta_vector ){
//      std::cout << theta << " ";
//    }
//    std::cout << std::endl;
    auto theta_min = *std::min_element( theta_vector.begin(), theta_vector.end() );
    auto theta_max = *std::max_element( theta_vector.begin(), theta_vector.end() );
    return {theta_min, theta_max};
  }
  int CalculateNumberOfChargedTracks( std::vector<double> );
private:
  bool is_mc_;

  ATI2::Branch* event_header_;
  ATI2::Branch* tracks_{nullptr};
  ATI2::Branch* sim_particles_{nullptr};

  TH1F* centrality_distribution_;
  TH3F* rec_y_pT_centrality_;
  TH3F* tru_y_pT_centrality_;
  TH2F* rec_pT_multiplicity_midtrapidity_;
  TH2F* tru_pT_multiplicity_midtrapidity_;

  double beta_cm_;
  double ref_mass_;
  std::vector<double> theta_range_;
  int reference_pdg_code_;

TASK_DEF(Yield, 0)
};

#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
