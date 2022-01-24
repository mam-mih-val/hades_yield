//
// Created by mikhail on 11/30/20.
//

#include "super_event.h"

#include "TLorentzVector.h"
#include "TDatabasePDG.h"

#include <AnalysisTree/DataHeader.hpp>

TASK_IMPL(SuperEvent)

boost::program_options::options_description SuperEvent::GetBoostOptions() {
  using namespace boost::program_options;
  options_description desc(GetName() + " options");
  desc.add_options()
    ("pdg-code", value(&reference_pdg_code_)->default_value(2212), "PDG-code of particle")
    ("group-events", value(&ref_n_events_)->default_value(5), "Number of events in a group");
  return desc;
}

void SuperEvent::PreInit() {
//  SetInputBranchNames({"mdc_vtx_tracks", "event_header", "sim_tracks"});
}

void SuperEvent::UserInit(std::map<std::string, void *> &Map) {
  // extracting event variables
  event_header_ = GetInBranch("event_header");
  sim_header_ = GetInBranch("sim_header");
  tracks_ = GetInBranch("mdc_vtx_tracks");
  sim_particles_ = GetInBranch("sim_tracks");

  h1_rec_all_theta_distribution_ = new TH1F( "h1_rec_all_theta_distribution_", ";theta",
                                            6, 0.2, 1.4 );

  h1_rec_pid_theta_distribution_ = new TH1F( "h1_rec_pid_theta_distribution_", ";theta",
                                            6, 0.2, 1.4 );
  h1_tru_pid_theta_distribution_ = new TH1F( "h1_tru_pid_theta_distribution_", ";theta",
                                            6, 0.2, 1.4 );

  out_file_->cd();
  for( int i=0; i<6; ++i  ) {
    linear_fitters_.push_back( new TLinearFitter(6) );
    linear_fitters_.back()->SetFormula("hyp6");
  }
  std::cout << "Initialized" << std::endl;
}

void SuperEvent::UserExec() {
  this->LoopRecTracks();
  this->LoopTruParticles();
  n_events_++;
  if( n_events_ == ref_n_events_ ){
    std::vector<double> n_particles_theta;
    for( int theta_bin=1; theta_bin<=h1_rec_all_theta_distribution_->GetNbinsX(); theta_bin++ ) {
      auto n_part = h1_rec_all_theta_distribution_->GetBinContent(theta_bin);
      n_part/= (double) n_events_;
      n_particles_theta.emplace_back(n_part);
    }
    std::vector<double> efficiency_theta;
    std::vector<double> err_efficiency_theta;
    std::vector<std::vector<double>> unbiased_nparicles_theta;
    for( int theta_bin=1; theta_bin<= h1_rec_pid_theta_distribution_->GetNbinsX(); theta_bin++ ){
      auto N_rec =
          h1_rec_pid_theta_distribution_->GetBinContent(theta_bin);
      auto rec_err = h1_rec_pid_theta_distribution_->GetBinError(theta_bin);
      auto N_tru =
          h1_tru_pid_theta_distribution_->GetBinContent(theta_bin);
      auto tru_err = h1_tru_pid_theta_distribution_->GetBinError(theta_bin);
      auto eff = N_tru > 0 ? N_rec / N_tru : 0;
      auto err = sqrt( pow(rec_err/N_tru, 2) + pow(N_rec/N_tru/N_tru*tru_err, 2) );
      efficiency_theta.push_back(eff);
      err_efficiency_theta.push_back( err );
      unbiased_nparicles_theta.push_back( n_particles_theta );
      unbiased_nparicles_theta.back().at( theta_bin-1 ) -= N_rec;
    }
    for( int i=0; i<6; ++i  ) {
      linear_fitters_.at(i)->AddPoint(unbiased_nparicles_theta.at(i).data(),
                                efficiency_theta.at(i)*100,
                                err_efficiency_theta.at(i)*100);
    }

    h1_rec_all_theta_distribution_->Clear("ICESM");
    h1_rec_pid_theta_distribution_->Clear("ICESM");
    h1_tru_pid_theta_distribution_->Clear("ICESM");
    n_events_=0;
  }
}

void SuperEvent::LoopRecTracks() {
  using AnalysisTree::Particle;
  auto rec_chi2_var = GetVar("mdc_vtx_tracks/chi2");
  auto rec_dca_xy_var = GetVar("mdc_vtx_tracks/dca_xy");
  auto rec_dca_z_var = GetVar("mdc_vtx_tracks/dca_z");

  for (auto track : tracks_->Loop()) {
    auto pid = track.DataT<Particle>()->GetPid();
    auto mass = track.DataT<Particle>()->GetMass();
    if( pid != 0 ) {
      if( TDatabasePDG::Instance()->GetParticle(pid) )
        mass = TDatabasePDG::Instance()->GetParticle(pid)->Mass();
    }
    auto mom4 = track.DataT<Particle>()->Get4MomentumByMass(mass);
    auto chi2 = track[rec_chi2_var].GetVal();
    auto dca_xy = track[rec_dca_xy_var].GetVal();
    auto dca_z = track[rec_dca_z_var].GetVal();
    if( chi2 > 100.0 )
      continue;
    h1_rec_all_theta_distribution_->Fill(mom4.Theta());
    if ( -10 > dca_xy || dca_xy > 10 )
      continue;
    if ( -10 > dca_z || dca_z > 10 )
      continue;
    if( pid != reference_pdg_code_ )
      continue;
    if( mom4.Pt() > 0.4 )
      h1_rec_pid_theta_distribution_->Fill( mom4.Theta() );
  }
}

void SuperEvent::LoopTruParticles() {
  using AnalysisTree::Particle;
  auto var_is_primary = GetVar("sim_tracks/is_primary");

  for( auto particle : sim_particles_->Loop() ){
    auto mass = particle.DataT<Particle>()->GetMass();
    auto pid = particle.DataT<Particle>()->GetPid();
    auto is_prim = particle[var_is_primary].GetBool();
    auto mom4 = particle.DataT<Particle>()->Get4MomentumByMass(mass);
    double charge=0.0;
    if( TDatabasePDG::Instance()->GetParticle( pid ) ){
      charge= TDatabasePDG::Instance()->GetParticle( pid )->Charge() / 3.0;
    }
    if( fabs(charge) < 0.01 )
      continue;
    if( !is_prim )
      continue;
    if( pid!=reference_pdg_code_ )
      continue;
    if( mom4.Pt() > 0.4 )
      h1_tru_pid_theta_distribution_->Fill( mom4.Theta() );
  }
}


void SuperEvent::UserFinish() {
  out_file_->cd();
  for( int i=0; i<6; ++i  ) {
    linear_fitters_.at(i)->Eval();
  }
  std::vector<double> v_pars;
  v_pars.reserve(6);
  std::cout << "Fit: " << std::endl;
  auto*p2_theta1_theta2 = new TH2F( "p2_theta1_theta2", ";#theta_1;#theta_2",
                                    6, 0.2, 1.4,
                                    6, 0.2, 1.4);
  auto*p1_efficiency_theta = new TH1F( "p1_efficiency_theta", ";#theta_1;#theta_2",
                                       6, 0.2, 1.4);
  for( int j=0; j<6; ++j ){
    for (int i = 0; i < 7; ++i) {
      auto par = linear_fitters_.at(j)->GetParameter(i);
      auto par_err = linear_fitters_.at(j)->GetParError(i);
      v_pars.push_back(linear_fitters_.at(j)->GetParameter(i));
      if( i == 0 ){
        p1_efficiency_theta->SetBinContent(j+1,par);
        p1_efficiency_theta->SetBinError( j+1, par_err);
      }
      if( i > 0 ) {
        p2_theta1_theta2->SetBinContent( i, j + 1, par);
        p2_theta1_theta2->SetBinError(i, j + 1, par_err);
      }
      std::cout << v_pars.back() << " ";
    }
    std::cout << std::endl;
  }
  p2_theta1_theta2->Write();
  p1_efficiency_theta->Write();
  std::cout << "Finished" << std::endl;
}
