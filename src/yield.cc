//
// Created by mikhail on 11/30/20.
//

#include "yield.h"

#include "TLorentzVector.h"
#include "TDatabasePDG.h"

#include <AnalysisTree/DataHeader.hpp>

TASK_IMPL(Yield)

boost::program_options::options_description Yield::GetBoostOptions() {
  using namespace boost::program_options;
  options_description desc(GetName() + " options");
  desc.add_options()
    ("pdg-code", value(&reference_pdg_code_)->default_value(2212), "PDG-code of particle");
  return desc;
}

void Yield::PreInit() {
//  SetInputBranchNames({"mdc_vtx_tracks", "event_header", "sim_tracks"});
}

void Yield::UserInit(std::map<std::string, void *> &Map) {
  // extracting event variables
  event_header_ = GetInBranch("event_header");
  tracks_ = GetInBranch("mdc_vtx_tracks");
  sim_particles_ = GetInBranch("sim_tracks");

  centrality_distribution_ = new TH1F( "centrality", ";TOF+RPC hits centrality (%)", 20, 0.0, 100.0 );

  rec_y_pT_centrality_ = new TH3F( "rec_y_pT_centrality", ";y_{cm};p_{T} [GeV/c];TOF+RPC hits centrality (%)",
                                  100, -0.85, 1.15,
                                  125, 0.0, 2.5,
                                  20, 0, 100);
  tru_y_pT_centrality_ = new TH3F( "tru_y_pT_centrality", ";y_{cm};p_{T} [GeV/c];TOF+RPC hits centrality (%)",
                                  100, -0.85, 1.15,
                                  125, 0.0, 2.5,
                                  20, 0, 100);
  std::vector<double> pt_axis;
  std::vector<double> m0_axis;

  if( abs(reference_pdg_code_) == 2212 ) {
    pt_axis = {0.0,     0.29375, 0.35625, 0.41875, 0.48125, 0.54375,
               0.61875, 0.70625, 0.81875, 1.01875, 2.0};
    for( int i=0; i<40; i+=2 ){m0_axis.push_back(i);}
  }
  if( abs(reference_pdg_code_) == 211 ) {
    pt_axis = {0,    0.08, 0.105, 0.13,  0.155, 0.18,
               0.21, 0.25, 0.315, 0.535, 1.0};
    for( int i=0; i<30; i+=2 ){m0_axis.push_back(i);}
  }

  rec_pT_multiplicity_midtrapidity_ = new TH1F( "rec_pT_multiplicity_midtrapidity",
                                               ";number of charged tracks",
                                               m0_axis.size()-1, m0_axis.data()
                                               );
  tru_pT_multiplicity_midtrapidity_ = new TH1F( "tru_pT_multiplicity_midtrapidity",
                                               ";p_{T} [GeV/c];number of charged tracks",
                                               m0_axis.size()-1, m0_axis.data()
                                               );
  out_file_->cd();
  auto y_cm = data_header_->GetBeamRapidity();
  std::vector<double> pT_midrapidity;
  double ref_mass = TDatabasePDG::Instance()->GetParticle( reference_pdg_code_ )->Mass();
  if( abs(reference_pdg_code_) == 2212 ){
    pT_midrapidity = {0.35625, 0.41875,};
  }
  if( abs(reference_pdg_code_) == 211 ){
    pT_midrapidity = {0.18, 0.21,};
  }
  theta_range_ = CalculateThetaRange( {0.69, 0.79}, pT_midrapidity, ref_mass  );
  beta_cm_ = tanh(y_cm);
  std::cout << "Initialized" << std::endl;
}

void Yield::UserExec() {
  using AnalysisTree::Particle;
  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();
  float y_beam = data_header_->GetBeamRapidity();
  double ref_mass = TDatabasePDG::Instance()->GetParticle( reference_pdg_code_ )->Mass();
  std::vector<double> pT_midrapidity;
  if( abs(reference_pdg_code_) == 2212 ){
    pT_midrapidity = {0.35625, 0.41875,};
  }
  if( abs(reference_pdg_code_) == 211 ){
    pT_midrapidity = {0.18, 0.21,};
  }
//  auto theta_range = CalculateThetaRange( {0.69, 0.79}, pT_midrapidity, ref_mass );
  auto rec_chi2_var = GetVar("mdc_vtx_tracks/chi2");
  auto rec_dca_xy_var = GetVar("mdc_vtx_tracks/dca_xy");
  auto rec_dca_z_var = GetVar("mdc_vtx_tracks/dca_z");
  auto tru_is_primary = GetVar("sim_tracks/is_primary");

  auto n_tracks_in_midrapidity=CalculateNumberOfChargedTracks(theta_range_);
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
    if( pid != reference_pdg_code_ )
      continue;
    if( chi2 > 100.0 )
      continue;
    if ( -10 > dca_xy || dca_xy > 10 )
      continue;
    if ( -10 > dca_z || dca_z > 10 )
      continue;
    auto y = mom4.Rapidity() - y_beam;
    rec_y_pT_centrality_->Fill( y, mom4.Pt(), centrality );
    if( -0.05 < y && y < 0.05 )
      if( pT_midrapidity.front() < mom4.Pt() && mom4.Pt() < pT_midrapidity.back()  )
        rec_pT_multiplicity_midtrapidity_->Fill( n_tracks_in_midrapidity );
  }
  for( auto particle : sim_particles_->Loop() ){
    auto mass = particle.DataT<Particle>()->GetMass();
    auto pid = particle.DataT<Particle>()->GetPid();
    auto is_prim = particle[tru_is_primary].GetBool();
    if( !is_prim )
      continue;
    if( pid!=reference_pdg_code_ )
      continue;
    auto mom4 = particle.DataT<Particle>()->Get4MomentumByMass(mass);
    tru_y_pT_centrality_->Fill( mom4.Rapidity() - y_beam, mom4.Pt(), centrality );
    if( -0.05 < mom4.Rapidity() - y_beam && mom4.Rapidity() - y_beam < 0.05 )
      if( pT_midrapidity.front() < mom4.Pt() && mom4.Pt() < pT_midrapidity.back()  )
        tru_pT_multiplicity_midtrapidity_->Fill( n_tracks_in_midrapidity );
  }

}
void Yield::UserFinish() {
  centrality_distribution_->Write();
  rec_pT_multiplicity_midtrapidity_->Write();
  tru_pT_multiplicity_midtrapidity_->Write();
  rec_y_pT_centrality_->Write();
  tru_y_pT_centrality_->Write();
  std::cout << "Finished" << std::endl;
}
int Yield::CalculateNumberOfChargedTracks(std::vector<double> theta_range) {
  using AnalysisTree::Particle;
  int n_particles=0;
  for( auto particle : sim_particles_->Loop() ){
    auto pid = particle.DataT<Particle>()->GetPid();
    double charge=0.0;
    if( TDatabasePDG::Instance()->GetParticle( pid ) ){
      charge= TDatabasePDG::Instance()->GetParticle( pid )->Charge() / 3.0;
    }
    if( fabs(charge) < 0.01 )
      continue;
    auto mass = particle.DataT<Particle>()->GetMass();
    auto mom4 = particle.DataT<Particle>()->Get4MomentumByMass(mass);
    if( theta_range.front() < mom4.Theta() && mom4.Theta() < theta_range.back() )
      n_particles++;
  }
  return n_particles;
}
