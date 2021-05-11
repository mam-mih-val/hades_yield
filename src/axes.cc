//
// Created by mikhail on 3/15/21.
//

#include <TFile.h>
#include <TH1F.h>
#include <stdexcept>
#include <string>
int main(int n_args, char** args){
  std::string in_file_name;
  if( n_args > 1 )
    in_file_name = args[1];
  else
    throw std::runtime_error( "minimum 1 argument expected: name of the input file" );
  auto in_file = TFile::Open(in_file_name.c_str());
  TH1F* pt_distribution;
  in_file->GetObject("pT_distribution", pt_distribution);
  auto uniform = dynamic_cast<TH1F*>(pt_distribution->Clone( "pT_uniform" ));
  auto non_uniform = dynamic_cast<TH1F*>(pt_distribution->Clone( "pT_non_uniform" ));
  auto out_file = TFile::Open( "output.root", "recreate" );

  std::vector<float> uniform_axis;
  float pT=0;
  while ( pT < 1.8 ){
    uniform_axis.push_back(pT);
    pT+=0.1;
  }
  std::vector<float> non_uniform_axis{0, 0.29375, 0.35625, 0.41875, 0.48125, 0.54375, 0.61875, 0.70625, 0.81875, 1.01875, 2.0};
  TH1F* histo;
  out_file->mkdir( "uniform" );
  out_file->cd( "uniform" );

  for( int i=0; i<uniform_axis.size()-1; ++i ){
    std::string name = "pT_uniform_"+std::to_string(uniform_axis.at(i))+"-"+std::to_string(uniform_axis.at(i+1));
    histo = dynamic_cast<TH1F*>(pt_distribution->Clone( name.c_str() ));
    histo->GetXaxis()->SetRangeUser(uniform_axis.at(i), uniform_axis.at(i+1));
    histo->Write();
  }
  out_file->cd( "/" );
  out_file->mkdir( "non-uniform" );
  out_file->cd( "non-uniform" );
  for( int i=0; i<non_uniform_axis.size()-1; ++i ){
    std::string name = "pT_non_uniform_"+std::to_string(non_uniform_axis.at(i))+"-"+std::to_string(non_uniform_axis.at(i+1));
    histo = dynamic_cast<TH1F*>(pt_distribution->Clone( name.c_str() ));
    histo->GetXaxis()->SetRangeUser(non_uniform_axis.at(i), non_uniform_axis.at(i+1));
    histo->Write();
  }

  out_file->Close();
}