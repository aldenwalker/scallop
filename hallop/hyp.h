#ifndef HALLOP_HYP_H
#define HALLOP_HYP_H

#include <vector>
#include <string>

//a hyperbolic group 
//it assumes that the given relators are a Dehn
//presentation

namespace HALLOP {

typedef int SignedInd;

struct HyperbolicGroup {
  int rank;
  vector<vector<HALLOP::SignedInd> > relators;
  
  HyperbolicGroup();
  HyperbolicGroup(int rank);
  HyperbolicGroup(int rank, std::string relator);
  HyperbolicGroup(int rank, std::vector<string>& realtors);
  
  SignedInd letter_to_signed_ind(char c);
  char signed_ind_to_letter(HALLOP::SignedInd);
  vector<HALLOP::SignedInd> word_to_vector(std::string& s);
  std::string vector_to_word(vector<HALLOP::SignedInd>& w);
  
  
  void add_relator(std::string r);
  void make_geodesic(std::vector<HALLOP::SignedInd>& w);
  std::string geodesic(std::string& w);
  
  //all_deformations
  
};

}
  
#endif