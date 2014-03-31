#include "hyp.h"

HALLOP::HyperbolicGroup::HyperbolicGroup() {
  rank = 0;
  relators.resize(0);
}

HALLOP::HyperbolicGroup::HyperbolicGroup(int rank) {
  this->rank = rank;
  relators.resize(0);
}

HALLOP::HyperbolicGroup::HyperbolicGroup(int rank, std::string relator) {
  this->rank = rank;
  relators.resize(1);
  relators[0] = this->word_to_vector(relator);
}

HALLOP::HyperbolicGroup::HyperbolicGroup(int rank, 
                                         std::vector<string>& realtors) {
  this->rank = rank;
  this->relators.resize(relators.size());
  for (int i=0; i<(int)relators.size(); ++i) {
    this->relators[i] = word_to_vector(relators[i]);
  }
}

  
  void add_relator(std::string r);
  void make_geodesic(std::vector<HALLOP::SignedInd>& w);
  std::string geodesic(std::string& w);
  
  
  
  SignedInd letter_to_signed_ind(char c);
  char signed_ind_to_letter(HALLOP::SignedInd);
  vector<HALLOP::SignedInd> word_to_vector(std::string& s);
  std::string vector_to_word(vector<HALLOP::SignedInd>& w);
  