#include <vector>
#include <string>
#include <iostream>
#include <cctype>
#include <cstdlib>



#include "free_group_chain.h"

#include "../word.h"

 
/****************************************************************************
 print a chain letter
 ***************************************************************************/
std::ostream& HALLOP::operator<<(std::ostream &os, HALLOP::ChainLetter &CL) {
  os << CL.letter << " (" << CL.word << "," << CL.index << ","
     << CL.index_in_group_reg_inv_list << ")";
  return os;
}
   
 
 
/***************************************************************************
 * Chain in a free group
 * *************************************************************************/
HALLOP::FreeGroupChain::FreeGroupChain(void) {
  words.resize(0);
  weights.resize(0);
}


HALLOP::FreeGroupChain::FreeGroupChain(char** input, int num_strings) {
  int i;
  int j;
  std::string word;
  std::string weight;
  
  std::vector<std::string> temp_words(0);
  std::vector<int> temp_weights(0);
  //input the raw words
  for (i=0; i<num_strings; i++) {
    j=0;
    while (isdigit(input[i][j]) || ispunct(input[i][j])) {
      j++;
    }
    weight = std::string(input[i]).substr(0,j);
    word = std::string(&input[i][j]);
    if (weight == "") {
      temp_weights.push_back(1);
    } else {
      temp_weights.push_back( atoi( weight.c_str() ));
    }
    temp_words.push_back(word);
  }
  
  chain_letters.resize(0);
  group_letters.resize(0);
  regular_letters.resize(0);
  inverse_letters.resize(0);
  words.resize(0);
  weights.resize(0);
  word_start_indices.resize(0);
  rank = 0;
  
  for (int i=0; i<(int)temp_words.size(); ++i) {
    add_word(temp_weights[i], temp_words[i]);
  }
  
}

void HALLOP::FreeGroupChain::add_word(int weight, std::string w) {
  
  words.push_back(w);
  weights.push_back(weight);
  word_start_indices.push_back(chain_letters.size());
  
  //compute the rank of the word
  int new_rank = chain_rank(w);
  
  if (new_rank > rank) {
    group_letters.resize(new_rank, std::vector<int>());
    regular_letters.resize(new_rank, std::vector<int>());
    inverse_letters.resize(new_rank, std::vector<int>());
    rank = new_rank;
  }
  
  //add the word to the chain
  HALLOP::ChainLetter temp_letter;
  temp_letter.word = (int)words.size()-1;
  temp_letter.inverse_relator_letter = -1;
  for (int j=0; j<(int)w.size(); j++) {
    temp_letter.index = j;
    temp_letter.letter = w[j];
    temp_letter.group = letter_index( temp_letter.letter );
    group_letters[ temp_letter.group ].push_back(chain_letters.size());
    if (isupper(temp_letter.letter)) {
      inverse_letters[temp_letter.group].push_back(chain_letters.size());
      temp_letter.index_in_group_reg_inv_list 
                            = inverse_letters[temp_letter.group].size()-1;
    } else {
      regular_letters[temp_letter.group].push_back(chain_letters.size());
      temp_letter.index_in_group_reg_inv_list 
                            = regular_letters[temp_letter.group].size()-1;
    }
    chain_letters.push_back(temp_letter);
  }
}


void HALLOP::FreeGroupChain::add_relator(std::string w) {
  
  //add the relator and the inverse
  add_word(-1, w);
  add_word(-1, inverse(w));

  //go back and fix which letters cannot be paired together
  int wi = num_words() - 2;
  int wIi = num_words() - 1;
  int w_chain_letter_start = word_start_index(wi);
  int wI_chain_letter_start = word_start_index(wIi);
  


  words.push_back(w);
  weights.push_back(-1);
  words.push_back(inverse(w));
  weights.push_back(-1);
  
  //compute the rank of the word
  int new_rank = chain_rank(w);
  
  if (new_rank > rank) {
    group_letters.resize(new_rank, std::vector<int>());
    regular_letters.resize(new_rank, std::vector<int>());
    inverse_letters.resize(new_rank, std::vector<int>());
    rank = new_rank;
  }
  
  //add the relator to the chain
  HALLOP::ChainLetter temp_letter;
  temp_letter.word = (int)words.size()-2;
  for (int j=0; j<(int)w.size(); j++) {
    temp_letter.index = j;
    temp_letter.letter = w[j];
    temp_letter.group = letter_index( temp_letter.letter );
    group_letters[ temp_letter.group ].push_back(chain_letters.size());
    if (isupper(temp_letter.letter)) {
      inverse_letters[temp_letter.group].push_back(chain_letters.size());
      temp_letter.index_in_group_reg_inv_list 
                            = inverse_letters[temp_letter.group].size()-1;
    } else {
      regular_letters[temp_letter.group].push_back(chain_letters.size());
      temp_letter.index_in_group_reg_inv_list 
                            = regular_letters[temp_letter.group].size()-1;
    }
    chain_letters.push_back(temp_letter);
  }

  //add the inverse relator to the chain

}



int HALLOP::FreeGroupChain::next_letter(int n) {
  int word = chain_letters[n].word;
  int index = chain_letters[n].index;
  if (index == (int)words[word].size()-1) {
    return n-words[word].size()+1;
  } else {
    return n+1;
  }
}

int HALLOP::FreeGroupChain::prev_letter(int n) {
  int word = chain_letters[n].word;
  int index = chain_letters[n].index;
  if (index == 0) {
    return n+words[word].size()-1;
  } else {
    return n-1;
  }
}

int HALLOP::FreeGroupChain::word_start_index(int n) {
  return word_start_indices[n];
}

int HALLOP::FreeGroupChain::num_words(void) {
  return words.size();
}

int HALLOP::FreeGroupChain::num_letters() {
  return chain_letters.size();
}

int HALLOP::FreeGroupChain::rnk() {
  return rank;
}

std::string HALLOP::FreeGroupChain::operator[](int index) {
  return words[index];
}


void HALLOP::FreeGroupChain::print_letters(std::ostream &os) {
  int i;
  for (i=0; i<(int)chain_letters.size(); i++) {
    os << i << ": " << chain_letters[i] << "\n";
  }
}

void HALLOP::FreeGroupChain::print_group_letters(std::ostream &os) {
  int i,j;
  for (i=0; i<rank; i++) {
    os << "Group " << i << ":\n";
    for (j=0; j<(int)group_letters[i].size(); j++) {
      os << "(" << chain_letters[group_letters[i][j]].word << "," 
         << chain_letters[group_letters[i][j]].index << "," 
         << chain_letters[group_letters[i][j]].letter << "),";
    }
    os << "\n";
  }
}


std::ostream& HALLOP::operator<<(std::ostream &os, FreeGroupChain &C) {
  int i;
  int len = (int)C.words.size();
  for (i=0; i<len; i++) {
    os << C.weights[i] << C.words[i] << " ";
  }
  return os;
}

