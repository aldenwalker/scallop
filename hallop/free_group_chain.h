#ifndef FREE_GROUP_CHAIN_H
#define FREE_GROUP_CHAIN_H

#include <vector>
#include <iostream>
#include <string>
 
 
 
namespace HALLOP {
  
/*****************************************************************************
 * a letter in a chain
 * ***************************************************************************/
struct ChainLetter {
  int word;
  int index;
  int index_in_group_reg_inv_list;  //index of the letter in the reg/inverse list
  char letter;
  int group;         //the index of the group in the free product
};

std::ostream &operator<<(std::ostream &os, HALLOP::ChainLetter &CL); 

 
 
 
 
 
/*****************************************************************************
 * A chain   
 * ***************************************************************************/
struct FreeGroupChain { 
  FreeGroupChain(void);
  FreeGroupChain(char** input, int num_strings);
  void add_word(int weight, std::string w);
  
  int next_letter(int n);
  int prev_letter(int n);
  int num_words(void);
  int num_letters();
  std::string operator[](int index);    //get a word
  void print_letters(std::ostream &os);
  void print_group_letters(std::ostream &os);

  int rank;
  std::vector<std::string> words;
  std::vector<int> weights;
  std::vector<std::vector<int> > group_letters;
  std::vector<ChainLetter> chain_letters;
  std::vector<std::vector<int> > regular_letters;
  std::vector<std::vector<int> > inverse_letters;
    
};

std::ostream &operator<<(std::ostream &os, HALLOP::FreeGroupChain &C);

}

#endif