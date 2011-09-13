#include <vector>
#include <iostream>
#include <string>

#include <ctype.h>
#include <stdlib.h>

#include "scyllop_classes.h"


/*****************************************************************************
* A free product of cyclic groups
* ****************************************************************************/
CyclicProduct::CyclicProduct(void) {
  gens.resize(0);
  orders.resize(0);
}

CyclicProduct::CyclicProduct(std::string input) {
  int i;
  int start;
  gens.resize(0);
  orders.resize(0);
  i=0;
  while (i < (int)input.size()) { 
    gens.push_back(input[i]);
    i++;
    start = i;
    while (i < (int)input.size() && isdigit(input[i])) {
      i++;
    }
    orders.push_back( atoi( input.substr(start, i-start).c_str() ));
  }
}

CyclicProduct::~CyclicProduct(void) {
  //I don't need to free the vectors, right?
}


std::vector<char> CyclicProduct::gen_list(void) {
  std::vector<char> ans = gens;
  return ans;
}

std::vector<int> CyclicProduct::order_list(void) {
  std::vector<int> ans = orders;
  return ans;
}

int CyclicProduct::num_groups(void) {
  return gens.size();
}

int CyclicProduct::gen_order(char gen) {
  int i;
  int len = gens.size();
  char gen_to_find = tolower(gen);
  for (i=0; i<len; i++) {
    if (gens[i] == gen_to_find) {
      return orders[i];
    }
  }
  return -1;
}

int CyclicProduct::gen_index(char gen) {
  int i;
  int len = gens.size();
  char gen_to_find = tolower(gen);
  for (i=0; i<len; i++) {
    if (gens[i] == gen_to_find) {
      return i;
    }
  }
  return -1;
}


void CyclicProduct::cyc_red(std::string* S) {
  
}    
    
std::ostream &operator<<(std::ostream &os, CyclicProduct &G) {
  int i;
  int len = G.orders.size();
  for (i=0; i<len-1; i++) {
    os << "<" << G.gens[i] << ">/<" << G.orders[i] << G.gens[i] << "> * ";
  }
  os << "<" << G.gens[len-1] << ">/<" << G.orders[len-1] << G.gens[len-1] << ">";
  return os;
}


/*****************************************************************************
 * A chain   
 * ***************************************************************************/
Chain::Chain(void) {
  words.resize(0);
  weights.resize(0);
}


Chain::Chain(CyclicProduct* G_in, char** input, int num_strings) {
  int i;
  int j;
  std::string word;
  std::string weight;
  ChainChunk temp;
  int wordLen;
  int order;
  std::string first, middle, last;
  
  G = G_in;                   // note this does NOT require a copy constructor
  
  //input the raw words
  for (i=0; i<num_strings; i++) {
    j=0;
    while (isdigit(input[i][j])) {
      j++;
    }
    weight = std::string(input[i]).substr(0,j);
    word = std::string(&input[i][j]);
    if (weight == "") {
      weights.push_back(1);
    } else {
      weights.push_back( atoi( weight.c_str() ));
    }
    words.push_back(word);
  }
  
  //simplify the words into all lower case (as much as possible)
  //first rotate the words so a chunk starts the word
  for (i=0; i<(int)words.size(); i++) {
    j=0;
    while (j<(int)words[i].size()-1 && words[i][j] == words[i][j+1]) {
      j++;
    }
    j = (j+1)%words[i].size();
    words[i] = words[i].substr(j, words[i].size()-j) + words[i].substr(0, j);
  }
  
  //first pass = all lower case
  for (i=0; i<(int)words.size(); i++) {
    j=0;
    while (j<(int)words[i].size()) {
      order = (*G).gen_order(words[i][j]);
      if (order > 0 && !islower(words[i][j])) {
        first = words[i].substr(0, j);
        middle = std::string(order-1, tolower(words[i][j]));
        last = (j < (int)words[i].size()-1 ? words[i].substr(j+1, words[i].size()-j) : "");
        words[i] = first + middle + last;
        j += order-1;
      } else {
        j++;
      }
    }
  }
  
  //now reduce them
  for (i=0; i<(int)words.size(); i++) {
    j=0;
    while ( true ) {
      order = (*G).gen_order(words[i][j]);
      if (j+order >= (int)words[i].size() ) {
        break;
      }
      if (order == 0) {
        j++;
        continue;
      }
      if (words[i].substr(j, order) == std::string(order, words[i][j])) {
        words[i] = words[i].substr(0, j) + words[i].substr(j+order, words[i].size()-j-order);
      } else {
        j++;
      }
    }
  }
  
  
  //now compute the chunks -- note we may assume that a chunk starts the word
  chunks.resize(words.size());
  for (i=0; i<(int)words.size(); i++) {
    chunks[i].resize(0);
    word = words[i];
    wordLen = word.size();
    temp.word = i;
    if (words[i][0] == words[i][wordLen-1]) { //if the word only has one letter
      temp.start_index = 0;
      temp.len = wordLen;
      temp.group = (*G).gen_index(words[i][0]);
      chunks[i].push_back(temp);
      continue;
    }
    j = 0;                                //j is at the beginning of a chunk
    while (j != (int)words[i].size()) {        //stop if we've looped around
      temp.start_index = j;
      if ((*G).gen_order(word[j]) == 0) {
        temp.len = 1;
        temp.group = (*G).gen_index(word[j]);
        chunks[i].push_back(temp);
        j++;
        continue;
      }
      while (word[j] == word[(j+1)%wordLen]) {  //go until the end of the chunk
        j++;
      }
      temp.len = j-temp.start_index;
      if (temp.len < 0) {
        temp.len += wordLen;
      }
      temp.len += 1;
      temp.group = (*G).gen_index(word[j]);
      chunks[i].push_back(temp);
      j++;
    }
  }
  
  //now compute the chain letters
  ChainLetter temp_letter;
  chain_letters.resize(0);
  group_letters.resize((*G).num_groups());
  for (i=0; i<(int)words.size(); i++) {
    temp_letter.word = i;
    for (j=0; j<(int)words[i].size(); j++) {
      temp_letter.index = j;
      temp_letter.letter = words[i][j];
      temp_letter.group = (*G).gen_index(words[i][j]);
      chain_letters.push_back(temp_letter);
      group_letters[ (*G).gen_index(words[i][j]) ].push_back(chain_letters[chain_letters.size()-1]);
    }
  }
}

Chain::~Chain(void) {
  //they should get destroyed
}

std::vector<std::string> Chain::word_list(void) {
  std::vector<std::string> ans = words;
  return ans;
}

std::vector<int> Chain::weights_list(void) {
  std::vector<int> ans = weights;
  return ans;
}


std::vector<std::vector<ChainChunk> > Chain::chunk_list(void) {
  std::vector<std::vector<ChainChunk> > ans = chunks;
  return ans;
}

std::vector<std::vector<int> > Chain::group_letter_list(void) {
  std::vector<std::vector<int> > ans = group_letters;
  return ans;
}

std::vector<std::vector<ChainLetter> > Chain::chain_letter_list(void) {
  std::vector<std::vector<ChainLetter> > ans = chain_letters;
  return ans;
}

int Chain::next_letter(int n) {
  int word = chain_letters[n].word;
  int index = chain_letters[n].index;
  if (index == (int)words[word].size()-1) {
    return n-words[word].size()+1;
  } else {
    return n+1;
  }
}

int Chain::num_words(void) {
  return words.size();
}

std::string Chain::operator[](int index) {
  return words[index];
}

void Chain::print_chunks(std::ostream &os) {
  int i,j;
  for (i=0; i<(int)chunks.size(); i++) {
    for (j=0; j<(int)chunks[i].size(); j++) {
      os << "(" << chunks[i][j].word << "," 
         << chunks[i][j].start_index << "," 
         << chunks[i][j].len << "," 
         << "(" << chunks[i][j].group << ")), ";
    }
    os << "\n";
  }
}

void Chain::print_group_letters(std::ostream &os) {
  int i,j;
  for (i=0; i<(*G).num_groups(); i++) {
    for (j=0; j<group_letters[i].size(); j++) {
      os << "(" << group_letters[i][j].word << "," 
         << group_letters[i][j].index << "," 
         << group_letters[i][j].letter << "),"
    }
    os << "\n";
  }
}


std::ostream &operator<<(std::ostream &os, Chain &C) {
  int i;
  int len = C.words.size();
  for (i=0; i<len; i++) {
    os << C.weights[i] << C.words[i] << " ";
  }
  return os;
}




Multiset::Multiset() {
  L = std::vector<int>(0);
  min = 0;
  max_plus_one = 0;
}

Multiset::Multiset(int len, int Min, int Max_plus_one) {
  L = std::vector<int>(len, Min);
  min = Min;
  max_plus_one = Max_plus_one;
}

std::vector<int>* Multiset::get_list(void) {
  return &L;
}

int Multiset::next(void) {
  return 1;
}
  







