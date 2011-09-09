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

int CyclicProduct::gen_order(char gen) {
  int i;
  int len = gens.size();
  for (i=0; i<len; i++) {
    if (gens[i] == gen) {
      return orders[i];
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

Chain::Chain(char** input, int num_strings) {
  int i;
  int j;
  std::string word;
  std::string weight;
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

int Chain::num_words(void) {
  return words.size();
}

std::string Chain::operator[](int index) {
  return words[index];
}

std::ostream &operator<<(std::ostream &os, Chain &C) {
  int i;
  int len = C.words.size();
  for (i=0; i<len; i++) {
    os << C.weights[i] << C.words[i] << " ";
  }
  return os;
}













