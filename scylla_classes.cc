#include <vector>
#include <iostream>
#include <string>

#include <ctype.h>
#include <stdlib.h>

#include "scylla_classes.h"



/****************************************************************************
 * make a list of all the central edges
 * note the central edges can be anything
 * the edge is the letter just before it, and the letter just after
 ****************************************************************************/
CentralEdgeList::CentralEdgeList() {
  edges.resize(0);
  edges_beginning_with.resize(0);
}

CentralEdgeList::CentralEdgeList(Chain &C) {
  int i,j;
  int num_letters = C.num_letters();
  CentralEdge temp_central_edge;
  edges.resize(0);
  edges_beginning_with.resize(num_letters);
  for (i=0; i<num_letters; i++) {
    edges_beginning_with[i].resize(0);
    temp_central_edge.first = i;
    for (j=0; j<num_letters; j++) {
      temp_central_edge.last = j;
      edges.push_back(temp_central_edge);
      edges_beginning_with[i].push_back(edges.size()-1);
    }
  }
}

int CentralEdgeList::get_index(int a, int b) {
  int i;
  for (i=0; i<(int)edges_beginning_with[a].size(); i++) {
    if (edges[edges_beginning_with[a][i]].last == b) {
      return edges_beginning_with[a][i];
    }
  }
  return -1;
}

CentralEdge CentralEdgeList::operator[](int index) {
  return edges[index];
}


/****************************************************************************
 * make a list of all the interface edges
 * note these are (1) any gen with any inverse (2) any (non)inverse gen 
 * with another (non) inverse gen (assuming order > 0)
 * these are from the POLYGON's PERSPECTIVE
 ****************************************************************************/
InterfaceEdgeList::InterfaceEdgeList() {
  edges.resize(0);
  edges_beginning_with.resize(0);
}

InterfaceEdgeList::InterfaceEdgeList(Chain &C) {
  int i,j;
  int num_groups = (C.G)->num_groups;
  InterfaceEdge temp_interface_edge;
  
  edges.resize(0);
  edges_beginning_with.resize(C.num_letters());
  
  for (i=0; i<num_groups; i++) {
    
    //just match every letter with every possible inverse
    for (j=0; j<(int)C.regular_letters[i].size(); j++) {
      for (k=0; k<(int)C.inverse_letters[i].size(); k++) {
        temp_interface_edge.first = C.regular_letters[i][j];
        temp_interface_edge.last = C.inverse_letters[i][k];
        edges.push_back(temp_interface_edge);
        edges_beginning_with[temp_interface_edge.first].push_back(edges.size()-1);
        temp_interface_edge.first = C.inverse_letters[i][k];
        temp_interface_edge.last = C.regular_letters[i][j];
        edges.push_back(temp_interface_edge);
        edges_beginning_with[temp_interface_edge.first].push_back(edges.size()-1);
      }
    }
    
    //now match every letter with every one of the same form (if order > 0)
    if (orders[i] == 0) {
      continue;
    }
    for (j=0; j<(int)C.regular_letters[i].size(); j++) {
      temp_interface_edge.first = C.regular_letters[i][j];
      for (k=0; k<(int)C.regular_letters[i].size(); k++) {
        temp_interface_edge.last = C.regular_letters[i][k];
        edges.push_back(temp_interface_edge);
        edges_beginning_with[temp_interface_edge.first].push_back(edges.size()-1);
      }
    }
    for (j=0; j<(int)C.inverse_letters[i].size(); j++) {
      temp_interface_edge.first = C.inverse_letters[i][j];
      for (k=0; k<(int)C.inverse_letters[i].size(); k++) {
        temp_interface_edge.last = C.inverse_letters[i][k];
        edges.push_back(temp_interface_edge);
        edges_beginning_with[temp_interface_edge.first].push_back(edges.size()-1);
      }
    }
  }
}

int InterfaceEdgeList::get_index_from_poly_side(int a, int b) {
  int i;
  for (i=0; i<(int)edges_beginning_with[a].size(); i++) {
    if (edges[edges_beginning_with[a][i]].last == b) {
      return edges_beginning_with[a][i];
    }
  }
  return -1;
}

int InterfaceEdgeList::get_index_from_group_side(int a, int b) {
  int i;
  for (i=0; i<(int)edges_beginning_with[b].size(); i++) {
    if (edges[edges_beginning_with[b][i]].last == a) {
      return edges_beginning_with[b][i];
    }
  }
  return -1;
}

InterfaceEdge InterfaceEdgeList::operator[](int index) {
  return edges[index];
}


/****************************************************************************
 * make a list of all the group edges
 ****************************************************************************/
GroupEdgeList::GroupEdgeList() {
  edges.resize(0);
  edges_beginning_with.resize(0);
}

GroupEdgeList::GroupEdgeList(Chain &C, int group_index) {
  int order = (C.G)->orders[group_index];
  edges.resize(0);
  edges_beginning_with.resize(C.num_letters());
  regular_edges.resize(0);
  inverse_edges.resize(0);
  GroupEdge temp_group_edge;
  if (order == 0) {
    return;
  }
  
  for (i=0; i<(int)C.regular_edges[group_index].size(); i++) {
    temp_group_edge.first = C.regular_edges[group_index][i];
    for (j=0; j<(int)C.regular_edges[group_index].size(); j++) {
      temp_group_edge.last = C.regular_edges[group_index][j];
      edges.push_back(temp_group_edge);
      edges_beginning_with.push_back(edges.size()-1);
      regular_edges.push_back(edges.size()-1);
    }
  }
  for (i=0; i<(int)C.inverse_edges[group_index].size(); i++) {
    temp_group_edge.first = C.inverse_edges[group_index][i];
    for (j=0; j<(int)C.inverse_edges[group_index].size(); j++) {
      temp_group_edge.last = C.inverse_edges[group_index][j];
      edges.push_back(temp_group_edge);
      edges_beginning_with.push_back(edges.size()-1);
      inverse_edges.push_back(edges.size()-1);
    }
  }
}

int GroupEdgeList::get_index(int a, int b) {
  int i;
  for (i=0; i<(int)edges_beginning_with[a].size(); i++) {
    if (edges[edges_beginning_with[a][i]].last == b) {
      return edges_beginning_with[a][i];
    }
  }
  return -1;
}


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


int CyclicProduct::gen_order(char gen) {
  int i;
  int len = gens.size();
  char gen_to_find = tolower(gen);
  for (i=0; i<len; i++) {
    if (gens[i] == gen_to_find) {
      return orders[i];
    }
  }
  std::cout << "Couldn't find " << gen_to_find << "\n";
  return -1;
}

int CyclicProduct::index_order(int index) {
  return orders[index];
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

int CyclicProduct::num_groups(void) {
  return gens.size();
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
      if (j >= (int)words[i].size()) {
        break;
      }
      order = (*G).gen_order(words[i][j]);
      if (j+order >= (int)words[i].size() ) {
        j++;
        continue;
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
  group_letters.resize(G->num_groups());
  regular_letters.resize(G->num_groups());
  inverse_letters.resize(G->num_groups());
  for (i=0; i<G->num_groups(); i++) {
    group_letters[i].resize(0);
    regular_letters[i].resize(0);
    inverse_letters[i].resize(0);
  }
  for (i=0; i<(int)words.size(); i++) {
    temp_letter.word = i;
    for (j=0; j<(int)words[i].size(); j++) {
      temp_letter.index = j;
      temp_letter.letter = words[i][j];
      temp_letter.group = (*G).gen_index(words[i][j]);
      chain_letters.push_back(temp_letter);
      group_letters[ temp_letter.group ].push_back(chain_letters.size()-1);
      if (isupper(temp_letter.letter)) {
        inverse_letters[temp_letter.group].push_back(chain_letters.size()-1);
      } else {
        regular_letters[temp_letter.group].push_back(chain_letters.size()-1);
      }
    }
  }
  
  
}

Chain::~Chain(void) {
  //they should get destroyed
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

int Chain::prev_letter(int n) {
  int word = chain_letters[n].word;
  int index = chain_letters[n].index;
  if (index == 0) {
    return n+words[word].size()-1;
  } else {
    return n-1;
  }
}


int Chain::num_words(void) {
  return words.size();
}

int Chain::num_letters() {
  return chain_letters.size();
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


void Chain::print_letters(std::ostream &os) {
  int i;
  for (i=0; i<(int)chain_letters.size(); i++) {
    os << i << ": (" << chain_letters[i].word << "," 
                     << chain_letters[i].index << ","
                     << chain_letters[i].letter << ")\n";
  }
}

void Chain::print_group_letters(std::ostream &os) {
  int i,j;
  for (i=0; i<(*G).num_groups(); i++) {
    os << "Group " << i << ":\n";
    for (j=0; j<(int)group_letters[i].size(); j++) {
      os << "(" << chain_letters[group_letters[i][j]].word << "," 
         << chain_letters[group_letters[i][j]].index << "," 
         << chain_letters[group_letters[i][j]].letter << "),";
    }
    os << "\n";
  }
}


std::ostream &operator<<(std::ostream &os, Chain &C) {
  int i;
  int len = (int)C.words.size();
  for (i=0; i<len; i++) {
    os << C.weights[i] << C.words[i] << " ";
  }
  return os;
}









Multiset::Multiset() {
  L = std::vector<int>(0);
  min = 0;
  max_plus_one = 0;
  len = 0;
}

Multiset::Multiset(int Len, int Min, int Max_plus_one) {
  L = std::vector<int>(Len, Min);
  min = Min;
  max_plus_one = Max_plus_one;
  len = Len;
}

int Multiset::operator[](int index) {
  return L[index];
}

std::vector<int>* Multiset::get_list(void) {
  return &L;
}

int Multiset::next(void) {
  int i = len-1;
  int j;
  i=len-1;
  while (i>=0 && L[i] == max_plus_one-1) {
    i--;
  }
  if (i==-1) {
    return 1;
  }
  L[i]++;
  for (j=i+1; j<len; j++) {
    L[j] = min;
  }
  return 0;
}













  







