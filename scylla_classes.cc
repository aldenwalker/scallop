#include <vector>
#include <iostream>
#include <string>

#include <ctype.h>
#include <stdlib.h>

#include "scylla_classes.h"






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

void CentralEdgeList::print(std::ostream &os) {
  int i;
  os << "Central Edges:\n";
  for (i=0; i<(int)edges.size(); i++) {
    os << i << ": " << edges[i].first << ", " << edges[i].last << "\n";
  }
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
  int i,j,k;
  int num_groups = (C.G)->num_groups();
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
    if ((C.G)->orders[i] == 0) {
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

void InterfaceEdgeList::print(std::ostream &os) {
  int i;
  os << "Interface Edges:\n";
  for (i=0; i<(int)edges.size(); i++) {
    os << i << ": " << edges[i].first << ", " << edges[i].last << "\n";
  }
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
  group = group_index;
  int i,j;
  if (order == 0) {
    return;
  }
  
  for (i=0; i<(int)C.regular_letters[group_index].size(); i++) {
    temp_group_edge.first = C.regular_letters[group_index][i];
    for (j=0; j<(int)C.regular_letters[group_index].size(); j++) {
      temp_group_edge.last = C.regular_letters[group_index][j];
      edges.push_back(temp_group_edge);
      edges_beginning_with[temp_group_edge.first].push_back(edges.size()-1);
      regular_edges.push_back(edges.size()-1);
    }
  }
  for (i=0; i<(int)C.inverse_letters[group_index].size(); i++) {
    temp_group_edge.first = C.inverse_letters[group_index][i];
    for (j=0; j<(int)C.inverse_letters[group_index].size(); j++) {
      temp_group_edge.last = C.inverse_letters[group_index][j];
      edges.push_back(temp_group_edge);
      edges_beginning_with[temp_group_edge.first].push_back(edges.size()-1);
      inverse_edges.push_back(edges.size()-1);
    }
  }
}

GroupEdge GroupEdgeList::operator[](int index) {
  return edges[index];
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

void GroupEdgeList::print(std::ostream &os) {
  int i;
  os << "Group " << group << " Edges:\n";
  for (i=0; i<(int)edges.size(); i++) {
    os << i << ": " << edges[i].first << ", " << edges[i].last << "\n";
  }
}



/*****************************************************************************
 Central polygons
 *****************************************************************************/
int CentralPolygon::chi_times_2(Chain &C, CentralEdgeList &CEL, InterfaceEdgeList &IEL) {
  int i;
  int sum=0;
  for (i=0; i<(int)edges.size(); i++) {
    if (interface[i]) {
      sum++;
    } else {
      if (C.next_letter( CEL[edges[i]].first ) != CEL[edges[i]].last) {
        sum++;
      }
    }
  }
  return 2 - sum;
}

void CentralPolygon::compute_ia_etc_for_edges(int col, 
                                              Chain &C,
                                              InterfaceEdgeList &IEL,
                                              CentralEdgeList &CEL,
                                              std::vector<EdgePair> &edge_pairs,
                                              std::vector<int> &central_edge_pairs,
                                              std::vector<int> &temp_ia,
                                              std::vector<int> &temp_ja,
                                              std::vector<int> &temp_ar) {
  int j,row,val;
  for (j=0; j<(int)edges.size(); j++) {
    if (interface[j]) {
      if (C.next_letter( IEL[edges[j]].last ) == IEL[edges[j]].first ) {         //don't restrict edges like this
        continue;
      }
      row = edges[j];
      val = 1;
    } else {
      if (C.next_letter( CEL[edges[j]].first ) == CEL[edges[j]].last) {   //nor these edges
        continue;
      }
      row = central_edge_pairs[edges[j]];
      if (edge_pairs[row].first == edges[j]) { 
        val = 1;
      } else {
        val = -1;
      }
    } 
    temp_ja.push_back(col+1);
    temp_ia.push_back(row+1);
    temp_ar.push_back(val);
  }
}  
  
  

std::ostream &operator<<(std::ostream &os, CentralPolygon &CP) {
  int j;
  os << "CP: ";
  for (j=0; j<(int)CP.edges.size(); j++) {
    os << CP.edges[j];
    if (CP.interface[j]) {
      os << "i";
    } else {
      os << "c";
    }
    os << " ";
  }
  return os;
}


/****************************************************************************
 group tooth
 ****************************************************************************/
int GroupTooth::chi_times_2(Chain &C) {
  int temp_letter_1, temp_letter_2;
  if (inverse) {
    temp_letter_1 = C.inverse_letters[group_index][first];
    temp_letter_2 = C.inverse_letters[group_index][last];
  } else {
    temp_letter_1 = C.regular_letters[group_index][first];
    temp_letter_2 = C.regular_letters[group_index][last];
  }
  if (C.next_letter(temp_letter_1) == temp_letter_2) {
    return 0;
  } else {
    return -1;
  }
}

void GroupTooth::compute_ia_etc_for_edges(int offset, 
                                          Chain &C,
                                          InterfaceEdgeList &IEL, 
                                          std::vector<std::vector<std::vector<int> > > &rows_for_letters_in_mouths, 
                                          std::vector<int> &ia, 
                                          std::vector<int> &ja, 
                                          std::vector<double> &ar) {
  int row;
  int col = offset;
  
  int temp_letter_1, temp_letter_2;
  if (inverse) {
    temp_letter_1 = C.inverse_letters[group_index][first];
    temp_letter_2 = C.inverse_letters[group_index][last];
  } else {
    temp_letter_1 = C.regular_letters[group_index][first];
    temp_letter_2 = C.regular_letters[group_index][last];
  }
  
  if (C.next_letter(temp_letter_1) != temp_letter_2) { 
    row = IEL.get_index_from_group_side(temp_letter_1, temp_letter_2);
    ia.push_back(row+1);
    ja.push_back(col+1);
    ar.push_back(-1);
    //std::cout << "Pushed " << row+1 << " " << col+1 << " " << -1 << "\n";
  }
  
  row = rows_for_letters_in_mouths[group_index][group_mouth_index][first] + position;
  ia.push_back(row+1);
  ja.push_back(col+1);
  ar.push_back(1);
  //std::cout << "Pushed " << row+1 << " " << col+1 << " " << 1 << "\n";
  
  row = rows_for_letters_in_mouths[group_index][group_mouth_index][last] + position + 1;
  ia.push_back(row+1);
  ja.push_back(col+1);
  ar.push_back(-1);
  //std::cout << "Pushed " << row+1 << " " << col+1 << " " << -1 << "\n";
}
  

void GroupTooth::compute_ia_etc_for_words(int offset, 
                                          int row_offset,
                                          Chain &C,
                                          std::vector<int> &ia, 
                                          std::vector<int> &ja, 
                                          std::vector<double> &ar) {
  int col = offset;
  int row;
  if (inverse) {
    row = row_offset + C.chain_letters[C.inverse_letters[group_index][first]].word;
  } else {
    row = row_offset + C.chain_letters[C.regular_letters[group_index][first]].word;
  }
  ia.push_back(row+1);
  ja.push_back(col+1);
  ar.push_back(1);
  //std::cout << "Pushed " << row+1 << " " << col+1 << " " << 1 << "\n";
}

std::ostream &operator<<(std::ostream &os, GroupTooth &GT) {
  os << "GT(position " << GT.position << ", mouth " << GT.group_mouth_index << "): " << GT.first << " " << GT.last;
  return os;
}

/****************************************************************************
 group mouth
 ****************************************************************************/
int GroupMouth::chi_times_2(Chain &C) {
  if (last == first) {
    return 2;
  } else {
    return 1;
  }
}

void GroupMouth::compute_ia_etc_for_edges(int offset,
                                          Chain &C,
                                          GroupEdgeList &GEL,
                                          int my_index,
                                          std::vector<std::vector<std::vector<int> > > &rows_for_letters_in_mouths, 
                                          std::vector<EdgePair> &edge_pairs,
                                          std::vector<int> &group_edge_pairs,
                                          std::vector<int> &ia,
                                          std::vector<int> &ja,
                                          std::vector<double> &ar) {
  int col = offset;
  int row;
  int edge_num;
  int val;
  row = rows_for_letters_in_mouths[group_index][my_index][first_letter_index] + 0;
  ia.push_back(row+1);
  ja.push_back(col+1);
  ar.push_back(-1);  
  //std::cout << "GM " << my_index << "\n";
  //std::cout << "Pushed " << row+1 << " " << col+1 << " " << -1 << "\n";
  
  row = rows_for_letters_in_mouths[group_index][my_index][last_letter_index] + (C.G)->orders[group_index];
  ia.push_back(row+1);
  ja.push_back(col+1);
  ar.push_back(1);  
  //std::cout << "Pushed " << row+1 << " " << col+1 << " " << 1 << "\n";
  
  if (last != first) {
    edge_num = GEL.get_index(last, first);
    row = group_edge_pairs[edge_num];
    if (edge_pairs[row].first == edge_num) {
      val = 1;
    } else {
      val = -1;
    }
    ia.push_back(row+1);
    ja.push_back(col+1);
    ar.push_back(val);
    //std::cout << "Pushed " << row+1 << " " << col+1 << " " << val << "\n";
  }
}

std::ostream &operator<<(std::ostream &os, GroupMouth &GM) {
  os << "GM: " << GM.first << "(" << GM.first_letter_index << ") " 
               << GM.last << "(" << GM.last_letter_index << ")";
  return os;
}



/****************************************************************************
 group polygon/square
 ****************************************************************************/
int GroupPolygon::chi_times_2(GroupEdgeList &GEL) {
  int i;
  int temp_letter_1, temp_letter_2;
  int sum=0;
  //std::cout << "computing chi for group polygon " << *this << " in group " << group << "\n";
  for (i=0; i<(int)edges.size(); i++) {
    //compute all the interface edges in the side
    temp_letter_1 = GEL[edges[i]].first;
    temp_letter_2 = GEL[edges[i]].last;
    if (temp_letter_1 != temp_letter_2) {
      sum++;
    }
  }
  //std::cout << "returning " << 2-sum << "\n";
  return 2 - sum;
}



void GroupPolygon::compute_ia_etc_for_edges(int offset, 
                                             Chain &C, 
                                             InterfaceEdgeList &IEL, 
                                             GroupEdgeList &GEL, 
                                             std::vector<EdgePair> &edge_pairs,
                                             std::vector<int> &group_edge_pairs, 
                                             std::vector<int> &temp_ia, 
                                             std::vector<int> &temp_ja, 
                                             std::vector<int> &temp_ar) {
  int i;
  int temp_letter_1, temp_letter_2;
  int temp_index;
  int row;
  for (i=0; i<(int)edges.size(); i++) {
    temp_letter_1 = GEL[edges[i]].first;
    temp_letter_2 = GEL[edges[i]].last;
    if (temp_letter_1 != temp_letter_2) {
      temp_index = edges[i];
      row = group_edge_pairs[temp_index];
      temp_ia.push_back(row+1);
      temp_ja.push_back(offset+1);
      if (edge_pairs[row].first == temp_index) {
        temp_ar.push_back(1);
      } else {
        temp_ar.push_back(-1);
      }
    }
  }
}

/*
void GroupPolygon::get_ia_etc_for_words(Chain &C, 
                                        InterfaceEdgeList &IEL, 
                                        GroupEdgeList &GEL, 
                                        std::vector<EdgePair> &edge_pairs,
                                        std::vector<int> &group_edge_pairs, 
                                        int &offset, 
                                        std::vector<int> &temp_ia, 
                                        std::vector<int> &temp_ja, 
                                        std::vector<int> &temp_ar) {
  int i,j;
  int temp_letter_1;
  int temp_index;
  for (i=0; i<(int)sides.size(); i++) {
    //get the letter from the group edge
    temp_letter_1 = GEL[edges[i]].first;
    temp_index = edge_pairs.size() + C.chain_letters[temp_letter_1].word;
    temp_ia.push_back(temp_index+1);  //temp_index is the actual row
    temp_ja.push_back(offset+1);
    temp_ar.push_back(1);
      
    for (j=0; j<(int)sides[i].letters.size(); j++) {
      temp_letter_1 = sides[i].letters[j];
      temp_index = edge_pairs.size() + C.chain_letters[temp_letter_1].word;
      temp_ia.push_back(temp_index+1);  //temp_index is the actual row
      temp_ja.push_back(offset+1);
      temp_ar.push_back(1);
    }
  }
}
 */

std::ostream &operator<<(std::ostream &os, GroupPolygon &GP) {
  int k;
  os << "GP(" << GP.group << "): ";
  for (k=0; k<(int)GP.edges.size(); k++) {
    os << GP.edges[k] << " ";
  }
  return os;
}



/****************************************************************************
 group rectangle
 ****************************************************************************/
void GroupRectangle::compute_ia_etc_for_edges(int offset, 
                                              std::vector<int> &ia, 
                                              std::vector<int> &ja, 
                                              std::vector<double> &ar) {
  int row,col;
  col = offset;
  row = first;
  ia.push_back(row+1);
  ja.push_back(col+1);
  ar.push_back(-1);
  //std::cout << "Pushed " << row+1 << " " << col+1 << " " << -1 << "\n";

  row = last;
  ia.push_back(row+1);
  ja.push_back(col+1);
  ar.push_back(-1);
  //std::cout << "Pushed " << row+1 << " " << col+1 << " " << -1 << "\n";
}

void GroupRectangle::compute_ia_etc_for_words(int offset, 
                                              int row_offset,
                                              Chain &C, 
                                              InterfaceEdgeList &IEL,
                                              std::vector<int> &ia,
                                              std::vector<int> &ja,
                                              std::vector<double> &ar) {
  int col = offset;
  int word_1 = C.chain_letters[IEL.edges[first].first].word;
  int word_2 = C.chain_letters[IEL.edges[last].first].word;
  if (word_1 == word_2) {
    ia.push_back(word_1 + row_offset + 1);
    ja.push_back(col + 1);
    ar.push_back(2);
    //std::cout << "Put " << word_1 + row_offset + 1 << ", " << col+1 << ", " << 2 << ".\n";
  } else {
    ia.push_back(word_1 + row_offset + 1);
    ja.push_back(col + 1);
    ar.push_back(1);
    //std::cout << "Put " << word_1 + row_offset + 1 << ", " << col+1 << ", " << 1 << ".\n";
    ia.push_back(word_2 + row_offset + 1);
    ja.push_back(col + 1);
    ar.push_back(1);
    //std::cout << "Put " << word_2 + row_offset + 1 << ", " << col+1 << ", " << 1 << ".\n";
  }
}


std::ostream &operator<<(std::ostream &os, GroupRectangle &GR) {
  os << "GR: " << GR.first << " " << GR.last;
  return os;
}



/****************************************************************************
 * a multiset
 * **************************************************************************/


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
  //std::cout << "next called on:";
  //for (j=0; j<L.size(); j++) {
  //  std::cout << L[j] << " ";
  //}
  while (i>=0 && L[i] == max_plus_one-1) {
    i--;
  }
  if (i==-1) {
    //std::cout << "returning 1\n";
    return 1;
  }
  L[i]++;
  for (j=i+1; j<len; j++) {
    L[j] = min;
  }
  //std::cout << "returning 0\n";
  return 0;
}


void Multiset::set_index(int index, int val) {
  L[index] = val;
}










  







