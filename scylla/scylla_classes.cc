#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>

#include <ctype.h>
#include <stdlib.h>

#include "scylla_classes.h"
#include "../word.h"

using namespace SCYLLA;

//this is useful; subtract one mod something
int sub_1_mod(int a, int m) {
  //std::cout << "Called with " << a << " and " << m << "\n";
  return (a+m-1)%m;
}


/****************************************************************************
 print a chain letter
 ***************************************************************************/
std::ostream& SCYLLA::operator<<(std::ostream &os, ChainLetter &CL) {
  os << CL.letter << " (" << CL.word << "," << CL.index << ","
     << CL.index_in_group_reg_inv_list << ")";
  return os;
}
   


/*****************************************************************************
* A free product of cyclic groups
* ****************************************************************************/
CyclicProduct::CyclicProduct(void) {
  gens.resize(0);
  orders.resize(0);
}

/* Create a cyclic group from an input of the form a0b2c3, giving Z*Z/2Z*Z/3Z
 */
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

/* Return the order of a factor, given the generator letter */
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

/* return the order of a factor, given the index of that factor */
int CyclicProduct::index_order(int index) {
  return orders[index];
}

/* return the index of a generator, given the generator letter */
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

/* Cyclically minimally reduce the word S.  This entails:
 * (1) cyclically reduce it
 * (2) replace strings of generators with inverses, if that makes it shorter
 * */
void CyclicProduct::cyc_red(std::string &S) {
  int i,j;
  char current_char;
  int current_start;
  int ord;
  std::vector<std::pair<char, int> > chunks(0);
  i=0;
  //reduce the interior of the word
  while (i < (int)S.size()-1) {
    if (S[i] == swapCaseChar(S[i+1])) {
      S.erase(i,2);
      if (i > 0) {
        i--;
      }
    } else {
      i++;
    }
  }
  
  //std::cout << "after reducing interior:" << S << "\n";
  
  //we'll cyclically reduce later, after reducing chunks
  //now make the chunks
  current_char = S[0];
  current_start = 0;
  i=1;
  while (i < (int)S.size()) {
    while (i < (int)S.size() && S[i] == current_char) {
      i++;
    }
    if (i==(int)S.size()) {
      break;
    }
    if (islower(current_char)) {
      chunks.push_back( std::make_pair( current_char, i-current_start ) );
    } else {
      chunks.push_back( std::make_pair( tolower(current_char), -(i-current_start)));
    }
    current_char = S[i];
    current_start = i;
    i++;
  }
  //finish the last chunk
  if (islower(current_char)) {
    chunks.push_back( std::make_pair( current_char, i-current_start ) );
  } else {
    chunks.push_back( std::make_pair( tolower(current_char), -(i-current_start)));
  }  

  //now reduce the chunks.  go around cyclically reducing the size 
  //and/or collecting chunks
  i=0;
  while (1) {
    //std::cout << "reducing chunk " << chunks[i].first << "," << chunks[i].second << "\n";
    //reduce the chunk
    ord = orders[ gen_index( chunks[i].first ) ];
    if (ord != 0) {
      chunks[i].second = chunks[i].second % ord;
      if (chunks[i].second < -(ord/2)) {
        chunks[i].second += ord;
      } else if (chunks[i].second > ord/2) {
        chunks[i].second -= ord;
      }
      if (chunks[i].second == 0) {
        chunks.erase(chunks.begin() + i);
        i=0;
        continue;
      }
    }
    
    if ((int)chunks.size() == 1) {
      break;
    }
    
    //combine the chunks
    j= (i+1) % chunks.size();
    if (chunks[i].first == chunks[j].first) {
      chunks[i].second += chunks[j].second;
      if (chunks[i].second == 0) {
        if (j < i) {
          chunks.erase(chunks.begin() + i);
          chunks.erase(chunks.begin() + j);
        } else {
          chunks.erase(chunks.begin() + j);
          chunks.erase(chunks.begin() + i);
        }
      } else {
        chunks.erase(chunks.begin() + j);
      }
      //start over
      i=0;
      continue;
    }
          
    i++;
    if (i==(int)chunks.size()) {
      break;
    }
  }
  
  //replace the word with the reduced one
  S = "";
  for (i=0; i<(int)chunks.size(); i++) {
    if (chunks[i].second < 0) {
      S += std::string(-chunks[i].second, toupper(chunks[i].first));
    } else {
      S += std::string(chunks[i].second, chunks[i].first);
    }
  }
}    


std::string SCYLLA::CyclicProduct::short_rep() {
  std::string rep = "";
  int len = (int)gens.size();
  if (len == 0) return rep;
  std::stringstream ss;
  if (orders[0] == 0) {
    rep += std::string(1,gens[0]);
  } else {
    ss << orders[0];
    rep += "(" + std::string(1,gens[0]) + "/" + ss.str() + std::string(1,gens[0]) + ")";
    ss.str("");
  }
  for (int i=1; i<(int)gens.size(); ++i) {
    if (orders[i] == 0) {
      rep += "*" + std::string(1,gens[i]);
    } else {
      ss << orders[i];
      rep += "*(" + std::string(1,gens[i]) + "/" + ss.str() + std::string(1,gens[i]) + ")";
      ss.str();
    }
  }
  return rep;
}
    
std::ostream& SCYLLA::operator<<(std::ostream &os, CyclicProduct &G) {
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
  
  G = G_in; //just a pointer assignment, not copy constructor
  
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
  
  //simplify the words
  for (i=0; i<(int)words.size(); i++) {
    (*G).cyc_red(words[i]);
    if ((int)words[i].size() == 0) {
      std::cout << "You gave me a trivial word\n";
      exit(1);
    }
  }
  
  //std::cout << "Done simplifying words\n";
  //std::cout.flush();
 
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

/*
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
*/

void Chain::print_letters(std::ostream &os) {
  int i;
  for (i=0; i<(int)chain_letters.size(); i++) {
    os << i << ": " << chain_letters[i] << "\n";
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


std::ostream& SCYLLA::operator<<(std::ostream &os, Chain &C) {
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
CentralEdgePairList::CentralEdgePairList() {
  edge_pairs.resize(0);
  edge_pairs_beginning_with.resize(0);
}

CentralEdgePairList::CentralEdgePairList(Chain &C) {
  int i,j;
  my_chain = &C;
  CentralEdgePair temp_central_edge_pair;
  num_letters = C.num_letters();
  edge_pairs.resize(0);
  edge_pairs_beginning_with.resize(num_letters);
  for (i=0; i<num_letters; i++) {
    edge_pairs_beginning_with[i].resize(0);
    temp_central_edge_pair.first = i;
    for (j=0; j<num_letters; j++) {
      if (C.prev_letter(j) <= i) {
        continue;
      }
      //std::cout << "pushing on " << i << " " << j << "\n";
      temp_central_edge_pair.last = j;
      edge_pairs.push_back(temp_central_edge_pair);
      edge_pairs_beginning_with[i].push_back(edge_pairs.size()-1);
    }
  }
}

//note b-1 and a+1 are computed in the WORD
//if (a,b) has minimal first entry between that and (b-1,a+1) (i.e. b-1>a),
//this returns ind+1
//otherwise, it returns -(ind+1)
int CentralEdgePairList::get_index(int a, int b) {
  int i;
  int bm1 = my_chain->prev_letter(b); // sub_1_mod(b, num_letters);
  int ap1 = my_chain->next_letter(a);
  //std::cout << "Getting index of " << a << ", " << b << "\n";
  if (bm1 > a) {
    for (i=0; i<(int)edge_pairs_beginning_with[a].size(); i++) {
      if (edge_pairs[edge_pairs_beginning_with[a][i]].last == b) {
        return edge_pairs_beginning_with[a][i]+1;
      }
    }
  } else {
    for (i=0; i<(int)edge_pairs_beginning_with[bm1].size(); i++) {
      if (edge_pairs[edge_pairs_beginning_with[bm1][i]].last == ap1) {
        return -(edge_pairs_beginning_with[bm1][i]+1);
      }
    }
  }
  return 0;  //this is bad
}

CentralEdgePair CentralEdgePairList::operator[](int index) {
  return edge_pairs[index];
}

void CentralEdgePairList::print(std::ostream &os) {
  int i;
  os << "Central Edge pairs:\n";
  for (i=0; i<(int)edge_pairs.size(); i++) {
    os << i << ": (" << edge_pairs[i].first << ", " << edge_pairs[i].last
       << "), (" << my_chain->prev_letter(edge_pairs[i].last) << ", " << my_chain->next_letter(edge_pairs[i].first) << ")\n";
  }
}


int CentralEdgePairList::size() {
  return edge_pairs.size();
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
    //don't do edge pairs of the form (i,i-1)
    if ((C.G)->orders[i] == 0) {
      continue;
    } 
    for (j=0; j<(int)C.regular_letters[i].size(); j++) {
      temp_interface_edge.first = C.regular_letters[i][j];
      for (k=0; k<(int)C.regular_letters[i].size(); k++) {
        temp_interface_edge.last = C.regular_letters[i][k];
        if (temp_interface_edge.last == C.prev_letter(temp_interface_edge.first)) {
          continue;
        }
        edges.push_back(temp_interface_edge);
        edges_beginning_with[temp_interface_edge.first].push_back(edges.size()-1);
      }
    }
    for (j=0; j<(int)C.inverse_letters[i].size(); j++) {
      temp_interface_edge.first = C.inverse_letters[i][j];
      for (k=0; k<(int)C.inverse_letters[i].size(); k++) {
        temp_interface_edge.last = C.inverse_letters[i][k];        
        if (temp_interface_edge.last == C.prev_letter(temp_interface_edge.first)) {
          continue;
        }
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

int InterfaceEdgeList::size() {
  return edges.size();
}



/*****************************************************************************
 Central polygons
 *****************************************************************************/
int CentralPolygon::chi_times_2() {
  return 2 - (int)edges.size();
}

void CentralPolygon::compute_ia_etc_for_edges(int col, 
                                              Chain &C,
                                              InterfaceEdgeList &IEL,
                                              CentralEdgePairList &CEL,
                                              SparseLP& LP) {
  int i,j,row,val;
  int num_sides = (int)edges.size();
  std::vector<int> temp_ia(0);
  std::vector<int> temp_ja(0);
  std::vector<int> temp_ar(0);
  for (j=0; j<num_sides; j++) {
    if (interface[j]) {
      row = IEL.get_index_from_poly_side(edges[j].first, edges[j].second);
      val = 1;
    } else {
      row = CEL.get_index(edges[j].first, edges[j].second);
      if (row < 0) {
        row = -row-1;
        row += IEL.size();
        val = -1;
      } else {
        row = row-1;
        row += IEL.size();
        val = 1;
      }
    }
    for (i=0; i<(int)temp_ia.size(); ++i) {
      if (temp_ia[i] == row && temp_ja[i] == col) {
        temp_ar[i] += val;
        break;
      }
    }
    if (i==(int)temp_ia.size()) {
      temp_ia.push_back(row);
      temp_ja.push_back(col);
      temp_ar.push_back(val);
    }
  }
  for (i=0; i<(int)temp_ia.size(); ++i) {
    LP.add_entry(temp_ia[i], temp_ja[i], temp_ar[i]);
  }
}  
  
  

std::ostream& SCYLLA::operator<<(std::ostream &os, CentralPolygon &CP) {
  int j;
  os << "CP: ";
  for (j=0; j<(int)CP.edges.size(); j++) {
    os << "(" << CP.edges[j].first << ", " << CP.edges[j].second << ")";
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
double GroupTooth::chi_times_2(Chain &C) {
  int ord = (C.G)->index_order( group_index );
  if (C.next_letter(first) == last) {
    return 2.0/(double)ord;
  } else {
    return 2.0/(double)ord - 1;
  }
}

void GroupTooth::compute_ia_etc_for_edges(int offset, 
                                          Chain &C,
                                          InterfaceEdgeList &IEL, 
                                          std::vector<std::vector<std::vector<int> > > &group_teeth_rows, 
                                          SparseLP& LP) {
  int row;
  int col = offset;
  int row1_offset, row2_offset;
  ChainLetter L1, L2, baseL;
  
  //do the interface edge
  if (C.next_letter(first) != last) { //only do it if it's not a dummy edge
    row = IEL.get_index_from_group_side(first, last);
    LP.add_entry(row, col, -1);
  }
  //std::cout << "Pushed interface " << row+1 << "," << col+1 << ",-1\n";

  //now do the group teeth edges
  L1 = C.chain_letters[first];
  L2 = C.chain_letters[last];
  baseL = C.chain_letters[base_letter];
  
  //std::cout << "letters: " << L1 << ", " << L2 << "\n";
  //std::cout << "based at: " << baseL << "\n";

  if (position > 0) {    
    row1_offset = group_teeth_rows[L1.group]
                                  [L1.index_in_group_reg_inv_list]
                                  [baseL.index_in_group_reg_inv_list];
    row = row1_offset + (position-1);
    LP.add_entry(row, col, -1);
    //ia.push_back(row+1);
    //ja.push_back(col+1);
    //ar.push_back(-1);
    //std::cout << "Compute row1_offset = " << row1_offset << 
    //"  so row+1 = " << row+1 << "\n";
    //std::cout << "Pushed tooth edge " << row+1 << "," << col+1 << ",-1\n";
    
  }
  if (position < (C.G)->index_order(group_index)-1) {
    row2_offset = group_teeth_rows[L2.group]
                                  [L2.index_in_group_reg_inv_list]
                                  [baseL.index_in_group_reg_inv_list];
    row = row2_offset + position;
    LP.add_entry(row, col, 1);
    //ia.push_back(row+1);
    //ja.push_back(col+1);
    //ar.push_back(1); 
    //std::cout << "Compute row2_offset = " << row2_offset << 
    //"  so row+1 = " << row+1 << "\n";
    //std::cout << "Pushed tooth edge " << row+1 << "," << col+1 << ",1\n";
  } 
}
  

void GroupTooth::compute_ia_etc_for_words(int offset, 
                                          Chain &C,
                                          int row_offset,
                                          SparseLP& LP) {
  int col = offset;
  int row;
  int word;
  if (C.chain_letters[first].index == 0) {
    word = C.chain_letters[first].word;
    row = row_offset + word;
    LP.add_entry(row, col, 1);
    //ia.push_back(row+1);
    //ja.push_back(col+1);
    //ar.push_back(1); 
  }
  //std::cout << "Pushed " << row+1 << " " << col+1 << " " << 1 << "\n";
}

std::ostream& SCYLLA::operator<<(std::ostream &os, GroupTooth &GT) {
  os << "GT: gp" << GT.group_index << " (" << GT.first << "," << GT.last << ") " << GT.position << " " << GT.base_letter;
  return os;
}



/****************************************************************************
 group rectangle
 ****************************************************************************/
void GroupRectangle::compute_ia_etc_for_edges(int col, 
                                              InterfaceEdgeList &IEL,
                                              SparseLP& LP) {
  int row;
  row = IEL.get_index_from_group_side(first, last);
  LP.add_entry(row, col, -1);
  //ia.push_back(row+1);
  //ja.push_back(col+1);
  //ar.push_back(-1);
  //std::cout << "Pushed " << row+1 << " " << col+1 << " " << -1 << "\n";

  row = IEL.get_index_from_group_side(last, first);
  LP.add_entry(row, col, -1);
  //ia.push_back(row+1);
  //ja.push_back(col+1);
  //ar.push_back(-1);
  //std::cout << "Pushed " << row+1 << " " << col+1 << " " << -1 << "\n";
}

void GroupRectangle::compute_ia_etc_for_words(int offset, 
                                              Chain &C, 
                                              int row_offset,
                                              InterfaceEdgeList &IEL,
                                              SparseLP& LP) {
  int col = offset;
  int row;
  int word;
  if (C.chain_letters[first].index == 0) { // note the other letter can't be index 0
    word = C.chain_letters[first].word;
    row = row_offset + word;
    LP.add_entry(row, col, 1);
    //ia.push_back(row+1);
    //ja.push_back(col+1);
    //ar.push_back(1);
    //std::cout << "Pushed " << row+1 << " " << col+1 << " " << 1 << "\n";
  }
  if (C.chain_letters[last].index == 0) { // note the other letter can't be index 0
    word = C.chain_letters[last].word;
    row = row_offset + word;
    LP.add_entry(row, col, 1);
    //ia.push_back(row+1);
    //ja.push_back(col+1);
    //ar.push_back(1);
    //std::cout << "Pushed " << row+1 << " " << col+1 << " " << 1 << "\n";
  }
}


std::ostream& SCYLLA::operator<<(std::ostream &os, GroupRectangle &GR) {
  os << "GR: gp" << GR.group_index << " (" << GR.first << "," << GR.last << ")";
  return os;
}









  







