#include <vector>
#include <iostream>
#include <string>
#include <utility>

#include <ctype.h>
#include <stdlib.h>

#include "trollop_classes.h"
#include "../word.h"

using namespace TROLLOP;

//this is useful; subtract one mod something
int sub_1_mod(int a, int m) {
  //std::cout << "Called with " << a << " and " << m << "\n";
  return (a+m-1)%m;
}


//raise an integer to an integer power -- this isn't in the 
//standard library?
int int_pow(int a, int b) {
  int e = b;
  int current_base = a;
  int collect = 1;
  while (e != 0) {
    if ((e&1) == 1) {
      collect *= current_base;
    }
    current_base *= current_base;
    e = e >> 1;
  }
  return collect;
}



//convert a character to a number (for base 2*rank or 2*rank-1)
int letter_to_number(char let, int rank, char letter_removed) {
  int letter_val;
  int removed_val;
  //std::cout << "Finding number for " << let << " with " << letter_removed 
  //          << " removed\n";
  if ((int)let < 97) {
    letter_val = (int)let - 65;
  } else {
    letter_val = rank + (int)let - 97;
  }
  if (letter_removed == '\0') {
    //working base 2*rank
    //std::cout << "got " << letter_val << "\n";
    return letter_val;
  } else {
    if ((int)letter_removed < 97) {
      removed_val = (int)letter_removed - 65;
    } else {
      removed_val = rank + (int)letter_removed - 97;
    }
    //std::cout << "Removed is " << removed_val << " so returning ";
    if (removed_val < letter_val) {
      //std::cout << letter_val-1 << "\n";
      return letter_val - 1;
    } else {
      //std::cout << letter_val << "\n";
      return letter_val;
    }
  }
}

//convert a number to a character 
char number_to_letter(int n, int rank, char letter_removed) {
  int looking_for;
  int val_removed;
  if (letter_removed != '\0') {
    if ((int)letter_removed < 97) {
      val_removed = (int)letter_removed - 65;
    } else {
      val_removed = rank + (int)letter_removed - 97;
    }
    if (n >= val_removed) {
      looking_for = n + 1;
    } else {
      looking_for = n;
    }
  } else {
    looking_for = n;
  }
  
  if (looking_for < rank) {
    return (char)(65 + looking_for);
  } else {
    return (char)(97 + (looking_for-rank));
  }
}



//convert a word to an integer, where the letters represent digits
//note which digits the letters represent changes with every letter!
//the leftmost letter is the most significant
//note we can't include the starting letter
int word_to_int(std::string& S, int start_index, int rank) {
  int val = 0;
  int i;
  int len = S.size();
  int base = 2*rank-1;
  int current_multiplier = 1;
  for (i=len-1; i>=start_index; i--) {
    val += current_multiplier * letter_to_number(S[i], rank, inverse_char(S[i-1]));
    current_multiplier *= base;
  }
  return val;
}

//here S is assumed to have the correct size
//also, a word of length 1 is done
void int_to_word(std::string &S, int n, int rank) {
  int Q, R;
  int i;
  int len = S.size();
  if (len == 1) {
    return;
  }
  int base = 2*rank-1;
  int current_multiplier = int_pow(base, len-2);
  for (i=1; i<len; i++) {
    Q = n/current_multiplier;
    R = n%current_multiplier;
    S[i] = number_to_letter(Q, rank, inverse_char(S[i-1]));
    n = R;
    current_multiplier /= base;
  }
}



/****************************************************************************
 WordTable methods
 *****************************************************************************/
WordTable::WordTable(int R, int L, bool DO_SUP) {
  rank = R;
  ell = L;
  do_sup = DO_SUP;
  have_assigned_indices = false;
  first_letter_offset_edges = int_pow(2*rank-1, ell-1);
  first_letter_offset_vertices = int_pow(2*rank-1, ell-2);
  num_real_edges = 2*rank * int_pow(2*rank-1, ell-1);
  num_real_verts = 2*rank * int_pow(2*rank-1, ell-2);
  if (!do_sup) {
    num_edges = -1;
    num_verts = -1;
  } else {
    num_verts = num_real_verts;
    num_edges = num_real_edges;
  }
  vertex_to_index = std::vector<int>(0);
  edge_to_index = std::vector<int>(0);
  index_to_vertex = std::vector<int>(0);
  index_to_edge = std::vector<int>(0);
  h_vector = std::vector<std::vector<int> >(0);
}


void WordTable::create_index_assignments(WVec &C) {
  int current_val;
  int i;
  int h_val;
  
  h_vector = std::vector<std::vector<int> > (2*rank, std::vector<int>(0));
  
  if (do_sup) {
    //do no assignments
    for (i=0; i<num_real_edges; i++) {
      h_val = letter_to_number(get_h_from_real_index(i), rank, '\0');
      h_vector[h_val].push_back(i);
    }
    num_edges = num_real_edges;
    num_verts = num_real_verts;
  } else {
    //go through the words which appear in C, and 
    //give them indices, and that's it
    edge_to_index = std::vector<int> (num_real_edges, 0);
    vertex_to_index = std::vector<int> (num_real_verts, 0);
    for (i=0; i<(int)C.word_list.size(); i++) {
      edge_to_index[ C.word_list[i].second ] = -1;
      vertex_to_index[ get_real_edge_dest(C.word_list[i].second) ] = -1;
      vertex_to_index[ get_real_edge_source(C.word_list[i].second) ] = -1;
    }
    current_val = 0;
    index_to_edge = std::vector<int>(0);
    for (i=0; i<num_real_edges; i++) {
      if (edge_to_index[i] == -1) {
        edge_to_index[i] = current_val;
        index_to_edge.push_back(i);
        h_val = letter_to_number(get_h_from_real_index(i), rank, '\0');
        h_vector[h_val].push_back(current_val);
        current_val++;
      }
    }
    num_edges = current_val;
    
    current_val = 0;
    index_to_vertex = std::vector<int>(0);
    for (i=0; i<num_real_verts; i++) {
      if (vertex_to_index[i] == -1) {
        vertex_to_index[i] = current_val;
        index_to_vertex.push_back(i);
        current_val++;
      }
    }
    num_verts = current_val;
    
  }
  have_assigned_indices = true;
}

void WordTable::print() {
  int i;
  std::string S;
  std::cout << "Word table:\nRank: " << rank << "\nell: " << ell << "\n";
  std::cout << "Number of edges: " << num_edges << "(full: " << num_real_edges << ")\n";
  std::cout << "Number of vertices: " << num_verts << "(full: " << num_real_verts << ")\n";
  std::cout << "Edges:\n";
  if (do_sup) {
    std::cout << "(there is no edge index assignment)\n";
    for (i=0; i<num_edges; i++) {
      get_real_word(S, ell, i);
      std::cout << i << " (" << get_index(S) << ") : " 
                << S << "; (" << get_edge_source(i) << "," << get_edge_dest(i) 
                << "); h = " << get_h_from_real_index(i) << "\n"; 
    }
  } else {
    for (i=0; i<num_edges; i++) {
      get_word(S, ell, i);
      std::cout << i << " (" << get_index(S) << ") : " 
                << S << "; real index: " << index_to_edge[i] 
                << "; (" << get_edge_source(i) << "," << get_edge_dest(i) 
                << "); h = " << get_h(i) << "\n"; 
    }    
  }
  std::cout << "Vertices: \n";
  if (do_sup) {
    std::cout << "(there is no vertex index assignment)\n";
    for (i=0; i<num_verts; i++) {
      get_real_word(S, ell-1, i);
      std::cout << i << ":  (" << get_index(S) << ") : "  << S 
      << " outgoing h = " << get_outgoing_h(i) << "\n"; 
    }
  } else {
    for (i=0; i<num_verts; i++) {
      get_word(S, ell-1, i);
      std::cout << i << ": (" << get_index(S) << ") : " << S << "; real index: " << index_to_vertex[i] 
      << " outgoing h = " << get_outgoing_h(i) << "\n";     
    }    
  }  
}

//this takes the first letter, computes the offset, 
//and then takes the rest of the word to be in base 
//(2*rank - 1), with the letters listed ABC...abc...
int WordTable::get_real_index(std::string &S) {
  int offset = -1;
  int value = -1;
  if ((int)S.size() == ell) { //it's an edge
    offset = first_letter_offset_edges * letter_to_number(S[0], rank, '\0');
  } else if ((int)S.size() == ell-1) { //it's a vertex
    offset = first_letter_offset_vertices * letter_to_number(S[0], rank, '\0');
  } else {
    std::cout << "Word " << S << " not the right size\n";
    return -1;
  }
  value = word_to_int(S, 1, rank);
  return offset + value;
}
  

void WordTable::get_real_word(std::string &S, int length, int index) {  
  int offset = -1;
  int ind = index;
  if (length == ell) {
    offset = ind / first_letter_offset_edges;
    ind -= (offset*first_letter_offset_edges);
  } else {
    offset = ind / first_letter_offset_vertices;
    ind -= (offset*first_letter_offset_vertices);
  }
  S.resize(length);
  S[0] = number_to_letter(offset, rank, '\0');
  int_to_word(S, ind, rank);
}


//only makes sense to call for edges
char WordTable::get_h_from_real_index(int index){
  int offset = index / first_letter_offset_edges;
  return number_to_letter(offset, rank, '\0');
}



int WordTable::get_real_edge_dest(int edge) {
  std::string S;
  get_real_word(S, ell, edge);
  S = S.substr(1, S.size()-1);
  return get_real_index(S);
}

int WordTable::get_real_edge_source(int edge){
  std::string S;
  get_real_word(S, ell, edge);
  S = S.substr(0, S.size()-1);
  return get_real_index(S);
}



char WordTable::get_real_outgoing_h(int vertex) {
  int offset = vertex / first_letter_offset_vertices;
  return number_to_letter(offset, rank, '\0');
}  



//given a word return the assigned index
int WordTable::get_index(std::string &S) {
  if (do_sup) {
    return get_real_index(S);
  } else {
    if ((int)S.size() == ell) { 
      return edge_to_index[get_real_index(S)];
    } else {
      return vertex_to_index[get_real_index(S)];
    }
  }
}

//given an index, return the assigned word
void WordTable::get_word(std::string &S, int length, int index) {
  if (do_sup) {
    get_real_word(S, length, index);
  } else {
    if (length == ell) {
      get_real_word(S, length, index_to_edge[index]);
    } else {
      get_real_word(S, length, index_to_vertex[index]);
    }
  }
}


char WordTable::get_h(int index) {
  return get_h_from_real_index(index_to_edge[index]);
}

int WordTable::get_edge_dest(int index) {
  if (do_sup) {
    return get_real_edge_dest(index);
  } else {
    return vertex_to_index[get_real_edge_dest(index_to_edge[index])];
  }
}

int WordTable::get_edge_source(int index) {
  if (do_sup) {
    return get_real_edge_source(index);
  } else {
    return vertex_to_index[get_real_edge_source(index_to_edge[index])];
  }
}
  

char WordTable::get_outgoing_h(int index) {
  if (do_sup) {
    return get_real_outgoing_h(index);
  } else {
    return get_real_outgoing_h(index_to_vertex[index]);
  }
}



/**************************************************************************
 WVec methods
 **************************************************************************/
void get_cyclic_subword(std::string &S, std::string &T, int start, int n) {
  int i;
  int TLen = T.size();
  S.resize(n);
  for (i=0; i<n; i++) {
    S[i] = T[(start+i)%TLen];
  }
}

WVec::WVec() {
  //do nothing
}

WVec::WVec (WordTable &WT, char** words, int num_words, bool USE_WORDS) {
  int i,j;
  std::string word;
  std::string coef_str;
  int coef;
  index_coefficients.resize(0);
  original_input.resize(0);
  std::string input_word;
  for (i=0; i<num_words; i++) {
    original_input.push_back(std::string(words[i]));
  }
  word_list.resize(0);
  wt = &WT;
  
  for (i=0; i<num_words; i++) {
    //std::cout << "Loading input word " << original_input[i] << "\n";
    //std::cout.flush();
    j=0;
    while (isdigit(original_input[i][j])) {
      j++;
    }
    if (j == 0) {
      coef = 1;
    } else {
      coef_str = original_input[i].substr(0, j);
      coef = atoi(coef_str.c_str());
    }
    input_word = original_input[i].substr(j,original_input[i].size()-j);
    
    //cyclically reduce!
    red(input_word);
    cyc_red(input_word);
    
    //std::cout << "Got input word: " << input_word << "\n";
    //std::cout.flush();
    
    //handle the word
    for (j=0; j<(int)input_word.size(); j++) {
      get_cyclic_subword(word, input_word, j, WT.ell);
      word_list.push_back(std::make_pair(coef, WT.get_real_index(word)));
      //std::cout << "Just input word: " << word <<"\n"; 
      //std::cout.flush();
    }
  }
  
  //collect the word list
  for (i=0; i<(int)word_list.size()-1; i++) {
    for (j=i+1; j<(int)word_list.size(); j++) {
      if (word_list[i].second == word_list[j].second) {
        word_list[i].first += word_list[j].first;
        word_list.erase(word_list.begin() + j);
        break;
      }
    }
  }
  i=0;
  while (i < (int)word_list.size()) {
    if (word_list[i].first == 0) {
      word_list.erase(word_list.begin() + i);
    } else {
      i++;
    }
  }
  
}


void WVec::fill_index_coefficients(WordTable &WT) {
  int i;
  index_coefficients.resize(WT.num_edges);
  for (i=0; i<WT.num_edges; i++) {
    index_coefficients[i] = 0;
  }
  for (i=0; i<(int)word_list.size(); i++) {
    if (WT.do_sup) {
      index_coefficients[ word_list[i].second ] += word_list[i].first;
    } else {
      index_coefficients[ WT.edge_to_index[word_list[i].second] ] += word_list[i].first;
    }
  }
  
}

void WVec::print_words(std::ostream& os) {
  int i;
  std::string temp;
  for (i=0; i<(int)word_list.size(); i++) {
    (*wt).get_real_word(temp, (*wt).ell, word_list[i].second);
    if (word_list[i].first == 1) {
      os << temp << " ";
    } else {
      os << word_list[i].first << temp << " ";
    }
  }
}

std::ostream& operator<<(std::ostream& os, WVec& C) {
  int i;
  for (i=0; i<(int)C.original_input.size(); i++) {
    std::cout  << C.original_input[i] << " ";
  }
  for (i=0; i<(int)C.word_list.size(); i++) {
    std::cout << "<" << C.word_list[i].first << ", " << C.word_list[i].second << "> ";
  }
  return os;
}


/***************************************************************************
 ArcPairList methods
 **************************************************************************/
 //simply go through all vertices and make all pairs
ArcPairList::ArcPairList(WordTable &WT, int num_copies) {
  int i,j,sidea,sideb;
  num_arcs = -500000;
  arc_list.resize(0);
  arcs_beginning_with.resize(num_copies);
  for (i=0; i<num_copies; i++) {
    arcs_beginning_with[i].resize(WT.num_verts);
  }
  for (sidea=0; sidea<num_copies; sidea++) {
    for (i=0; i<WT.num_verts; i++) {
      arcs_beginning_with[sidea][i].resize(0);
      for (sideb=sidea; sideb<num_copies; sideb++) {
        for (j=(sidea<sideb? 0 : i+1); j<WT.num_verts; j++) {
          arc_list.push_back( std::make_pair( std::make_pair(sidea, i), 
                                              std::make_pair(sideb, j) ) );
          arcs_beginning_with[sidea][i].push_back(arc_list.size()-1);
        }
      }
    }
  }
  num_arcs = (int)arc_list.size();
}

//if a<b, this returns 1+index
//if b<a, this returns -(1+index)
int ArcPairList::index_of_arc(int sidea, int a, int sideb, int b) {
  int i;
  if (sidea < sideb || (sidea==sideb && a<b)) {
    for (i=0; i<(int)arcs_beginning_with[sidea][a].size(); i++) {
      if (arc_list[arcs_beginning_with[sidea][a][i]].second == std::make_pair(sideb, b) ) {
        return 1 + arcs_beginning_with[sidea][a][i];
      }
    }
  } else {
    for (i=0; i<(int)arcs_beginning_with[sideb][b].size(); i++) {
      if (arc_list[arcs_beginning_with[sideb][b][i]].second == std::make_pair(sidea, a) ) {
        return -(1 + arcs_beginning_with[sideb][b][i]);
      }
    }    
  }
  return -1;
}

void ArcPairList::print(std::ostream& os) {
  int i;
  os << "Arc pair list: (" << num_arcs << " pairs)\n";
  for (i=0; i<(int)arc_list.size(); i++) {
    os << i << ": (" << arc_list[i].first.first << "-" << arc_list[i].first.second << ", " 
                     << arc_list[i].second.first << "-" << arc_list[i].second.second << ")\n";
  }
}



std::ostream& operator<<(std::ostream& os, Triangle& T) {
  os << "T verts: (" << T.v0 << "," << T.v1 << "," << T.v2 
     << ") arcs: (";
  if (T.a0 < 0) { os << "-" << -T.a0-1; } else { os<< T.a0-1;} os << ",";
  if (T.a1 < 0) { os << "-" << -T.a1-1; } else { os<< T.a1-1;} os << ",";
  if (T.a2 < 0) { os << "-" << -T.a2-1; } else { os<< T.a2-1;} os << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, Rectangle& R) {
  int sign, index;
  extract_signed_index(&sign, &index, R.a0);
  if (sign < 0) {
    os << "R(-" << index << ", ";
  } else {
    os << "R(" << index << ", ";
  }
  os << R.side0 << "-" << R.e0 << ", ";
  extract_signed_index(&sign, &index, R.a1);
  if (sign < 0) {
    os << "-" << index << ", ";
  } else {
    os << index << ", ";
  }
  os << R.side1 << "-" << R.e1 << ")";
  return os;
}




void extract_signed_index(int* sign, int* index, int signed_index) {
  if (signed_index < 0) {
    *sign = -1;
    *index = -signed_index-1;
  } else {
    *sign = 1;
    *index = signed_index-1;
  }
}




































