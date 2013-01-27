#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <utility>
#include <cmath>

#include "../rational.h"
#include "../word.h"
#include "../lp.h"

#include "scabble.h"



/****************************************************************************
 a n-dimensional point class
 ****************************************************************************/

SCABBLE::Pt::Pt() {
  x.resize(0);
}

SCABBLE::Pt::Pt(int dim) {
  x.resize(dim);
}

SCABBLE::Pt::Pt(const SCABBLE::Pt& p) {
  int dim = (int)p.x.size();
  x.resize(dim);
  for (int i=0; i<dim; ++i) {
    x[i] = p.x[i];
  }
}

SCABBLE::Pt::Pt(int dim, Rational& init) {
  x.resize(dim);
  for (int i=0; i<dim; ++i) {
    x[i] = init;
  }
}

SCABBLE::Pt::Pt(int dim, int init) {
  x.resize(dim);
  for (int i=0; i<dim; ++i) {
    x[i] = init;
  }
}

int SCABBLE::Pt::dim() {
  return (int)x.size();
}

Rational& SCABBLE::Pt::operator[](int i) {
  return x[i];
}

SCABBLE::Pt SCABBLE::Pt::operator-(const SCABBLE::Pt& other) {
  int dim = (int)x.size();
  SCABBLE::Pt p(dim);
  for (int i=0; i<dim; ++i) {
    p[i] = x[i] - other.x[i];
  }
  return p;
}

SCABBLE::Pt SCABBLE::Pt::operator+(const SCABBLE::Pt& other) {
  int dim = (int)x.size();
  SCABBLE::Pt p(dim);
  for (int i=0; i<dim; ++i) {
    p[i] = x[i] + other.x[i];
  }
  return p;
}

SCABBLE::Pt SCABBLE::Pt::operator-() {
  int dim = (int)x.size();
  SCABBLE::Pt p(dim);
  for (int i=0; i<dim; ++i) {
    p[i] = -x[i];
  }
  return p;
}

Rational SCABBLE::Pt::dot(const SCABBLE::Pt& other) {
  Rational r = 0;
  int dim = (int)x.size();
  for (int i=0; i<dim; ++i) {
    r += (x[i] * other.x[i]);
  }
  return r;
}

SCABBLE::Pt SCABBLE::Pt::cross(const SCABBLE::Pt& other) {
  if ((int)x.size() != 3) {
    std::cout << "Can't cross non-3-dimensional vectors!\n";
    return SCABBLE::Pt();
  }
  SCABBLE::Pt p(3);
  p[0] = x[1]*other.x[2] - x[2]*other.x[1];
  p[1] = x[2]*other.x[0] - x[0]*other.x[2];
  p[2] = x[0]*other.x[1] - x[1]*other.x[0];
  return p;
}

void SCABBLE::Pt::rescale_to_integer() {
  int dim = (int)x.size();
  int L = 1;
  for (int i=0; i<dim; ++i) {
    L = lcm(L, x[i].d());
  }
  for (int i=0; i<dim; ++i) {
    x[i] *= L;
  }
}

SCABBLE::Pt& SCABBLE::Pt::operator*=(const Rational& m) {
  int dim = (int)x.size();
  for (int i=0; i<dim; ++i) {
    x[i] *= m;
  }
  return *this;
}
  
SCABBLE::Pt& SCABBLE::Pt::operator*=(int m) {
  int dim = (int)x.size();
  Rational r(m);
  for (int i=0; i<dim; ++i) {
    x[i] *= r;
  }
  return *this;
}

SCABBLE::Pt& SCABBLE::Pt::operator/=(const Rational& m) {
  int dim = (int)x.size();
  for (int i=0; i<dim; ++i) {
    x[i] /= m;
  }
  return *this;
}

SCABBLE::Pt& SCABBLE::Pt::operator/=(int m) {
  int dim = (int)x.size();
  Rational r(m);
  for (int i=0; i<dim; ++i) {
    x[i] /= r;
  }
  return *this;
}

//negate the coordinates which are masked in i
SCABBLE::Pt SCABBLE::Pt::negate_coords(int i) {
  int dim = (int)x.size();
  SCABBLE::Pt p(dim);
  for (int j=0; j<dim; ++j) {
    if ( ((i>>j)&1) == 1 ) {
      p[j] = -x[j];
    } else {
      p[j] = x[j];
    }
  }
  return p;
}

std::ostream& SCABBLE::operator<<(std::ostream& os, const SCABBLE::Pt& p) {
  int dim = p.x.size();
  os << "[";
  if (dim > 0) {
    os << p.x[0];
  }
  for (int i=1; i<dim; ++i) {
    os << ", " << p.x[i];
  }
  os << "]";
  return os;
}



    


//compute a normal and normal value 
void SCABBLE::affine_hyperplane(std::vector<SCABBLE::Pt>& face, 
                                SCABBLE::Pt& normal, 
                                Rational& normal_value) {
  if ((int)face.size() == 2) {
    SCABBLE::Pt parallel = face[1]-face[0];
    normal[0] = -parallel[1];
    normal[1] = parallel[0];
    normal.rescale_to_integer();
    normal_value = normal.dot(face[0]);
    normal *= normal_value.d();
    normal_value = normal.dot(face[0]);
  } else if ((int)face.size() == 3) {
    SCABBLE::Pt parallel1 = face[1]-face[0];
    SCABBLE::Pt parallel2 = face[2]-face[0];
    normal = parallel1.cross(parallel2);
    normal.rescale_to_integer();
    normal_value = normal.dot(face[0]);
    normal *= normal_value.d();
    normal_value = normal.dot(face[0]);
  } else {
   std::cout << "Not implemented\n"; 
  } 
}




/****************************************************************************
 print a chain letter
 ***************************************************************************/
std::ostream& SCABBLE::operator<<(std::ostream &os, SCABBLE::ChainLetter &CL) {
  os << CL.letter << " (" << CL.word << "," << CL.index << ","
     << CL.index_in_group_reg_inv_list << ")";
  return os;
}

/*****************************************************************************
* A free product of cyclic groups
* ****************************************************************************/
SCABBLE::CyclicProduct::CyclicProduct(void) {
  gens.resize(0);
  orders.resize(0);
}

/* Create a cyclic group from an input of the form a0b2c3, giving Z*Z/2Z*Z/3Z
 */
SCABBLE::CyclicProduct::CyclicProduct(std::string input) {
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

SCABBLE::CyclicProduct::~CyclicProduct(void) {
  //I don't need to free the vectors, right?
}

/* Return the order of a factor, given the generator letter */
int SCABBLE::CyclicProduct::gen_order(char gen) {
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
int SCABBLE::CyclicProduct::index_order(int index) {
  return orders[index];
}

/* return the index of a generator, given the generator letter */
int SCABBLE::CyclicProduct::gen_index(char gen) {
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

int SCABBLE::CyclicProduct::num_groups(void) {
  return gens.size();
}

/* Cyclically minimally reduce the word S.  This entails:
 * (1) cyclically reduce it
 * (2) replace strings of generators with inverses, if that makes it shorter
 * */
void SCABBLE::CyclicProduct::cyc_red(std::string &S) {
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
    
std::ostream& SCABBLE::operator<<(std::ostream &os, SCABBLE::CyclicProduct &G) {
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
SCABBLE::Chain::Chain(void) {
  words.resize(0);
  weights.resize(0);
}



SCABBLE::Chain::Chain(CyclicProduct* G, std::vector<std::string>& words) {
  std::vector<const char*> h(words.size());
  for (int i=0; i<(int)h.size(); ++i) {
    h[i] = words[i].c_str();
  }
  //Chain(G, &h[0], (int)h.size());
}


SCABBLE::Chain::Chain(CyclicProduct* G_in, const char** input, int num_strings) {
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
    //std::cout << "Simplifying " << words[i] << " to:\n";
    (*G).cyc_red(words[i]);
    if ((int)words[i].size() == 0) {
      std::cout << "You gave me a trivial word\n";
      exit(1);
    }
    //std::cout << words[i] << "\n";
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

int SCABBLE::Chain::next_letter(int n) {
  int word = chain_letters[n].word;
  int index = chain_letters[n].index;
  if (index == (int)words[word].size()-1) {
    return n-words[word].size()+1;
  } else {
    return n+1;
  }
}

int SCABBLE::Chain::prev_letter(int n) {
  int word = chain_letters[n].word;
  int index = chain_letters[n].index;
  if (index == 0) {
    return n+words[word].size()-1;
  } else {
    return n-1;
  }
}


int SCABBLE::Chain::num_words(void) {
  return words.size();
}

int SCABBLE::Chain::num_letters() {
  return chain_letters.size();
}

std::string SCABBLE::Chain::operator[](int index) {
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

void SCABBLE::Chain::print_letters(std::ostream &os) {
  int i;
  for (i=0; i<(int)chain_letters.size(); i++) {
    os << i << ": " << chain_letters[i] << "\n";
  }
}

void SCABBLE::Chain::print_group_letters(std::ostream &os) {
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


std::ostream& SCABBLE::operator<<(std::ostream &os, Chain &C) {
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
SCABBLE::CentralEdgePairList::CentralEdgePairList() {
  edge_pairs.resize(0);
  edge_pairs_beginning_with.resize(0);
}

SCABBLE::CentralEdgePairList::CentralEdgePairList(Chain &C) {
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
int SCABBLE::CentralEdgePairList::get_index(int a, int b) {
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

SCABBLE::CentralEdgePair SCABBLE::CentralEdgePairList::operator[](int index) {
  return edge_pairs[index];
}

void SCABBLE::CentralEdgePairList::print(std::ostream &os) {
  int i;
  os << "Central Edge pairs:\n";
  for (i=0; i<(int)edge_pairs.size(); i++) {
    os << i << ": (" << edge_pairs[i].first << ", " << edge_pairs[i].last
       << "), (" << my_chain->prev_letter(edge_pairs[i].last) << ", " << my_chain->next_letter(edge_pairs[i].first) << ")\n";
  }
}


int SCABBLE::CentralEdgePairList::size() {
  return edge_pairs.size();
}
  


/****************************************************************************
 * make a list of all the interface edges
 * note these are (1) any gen with any inverse (2) any (non)inverse gen 
 * with another (non) inverse gen (assuming order > 0)
 * these are from the POLYGON's PERSPECTIVE
 ****************************************************************************/
SCABBLE::InterfaceEdgeList::InterfaceEdgeList() {
  edges.resize(0);
  edges_beginning_with.resize(0);
}

SCABBLE::InterfaceEdgeList::InterfaceEdgeList(Chain &C) {
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

int SCABBLE::InterfaceEdgeList::get_index_from_poly_side(int a, int b) {
  int i;
  for (i=0; i<(int)edges_beginning_with[a].size(); i++) {
    if (edges[edges_beginning_with[a][i]].last == b) {
      return edges_beginning_with[a][i];
    }
  }
  return -1;
}

int SCABBLE::InterfaceEdgeList::get_index_from_group_side(int a, int b) {
  int i;
  for (i=0; i<(int)edges_beginning_with[b].size(); i++) {
    if (edges[edges_beginning_with[b][i]].last == a) {
      return edges_beginning_with[b][i];
    }
  }
  return -1;
}

SCABBLE::InterfaceEdge SCABBLE::InterfaceEdgeList::operator[](int index) {
  return edges[index];
}

void SCABBLE::InterfaceEdgeList::print(std::ostream &os) {
  int i;
  os << "Interface Edges:\n";
  for (i=0; i<(int)edges.size(); i++) {
    os << i << ": " << edges[i].first << ", " << edges[i].last << "\n";
  }
}

int SCABBLE::InterfaceEdgeList::size() {
  return edges.size();
}



/*****************************************************************************
 Central polygons
 *****************************************************************************/
int SCABBLE::CentralPolygon::chi_times_2() {
  return 2 - (int)edges.size();
}

void SCABBLE::CentralPolygon::compute_ia_etc_for_edges(int col, 
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
  
  

std::ostream& SCABBLE::operator<<(std::ostream &os, CentralPolygon &CP) {
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
double SCABBLE::GroupTooth::chi_times_2(Chain &C) {
  int ord = (C.G)->index_order( group_index );
  if (C.next_letter(first) == last) {
    return 2.0/(double)ord;
  } else {
    return 2.0/(double)ord - 1;
  }
}

void SCABBLE::GroupTooth::compute_ia_etc_for_edges(int offset, 
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
  

void SCABBLE::GroupTooth::compute_ia_etc_for_words(int offset, 
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




std::ostream& SCABBLE::operator<<(std::ostream &os, GroupTooth &GT) {
  os << "GT: gp" << GT.group_index << " (" << GT.first << "," << GT.last << ") " << GT.position << " " << GT.base_letter;
  return os;
}



/****************************************************************************
 group rectangle
 ****************************************************************************/
void SCABBLE::GroupRectangle::compute_ia_etc_for_edges(int col, 
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

void SCABBLE::GroupRectangle::compute_ia_etc_for_words(int offset, 
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


std::ostream& SCABBLE::operator<<(std::ostream &os, GroupRectangle &GR) {
  os << "GR: gp" << GR.group_index << " (" << GR.first << "," << GR.last << ")";
  return os;
}


void SCABBLE::print_central_polys(std::vector<SCABBLE::CentralPolygon> &CP, 
                         std::ostream &os, 
                         int level) {
  int i;
  os << "Central polygons: (" << CP.size() << "):\n"; 
  if (level < 3 && (int)CP.size() > 30) {
    os << "(" << CP.size() << " polygons hidden)\n";
  } else {
    for (i=0; i<(int)CP.size(); i++) {
      os << i << ": " << CP[i] << "\n";
    }
  }
}


void SCABBLE::print_group_teeth_and_rectangles(std::vector<SCABBLE::GroupTooth> &GT,
                                  std::vector<SCABBLE::GroupRectangle> &GR,
                                  std::ostream &os,
                                  int level) {
  int i;
  os << "Group teeth: (" << GT.size() << ")\n";
  if (level < 3 && (int)GT.size() > 30) {
    os << "(" << GT.size() << " hidden)\n";
  } else {
    for (i=0; i<(int)GT.size(); i++) {
      os << i << ": " << GT[i] << "\n";
    }  
    os << "Group rectangles: (" << GR.size() << ")\n";
     if (level < 3 && (int)GR.size() > 30) {
      os << "(" << GR.size() << " rectangles hidden)\n";
    } else {
      for (i=0; i<(int)GR.size(); i++) {
        os << i << ": " << GR[i] << "\n";
      }
    }
  }
}


/*****************************************************************************
 * Make the list of group polygons and rectangles. 
 * ***************************************************************************/
void SCABBLE::compute_group_teeth_and_rectangles(SCABBLE::Chain &C, 
                                       std::vector<SCABBLE::GroupTooth > &GT,
                                       std::vector<SCABBLE::GroupRectangle > &GR) {
  int i,j,k,l,m;
  int ord;
  int num_groups = (C.G)->num_groups();
  SCABBLE::GroupTooth temp_group_tooth;
  SCABBLE::GroupRectangle temp_group_rect;
  
  GT.resize(0);
  GR.resize(0);
  
  for (i=0; i<num_groups; i++) {
    ord = (C.G)->index_order(i);
    
    //compute the group rectangles
    temp_group_rect.group_index = i;
    for (j=0; j<(int)C.regular_letters[i].size(); j++) {
      temp_group_rect.first = C.regular_letters[i][j];
      for (k=0; k<(int)C.inverse_letters[i].size(); k++) {
        temp_group_rect.last = C.inverse_letters[i][k];
        GR.push_back(temp_group_rect);
      }
    }
    
    //compute the group teeth
    if (ord == 0) {
      continue;
    }
    //non-inverses 
    temp_group_tooth.inverse = false;
    temp_group_tooth.group_index = i;
    for (j=0; j<(int)C.regular_letters[i].size(); j++) {   //this is the basepoint
      temp_group_tooth.base_letter = C.regular_letters[i][j];
      //do the group teeth at position 0
      temp_group_tooth.position = 0;
      temp_group_tooth.first = C.regular_letters[i][j];
      for (k=0; k<(int)C.regular_letters[i].size(); k++) {
        temp_group_tooth.last = C.regular_letters[i][k];
        GT.push_back(temp_group_tooth);
      }
      //do the rest of the positions (=k)
      for (k=1; k<ord-1; k++) {                              
        temp_group_tooth.position = k;
        for (l=0; l<(int)C.regular_letters[i].size(); l++) {
          temp_group_tooth.first = C.regular_letters[i][l];
          for (m=0; m<(int)C.regular_letters[i].size(); m++) {
            temp_group_tooth.last = C.regular_letters[i][m];
            GT.push_back(temp_group_tooth);
          }
        }
      }
      //do position ord-1
      temp_group_tooth.position = ord-1;
      temp_group_tooth.last = C.regular_letters[i][j];
      for (k=0; k<(int)C.regular_letters[i].size(); k++) {
        temp_group_tooth.first = C.regular_letters[i][k];
        GT.push_back(temp_group_tooth);
      }
    }
    
    //now inverses
    temp_group_tooth.inverse = true;
    temp_group_tooth.group_index = i;
    for (j=0; j<(int)C.inverse_letters[i].size(); j++) {   //this is the basepoint
      temp_group_tooth.base_letter = C.inverse_letters[i][j];
      //do the group teeth at position 0
      temp_group_tooth.position = 0;
      temp_group_tooth.first = C.inverse_letters[i][j];
      for (k=0; k<(int)C.inverse_letters[i].size(); k++) {
        temp_group_tooth.last = C.inverse_letters[i][k];
        GT.push_back(temp_group_tooth);
      }
      //do the rest of the positions (=k)
      for (k=1; k<ord-1; k++) {                              
        temp_group_tooth.position = k;
        for (l=0; l<(int)C.inverse_letters[i].size(); l++) {
          temp_group_tooth.first = C.inverse_letters[i][l];
          for (m=0; m<(int)C.inverse_letters[i].size(); m++) {
            temp_group_tooth.last = C.inverse_letters[i][m];
            GT.push_back(temp_group_tooth);
          }
        }
      }
      //do position ord-1
      temp_group_tooth.position = ord-1;
      temp_group_tooth.last = C.inverse_letters[i][j];
      for (k=0; k<(int)C.inverse_letters[i].size(); k++) {
        temp_group_tooth.first = C.inverse_letters[i][k];
        GT.push_back(temp_group_tooth);
      }
    }
  }
  
}




/*****************************************************************************
 * compute the list of polys. 
 Setting limit_central_sides limits the central sides to only a single 
 one in a triangle.  This effectively duplicates the nonrigorous scallop
 (the intention is that it should be faster)
 *****************************************************************************/
void SCABBLE::compute_central_polys(SCABBLE::Chain &C, 
                                    SCABBLE::InterfaceEdgeList &IEL, 
                                    std::vector<SCABBLE::CentralPolygon> &CP) {
  int i,j,k,l;
  int e1L1, e1L2, e2L1, e2L2, e3L1, e3L2;  //edge1, letter1, etc
  SCABBLE::CentralPolygon temp_central_poly;
  
  CP.resize(0);
  
  //first, enumerate all polys with two sides.  These are 
  //always both interface edges.  We may always assume that 
  //the smallest letter is at position 0
  temp_central_poly.edges.resize(2);
  temp_central_poly.interface = std::vector<bool>(2, true);
  for (i=0; i<(int)C.chain_letters.size(); i++) {               //the first letter
    e1L1 = i;
    for (j=0; j<(int)IEL.edges_beginning_with[i].size(); j++) { //which edge we are thinking about
      e1L2 = IEL[IEL.edges_beginning_with[i][j]].last;
      e2L1 = C.next_letter(e1L2);
      if (e2L1 < e1L1) {
        continue;
      }
      for (k=0; k<(int)IEL.edges_beginning_with[e2L1].size(); k++) {
        e2L2 = IEL[IEL.edges_beginning_with[e2L1][k]].last;
        if (e2L2 == C.prev_letter(e1L1)) {
          temp_central_poly.edges[0] = std::make_pair( e1L1, e1L2 );
          temp_central_poly.edges[1] = std::make_pair( e2L1, e2L2 );
          CP.push_back(temp_central_poly);
          break;
        }
      }
    }
  }
  
  //enumerate all polys with three sides, all interface
  temp_central_poly.edges.resize(3);
  temp_central_poly.interface = std::vector<bool>(3, true);
  for (i=0; i<(int)C.chain_letters.size(); i++) {               //the first letter
    e1L1 = i;
    for (j=0; j<(int)IEL.edges_beginning_with[i].size(); j++) { //which edge we are thinking about
      e1L2 = IEL[IEL.edges_beginning_with[i][j]].last;
      e2L1 = C.next_letter(e1L2);
      if (e2L1 < e1L1) {  //if first letter isn't smallest
        continue;
      }
      for (k=0; k<(int)IEL.edges_beginning_with[e2L1].size(); k++) {
        e2L2 = IEL[IEL.edges_beginning_with[e2L1][k]].last;
        e3L1 = C.next_letter(e2L2);      
        if (e3L1 < e1L1) {  //if first letter isn't smallest
          continue;
        }
        for (l=0; l<(int)IEL.edges_beginning_with[e3L1].size(); l++) {
          e3L2 = IEL[IEL.edges_beginning_with[e3L1][l]].last;
          if (e3L2 == C.prev_letter(e1L1)) {
            temp_central_poly.edges[0] = std::make_pair( e1L1, e1L2 );
            temp_central_poly.edges[1] = std::make_pair( e2L1, e2L2 );
            temp_central_poly.edges[2] = std::make_pair( e3L1, e3L2 );
            CP.push_back(temp_central_poly);
            break;
          }
        }
      }
    }
  }
  
  //enumerate all polys with three sides, 2 interface
  //we can no longer assume the first one is the smallest letter
  temp_central_poly.edges.resize(3);
  temp_central_poly.interface = std::vector<bool>(3, true);
  temp_central_poly.interface[2] = false;
  for (i=0; i<(int)C.chain_letters.size(); i++) {               //the first letter
    e1L1 = i;
    for (j=0; j<(int)IEL.edges_beginning_with[i].size(); j++) { //which edge we are thinking about
      e1L2 = IEL[IEL.edges_beginning_with[i][j]].last;
      e2L1 = C.next_letter(e1L2);   
      for (k=0; k<(int)IEL.edges_beginning_with[e2L1].size(); k++) {
        e2L2 = IEL[IEL.edges_beginning_with[e2L1][k]].last;
        if (C.next_letter(e2L2) == e1L1) { //this is really a bigon
          continue;
        }
        temp_central_poly.edges[0] = std::make_pair( e1L1, e1L2 );
        temp_central_poly.edges[1] = std::make_pair( e2L1, e2L2 );
        temp_central_poly.edges[2] = std::make_pair( e2L2, e1L1 );
        CP.push_back(temp_central_poly);
      }
    }
  }
  
  //if (limit_central_sides) {
  //  return;
  //}
  
  //enumerate all polys with three sides, 1 interface
  //we can no longer assume the first one is the smallest letter
  temp_central_poly.edges.resize(3);
  temp_central_poly.interface = std::vector<bool>(3, false);
  temp_central_poly.interface[0] = true;
  for (i=0; i<(int)C.chain_letters.size(); i++) {               //the first letter
    e1L1 = i;
    for (j=0; j<(int)IEL.edges_beginning_with[i].size(); j++) { //which edge we are thinking about
      e1L2 = IEL[IEL.edges_beginning_with[i][j]].last;
      e2L1 = e1L2;     
      for (e2L2 = 0; e2L2<(int)C.chain_letters.size(); e2L2++) {
        if (e2L2 == C.next_letter(e2L1)) {
          continue;
        }
        e3L1 = C.prev_letter(e2L2); //this is really weird
        e3L2 = e1L1;
        if (e3L2 == C.next_letter(e3L1)) {
          continue;
        }
        temp_central_poly.edges[0] = std::make_pair( e1L1, e1L2 );
        temp_central_poly.edges[1] = std::make_pair( e2L1, e2L2 );
        temp_central_poly.edges[2] = std::make_pair( e3L1, e3L2 );
        CP.push_back(temp_central_poly);       
      }
    }
  }
  
}




//this create the LP -- everything except the final row
void SCABBLE::init_orthant_lp(std::vector<std::pair<int, int> >& chain_locs,
                              SCABBLE::Chain& C, 
                              SCABBLE::InterfaceEdgeList& IEL, 
                              SCABBLE::CentralEdgePairList& CEL,
                              std::vector<SCABBLE::CentralPolygon>& CP, 
                              std::vector<SCABBLE::GroupTooth>& GT, 
                              std::vector<SCABBLE::GroupRectangle>& GR,
                              SparseLP& LP,
                              SparseLPSolver solver,
                              std::vector<std::vector<int> >& chain_cols,
                              std::vector<int>& chain_rows,
                              int verbose) {
  int i,j,k,m;
  int ord;
  int num_cols, offset, num_rows;
  int num_equality_rows;
  int i_edge_pairs = IEL.size();
  int c_edge_pairs = CEL.size();
  int num_words = C.num_words();
  //ROWS (see above)
  //we need to construct something to tell us the row 
  //number for the group teeth rows
  // group_teeth_rows_reg[i][j][k] gives the row of 
  // group i, letter j base letter k, position 1.
  // all other positions must offset from there.
  std::vector<std::vector<std::vector<int> > > group_teeth_rows_reg;
  std::vector<std::vector<std::vector<int> > > group_teeth_rows_inv;
  num_rows = i_edge_pairs + c_edge_pairs;
  group_teeth_rows_reg.resize((C.G)->num_groups());
  group_teeth_rows_inv.resize((C.G)->num_groups());
  for (i=0; i<(int)(C.G)->num_groups(); i++) {
    ord = (C.G)->index_order(i);
    if (ord == 0) {
      group_teeth_rows_reg[i].resize(0);
      group_teeth_rows_inv[i].resize(0);
      continue;
    }
    group_teeth_rows_reg[i].resize(C.regular_letters[i].size());
    for (j=0; j<(int)C.regular_letters[i].size(); j++) {
      group_teeth_rows_reg[i][j].resize(C.regular_letters[i].size());
      for (k=0; k<(int)C.regular_letters[i].size(); k++) {
        group_teeth_rows_reg[i][j][k] = num_rows;
        num_rows += ord-1;
      }
    }
    group_teeth_rows_inv[i].resize(C.inverse_letters[i].size());
    for (j=0; j<(int)C.inverse_letters[i].size(); j++) {
      group_teeth_rows_inv[i][j].resize(C.inverse_letters[i].size());
      for (k=0; k<(int)C.inverse_letters[i].size(); k++) {
        group_teeth_rows_inv[i][j][k] = num_rows;
        num_rows += ord-1;
      }
    }    
  }   
  
  num_equality_rows = num_rows;
  num_rows += (num_words-(int)chain_locs.size());
  num_rows += chain_locs.size();
  
  num_cols = CP.size() + GT.size() + GR.size();
  
  if (verbose > 2) {
    std::cout << "Num equality rows: " << num_equality_rows
    << "Num rows: " << num_rows << "\n";
  }
  
  if (verbose>1) {
    std::cout << "Started linear programming setup\n";
  }
  
  
  //Create the LP problem
  LP = SparseLP(solver, num_rows, num_cols);
  
  for(i=0; i<(int)num_equality_rows; i++){
    LP.set_RHS(i, 0); // glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
    LP.set_equality_type(i, EQ);
    if (verbose > 3) {
      std::cout << "Set row " << i << " bounded to " << 0 << "\n";
    }
  }
  //we need to set all words as appearing as many times as the first word
  for(i=0; i<(int)num_words; i++){
    //RHS[num_equality_rows+i] = C.weights[i];
    LP.set_RHS(num_equality_rows+i, 0);
    LP.set_equality_type(num_equality_rows+i, EQ);
    //glp_set_row_bnds(lp, 
    //                  num_equality_rows+i+1, 
    //                  GLP_FX, 
    //                  C.weights[i], 
    //                  C.weights[i]);  
    if (verbose > 3) {
      std::cout << "Set row " << num_equality_rows+i << " bounded to " << 0 << "\n";
    }
  }
  
  
  //COLS
  for(i=0; i<(int)CP.size(); i++){
    //glp_set_col_bnds(lp, i+1, GLP_LO, 0.0, 0.0);
    //objective[i] = -CP[i].chi_times_2(); //glp_set_obj_coef(lp, i+1, -CP[i].chi_times_2());
    LP.set_obj(i, -CP[i].chi_times_2());
    if (verbose>3) {
      std::cout << "Set objective " << i << " to " << -CP[i].chi_times_2() << "\n";
    }
  }
  offset = CP.size();
  for (i=0; i<(int)GT.size(); i++) {
    //glp_set_col_bnds(lp, offset+i+1, GLP_LO, 0.0, 0.0);
    //objective[offset + i] = -GT[i].chi_times_2(C); //glp_set_obj_coef(lp, offset+i+1, -GT[i].chi_times_2(C));
    LP.set_obj(offset+i, -GT[i].chi_times_2(C));
    if (verbose>3) {
      std::cout << "GT Set objective " << offset+i << " to " << -GT[i].chi_times_2(C) << "\n";
    }
  }
  offset = CP.size() + GT.size();
  for (i=0; i<(int)GR.size(); i++) {
    //glp_set_col_bnds(lp, offset+i+1, GLP_LO, 0.0, 0.0);
    //objective[offset+i] = 0; //glp_set_obj_coef(lp, offset+i+1, 0);
    LP.set_obj(offset+i, 0);
    if (verbose>3) {
      std::cout << "GR Set objective " << offset+i << " to " << 0 << "\n";
    }
  }    
  
  //CENTRAL POLYGONS
  for (i=0; i<(int)CP.size(); i++) {
    CP[i].compute_ia_etc_for_edges(i,
                                   C,
                                   IEL, 
                                   CEL, 
                                   LP);
    if (verbose > 3) {
      std::cout << "CP number " << i << ":\n";
    }
  }
  
  if (verbose>2) { 
    std::cout << "Loaded central polygon edge constraints\n";
  }     
  
  //GROUP TEETH and RECTANGLES
  offset = CP.size();
  for (m=0; m<(int)GT.size(); m++) {
    if (GT[m].inverse) {
      if (verbose > 2) {
        std::cout << GT[m] << "\n";
      }
      GT[m].compute_ia_etc_for_edges(offset + m, 
                                     C, 
                                     IEL, 
                                     group_teeth_rows_inv, 
                                     LP);
    } else {        
      if (verbose > 3) {
        std::cout << GT[m] << "\n";
      }
      GT[m].compute_ia_etc_for_edges(offset + m, 
                                     C, 
                                     IEL, 
                                     group_teeth_rows_reg, 
                                     LP);
    }
  }
  offset = CP.size() + GT.size();
  for (m=0; m<(int)GR.size(); m++) {
    if (verbose>3) {
      std::cout << GR[m] << "\n";
    }
    GR[m].compute_ia_etc_for_edges(offset + m, IEL, LP);
  }
  
  
  if (verbose>2) {
    std::cout << "Loaded group constraints\n";
  }
  
  //word constraints: for every word, ensure that it appears as many times
  //as the first word in that chain
  //first, record which columns are the special ones
  //(are the first letter of the first word in a chain
  chain_cols.resize(chain_locs.size());
  for (i=0; i<(int)chain_locs.size(); ++i) {
    chain_cols[i].resize(0);
  }
  chain_rows.resize(chain_locs.size());
  offset = CP.size();
  for (i=0; i<(int)GT.size(); ++i) {
    for (j=0; j<(int)chain_locs.size(); ++j) {
      if (C.chain_letters[GT[i].first].word == chain_locs[j].first &&
          C.chain_letters[GT[i].first].index == 0) {
        chain_cols[j].push_back( offset + i );
        break;
      }
    }
  }
  offset = CP.size() + GT.size();
  for (i=0; i<(int)GR.size(); ++i) {
    for (j=0; j<(int)chain_locs.size(); ++j) {
      //std::cout << "Looking at GR " << i << " with chain " << j <<"\n";
      //std::cout << "GR: " << GR[i] << "\n";
      //std::cout << "Letter: " << C.chain_letters[GR[i].first] << "\n";
      if (C.chain_letters[GR[i].first].word == chain_locs[j].first &&
          C.chain_letters[GR[i].first].index == 0) {
        //std::cout << "Yep this is a chain column\n";
        chain_cols[j].push_back( offset + i );
        break;
      }
    }
  }  
  
  if (verbose > 3) {
    std::cout << "Chain columns:\n";
    for (i=0; i<(int)chain_cols.size(); ++i) {
      std::cout << i << ": ";
      for (j=0; j<(int)chain_cols[i].size(); ++j) {
        std::cout << chain_cols[i][j] << " ";
      }
      std::cout << "\n";
    }
  }
    
  //now make the constraints -- for all words, they must appear the same
  //number of times as the first word in that chain
  //only put the -1 for non-first words
  
  //get a list of which words are main words
  std::vector<bool> is_main_word(num_words, false);
  for (i=0; i<(int)chain_locs.size(); ++i) {
    is_main_word[chain_locs[i].first] = true;
  }
  is_main_word.push_back(true); //<-- this makes all the loops below stop
  
  int row_offset = num_equality_rows;
  offset = CP.size();
  for (i=0; i<(int)GT.size(); ++i) {
    int wi1 = C.chain_letters[GT[i].first].index;
    if (wi1 != 0) continue;
    int w1 = C.chain_letters[GT[i].first].word;
    if (is_main_word[w1]) continue;
    LP.add_entry( row_offset + w1,
                  offset + i,
                 -1 );
  }
  
  offset = CP.size() + GT.size();
  //for (i=0; i<num_words; ++i) {
  //  std::cout << is_main_word[i] << " ";
  //}
  //std::cout << "\n";
  //now enter the word columns
  //this is just horrible horrible
  for (i=0; i<(int)GR.size(); ++i) {
    int w1 = C.chain_letters[GR[i].first].word;
    int wi1 = C.chain_letters[GR[i].first].index;
    int w2 = C.chain_letters[GR[i].last].word;
    int wi2 = C.chain_letters [GR[i].last].index;
    if (wi1 != 0 && wi2 != 0) continue;
    if (wi1 == 0 && wi2 != 0) { //only w1 is active
      if (is_main_word[w1]) {
        for (j=w1+1; !is_main_word[j]; ++j) {
          LP.add_entry(row_offset + j,
                       offset + i,
                       1);
        }
      } else {
        LP.add_entry(row_offset + w1,
                     offset + i,
                     -1);
      }
      
    } else if (wi1 != 0 && wi2 == 0) {
      if (is_main_word[w2]) {
        for (j=w2+1; !is_main_word[j]; ++j) {
          LP.add_entry(row_offset + j,
                       offset + i,
                       1);
        }
      } else {
        LP.add_entry(row_offset + w2,
                     offset + i,
                     -1);
      }
      
    } else if (wi1 == 0 && wi2 == 0) {
      if (is_main_word[w1] && is_main_word[w2]) {
        if (w1 == w2) {
          for (j=w1; !is_main_word[j]; ++j) {
            LP.add_entry(row_offset + j,
                         offset + i,
                         2);
          }
        } else {
          for (j=w1; !is_main_word[j]; ++j) {
            LP.add_entry(row_offset + j,
                         offset + i,
                         1);
          }
          for (j=w2; !is_main_word[j]; ++j) {
            LP.add_entry(row_offset + j,
                         offset + i,
                         1);
          }
        }
        
      } else if (is_main_word[w1]) { //w1 is main, w2 not
        bool did_w2 = false;
        for (j=w1+1; !is_main_word[j]; ++j) {
          if (j==w2) {
            did_w2 = true;
            continue;
          }
          LP.add_entry(row_offset + j,
                       offset + i,
                       1);
        }
        if (!did_w2) LP.add_entry(row_offset + w2, offset + i, -1);
      
      } else if (is_main_word[w2]) { //w2 is main, w1 not
        bool did_w1 = false;
        for (j=w2+1; !is_main_word[j]; ++j) {
          if (j==w1) {
            did_w1 = true;
            continue;
          }
          LP.add_entry(row_offset + j,
                       offset + i,
                       1);
        }
        if (!did_w1) LP.add_entry(row_offset + w1, offset + i, -1);
      
      } else { //neither are main
        if (w1 == w2) {
          LP.add_entry(row_offset + w1, offset + i, -2);
        } else {
          LP.add_entry(row_offset + w1, offset + i, -1);
          LP.add_entry(row_offset + w2, offset + i, -1);
        }
      }
    }
  }
  
  //find the chain rows
  chain_rows.resize(chain_locs.size());
  for (i=0; i<(int)chain_locs.size(); ++i) {
    chain_rows[i] = row_offset + chain_locs[i].first;
  }
  
  if (verbose > 3) {
    std::cout << "Chain rows: \n";
    for (i=0; i<(int)chain_rows.size(); ++i) {
      std::cout << chain_rows[i] << " ";
    }
    std::cout << "\n";
  }
  
  if (verbose > 2) {
    std::cout << "Loaded word constraints\n";
  }
  
  if (verbose > 3) {
    std::cout << "LP problem:\n";
    LP.print_LP();
  }
  
}


//compute the scl at a particular point
void SCABBLE::point_scl(std::vector<std::vector<int> >& chain_cols,
                        std::vector<int>& chain_rows,
                        SparseLP& LP,
                        SCABBLE::Pt& p,
                        Rational& scl,
                        int verbose) {
  int init_matrix_entries = LP.get_num_entries();
  
  if (verbose > 2) {
    std::cout << "Computing the scl of the point: " << p << "\n";
  }
  
  //just set the key rows
  for (int i=0; i<(int)chain_cols.size(); ++i) {
    for (int j=0; j<(int)chain_cols[i].size(); ++j) {
      LP.add_entry( chain_rows[i], chain_cols[i][j], 1 );
    }
    LP.set_RHS(chain_rows[i], p[i]);
  }
  
  if (verbose > 3) {
    std::cout << "Running the LP problem:\n";
    LP.print_LP();
  }
  
  LP.solve(verbose);
  
  LP.get_optimal_value(scl);
  
  if (verbose > 2) {
    std::cout << "Got scl " << scl << "\n";
  }
  
  LP.reset_num_entries(init_matrix_entries);
  
}

//compute the min scl over a face
void SCABBLE::face_scl(std::vector<std::vector<int> >& chain_cols,
                       std::vector<int>& chain_rows,
                       SparseLP& LP,
                       std::vector<SCABBLE::Pt>& face,
                       Rational& scl,
                       SCABBLE::Pt& min_p,
                       int verbose) {
  int init_matrix_entries = LP.get_num_entries();
  
  if (verbose > 2) {
    std::cout << "Computing the scl of the face: \n";
    for (int i=0; i<(int)face.size(); ++i) {
      std::cout << face[i] << "\n";
    }
  }
  
  SCABBLE::Pt normal((int)chain_cols.size());
  Rational normal_value;
  affine_hyperplane(face, normal, normal_value);  
  
  if (verbose > 2) {
    std::cout << "Got the normal vector: " << normal << "\n";
    std::cout << "With normal value: " << normal_value << "\n";
  }
  
  //set the first row, and make everything else zero
  for (int i=0; i<(int)chain_cols.size(); ++i) {
    for (int j=0; j<(int)chain_cols[i].size(); ++j) {
      LP.add_entry( chain_rows[0], chain_cols[i][j], normal[i] );
    }
  }
  LP.set_RHS(chain_rows[0], normal_value);
  
  for (int i=1; i<(int)chain_cols.size(); ++i) {
    LP.set_RHS(chain_rows[i], 0);
  }
  
  if (verbose > 3) {
    std::cout << "Running the LP problem:\n";
    LP.print_LP();
  }
    
  LP.solve(verbose);
  
  LP.get_optimal_value(scl);
  
  if (verbose > 2) {
    std::cout << "Got scl " << scl << "\n";
  }
  
  std::vector<Rational> soln_vector(0);
  LP.get_soln_vector(soln_vector);
  for (int i=0; i<(int)chain_cols.size(); ++i) {
    Rational sum = 0;
    for (int j=0; j<(int)chain_cols[i].size(); ++j) {
      sum += soln_vector[chain_cols[i][j]];
    }
    min_p[i] = sum;
  }
  
  if (verbose > 2) {
    std::cout << "Got min point: " << min_p << "\n";
  }
  
  LP.reset_num_entries(init_matrix_entries);
    
}

void SCABBLE::subdivide_face(std::vector<SCABBLE::Pt>& working_face, 
                             SCABBLE::Pt& new_vert, 
                             std::vector<std::vector<SCABBLE::Pt> >& face_stack, 
                             int verbose) {
  int dim = new_vert.dim();
  std::vector<SCABBLE::Pt> new_face(dim);
  for (int i=0; i<dim; ++i) {
    for (int j=0; j<dim; ++j) {
      if (j==i) {
        new_face[j] = new_vert;
      } else {
        new_face[j] = working_face[j];
      }
    }
    face_stack.push_back(new_face);
  }
}


//computes an (arbitrary dimenion) positive orthant of the ball
void SCABBLE::compute_ball_ant(std::vector<std::pair<int, int> >& chain_locs,
                               SCABBLE::Chain& C, 
                               std::vector<SCABBLE::Pt>& orthant_verts, 
                               std::vector<std::vector<SCABBLE::Pt> >& orthant_faces, 
                               SparseLPSolver solver,
                               int verbose) {
  int dim = (int)chain_locs.size();
  
  //now we build the rectangles and whatever for the *whole* chain
  //we'll restrict later
  SCABBLE::InterfaceEdgeList IEL(C);
  if (verbose>3) IEL.print(std::cout);
  
  SCABBLE::CentralEdgePairList CEL(C);
  if (verbose>3) CEL.print(std::cout);
  
  std::vector<SCABBLE::CentralPolygon> CP;
  SCABBLE::compute_central_polys(C, IEL, CP);
  if (verbose > 3) {
    std::cout << "computed polys (" << CP.size() << ")\n"; std::cout.flush();
    SCABBLE::print_central_polys(CP, std::cout, verbose);
  }
  
  std::vector<SCABBLE::GroupTooth> GT;
  std::vector<SCABBLE::GroupRectangle> GR;
  SCABBLE::compute_group_teeth_and_rectangles(C, GT, GR);
  if (verbose > 3) {
    std::cout << "computed group teeth and rectangles\n"; std::cout.flush();
    SCABBLE::print_group_teeth_and_rectangles(GT, GR, std::cout, verbose);
  }
  
  SparseLP LP(solver);
  std::vector<std::vector<int> > chain_cols(0);
  std::vector<int> chain_rows(0);
  SCABBLE::init_orthant_lp(chain_locs, 
                           C, IEL, CEL, CP, GT, GR, 
                           LP, solver, 
                           chain_cols, chain_rows, 
                           verbose);
  
  //find the scls of the basis vectors
  std::vector<SCABBLE::Pt> initial_face(dim);
  for (int i=0; i<dim; ++i) {
    initial_face[i] = SCABBLE::Pt(dim, 0);
    initial_face[i][i] = 1;
    Rational s;
    point_scl(chain_cols, chain_rows, LP, initial_face[i], s, verbose);
    initial_face[i][i] /= s;
  }
  
  if (verbose > 2) {
    std::cout << "Computed scl of basis vectors\n";
    std::cout << "Initial face:\n";
    for (int i=0; i<(int)initial_face.size(); ++i) {
      std::cout << i << ": " << initial_face[i] << "\n";
    }
  }
  
  //create the stack
  std::vector<std::vector<SCABBLE::Pt> > face_stack(0);
  face_stack.push_back(initial_face);
  
  //enter the main loop
  orthant_verts.resize(0);
  orthant_faces.resize(0);
  while (face_stack.size() > 0) {    
    std::vector<SCABBLE::Pt> working_face = face_stack.back();
    face_stack.pop_back();
    
    if (verbose > 2) {
      std::cout << "Remaining stack size: " << face_stack.size() << "\n";
      std::cout << "Working face: \n";
      for (int i=0; i<(int)working_face.size(); ++i) {
        std::cout << working_face[i] << "\n";
      }
    }
    
    Rational s;
    SCABBLE::Pt new_vert(dim,0);
    
    face_scl(chain_cols, chain_rows, LP, working_face, s, new_vert, verbose);
    
    if (s == 1) {
      orthant_faces.push_back(working_face);
      if (verbose > 2) {
        std::cout << "This face is good\n";
      }
      continue;
    }
   
    new_vert /= s;
   
    if (verbose > 2) {
      std::cout << "new vertex: " << new_vert << "\n";
    }
   
    //now we add the new simplices to the face.
    subdivide_face(working_face, new_vert, face_stack, verbose);
  }
  
}



void SCABBLE::compute_ball( SCABBLE::CyclicProduct& G,
                            std::vector<std::string>& words,
                            std::vector<std::pair<int,int> >& chain_locs,
                            std::vector<SCABBLE::Pt>& verts,
                            std::vector<std::vector<SCABBLE::Pt> >& faces,
                            SparseLPSolver solver,
                            int verbose ) {
  
  int dim = (int)chain_locs.size();
  
  verts.resize(0);
  faces.resize(0);
  
  //we run through all the orthants that we need to
  int num_orthants = (1<<dim);
  int orthant_mask = 0;
  for (int i=0; i<dim; ++i) {
    orthant_mask |= (1<<i);
  }
  std::vector<int> orthant_done(num_orthants, 0);
  
  std::vector<SCABBLE::Pt> orthant_verts(0);
  std::vector<std::vector<SCABBLE::Pt> > orthant_faces(0);
  std::vector<std::string> orthant_words(words.size());
  
  for (int i=0; i<num_orthants; ++i) {
    
    if (verbose > 1) {
      std::cout << "Running orthant " << i << "\n";
    }
    
    if (orthant_done[i] != 0) continue;
    
    //invert the chains as necessary
    for (int j=0; j<dim; ++j) {
      bool invert = ( ((i>>j)&1) == 1 ? true : false );  //invert this chain?
      if (invert) {
        for (int k=chain_locs[j].first; k<chain_locs[j].first + chain_locs[j].second; ++k) {
          orthant_words[k] = inverse(words[k]);
        }
      } else {
        for (int k=chain_locs[j].first; k<chain_locs[j].first + chain_locs[j].second; ++k) {
          orthant_words[k] = words[k];
        }
      }
    }
    
    //load the chain from the orthant words
    std::vector<const char*> c_words(words.size());
    for (int j=0; j<(int)words.size(); ++j) {
      c_words[j] = orthant_words[j].c_str();
    }
    SCABBLE::Chain C(&G, &c_words[0], (int)c_words.size());   //process the chain argument
    
    //print if verbose
    if (verbose>1) {
      std::cout << "Group: " << G << "\n";
      std::cout << "Chain: " << C << "\n";
      if (verbose>2) {
        std::cout << "Letters:\n";
        C.print_letters(std::cout);
        std::cout << "Group letters:\n";
        C.print_group_letters(std::cout);
      }
    }
    
    //now produce the positive orthant
    SCABBLE::compute_ball_ant(chain_locs,
                              C, 
                              orthant_verts, 
                              orthant_faces, 
                              solver,
                              verbose);
    
    //append the verts and faces to the list
    //we have to remember to negate the coordinates which were inverted
    for (int j=0; j<(int)orthant_verts.size(); ++j) {
      verts.push_back(orthant_verts[j].negate_coords(i));
    }
    if (verbose > 1) {
      std::cout << "Finished orthant " << i <<"\n";
      std::cout << "With faces:\n";
    }
    for (int j=0; j<(int)orthant_faces.size(); ++j) {
      std::vector<SCABBLE::Pt> temp_face(0);
      for (int k=0; k<(int)orthant_faces[j].size(); ++k) {
        temp_face.push_back(orthant_faces[j][k].negate_coords(i));
      }
      if (verbose > 1) {
        std::cout << "Face " << j << ":\n";
        for (int k=0; k<(int)temp_face.size(); ++k) {
          std::cout << temp_face[k] << "\n";
        }
      }
      faces.push_back(temp_face);
    }
    
    //mark this orthant, and its reflection, as done
    orthant_done[i] = 1;
    orthant_done[ (~i)&orthant_mask ] = 1;
    
  }
  
  //now we must reflect every face and vertex through the origin
  int old_num_verts = (int)verts.size();
  for (int i=0; i<old_num_verts; ++i) {
    verts.push_back(-verts[i]);
  }
  int old_num_faces = (int)faces.size();
  for (int i=0; i<old_num_faces; ++i) {
    std::vector<SCABBLE::Pt> temp_face(faces[i].size());
    for (int j=0; j<(int)faces[i].size(); ++j) {
      temp_face[j] = -faces[i][j];
    }
    faces.push_back(temp_face);
  }
  
}
    
  
  
  


void SCABBLE::write_ball_to_file(std::string output_filename, 
                                  std::vector<SCABBLE::Pt>& verts, 
                                  std::vector<std::vector<SCABBLE::Pt> >& faces, 
                                  int verbose) {
  std::cout << "Sorry, writing the ball to a file is currently not implemented\n";
}

void SCABBLE::draw_ball_to_file(std::string output_filename, 
                                std::vector<SCABBLE::Pt>& verts, 
                                std::vector<std::vector<SCABBLE::Pt> >& faces, 
                                int verbose) {
    
  if (verbose > 1) {
    std::cout << "Drawing ball to file: " << output_filename << "\n";
  }
  
  std::fstream fs;
  fs.open(output_filename.c_str(), std::fstream::out);
  fs << "%!PS-Adobe-2.0 EPSF-2.0\n";
  fs << "%%BoundingBox: 0 0 800 800\n";
  fs << "400 400 translate\n";
  fs << "0 0 0 setcolor\n";
  fs << "/Arial findfont\n";
  fs << "14 scalefont\n";
  fs << "setfont\n";
  
  //first, find the x and y scaling factors
  double max_x = 0;
  double max_y = 0;
  for (int i=0; i<(int)faces.size(); ++i) {
    for (int j=0; j<(int)faces[i].size(); ++j) {
      double x = fabs(faces[i][j][0].get_d());
      double y = fabs(faces[i][j][1].get_d());
      if (x > max_x) max_x = x;
      if (y > max_y) max_y = y;
    }
  }
  //make the scaling factor so that max_x and max_y would be about 90
  double scaling_x, scaling_y;
  scaling_x = 350.0/max_x;
  scaling_y = 350.0/max_y;
  //now the drawing coords of a point are (100,100) + (scaling_x*x, scaling_y*y)
  for (int i=0; i<(int)faces.size(); ++i) {
    double start_x = scaling_x*faces[i][0][0].get_d();
    double start_y = scaling_y*faces[i][0][1].get_d();
    double end_x = scaling_x*faces[i][1][0].get_d();
    double end_y = scaling_y*faces[i][1][1].get_d();
    //draw the line
    fs << "1 setlinewidth\n";
    fs << start_x << " " << start_y << " moveto\n";
    fs << end_x << " " << end_y << " lineto\n";
    fs << "stroke\n";
    
    //draw the vertices (they'll get drawn twice)    
    fs << "3 setlinewidth\n";
    fs << start_x << " " << start_y << " moveto\n";
    fs << start_x << " " << start_y << " 1.5 0 360 arc\n";
    fs << "stroke\n";
    fs << end_x << " " << end_y << " moveto\n";
    fs << end_x << " " << end_y << " 1.5 0 360 arc\n";
    fs << "stroke\n\n";
    
    //put the text labels on
    fs << "1 setlinewidth\n";
    fs << start_x << " " << start_y << " moveto\n";
    fs << "(  (" << faces[i][0][0] << ", " << faces[i][0][1] << ")) show\n";
    fs << end_x << " " << end_y << " moveto\n";
    fs << "(  (" << faces[i][1][0] << ", " << faces[i][1][1] << ")) show\n";
  }
  
  //now draw the axes
  fs << "0.5 setlinewidth\n";
  fs << "-398 0 moveto\n";
  fs << "398 0 lineto\n";
  fs << "stroke\n";
  fs << "0 -398 moveto\n";
  fs << "0 398 lineto\n";
  fs << "stroke\n";  
  
  fs.close();
  
  if (verbose > 1) {
    std::cout << "Done drawing\n";
  }
  
}







int SCABBLE::scabble(int argc, char** argv) {
  
  SparseLPSolver solver = EXLP;
  bool output_polyhedron = false;
  std::string output_filename = "";
  int verbose = 1;
  int current_arg = 0;
  
  if (argc < 1 || std::string(argv[0]) == "-h") {
    std::cout << "usage: ./scallop -ball [-h] [-v[n]] [-P] [-m<GLPK,GIPT,EXLP,GUROBI>] <filename> [gen string] <chain1> , <chain2> , ...\n";
    std::cout << "\twhere [gen string] is of the form <gen1><order1><gen2><order2>...\n";
    std::cout << "\te.g. a5b0 computes in Z/5Z * Z\n";
    std::cout << "\tand the <chain>s are integer linear combinations of words in the generators\n";
    std::cout << "\tthe commas are necessary.  There is a maximum of 2 chains unless -P is specified\n";
    std::cout << "\te.g. ./scallop -cyclic a5b0 aabaaaB\n";
    std::cout << "\tIf the gen string is omitted, the group is assumed to be free\n";
    std::cout << "\t-h: print this message\n";
    std::cout << "\t-v[n]: verbose output (n=0,1,2,3); 0 gives quiet output\n";
    std::cout << "\t-P: output the polygon in CDD file format\n";
    std::cout << "\t-m<format>: use the LP solver specified (EXLP uses GMP for exact output)\n";
    exit(0);
  }
  
  while (argv[current_arg][0] == '-') {
    if (argv[current_arg][1] == 'm') {
      if (argv[current_arg][2] == 'G' && argv[current_arg][3] == 'L') {
        solver = GLPK;
      } else if (argv[current_arg][2] == 'G' && argv[current_arg][3] == 'I') {
        solver = GLPK_IPT;
      } else if (argv[current_arg][2] == 'G') {
        solver = GUROBI;
      } else if (argv[current_arg][2] == 'E') {
        solver = EXLP;
      }
      
    } else if (argv[current_arg][1] == 'v') {
      if (argv[current_arg][2] == '\0') {
        verbose = 2;
      } else {
        verbose = atoi(&argv[current_arg][2]);
      }
      
    } else if (argv[current_arg][1] == 'P') {
      output_polyhedron = true;
    }
    current_arg++;
  }
  
  output_filename = std::string(argv[current_arg]);
  current_arg++;
  
  //if the first argument is a group string, then good
  //otherwise, assume it's a free group
  std::string first_arg = std::string(argv[current_arg]);
  std::string G_in = "";
  if (first_arg.size() < 2 || isalpha(first_arg[1])) {
    //it's not a group string
    int r = chain_rank(argc-current_arg, &argv[current_arg]);
    for (int i=0; i<r; ++i) {
      G_in += (char)(97+i);
      G_in += "0";
    }
  } else {  
    G_in = std::string(argv[current_arg]);                                                        //create the group
    current_arg++;
  }
  SCABBLE::CyclicProduct G(G_in); 
  
  if (verbose > 1) {
    std::cout << "Group: " << G <<"\n";
  }
  
  //record where the commas are and stuff
  std::vector<std::string> words(0);
  std::vector<std::pair<int,int> > chain_locs(0); //it's (start, len)
  int word_start = 0;
  for (int i=current_arg; i<argc; ++i) {
    std::string w = std::string(argv[i]);
    if (w == ",") {
      chain_locs.push_back(std::make_pair(word_start, words.size()-word_start));
      word_start = words.size();
      continue;
    }
    words.push_back(w);
  }
  //push on the last chain data
  chain_locs.push_back(std::make_pair(word_start, words.size()-word_start));
  
  if (verbose > 1) {
    for (int i=0; i<(int)chain_locs.size(); ++i) {
      std::cout << "Chain " << i << ": \n";
      for (int j=chain_locs[i].first; j<chain_locs[i].first + chain_locs[i].second; ++j) {
        std::cout << words[j] << " ";
      }
      std::cout << "\n";
    }
  }
  
  std::vector<SCABBLE::Pt> verts(0);
  std::vector<std::vector<SCABBLE::Pt> > faces(0);
  
  if (verbose > 1) {
    std::cout << "Computing ball\n";
  }    
  
  if ((int)chain_locs.size() > 2) {
    std::cout << "Sorry, only 2-d balls are currently supported\n";
  }
  
  SCABBLE::compute_ball(G, 
                        words, 
                        chain_locs, 
                        verts, 
                        faces, 
                        solver, 
                        verbose);
  
  if (verbose > 1) {
    std::cout << "Writing ball output\n";
  }
  
  if (output_polyhedron) {
    SCABBLE::write_ball_to_file(output_filename, verts, faces, verbose);
  } else {
    SCABBLE::draw_ball_to_file(output_filename, verts, faces, verbose);
    std::cout << "Drew ball to file\n";
  }
  
  return 0;
}

  
 
  

  
  
  