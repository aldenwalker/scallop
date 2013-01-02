/*****************************************************************************
* Bad form, but all the classes are going in here, because they are small
******************************************************************************/

#ifndef trollop_classes_H
#define trollop_classes_H

#include <iostream>
#include <string>
#include <vector>
#include <utility>


namespace TROLLOP {

struct WVec;


int letter_to_number(char let, int rank, char letter_removed);
char number_to_letter(int n, int rank, char letter_removed);
  void int_to_word(std::string &S, int n, int rank);
  int word_to_int(std::string& S, int start_index, int rank);
  int sub_1_mod(int a, int m);
  int int_pow(int a, int b);
  void get_cyclic_subword(std::string &S, std::string &T, int start, int n);
    


/*****************************************************************************
 This class creates a table of the words of length ell and ell-1
 (it doesn't actually create a list of all of them).  These are the vertices 
 and edges of a train track.  It also computes unique indices for the vertices
 and edges so that if we are asked for an scl computation, we can 
 do it using a small amount of space
 ****************************************************************************/
struct WordTable {
  int rank;
  int ell;
  int first_letter_offset_edges;
  int first_letter_offset_vertices;
  int num_real_edges;
  int num_real_verts;
  int num_edges;
  int num_verts;
  bool do_sup;
  bool have_assigned_indices;
  std::vector<std::vector<int> > h_vector;
  std::vector<int> vertex_to_index;
  std::vector<int> edge_to_index;
  std::vector<int> index_to_vertex;
  std::vector<int> index_to_edge;
  
  WordTable(int R, int L, bool DO_SUP);
  void create_index_assignments(WVec &C);
  void print();
  
  int get_real_index(std::string &S);
  void get_real_word(std::string &S, int length, int index);
  char get_h_from_real_index(int index);
  int get_real_edge_dest(int edge);
  int get_real_edge_source(int edge);
  char get_real_outgoing_h(int vertex);
  
  int get_index(std::string &S);
  void get_word(std::string &S, int length, int index);
  char get_h(int index);
  int get_edge_dest(int index);
  int get_edge_source(int index);
  char get_outgoing_h(int index);
};


/*****************************************************************************
 This class gives a vector of words in W_ell (or, more precisely, 
 R[E], since it doesn't check incorrect input)
 It gives a list of real indices in the word table.  If fill_index_assignments 
 is called, it also stores the indices of the words.
 Note this is basically just a list of edges
 *****************************************************************************/
struct WVec {
  std::vector<std::string> original_input;
  std::vector<std::pair<int, int> > word_list;
  std::vector<int> index_coefficients;
  WordTable* wt;
  
  WVec( WordTable &WT, char** words, int num_words, bool USE_WORDS);
  WVec();
  void fill_index_coefficients(WordTable &WT);
  void print_words(std::ostream& os);
};

std::ostream& operator<<(std::ostream& os, WVec& C);

/****************************************************************************
 This class gives a list of pairs of arcs.  Each arc gives a pair of 
 vertices.  We have (a,b) = -(b,a), so we will only record whichever 
 is lexicographically first
 ****************************************************************************/
struct ArcPairList {
  int num_arcs;
  std::vector<std::pair<std::pair<int, int>, std::pair<int, int> > > arc_list;
  std::vector<std::vector<std::vector<int> > > arcs_beginning_with;
  
  ArcPairList(WordTable &WT, int num_copies);
  int index_of_arc(int side_a, int a, int side_b, int b);
  void print(std::ostream& os);
};


/***************************************************************************
 This class gives a triangle
 it's a list of three vertices, by index in the word table
 * the arcs go v0 --a0--> v1, etc
*****************************************************************************/
struct Triangle {
  int v0,v1,v2;
  int a0,a1,a2;
};
std::ostream& operator<<(std::ostream& os, Triangle& C);

/*****************************************************************************
 This class gives a rectangle, as given by a list of two edges e_0, e_1 and 
 arcs a_0, a_1.  They go (a_0, e_0, a_1, e_1) such that:
 -  destination(a_i) = e_1 and origin(a_i) = e_{i-1}
 -  h(e_0) = inverse(h(e_1)), where h gives the labels.  
 Recall our labels are simply the last letter of the edges
 everything is listed by index in the word table
 *****************************************************************************/
struct Rectangle {
  int side0, side1;
  int a0,e0,a1,e1;
};
std::ostream& operator<<(std::ostream& os, Rectangle& C);
  
  
void extract_signed_index(int* sign, int* index, int signed_index);
  
  
  
  
  
  
  
  
}
  
  
  
  
  

#endif 
