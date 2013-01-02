#ifndef gallop_H
#define gallop_H


#include <vector>
#include <string>
#include <iostream>
#include "../rational.h"
#include "graph.h"

namespace GALLOP {

struct ChainLetter {
  std::pair<int,int> position;  //the word index, then the index in the word (zero based) 
  bool marked;                  //whether the chain letter is marked (in a.b, b is marked)
  int in_tag;                   //0 means not in a tag, 1 means in one
  int edge_index;               //the graph index of the edge, -1 for opposite direction, and indices start at 1
  bool in_delta_minus;          //is it in delta minus?
};

struct Chain {
  std::vector<std::string> words;
  std::vector<int> weights;
  std::vector<ChainLetter> CL;
  std::vector<std::vector<int> > CLs_from_gen;
  Graph* G;
  
  int next_letter(int a);  
  int prev_letter(int a);
  
  Chain(Graph& G, 
        char** input, 
        int num_inputs,
        int require_f_folded,
        int verbose);
};


struct Rect {
  int first;  //the chainletters
  int second;
  int graph_edge;
  int source_edge_idx;  //where in the vertex edge list it is
  int dest_edge_idx;
  bool is_tag;
  bool in_delta_minus;
  
  friend std::ostream& operator<<(std::ostream& os, Rect R);
};


struct RectList {
  std::vector<Rect> r;
  std::vector<std::vector<int> > rects_starting_with;
  
  RectList(Chain& C, int require_f_folded, int verbose);
  RectList();
  int find_index_from_pair(int a, int b);
};


struct Poly {
  std::vector<int> rects;
  
  friend std::ostream& operator<<(std::ostream& os, Poly P);
};


Graph read_branched_surface(RectList& RL, 
                            std::vector<Poly>& P, 
                            std::vector<Rational>& solution,
                            int verbose);


int vector_clear_dens(std::vector<Rational>& a);

void gallop(int argc, char** argv);

}
#endif