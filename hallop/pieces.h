#ifndef HALLOP_PIECES_H
#define HALLOP_PIECES_H

#include <vector>
#include <map>
#include <utility>
#include <iostream>

#include "free_group_chain.h"

namespace HALLOP {
  
  typedef int SignedInd;
  
  void extract_signed_index(SignedInd si, int& ind, int& s);
  
  struct GluingEdge {
    int in1; //the ChainLetter (indices) which are incoming; in1 < in2 always
    int in2;
  };

  struct Rectangle {
    int let1; //the letters on either side of the rectangle
    int let2;
    SignedInd b1; //the boundary; this is the edge AFTER let1, and so on
    SignedInd b2; 
  };

  struct Triangle {
    int let1; //the three incoming letters, counterclockwise
    int let2;
    int let3;
    SignedInd b1; //the boundary; this is the edge AFTER let1, and so on
    SignedInd b2;
    SignedInd b3;
  };
  
  struct Pieces {
    
    Pieces(FreeGroupChain& C);
    
    std::vector<GluingEdge> edges;
    std::vector<Rectangle> rects;
    std::vector<Triangle> tris;
    
    std::map<std::pair<int,int>, int> edge_index_from_pair; //the reverse-lookup map (only includes edges, not negatives)
    SignedInd GE_index(int let1, int let2);               //uses the map and returns a signed index of the GluingEdge for any pair input
    
    std::vector<std::vector<SignedInd> > rects_containing_edge; //for every edge, which rectangles contain it (positively or negatively)
    std::vector<std::vector<SignedInd> > tris_containing_edge;  //similarly for triangles
    
    void print(std::ostream& os);
  };
  
  std::ostream& operator<<(std::ostream& os, GluingEdge& e);
  std::ostream& operator<<(std::ostream& os, Rectangle& r);
  std::ostream& operator<<(std::ostream& os, Triangle& t);
  
}


#endif