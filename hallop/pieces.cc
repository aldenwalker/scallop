#include <vector>
#include <iostream>

#include "free_group_chain.h"
#include "pieces.h" 

void HALLOP::extract_signed_index(SignedInd si, int& ind, int& s) {
  if (si < 0) {
    ind = (-si)-1;
    s = -1;
  } else {
    ind = si-1;
    s = 1;
  }
}

HALLOP::Pieces::Pieces(HALLOP::FreeGroupChain& C) {
  
  //construct the GluingEdges; this is all pairs of letters, as long as 
  //in1 < in2
  int NL = C.num_letters();
  edges.resize((NL*(NL-1))/2); //this is how many edges
  edge_index_from_pair.clear();
  GluingEdge temp_edge;
  int num_edges_recorded = 0;
  for (int i=0; i<NL; ++i) {
    temp_edge.in1 = i;
    for (int j=i+1; j<NL; ++j) {
      temp_edge.in2 = j;
      edges[num_edges_recorded] = temp_edge;
      edge_index_from_pair[ std::make_pair(i,j) ] = num_edges_recorded;
      num_edges_recorded++;
    }
  }
  
  //construct the rectangles; this is all pairs of letters 
  //which are a matching letter-inverse pair
  rects.resize(0);
  rects_containing_edge = std::vector<std::vector<int> >( edges.size(), std::vector<int>(0) );
  Rectangle temp_rect;
  for (int i=0; i<C.rank; ++i) {
    for (int j=0; j<(int)C.regular_letters[i].size(); ++j) {
      temp_rect.let1 = C.regular_letters[i][j];
      int ind,s;
      for (int k=0; k<(int)C.inverse_letters[i].size(); ++k) {
        temp_rect.let2 = C.inverse_letters[i][k];
        temp_rect.b1 = GE_index(temp_rect.let1, C.prev_letter(temp_rect.let2));
        temp_rect.b2 = GE_index(temp_rect.let2, C.prev_letter(temp_rect.let1));
        extract_signed_index(temp_rect.b1, ind, s);
        rects_containing_edge[ind].push_back(s*(rects.size()+1));
        extract_signed_index(temp_rect.b2, ind, s);
        rects_containing_edge[ind].push_back(s*(rects.size()+1));
        rects.push_back(temp_rect);
      }
    }
  }
  
  //construct the triangles; this is all triples of letters
  //no pair of letters can be the same
  //we may assume that the first letter has the smallest index
  tris.resize(0);
  tris_containing_edge = std::vector<std::vector<int> >( edges.size(), std::vector<int>(0) );
  Triangle temp_tri;
  for (int i=0; i<NL; ++i) {
    temp_tri.let1 = i;
    for (int j=i+1; j<NL; ++j) {
      temp_tri.let2 = j;
      temp_tri.b1 = GE_index(i,j);
      int ind,s;
      for (int k=i+1; k<NL; ++k) {
        if (j==k) continue;
        temp_tri.let3 = k;
        temp_tri.b2 = GE_index(j,k);
        temp_tri.b3 = GE_index(k,i);
        extract_signed_index(temp_tri.b1, ind, s);
        tris_containing_edge[ind].push_back(s*(tris.size()+1));
        extract_signed_index(temp_tri.b2, ind, s);
        tris_containing_edge[ind].push_back(s*(tris.size()+1));
        extract_signed_index(temp_tri.b3, ind, s);
        tris_containing_edge[ind].push_back(s*(tris.size()+1));
        tris.push_back(temp_tri);
      }
    }
  } 
}


HALLOP::SignedInd HALLOP::Pieces::GE_index(int let1, int let2) {
  int L1, L2, sign;
  if (let1 < let2) {
    L1 = let1; L2 = let2; sign = 1;
  } else {
    L1 = let2; L2 = let1; sign = -1;
  }
  
  return sign * (edge_index_from_pair[ std::make_pair(L1, L2) ] + 1);

}


void HALLOP::Pieces::print(std::ostream& os) {
  os << "Pieces:\n";
  os << "Gluing Edges (" << edges.size() << "): \n";
  for (int i=0; i<(int)edges.size(); ++i) {
    os << i << ": " << edges[i] << "\n";
  }  
  os << "Rectangles (" << rects.size() << "): \n";
  for (int i=0; i<(int)rects.size(); ++i) {
    os << i << ": " << rects[i] << "\n";
  }  
  os << "Triangles (" << tris.size() << "): \n";
  for (int i=0; i<(int)tris.size(); ++i) {
    os << i << ": " << tris[i] << "\n";
  }
  
}


std::ostream& HALLOP::operator<<(std::ostream& os, HALLOP::GluingEdge& e) {
  return os << "GE(" << e.in1 << ", " << e.in2 << ")";
}

std::ostream& HALLOP::operator<<(std::ostream& os, HALLOP::Rectangle& r) {
  return os << "R(" << r.let1 << ", " << r.let2 << "); bd(" << r.b1 << ", " << r.b2 << ")";
}

std::ostream& HALLOP::operator<<(std::ostream& os, HALLOP::Triangle& t) {
  return os << "T(" << t.let1 << ", " << t.let2 << ", " << t.let3 << "); bd(" << t.b1 << ", " << t.b2 << ", " << t.b3 << ")";
}