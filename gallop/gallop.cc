#include <vector>
#include <iostream>
#include <string>
#include <cctype>
#include <cstdlib>
#include <cstdio>

#include "gallop.h"
#include "graph.h"
#include "../rational.h"
#include "../lp.h"
#include "../word.h"


using namespace GALLOP;

namespace GALLOP {

void extend_vector(std::vector<int>& a, std::vector<int>& b) {
  int old_size = a.size();
  int i;
  a.resize(a.size() + b.size());
  for (i=0; i<(int)b.size(); i++) {
    a[old_size + i] = b[i];
  }
}


int vector_clear_dens(std::vector<Rational>& a) {
  int L = 1;
  int i;
  for (i=0; i<(int)a.size(); i++) {
    L = lcm(L, a[i].d());
  }
  for (i=0; i<(int)a.size(); i++) {
    a[i] = a[i]*L;
  }
  return L;  
}
}


//this data structure is supposed to be a transparent way of 
//doing folded, f-folded, etc
namespace GALLOP {
struct local_condition_test {
  std::vector<bool> used_edges; //record which edges incident to a vertex are used (for foldedness)
  bool has_marked_point; //record if we've seen a marked point (for f-folded)
  bool has_delta_minus;  //record if we've seen a point on delta minus
  bool folded, f_folded;
  
  local_condition_test(bool f, bool ff);
  void new_vertex(int num_incident_edges);
  bool test_new_rect(Chain& C, RectList& RL, Poly& temp_poly, int new_rect, int verbose);
  void add_rect(Chain& C, RectList& RL, Poly& temp_poly, int new_rect);
  void pop_rect(Chain& C, RectList& RL, Poly& temp_poly);
};
}

GALLOP::local_condition_test::local_condition_test(bool f, bool ff) {
  folded = f;
  f_folded = ff;
  used_edges.resize(0);
  has_marked_point = false;
  has_delta_minus = false;
}

void GALLOP::local_condition_test::new_vertex(int num_incident_edges) {
  int j;
  if (folded) {
    used_edges.resize(num_incident_edges);
    for (j=0; j<num_incident_edges; j++) {
      used_edges[j] = false;
    }
  }
  if (f_folded) {
    has_marked_point = false;
    has_delta_minus = false;
  }
}

bool GALLOP::local_condition_test::test_new_rect(Chain& C, RectList& RL, Poly& temp_poly, int new_rect, int verbose) {
  int ind, sign;
  int edge_loc;
  bool fok=true;
  bool ffok=true;
  int j;
  int potential_chainletter;
  
  if (folded) {
    extract_signed_index(new_rect, &ind, &sign);
    edge_loc = (sign < 0 ? RL.r[ind].dest_edge_idx : RL.r[ind].source_edge_idx);
    fok = !used_edges[edge_loc];
    if (verbose > 3) {
      std::cout << "Used location vector: ";
      for (j=0; j<(int)used_edges.size(); j++) {
        std::cout << used_edges[j] << " ";
      }
      std::cout << "\n";
    }
  }
  
  if (f_folded) {
    if (has_marked_point && temp_poly.rects.size() > 1) {
      ffok = false;
    } else {
      extract_signed_index(new_rect, &ind, &sign);
      potential_chainletter = (sign < 0 ? RL.r[ind].second : RL.r[ind].first );
      if (C.CL[potential_chainletter].marked) {
        if (has_marked_point || temp_poly.rects.size() > 1) ffok = false;
      }
      if (C.CL[potential_chainletter].in_delta_minus) {
        if (has_delta_minus) {
          ffok = false;
        }
      }
    }
  }
  
  return fok && ffok;
}


void GALLOP::local_condition_test::add_rect(Chain& C, RectList& RL, Poly& temp_poly, int new_rect) {
  int ind, sign;
  int edge_loc;
  int new_chainletter;
  extract_signed_index(new_rect, &ind, &sign);
  if (folded) {
    edge_loc = (sign < 0 ? RL.r[ind].dest_edge_idx : RL.r[ind].source_edge_idx);
    used_edges[edge_loc] = true;
  }
  if (f_folded) {
    new_chainletter = (sign < 0 ? RL.r[ind].second : RL.r[ind].first);
    if (C.CL[new_chainletter].marked) {
      has_marked_point = true;
    }
    if (C.CL[new_chainletter].in_delta_minus) {
      has_delta_minus = true;
    }
  }
}

void GALLOP::local_condition_test::pop_rect(Chain& C, RectList& RL, Poly& temp_poly) {
  int ind, sign;
  int edge_loc;
  int old_chainletter;
  extract_signed_index(temp_poly.rects.back(), &ind, &sign);
  if (folded) {
    edge_loc = (sign < 0 ? RL.r[ind].dest_edge_idx : RL.r[ind].source_edge_idx);
    used_edges[edge_loc] = false;
  }
  if (f_folded) {
    old_chainletter = (sign < 0 ? RL.r[ind].second : RL.r[ind].first);
    if (C.CL[old_chainletter].marked) {
      has_marked_point = false;
    }
    if (C.CL[old_chainletter].in_delta_minus) {
      has_delta_minus = false;
    }
  }
}







void compute_polys(std::vector<Poly>& P, 
                   RectList& RL, 
                   Chain& C, 
                   bool require_folded,
                   int require_f_folded,
                   int max_polygon_valence,
                   int verbose) {
  P.resize(0);
  Poly temp_poly;
  Graph* G = C.G;
  std::vector<std::vector<int> > rect_choices(1);
  int gen;
  int ind;
  int sign;
  int current_rect; //the current polygon length
  int first_letter;
  int next_letter, second_letter;
  int i,j,k,m,n;
  int max_rects;
  int valence;
  bool ffolded = (bool)require_f_folded;
  bool do_tester = require_folded || ffolded;
  
  local_condition_test LCT(require_folded, ffolded);
  
  for (i=0; i<G->num_verts; i++) {
    
    if (verbose > 3) std::cout << "Running vertex " << i << "\n";
    
    valence = G->verts[i].incident_edges.size();
    
    if (max_polygon_valence > 0) {
      max_rects = max_polygon_valence;
    } else if (require_f_folded) {
      max_rects = valence;
    } else {
      max_rects = valence;
    }
    
    //initialize the local tester
    if (do_tester) LCT.new_vertex(G->verts[i].incident_edges.size());
  
    //build all polys at this vertex
    //we may assume that the signed index of the starting rectangle is
    //the smallest
    rect_choices.resize(1);
    rect_choices[0].resize(0);
    for (j=0; j<(int)G->verts[i].incident_edges.size(); j++) {
      gen = G->verts[i].incident_edges[j];
      sign = (G->verts[i].is_outgoing[j] ? 0 : 1);
      for (k=0; k<(int)C.CLs_from_gen[2*gen+sign].size(); k++) {
        extend_vector(rect_choices[0], 
                      RL.rects_starting_with[C.CLs_from_gen[2*gen+sign][k]] );
      }
    }
    
    current_rect = 0;
    temp_poly.rects.resize(0);
    while (rect_choices[0].size() != 0 || rect_choices.size() > 1) {
      
      if (verbose > 3) {
        std::cout << "Current rect: " << current_rect << "\n";
        std::cout << "Current polygon stack:\n";
        for (m=0; m<(int)rect_choices.size(); m++) {
          std::cout << m << ": ";
          for (n=0; n<(int)rect_choices[m].size(); n++) {
            std::cout << rect_choices[m][n] << " ";
          }
          std::cout << "\n";
        }
        std::cout << "Current temp_poly: " << temp_poly << "\n";
      }      
      
      //are we dealing with an empty choice?  If so, back up
      if (rect_choices[current_rect].size() == 0) {
        current_rect--;
        rect_choices.pop_back();
        if (do_tester) LCT.pop_rect(C, RL, temp_poly);      //we're about to back up -- if require_folded, un-check this edge
        temp_poly.rects.pop_back();
        if (verbose>3) std::cout << "backing up\n";
        continue;
      } else {
        //if the potential choice has a smaller index, forget it
        if (current_rect > 0 && rect_choices[current_rect].back() < temp_poly.rects[0]) {
          rect_choices[current_rect].pop_back();
          if (verbose>3) std::cout << "next rect too small\n";
          continue;
        }
        //do the local test
        if (do_tester) {
          if (LCT.test_new_rect(C, RL, temp_poly, rect_choices[current_rect].back(), verbose)) {
            LCT.add_rect(C, RL, temp_poly, rect_choices[current_rect].back());
          } else {
            if (verbose > 3) std::cout << "next rect would violate tester\n";
            rect_choices[current_rect].pop_back();
            continue;
          }
        } 
        temp_poly.rects.push_back( rect_choices[current_rect].back() );
        rect_choices[current_rect].pop_back();
      }
      //can we add a rectangle?
      extract_signed_index(temp_poly.rects[current_rect], &ind, &sign);
      second_letter = ( sign < 0 ? RL.r[ind].first : RL.r[ind].second);
      next_letter = C.next_letter(second_letter);
      //if this rectangle is the 0th one, then we need to reset the first_letter
      if (current_rect == 0) {
        first_letter = (sign < 0 ? RL.r[ind].second : RL.r[ind].first);
      }
      if (next_letter == first_letter) {               //the polygon can close up
        //add this polygon
        //but only if the valence is correct
        if (verbose > 3) std::cout << "This one is good: " << temp_poly << "\n";
        //if (require_trivalent) {
        //  if (temp_poly.rects.size() == 3 || valence < 3) {
        //    P.push_back(temp_poly);
        //  }
        //} else {
          P.push_back(temp_poly);
        //}
        if (do_tester) LCT.pop_rect(C, RL, temp_poly); //we're about to back up -- if require_folded, un-check this edge
        temp_poly.rects.pop_back();
        continue;
      } 
      //the polygon doesn't close up, so add on the next set of rectangle options
      //if the length isn't too big
      if (current_rect+1 == max_rects) {        
        //we're about to back up -- if require_folded, un-check this edge
        if (do_tester) LCT.pop_rect(C, RL, temp_poly);
        temp_poly.rects.pop_back();
        if (verbose>3) std::cout << "sized maxed out -- backing up\n";
        continue;
      } else {    
        rect_choices.resize(current_rect+2);
        rect_choices[current_rect+1] = RL.rects_starting_with[next_letter];
        current_rect++;
        if (verbose>3) std::cout << "added a new option\n";
        continue;
      }
    }
  }    
}


std::ostream& GALLOP::operator<<(std::ostream& os, Poly P) {
  int i;
  for (i=0; i<(int)P.rects.size(); i++) {
    std::cout << P.rects[i] << " ";
  }
  return os;
}

RectList::RectList() {
  r.resize(0);
  rects_starting_with.resize(0);
}

RectList::RectList(Chain& C,
                   int require_f_folded,
                   bool check_no_inverse_pairings,
                   int verbose) {
  int a,b;
  Rect temp_rect;
  int temp_vert;
  r.resize(0);
  rects_starting_with.resize((int)C.CL.size());
  int i,j,k,m;
  for (i=0; i<(int)C.CL.size(); i++) {
    rects_starting_with[i].resize(0);
  }
  
  for (i=0; i<(int)C.G->edges.size(); i++) {
    for (j=0; j<(int)C.CLs_from_gen[2*i].size(); j++) {
      for (k=0; k<(int)C.CLs_from_gen[2*i+1].size(); k++) {
        a = C.CLs_from_gen[2*i][j];
        b = C.CLs_from_gen[2*i+1][k];
        if (check_no_inverse_pairings && C.CL[a].inverse_letter == b) {
          continue;
        }
        if (C.CL[a].in_tag ^ C.CL[b].in_tag) {
          continue;
        }
        if (C.CL[a].in_tag && !(a == b+1 || b == a+1)) {
          continue;
        }
        if (C.CL[a].in_tag) {
          temp_rect.is_tag = true;
        } else {
          temp_rect.is_tag = false;
        }
        if (C.CL[a].in_delta_minus && C.CL[b].in_delta_minus) {
          continue;
        }
        if (C.CL[a].in_delta_minus || C.CL[b].in_delta_minus) {
          temp_rect.in_delta_minus = true;
        } else {
          temp_rect.in_delta_minus = false;
        }
        
        if (a < b) {
          rects_starting_with[a].push_back((int)r.size()+1);
          rects_starting_with[b].push_back(-((int)r.size()+1));
          temp_rect.first = a;
          temp_rect.second = b;
          temp_rect.graph_edge = i+1;
          temp_vert = C.G->edges[i].source;
          for (m=0; m<(int)C.G->verts[temp_vert].incident_edges.size(); m++) { 
            if (C.G->verts[temp_vert].incident_edges[m] == i && 
                C.G->verts[temp_vert].is_outgoing[m]) {
              break;
            }
          }
          temp_rect.source_edge_idx = m;
          temp_vert = C.G->edges[i].dest;
          for (m=0; m<(int)C.G->verts[temp_vert].incident_edges.size(); m++) { 
            if (C.G->verts[temp_vert].incident_edges[m] == i && 
                !C.G->verts[temp_vert].is_outgoing[m]) {
              break;
            }
          }          
          temp_rect.dest_edge_idx = m;
          r.push_back( temp_rect );
          
        } else {
          rects_starting_with[a].push_back(-((int)r.size()+1));
          rects_starting_with[b].push_back((int)r.size()+1);
          temp_rect.first = b;
          temp_rect.second = a;
          temp_rect.graph_edge = -(i+1);
          temp_vert = C.G->edges[i].dest;
          for (m=0; m<(int)C.G->verts[temp_vert].incident_edges.size(); m++) { 
            if (C.G->verts[temp_vert].incident_edges[m] == i && 
                !C.G->verts[temp_vert].is_outgoing[m]) {
              break;
            }
          }
          temp_rect.source_edge_idx = m;
          temp_vert = C.G->edges[i].source;
          for (m=0; m<(int)C.G->verts[temp_vert].incident_edges.size(); m++) { 
            if (C.G->verts[temp_vert].incident_edges[m] == i && 
                C.G->verts[temp_vert].is_outgoing[m]) {
              break;
            }
          }          
          temp_rect.dest_edge_idx = m;
          r.push_back( temp_rect );
        }
      }
    }
  }  
}


int RectList::find_index_from_pair(int a, int b) {
  int i,ind,sign;
  for (i=0; i<(int)rects_starting_with[a].size(); i++) {
    extract_signed_index(rects_starting_with[a][i], &ind, &sign);
    if (sign == 1 && r[ind].second == b) {
      return sign*ind;
    } else if (sign == -1 && r[ind].first == b) {
      return sign*ind;
    }
  }
  return 0; 
}
  
std::ostream& GALLOP::operator<<(std::ostream& os, Rect R) {
  os << "R(" << R.first << ", " << R.second << "); edge: " << R.graph_edge << 
        "; source_idx: " << R.source_edge_idx << "; dest_idx: " << R.dest_edge_idx;
  return os;
}

void cyc_red_marked(std::string& s) {
  while ( (s.size() > 0) && (s[0] == swapCaseChar(s[s.size()-1])) ) {
    s.erase(0,1);
    s.erase(s.size()-1, 1);
  }
  if (s[0] == '.') {
    while ( (s.size() > 0) && (s[1] == swapCaseChar(s[s.size()-1])) ) {
      s.erase(1,1);
      s.erase(s.size()-1, 1);
    }
  }
}

Chain::Chain(Graph& input_G, 
              char** input, 
              int num_inputs, 
              int require_f_folded,
              int verbose) {
  int i,j;
  int ind, sign;
  
  G = &input_G;
  
  //load in the words
  words.resize(num_inputs);
  weights.resize(num_inputs);
  for (i=0; i<num_inputs; i++) {
    words[i] = std::string(input[i]);
    int coef = 1;
    if (isdigit(words[i][0])) {
      j=0;
      while (isdigit(words[i][j])) j++;
      coef = atoi(words[i].substr(0,j).c_str());
      words[i] = words[i].substr(j, words[i].size());
    }
    cyc_red_marked(words[i]);
    weights[i] = coef;
  }
  
  //load in the letters which we actually use
  //scan for marked letters
  ChainLetter temp;
  temp.inverse_letter = -1;
  CL.resize(0);
  for (i=0; i<(int)words.size(); i++) {
    temp.in_delta_minus = (i < require_f_folded);
    j=0;
    while (j < (int)words[i].size()) {
      if (words[i][j] == '.') {
        words[i].erase(j,1);
        temp.marked = true;
      } else {
        temp.marked = false;
      }
      
      if (words[i][j] == '/') { //make a tag
        words[i].erase(j,1);
        temp.in_tag = 1;
        temp.position.first = i;
        temp.position.second = j;
        temp.edge_index = G->get_edge_from_label(words[i][j]);
        CL.push_back(temp);
        temp.position.second = j+1;
        temp.edge_index = G->get_edge_from_label(words[i][j+1]);
        temp.in_tag = 1;
        CL.push_back(temp);
        words[i].erase(j+2,1);
        j+=2;
      
      } else {      
        temp.in_tag = 0;
        temp.position.first = i;
        temp.position.second = j;
        temp.edge_index = G->get_edge_from_label(words[i][j]);
        CL.push_back(temp);
        j++;
      }
    }
  }
  
  CLs_from_gen.resize( 2*(int)G->edges.size() );
  for (i=0; i<2*(int)G->edges.size(); i++) {
    CLs_from_gen[i].resize(0);
  }
  for (i=0; i<(int)CL.size(); i++) {
    extract_signed_index( CL[i].edge_index, &ind, &sign);
    sign = (sign < 0 ? 1 : 0);
    CLs_from_gen[2*ind + sign].push_back(i);
  }
  
}


int Chain::next_letter(int a) {
  int word = CL[a].position.first;
  int pos = CL[a].position.second;
  int word_len = words[word].size();
  if (pos == word_len-1) {
    return a-(word_len-1);    
  } else {
    return a+1;
  }
}

int Chain::prev_letter(int a) {
  int word = CL[a].position.first;
  int pos = CL[a].position.second;
  int word_len = words[word].size();
  if (pos == 0) {
    return a+(word_len-1);    
  } else {
    return a-1;
  }
}


int Chain::CL_from_inds(int w, int ell) {
  std::pair<int, int> pos_to_find = std::make_pair(w, ell);
  for (int i=0; i<(int)CL.size(); ++i) {
    if (pos_to_find == CL[i].position) {
      return i;
    }
  }
  return -1;
}

void Chain::add_inverse_word(int i) {
  std::string w_inv = inverse(words[i]);
  int wL = w_inv.size();
  int new_ind = words.size();
  
  words.push_back(w_inv);
  weights.push_back(weights[i]);
  
  for (int j=0; j<wL; ++j) {
    int let_inv_ind = CL_from_inds( i, wL-1-j );
    ChainLetter temp;
    temp = CL[let_inv_ind];
    temp.position = std::make_pair( new_ind, j );
    temp.edge_index = -temp.edge_index;
    temp.inverse_letter = let_inv_ind;
    
    CL[let_inv_ind].inverse_letter = CL.size();
    int ind, sign;
    extract_signed_index(temp.edge_index, &ind, &sign);
    sign = (sign < 0 ? 1 : 0);
    CLs_from_gen[2*ind+sign].push_back(CL.size());
    CL.push_back(temp);
  } 
  
}

std::ostream& GALLOP::operator<<(std::ostream &os, Chain &C) {
  int i;
  int len = (int)C.words.size();
  for (i=0; i<len; i++) {
    if (C.weights[i] != 1) {
      os << C.weights[i] << C.words[i] << " ";
    } else {
      os << C.words[i] << " ";
    }
  }
  return os;
}





Graph read_branched_surface(Chain& C,
                             RectList& RL, 
                             std::vector<Poly>& P, 
                             std::vector<Rational>& solution,
                             int verbose) {
  Graph G;
  char temp[100];
  int i,j;
  std::vector<Poly> available_polys(0);
  
  vector_clear_dens(solution);  //we don't return this
  
  for (i=0; i<(int)solution.size(); i++) {
    for (j=0; j<solution[i].n(); j++) {
      available_polys.push_back(P[i]);
    }
  }
  G.num_verts = available_polys.size();
  G.verts.resize(G.num_verts);
  for (i=0; i<G.num_verts; i++) {
    sprintf(temp, "VERT%d", i);
    G.verts[i].name = std::string(temp);
    G.verts[i].incident_edges.resize(available_polys[i].rects.size());
    G.verts[i].is_outgoing.resize(available_polys[i].rects.size());
    for (j=0; j<(int)available_polys[i].rects.size(); j++) {
      G.verts[i].incident_edges[j] = -1;
      G.verts[i].is_outgoing[j] = (available_polys[i].rects[j] > 0);
    }
  }
  
  bool we_did_something = true;
  int first_edge_vert, first_edge_index, first_edge_rect;
  int second_edge_vert, second_edge_index, second_edge_rect;
  ChainLetter first_chainletter;
  ChainLetter second_chainletter;
  int ind1, sign1, ind2, sign2;
  int ind,sign;
  
  G.num_edges = 0;
  G.edges.resize(0);
  
  while (we_did_something) {
    we_did_something = false;
    
    //find a blank edge
    first_edge_vert = -1;
    first_edge_index = -1;
    first_edge_rect = -1;
    for (i=0; i<G.num_verts; i++) {
      for (j=0; j<(int)G.verts[i].incident_edges.size(); j++) {
        if (G.verts[i].incident_edges[j] == -1) {
          first_edge_vert = i;
          first_edge_index = j;
          first_edge_rect = available_polys[i].rects[j];
          we_did_something = true;
          break;
        }
      }
      if (first_edge_index != -1) {
        break;
      }
    }
    if (first_edge_vert == -1) {
      if (verbose > 1) std::cout << "Unable to find a first edge -- we're done\n";
      continue;
    }
    
    //now find a matching edge
    second_edge_index = -1;
    second_edge_vert = -1;
    second_edge_rect = -first_edge_rect;
    for (i=0; i<G.num_verts; i++) {
      for (j=0; j<(int)G.verts[i].incident_edges.size(); j++) {
        if (G.verts[i].incident_edges[j] == -1 &&
            available_polys[i].rects[j] == second_edge_rect) {
          second_edge_vert = i;
          second_edge_index = j;
          break;
        }
      }
      if (second_edge_vert != -1) {
        break;
      }
    }
    if (second_edge_vert == -1) {
      std::cout << "Unable to find a matching edge\n";
    }
    
    //attach a new edge
    Edge e;
    sprintf(temp, "EDGE%d", G.num_edges);
    e.name = std::string(temp);
    e.source = (first_edge_rect < 0 ? second_edge_vert : first_edge_vert);
    e.dest = (first_edge_rect < 0 ? first_edge_vert : second_edge_vert);
    extract_signed_index(first_edge_rect, &ind, &sign);
    first_chainletter = C.CL[RL.r[ind].first];
    second_chainletter = C.CL[RL.r[ind].second];
    extract_signed_index(first_chainletter.edge_index, &ind1, &sign1);
    extract_signed_index(second_chainletter.edge_index, &ind2, &sign2);
    e.label_forward = (sign1 < 0 ? C.G->edges[ind1].label_backward : C.G->edges[ind1].label_forward );
    //sprintf(temp, "%d", RL.r[ind].first);
    //e.label_forward += std::string(temp);
    e.label_backward = (sign2 < 0 ? C.G->edges[ind2].label_backward : C.G->edges[ind2].label_forward );
    //sprintf(temp, "%d", RL.r[ind].second);
    //e.label_backward += std::string(temp);
    G.edges.push_back(e);
    G.verts[first_edge_vert].incident_edges[first_edge_index] = G.num_edges;
    G.verts[second_edge_vert].incident_edges[second_edge_index] = G.num_edges;
    
    if (verbose > 2) {
      std::cout << "Joined vertices " << first_edge_vert << " and " << second_edge_vert 
      << " with edge " << G.num_edges << "\n";
    }
    
    G.num_edges ++;
  }
  
  
  return G;
}






void gallop_lp(Chain& C,
               RectList& RL, 
               std::vector<Poly>& P, 
               std::vector<Rational>& solution_vector, 
               Rational& scl, 
               bool only_check_exists,
               bool check_polygonal,
               SparseLPSolver solver,
               int time_limit,
               int verbose) {
  int i,j;
  int ind, sign;
  int letter;
  //build the sparse matrix
  //there is a row for each rectangle, plus one for each word
  //there is a column for each polygon
  int num_extra_rows = (check_polygonal ? 1 : C.words.size());
  int num_rows = RL.r.size() + num_extra_rows;
  int num_cols = P.size();
  
  std::vector<int> temp_ia(0);
  std::vector<int> temp_ja(0);
  std::vector<int> temp_ar(0);
  
  SparseLP LP(solver, num_rows, num_cols);
  
  //run through the columns
  for (i=0; i<num_cols; i++) {
    temp_ia.resize(0); temp_ja.resize(0); temp_ar.resize(0);
    LP.set_obj(i, (int)(only_check_exists ? 0 : P[i].rects.size()-2));
    //objective[i] = (only_check_exists ? 0 : ((double)P[i].rects.size())-2);
    for (j=0; j<(int)P[i].rects.size(); j++) {
      //this is the rectangle rows
      extract_signed_index(P[i].rects[j], &ind, &sign);
      temp_ia.push_back(ind);
      temp_ja.push_back(i);
      temp_ar.push_back(sign);
      //this is the rows for the words -- we count 1 whenever there is a first letter
      letter = ( sign < 0 ? RL.r[ind].second : RL.r[ind].first );
      if (C.CL[letter].position.second == 0 && 
          (!check_polygonal || C.CL[letter].position.first == 0)) {
        temp_ia.push_back( RL.r.size() + C.CL[letter].position.first );
        temp_ja.push_back( i );
        temp_ar.push_back(1);
      }        
    }
    LP.extend_entries_no_dups(temp_ia, temp_ja, temp_ar);
  }
  
  //run through the rows
  for (i=0; i<(int)RL.r.size(); i++) {
    LP.set_RHS(i,0);
    LP.set_equality_type(i, EQ);
  }
  for (i=0; i<num_extra_rows; i++) {
    LP.set_RHS(RL.r.size() + i, C.weights[i]);
    LP.set_equality_type(RL.r.size() + i, EQ);
    //RHS[RL.r.size() + i] = C.weights[i];
  }
  
  SparseLPSolveCode code = LP.solve(verbose);
  
  if (code != LP_OPTIMAL) {
    scl = Rational(-1,1);
  } else {
    LP.get_optimal_value(scl);
    LP.get_soln_vector(solution_vector);
  }
}





void GALLOP::gallop(int argc, char** argv) {
  
  int verbose = 1;
  int current_arg = 0;
  std::string input_file_name="";
  std::string output_file_name;
  Graph G;
  RectList RL;
  int i,j;
  bool require_folded = false;
  int require_f_folded = 0;
  bool only_check_exists = false;
  int max_polygon_valence = -1;
  bool do_output = false;
  bool check_polygonal = false;
  bool check_polygonal_relaxed = false;
  SparseLPSolver solver = GLPK_SIMPLEX;
  int time_limit=0;
  
  if (argc < 1 || std::string(argv[0]) == "-h") {
    std::cout << "usage: ./scallop -local [-v[n]] [-f] [-ff[n=1]] [-e] [-tn] [-pn] [-y,Y] [-m<GLPK,GIPT,GUROBI,EXLP>] [-G<graph input file>] [-o <surface (graph) output file>] <chain>\n";
    std::cout << "\t-v[n]: verbose output (level n)\n";
    std::cout << "\t-y: check if the chain is polygonal (overrides -f,-ff,-p)\n";
    std::cout << "\t-Y: check polygonal without folded restriction\n";
    std::cout << "\t-f: only search for surfaces that are folded\n";
    std::cout << "\t-ff[n=1]: only search for surfaces that are f-folded:\n";
    std::cout << "\t\tn is the number of words in delta^- (these must come first)\n";
    std::cout << "\t\tmark verts in delta^+ with periods\n";
    std::cout << "\t-e: only check the existence of a feasible surface\n";
    std::cout << "\t-tn: set an LP time limit of n seconds (only works with Gurobi)\n";
    std::cout << "\t-pn: only use polygons which have at most n sides\n";
    std::cout << "\t-m<method>: specify which LP solver to use\n";
    std::cout << "\t-o filename: write out the solution fatgraph\n";
    return;
  }
  
  while (argv[current_arg][0] == '-') {
    if (argv[current_arg][1] == 'v') {
      if (argv[current_arg][2] == '\0') {
        verbose = 2;
      } else {
        verbose = atoi(&argv[current_arg][2]);
      }
    }
    
    else if (argv[current_arg][1] == 'f') {
      if (argv[current_arg][2] == 'f') {
        if (argv[current_arg][3] == '\0') {
          require_f_folded = 1;
        } else {
          require_f_folded = atoi(&argv[current_arg][3]);
        }
      } else {
        require_folded = true;
      }
    }
    
    else if (argv[current_arg][1] == 'e') {
      only_check_exists = true;
    }
    
    else if (argv[current_arg][1] == 'o') {
      do_output = true;
      current_arg++;
      output_file_name = std::string(argv[current_arg]);
    }
    
    else if (argv[current_arg][1] == 't') {
      time_limit = atoi(&argv[current_arg][2]);
    }
    
    else if (argv[current_arg][1] == 'm') {
      if (argv[current_arg][2] == 'E') {
        solver = EXLP;
      } else {
        if (argv[current_arg][3] == 'L') {
          solver = GLPK;
        } else if (argv[current_arg][3] == 'I') {
          solver = GLPK_IPT;
        } else {
          solver = GUROBI;
        }
      }
    }
    
    else if (argv[current_arg][1] == 'p') {
      max_polygon_valence = atoi(&argv[current_arg][2]);
    }
    
    else if (argv[current_arg][1] == 'G') {
      input_file_name = std::string(&argv[current_arg][2]);
    }
    
    else if (argv[current_arg][1] == 'y') {
      check_polygonal = true;
    }
    
    else if (argv[current_arg][1] == 'Y') {
      check_polygonal_relaxed = true;
    }
    
    current_arg++;
  }
  
  if (check_polygonal) {
    require_f_folded = false;
    require_folded = true;
    max_polygon_valence = -1;
  } else if (check_polygonal_relaxed) {
    check_polygonal = true;
    require_f_folded = false;
    require_folded = false;
    max_polygon_valence = -1;
  }
    
  
  if (input_file_name != "") {
    //read in the file
    G.read_file(input_file_name, verbose);
  } else {
    G.set_labeled_rose(chain_gens(argc-current_arg, &argv[current_arg]), verbose);
  }
  
  if (verbose>1) std::cout << "Set labeled rose\n";
  
  //and the chain
  Chain C(G, &argv[current_arg], argc-current_arg, require_f_folded, verbose);
  if (check_polygonal) {
    int num_words = C.words.size();
    for (int i=0; i<(int)num_words; ++i) {
      C.add_inverse_word(i);
    }
  }
  
  if (verbose > 1) {
    if (require_folded) std::cout << "Requiring a folded surface\n";
    if (max_polygon_valence > 0) std::cout << "Requiring valence at most " << max_polygon_valence << "\n";
    std::cout << "Input chain:\n";
    for (i=0; i<(int)C.words.size(); i++) {
      std::cout << C.weights[i] << "*" << C.words[i] << " ";
    }
    std::cout << "\n";
    std::cout << "Letters for each generator:\n";
    for (i=0; i<(int)C.G->edges.size(); i++) {
      std::cout << C.G->edges[i].label_forward << ": ";
      for (j=0; j<(int)C.CLs_from_gen[2*i].size(); j++) {
        if (C.CL[C.CLs_from_gen[2*i][j]].marked) {
          std::cout << C.CLs_from_gen[2*i][j] << "* ";
        } else {
          std::cout << C.CLs_from_gen[2*i][j] << " ";
        }
      }
      std::cout << "\n";
      std::cout << C.G->edges[i].label_backward << ": ";
      for (j=0; j<(int)C.CLs_from_gen[2*i+1].size(); j++) {        
        if (C.CL[C.CLs_from_gen[2*i+1][j]].marked) {
          std::cout << C.CLs_from_gen[2*i+1][j] << "* ";
        } else {
          std::cout << C.CLs_from_gen[2*i+1][j] << " ";
        }
      }
      std::cout << "\n";
    }
  }
  
  
  //compute the polygons and rectangles
  RL = RectList(C, require_f_folded, check_polygonal, verbose);
  if (verbose > 1) {
    std::cout << "Generated " << RL.r.size() << " rectangles\n";
    if (verbose > 2) {
      for (i=0; i<(int)RL.r.size(); i++) {
        std::cout << i << ": " << RL.r[i] << "\n";
      }
    }
  }
  
  std::vector<Poly> P(0);
  compute_polys(P, RL, C, require_folded, require_f_folded, max_polygon_valence, verbose);
   if (verbose > 1) {
    std::cout << "Generated " << P.size() << " polygons\n";
    if (verbose > 2) {
      for (i=0; i<(int)P.size(); i++) {
        std::cout << i << ": ";
        for (j=0; j<(int)P[i].rects.size(); j++) {
          std::cout << P[i].rects[j] << " ";
        }
        std::cout << "\n";
      }
    }
  } 
  
  if (P.size() == 0) {
    if (verbose > 0) {
      std::cout << "No feasible solution found\n";
    }
    return;
  }
  
  //do the linear programming
  std::vector<Rational> solution_vector(0);
  Rational scl;
  gallop_lp(C, 
            RL, 
            P, 
            solution_vector, 
            scl, 
            only_check_exists, 
            check_polygonal,
            solver,
            time_limit,
            verbose);
  
  if (verbose > 0) {
    if (scl < Rational(0,1)) {
      std::cout << "No feasible solution found\n";
    } else {
      if (only_check_exists) {
        std::cout << "Feasible solution found\n";
      } else {
        if (check_polygonal) {
          std::cout << C << "is polygonal with min -chi/2n = " << scl << " = " << scl.get_d() << "\n";
        } else {
          std::cout << "scl( " << C << ") = " << scl << " = " << scl.get_d() << "\n";
        }
      }
    }
  }
  
  //build a solution graph
  if (do_output) {
    Graph surface = read_branched_surface(C, RL, P, solution_vector, verbose);
    surface.write_file(output_file_name, verbose);
  }
  
}
