/****************************************************************************
* compute scl in free products of cyclic groups (finite and infinite)
* By Alden Walker
* Implements a generalization of the algorithm described in "scl"
*****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <ctype.h>

#include "scylla_classes.h"
#include "rational.h"
#include "scylla_lp.h"










/*****************************************************************************
 * Make the list of group polygons and rectangles. 
 * ***************************************************************************/
void compute_group_polygons_and_rectangles(Chain &C, 
                                           InterfaceEdgeList &IEL,
                                           std::vector<GroupEdgeList> &GEL,
                                           std::vector<std::vector<GroupPolygon> > &GP,
                                           std::vector<std::vector<GroupRectangle> > &GR) {
  int i;
  CyclicProduct* G = C.group();
  int num_groups = (*G).num_groups();
  std::vector<std::vector<int> > group_letters = C.group_letter_list();
  std::vector<int> orders = (*G).order_list();
  Multiset letter_selection;
  std::vector<int> regular_letters;
  std::vecctor<int> inverse_letters;
  std::vector<Multiarc> regular_multiarcs;
  std::vector<Multiarc> inverse_multiarcs;
  Multiarc temp_marc;
  GroupRectangle temp_group_rect;
  GroupPolygon temp_group_poly;
  
  arcs.resize(num_groups);
  GP.resize(num_groups);
  GR.resize(num_groups);
  
  for (i=0; i<num_groups; i++) {
    regular_letters.resize(0);
    inverse_letters.resize(0);
    for (j=0; j<(int)group_letters[i].size(); j++) {
      if ( isupper(chain_letters[group_letters[i][j]].letter) ) {
        inverse_letters.push_back(group_letters[i][j]);
      } else {
        regular_letters.push_back(group_letters[i][j]);
      }
    }
    
    //get the regulars (non-inverse)
    letter_selection = Multiset(orders[i]-1, 0, regular_letters.size());
    regular_multiarcs.resize(0);
    do {
      for (j=0; j<orders[i]-1; j++) {
        temp_marc.letters[j] = regular_letters[letter_selection[j]];
      }
      regular_multiarcs.push_back(temp_marc);
    } while (1 != letter_selection.next());
    
    //get the inverse multiarcs
    letter_selection = Multiset(orders[i]-1, 0, inverse_letters.size());
    inverse_multiarcs.resize(0);
    do {
      for (j=0; j<orders[i]-1; j++) {
        temp_marc.letters[j] = inverse_letters[letter_selection[j]];
      }
      inverse_multiarcs.push_back(temp_marc);
    } while (1 != letter_selection.next());    
    
    //now assemble the group rectangles and group polygons
    //first, let's do the rectangles
    temp_group_rect.group = i;
    temp_group_rect.edges.resize(2);
    GR[i].resize(0);
    for (j=0; j<(int)regular_letters.size(); j++) {
      for (k=0; k<(int)inverse_letters.size(); k++) {
        temp_group_rect.edges[0] = IEL.get_index_from_group_side(regular_letters[j], inverse_letters[k]);
        temp_group_rect.edges[1] = IEL.get_index_from_group_side(inverse_letters[k], regular_letters[j]);
        GR[i].push_back(temp_group_rect);
      }
    }

    //now we do the group rectangles 
    //these either have one or two multiarc sides.  first, the ones with a 
    //single multiarc side
    GP[i].resize(0);
    temp_group_poly.group = i;
    temp_group_poly.sides.resize(1);
    temp_group_poly.edges.resize(1);
    for (j=0; j<(int)regular_multiarcs.size(); j++) {
      temp_group_poly.sides[0] = regular_multiarcs[j];
      for (k=0; k<(int)GEL.regular_edges.size(); k++) {
        temp_group_poly.edges[0] = GEL[i].regular_edges[k];
        GP[i].push_back(temp_group_poly);
      }
    }
    for (j=0; j<(int)inverse_multiarcs.size(); j++) {
      temp_group_poly.sides[0] = inverse_multiarcs[j];
      for (k=0; k<(int)GEL.inverse_edges.size(); k++) {
        temp_group_poly.edges[0] = GEL[i].inverse_edges[k];
        GP[i].push_back(temp_group_poly);
      }
    }
    
    //now the ones with two multiarc sides
    temp_group_poly.sides.resize(2);
    temp_group_poly.edges.resize(2);
    for (j=0; j<(int)regular_multiarcs.size(); j++) { //first multiarc
      temp_group_poly.sides[0] = regular_multiarcs[j];
      for (k=0; k<(int)regular_multiarcs.size(); k++) { //second multiarc
        temp_group_poly.sides[0] = regular_multiarcs[k];
        for (m=0; m<(int)GEL.regular_edges.size(); m++) { //first side (between multiarc 0 and 1)
          temp_group_poly.edges[0] = GEL[i].regular_edges[m];
          for (n=0; n<(int)GEL.regular_edges.size(); n++) {
            temp_group_poly.edges[1] = GEL[i].regular_edges[n];
            GP[i].push_back(temp_group_poly);
          }
        }
      }
    }
    for (j=0; j<(int)inverse_multiarcs.size(); j++) { //first multiarc
      temp_group_poly.sides[0] = inverse_multiarcs[j];
      for (k=0; k<(int)inverse_multiarcs.size(); k++) { //second multiarc
        temp_group_poly.sides[0] = inverse_multiarcs[k];
        for (m=0; m<(int)GEL.inverse_edges.size(); m++) { //first side (between multiarc 0 and 1)
          temp_group_poly.edges[0] = GEL[i].inverse_edges[m];
          for (n=0; n<(int)GEL.inverse_edges.size(); n++) {
            temp_group_poly.edges[1] = GEL[i].inverse_edges[n];
            GP[i].push_back(temp_group_poly);
          }
        }
      }
    }
    
  }  
  
}




/*****************************************************************************
 * compute the list of polys. 
 *****************************************************************************/
void compute_polys(Chain &C, 
                   InterfaceEdgeList &IEL, 
                   CentralEdgeList &CEL,
                   std::vector<CentralPolygon> &CP) {
  int i,j,k;
  
  CentralPolygon temp_central_poly;
  
  CP.resize(0);
  temp_central_poly.edges.resize(4);
  temp_central_poly.interface.resize(4);
  
  //first, we go through everything with two central edges
  //for this, all we have to do is enumerate everything with 
  //2 interface edges, and that's it
  for (i=0; i<4; i++) {
    temp_central_poly.interface[i] = (i%2 == 0 ? true : false);
  }
  for (i=0; i<(int)IEL.edges.size(); i++) {
    temp_central_poly.edges[0] = i;
    for (j=0; j<(int)IEL.edges.size(); j++) {
      temp_central_poly.edges[2] = j;
      temp_central_poly.edges[1] = CEL.get_index( IEL[i].last, IEL[j].first );
      temp_central_poly.edges[3] = CEL.get_index( IEL[j].last, IEL[i].first );
      CP.push_back(temp_central_poly);
    }
  }
  
  //now, we go through all guys with 1 central edge
  //build chains of interface edges; whenever the length is 3, we tack on 
  //a central edge, and we're done
  std::vector<int> current_beginning_letters(3);    //this records the first letters of the interface edge choices
  std::vector<int> current_edges(3);                //this records where we are in the lists real_edges_beginning_with
  int current_len;                                  //this records the current length
  int temp_CE_1, temp_CE_2, temp_CE_3, temp_letter;
  current_len = 1;
  current_beginning_letters[0] = 0;
  current_edges[0] = 0;
  temp_central_poly.interface[0] = true;
  temp_central_poly.interface[1] = true;
  temp_central_poly.interface[2] = true;
  temp_central_poly.interface[3] = false;
  while (true) {   
    if (current_len == 3) {
      temp_CE_1 = IEL.edges_beginning_with[current_beginning_letters[current_len-1]]
                                          [current_edges[current_len-1]];
      temp_CE_2 = IEL.edges_beginning_with[current_beginning_letters[0]]
                                          [current_edges[0]];
      for (i=0; i<current_len; i++) {
        temp_central_poly.edges[i] = IEL.edges_beginning_with[current_beginning_letters[i]]
                                                             [current_edges[i]];
      }
      temp_central_poly.edges[3] = CEL.get_index(IEL[temp_CE_1].last, IEL[temp_CE_2].first);
      CP.push_back(temp_central_poly);
    }
    //now we advance it: if the list is shorter than maxmal (3), then add
    //one on.  otherwise, step back and leave it short
    if (current_len < 3) { 
      temp_CE_1 = IEL.edges_beginning_with[current_beginning_letters[current_len-1]]
                                          [current_edges[current_len-1]];
      temp_letter = C.next_letter( IEL[temp_CE_1].last );
      current_beginning_letters[current_len] = temp_letter;
      current_edges[current_len] = 0;
      current_len++;
      continue;
    }
    //if we're here, then the length is 3, so advance to the next one
    i = current_len-1;
    while (i>=0 && current_edges[i] == (int)IEL.edges_beginning_with[current_beginning_letters[i]].size()) {
      i--;
    }
    if (i==-1) {
      if (current_beginning_letters[0] == (int)chain_letters.size()-1) {
        break;
      } else {
        current_beginning_letters[0]++;
        current_edges[0] = 0;
        current_len = 1;
        continue;
      }
    } 
    current_edges[i]++;
    current_len = i+1;
  }
                 
  
  //now we search for all interface-edge polygons
  current_len = 1;
  current_beginning_letters[0] = 0;
  current_edges[0] = 0;
  temp_central_poly.interface[0] = true;
  temp_central_poly.interface[1] = true;
  temp_central_poly.interface[2] = true;
  temp_central_poly.interface[3] = true;
  while (true) {   
    if (current_len == 4) {
      temp_CE_1 = IEL.edges_beginning_with[current_beginning_letters[current_len-1]]
                                          [current_edges[current_len-1]];
      temp_CE_2 = IEL.edges_beginning_with[current_beginning_letters[0]]
                                          [current_edges[0]];
      temp_CE_3 = IEL.get_index(C.next_letter( IEL[temp_CE_1].last ),
                                C.prev_letter( IEL[temp_CE_2].first ) );
      if (temp_CE_3 != -1) {
        for (i=0; i<current_len; i++) {
          temp_central_poly.edges[i] = IEL.edges_beginning_with[current_beginning_letters[i]]
                                                               [current_edges[i]];
        }
        temp_central_poly.edges[3] = temp_CE_3;
        CP.push_back(temp_central_poly);
      }
    }
    //now we advance it: if the list is shorter than maxmal (4), then add
    //one on.  otherwise, step back and leave it short
    if (current_len < 4) { 
      temp_CE_1 = IEL.edges_beginning_with[current_beginning_letters[current_len-1]]
                                          [current_edges[current_len-1]];
      temp_letter = C.next_letter( IEL[temp_CE_1].last );
      current_beginning_letters[current_len] = temp_letter;
      current_edges[current_len] = 0;
      current_len++;
      continue;
    }
    //if we're here, then the length is 3, so advance to the next one
    i = current_len-1;
    while (i>=0 && current_edges[i] == (int)IEL.edges_beginning_with[current_beginning_letters[i]].size()) {
      i--;
    }
    if (i==-1) {
      if (current_beginning_letters[0] == (int)chain_letters.size()-1) {
        break;
      } else {
        current_beginning_letters[0]++;
        current_edges[0] = 0;
        current_len = 1;
        continue;
      }
    } 
    current_edges[i]++;
    current_len = i+1;
  }
  
}








void print_multiarcs(std::ostream &os, std::vector<std::vector<Multiarc> > &arcs) {
  int i,j,k;
  os << "Multiarcs:\n";
  for (i=0; i<(int)arcs.size(); i++) {
    os << "Group " << i << ":\n";
    for (j=0; j<(int)arcs[i].size(); j++) {
      os << "(";
      for (k=0; k<(int)arcs[i][j].letters.size(); k++) {
        os << arcs[i][j].letters[k] << ",";
      }
      os << ")\n";
    }
  }
}

void print_edges(std::ostream &os, std::vector<Edge> &edges) {
  int i;
  os << "Edges: \n";
  for (i=0; i<(int)edges.size(); i++) {
    os << i << ": (" << edges[i].first << "," << edges[i].last << ")";
    if (edges[i].blank) {
      os << "b";
    }
    os << "\n";
  }
}

void print_polys(std::ostream &os, std::vector<Edge> &edges, std::vector<Polygon> &polys) {
  int i,j;
  os << "Polygons: \n";
  for (i=0; i<(int)polys.size(); i++) {
    os << "(";
    for (j=0; j<(int)polys[i].edges.size(); j++) {
      os << polys[i].edges[j];
      if (edges[polys[i].edges[j]].blank) {
        os << "!";
      }
      os << ",";
    }
    os << ")\n";
  }
}

int main(int argc, char* argv[]) {
  int current_arg = 1;
  if (argc < 3 || std::string(argv[1]) == "-h") {
    std::cout << "usage: ./scyllop [-h] <gen string> <chain>\n";
    std::cout << "\twhere <gen string> is of the form <gen1><order1><gen2><order2>...\n";
    std::cout << "\te.g. a5b0 computes in Z/5Z * Z\n";
    std::cout << "\tand <chain> is an integer linear combination of words in the generators\n";
    std::cout << "\te.g. ./scyllop a5b0 aabaaaB\n";
    std::cout << "\t-h: print this message\n";
    exit(0);
  }
  while (argv[current_arg][0] == '-') {
    current_arg++;
    //handle arguments (none yet)
  }
  
  std::string G_in = std::string(argv[current_arg]);
  CyclicProduct G(G_in);                                                         //create the group
  current_arg++;
  
  Chain C(&G, &argv[current_arg], argc-current_arg);                              //process the chain arguments
  
  std::vector<std::vector<Multiarc> > arcs(0);
  std::vector<Edge> edges(0);
  std::vector<Polygon> polys(0);

  std::cout << "Group: " << G << "\n";
  std::cout << "Chain: " << C << "\n";
  C.print_chunks(std::cout);
  
  std::cout << "Letters:\n";
  C.print_letters(std::cout);
  
  std::cout << "Group letters:\n";
  C.print_group_letters(std::cout);
  
  compute_multiarcs(G, C, arcs);                                                 //calls arcs and polys by reference
  print_multiarcs(std::cout, arcs);
  compute_edges(G, C, edges);
  print_edges(std::cout, edges);
  compute_polys(G, C, edges, polys);
  print_polys(std::cout, edges, polys);
  
  rational scl;
  std::vector<rational> solution_vector(0);                           //run the LP
  scyllop_lp(G, C, arcs, edges, polys, &scl, &solution_vector, GLPK_DOUBLE, true); 
  
  std::cout << "scl( " << C << ") = " << scl << " = " << scl.get_d() << "\n";    //output the answer
  
  return 0;
}
