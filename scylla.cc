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
void compute_group_teeth_and_rectangles(Chain &C, 
                                         InterfaceEdgeList &IEL,
                                         std::vector<GroupTooth > &GT,
                                         std::vector<GroupRectangle > &GR) {
  int i,j,k,l,m,n;
  int ord;
  int num_groups = (C.G)->num_groups();
  GroupTooth temp_group_tooth;
  GroupRectangle temp_group_rect;
  
  GT.resize(0);
  GR.resize(0);
  
  for (i=0; i<num_groups; i++) {
    ord = (C.G)->index_order(i);
    
    //compute the group rectangles
    temp_group_rect.group_index = i;
    for (j=0; j<(int)C.regular_letters[i].size(); j++) {
      temp_group_rect.first = C.regular_letters[i][j];
      for (k=0; k<(int)C.inverse_letters[i].size(); k++) {
        temp_group_rect.last = C.inverse_letters[i][j];
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
 *****************************************************************************/
void compute_central_polys(Chain &C, 
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
    //for (i=0; i<current_len; i++) {
    //  std::cout << "(" << current_beginning_letters[i] << ", " << current_edges[i] << ") ";
    //}
    //std::cout << "\n";
    if (current_len == 3) {
      temp_CE_1 = IEL.edges_beginning_with[current_beginning_letters[current_len-1]][current_edges[current_len-1]];
      temp_CE_2 = IEL.edges_beginning_with[current_beginning_letters[0]][current_edges[0]];
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
    while (i>=0 && current_edges[i] == (int)IEL.edges_beginning_with[current_beginning_letters[i]].size()-1) {
      i--;
    }
    if (i==-1) {
      if (current_beginning_letters[0] == (int)C.chain_letters.size()-1) {
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
  current_beginning_letters.resize(3);
  current_edges.resize(3);
  current_beginning_letters[0] = 0;
  current_edges[0] = 0;
  temp_central_poly.interface[0] = true;
  temp_central_poly.interface[1] = true;
  temp_central_poly.interface[2] = true;
  temp_central_poly.interface[3] = true;
  while (true) {     
    //for (i=0; i<current_len; i++) {
    //  std::cout << "(" << current_beginning_letters[i] << ", " << current_edges[i] << ") ";
    //}
    //std::cout << "\n";
    if (current_len == 3) {
      temp_CE_1 = IEL.edges_beginning_with[current_beginning_letters[current_len-1]]
                                          [current_edges[current_len-1]];
      temp_CE_2 = IEL.edges_beginning_with[current_beginning_letters[0]]
                                          [current_edges[0]];
      temp_CE_3 = IEL.get_index_from_poly_side(C.next_letter( IEL[temp_CE_1].last ),
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
    if (current_len < 3) { 
      for (k=current_edges[current_len-1];
           k<(int)IEL.edges_beginning_with[current_beginning_letters[current_len-1]].size();
           k++) {
        temp_CE_1 = IEL.edges_beginning_with[current_beginning_letters[current_len-1]][k];
        temp_letter = C.next_letter( IEL[temp_CE_1].last );
        if (temp_letter > current_beginning_letters[0]) {
          break;
        }
      }
      if (k<(int)IEL.edges_beginning_with[current_beginning_letters[current_len-1]].size()) {
        current_beginning_letters[current_len] = temp_letter;
        current_edges[current_len] = 0;
        current_edges[current_len-1] = k;
        current_len++;
        continue;
      } else {
        current_edges[current_len-1] = k-1;
        //we need to advance to the next letter; the following code will do that
        //because we just set current_edges[current_len-1] to the last one
      }
    }
    //if we're here, then the length is 3, so advance to the next one
    i = current_len-1;
    while (i>=0 && current_edges[i] == (int)IEL.edges_beginning_with[current_beginning_letters[i]].size()-1) {
      i--;
    }
    if (i==-1) {
      if (current_beginning_letters[0] == (int)C.chain_letters.size()-1) {
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



void print_central_polys(std::vector<CentralPolygon> &CP, std::ostream &os, int level) {
  int i,j;
  os << "Central polygons: (" << CP.size() << "):\n"; 
  if (level == 0) {
    os << "(" << CP.size() << " polygons hidden)\n";
  } else {
    for (i=0; i<(int)CP.size(); i++) {
      os << i << ": ";
      for (j=0; j<(int)CP[i].edges.size(); j++) {
        os << CP[i].edges[j];
        if (CP[i].interface[j]) {
          os << "i";
        } else {
          os << "c";
        }
        os << " ";
      }
      os << "\n";
    }
  }
}


void print_group_teeth_mouths_polys_and_rectangles(std::vector<std::vector<GroupTooth> > &GT,
                                                  std::vector<std::vector<GroupMouth> > &GM,
                                                  std::vector<std::vector<GroupPolygon> > &GP,
                                                  std::vector<std::vector<GroupRectangle> > &GR,
                                                  std::ostream &os,
                                                  int level) {
  int i,j;
  for (i=0; i<(int)GP.size(); i++) {
    os << "Group " << i << " teeth:\n";
    if (level == 0) {
      os << "(" << GT[i].size() << " hidden)\n";
    } else {
      for (j=0; j<(int)GT[i].size(); j++) {
        os << j << ": " << GT[i][j] << "\n";
      }
    }
    os << "Group " << i << " mouths:\n";
    if (level == 0) {
      os << "(" << GM[i].size() << " hidden)\n";
    } else {
      for (j=0; j<(int)GM[i].size(); j++) {
        os << j << ": " << GM[i][j] << "\n";
      }
    }
    os << "Group " << i << " polygons:\n";
    if (level == 0) {
      os << "(" << GP[i].size() << " polygons hidden)\n";
    } else {
      for (j=0; j<(int)GP[i].size(); j++) {
        os << j << ": " << GP[i][j] << "\n";
      }
    }    
    os << "Group " << i << " rectangles:\n";
    if (level == 0) {
      os << "(" << GR[i].size() << " rectangles hidden)\n";
    } else {
      for (j=0; j<(int)GR[i].size(); j++) {
        os << j << ": " << GR[i][j] << "\n";
      }
    }
  }
}
      




int main(int argc, char* argv[]) {
  int current_arg = 1;
  //int i;
  bool VERBOSE = false;
  bool IPT = false;
  
  if (argc < 3 || std::string(argv[1]) == "-h") {
    std::cout << "usage: ./scyllop [-h] [-v] [-i] <gen string> <chain>\n";
    std::cout << "\twhere <gen string> is of the form <gen1><order1><gen2><order2>...\n";
    std::cout << "\te.g. a5b0 computes in Z/5Z * Z\n";
    std::cout << "\tand <chain> is an integer linear combination of words in the generators\n";
    std::cout << "\te.g. ./scyllop a5b0 aabaaaB\n";
    std::cout << "\t-h: print this message\n";
    std::cout << "\t-v: verbose output\n";
    std::cout << "\t-i: use the interior point LP method (faster but rational output is sometimes wrong)\n";
    exit(0);
  }
  while (argv[current_arg][0] == '-') {
    if (argv[current_arg][1] == 'i') {
      IPT = true;
    } else if (argv[current_arg][1] == 'v') {
      VERBOSE = true;
    }
    current_arg++;
  }
  
  std::string G_in = std::string(argv[current_arg]);
  CyclicProduct G(G_in);                                                         //create the group
  current_arg++;
  
  Chain C(&G, &argv[current_arg], argc-current_arg);                              //process the chain argument

  if (VERBOSE) {
    std::cout << "Group: " << G << "\n";
    std::cout << "Chain: " << C << "\n";
    std::cout << "Letters:\n";
    C.print_letters(std::cout);
    std::cout << "Group letters:\n";
    C.print_group_letters(std::cout);
  }
  
  CentralEdgeList CEL(C);
  if (VERBOSE) CEL.print(std::cout);
  
  InterfaceEdgeList IEL(C);
  if (VERBOSE) IEL.print(std::cout);
  
  std::vector<CentralPolygon> CP;
  compute_central_polys(C, IEL, CEL, CP);
  std::cout << "computed polys (" << CP.size() << ")\n"; std::cout.flush();
  if (VERBOSE) print_central_polys(CP, std::cout, VERBOSE);
  
  std::vector<std::vector<GroupTooth> > GT;
  std::vector<std::vector<GroupMouth> > GM;
  std::vector<std::vector<GroupRectangle> > GR;
  std::vector<std::vector<GroupPolygon> > GP;
  compute_group_teeth_mouths_polygons_and_rectangles(C, IEL, GEL, GT, GM, GP, GR);
  std::cout << "computed group stuff\n"; std::cout.flush();
  if (VERBOSE) print_group_teeth_mouths_polys_and_rectangles(GT, GM, GP, GR, std::cout, VERBOSE);
  
   
  rational scl;
  std::vector<rational> solution_vector(0);                           //run the LP
  scylla_lp(C, GEL, IEL, CEL, CP, GT, GM, GP, GR, 
            &scl, 
            &solution_vector, 
            (IPT ? GLPK_IPT : GLPK_SIMPLEX),
            true); 
  
  std::cout << "scl_{" << G << "}( " << C << ") = " << scl << " = " << scl.get_d() << "\n";    //output the answer
  
  return 0;
}
