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
                                       std::vector<GroupTooth > &GT,
                                       std::vector<GroupRectangle > &GR) {
  int i,j,k,l,m;
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
void compute_central_polys(Chain &C, 
                           InterfaceEdgeList &IEL, 
                           std::vector<CentralPolygon> &CP,
                           bool limit_central_sides) {
  int i,j,k,l;
  int e1L1, e1L2, e2L1, e2L2, e3L1, e3L2;  //edge1, letter1, etc
  CentralPolygon temp_central_poly;
  
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
  
  if (limit_central_sides) {
    return;
  }
  
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



void print_central_polys(std::vector<CentralPolygon> &CP, 
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


void print_group_teeth_and_rectangles(std::vector<GroupTooth> &GT,
                                  std::vector<GroupRectangle> &GR,
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
      




int main(int argc, char* argv[]) {
  int current_arg = 1;
  //int i;
  int VERBOSE = 1;
  bool IPT = false;
  bool LIMIT_CENTRAL_SIDES = false;
  
  if (argc < 3 || std::string(argv[1]) == "-h") {
    std::cout << "usage: ./scyllop [-h] [-v[n]] [-l] [-i] <gen string> <chain>\n";
    std::cout << "\twhere <gen string> is of the form <gen1><order1><gen2><order2>...\n";
    std::cout << "\te.g. a5b0 computes in Z/5Z * Z\n";
    std::cout << "\tand <chain> is an integer linear combination of words in the generators\n";
    std::cout << "\te.g. ./scyllop a5b0 aabaaaB\n";
    std::cout << "\t-h: print this message\n";
    std::cout << "\t-v[n]: verbose output (n=0,1,2,3); 0 gives quiet output\n";
    std::cout << "\t-l: limit the number of central sides to 1 (see readme)\n";
    std::cout << "\t-i: use the interior point LP method (faster but rational output is sometimes wrong)\n";
    exit(0);
  }
  while (argv[current_arg][0] == '-') {
    if (argv[current_arg][1] == 'i') {
      IPT = true;
    } else if (argv[current_arg][1] == 'v') {
      if (argv[current_arg][2] == '\0') {
        VERBOSE = 2;
      } else {
        VERBOSE = atoi(&argv[current_arg][2]);
      }
    } else if (argv[current_arg][1] == 'l') {
      LIMIT_CENTRAL_SIDES = true;
    }
    current_arg++;
  }
  
  std::string G_in = std::string(argv[current_arg]);
  CyclicProduct G(G_in);                                                         //create the group
  current_arg++;
  
  Chain C(&G, &argv[current_arg], argc-current_arg);                              //process the chain argument

  if (VERBOSE>1) {
    std::cout << "Group: " << G << "\n";
    std::cout << "Chain: " << C << "\n";
    if (VERBOSE>2) {
      std::cout << "Letters:\n";
      C.print_letters(std::cout);
      std::cout << "Group letters:\n";
      C.print_group_letters(std::cout);
    }
  }
  
  InterfaceEdgeList IEL(C);
  if (VERBOSE>1) IEL.print(std::cout);
  
  CentralEdgePairList CEL(C);
  if (VERBOSE>1) CEL.print(std::cout);
  
  std::vector<CentralPolygon> CP;
  compute_central_polys(C, IEL, CP, LIMIT_CENTRAL_SIDES);
  if (VERBOSE > 1) {
    std::cout << "computed polys (" << CP.size() << ")\n"; std::cout.flush();
    print_central_polys(CP, std::cout, VERBOSE);
  }
  
  std::vector<GroupTooth> GT;
  std::vector<GroupRectangle> GR;
  compute_group_teeth_and_rectangles(C, GT, GR);
  if (VERBOSE > 1) {
    std::cout << "computed group teeth and rectangles\n"; std::cout.flush();
    print_group_teeth_and_rectangles(GT, GR, std::cout, VERBOSE);
  }
  
  rational scl;
  std::vector<rational> solution_vector(0);                           //run the LP
  scylla_lp(C, IEL, CEL, CP, GT, GR, 
            &scl, 
            &solution_vector, 
            (IPT ? GLPK_IPT : GLPK_SIMPLEX),
            VERBOSE); 
  
  if (VERBOSE>0) {
    std::cout << "scl_{" << G << "}( " << C << ") = " << scl << " = " << scl.get_d() << "\n";    //output the answer
  } else {
    std::cout << scl.get_d() << "\n";
  }
  return 0;
}
