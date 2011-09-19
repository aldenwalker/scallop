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
void compute_group_teeth_mouths_polygons_and_rectangles(Chain &C, 
                                                         InterfaceEdgeList &IEL,
                                                         std::vector<GroupEdgeList> &GEL,
                                                         std::vector<std::vector<GroupTooth> > &GT,
                                                         std::vector<std::vector<GroupMouth> > &GM,
                                                         std::vector<std::vector<GroupPolygon> > &GP,
                                                         std::vector<std::vector<GroupRectangle> > &GR) {
  int i,j,k,l,m,n;
  int num_groups = (C.G)->num_groups();
  GroupRectangle temp_group_rect;
  GroupPolygon temp_group_poly;
  GroupTooth temp_group_tooth;
  GroupMouth temp_group_mouth;
  Multiset letter_selection;
  int temp_letter_1, temp_letter_2;
  
  GP.resize(num_groups);
  GR.resize(num_groups);
  GT.resize(num_groups);
  GM.resize(num_groups);
  
  for (i=0; i<num_groups; i++) {
    
    //first, let's do the rectangles
    GR[i].resize(0);
    for (j=0; j<(int)C.regular_letters[i].size(); j++) {
      for (k=0; k<(int)C.inverse_letters[i].size(); k++) {
        temp_group_rect.first = IEL.get_index_from_group_side(C.regular_letters[i][j], 
                                                              C.inverse_letters[i][k]);
        temp_group_rect.last = IEL.get_index_from_group_side(C.inverse_letters[i][k], 
                                                             C.regular_letters[i][j]);
        GR[i].push_back(temp_group_rect);
      }
    }
  
    if ((C.G)->orders[i] == 0) {
      GT[i].resize(0);
      GM[i].resize(0);
      GP[i].resize(0);
      continue;
    }

    //now we make the group teeth and mouths
    GT[i].resize(0);   
    GM[i].resize(0);
    temp_group_mouth.inverse = false;
    temp_group_mouth.group_index = i;
    temp_group_tooth.inverse = false;
    temp_group_tooth.group_index = i;
    for (j=0; j<(int)C.regular_letters[i].size(); j++) {
      temp_group_mouth.first = C.regular_letters[i][j];
      temp_group_mouth.first_letter_index = j;
      for (k=0; k<(int)C.regular_letters[i].size(); k++) {
        temp_group_mouth.last = C.regular_letters[i][k];
        temp_group_mouth.last_letter_index = k;
        GM[i].push_back(temp_group_mouth);
        
        //teeth for the 0th position
        temp_group_tooth.first = j;
        temp_group_tooth.position = 0;
        temp_group_tooth.group_mouth_index = GM[i].size()-1;
        for (m=0; m<(int)C.regular_letters[i].size(); m++) {
          temp_group_tooth.last = m;
          GT[i].push_back(temp_group_tooth);
        }
        //teeth for the orderth position
        temp_group_tooth.last = k;
        temp_group_tooth.position = (C.G)->orders[i]-1;
        temp_group_tooth.group_mouth_index = GM[i].size()-1;
        for (m=0; m<(int)C.regular_letters[i].size(); m++) {
          temp_group_tooth.first = m;
          GT[i].push_back(temp_group_tooth);
        }     
        //teeth for the middle positions
        for (l=1; l<(int)(C.G)->orders[i]-1; l++) {
          temp_group_tooth.position = l;
          temp_group_tooth.group_mouth_index = GM[i].size()-1;
          for (m=0; m<(int)C.regular_letters[i].size(); m++) {
            temp_group_tooth.first = m;
            for (n=0; n<(int)C.regular_letters[i].size(); n++) {
              temp_group_tooth.last = n;
              GT[i].push_back(temp_group_tooth);
            }
          }
        }
      }
    }
    temp_group_mouth.inverse = true;
    temp_group_tooth.inverse = true;
    for (j=0; j<(int)C.inverse_letters[i].size(); j++) {
      temp_group_mouth.first = C.inverse_letters[i][j];
      temp_group_mouth.last_letter_index = j;
      for (k=0; k<(int)C.inverse_letters[i].size(); k++) {
        temp_group_mouth.last = C.inverse_letters[i][k];
        temp_group_mouth.last_letter_index = k;
        GM[i].push_back(temp_group_mouth);
        
        //teeth for the 0th position
        temp_group_tooth.first = j;
        temp_group_tooth.position = 0;
        temp_group_tooth.group_mouth_index = GM[i].size()-1;
        for (m=0; m<(int)C.inverse_letters[i].size(); m++) {
          temp_group_tooth.last = m;
          GT[i].push_back(temp_group_tooth);
        }
        //teeth for the orderth position
        temp_group_tooth.last = k;
        temp_group_tooth.position = (C.G)->orders[i]-1;
        temp_group_tooth.group_mouth_index = GM[i].size()-1;
        for (m=0; m<(int)C.inverse_letters[i].size(); m++) {
          temp_group_tooth.first = m;
          GT[i].push_back(temp_group_tooth);
        }     
        //teeth for the middle positions
        for (l=1; l<(int)(C.G)->orders[i]-1; l++) {
          temp_group_tooth.position = l;
          temp_group_tooth.group_mouth_index = GM[i].size()-1;
          for (m=0; m<(int)C.inverse_letters[i].size(); m++) {
            temp_group_tooth.first = m;
            for (n=0; n<(int)C.inverse_letters[i].size(); n++) {
              temp_group_tooth.last = n;
              GT[i].push_back(temp_group_tooth);
            }
          }
        }
      }
    }
      
    


    
    //now the group internal polygons -- note there can be just any four letters
    //on the corners, though we require there to be at least one nontrivial edge
    GP[i].resize(0);
    temp_group_poly.edges.resize(4);
    temp_group_poly.group = i;
    if (C.regular_letters[i].size() > 0) {
      letter_selection = Multiset(4, 0, C.regular_letters[i].size());
      do {
        if (letter_selection[0] == letter_selection[1]) {
          continue;                                       //this is a waste of time
        }
        for (j=0; j<4; j++) {
          temp_letter_1 = C.regular_letters[i][letter_selection[j]];
          temp_letter_2 = C.regular_letters[i][letter_selection[(j+1)%4]];
          temp_group_poly.edges[j] = GEL[i].get_index(temp_letter_1, temp_letter_2);
        }
        GP[i].push_back(temp_group_poly);
      } while (1 != letter_selection.next());
    }
    
    if (C.inverse_letters[i].size() > 0) {
      letter_selection = Multiset(4, 0, C.inverse_letters[i].size());
      do {
        if (letter_selection[0] == letter_selection[1]) {
          continue;                                       //this is a waste of time
        }
        for (j=0; j<4; j++) {
          temp_letter_1 = C.inverse_letters[i][letter_selection[j]];
          temp_letter_2 = C.inverse_letters[i][letter_selection[(j+1)%4]];
          temp_group_poly.edges[j] = GEL[i].get_index(temp_letter_1, temp_letter_2);
        }
        GP[i].push_back(temp_group_poly);
      } while (1 != letter_selection.next());
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
  int i;
  bool VERBOSE = false;
  
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
  
  std::vector<GroupEdgeList> GEL((C.G)->num_groups());
  for (i=0; i<(C.G)->num_groups(); i++) {
    GEL[i] = GroupEdgeList(C, i);
    if (VERBOSE) GEL[i].print(std::cout);
  }
  
  std::vector<CentralPolygon> CP;
  compute_central_polys(C, IEL, CEL, CP);
  if (VERBOSE) print_central_polys(CP, std::cout, 0);
  
  std::vector<std::vector<GroupTooth> > GT;
  std::vector<std::vector<GroupMouth> > GM;
  std::vector<std::vector<GroupRectangle> > GR;
  std::vector<std::vector<GroupPolygon> > GP;
  compute_group_teeth_mouths_polygons_and_rectangles(C, IEL, GEL, GT, GM, GP, GR);
  if (VERBOSE) print_group_teeth_mouths_polys_and_rectangles(GT, GM, GP, GR, std::cout, 0);
  
   
  rational scl;
  std::vector<rational> solution_vector(0);                           //run the LP
  scylla_lp(C, GEL, IEL, CEL, CP, GT, GM, GP, GR, &scl, &solution_vector, GLPK_DOUBLE, VERBOSE); 
  
  std::cout << "scl_{" << G << "}( " << C << ") = " << scl << " = " << scl.get_d() << "\n";    //output the answer
  
  return 0;
}
