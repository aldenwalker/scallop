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

#include "scylla.h"
#include "scylla_classes.h"
#include "../rational.h"
#include "../lp.h"
#include "../word.h"


using namespace SCYLLA;


/*****************************************************************************
 * Make the list of group polygons and rectangles. 
 * ***************************************************************************/
void SCYLLA::compute_group_teeth_and_rectangles(Chain &C, 
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
void SCYLLA::compute_central_polys(Chain &C, 
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



void SCYLLA::print_central_polys(std::vector<CentralPolygon> &CP, 
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


void SCYLLA::print_group_teeth_and_rectangles(std::vector<GroupTooth> &GT,
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
      


/****************************************************************************
 * The columns are the central polygons, followed by the group teeth, 
 * followed by the group rectangles
 * 
 * the rows are the following:
 * (1) the interface edges
 * (2) the central edge *pairs*
 * (3) the group teeth restrictions
 * (4) the word constraints (N=1)
 * 
 * it's this order for simplicity later
 * ***************************************************************************/
void SCYLLA::scylla_lp(Chain& C, 
               InterfaceEdgeList &IEL,
               CentralEdgePairList &CEL, 
               std::vector<CentralPolygon> &CP,
               std::vector<GroupTooth> &GT,
               std::vector<GroupRectangle> &GR,
               Rational* scl, 
               std::vector<Rational>* solution_vector, 
               SparseLPSolver solver, 
               bool WRITE_LP,
               std::string LP_filename,
               bool cl_not_scl,
               int VERBOSE,
               int LP_VERBOSE) {
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
  num_rows += num_words;
  
  num_cols = CP.size() + GT.size() + GR.size();
  
  if (VERBOSE > 2) {
    std::cout << "Num equality rows: " << num_equality_rows
    << "Num rows: " << num_rows << "\n";
  }
  
  if (VERBOSE>1) {
    std::cout << "Started linear programming setup\n";
  }
  
  
  //Create the LP problem
  SparseLP LP(solver, num_rows, num_cols);
  if (cl_not_scl) {
    for (i=0; i<num_cols; ++i) {
      LP.set_col_type(i, INT);
    }
  }
  
  for(i=0; i<(int)num_equality_rows; i++){
    LP.set_RHS(i, 0); // glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
    LP.set_equality_type(i, EQ);
    if (VERBOSE > 2) {
      std::cout << "Set row " << i << " bounded to " << 0 << "\n";
    }
  }
  for(i=0; i<(int)num_words; i++){
    //RHS[num_equality_rows+i] = C.weights[i];
    LP.set_RHS(num_equality_rows+i, C.weights[i]);
    LP.set_equality_type(num_equality_rows+i, EQ);
    //glp_set_row_bnds(lp, 
    //                  num_equality_rows+i+1, 
    //                  GLP_FX, 
    //                  C.weights[i], 
    //                  C.weights[i]);	
    if (VERBOSE > 2) {
      std::cout << "Set row " << num_equality_rows+i << " bounded to " << C.weights[i] << "\n";
    }
  }
  
  
  //COLS
  for(i=0; i<(int)CP.size(); i++){
    //glp_set_col_bnds(lp, i+1, GLP_LO, 0.0, 0.0);
    //objective[i] = -CP[i].chi_times_2(); //glp_set_obj_coef(lp, i+1, -CP[i].chi_times_2());
    LP.set_obj(i, -CP[i].chi_times_2());
    if (VERBOSE>2) {
      std::cout << "Set objective " << i << " to " << -CP[i].chi_times_2() << "\n";
    }
  }
  offset = CP.size();
  for (i=0; i<(int)GT.size(); i++) {
    //glp_set_col_bnds(lp, offset+i+1, GLP_LO, 0.0, 0.0);
    //objective[offset + i] = -GT[i].chi_times_2(C); //glp_set_obj_coef(lp, offset+i+1, -GT[i].chi_times_2(C));
    LP.set_obj(offset+i, -GT[i].chi_times_2(C));
    if (VERBOSE>2) {
      std::cout << "GT Set objective " << offset+i << " to " << -GT[i].chi_times_2(C) << "\n";
    }
  }
  offset = CP.size() + GT.size();
  for (i=0; i<(int)GR.size(); i++) {
    //glp_set_col_bnds(lp, offset+i+1, GLP_LO, 0.0, 0.0);
    //objective[offset+i] = 0; //glp_set_obj_coef(lp, offset+i+1, 0);
    LP.set_obj(offset+i, 0);
    if (VERBOSE>2) {
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
    if (VERBOSE > 2) {
      std::cout << "CP number " << i << ":\n";
    }
  }
  
  if (VERBOSE>1) { 
    std::cout << "Loaded central polygon edge constraints\n";
  }     
  
  //GROUP TEETH and RECTANGLES
  offset = CP.size();
  for (m=0; m<(int)GT.size(); m++) {
    if (GT[m].inverse) {
      if (VERBOSE > 2) {
        std::cout << GT[m] << "\n";
      }
      GT[m].compute_ia_etc_for_edges(offset + m, 
                                     C, 
                                     IEL, 
                                     group_teeth_rows_inv, 
                                     LP);
    } else {        
      if (VERBOSE > 2) {
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
    if (VERBOSE>2) {
      std::cout << GR[m] << "\n";
    }
    GR[m].compute_ia_etc_for_edges(offset + m, IEL, LP);
  }
  
  
  if (VERBOSE>1) {
    std::cout << "Loaded group constraints\n";
  }
  
  //word constraints: for every group rectangle and group polygon, for every edge, put a 1 in the 
  //row corresponding to the word for the first letter
  offset = CP.size();
  for (j=0; j<(int)GT.size(); j++) {
    GT[j].compute_ia_etc_for_words(offset + j, C, num_equality_rows, LP);
  }
  offset = CP.size() + GT.size();
  for (j=0; j<(int)GR.size(); j++) {
    if (VERBOSE > 2) {
      std::cout << "word cons GR " << GR[j] << "\n";
    }
    GR[j].compute_ia_etc_for_words(offset + j, 
                                   C, 
                                   num_equality_rows, 
                                   IEL,
                                   LP);
  }
  
  if (VERBOSE > 1) {
    std::cout << "Loaded word constraints\n";
  }
  
  
  if (WRITE_LP) {
    LP.write_to_file(LP_filename);    
    return;
  }
  
  if (VERBOSE > 2) {
    std::cout << "LP problem:\n";
    LP.print_LP();
  }
  
  LP.solve(VERBOSE);
  
  LP.get_optimal_value(*scl);
  LP.get_soln_vector(*solution_vector);
  
  
}



void SCYLLA::write_solution_to_fatgraph(std::string& output_filename,
                                        Chain& C,
                                        InterfaceEdgeList& IEL,
                                        CentralEdgePairList &CEL, 
                                        std::vector<CentralPolygon> &CP,
                                        std::vector<GroupTooth> &GT,
                                        std::vector<GroupRectangle> &GR,
                                        std::vector<Rational>& solution_vector,
                                        int verbose ) {}

void SCYLLA::scylla(int argc, char** argv) {
  int current_arg = 0;
  //int i;
  int VERBOSE = 1;
  int LP_VERBOSE = 0;
  SparseLPSolver solver = GLPK_SIMPLEX;
  bool LIMIT_CENTRAL_SIDES = false;
  bool WRITE_LP = false;
  bool CL = false;
  std::string LP_filename;
  bool WRITE_FATGRAPH = false;
  std::string fatgraph_file = "";
  
  if (argc < 1 || std::string(argv[0]) == "-h") {
    std::cout << "usage: ./scallop -cyclic [-h] [-v[n]] [-o <filename>] [-L <filename>] [-l] [-C] [-m<GLPK,GIPT,EXLP,GUROBI>] <gen string> <chain>\n";
    std::cout << "\twhere <gen string> is of the form <gen1><order1><gen2><order2>...\n";
    std::cout << "\te.g. a5b0 computes in Z/5Z * Z\n";
    std::cout << "\tand <chain> is an integer linear combination of words in the generators\n";
    std::cout << "\te.g. ./scallop -cyclic a5b0 aabaaaB\n";
    std::cout << "\tIf the gen string is omitted, the group is assumed to be free\n";
    std::cout << "\t-h: print this message\n";
    std::cout << "\t-v[n]: verbose output (n=0,1,2,3); 0 gives quiet output\n";
    std::cout << "\t-o <filename>: write out *a* (not necessarily *the*) minimal surface as a fatgraph\n";
    std::cout << "\t-L <filename>: write out a sparse lp to the filename .A, .b, and .c\n";
    std::cout << "\t-C compute commutator length (not scl)\n";
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
        if (argv[current_arg][5] == 'O') {
          solver = GUROBI_SIMPLEX;
        } else {
          solver = GUROBI_IPT;
        }
      } else if (argv[current_arg][2] == 'E') {
        solver = EXLP;
      }
      
    } else if (argv[current_arg][1] == 'v') {
      if (argv[current_arg][2] == '\0') {
        VERBOSE = 2;
      } else {
        VERBOSE = atoi(&argv[current_arg][2]);
      }
      
    } else if (argv[current_arg][1] == 'V') {
      LP_VERBOSE = 1;
      
    } else if (argv[current_arg][1] == 'L') {
      WRITE_LP = true;
      LP_filename = std::string(argv[current_arg+1]);
      current_arg++;
      
    } else if (argv[current_arg][1] == 'l') {
      LIMIT_CENTRAL_SIDES = true;
    
    } else if (argv[current_arg][1] == 'C') {
      CL = true;
    
      
    } else if (argv[current_arg][1] == 'o') {
      WRITE_FATGRAPH = true;
      fatgraph_file = std::string(argv[current_arg+1]);
      current_arg++;
    }
    
    current_arg++;
  }
  
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
  CyclicProduct G(G_in); 
  
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
  
  Rational scl;
  std::vector<Rational> solution_vector(0);                           //run the LP
  scylla_lp(C, IEL, CEL, CP, GT, GR, 
            &scl, 
            &solution_vector, 
            solver,
            WRITE_LP, LP_filename,
            CL,
            VERBOSE,
            LP_VERBOSE); 
  
  if (WRITE_LP) {
    std::cout << "Wrote linear program\n";
    return;
  }
  
  if (scl == -1) {
    std::cout << "There was some linear programming error\n";
    return;
  }
  if (CL) {
    scl = scl + Rational(1,2);
  }
  if (VERBOSE>0) {
    if (CL) {
      std::cout << "cl_{" << G.short_rep() << "}(" << C << ") = " << scl << " = " << scl.get_d() << "\n";
    } else {
      std::cout << "scl_{" << G.short_rep() << "}( " << C << ") = " << scl << " = " << scl.get_d() << "\n";    //output the answer
    } 
  } else {
    std::cout << scl.get_d() << "\n";
  }
  
  if (WRITE_FATGRAPH) {
    write_solution_to_fatgraph(fatgraph_file,
                               C, IEL, CEL, CP, GT, GR, solution_vector, VERBOSE);
  }
  
  return;
}
