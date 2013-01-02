/****************************************************************************
* compute scl(w), where w is the collection of words of length ell
* in Gamma.  Also compute sup phi(Gamma)/2D(phi) for counting 
* functions phi of a certain given length ell.  
* See the included traintracks.pdf for details.
*
* By Alden Walker.  Implements the algorithm described in 
* "Traintracks and scl" by Danny Calegari
*****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <ctype.h>

#include "trollop_classes.h"
#include "trollop.h"
#include "../rational.h"
#include "../word.h"
#include "../lp.h"



using namespace TROLLOP;



/***************************************************************************
 a helper function to read in a matrix. The matrix file is formatted as follows:
 <rows> <columns>
 <row which has all the words in the order used in the matrix>
 followed by the matrix, whitespace delimited, row by row 
 ***************************************************************************/
void TROLLOP::read_matrix(std::vector<std::vector<int> >& matrix, 
                 std::string matrix_filename, 
                 WordTable &WT,
                 bool rearrange_rows,
                 int VERBOSE) {
  std::fstream mat_file;
  int rows, cols;
  int i,j;
  std::string word;
  
  mat_file.open(matrix_filename.c_str(), std::fstream::in);
  
  if (mat_file.fail()) {
    std::cout << "File open failure\n";
    return;
  }
  
  mat_file >> rows >> cols;
  
  if (VERBOSE>1) {
    std::cout << "Reading matrix from file " << matrix_filename << "\n";
    std::cout << "Matrix has " << rows << " rows and " << cols << " columns\n";
  }
  
  std::vector<int> file_index_to_wordtable_index(cols);
  for (i=0; i<cols; i++) {
    mat_file >> word;
    if (word == "RHS") {
      file_index_to_wordtable_index[i] = cols-1;
    } else {
      file_index_to_wordtable_index[i] = WT.get_index(word);
    }
    if (VERBOSE>3) {
      std::cout << "Read translation table: " << word << ": file: " << i << " WT: " << file_index_to_wordtable_index[i] << "\n";
    }
  }
  matrix.resize(rows);
  for (i=0; i<rows; i++) {
    matrix[i].resize(cols);
  }
  for (i=0; i<rows; i++) {
    for (j=0; j<cols; j++) {
      if (rearrange_rows) {
        mat_file >> matrix[file_index_to_wordtable_index[i]]
                          [file_index_to_wordtable_index[j]];
      } else {
        mat_file >> matrix[i][file_index_to_wordtable_index[j]];
      }
    }
  }
  mat_file.close();
}


/***************************************************************************
 a helper function to read in a vector. The vector file is formatted as follows:
 <length>
 followed by the vector, whitespace delimited
 ***************************************************************************/
void TROLLOP::read_vector(std::vector<int>& vector, 
                 std::string vector_filename,
                 int VERBOSE) {
  std::fstream vec_file;
  int len;
  int i;//,j;
  
  vec_file.open(vector_filename.c_str(), std::fstream::in);
  
  if (vec_file.fail()) {
    std::cout << "File open failure\n";
    return;
  }
  
  vec_file >> len;
  
  if (VERBOSE > 1) {
    std::cout << "Reading matrix from file " << vector_filename << "\n";
    std::cout << "Vector has length" << len << "\n";
  }
  
  vector.resize(len);
  for (i=0; i<len; i++) {
    vec_file >> vector[i];
  }
  vec_file.close();
}





void TROLLOP::compute_triangles(WordTable& WT, 
                       ArcPairList& AL, 
                       int num_copies,
                       std::vector<Triangle>& TR) {
  int v0,v1,v2,a0,a1,a2;
  int side0, side1, side2;
  TR.resize(0);
  Triangle temp_t;
  for (side0=0; side0<num_copies; side0++) {
    for (v0=0; v0<WT.num_verts; v0++) {
      
      for (side1=0; side1<num_copies; side1++) {
        for (v1 = v0+1; v1<WT.num_verts; v1++) {
          a0 = AL.index_of_arc(side0, v0, side1, v1);
          
          for (side2=0; side2<num_copies; side2++) {
            for (v2 = v0+1; v2<WT.num_verts; v2++) {
              if (v2 == v1) {
                continue;
              }
              a1 = AL.index_of_arc(side1, v1, side2, v2);
              a2 = AL.index_of_arc(side2, v2, side0, v0);
              temp_t.v0 = v0; temp_t.v1 = v1; temp_t.v2 = v2;
              temp_t.a0 = a0; temp_t.a1 = a1; temp_t.a2 = a2;
              TR.push_back(temp_t);
            }
          }
        }
      }
    }
  }
}


void TROLLOP::print_triangles(std::vector<Triangle>& TR, std::ostream& os) {
  int i;
  os << "Triangles (" << TR.size() << "):\n";
  for (i=0; i<(int)TR.size(); i++) {
    os << i << ": " << TR[i] << "\n";
  }
}


//for every pair of edges with inverse h values, make a rectangle
//recall it's (a0,e0,a1,e1)
void TROLLOP::compute_rectangles(WordTable& WT, 
                        ArcPairList& AL, 
                        int num_copies, 
                        std::vector<Rectangle>& RE) {
  int i,j;
  RE.resize(0);
  Rectangle temp_r;
  int a0,a1,e0,e1;
  int a0v0, a0v1, a1v0, a1v1;
  int e0_letter, e1_letter;
  int side0, side1;
  for (e0_letter=0; e0_letter<WT.rank; e0_letter++) {
    e1_letter = letter_to_number( 
                    inverse_char(
                        number_to_letter(e0_letter, WT.rank, '\0')
                    ), 
                    WT.rank, '\0');
    for (i=0; i<(int)WT.h_vector[e0_letter].size(); i++) {
      e0 = WT.h_vector[e0_letter][i];
      for (j=0; j<(int)WT.h_vector[e1_letter].size(); j++) {
        e1 = WT.h_vector[e1_letter][j];
        a0v0 = WT.get_edge_dest(e1);
        a0v1 = WT.get_edge_source(e0);
        a1v0 = WT.get_edge_dest(e0);
        a1v1 = WT.get_edge_source(e1);
        for (side0=0; side0<num_copies; side0++) {
          temp_r.side0 = side0;
          for (side1=0; side1<num_copies; side1++) {
            temp_r.side1 = side1;
            a0 = AL.index_of_arc(side1, a0v0, side0, a0v1);
            a1 = AL.index_of_arc(side0, a1v0, side1, a1v1);
            
            temp_r.a0 = a0; temp_r.a1 = a1;
            temp_r.e0 = e0; temp_r.e1 = e1;
            RE.push_back(temp_r);
          }
        } 
      }
    }
  }    
}

void TROLLOP::print_rectangles(std::vector<Rectangle>& RE, std::ostream& os) {
  int i;
  //int sign, index;
  os << "Rectangles (" << RE.size() << "):\n";
  for (i=0; i<(int)RE.size(); i++) {
    os << i << ": " << RE[i] << "\n";
  }
}







template <typename T>
void free_vector(T & t) {
  T tmp_vector;
  t.swap( tmp_vector ); 
  //now the scope of tmp_vector will end
}


/***************************************************************************
 a helper function to collect duplicates
 ***************************************************************************/
void TROLLOP::collect_dups_and_push(std::vector<int> &temp_ia,
                           std::vector<int> &temp_ja,
                           std::vector<int> &temp_ar,
                           std::vector<int> &ia,
                           std::vector<int> &ja,
                           std::vector<int> &ar,
                           int VERBOSE) {
  int j, k, temp;
  if (VERBOSE < 4) {
    for (j=0; j<(int)temp_ia.size(); j++) {
      if (temp_ar[j] == 0) {
        continue;
      }
      temp = 0;
      for (k=0; k<(int)temp_ia.size(); k++) {
        if (temp_ia[k] == temp_ia[j]) {
          temp += temp_ar[k];
          temp_ar[k] = 0;
        }
      }
      ja.push_back(temp_ja[j]);
      ia.push_back(temp_ia[j]);
      ar.push_back(temp);
    }
  } else {
    for (j=0; j<(int)temp_ia.size(); j++) {
      if (temp_ar[j] == 0) {
        continue;
      }
      temp = 0;
      for (k=0; k<(int)temp_ia.size(); k++) {
        if (temp_ia[k] == temp_ia[j]) {
          temp += temp_ar[k];
          temp_ar[k] = 0;
        }
      }
      ja.push_back(temp_ja[j]);
      ia.push_back(temp_ia[j]);
      ar.push_back(temp);
      std::cout << "Put " << temp_ia[j] << ", " << temp_ja[j] << ", " << temp << ".\n";
    }
  }    
}


/***************************************************************************
 * Write the lp to a file
 * *************************************************************************/
void TROLLOP::write_lp(WordTable& WT, 
              WVec& C, 
              ArcPairList& AL, 
              int num_copies,
              std::vector<Triangle>& TR, 
              std::vector<Rectangle>& RE, 
              bool DO_SUP,
              bool MAT_COMP,
              bool MAT_SEPARATE_DOMAIN,
              std::vector<std::vector<int> >& M,
              std::vector<std::vector<int> >& N,
              std::vector<int>& b,
              int VERBOSE,
              int LP_VERBOSE,
              std::string& programFile) {
  int i,j;
  std::string AFile = programFile + ".A";
  std::string bFile = programFile + ".b";
  std::string cFile = programFile + ".c";
  std::fstream outAFile;
  std::fstream outbFile;
  std::fstream outcFile;
  //these are all just for removing duplicates:
  std::vector<int> temp_ia2(0);
	std::vector<int> temp_ja2(0);
	std::vector<int> temp_ar2(0); 
  std::vector<int> temp_ia(0);
	std::vector<int> temp_ja(0);
	std::vector<int> temp_ar(0); 
  int row, col, val;
  int sign, index;  
  int M_rows;
  int M_cols;
  int N_rows, N_cols;
  int b_len;
  int num_rows;
  int num_cols;
  int offset;
  int col_offset, row_offset, row1, row2;
  
  //the number of rows is the number of arcs plus the number of 
  //edges
  if (MAT_COMP) {
    M_rows = M.size(); 
    M_cols = M[0].size();
    N_rows = N.size();
    N_cols = N[0].size();
    b_len = b.size(); 
    if (b_len != N_rows) {
      std::cout << "N, b size mismatch\n";
    }
    if (N_cols != WT.num_edges || M_cols != WT.num_edges || M_cols != M_rows) {
      std::cout << "Matrix size mismatch\n";
    }
    if (MAT_SEPARATE_DOMAIN) {
      //a row for each of the arcs, rows for the edges, plus rows for the rows in N
      //side0 and side1 are forced to be weights because they bound a side of the surface
      //plus we need rows to specify that side1 = M*side0
      num_rows = AL.num_arcs + 2*WT.num_edges + b_len + M_rows;
      //a column for each rectangle, each triangle, and each edge (for sides 0 and 1)
      num_cols = RE.size() + TR.size() + 2*WT.num_edges;
    } else {
      //a row for each of the arcs, rows for the edges, plus rows for the rows in N, plus rows to make sure x is a weight
      num_rows = AL.num_arcs + WT.num_edges + b_len + WT.num_verts;
      //a column for each rectangle, each triangle, and each edge (for the M matrix)
      num_cols = RE.size() + TR.size() + M_rows;
    }
  } else {
    num_rows = AL.num_arcs + WT.num_edges;
    //the number of columns is the number of rectangles and triangles, 
    //plus a single column if DO_SUP
    num_cols = RE.size() + TR.size() + (DO_SUP ? 1 : 0);
  }
  
  
  
  /////  A matrix
  outAFile.open(AFile.c_str(), std::fstream::out);
  outAFile << num_rows << " " << num_cols << " 0\n";
  for (i=0; i<(int)RE.size(); i++) {
    //the arcs
    temp_ia.resize(0); temp_ja.resize(0); temp_ar.resize(0);
    temp_ia2.resize(0); temp_ja2.resize(0); temp_ar2.resize(0);
    extract_signed_index(&sign, &index, RE[i].a0);
    temp_ia.push_back(index+1);
    temp_ja.push_back(i+1);
    if (sign < 0) {
      temp_ar.push_back(-1);
    } else {
      temp_ar.push_back(1);
    }
    extract_signed_index(&sign, &index, RE[i].a1);
    temp_ia.push_back(index+1);
    temp_ja.push_back(i+1);
    if (sign < 0) {
      temp_ar.push_back(-1);
    } else {
      temp_ar.push_back(1);
    }
    
    //the edges -- still for MAT_COMP
    temp_ia.push_back(AL.num_arcs+RE[i].e0 + 1);
    temp_ja.push_back(i+1);
    temp_ar.push_back(1);
    temp_ia.push_back(AL.num_arcs+RE[i].e1 + 1);
    temp_ja.push_back(i+1);
    temp_ar.push_back(1);
    collect_dups_and_push(temp_ia, temp_ja, temp_ar, 
                          temp_ia2, temp_ja2, temp_ar2, VERBOSE);
    for (j=0; j<(int)temp_ia2.size(); j++) {
      outAFile << temp_ia2[j] << " " << temp_ja2[j] << " " << temp_ar2[j] << "\n";
    }
    
  }
  
  if (VERBOSE>1) {
    std::cout << "wrote rectange columns to file\n";
  }
  
  offset = RE.size();
  for (i=0; i<(int)TR.size(); i++) {
    for (j=0; j<3; j++) {
      extract_signed_index(&sign, &index, (j==0 ? TR[i].a0 : (j==1 ? TR[i].a1 : TR[i].a2) ) ) ;
      row = index + 1; //ia.push_back(index+1);
      col = offset + i + 1; //ja.push_back(offset+i+1);
      if (sign < 0) {
        val = -1; //ar.push_back(-1);
      } else {
        val = 1;  //ar.push_back(1);
      }
      outAFile << row << " " << col << " " << val << "\n";
    }
  }
  
  if (VERBOSE>1) {
    std::cout << "wrote triangle columns to file\n";
  }
  
  
  if (MAT_COMP) {
    //load in the M matrix
    col_offset = RE.size() + TR.size();
    row_offset = AL.num_arcs;
    for (i=0; i<WT.num_edges; i++) { //these are the COLUMNS (the weights which form the input weight)
      for (j=0; j<WT.num_edges; j++) { //these are the ROWS  (the edge rows
        //for each column, compute how much it contributes to each row, 
        //and add these as negatives.
        //that is, we are simply putting the ith column of M in, with minus signs
        if (M[j][i] == 0) continue;
        outAFile << row_offset + j + 1 << " " << col_offset + i + 1 << " " <<  -M[j][i] << "\n";
      }
    }
    
    //load in the N matrix (the b's are already set on the RHS)
    row_offset = AL.num_arcs + WT.num_edges;
    col_offset = RE.size() + TR.size();
    for (i=0; i<WT.num_edges; i++) { //these are the columns (weights form the input x weight)
      for (j=0; j<b_len; j++) { //these are the ROWS, the output rows of the N matrix
        if (N[j][i] == 0) continue;
        outAFile << row_offset + j + 1 << " " << col_offset + i + 1 << " " <<  N[j][i] << "\n";
      }
    }
    if (VERBOSE>1) { 
      std::cout << "Loaded the M and N matrices\n";
    }
    
    //load in the constraints that make sure the input is a weight
    row_offset = AL.num_arcs + WT.num_edges + b_len;
    col_offset = RE.size() + TR.size();
    for (i=0; i<WT.num_edges; i++) {  //these are the columns (input x weights)
      temp_ia.resize(0); temp_ja.resize(0); temp_ar.resize(0);
      temp_ia2.resize(0); temp_ja2.resize(0); temp_ar2.resize(0);
      row1 = WT.get_edge_source(i);
      temp_ia.push_back(row_offset + row1 + 1);
      temp_ja.push_back(col_offset + i + 1);
      temp_ar.push_back(1);
      row2 = WT.get_edge_dest(i);
      temp_ia.push_back(row_offset + row2 + 1);
      temp_ja.push_back(col_offset + i + 1);
      temp_ar.push_back(-1);
      collect_dups_and_push(temp_ia, temp_ja, temp_ar, temp_ia2, temp_ja2, temp_ar2, VERBOSE);
      for (j=0; j<(int)temp_ia2.size(); j++) {
        outAFile << temp_ia2[j] << " " << temp_ja2[j] << " " << temp_ar2[j] << "\n";
      }
    }
  }
  
  
  if (DO_SUP) {
    offset = AL.num_arcs;
    for (i=0; i<WT.num_edges; i++) {
      row = offset + i + 1; //ia.push_back(offset + i + 1);
      col = num_cols; //ja.push_back(num_cols);
      val = -1; //ar.push_back(-1);
      outAFile << row << " " << col << " " << val << "\n";
    }
  }
  
  outAFile.close(); 
  if (VERBOSE>1) {
    std::cout << "Finished with A file\n";
  }
  
  /////// b vector
  outbFile.open(bFile.c_str(), std::fstream::out);
  
  for (i=0; i<AL.num_arcs; i++) {          //all the arc constraints are zero
    outbFile << "0\n";
  }
  
  if (MAT_COMP) {
    offset = AL.num_arcs;
    for (i=0; i<WT.num_edges; i++) {                 //the weight rows all equal 0
      outbFile << 0 << "\n";
    }
    
    offset = AL.num_arcs + WT.num_edges;    //the Nx = b rows
    for (i=0; i<b_len; i++) {
      outbFile << b[i] << "\n";
    }
    
    offset = AL.num_arcs + WT.num_edges + b_len;  //the rows making sure it's a traintrack
    for (i=0; i<WT.num_verts; i++) {
      outbFile << 0 <<"\n";
    }
    
  } else {
    offset = AL.num_arcs;
    for (i=0; i<WT.num_edges; i++) {                 //the word RHS is the weight
      outbFile << C.index_coefficients[i] << "\n";
    }
  }
  
  outbFile.close();
  if (VERBOSE>1) {
    std::cout << "Finished with b vector\n";
  }  
  
  /////// c objective function
  outcFile.open(cFile.c_str(), std::fstream::out);
  
  for (i=0; i<(int)RE.size(); i++) {   //chi is zero
    outcFile << "0\n";
  }
  offset = RE.size();
  for (i=0; i<(int)TR.size(); i++) {   //chi is -1/2  (we minimize -2chi)
    outcFile << "1\n";
  }  
  
  if (MAT_COMP) {  //add the columns for the x vector (all zero coefficients)
    for (i=0; i<WT.num_edges; i++) {
      outcFile << "0\n";
    }
  }
  
  if (DO_SUP) {
    outcFile << "0\n";
  }
  outcFile.close();
  if (VERBOSE>1) {
    std::cout << "Finished with c objective vector\n";
  }  
}



/****************************************************************************
 run the lp
 * ***************************************************************************/
void TROLLOP::trollop_lp(WordTable& WT, 
                WVec& C, 
                ArcPairList& AL, 
                int num_copies,
                std::vector<Triangle>& TR, 
                std::vector<Rectangle>& RE, 
                bool DO_SUP,
                bool MAT_COMP,
                bool MAT_SEPARATE_DOMAIN,
                std::vector<std::vector<int> >& M,
                std::vector<std::vector<int> >& N,
                std::vector<int>& b,
                Rational& ans,
                std::vector<Rational>& solution_vector,
                SparseLPSolver solver,
                int VERBOSE,
                int LP_VERBOSE) {
  std::vector<int> temp_ia(0);
	std::vector<int> temp_ja(0);
	std::vector<int> temp_ar(0); 
  int sign, index;
  int i,j;
  int num_cols, offset, num_rows;
  int M_rows;
  int M_cols;
  int N_rows, N_cols;
  int b_len;
  int row_offset, col_offset;
  int row1, row2;
  
  //the number of rows is the number of arcs plus the number of 
  //edges
  if (MAT_COMP) {
    M_rows = M.size();
    M_cols = M[0].size();
    N_rows = N.size();
    N_cols = N[0].size();
    b_len = b.size();
    if (b_len != N_rows) {
      std::cout << "N, b size mismatch\n";
    }
    if (N_cols != WT.num_edges || M_cols != WT.num_edges || M_cols != M_rows) {
      std::cout << "Matrix size mismatch\n";
    }
    if (MAT_SEPARATE_DOMAIN) {
      //a row for each of the arcs, rows for the edges, plus rows for the rows in N
      //side0 and side1 are forced to be weights because they bound a side of the surface
      //plus we need rows to specify that side1 = M*side0
      num_rows = AL.num_arcs + 2*WT.num_edges + b_len + M_rows;
      //a column for each rectangle, each triangle, and each edge (for sides 0 and 1)
      num_cols = RE.size() + TR.size() + 2*WT.num_edges;
    } else {
      //a row for each of the arcs, rows for the edges, plus rows for the rows in N, plus rows to make sure x is a weight
      num_rows = AL.num_arcs + WT.num_edges + b_len + WT.num_verts;
      //a column for each rectangle, each triangle, and each edge (for the M matrix)
      num_cols = RE.size() + TR.size() + M_rows;
    }
  } else {
    num_rows = AL.num_arcs + WT.num_edges;
    //the number of columns is the number of rectangles and triangles, 
    //plus a single column if DO_SUP
    num_cols = RE.size() + TR.size() + (DO_SUP ? 1 : 0);
  }
  
  
  SparseLP LP(solver, num_rows, num_cols);
  
  if (VERBOSE>1) {
    std::cout << "Started linear programming setup\n";
  }
  
  
  //now we actually do the LP setup  
  
  //ROWS
  
  //we always have the arc rows
  for(i=0; i<(int)AL.num_arcs; i++){
    LP.set_RHS(i, 0);
    LP.set_equality_type(i , EQ);
    //RHS[i+1] = 0; //glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
    if (VERBOSE > 3) {
      std::cout << "Set row " << i << " fixed to " << 0 << "\n";
    }
  }  
  
  if (MAT_COMP) {
    if (MAT_SEPARATE_DOMAIN) {
      offset = AL.num_arcs;
      for (i=0; i<WT.num_edges; i++) {    //set the edge rows for side 0
        LP.set_RHS(offset+i, 0);
        //RHS[offset+i+1] = 0; //glp_set_row_bnds(lp, offset+i+1, GLP_FX, 0, 0);
        if (VERBOSE > 3) {
          std::cout << "Set (edge) row " << offset+i << " fixed to " << 0 << "\n";
        }
      }
      offset = AL.num_arcs + WT.num_edges;
      for (i=0; i<WT.num_edges; i++) {    //set the edge rows for side 1
        LP.set_RHS(offset+i, 0);
        //RHS[offset+i+1] = 0; //glp_set_row_bnds(lp, offset+i+1, GLP_FX, 0, 0);
        if (VERBOSE > 3) {
          std::cout << "Set (edge) row " << offset+i << " fixed to " << 0 << "\n";
        }
      }
      offset = AL.num_arcs + 2*WT.num_edges; 
      for (i=0; i<b_len; i++) {           //set the b rows
        LP.set_RHS(offset+i, b[i]);
        //RHS[offset+i+1] = b[i]; //glp_set_row_bnds(lp, offset+i+1, GLP_FX, b[i], b[i]);        
        if (VERBOSE > 3) {
          std::cout << "Set (b) row " << offset+i << " fixed to " << b[i] << "\n";
        }
      }
      offset = AL.num_arcs + 2*WT.num_edges + b_len;
      for (i=0; i<WT.num_edges; i++) {    //set the matrix constraint rows
        LP.set_RHS(offset+i, 0);
        //RHS[offset+i+1] = 0; //glp_set_row_bnds(lp, offset+i+1, GLP_FX, 0, 0);
        if (VERBOSE > 3) {
          std::cout << "Set (matrix) row " << offset+i << " fixed to " << 0 << "\n";
        }
      }
      
      
    } else {
      offset = AL.num_arcs;      
      for (i=0; i<WT.num_edges; i++) {    //set the edge rows
        LP.set_RHS(offset+i, 0);
        //RHS[offset+i+1] = 0; //glp_set_row_bnds(lp, offset+i+1, GLP_FX, 0, 0);
        if (VERBOSE > 3) {
          std::cout << "Set (edge) row " << offset+i+1 << " fixed to " << 0 << "\n";
        }
      }
      offset = AL.num_arcs + WT.num_edges;  //set the b rows
      for (i=0; i<b_len; i++) {
        LP.set_RHS(offset+i, b[i]);
        //RHS[offset+i+1] = b[i]; //glp_set_row_bnds(lp, offset+i+1, GLP_FX, b[i], b[i]);        
        if (VERBOSE > 3) {
          std::cout << "Set (b) row " << offset+i+1 << " fixed to " << b[i] << "\n";
        }
      }
      offset = AL.num_arcs + WT.num_edges + b_len;  //set the make-sure-x-is-weight rows
      for (i=0; i<WT.num_verts; i++) {
        LP.set_RHS(offset+i, 0);
        //RHS[offset+i+1] = 0; //glp_set_row_bnds(lp, offset+i+1, GLP_FX, 0,0);        
        if (VERBOSE > 3) {
          std::cout << "Set (traintrack) row " << offset+i+1 << " fixed to " << 0 << "\n";
        }
      }
    }
    
  } else {  //usual computation
    offset = AL.num_arcs;
    for (i=0; i<WT.num_edges; i++) {
      LP.set_RHS(offset+i, C.index_coefficients[i]);
      //RHS[offset+i+1] = C.index_coefficients[i]; //glp_set_row_bnds(lp, offset+i+1, GLP_FX, C.index_coefficients[i], C.index_coefficients[i]);
      if (VERBOSE > 3) {
        std::cout << "Set (edge) row " << offset+i+1 << " fixed to " << C.index_coefficients[i] << "\n";
      }
    }
  }  
  
  //COLS
  for(i=0; i<(int)RE.size(); i++){
    //glp_set_col_bnds(lp, i+1, GLP_LO, 0.0, 0.0);
    LP.set_obj(i, 0);
    //objective[i+1] = 0; //glp_set_obj_coef(lp, i+1, 0);
    if (VERBOSE>3) {
      std::cout << "Set objective rectangle" << i+1 << " to " << 0 << "\n";
    }
  }
  offset = RE.size();
  for (i=0; i<(int)TR.size(); i++) {
    //glp_set_col_bnds(lp, offset+i+1, GLP_LO, 0.0, 0.0);
    LP.set_obj(offset+i, 1);
    //objective[offset+i+1] = 1; //glp_set_obj_coef(lp, offset+i+1, 1);
    if (VERBOSE>3) {
      std::cout << "Set objective triangle" << offset+i+1 << " to " << 1 << "\n";
    }
  }
  
  if (MAT_COMP) {
    if (MAT_SEPARATE_DOMAIN) {
      offset = RE.size() + TR.size();
      for (i=0; i<2*WT.num_edges; i++) {
        //glp_set_col_bnds(lp, offset+i+1, GLP_LO, 0, 0);
        LP.set_obj(offset+i, 0);
        //objective[offset+i+1] = 0; //glp_set_obj_coef(lp, offset+i+1, 0); 
        if (VERBOSE>3) {
          std::cout << "Set objective side value" << offset+i+1 << " to " << 0 << "\n";
        }
      }
      
    } else { 
      offset = RE.size() + TR.size();
      for (i=0; i<M_cols; i++) {
        //glp_set_col_bnds(lp, offset+i+1, GLP_LO, 0, 0);
        LP.set_obj(offset+i, 0);
        //objective[offset+i+1] = 0; //glp_set_obj_coef(lp, offset+i+1, 0);      
        if (VERBOSE>3) {
          std::cout << "Set objective M value" << offset+i+1 << " to " << 0 << "\n";
        }
      }
    }
  }
  
  if (DO_SUP) {
    //glp_set_col_bnds(lp, num_cols, GLP_LO, 0.0, 0.0);
    LP.set_obj(num_cols-1, 0);
    //objective[num_cols] = 0;  //glp_set_obj_coef(lp, num_cols, 0);
  }
  
  for (i=0; i<(int)RE.size(); i++) {
    if (VERBOSE > 3) {
      std::cout << "Rectangle " << i << " " << RE[i] << "\n";
    }
    
    //the arcs
    temp_ia.resize(0); temp_ja.resize(0); temp_ar.resize(0);
    extract_signed_index(&sign, &index, RE[i].a0);
    temp_ia.push_back(index);
    temp_ja.push_back(i);
    if (sign < 0) {
      temp_ar.push_back(-1);
    } else {
      temp_ar.push_back(1);
    }
    extract_signed_index(&sign, &index, RE[i].a1);
    temp_ia.push_back(index);
    temp_ja.push_back(i);
    if (sign < 0) {
      temp_ar.push_back(-1);
    } else {
      temp_ar.push_back(1);
    }
    
    //the edges; we still do this even for MAT_COMP (MAT_SEPARATE_DOMAIN is also ok in here)
    temp_ia.push_back(AL.num_arcs + (RE[i].side0 * WT.num_edges) +  RE[i].e0 );
    temp_ja.push_back(i);
    temp_ar.push_back(1);
    temp_ia.push_back(AL.num_arcs + (RE[i].side1 * WT.num_edges) + RE[i].e1 );
    temp_ja.push_back(i);
    temp_ar.push_back(1);
    
    //ACTUALLY push the rectangle stuff into the matrix
    LP.extend_entries_no_dups(temp_ia, temp_ja, temp_ar);
    //collect_dups_and_push(temp_ia, temp_ja, temp_ar, ia, ja, ar, VERBOSE);
  }
  
  if (VERBOSE>1) {
    std::cout << "loaded rectange columns\n";
  }
  
  //loading the triangles shouldn't affect MAT_COMP
  offset = RE.size();
  for (i=0; i<(int)TR.size(); i++) {
    if (VERBOSE>3) {
      std::cout << "Triangle " << i << " " << TR[i] << "\n";
    }
    for (j=0; j<3; j++) {
      extract_signed_index(&sign, &index, (j==0 ? TR[i].a0 : (j==1 ? TR[i].a1 : TR[i].a2) ) ) ;
      LP.add_entry(index, offset+i, (sign<0 ? -1 : 1));
      //ia.push_back(index+1);
      //ja.push_back(offset+i+1);
      //if (sign < 0) {
      //  ar.push_back(-1);
      //} else {
      //  ar.push_back(1);
      //}
      if (VERBOSE>3) {
        //std::cout << "put " << ia[ia.size()-1] << " " << ja[ja.size()-1] << " " << ar[ar.size()-1] <<"\n";
      }
    }
  }
  if (VERBOSE>1) {
    std::cout << "loaded triangle columns\n";
  }  
  
  
  //if we are doing MAT_COMP, then load in the M and N matrix entries
  if (MAT_COMP) {
    if (MAT_SEPARATE_DOMAIN) {
      //load in the negative identity matrices
      row_offset = AL.num_arcs;
      col_offset = RE.size() + TR.size();
      for (i=0; i<2*WT.num_edges; i++) {
        LP.add_entry(row_offset+i, col_offset+i, -1); 
        //ia.push_back(row_offset + i + 1);
        //ja.push_back(col_offset + i + 1);
        //ar.push_back(-1);
      }
      //load in the N matrix (the b's are already set on the RHS)
      row_offset = AL.num_arcs + 2*WT.num_edges;
      col_offset = RE.size() + TR.size();
      for (i=0; i<WT.num_edges; i++) { //these are the columns (weights form the input x weight)
        for (j=0; j<b_len; j++) { //these are the ROWS, the output rows of the N matrix
          if (N[j][i] == 0) continue;
          LP.add_entry(row_offset+j, col_offset+i, N[j][i]);
          //ia.push_back(row_offset + j + 1);
          //ja.push_back(col_offset + i + 1);
          //ar.push_back(N[j][i]);
        }
      }
      //load in the matrix constraint
      col_offset = RE.size() + TR.size();
      row_offset = AL.num_arcs + 2*WT.num_edges + b_len;
      for (i=0; i<WT.num_edges; i++) { //these are the COLUMNS (the weights which form the input weight)
        for (j=0; j<WT.num_edges; j++) { //these are the ROWS  (the edge rows
          //for each column, compute how much it contributes to each row, 
          //that is, we are simply putting the ith column of M in
          if (M[j][i] == 0) continue;
          LP.add_entry(row_offset+j, col_offset+i, M[j][i]);
          //ia.push_back(row_offset + j + 1);
          //ja.push_back(col_offset + i + 1);
          //ar.push_back(M[j][i]);
        }
      }
      //we also need the negative identity on the right (for the side1)
      row_offset = AL.num_arcs+2*WT.num_edges + b_len;
      col_offset = RE.size() + TR.size() + WT.num_edges;
      for (i=0; i<WT.num_edges; i++) {
        LP.add_entry(row_offset+i, col_offset+i, -1);
        //ia.push_back(row_offset + i + 1);
        //ja.push_back(col_offset + i + 1);
        //ar.push_back(-1);
      }        
      
      if (VERBOSE>1) { 
        std::cout << "Loaded the M and N matrices\n";
      }
      
    } else {
      
      //load in the M matrix
      col_offset = RE.size() + TR.size();
      row_offset = AL.num_arcs;
      for (i=0; i<WT.num_edges; i++) { //these are the COLUMNS (the weights which form the input weight)
        for (j=0; j<WT.num_edges; j++) { //these are the ROWS  (the edge rows
          //for each column, compute how much it contributes to each row, 
          //and add these as negatives.
          //that is, we are simply putting the ith column of M in, with minus signs
          if (M[j][i] == 0) continue;
          LP.add_entry(row_offset+j, col_offset+i, -M[j][i]);
          //ia.push_back(row_offset + j + 1);
          //ja.push_back(col_offset + i + 1);
          //ar.push_back(-M[j][i]);
        }
      }
      
      //load in the N matrix (the b's are already set on the RHS)
      row_offset = AL.num_arcs + WT.num_edges;
      col_offset = RE.size() + TR.size();
      for (i=0; i<WT.num_edges; i++) { //these are the columns (weights form the input x weight)
        for (j=0; j<b_len; j++) { //these are the ROWS, the output rows of the N matrix
          if (N[j][i] == 0) continue;
          LP.add_entry(row_offset+j, col_offset+i, N[j][i]);
          //ia.push_back(row_offset + j + 1);
          //ja.push_back(col_offset + i + 1);
          //ar.push_back(N[j][i]);
        }
      }
      if (VERBOSE>1) { 
        std::cout << "Loaded the M and N matrices\n";
      }
      
      //load in the constraints that make sure the input is a weight
      row_offset = AL.num_arcs + WT.num_edges + b_len;
      col_offset = RE.size() + TR.size();
      for (i=0; i<WT.num_edges; i++) {  //these are the columns (input x weights)
        temp_ia.resize(0); temp_ja.resize(0); temp_ar.resize(0);
        row1 = WT.get_edge_source(i);
        temp_ia.push_back(row_offset + row1 + 1);
        temp_ja.push_back(col_offset + i + 1);
        temp_ar.push_back(1);
        row2 = WT.get_edge_dest(i);
        temp_ia.push_back(row_offset + row2 + 1);
        temp_ja.push_back(col_offset + i + 1);
        temp_ar.push_back(-1);
        LP.extend_entries_no_dups(temp_ia, temp_ja, temp_ar);
        //collect_dups_and_push(temp_ia, temp_ja, temp_ar, ia, ja, ar, VERBOSE);
      }
    }
  }
  
  
  if (DO_SUP) {
    offset = AL.num_arcs;
    for (i=0; i<WT.num_edges; i++) {
      LP.add_entry(offset+i, num_cols-1, -1);
      //ia.push_back(offset + i + 1);
      //ja.push_back(num_cols);
      //ar.push_back(-1);
    }
  }
  
  if (VERBOSE>1 && DO_SUP) {
    std::cout << "loaded supremum (t) column\n";
  }  
  
  if (VERBOSE > 1) {
    //std::cout << "Created " << ia.size() << " nonzeroes on " << num_rows << " rows and " << num_cols << " columns\n";
  }
  
  LP.solve(VERBOSE);
  
  LP.get_optimal_value(ans);
  LP.get_soln_vector(solution_vector);
  
  
}


int TROLLOP::trollop(int argc, char* argv[]) {
  int current_arg = 0;
  int i;
  int rank;
  int VERBOSE = 1;
  int LP_VERBOSE = 0;
  bool DO_SUP = false;
  bool USE_WORDS = false;
  bool OUTPUT_PROGRAM = false;
  bool MAT_COMP = false;
  bool MAT_SEPARATE_DOMAIN = false;
  SparseLPSolver solver = GLPK;
  int num_copies=1;
  
  std::string M_filename;
  std::string N_filename;
  std::string b_filename;
  std::string filename;
  
  if (argc < 2 || std::string(argv[1]) == "-h") {
    std::cout << "usage: ./scallop -train [-h] [-v[n]] [-V] [-m<GLPK,GIPT,GUROBI,EXLP>] [-L filename] [-sup,-scl, [-dom] -mat M_file N_file b_file rank length] [-w] <length> <chain or list of words>\n";
    std::cout << "\twhere <length gives the length of the words we want\n";
    std::cout << "\tand <chain...> is a chain OR, if -w, a list of words\n";
    std::cout << "\te.g. ./trollop 3 abABAbaB\n";
    std::cout << "\t-h: print this message\n";
    std::cout << "\t-v[n]: verbose output (n=0,1,2,3); 0 gives quiet output\n";
    std::cout << "\t-V: verbose LP output\n";
    std::cout << "\t-m<method>: use the LP solver specified\n";
    std::cout << "\t-L: output the linear program as a sparse matrix and two vectors to filename.A, .b, and .c\n";
    std::cout << "\t-sup: computes sup phi(C)/2D(phi) (this is the default behavior\n";
    std::cout << "\t-scl: computes scl(Psi_ell(C))\n";
    std::cout << "\t-w: use a list of words rather than a chain\n";
    std::cout << "\t-mat: minimize the scl over all weights of the form Mx, where Nx=b\n";
    std::cout << "\t-dom: separate the domain for the matrix computation (better lower bound)\n";
    exit(0);
  }
  while (argv[current_arg][0] == '-') {
    if (argv[current_arg][1] == 'm' && argv[current_arg][2] != 'a') {
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
    
    } else if (argv[current_arg][1] == 'v') {
      if (argv[current_arg][2] == '\0') {
        VERBOSE = 2;
      } else {
        VERBOSE = atoi(&argv[current_arg][2]);
      }
    
    } else if (argv[current_arg][1] == 'V') {
      LP_VERBOSE = 1;
      
    } else if (argv[current_arg][1] == 's') {
      if (argv[current_arg][2] == 'u') {
        DO_SUP = true;
      } else {
        DO_SUP = false;
      }
    } else if (argv[current_arg][1] == 'd') {
      MAT_SEPARATE_DOMAIN = true;
      num_copies=2;
      
    } else if (argv[current_arg][1] == 'm' && argv[current_arg][2] == 'a') {
      MAT_COMP = true;
      M_filename = std::string(argv[current_arg+1]);
      N_filename = std::string(argv[current_arg+2]);
      b_filename = std::string(argv[current_arg+3]);
      current_arg += 3;
    
    } else if (argv[current_arg][1] == 'w') {
      USE_WORDS = true;
      
    } else if (argv[current_arg][1] == 'L') {
      OUTPUT_PROGRAM = true;
      filename = std::string(argv[current_arg+1]);
      current_arg++;
    }
    
    current_arg++;
  }
  
  if (MAT_COMP) {
    rank = atoi(argv[current_arg]);
    current_arg++;
  } else {
    rank = chain_rank(argc-current_arg-1, &argv[current_arg+1]);
  }
  WordTable WT(rank, atoi(argv[current_arg]), DO_SUP || MAT_COMP);
  WVec C;
  std::vector<std::vector<int> > M;
  std::vector<std::vector<int> > N;
  std::vector<int> b;
  
  current_arg++;
  
  if (VERBOSE > 1) {
    std::cout << "Working word length: " << WT.ell << "\n";
    if (VERBOSE > 2) {
      WT.print();
    }
  }
  
  
  if (MAT_COMP) {
    //we still need to do this, just to create the homology info
    WT.create_index_assignments(C);
    //load the matrices
    read_matrix(M, M_filename, WT, true, VERBOSE);
    read_matrix(N, N_filename, WT, false, VERBOSE);
    read_vector(b, b_filename, VERBOSE);
    
  } else {
    //if we're computing null scl, then we don't need all this
    
    C =  WVec(WT, &argv[current_arg], argc-current_arg, USE_WORDS);   //process the chain argument
    
    //exit(0);
    
    if (VERBOSE>1) {
      std::cout << "Input vector: " << C << "\n";
      std::cout.flush();
    }
    
    //exit(0);
    
    WT.create_index_assignments(C); //assign working indices to words (edges and vertices)
    
    if (VERBOSE>1) {
      std::cout << "Made index assignments:\n";
      WT.print();
    }  
    
    //exit(0);
    
    C.fill_index_coefficients(WT);   //create the vector of edges from C
    
    if (VERBOSE>2) {
      std::cout << "Got the vector of edges from the input vector:\n";
      for (i=0; i<WT.num_edges; i++) {
        std::cout << C.index_coefficients[i] << " ";
      }
      std::cout << "\n";
    }
  }
  
  //exit(0);
  
  ArcPairList AL(WT, num_copies);
  if (VERBOSE > 1) {
    std::cout << "computed arcs (" << AL.num_arcs << ")\n";
    if (VERBOSE > 2) {
      AL.print(std::cout);
    }
  }
  
  //exit(0);
  
  std::vector<Triangle> TR(0);
  compute_triangles(WT, AL, num_copies, TR);
  if (VERBOSE > 1) {
    std::cout << "computed triangles (" << TR.size() << ")\n"; std::cout.flush();
    if (VERBOSE > 2) {
      print_triangles(TR, std::cout);
    }
  }
  
  //exit(0);
  
  std::vector<Rectangle> RE;
  compute_rectangles(WT, AL, num_copies, RE);
  if (VERBOSE > 1) {
    std::cout << "computed rectangles (" << RE.size() << ")\n"; std::cout.flush();
    if (VERBOSE > 2) {
      print_rectangles(RE, std::cout);
    }
  }
  
  //exit(0);
  
  
  if (OUTPUT_PROGRAM) {
    write_lp(WT, C, AL, num_copies, TR, RE, DO_SUP, MAT_COMP, MAT_SEPARATE_DOMAIN, M, N, b, VERBOSE, LP_VERBOSE, filename);
    if (VERBOSE>0) {
      std::cout << "Wrote linear program\n";
    }
  } else {

    Rational ans;
    std::vector<Rational> solution_vector(0);                           //run the LP
    
    trollop_lp(WT, C, AL, num_copies, TR, RE, DO_SUP,
               MAT_COMP,
               MAT_SEPARATE_DOMAIN,
               M, N, b,
              ans, 
              solution_vector, 
              solver,
              VERBOSE,
              LP_VERBOSE); 
    
    if (VERBOSE>0) {
        //output for supremum
      if (DO_SUP) {
        std::cout << "sup_{Q_" << WT.ell << "} phi(C)/2D(phi) = " 
                  << "inf_t(w + tE) = (t->" <<  solution_vector[solution_vector.size()-1] 
                  << "); " << ans << " = " << ans.get_d() << "\n";
        if (VERBOSE > 1) {
          std::cout << "t = " << solution_vector[solution_vector.size()-1] << "\n"; 
        }
        
        //output for MAT_COMP
      } else if (MAT_COMP) {
        if (MAT_SEPARATE_DOMAIN) {
          std::vector<std::vector<Rational> > weight_vector;
          weight_vector.resize(num_copies);
          for (i=0; i<num_copies; i++) {
            weight_vector[i].resize(WT.num_edges);
          }
          for (i=0; i<WT.num_edges; i++) {
            weight_vector[0][i] = solution_vector[RE.size() + TR.size() + i];
            weight_vector[1][i] = solution_vector[RE.size() + TR.size() + WT.num_edges + i];
          }          
          std::cout << "Min side 0 (input weight x): ";
          for (i=0; i<(int)weight_vector[0].size(); i++) {
            std::string temp;
            if (weight_vector[0][i] > 0) {
              WT.get_word(temp, WT.ell, i);
              std::cout << weight_vector[0][i] << temp << " ";
            }
          }
          std::cout << "\n" << "Min side 1 (phi(x)): ";
          for (i=0; i<(int)weight_vector[1].size(); i++) {
            std::string temp;
            if (weight_vector[1][i] > 0) {
              WT.get_word(temp, WT.ell, i);
              std::cout << weight_vector[1][i] << temp << " ";
            }
          }
          std::cout << "\nwith scl = " << ans << " = " << ans.get_d() << "\n";
          
        } else {
          std::vector<Rational> weight_vector(WT.num_edges);
          std::vector<Rational> input_weight_vector(WT.num_edges);
          for (i=0; i<(int)weight_vector.size(); i++) {
            weight_vector[i] = Rational(0,1);
            input_weight_vector[i] = Rational(0,1);
          }
          for (i=0; i<(int)RE.size(); i++) {
            if (solution_vector[i] > 0) {
              weight_vector[RE[i].e0] = weight_vector[RE[i].e0] + solution_vector[i];
              weight_vector[RE[i].e1] = weight_vector[RE[i].e1] + solution_vector[i];
            }
          }
          for (i=0; i<WT.num_edges; i++) {
            input_weight_vector[i] = solution_vector[RE.size() + TR.size() + i];
          }
          std::cout << "Min weight (i.e. partial S): ";
          for (i=0; i<(int)weight_vector.size(); i++) {
            std::string temp;
            if (weight_vector[i] > 0) {
              WT.get_word(temp, WT.ell, i);
              std::cout << weight_vector[i] << temp << " ";
            }
          }
          std::cout << "\nwith scl = " << ans << " = " << ans.get_d() << "\n";
          std::cout << "Min input weight (the x in partial S = Mx): ";
          for (i=0; i<(int)input_weight_vector.size(); i++) {
            std::string temp;
            if (input_weight_vector[i] > 0) {
              WT.get_word(temp, WT.ell, i);
              std::cout << input_weight_vector[i] << temp << " ";
            }
          }
          std::cout << "\n";
        }
        //output for scl 
      } else {
        std::cout << "scl( ";
        C.print_words(std::cout);
        std::cout << ") = " << ans << " = " << ans.get_d() << "\n";    //output the answer
      }
      
      if (VERBOSE > 1){// && (int)solution_vector.size() < 1000) {
        for (i=0; i<(int)solution_vector.size()-1; i++) {
          if (solution_vector[i] > 0 && i < (int)RE.size() + (int)TR.size()) { 
            if (i < (int)RE.size()) {
              std::cout << solution_vector[i] << " * " << RE[i] << "\n";
            } else {
              std::cout << solution_vector[i] << " * " << TR[i-(int)RE.size()] << "\n";
            }
          }
        }
      }
    } else {
      std::cout << ans.get_d() << "\n";
    }
    
  }
    
  return 0;
}
