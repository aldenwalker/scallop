

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "scallop.h"
#include "rational.h"
#include "io.h"

using namespace std;

void write_linear_program(vector<string>& w, 
                          vector<int>& weight, 
                          vector<arc>& arc_list, 
                          vector<polygon>& polygon_list, 
                          string& programFile) {
  int i,j;
  string AFile = programFile + ".A";
  string bFile = programFile + ".b";
  string cFile = programFile + ".c";
  fstream outAFile;
  fstream outbFile;
  fstream outcFile;
  
  /////  A matrix
  /*
  outAFile.open(AFile.c_str(), fstream::out);
  //write out the arc constraints matrix
  outAFile << (int)arc_list.size()/2 + (int)w.size() << " " << (int)polygon_list.size() << " 0\n";
  for (i=0; i<(int)arc_list.size()/2; i++) {          //for each arc (row)
    for (j=0; j<(int)polygon_list.size(); j++) {     //for each polygon (column)
      coef = 0;
      for (k=0; k<(int)polygon_list[j].size; k++) {  //for each arc in the poly
        if ( polygon_list[j].arc[k] == 2*i ) {         //check if it's our arc
          coef++;
        } else if (polygon_list[j].arc[k] == 2*i+1) {   //it's the inverse
          coef--;
        }
      }
      //outAFile << coef << " ";                           //write the coefficient
      if (coef != 0) {
        outAFile << i+1 << " " << j+1 << " " << coef << "\n";   //sparse output
      }
    }
    //outAFile << "\n";                              //this ends the row
  }
  
  //write out the chain constraint
  for (i=0; i<(int)w.size(); i++) {               //for each word (row here)
    for (j=0; j<(int)polygon_list.size(); j++) {  //for each polygon (column)
      coef = 0;
      for (k=0; k<(int)polygon_list[j].size; k++) {
        if (arc_list[polygon_list[j].arc[k]].first == 0 &&
            arc_list[polygon_list[j].arc[k]].first_word == i) {
          coef++;
        }
      }
      //outAFile << coef << " ";
      if (coef != 0) {
        outAFile << (int)arc_list.size()/2 + i + 1 << " " << j+1 << " " << coef << "\n";
      }
    }
    //outAFile << "\n";
  }
  
  outAFile.close();
  
  */
  /////  A matrix
  
  outAFile.open(AFile.c_str(), fstream::out);
  //write out the dimensions
  outAFile << (int)arc_list.size()/2 + (int)w.size() << " " << (int)polygon_list.size() << " 0\n";
  //we need to do this intelligently, so we do each column (polygon)
  //every arc can appear at most once, so we just read off the arcs
  for (i=0; i<(int)polygon_list.size(); i++) {
    for (j=0; j<(int)polygon_list[i].size; j++) {
      if (polygon_list[i].arc[j] % 2 == 0) {
        outAFile << (polygon_list[i].arc[j]/2)+1 << " " << i+1 << " 1\n";
      } else {
        outAFile << ((polygon_list[i].arc[j]-1)/2)+1 << " " << i+1 << " -1\n";
      }
    }
  }
  
  //now write the chain constraint
  for (i=0; i<(int)polygon_list.size(); i++) {
    for (j=0; j<(int)polygon_list[i].size; j++) {
      if (arc_list[polygon_list[i].arc[j]].first == 0) {
        outAFile << (int)arc_list.size()/2 + arc_list[polygon_list[i].arc[j]].first_word + 1
                 << " " << i+1 << " 1\n";
      }
    }
  }
  
  outAFile.close(); 
  
  
  
  /////// b vector
  outbFile.open(bFile.c_str(), fstream::out);
  
  for (i=0; i<(int)arc_list.size()/2; i++) {          //all the arc constraints are zero
    outbFile << "0\n";
  }
  
  for (i=0; i<(int)w.size(); i++) {                 //the word RHS is the weight
    outbFile << weight[i] << "\n";
  }
  
  outbFile.close();
  
  /////// c objective function
  outcFile.open(cFile.c_str(), fstream::out);
  
  for (i=0; i<(int)polygon_list.size(); i++) {   //chi is just size -2
    outcFile << polygon_list[i].size-2 << "\n";
  }
  outcFile.close();
  
}


/*****************************************************************************/
/* this function writes out the polyhedron of all (optimal) solutions, by    */
/* simply adding the constraint that the objective function is 4*scl         */
/* perhaps in the future this will actually do the computation?  right now   */
/* it only outputs the description for cddlib and/or lrslib                  */
/*****************************************************************************/
void write_solution_polyhedron(vector<string>& w, 
                               vector<int>& weight, 
                               vector<arc>& arc_list, 
                               vector<polygon>& polygon_list, 
                               string& polyFileName,
                               rational scl) {
  int i,j,k, coef;
  int numLins = (int)arc_list.size()/2 + (int)w.size() + 1;  //+1 is optimality
  int dim = (int)polygon_list.size();
  int rows = numLins + dim; //lins, plus >=0
  fstream outPolyFile;
  outPolyFile.open(polyFileName.c_str(), fstream::out);
  
  //we can't output a sparse matrix or anything like that
  outPolyFile << "H-representation\n";
  outPolyFile << "linearity " << numLins;
  for (i=0; i<numLins; i++) {
    outPolyFile << " " << i+1;
  }
  outPolyFile << "\n";
  outPolyFile << "begin\n";
  outPolyFile << rows << " " << dim+1 << " rational\n";
  
  
  //write the linearities:
  //the arc constraints:
  for (i=0; i<(int)arc_list.size()/2; i++) {          //for each arc (row)
    outPolyFile << "0";                               //the RHS is zero
    for (j=0; j<(int)polygon_list.size(); j++) {     //for each polygon (column)
      coef = 0;
      for (k=0; k<(int)polygon_list[j].size; k++) {  //for each arc in the poly
        if ( polygon_list[j].arc[k] == 2*i ) {         //check if it's our arc
          coef++;
        } else if (polygon_list[j].arc[k] == 2*i+1) {   //it's the inverse
          coef--;
        }
      }
      outPolyFile << " " << coef;                           //write the coefficient
    }
    outPolyFile << "\n";                              //this ends the row
  }
  
  //the chain constraints (word i appears weight[i] times)
  for (i=0; i<(int)w.size(); i++) {               //for each word (row here)
    outPolyFile << -weight[i];             //the RHS is the weight
                                           //(but negative because of the format)
    for (j=0; j<(int)polygon_list.size(); j++) {  //for each polygon (column)
      coef = 0;
      for (k=0; k<(int)polygon_list[j].size; k++) {
        if (arc_list[polygon_list[j].arc[k]].first == 0 &&
            arc_list[polygon_list[j].arc[k]].first_word == i) {
          coef++;
        }
      }
      outPolyFile << " "  << coef;
    }
    outPolyFile << "\n";
  }
  
  //the constraint that says we have an optimal solution
  //(maybe this should be first?)
  outPolyFile << rational(-4,1)*scl;                          //the RHS is 4*scl
                                                  //(negative because of format)
  for (i=0; i<(int)polygon_list.size(); i++) {
    outPolyFile << " " << polygon_list[i].size-2;         //-chi is the sum of these
  }
  outPolyFile << "\n";
  
  
  //the inequalities:
  //these just require that the answer is nonnegative
  for (i=0; i<(int)polygon_list.size(); i++) {
    outPolyFile << "0";
    for (j=0; j<(int)polygon_list.size(); j++) {
      outPolyFile << " " << (i==j ? 1 : 0);
    }
    outPolyFile << "\n";
  }
  
  outPolyFile << "end\n";
  
  outPolyFile.close();
  
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  















