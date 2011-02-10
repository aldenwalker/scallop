

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "scallop.h"
#include "io.h"

using namespace std;

void write_linear_program(vector<string>& w, 
                          vector<int>& weight, 
                          vector<arc>& arc_list, 
                          vector<polygon>& polygon_list, 
                          string& programFile) {
  int i,j,k;
  int coef;
  string AFile = programFile + ".A";
  string bFile = programFile + ".b";
  string cFile = programFile + ".c";
  fstream outAFile;
  fstream outbFile;
  fstream outcFile;
  
  /////  A matrix
  
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
