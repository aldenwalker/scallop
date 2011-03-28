/************************************
*                                   *
*	scallop.cc             	        	*
*								                  	*	
*	calculates scl on finite	      	*
*	linear sums of elements in	    	*
*	free groups					            	*
*								                  	*
*	Copyright Danny Calegari 2008	    *
*   Copyright Danny Calegari and	  *
*		Alden Walker 2009, 2010, 2011  	*
*									                  *
*	Includes modifications by	      	*
*		Alden Walker				            *
*									                  *
*	Released under the GPL license	  *
*									                  *
************************************/
	

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <ctype.h>

#include "word.h"
#include "scallop.h"
#include "draw.h"
#include "rational.h"
#include "lp.h"
#include "io.h"

using namespace std;



/******************************************************************************/
/*  generate a list of all possible arcs                                      */
/******************************************************************************/
void generate_arcs(vector<arc>* arc_list, vector<string>& w, int WORD) {
	int g,h,i,j;
	int a,b;
	arc tempArc;
	arc_list->resize(0);
	for(g=0;g<WORD;g++){	
		for(i=0;i< (int) w[g].length();i++){
			a=(int) w[g][i];		
			for(h=g;h<WORD;h++){
				if(h==g){
					for(j=i+1;j< (int) w[h].length();j++){
						b=(int) w[h][j];
						if((32+a-b)%64==0){		// are letters inverse?
						  tempArc.first = i;
						  tempArc.last = j;
						  tempArc.first_word=g;
						  tempArc.last_word=h;
						  arc_list->push_back(tempArc);
						  tempArc.last = i;
						  tempArc.first = j;
						  tempArc.first_word = h;
						  tempArc.last_word = g;
						  arc_list->push_back(tempArc);
						};		
					};
				} else {
					for(j=0;j< (int) w[h].length();j++){
						b=(int) w[h][j];
						if((32+a-b)%64==0){		// are letters inverse?
						  tempArc.first = i;
						  tempArc.last = j;
						  tempArc.first_word=g;
						  tempArc.last_word=h;
						  arc_list->push_back(tempArc);
						  tempArc.last = i;
						  tempArc.first = j;
						  tempArc.first_word = h;
						  tempArc.last_word = g;
						  arc_list->push_back(tempArc);
						};							
					};
				};
			};			
		};
	};

};


/******************************************************************************/
/* generate the polygon list                                                  */
/******************************************************************************/
void generate_polygons(vector<string> w,
                       vector<polygon> &polygon_list,  
                       vector<arc>  &arc_list, 
                       int maxjun) {
  int i,j;
  
  polygon_list.resize(0);

  //we need to make a structure saying which arcs start where
  //the array structure is arc_locs[word][index] is a list
  //of the arcs which start at index index in word word
  vector<vector<vector<int> > > arc_locs(0);
  arc_locs.resize(w.size());

  for (i=0; i<(int)w.size(); i++) {
    arc_locs[i].resize(w[i].size());
    for (j=0; j<(int)w[i].size(); j++) {
      arc_locs[i][j].resize(0);
    }
  }

  //go through the arcs and add them all
  for(i=0; i<(int)arc_list.size(); i++) {
    arc_locs[arc_list[i].first_word][arc_list[i].first].push_back(i);
  }

  //now we build the polygon
  vector<int> temp_arc_list(maxjun);
  int temp_arc_list_size = 0; // this is probably redundant
  vector<int> arc_end_indices(maxjun); //this keeps track of where the current
                                       //arc list selection has ending indices
                                       //in the data structure
  polygon temp_poly;
  temp_poly.arc = vector<int>(0);
  int next_arc_ind;
  arc first_arc;
  arc cur_arc;
  int cur_target_word;
  int cur_target_ind;
  int we_pushed;
  for (i=0; i<(int)arc_list.size(); i++) {  //loop through which arc is the first
                                            //note we may assume it is minimal
    temp_arc_list[0] = i;
    temp_arc_list_size = 1;
    arc_end_indices[0] = 0;
    first_arc = arc_list[i];
    while (1) {
      if (temp_arc_list_size==0) {  //we are really done
        break;
      }
      //first, load the arc that we're talking about
      cur_arc = arc_list[temp_arc_list[temp_arc_list_size-1]];
      cur_target_word = cur_arc.last_word;
      cur_target_ind = (cur_arc.last+1)%(int)w[cur_target_word].size();

      //if we have reached this point, we know that 
      //(A) our arc doesn't glue to the original arc and
      //(B) our arc has index bigger than the original arc
      //therefore, all we need to do is to find the next
      //ok arc to push on, and then continue
      //the next candidate is contained in arc_locs at index
      //arc_end_indices[temp_arc_list_size-1]
      we_pushed = 0;
      while (arc_end_indices[temp_arc_list_size-1] < 
             (int)arc_locs[cur_target_word][cur_target_ind].size()) {
        next_arc_ind = arc_locs[cur_target_word]
                               [cur_target_ind]
                               [arc_end_indices[temp_arc_list_size-1]];
        //check if the next arc glues up or has a bad index or we've got maxjun
        if (next_arc_ind < i) {
          //index is too small, skip it
        
        } else if (arc_list[next_arc_ind].last_word == first_arc.first_word
             && (arc_list[next_arc_ind].last+1)%(int)w[first_arc.first_word].size() == first_arc.first) {
          //it glues!
          temp_poly.arc.resize(temp_arc_list_size+1);
          temp_poly.size = temp_arc_list_size+1; 
          for (j=0; j<temp_arc_list_size; j++) {
            temp_poly.arc[j] = temp_arc_list[j];
          }
          temp_poly.arc[j] = next_arc_ind;
          polygon_list.push_back(temp_poly);
        } else if (temp_arc_list_size == maxjun-1) {
          //we've reached maxjun -- can't push anything

        } else {
          //it doesn't glue, its index is ok, and we haven't reached maxjun,
          //then we add it on, step forward, and continue
          temp_arc_list[temp_arc_list_size] = next_arc_ind;
          arc_end_indices[temp_arc_list_size] = 0;
          temp_arc_list_size++;
          we_pushed = 1;
          break;
        }
        arc_end_indices[temp_arc_list_size-1]++;
      }
      
      //if we pushed a new arc, ok just loop
      //if we didn't push a new arc, we need to step back
      if (we_pushed==1) {
        //nothing
      } else {
        temp_arc_list_size--;
        if (temp_arc_list_size > 0) {
          arc_end_indices[temp_arc_list_size-1]++;
        }
      }
    }
  }	
  
}










int main(int argc, char* argv[]){


	int i,j,k;
	int chainStart = 1;
	int WORD;
	int VERBOSE = 0;
	int DRAW = 0;
	int RIGOROUS = 0;
	int WRITE_PROGRAM = 0;
	int ONLY_WRITE_PROGRAM = 0;
  int WRITE_POLYHEDRON = 0;
  string polyFileName = "";
	string programFile = "";
	int RAT = 0;
	int overrideMaxjun = 0;
	int overridedMaxjun = -1;
	string drawFile = "";

	if(argc < 2 || strcmp(argv[1],"-h")==0){
	  cout << "scallop\n";
	  cout << "version 2.6 - March 28, 2011\n";
	  cout << "by Danny Calegari and Alden Walker\n";
	  cout << "see the README for details\n";
		cout << "usage: scallop [-revh, -mn, -s filename, -L[!] filename, -P filename] [i_1]w_1 [i_1]w_2 . . [i_n]w_n for scl(i_1*w_1 + i_2*w_2 + . . + i_n*w_n) \n";
		cout << "options:\n\t-s filename :\tdraw the components of an extremal surface, each to \n";
		cout <<                      "\t\t\ta different file and print generators for the \n";
		cout <<                      "\t\t\timage of the fundamental group\n";
		cout << "\t-r : use rational arithmetic (GMP) internally\n";
		cout << "\t-e : rigorous (slower) calculation (see README)\n";
		cout << "\t-L[!] : (!=no solving) output the linear program to filename.A, .b, and .c\n";
		cout << "\t-P : output a polyhedron representation of the optimal solutions to filename\n";
    cout << "\t-mn : advanced: use polygons with up to n edges\n";
		cout << "\t-v : verbose output\n";
		cout << "\t-h : print this message\n";
		return(0);
	};
	
	//get the options
	WORD=argc-1;
  chainStart = 1;
  while (argv[chainStart][0] == '-') {
    WORD--;
    switch (argv[chainStart][1]) {
      case 'v':
        VERBOSE = 1;
        break;
      case 's':
        DRAW = 1;
        drawFile = argv[chainStart+1];
        chainStart++;
        WORD--;
        break;
      case 'L':
        WRITE_PROGRAM = 1;
        if (argv[chainStart][2] == '!') {
          ONLY_WRITE_PROGRAM = 1;
         }
        programFile = argv[chainStart+1];
        chainStart++;
        WORD--;
        break;
      case 'P':
        WRITE_POLYHEDRON = 1;
        polyFileName = argv[chainStart+1];
        chainStart++;
        WORD--;
        break;
      case 'r':
        RAT = 1;
        break;
      case 'e':
        RIGOROUS = 1;
        break;
      case 'm':
        overrideMaxjun = 1;
        overridedMaxjun = atoi(&(argv[chainStart][2]));
        break;
      default:
        cout << "Unexpected option\n";
        exit(1);
        break;
    }
    chainStart++;
  }

  //chainStart is now the index of the first word in the chain

	vector<string> w(WORD);
	vector<string> oldW(WORD);
	int arc_list_length;
	int polygon_list_length;
	vector<arc> arc_list(0);
	vector<polygon> polygon_list(0);
	vector<char> lettersSeen(0);
	int maxjun;
	char c;
	rational scl;
	vector<int> weight(WORD);
	vector<rational> solutionVector(0);

	

	//cout << "starting...\n";
	//fflush(stdout);

	//preprocess the words to remove the weights
	int weightLen;
	for(i=0;i<WORD;i++){		// or should that be WORDS?
		w[i]=argv[chainStart + i];
		weightLen = 0;
		j=0;
		while (isdigit(w[i][j])) {
		  j++;
		}
		if (j>0) {
		  weight[i] = atoi(w[i].substr(0,j).c_str());
		} else {
		  weight[i] = 1;
		}
		w[i] = w[i].substr(j,w[i].length()-j);
	};
	
	
	//find the rank (count the letters)
	for (i=0; i<WORD; i++) {
	  for (j=0; j<(int)w[i].size(); j++) {
	    c = tolower(w[i][j]);
	    for (k=0; k<(int)lettersSeen.size(); k++) { 
	      if (lettersSeen[k] == c) {
	        break;
	      }
	    }
	    if (k==(int)lettersSeen.size()) { //did we not break?
	      lettersSeen.push_back(c);
	    }
	  }
	}
	
	//possibly override the number of polygon sides
	if (overrideMaxjun == 0) {
	  maxjun = 2*lettersSeen.size(); //i.e. maxjun is twice the rank
	} else {
	  maxjun = overridedMaxjun;
	  cout << "Overriding max polygon sides to " << maxjun << "\n";
	}
	
	//make sure the words are cyclically reduced
  for (i=0; i<(int)w.size(); i++) {
    cyc_red(w[i]);
  }
  if (VERBOSE==1) {
    cout << "Words after cyclic reduction:\n";
    for (i=0; i<(int)w.size(); i++) {
      cout << w[i] << "\n";
    }
  }
  
  //if we're doing a rigorous calculation, triple every letter
  if (RIGOROUS == 1) {
    vector<string> wTripled(WORD);
    for (i=0; i<WORD; i++) {
      wTripled[i] = "";
      for (j=0; j<(int)w[i].size(); j++) {
        wTripled[i] += w[i][j];
        wTripled[i] += w[i][j];
        wTripled[i] += w[i][j];
      }
      oldW[i] = w[i];
      w[i] = wTripled[i];
    }
    if (VERBOSE == 1) {
      cout << "Words tripled for rigorous calculation:\n";
      for (i=0; i<WORD; i++) {
        cout << w[i] << "\n";
      }
    }
  }
        
  



  //generate the arcs!
	generate_arcs(&arc_list, w, WORD);
	arc_list_length = arc_list.size();
	
	if(VERBOSE==1){
		cout << arc_list_length << " arcs (start letter, end letter, start word, end word) \n";
		for(i=0;i<arc_list_length;i++){
			cout << "arc " << i << " : ";
			cout << arc_list[i].first << " " << arc_list[i].last << " " << arc_list[i].first_word << " " << arc_list[i].last_word << "\n";
	
		};
	};
	
	//generate the polygons!
	//generate_polygons(w, polygon_list, arc_list, maxjun);
	generate_polygons(w, polygon_list, arc_list, maxjun);
	if (VERBOSE==1) {
	  cout << "generated polys\n";
	}

	polygon_list_length = polygon_list.size();
	
	if(VERBOSE==1){
		cout << polygon_list_length << " polygons (cyclic list of arcs) \n";
		for(i=0;i<polygon_list_length;i++){
			cout << "polygon " << i << " : ";
			for(j=0;j<polygon_list[i].size;j++){
				cout << polygon_list[i].arc[j] << " ";
			};
			cout << '\n';
		}
	}

	if (WRITE_PROGRAM == 1) {
	  write_linear_program(w, weight, arc_list, polygon_list, programFile);
	  if (ONLY_WRITE_PROGRAM == 1) {
	    return 0;
	  }
  }

	solutionVector.resize(polygon_list_length);
	
	//Run the linear program!
  do_linear_program(w, weight, arc_list, polygon_list, scl, 
                                              solutionVector, 
                                              (RAT == 0 ? GLPK_DOUBLE : EXLP),
                                              VERBOSE);	
	
  if(VERBOSE==1){
		cout << "optimal solution vector (weight by polygon type): \n";
		for(i=0;i<polygon_list_length;i++){
			cout << "polygon " << i << " : " << solutionVector[i] << "\n";
		}
	}
  
  //write the optimal polyhedron if we're supposed to 
  if (WRITE_POLYHEDRON == 1) {
    write_solution_polyhedron(w, weight, arc_list, polygon_list, polyFileName, scl);
  }
  
	
	//draw the surface if we're supposed to
	if (DRAW == 1) {
	  make_graph_and_print(drawFile, w, arc_list, polygon_list, solutionVector);
	}

  //output the answer
  if (RIGOROUS==1) { //bring back the original words
    for (i=0; i<WORD; i++) {
      w[i] = oldW[i];
    }
  }  
  ostringstream weightString;
	cout << "scl ( ";
	for(i=0;i<WORD;i++){
	  weightString.str("");
	  if (weight[i] != 1) {
	    weightString << weight[i];
	  }
		if(i==0){
			cout << weightString.str() << w[0] << " ";
		} else {
			cout << "+ " << weightString.str() << w[i] << " "; 
		};
	};
	cout << ") = ";
	cout << scl << " = "  << scl.get_d() << "\n";
	
	//The vector<>'s should get deconstructed when their scope ends?
	
	return 0;
}
