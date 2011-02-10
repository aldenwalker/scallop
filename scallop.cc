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
*		Alden Walker 2009, 2010        	*
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
                       
	int i,j,k,size,add,close,appeared;
	polygon testpoly;
	testpoly.arc.resize(maxjun);
	int arc_list_length = arc_list.size();

  vector<int> word_length(w.size());
  
  for (i=0; i<(int)w.size(); i++) {
    word_length[i] = w[i].size();
  }

	polygon_list.resize(0);

	for(i=0;i<arc_list_length-1;i++){	// i is initial index of hypothetical polygon
		testpoly.arc[0]=i;
		j=(i+1);
		size=1;
		while(size>0){

			add=1;		// haven't yet added
			close=1;	// haven't yet closed
			
			// does arc j glue up to end of last arc of testpoly?
			if(arc_list[j].first_word == arc_list[testpoly.arc[size-1]].last_word && 
			  (arc_list[j].first - arc_list[testpoly.arc[size-1]].last-1)%word_length[arc_list[j].first_word]==0){		
				appeared=1;
				for(k=0;k<size;k++){
					if(j==testpoly.arc[k]){		//has arc j already appeared in testpoly?
						appeared=0;
					}
				}
			
				if(appeared==1){	        // if j has not appeared
					testpoly.arc[size]=j;		// then add it to testpoly
					add=0;						      // note that we have added an arc
					
					// does it close up?
					if(arc_list[j].last_word == arc_list[testpoly.arc[0]].first_word && 
					  (arc_list[testpoly.arc[0]].first-arc_list[j].last-1)%word_length[arc_list[j].last_word]==0){		
						testpoly.size = size+1;
						polygon_list.push_back(testpoly);
						close=0;				// note that we have closed a polygon
					}
				}
			}
	
	
			if(add==0 && close==1){
				size++;
				j=i+1;
				if(size>=maxjun){
					j=arc_list_length;
				}
			} else {
				j++;		// increment j
			}
			while(j>=arc_list_length){
				size--;
				j=testpoly.arc[size]+1;
			}
		}
	}

}


void generate_polygons_new(vector<string> w,
                          vector<polygon> &polygon_list,  
                          vector<arc>  &arc_list, 
                          int maxjun) {
  int i,j,k;
  int ALL = (int)arc_list.size();
  int len_bound = (maxjun%2 == 0 ? maxjun/2 : maxjun/2 + 1);
  int total_len_bound = maxjun;
  int old_old_arc_seq_len;
  int old_arc_seq_len;
  int current_poly_len;
  arc temp_arc;
  arc temp_arc_2;
  polygon_list.resize(0);
  
  vector<int> templist(1);
  
  vector<vector<int> > arc_sequences(0);
  //build all the arc sequences of up to half the necessary length
  //hopefully, there won't be a ridiculous number of these
  
  //start with length 1 (just all arcs)
  for (i=0; i<ALL; i++) {
    templist[0] = i;
    arc_sequences.push_back(templist);
  }
  
  //now build up
  old_old_arc_seq_len = 0;
  for (j=0; j<len_bound-1; j++) {     //we need to add on len_bound-1 arcs
    old_arc_seq_len = arc_sequences.size();
    current_poly_len = arc_sequences[old_old_arc_seq_len].size();
    for (i=old_old_arc_seq_len; i<old_arc_seq_len; i++) {     //for each arc sequence
      for (k=0; k<ALL; k++) {                               //for each arc
        temp_arc = arc_list[arc_sequences[i][current_poly_len-1]];
        if (temp_arc.last_word == arc_list[k].first_word
            &&
            (temp_arc.last + 1)%w[temp_arc.last_word].size() 
                              == arc_list[k].first) {
        //so this arc can be glued to the end of this arc sequence
        //there are two options -- either the polygon (sequence) closes up, in 
        //which case we put it on polygon_list and *not* into the larger 
        //arc sequence list
        temp_arc_2 = arc_list[arc_sequences[i][0]];
        if (arc_list[k].last_word == temp_arc_2.first_word
        ) {
          
        } else {  //or it doesn't close up, so we just append this to the
                  //list of arc sequences
        
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
	string programFile = "";
	int RAT = 0;
	int overrideMaxjun = 0;
	int overridedMaxjun = -1;
	string drawFile = "";

	if(argc < 2 || strcmp(argv[1],"-h")==0){
	  cout << "scallop\n";
	  cout << "version 2.11 - December 10, 2010\n";
	  cout << "by Danny Calegari and Alden Walker\n";
	  cout << "see the README for details\n";
		cout << "usage: scallop [-revh, -mn, -s filename, -L[!] filename] [i_1]w_1 [i_1]w_2 . . [i_n]w_n for scl(i_1*w_1 + i_2*w_2 + . . + i_n*w_n) \n";
		cout << "options:\n\t-s filename :\tdraw the components of an extremal surface, each to \n";
		cout <<                      "\t\t\ta different file and print generators for the \n";
		cout <<                      "\t\t\timage of the fundamental group\n";
		cout << "\t-r : use rational arithmetic (GMP) internally\n";
		cout << "\t-e : rigorous (slower) calculation (see README)\n";
		cout << "\t-L[!] : (!=no solving) output the linear program to filename.A, .b, and .c\n";
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


  //NOTE: it would be really awesome if we knew how big to make the polygon
  //list, since right now it probably does tons of single-element reallocs, 
  //which is just horrible

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
	generate_polygons(w, polygon_list, arc_list, maxjun);
	polygon_list_length = polygon_list.size();
	
	if(VERBOSE==1){
		cout << polygon_list_length << " polygons (cyclic list of arcs) \n";
		for(i=0;i<polygon_list_length;i++){
			cout << "polygon " << i << " : ";
			for(j=0;j<polygon_list[i].size;j++){
				cout << polygon_list[i].arc[j] << " ";
			};
			cout << '\n';
		};
	};

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
