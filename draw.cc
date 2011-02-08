#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include <cmath>

#include "word.h"
#include "rational.h"
#include "draw.h"


using namespace std;
  

struct vert {
  int polyRep;
  vector<int> e;  //in order!
  vector<int> end;   //1 means it's starting here, -1 that it's leaving
};

struct edge {
  string label;
  int start;
  int end;
};



void print_tree_ps(ostream& OS, vector<vert>& verts, vector<edge>& edges, vector<int>& shortEdges, 
                                                               int root, 
                                                               int depth,
                                                               int rootX,
                                                               int rootY,
                                                               int xwindow,
                                                               int levelIncrement) {
  int i;
  int xmin = rootX - (int)( 0.5 * xwindow);
  int xmax = rootX + (int)( 0.5 * xwindow);
  //put the root text
  //cout << "/Arial findfont\n";
  //cout << "20 scalefont setfont\n";
  //cout << rootX+10 << " " << rootY-10 << " moveto\n";
  //cout << "(" << root << ")" << "show\n";
  OS << "5 setlinewidth\n";
  OS << rootX << " " << rootY << " moveto\n";
  OS << rootX << " " << rootY << " 2  0 360 arc\n"; 
  OS << "stroke\n";
  OS << "2 setlinewidth\n";
  
  //draw the lines to the children
  int numChildren = verts[root].e.size();
  int destinationX;
  int destinationY;
  //int destinationRoot;
  string tempLabel;
  //the first one goes to xmin + (1/2numChildren)(xmax-xmin), then we add (1/numChildren)(xmax-xmin)
  destinationX = xmin + (int)( (1.0/(2.0*numChildren)) * (xmax-xmin) );
  destinationX -= (int)( (1.0/numChildren) * (xmax-xmin) );
  destinationY = rootY - levelIncrement;
  //cout << "\tI'm doing root: " << root << "\n";
  for (i=0; i<numChildren; i++) {
    if (shortEdges[ verts[root].e[i] ] != 1) {
      destinationX = destinationX + (int)( (1.0/numChildren) * (xmax-xmin) );
      OS << rootX << " " << rootY << " moveto\n";
      OS << destinationX << " " << destinationY << " lineto\n";
      OS << "stroke\n";
      //draw the edge label
      OS << rootX + (int)((destinationX - rootX)/2.0) - 21 << " " << rootY + (int)( (destinationY - rootY)/2.0) << " moveto\n";
      OS << "(" << edges[verts[root].e[i]].label << ")" << "show\n";     
      print_tree_ps(OS, verts, edges, shortEdges, edges[verts[root].e[i]].end, depth, 
                                                           destinationX, 
                                                           destinationY, 
                                                           (int)( (1.0/(numChildren)) * (xmax-xmin) ), 
                                                           levelIncrement);
    } else { //it's a short edge
      destinationX = destinationX + (int)( (1.0/numChildren) * (xmax-xmin) );
      destinationY += (int)(levelIncrement/2.0);
      OS << rootX << " " << rootY << " moveto\n";
      OS << destinationX << " " << destinationY << " lineto\n";
      OS << "stroke\n";
      //draw the edge label
      OS << rootX + (int)((destinationX - rootX)/2.0) - 17 << " " << rootY + (int)( (destinationY - rootY)/2.0) << " moveto\n";
      tempLabel = edges[verts[root].e[i]].label;
      if (verts[root].end[i] == -1) {
        swapCase(tempLabel);
      }
      OS << "(" << tempLabel << ")" << "show\n";  
      
      //draw the edge number so it can be matched
      OS << destinationX << " " << destinationY-10 << " moveto\n";
      OS << "(" << verts[root].e[i] << ")show\n";
      
        
      destinationY -= (int)(levelIncrement/2.0);      
    }
  }
}
  
//print the root and all edges (including short edges) then call overloaded print_tree_ps
//with the new root and current x,y coords and scaling factor
void print_tree_ps(ostream& OS, vector<vert>& verts, vector<edge>& edges, vector<int>& shortEdges, int root, int depth) {
  int xmin = 0;
  int xmax = 1024;
  int rootX = 512;
  int rootY = 1000;
  int levelIncrement = 1024/(depth+1);  
  OS << "\%!PS-Adobe-2.0 EPSF-2.0\n";
  OS << "\%\%BoundingBox: 0 0 1024 1024\n\n";      
  OS << "/Arial findfont\n";
  OS << "20 scalefont setfont\n";
  OS << "1 setlinejoin\n";
  OS << "2 setlinewidth\n";
  print_tree_ps(OS, verts, edges, shortEdges, root, depth, rootX, rootY, xmax-xmin, levelIncrement);
  OS << "\%eof\n";
}
    
  




void find_connected_components_and_min_depth_roots(vector<vert>& verts, 
                                                   vector<edge>& edges,
                                                   vector<int>& comps,
                                                   vector<int>& roots,
                                                   vector<int>& depths) {
  //First, find the connected components
  //this is redundant, but I'm just being careful
  comps.assign(verts.size(), -1);
  int root = 0;
  int i,j;
  int numVerts = verts.size();
  int numVertsFound = 0;
  queue<int> Q;
  int currentVertex;
  int destinationVertex;
  int currentComponentNum = 0;
  int numComps = 0;
  
  do {
    
    //find an unfound vertex to be the root
    for (i=0; i<(int)verts.size(); i++) {
      if (comps[i] == -1) {
        root = i;
        break;
      }
    }
    
    Q.push(root);
    
    //find the rest of the component 
    while (Q.size() > 0) {
      currentVertex = Q.front();
      Q.pop();
      comps[currentVertex] = currentComponentNum;
      numVertsFound++;
      
      for (j=0; j<(int)verts[currentVertex].e.size(); j++) {
        destinationVertex = ( verts[currentVertex].end[j] == 1 ?
                                edges[ verts[currentVertex].e[j] ].end :
                                edges[ verts[currentVertex].e[j] ].start );
        if (comps[destinationVertex] != -1) {
          continue;
        }
        Q.push(destinationVertex);
      }
    }
    
    //the queue is empty, so we're done with the component
    currentComponentNum++;
    
  } while (numVertsFound < (int)verts.size());
  
  //now each vertex has a component label, and the number of component is
  //currentComponentNum 
  
  numComps = currentComponentNum;
  
  //cout << "Found the number of components to be " << numComps << "\n";
  
  //now we go through each component and find the depth of the tree with that
  //vertex as root   
  vector<int> depthOfTreeWithRoot(verts.size());
  vector<int> tempDepths(verts.size());
  int currentDepth;
  for (i=0; i<numVerts; i++) {
    //find the depth for this vertex
    for (j=0; j<numVerts; j++) {
      tempDepths[j] = -1;
    }
    currentDepth = 0;
    Q.push(i);
    tempDepths[i] = 1;
    while (Q.size() > 0) {
      currentVertex = Q.front();
      Q.pop();
      for (j=0; j<(int)verts[currentVertex].e.size(); j++) {
        destinationVertex = ( verts[currentVertex].end[j] == 1 ?
                              edges[ verts[currentVertex].e[j] ].end :
                              edges[ verts[currentVertex].e[j] ].start );
        if (tempDepths[destinationVertex] != -1) {
          continue;
        }
        tempDepths[destinationVertex] = tempDepths[currentVertex] + 1;
        Q.push(destinationVertex);
      }
    }
    currentDepth = 1;
    for (j=0; j<numVerts; j++) {
      if (tempDepths[j] > currentDepth) {
        currentDepth = tempDepths[j];
      }
    }
    depthOfTreeWithRoot[i] = currentDepth;
  }
  
  //now for each connected component, find the minimum depth
  int whichComp = 0;
  roots.assign(numComps, -1);
  depths.assign(numComps, -1);
  for (whichComp = 0; whichComp < numComps; whichComp++) {
    for (i=0; i<numVerts; i++) {
      if (comps[i] == whichComp && 
           ( roots[whichComp] == -1 || 
             depthOfTreeWithRoot[i] < depthOfTreeWithRoot[roots[whichComp]])){
        roots[whichComp] = i;
      }
    }
    depths[whichComp] = depthOfTreeWithRoot[roots[whichComp]];
  }
 
}

//this rotates the list so element k is on the end, and then it deletes it
void rotate_list_and_truncate(vector<int>& L, int k) {
  vector<int> newL(L.size()-1);
  int i;
  int ls = L.size();
  for (i=0; i<ls-1; i++) {
    newL[i] = L[(k+i+1)%ls];
  }
  L = newL;
}



/* this creates a spanning tree out of a graph
It makes all the edges point down
roots should be a vector of vertices which are the roots of the connected 
components
*/
vector<int> make_trees_from_graphs(vector<vert>& verts, vector<edge>& edges, vector<int>& roots) {
  //do a BFS with the order done carefully
  vector<int> depths(verts.size());
  vector<int> shortEdges(edges.size());
  queue<int> Q;
  int i,j,k,r;
  int numVerts = verts.size();
  int currentVertex;
  int destinationVertex;
  int temp;
  int root;

  for (i=0; i<numVerts; i++) {
    depths[i] = -1;
  }
  for (i=0; i<(int)edges.size(); i++) {
    shortEdges[i] = 0;
  }
  
  for (r=0; r<(int)roots.size(); r++) {
    root = roots[r];
    Q.push(root);
    depths[root] = 1;
    while ((int)Q.size() > 0) {
      currentVertex = Q.front();
      Q.pop();
      for (j=0; j<(int)verts[currentVertex].e.size(); j++) {
        destinationVertex = ( verts[currentVertex].end[j] == 1 ?
                              edges[ verts[currentVertex].e[j] ].end :
                              edges[ verts[currentVertex].e[j] ].start );
        if (depths[destinationVertex] != -1) {
          if (shortEdges[ verts[currentVertex].e[j] ] == 1) { //ok it's already a short edge
          } else {
            shortEdges[ verts[currentVertex].e[j] ] = 1;
          }
          continue;
        }
        depths[destinationVertex] = depths[currentVertex] + 1;
        Q.push(destinationVertex);
        //change the order of the edges in destinationVertex, and change
        //this edge to point down
        for (k=0; k<(int)verts[destinationVertex].e.size(); k++) {
          if (verts[destinationVertex].e[k] == verts[currentVertex].e[j] && 
              verts[destinationVertex].end[k] == (-1) * verts[currentVertex].end[j]) {
            break;
          }
        }
        rotate_list_and_truncate(verts[destinationVertex].e, k);
        rotate_list_and_truncate(verts[destinationVertex].end, k);
        //change the edge
        if (verts[currentVertex].end[j] == -1) {
          verts[currentVertex].end[j] = 1;
          temp = edges[verts[currentVertex].e[j]].start;
          edges[verts[currentVertex].e[j]].start = edges[verts[currentVertex].e[j]].end;
          edges[verts[currentVertex].e[j]].end = temp;
          swapCase(edges[verts[currentVertex].e[j]].label);
        }
      }
    }
  }
  
  return shortEdges;
  
}







void print_graph_sage(vector<vert>& verts, vector<edge>& edges) {
  int i,j,k;
  int numVerts = verts.size();
  
 //ok now we should be done -- output a dictionary of dictionaries
  //for each vertex, output the edges starting there and with labels
  cout << "Sage graph: \n";
  cout << "{ ";
  bool first = true;
  for (i=0; i<numVerts; i++) {
    first = true;
    cout << i << ":{ ";
    for (j=0; j<numVerts; j++) {
      for (k=0; k<(int)verts[i].e.size(); k++) {
        if (edges[ verts[i].e[k] ].end == j && verts[i].end[k] == 1) {
          break;
        }
      }
      if (k == (int)verts[i].e.size()) {
        continue;
      }
      if (!first) {
        cout << ", ";
      } else {
        first = false;
      }
      cout << j << ":[";
      bool firstL = true;
      for (k=0; k<(int)verts[i].e.size(); k++) {
        if (edges[ verts[i].e[k] ].end == j && verts[i].end[k] == 1) {
          if (!firstL) {
            cout << ", ";
          } else {
            firstL = false;
          }
          cout << "'" << verts[i].e[k] << " - " << edges[ verts[i].e[k] ].label << "'";
        }
      }
      cout << "]";
    }
    cout << "}, ";
  }
  cout << "}\n";  
  
  cout << "Vertex edge orders:\n";
  cout << "{ ";
  for (i=0; i<numVerts; i++) {
    cout << i << ":(";
    cout << verts[i].e[0];
    if (verts[i].end[0] == -1) {
      cout << "*";
    }
    for (j=1; j<(int)verts[i].e.size(); j++) {
      cout << "," << verts[i].e[j];
      if (verts[i].end[j] == -1) {
        cout << "*";
      }
    }
    cout << "), ";
  }
  cout << "}\n"; 
}



void write_fatgraph_file(string filename, vector<vert>& verts, 
                                          vector<edge>& edges) {
  int i,j;
  fstream fgFile;
  fgFile.open(filename.c_str(), fstream::out);
  fgFile << "vertices " << verts.size() << "\n";
  for (i=0; i<verts.size(); i++) {
    fgFile << "VERT" << i << " " << verts[i].e.size() << "\n";
    for (j=0; j<verts[i].e.size(); j++) {
      fgFile << "EDGE" << verts[i].e[j] << " ";
    }
    fgFile << "\n";
    for (j=0; j<verts[i].e.size(); j++) {
      fgFile << (verts[i].end[j] == 1 ? 1 : 0) << " ";
    }
    fgFile << "\n";
  }
  fgFile << "edges " << edges.size() << "\n";
  for (i=0; i<edges.size(); i++) {
    fgFile << "EDGE" << i << " " << edges[i].label 
                            << " " << inverse(edges[i].label)
                            << " " << "VERT" << edges[i].start
                            << " " << "VERT" << edges[i].end 
                            << "\n";
  }
  fgFile.close();
}





void print_fund_group_gens(ostream& OS, vector<int>& shortEdges, vector<vert>& verts, 
                                                vector<edge>& edges, vector<int>& comps, vector<int>& roots) {
  //do a depth-first search of the tree, keeping the current path
  //when we reach a short edge, record the path we used to get here
  //When we have done this, match up the short edges to get the gens
  int i,r;
  vector<string> pathToVertex(verts.size(), "");
  vector<int> stack(0);
  int currentVertex,nextVertex;
  string word;
  
  for (r=0; r<(int)roots.size(); r++) {
  
    stack.push_back(roots[r]);
    
    //note that it is (better be) a tree, so we don't
    //need to worry about checking whether we've already visited
    while ((int)stack.size() > 0) {  
      currentVertex = stack.back();
      stack.pop_back();
      for (i=0; i<(int)verts[currentVertex].e.size(); i++) {
        if (shortEdges[ verts[currentVertex].e[i] ] == 1) {
          continue;
        }
        nextVertex = edges[verts[currentVertex].e[i]].end;
        if (nextVertex == currentVertex) { //this shouldn't happen
          cout << "Tree problems\n";
          return;
        }
        stack.push_back(nextVertex);
        pathToVertex[nextVertex] = pathToVertex[currentVertex] + edges[verts[currentVertex].e[i]].label;
      }
    }
  }
  
  //OS << "Paths to vertices:\n";
  //for (i=0; i<verts.size(); i++) {
  //  OS << pathToVertex[i] << "\n";
  //}
  //OS << "\n";
  
  
  //ok, now we just look at the short edges and give a generator for each one
  OS << "Generators for fundamental group:\n";
  for (r = 0; r<(int)roots.size(); r++) {
    OS << "For connected component " << r << "\n";
    for (i=0; i<(int)edges.size(); i++) {
      if (shortEdges[i] == 1 && comps[edges[i].end] == r) {
        currentVertex = edges[i].start;
        nextVertex = edges[i].end;
        word = pathToVertex[currentVertex] + edges[i].label + inverse(pathToVertex[nextVertex]);
        red(word);
        OS << word << "\n";
      }
    }
  }
  
  
}
      
  
  
  
  
void print_ends(ostream& endsOut, vector<int>& shortEdgeOrder, vector<string>& shortEdgeGenList) {
  int i,j;
  string tempEnd = "";
  int numEnds = shortEdgeOrder.size();
  int currentIndex;
  string gen;
  int inverseIndex;
  
  for (i=0; i<numEnds; i++) {
    tempEnd = "";
    currentIndex = i;
    while (tempEnd.size() < 30) {
      //cout << "currentIndex = " << currentIndex << "\n";
      if (shortEdgeOrder[currentIndex] < 0) {
        for (inverseIndex = 0; inverseIndex<numEnds; inverseIndex++) {
          if (shortEdgeOrder[inverseIndex] == -shortEdgeOrder[currentIndex]) {
            break;
          }
        }
        gen = inverse(shortEdgeGenList[inverseIndex]);
      } else {
        gen = shortEdgeGenList[currentIndex];
      }
      //cout << "gen = "<< gen << "\n";
      tempEnd = multiply_words(tempEnd, gen);
      for (j=0; j<numEnds; j++) {
        if (shortEdgeOrder[j] == -shortEdgeOrder[currentIndex]) {
          currentIndex = (j + (numEnds/2)) % numEnds;
          break;
        }
      }
      //cout << tempEnd << "\n";
    }
    cout << tempEnd << "\n";
  }
}
  
void print_polygon(ostream& OS, vector<int> shortEdgeOrder, 
                                         vector<string> shortEdgeGenList){
  int i;
  //print a polygon
  //for each edge:
  //if it is negative, it gets the number*
  //if positive, the number
  //on the key, the number is associated with the generator
  OS << "\%!PS-Adobe-2.0 EPSF-2.0\n";
  OS << "\%\%BoundingBox: 0 0 1024 1024\n\n";      
  OS << "/Arial findfont\n";
  OS << "20 scalefont setfont\n";
  OS << "1 setlinejoin\n";
  OS << "2 setlinewidth\n";
  
  //draw the polygon
  int numFaces = shortEdgeOrder.size();
  int centerX = 512;
  int centerY = 512;
  double rad = 256;
  double pi = 3.141592654;
  double currentX, currentY;
  
  OS << centerX+rad << " " << centerY << " moveto\n";  
  
  for (i=1; i<numFaces; i++) {
    currentX = centerX + rad*cos( (double)(2*pi*i)/(double)numFaces);
    currentY = centerY + rad*sin( (double)(2*pi*i)/(double)numFaces);
    OS << currentX << " " << currentY << " lineto\n";
  }
  OS << " closepath\nstroke\n";
  
  //now the labels
  for (i=0; i<numFaces; i++) {
    currentX = centerX + 1.1*rad*cos( (double)(2*pi*(i + 0.5))/(double)numFaces);
    currentY = centerY + 1.1*rad*sin( (double)(2*pi*(i + 0.5))/(double)numFaces);
    OS << currentX << " " << currentY << " moveto\n";
    if (shortEdgeOrder[i] < 0) {
      OS << "(" << -shortEdgeOrder[i]-1 << "*) dup stringwidth pop 2 div neg 0 rmoveto show\n";
    } else {
    OS << "(" << shortEdgeOrder[i]-1 << ") dup stringwidth pop 2 div neg 0 rmoveto show\n";
    }
  } 
  
  //now the key
  currentX = 100;
  currentY = 924;
  for (i=0; i<numFaces; i++) {
    if (shortEdgeOrder[i] < 0) {
      continue;
    }
    OS << currentX << " " << currentY << " moveto\n";
    OS << "(" << shortEdgeOrder[i]-1 << ": " << shortEdgeGenList[i] << ") show\n";
    currentY -= 30;
  }
   
  OS << "\%eof";
}
  
  
  
void print_polygon_gluing_to_file_and_ends(string drawFile, ostream& endsOut, vector<int>& shortEdges, 
                          vector<vert>& verts, 
                          vector<edge>& edges, vector<int>& comps, vector<int>& roots) {
                          
  //do a depth-first search of the tree, keeping the current path
  //when we reach a short edge, record the path we used to get here
  //When we have done this, match up the short edges to get the gens
  int i,r;
  vector<string> pathToVertex(verts.size(), "");
  vector<int> stack(0);
  int currentVertex,nextVertex;
  string word;
  
  
  fstream imageFile;
  string fileName = "";
  stringstream tempNum;
  
  
  for (r=0; r<(int)roots.size(); r++) {
  
    //keep a list of the order we see short edges -- if the short edge is 
    //starting at our vertex, it's positive, if it's ending, it's negative
    //(add 1 so this always makes sense and is confusing)
    vector<int> shortEdgeOrder;
    shortEdgeOrder.resize(0);
  
    stack.push_back(roots[r]);
    
    //note that it is (better be) a tree, so we don't
    //need to worry about checking whether we've already visited
    while ((int)stack.size() > 0) {  
      currentVertex = stack.back();
      stack.pop_back();
      for (i=0; i<(int)verts[currentVertex].e.size(); i++) {
        if (shortEdges[ verts[currentVertex].e[i] ] == 1) {
          if (edges[ verts[currentVertex].e[i] ].start == currentVertex &&
              verts[currentVertex].end[i] == 1) {
            shortEdgeOrder.push_back(verts[currentVertex].e[i]+1);
          } else {
            shortEdgeOrder.push_back(-(verts[currentVertex].e[i]+1));
          }
          continue;
        }
        nextVertex = edges[verts[currentVertex].e[i]].end;
        if (nextVertex == currentVertex) { //this shouldn't happen
          cout << "Tree problems\n";
          return;
        }
        stack.push_back(nextVertex);
        pathToVertex[nextVertex] = pathToVertex[currentVertex] + edges[verts[currentVertex].e[i]].label;
      }
    }
    
    //now we get a polygon with shortEdgeOrder.size() sides
    //first, we get the list of generators
    vector<string> shortEdgeGenList;
    shortEdgeGenList.resize(0);
    for (i=0; i<(int)shortEdgeOrder.size(); i++) {
      if (shortEdgeOrder[i] < 0) {
        shortEdgeGenList.push_back(""); // just a placeholder
        continue;
      }
      currentVertex = edges[shortEdgeOrder[i]-1].start;
      nextVertex = edges[shortEdgeOrder[i]-1].end;
      word = pathToVertex[currentVertex] + edges[shortEdgeOrder[i]-1].label + inverse(pathToVertex[nextVertex]);
      red(word);
      shortEdgeGenList.push_back(word);
    }
      
 
    fileName = drawFile;
    tempNum.str("");
    tempNum << r;
    fileName += tempNum.str();
    fileName += ".eps";
    imageFile.open(fileName.c_str(), fstream::out);
    
    print_polygon(imageFile, shortEdgeOrder, shortEdgeGenList);
    
    imageFile.close();
    
    cout << "Ends for component " << r << ":\n";
    print_ends(endsOut, shortEdgeOrder, shortEdgeGenList);
    
  }
  
}
  
  
  
  
  
  
  
  
  
  




void make_graph_and_print(string drawFile, vector<string> chain, 
                          vector<arc>& arc_list, vector<polygon>& poly_list, 
                          vector<rational>& weights) {
  vector<int> IWeights(weights.size());
  int i;
  int j;
  int k,l;
  int commonClear = 1;
  for (i=0; i<(int)weights.size(); i++) {
    //cout << "approxmated " << weights[i] << " by " << RWeights[i].n() << "/" << RWeights[i].d()<< "\n";
    commonClear = lcm(commonClear, weights[i].d());
  }
  for (i=0; i<(int)weights.size(); i++) {
    IWeights[i] = (commonClear/weights[i].d()) * weights[i].n();
  }
  //create the list of vertices (nonzero polys with multiplicity)
  int numVerts = 0;
  for (i=0; i<(int)weights.size(); i++) {
    numVerts += IWeights[i];
  }
  vector<vert> verts(numVerts);
  numVerts = 0;
  for (i=0; i<(int)weights.size(); i++) {
    if (IWeights[i] >0) {
      for (j=0; j<IWeights[i]; j++) {
        verts[numVerts].polyRep = i;
        verts[numVerts].e.assign(poly_list[i].size, -1);
        verts[numVerts].end.assign(poly_list[i].size, -1);
        numVerts++;
      }
    }
  }
  //cout << "I have created the vertex list:\n";
  //for (i=0; i<numVerts; i++) {
  //  cout << i << ": " << verts[i].polyRep << "\n";
  //  cout << "\t";
  //  for (j=0; j<verts[i].e.size(); j++) {
  //    cout << verts[i].e[j] << " " << verts[i].end[j] << "   ";
  //  }
  //  cout << "\n";
  //}
  //now we have all the vertices -- go back and add in the edges
  vector<edge> edges(0);
  int newEdgeNum;
  int startWordNum;
  int startCharNum;
  int endWordNum;
  int endCharNum;
  edge newE;
  while (1) {
    //find a vert with a missing edge
    for (i=0; i<numVerts; i++) {
      for (j=0; j<(int)verts[i].e.size(); j++) {
        if (verts[i].e[j] == -1) {
          break;
        }
      }
      if (j < (int)verts[i].e.size()) {
        break;
      }
    }
    if (i == numVerts) { // we're done
      break;
    }
    //cout << "I'm trying to fill in vertex " << i << " position " << j << "\n";
    //now vert i, edge j needs to be filled in
    startWordNum = arc_list[ poly_list[verts[i].polyRep].arc[j] ].first_word;
    startCharNum = arc_list[ poly_list[verts[i].polyRep].arc[j] ].first;
    endWordNum = arc_list[ poly_list[verts[i].polyRep].arc[j] ].last_word;
    endCharNum = arc_list[ poly_list[verts[i].polyRep].arc[j] ].last;
    newEdgeNum = edges.size();
    newE.start = i;
    newE.label = chain[ startWordNum ].substr( startCharNum, 1  );
    //cout << "Ok I'm trying to find the edge " << endCharNum << "," << startCharNum << "," << endWordNum << "," << startWordNum << "\n";
    //find a place to fill it in
    for (k = 0; k<numVerts; k++) {
      for (l = 0; l<(int)verts[k].e.size(); l++) {
        if ( verts[k].e[l] == -1 &&  //we only fill in empty edges
             arc_list[ poly_list[verts[k].polyRep].arc[l] ].first_word == endWordNum &&
             arc_list[ poly_list[verts[k].polyRep].arc[l] ].first == endCharNum &&
             arc_list[ poly_list[verts[k].polyRep].arc[l] ].last_word == startWordNum &&
             arc_list[ poly_list[verts[k].polyRep].arc[l] ].last == startCharNum ) {
          break;
        }
      }
      if (l < (int)verts[k].e.size()) {
        break;
      }
    }
    if (k == numVerts) {
      cout << "I can't find a matching edge somehow\n";
      return;
    }
    
    //cout << "Found it in vertex " << k << " at position " << l << "\n";
    
    newE.end = k;
    //fill in which edge it is
    verts[i].e[j] = newEdgeNum;
    verts[i].end[j] = 1;
    verts[k].e[l] = newEdgeNum;
    verts[k].end[l] = -1;
    edges.push_back(newE);
  }
  
  
  //print it
  print_graph_sage(verts, edges);
 
 
  //write it to a fatgraph file
  write_fatgraph_file(drawFile+".fg", verts, edges);
  
 
 
  
  //eps it
  vector<int> roots;
  vector<int> comps;
  vector<int> depths;
  find_connected_components_and_min_depth_roots(verts, edges, comps, roots, depths);
  
  cout << "There are " << roots.size() << " components\n";
  
  vector<int> shortEdges = make_trees_from_graphs(verts, edges, roots);
  
  //cout << "Created spanning trees\n";
  
  //opening the files

  fstream imageFile;
  string fileName = "";
  stringstream tempNum;
  for (i=0; i<(int)roots.size(); i++) {
    fileName = drawFile;
    tempNum.str("");
    tempNum << i;
    fileName += tempNum.str();
    fileName += ".eps";
    imageFile.open(fileName.c_str(), fstream::out);
    print_tree_ps(imageFile, verts, edges, shortEdges, roots[i], depths[i]);
    imageFile.close();
  }


  print_fund_group_gens(cout, shortEdges, verts, edges, comps, roots); 

  print_polygon_gluing_to_file_and_ends(drawFile+"P", cout, shortEdges, verts, edges, comps, roots);  

   
}
        
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
