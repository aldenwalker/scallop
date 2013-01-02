#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdio>

#include "graph.h"

using namespace GALLOP;

namespace GALLOP {
void extract_signed_index(int si, int* ind, int* sign) {
  if (si < 0) {
    *ind = -(si+1);
    *sign = -1;
  } else {
    *ind = si-1;
    *sign = 1;
  }
}
}


int Graph::get_edge_from_label(char label) {
  int i;
  for (i=0; i<(int)edges.size(); i++) {
    if (edges[i].label_forward[0] == label) {
      return i+1;
    } else if (edges[i].label_backward[0] == label) {
      return -(i+1);
    }
  }
  return 0;
}



void Graph::set_standard_rose(int rank, int verbose) {
  num_verts = 1;
  num_edges = rank;
  verts.resize(1);
  edges.resize(0);
  Vert V;
  V.name = "VERT0";
  V.incident_edges.resize(2*rank);
  V.is_outgoing.resize(2*rank);
  Edge E;
  char c[2];
  c[1] = '\0';
  for (int i=0; i<rank; ++i) {
    V.incident_edges[2*i] = i;
    V.incident_edges[2*i+1] = i;
    V.is_outgoing[2*i] = true;
    V.is_outgoing[2*i+1] = false;
    c[0] = (char)(48+i);
    E.name = "EDGE" + std::string(c);
    c[0] = (char)(97+i);
    E.label_forward = std::string(c);
    c[0] = (char)(65+i);
    E.label_backward = std::string(c);
    E.source = 0;
    E.dest = 0;
    edges.push_back(E);
  }
  verts[0] = V;
}



int Graph::write_file(std::string filename, int verbose) {
  int i,j;
  std::fstream ofile;
  ofile.open(filename.c_str(), std::fstream::out);

  if (verbose>1) std::cout << "Writing fatgraph to file " << filename << "\n";
  
  ofile << "#output produced from fatgraph program\n";
  ofile << "vertices " << num_verts << "\n";
  for (i=0; i<num_verts; i++) {
    if (verts[i].name == "") {
      ofile << "VERT" << i << " " << verts[i].incident_edges.size() << "\n";
    } else {
      ofile << verts[i].name << " " << verts[i].incident_edges.size() << "\n";
    }
    for (j=0; j<(int)verts[i].incident_edges.size(); j++) {
      if (edges[verts[i].incident_edges[j]].name == "") {
        ofile << "EDGE" << verts[i].incident_edges[j] << " ";
      } else {
        ofile << edges[verts[i].incident_edges[j]].name << " ";
      }
    }
    ofile << "\n";
    for (j=0; j<(int)verts[i].incident_edges.size(); j++) {
      ofile << (verts[i].is_outgoing[j] ? 1 : 0) << " ";
    }
    ofile << "\n";
  }
  ofile << "edges " << num_edges << "\n";
  for (i=0; i<num_edges; i++) {
    if (edges[i].name == "") {
      ofile << "EDGE" << i << " ";
    } else {
      ofile << edges[i].name << " ";
    }
    ofile << edges[i].label_forward << " " << edges[i].label_backward << " ";
    if (verts[edges[i].source].name == "") {
      ofile << "VERT" << edges[i].source << " ";
    } else {
      ofile << verts[edges[i].source].name << " ";
    }
    if (verts[edges[i].dest].name == "") {
      ofile << "VERT" << edges[i].dest << "\n";
    } else {
      ofile << verts[edges[i].dest].name << "\n";
    }
  }
  ofile.close();
  return 0;
}
  
  
  
int Graph::read_file(std::string filename, int verbose) {
  FILE* ifile = fopen(filename.c_str(), "r");

  int i,j,k;
  char tempc;
  char temps[1000];
  char temp_edge_name[50];
  char temp_label_forward[50];
  char temp_label_backward[50];
  char temp_vert_start[50];
  char temp_vert_end[50];
  int tempi;
  int tempe;
  
  
  while (fgetc(ifile) == '#') {  
    do {
      tempc = fgetc(ifile);
    } while (tempc != EOF && tempc != '\n');
  }
  fseek(ifile, -1, SEEK_CUR);
  fscanf(ifile, "vertices %d\n", &num_verts);
  verts.resize(num_verts);
  for (i=0; i<num_verts; i++) {
    fscanf(ifile, "%s %d\n", temps, &tempe);
    if (verbose > 1) std::cout << "new vertex " << temps << " with " << tempe << " edges\n"; 
    verts[i].name = std::string(temps);  
    if (tempe > 0 ) {
      verts[i].incident_edges.resize(tempe);
      verts[i].is_outgoing.resize(tempe);
    } else {
      verts[i].incident_edges.resize(0);
      verts[i].is_outgoing.resize(0);
    }
    
    if (verts[i].incident_edges.size() > 0) {
      //skip the edge line for now
      fgets(temps, 1000, ifile);
      
      //read in the edges direction(even though we don't know what they are!)
      for (j=0; j<(int)verts[i].incident_edges.size(); j++) {
        fscanf(ifile, "%d", &tempi);
        //std::cout << "Read direction: " << tempi << "\n";
        verts[i].is_outgoing[j] = (bool)tempi;
      }
    }
    
    //read the bezier and location, if we've got them
    fscanf(ifile, "\n%s", temps);
    if (verbose > 1) std::cout << "next word: " << temps << "\n";
    if (strcmp(temps, "bezier") == 0) {
      if (verbose > 1) std::cout << "It's got bezier stuff\n";
      for (j=0; j<(int)verts[i].incident_edges.size(); j++) {
        fscanf(ifile, "%*f %*f");
      }
      fscanf(ifile, "%s", temps);
    } else {
      //there's no bezier, so deal with that
      if (verbose > 1) std::cout << "There's no bezier\n";
    }
    
    if (strcmp(temps, "loc") == 0) {
      fscanf(ifile, "%*f %*f");
    } else {
      //there's no loc, so deal with that
      fseek(ifile, -strlen(temps)-1, SEEK_CUR);
    } 
  }
  
  if (verbose > 1) std::cout << "I've read in the vertices\n";  
  
  fscanf(ifile, " edges %d\n", &num_edges);
  if (verbose > 1) std::cout << "There should be " << num_edges << " edges\n";
  edges.resize(num_edges);
  for (i=0; i<num_edges; i++) {
    fscanf(ifile, " %s %s %s %s %s\n", temp_edge_name, 
                                       temp_label_forward,
                                       temp_label_backward,
                                       temp_vert_start,
                                       temp_vert_end);
    if (verbose > 1) {
      std::cout << "Read edge: " << temp_edge_name << " " 
                               << temp_label_forward << " "
                               << temp_label_backward << " "
                               << temp_vert_start << " "
                               << temp_vert_end << "\n";
    }
    edges[i].name = std::string(temp_edge_name);
    edges[i].label_forward = std::string(temp_label_forward);
    edges[i].label_backward = std::string(temp_label_backward);
    for (j=0; j<num_verts; j++) {
      if (std::string(temp_vert_start) == verts[j].name) {
        edges[i].source = j;
        break;
      }
    }
    for (j=0; j<num_verts; j++) {
      if (std::string(temp_vert_end) == verts[j].name) {
        edges[i].dest = j;
        break;
      }
    }
  }
  
  //reset to the beginning
  fseek(ifile, 0, SEEK_SET);
  
  if (fgetc(ifile) == '#') {
    do {
      tempc = fgetc(ifile);
    } while (tempc != EOF && tempc != '\n');
  }
  
  fscanf(ifile, "vertices %*d\n");
  for (i=0; i<num_verts; i++) {
    fscanf(ifile, " %s ", temps);
    while (std::string(temps) != verts[i].name) {
      fscanf(ifile, " %s ", temps);
    }
    fscanf(ifile, "%*d"); //get rid of the number of edges
    for (j=0; j<(int)verts[i].incident_edges.size(); j++) {
      fscanf(ifile, " %s ", temp_edge_name);
      for (k=0; k<num_edges; k++) {
        if (std::string(temp_edge_name) == edges[k].name) {
          break;
        }
      }
      verts[i].incident_edges[j] = k;
    } 
  }
  
  return 0;
} 


