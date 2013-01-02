#ifndef graph_H
#define graph_H

#include <string>
#include <vector>

namespace GALLOP {

struct Edge {
  std::string name;
  std::string label_forward;
  std::string label_backward;
  int source;
  int dest;
};

struct Vert {
  std::string name;
  std::vector<int> incident_edges;
  std::vector<bool> is_outgoing;
};
  
struct Graph {
  int num_verts;
  int num_edges;
  std::vector<Edge> edges; //the edges 
  std::vector<Vert> verts; //the verts
  
  int get_edge_from_label(char label);
  void set_standard_rose(int rank, int verbose);
  int read_file(std::string file_name, int verbose);
  int write_file(std::string file_name, int verbose);
};
  
void extract_signed_index(int si, int* ind, int* sign);
  
}
#endif