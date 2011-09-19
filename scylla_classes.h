/*****************************************************************************
* Bad form, but all the classes are going in here, because they are small
******************************************************************************/

#ifndef SCYLLA_CLASSES_H
#define SCYLLA_CLASSES_H

#include <iostream>
#include <string>
#include <vector>


/*****************************************************************************
 * a chunk in a chain
 * ***************************************************************************/
struct ChainChunk {
  int word;
  int group;
  int start_index;
  int len;
};

/*****************************************************************************
 * a letter in a chain
 * ***************************************************************************/
struct ChainLetter {
  int word;
  int index;
  char letter;
  int group;
};



/*****************************************************************************
* A free product of cyclic groups
* ****************************************************************************/
struct CyclicProduct {
  CyclicProduct(void);
  CyclicProduct(std::string input);
  ~CyclicProduct(void);

  int gen_order(char gen);                 //return the order of the given generator
  int index_order(int index);
  int gen_index(char gen);                 //return the number of the gen
  int num_groups(void);
  
  void cyc_red(std::string* S);                 //cyclically reduce a string
    
  std::vector<char> gens;
  std::vector<int> orders;
    
};
std::ostream &operator<<(std::ostream &os, CyclicProduct &G);

/*****************************************************************************
 * A chain   
 * ***************************************************************************/
struct Chain { 
  Chain(void);
  Chain(CyclicProduct* G, char** input, int num_strings);

  int next_letter(int n);
  int prev_letter(int n);
  int num_words(void);
  int num_letters();
  std::string operator[](int index);    //get a word
  void print_chunks(std::ostream &os);
  void print_letters(std::ostream &os);
  void print_group_letters(std::ostream &os);
  
  CyclicProduct* G;
  std::vector<std::string> words;
  std::vector<int> weights;
  std::vector<std::vector<ChainChunk> > chunks;
  std::vector<std::vector<int> > group_letters;
  std::vector<ChainLetter> chain_letters;
  std::vector<std::vector<int> > regular_letters;
  std::vector<std::vector<int> > inverse_letters;
    
};

std::ostream &operator<<(std::ostream &os, Chain &C);


/****************************************************************************
 * an edge joining two central polygons
 ****************************************************************************/
struct CentralEdge {
  int first;
  int last;
};

/****************************************************************************
 * A list of central edges
 * **************************************************************************/
struct CentralEdgeList {
  CentralEdgeList();
  CentralEdgeList(Chain &C);
  
  int get_index(int a, int b);
  CentralEdge operator[](int index);
  void print(std::ostream &os);
  
  std::vector<CentralEdge> edges;
  std::vector<std::vector<int> > edges_beginning_with;
};


/*****************************************************************************
 * an edge joining a central polygon to a group polygon
 * these are ALWAYS LISTED FROM THE POLYGON's PERSPECTIVE
 * ***************************************************************************/
struct InterfaceEdge {
  int first;
  int last;
};

/****************************************************************************
 * a list of interface edges
 * **************************************************************************/
struct InterfaceEdgeList {
  InterfaceEdgeList();
  InterfaceEdgeList(Chain &C);
  int get_index_from_poly_side(int a, int b);
  int get_index_from_group_side(int a, int b);
  InterfaceEdge operator[](int index);
  void print(std::ostream &os);
  
  std::vector<InterfaceEdge> edges;
  std::vector<std::vector<int> > edges_beginning_with;
};

/*****************************************************************************
 * an edge between group polygons
 * ***************************************************************************/
struct GroupEdge {
  int first;
  int last;
};

/****************************************************************************
 * a list of group edges
 * **************************************************************************/
struct GroupEdgeList {
  GroupEdgeList();
  GroupEdgeList(Chain &C, int group_index);
  int get_index(int a, int b);
  void print(std::ostream &os);
  GroupEdge operator[](int index);
  
  int group;
  std::vector<GroupEdge> edges;
  std::vector<std::vector<int> > edges_beginning_with;
  std::vector<int> regular_edges;
  std::vector<int> inverse_edges;
};

/****************************************************************************
 * an edge pair
 * **************************************************************************/
enum EDGE_TYPE {SCYLLA_EDGE_CENTRAL, SCYLLA_EDGE_INTERFACE, SCYLLA_EDGE_GROUP};
struct EdgePair {
  EDGE_TYPE edge_type;
  int group;
  int first;
  int second;
};


/*****************************************************************************
 * a central polygon  (this is a list of interface and polygon edges)
 * ***************************************************************************/
struct CentralPolygon {
  std::vector<int> edges;
  std::vector<bool> interface;
  int chi_times_2(Chain &C, CentralEdgeList &CEL, InterfaceEdgeList &IEL);
  void compute_ia_etc_for_edges(int i,
                                Chain &C,
                                InterfaceEdgeList &IEL,
                                CentralEdgeList &CEL,
                                std::vector<EdgePair> &edge_pairs,
                                std::vector<int> &central_edge_pairs,
                                std::vector<int> &temp_ia,
                                std::vector<int> &temp_ja,
                                std::vector<ing> &temp_ar);
};

std::ostream &operator<<(std::ostream &os, CentralPolygon &CP);

/****************************************************************************
 * a group outside edge (these pair with the interface edges
 * **************************************************************************/
struct GroupTooth {
  int position;
  int first;
  int last;
  int chi_times_2(Chain &C);
  void compute_ia_etc_for_edges(int offset, 
                                Chain &C,
                                InterfaceEdgeList &IEL, 
                                std::vector<int> &rows_for_letters, 
                                std::vector<int> &ia, 
                                std::vector<int> &ja, 
                                std::vector<int> &ar);
  void compute_ia_etc_for_words(int offset, 
                                int row_offset,
                                Chain &C,
                                std::vector<int> &ia, 
                                std::vector<int> &ja, 
                                std::vector<int> &ar);
                                
};

std::ostream &operator<<(std::ostream &os, GroupTooth &GT);

/****************************************************************************
 * These are really just special group polygons designated to face the teeth
 * They join the group edge (first, last) (it's not reversed)
 * **************************************************************************/
struct GroupMouth {
  int first;
  int last;
  int chi_times_2(Chain &C);
  void compute_ia_etc_for_edges(int offset, 
                                int group_index,
                                Chain &C,
                                GroupEdgeList &GEL,
                                std::vector<int> &rows_for_letters,
                                std::vector<EdgePair> &edge_pairs,
                                std::vector<int> &group_edge_pairs,
                                std::vector<int> &ia,
                                std::vector<int> &ja,
                                std::vector<int> &ar);
};

std::ostream &operator<<(std::ostream &os, GroupMouth &GM);

/****************************************************************************
 * a group polygon -- these are really just squares with all incident 
 * edges as group edges -- when computing chi, we must ignore the joined edges
 * **************************************************************************/
struct GroupPolygon {
  int group;
  std::vector<int> edges;
  int chi_times_2(GroupEdgeList &GEL);
  void get_ia_etc_for_edges(int offset, 
                            Chain &C, 
                            InterfaceEdgeList &IEL, 
                            GroupEdgeList &GEL, 
                            std::vector<EdgePair> &edge_pairs,
                            std::vector<int> &group_edge_pairs, 
                            std::vector<int> &temp_ia, 
                            std::vector<int> &temp_ja, 
                            std::vector<int> &temp_ar);
};

std::ostream &operator<<(std::ostream &os, GroupPolygon &GP);

/****************************************************************************
 * a group rectangle
 * **************************************************************************/
struct GroupRectangle {
  int first;
  int last;
  void compute_ia_etc_for_edges(int offset, 
                                std::vector<int> &temp_ia, 
                                std::vector<int> &temp_ja, 
                                std::vector<int> &temp_ar);
  void compute_ia_etc_for_words(int offset, 
                                int row_offset,
                                Chain &C, 
                                InterfaceEdgeList &IEL,
                                std::vector<int> &ia,
                                std::vector<int> &ja,
                                std::vector<int> &ar);
                                
};
 
std::ostream &operator<<(std::ostream &os, GroupRectangle &GR);


/*****************************************************************************
 * a multiset
 * ***************************************************************************/
class Multiset {
  public:
    Multiset();
    Multiset(int len, int Min, int Max_plus_one);
    int next(void);
    std::vector<int>* get_list(void);
    int operator[](int index);
    void set_index(int index, int val);
  private:
    std::vector<int> L;
    int min;
    int len;
    int max_plus_one;
};



#endif 
