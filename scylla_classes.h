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
  void get_index_from_poly_side(int a, int b);
  void get_index_from_group_side(int a, int b);
  InterfaceEdge operator[](int index);
  
  std::vector<InterfaceEdge> edges;
  std::vector<std:vector<int> > edges_beginning_with;
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
  GroupEdgeList(Chain &C);
  void get_index(int a, int b);
  
  std::vector<GroupEdge> edges;
  std::vector<std::vector<int> > edges_beginning_with;
  std::vector<int> regular_edges;
  std::vector<int> inverse_edges;
};



/*****************************************************************************
 * a central polygon  (this is a list of interface and polygon edges)
 * ***************************************************************************/
struct CentralPolygon {
  std::vector<int> edges;
  std::vector<bool> interface;
  int chi_times_2(CentralEdgeList &CEL, InterfaceEdgeList &IEL);
};


/****************************************************************************
 * a group multiarc
 * **************************************************************************/
struct Multiarc {
  std::vector<int> letters;
};

/****************************************************************************
 * a group polygon
 * **************************************************************************/
struct GroupPolygon {
  int group;
  std::vector<Multiarc> sides;
  std::vector<int> edges;
  int chi_times_2(GroupEdgeList &GEL, InterfaceEdgeList &IEL);
};

/****************************************************************************
 * a group rectangle
 * **************************************************************************/
struct GroupRectangle {
  int group;
  std::vector<int> edges;
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
  private:
    std::vector<int> L;
    int min;
    int len;
    int max_plus_one;
};



#endif 
