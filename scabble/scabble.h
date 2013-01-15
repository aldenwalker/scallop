


namespace SCABBLE {

  /*****************************************************************************
  * a letter in a chain
  * ***************************************************************************/
  struct ChainLetter {
    int word;
    int index;
    int index_in_group_reg_inv_list;
    char letter;
    int group;
  };

  std::ostream &operator<<(std::ostream &os, ChainLetter &CL);

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
    
    void cyc_red(std::string &S);                 //cyclically reduce a string
      
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
    std::vector<std::vector<int> > group_letters;
    std::vector<ChainLetter> chain_letters;
    std::vector<std::vector<int> > regular_letters;
    std::vector<std::vector<int> > inverse_letters;
      
  };

  std::ostream &operator<<(std::ostream &os, Chain &C);
    
  /****************************************************************************
 * an edge pair which can join central polygons 
 * Note this pair is edges (i,j) and (j-1,i+1)
 ****************************************************************************/
  struct CentralEdgePair {
    int first;
    int last;
  };
  /****************************************************************************
  * A list of central edges
  * (we will record only the one with minimal first letter.)
  * Also, we don't allow the edge pair (i,i+1), (i, i+1)
  * **************************************************************************/
  struct CentralEdgePairList {
    CentralEdgePairList();
    CentralEdgePairList(Chain &C);
    
    int get_index(int a, int b);      //this returns ind+1 for the first edge and -(ind+1) for the second edge
    CentralEdgePair operator[](int index);
    void print(std::ostream &os);
    int size();
    
    int num_letters;
    Chain* my_chain;;
    std::vector<CentralEdgePair> edge_pairs;
    std::vector<std::vector<int> > edge_pairs_beginning_with;
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
    int size();
    
    std::vector<InterfaceEdge> edges;
    std::vector<std::vector<int> > edges_beginning_with;
  };



  /*****************************************************************************
  * a central polygon  (this is a list of interface and polygon edges)
  * ***************************************************************************/
  struct CentralPolygon {
    std::vector<std::pair<int, int> > edges; //these record the first and second letters in the edge pair
    std::vector<bool> interface;   //this records whether the edge is an interface edge
    int chi_times_2();
    void compute_ia_etc_for_edges(int col,
                                  Chain &C,
                                  InterfaceEdgeList &IEL,
                                  CentralEdgePairList &CEL,
                                  SparseLP& LP);
  };

  std::ostream &operator<<(std::ostream &os, CentralPolygon &CP);

  /****************************************************************************
  * a group outside edge (these pair with the interface edges)
  * **************************************************************************/
  struct GroupTooth {
    int position;     //the position of this tooth around the circle
    int first;        //the first letter (as position in the chain)
    int last;         //the last letter 
    bool inverse;     //is it made up of inverse letters?
    int group_index;  //index of the group
    int base_letter;  //the base letter (as position in the chain)
    double chi_times_2(Chain &C);
    void compute_ia_etc_for_edges(int offset, 
                                  Chain &C,
                                  InterfaceEdgeList &IEL, 
                                  std::vector<std::vector<std::vector<int> > > &group_teeth_rows, 
                                  SparseLP& LP);
    void compute_ia_etc_for_words(int offset, 
                                  Chain &C,
                                  int row_offset,
                                  SparseLP& LP);
                                  
  };

  std::ostream &operator<<(std::ostream &os, GroupTooth &GT);


  /****************************************************************************
  * a group rectangle
  * **************************************************************************/
  struct GroupRectangle {
    int group_index;
    int first;            //the first letter of the rectangle (position in the chain)
    int last;             //the second letter of the rectangle
    void compute_ia_etc_for_edges(int offset, 
                                  InterfaceEdgeList &IEL, 
                                  SparseLP& LP);
    void compute_ia_etc_for_words(int offset, 
                                  Chain &C, 
                                  int row_offset,
                                  InterfaceEdgeList &IEL,
                                  SparseLP& LP);
                                  
  };
  
  std::ostream &operator<<(std::ostream &os, GroupRectangle &GR);
    
  int scabble(int argc, char** argv);
  
}