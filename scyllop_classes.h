/*****************************************************************************
* Bad form, but all the classes are going in here, because they are small
******************************************************************************/

#ifndef SCYLLOP_CLASSES_H
#define SCYLLOP_CLASSES_H

#include <iostream>
#include <string>
#include <vector>


/*****************************************************************************
 * a letter in a chain
 * ***************************************************************************/
struct ChainLetter {
  int word;
  int index;
};

/*****************************************************************************
 * A multiarc
 * ***************************************************************************/
struct Multiarc {
  int group;
  std::vector<ChainLetter> letters;
};

/*****************************************************************************
 * a poylgon
 * ***************************************************************************/
struct Polygon {
  std::vector<int> arcs;
};



/*****************************************************************************
* A free product of cyclic groups
* ****************************************************************************/
class CyclicProduct {
  
  public:
    CyclicProduct(void);
    CyclicProduct(std::string input);
    ~CyclicProduct(void);
  
    std::vector<char> gen_list(void);        //return the generator list
    std::vector<int> order_list(void);       //return the list of orders (0=inf)
    int gen_order(char gen);                 //return the order of the given generator
    
    void cyc_red(std::string* S);                 //cyclically reduce a string
    
    friend std::ostream &operator<<(std::ostream &os, CyclicProduct &G);
  
  private:
    std::vector<char> gens;
    std::vector<int> orders;
    
};


/*****************************************************************************
 * A chain   
 * ***************************************************************************/
class Chain {
  
  public:
    Chain(void);
    Chain(CyclicProduct* G, char** input, int num_strings);
    ~Chain(void);
    
    std::vector<std::string> word_list(void);
    std::vector<int> weights_list(void);
    //std::vector<std::vector<ChainChunk> > chunk_list(void);
    int num_words(void);
    std::string operator[](int index);    //get a word
    //void print_chunks(std::ostream &os);
    
    friend std::ostream &operator<<(std::ostream &os, Chain &C);
    
  private:
    CyclicProduct* G;
    std::vector<std::string> words;
    std::vector<int> weights;
    std::vector<std::vector<ChainChunk> > chunks;
    
};




/*****************************************************************************
 * a multiset
 * ***************************************************************************/
class Multiset {
  public:
    Multiset();
    Multiset(int len, int Min, int Max_plus_one);
    int next(void);
    std::vector<int>* get_list(void);
  private:
    std::vector<int> L;
    int min;
    int max_plus_one;
};



#endif 
