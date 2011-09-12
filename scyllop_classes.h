/*****************************************************************************
* Bad form, but all the classes are going in here, because they are small
******************************************************************************/

#ifndef SCYLLOP_CLASSES_H
#define SCYLLOP_CLASSES_H

#include <iostream>
#include <string>
#include <vector>

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
    Chain(char** input, int num_strings);
    ~Chain(void);
    
    std::vector<std::string> word_list(void);
    std::vector<int> weights_list(void);
    int num_words(void);
    std::string operator[](int index);    //get a word
    
    friend std::ostream &operator<<(std::ostream &os, Chain &C);
    
  private:
    std::vector<std::string> words;
    std::vector<int> weights;
    
};


/*****************************************************************************
 * a letter in a chain
 * ***************************************************************************/
struct ChainChunk {
  int word;
  int start_index;
  int len;
};

/*****************************************************************************
 * A multiarc
 * ***************************************************************************/
struct Multiarc {
  int group;
  std::vector<ChainChunk> letters;
};

/*****************************************************************************
 * a poylgon
 * ***************************************************************************/
struct Polygon {
  std::vector<int> arcs;
};




#endif 
