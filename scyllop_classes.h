/*****************************************************************************
* Bad form, but all the classes are going in here, because they are small
******************************************************************************/

#ifndef SCYLLOP_CLASSES_H
#define SCYLLOP_CLASSES_H

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
 * A multiarc
 * ***************************************************************************/
struct Multiarc {
  int group;
  std::vector<int> letters;
};

/****************************************************************************
 * an edge
 ****************************************************************************/
struct Edge {
  bool blank;
  int first;
  int last;
};


/*****************************************************************************
 * a poylgon
 * ***************************************************************************/
struct Polygon {
  std::vector<int> edges;
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
    int gen_index(char gen);                 //return the number of the gen
    int num_groups(void);
    
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
    std::vector<std::vector<ChainChunk> > chunk_list(void);
    std::vector<std::vector<int> > group_letter_list(void);
    std::vector<ChainLetter> chain_letter_list(void);
    int next_letter(int n);
    int prev_letter(int n);
    int num_words(void);
    std::string operator[](int index);    //get a word
    void print_chunks(std::ostream &os);
    void print_letters(std::ostream &os);
    void print_group_letters(std::ostream &os);
    
    friend std::ostream &operator<<(std::ostream &os, Chain &C);
    
  private:
    CyclicProduct* G;
    std::vector<std::string> words;
    std::vector<int> weights;
    std::vector<std::vector<ChainChunk> > chunks;
    std::vector<std::vector<int> > group_letters;
    std::vector<ChainLetter> chain_letters;
    
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
    int operator[](int index);
  private:
    std::vector<int> L;
    int min;
    int len;
    int max_plus_one;
};



#endif 
