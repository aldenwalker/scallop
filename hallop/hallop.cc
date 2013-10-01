#include <string>
#include <iostream>

#include "free_group_chain.h"
#include "hallop.h"

void HALLOP::hallop(int argc, char** argv) {
  if (argc < 1 || std::string(argv[0]) == "-h") {
    std::cout << "usage: ./scallop -hyp [-R<relator>] <chain>\n";
    std::cout << "\twhere <chain> allows integral weights on the words\n";
    std::cout << "\t-h: print this message\n";
    std::cout << "\t-R relator: add a relator\n";
    std::cout << "\tExample: ./scallop -hyp -RabABcdCD abAB\n";
  }
  
  std::vector<std::string> relators(0);
  
  int current_arg = 0;
  while (argv[current_arg][0] == '-') {
    if (argv[current_arg][1] == 'R') {
      relators.push_back( std::string(&argv[current_arg][2]) );
    }
    current_arg++;
  }
  
  FreeGroupChain C(&argv[current_arg], argc-current_arg);
  
  std::cout << "Got relators:\n";
  for (int i=0; i<(int)relators.size(); ++i) {
    std::cout << i << ": " << relators[i] << "\n";
  }
  std::cout << "Got chain: " << C << "\n";
  
  FreeGroupChain CR = C; 
  
  for (int i=0; i<(int)relators.size(); ++i) {
    CR.add_word(-1, relators[i]);
  }
  
  std::cout << "Chain with relators: " << CR << "\n";
  
}