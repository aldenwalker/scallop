#include <iostream>

#include "scylla/scylla.h"
#include "gallop/gallop.h"
#include "trollop/trollop.h"
#include "scabble/scabble.h"

int main(int argc, char* argv[]) {

  enum {CYCLIC, LOCAL, TRAIN, BALL} comp_func;
  char** arg_array = NULL;
  int num_args;
  
  if (argc < 2 || std::string(argv[1]) == "-h") {
    std::cout << "Scallop\n";
    std::cout << "Usage: ./scallop [-cyclic, -local, -train, -ball] <option-specific arguments>\n";
    std::cout << "Enter ./scallop <option> -h (e.g. ./scallop -cyclic -h) for specific info\n";
    std::cout << "\t-cyclic (default if no option is given):\n";
    std::cout << "\t\tCompute scl in free products of cyclic groups\n";
    std::cout << "\n";
    std::cout << "\t-local:\n";
    std::cout << "\t\tFind a minimal surface satisfying local contraints.\n";
    std::cout << "\n";
    std::cout << "\t-train:\n";
    std::cout << "\t\tCompute with traintracks\n";
    std::cout << "\t-ball:\n";
    std::cout << "\t\tCompute a ball in the scl norm in a free product of cyclic groups\n";
    return 0;
  }
  
  if (argv[1][0] != '-') {
    comp_func = LOCAL;
    arg_array = &argv[1];
    num_args = argc-1;
  } else if (argv[1][1] == 't') {
    comp_func = TRAIN;
    arg_array = &argv[2];
    num_args = argc-2;
  } else if (argv[1][1] == 'l') {
    comp_func = LOCAL;
    arg_array = &argv[2];
    num_args = argc-2;
  } else if (argv[1][1] == 'c') {
    comp_func = CYCLIC;
    arg_array = &argv[2];
    num_args = argc-2;
  } else if (argv[1][1] == 'b') {
    comp_func = BALL;
    arg_array = &argv[2];
    num_args = argc-2;
  } else {
    comp_func = LOCAL;
    arg_array = &argv[1];
    num_args = argc-1;
  }
  
  switch (comp_func) {
    case CYCLIC:
      SCYLLA::scylla(num_args, arg_array);
      break;
    case LOCAL:
      GALLOP::gallop(num_args, arg_array);
      break;
    case TRAIN:
      TROLLOP::trollop(num_args, arg_array);
      break;
    case BALL:
      SCABBLE::scabble(num_args, arg_array);
      break;
  }
  
  return 0;
}
    
