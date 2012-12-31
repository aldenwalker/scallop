#include <iostream>

int main(int argc, char* argv[]) {

  enum {CYCLIC, LOCAL, TRAIN} comp_func;
  char** arg_array = NULL;
  
  if (argc < 2 || std::string(argv[1]) == "-h") {
    std::cout << "Scallop\n";
    std::cout << "Usage: ./scallop [-cyclic, -local, -train] <option-specific arguments>\n"; 
    std::cout << "\t-cyclic (default if no option is given):\n";
    std::cout << "\t\tUsage: ./scallop -cyclic [factor orders] <chain>\n";
    std::cout << "\t\te.g.: ./scallop -cyclic a0b5 abbAbbb\n";
    std::cout << "\t\tcomputes scl(abbAbbb) in Z*Z/5Z\n";
    std::cout << "\t\tIf omitted, the factor order string is taken\n";
    std::cout << "\t\tto be a free group (i.e. a0b0c0...)\n";
    std::cout << "\n";
    std::cout << "\t-local:\n";
    std::cout << "\t\tUsage: ./scallop -local [-f -ff[n] -e] <chain>\n";
    std::cout << "\t\tFind a minimal surface satisfying the contraints.\n";
    std::cout << "\t\t-f: require that the surface be folded\n";
    std::cout << "\t\t-ff[n]: require that the surface be f-folded\n";
    std::cout << "\t\t\tn is the length of C in C-phi(C).  phi(C) is marked\n";
    std::cout << "\t\t\twith periods e.g. -ff1 abAB .abba.baab.ABBA.BAAB\n";
    std::cout << "\t\t-e: only check if the surface exists (no minimization -- faster)\n";
    std::cout << "\n";
    std::cout << "\t-train:\n";
    std::cout << "\t\tUsage: ./scallop -train <chain>\n";
    std::cout << "\t\tCompute with traintracks (not implemented)\n";
    exit(0);
  }
  
  if (argv[1][0] != '-') {
    comp_func = CYCLIC;
    arg_array = &argv[1];
  } else if (argv[1][1] == 'c') {
    comp_func = CYCLIC;
    arg_array = &argv[2];
  } else if (argv[1][1] == 'l') {
    comp_func = LOCAL;
    arg_array = &argv[2];
  } else {
    comp_func = TRAIN;
    arg_array = &argv[2];
  }
  
  return 0;
}
    