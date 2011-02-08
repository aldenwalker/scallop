#ifndef __word__
#define __word__

#include <string>


void red(std::string& s);			// reduce by cancelling adjacent inverse letters
void swapCase(std::string& s);
void cyc_red(std::string& s);
std::string inverse(std::string& s);
std::string reverse(std::string s);
std::string multiply_words(std::string& s1, std::string& s2);

#endif
