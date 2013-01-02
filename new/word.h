#ifndef __word__
#define __word__

#include <string>


void red(std::string& s);			// reduce by cancelling adjacent inverse letters
void swapCase(std::string& s);
char swapCaseChar(char c);
char inverse_char(char c);
void cyc_red(std::string& s);
std::string inverse(std::string& s);
std::string reverse(std::string s);
std::string multiply_words(std::string& s1, std::string& s2);
int chain_rank(int num_words, char** words); 

#endif
