#include <set>
#include <cctype>

#include "word.h"


using namespace std;

/*
struct pair{
	int loc;
	int siz;
};
*/

char swapCaseChar(char c) {
  return (char)((int)c > 96 ? (int)c-32 : (int)c+32);
}

char inverse_char(char c) {
  return (char)((int)c > 96 ? (int)c-32 : (int)c+32);
}

void swapCase(string& s) {
  int i;
  int sl = s.size();
  for (i=0; i<sl; i++) {
    s[i] = 	(char)(s[i] > 96 ? (int)s[i]-32 : (int)s[i]+32);
  }
}

string inverse(string& s) {
  string ans;
  int i;
  int ss = s.size();
  ans.resize(ss);
  for (i=0; i<ss; i++) {
    ans[ss-i-1] = s[i];
  }
  swapCase(ans);
  return ans;
}

void red(string& s) {			// reduce by cancelling adjacent inverse letters
	int i,a,b;
	i=0;
	while(i<=(int) s.length()){
		if(i>0){
			a = (int) s[i-1];
			b = (int) s[i];
			if((32+a-b)%64==0){
	//			cout << (*s) << " " << i << "\n";
				s.erase(i-1,2);
	//			cout << (*s) << "\n";
				i=i-2;
			};
		};
		i++;
	};
	return;
};


void cyc_red(string& s) {
  while ( (s.size() > 0) && (s[0] == swapCaseChar(s[s.size()-1])) ) {
    s.erase(0,1);
    s.erase(s.size()-1, 1);
  }
}




string multiply_words(string& s1, string& s2) {
  string ans = s1+s2;
  red(ans);
  return ans;
}


//returns 0 for a/A, etc
int letter_index(char let) {
  int val = (int)let;
  if (val < 97) val += 32;
  val -= 97;
  return val;
}


std::vector<char> chain_gens(int num_words, char** words) {
  std::set<char> found_letters;
  for (int i=0; i<num_words; ++i) {
    int j=0;
    while (words[i][j] != '\0') {
      found_letters.insert(tolower(words[i][j]));
      ++j;
    }
  }
  std::vector<char> gens(found_letters.begin(), found_letters.end());
  return gens;
}

int chain_rank(int num_words, char** words) {
  int max_found = 0;
  for (int i=0; i<num_words; ++i) {
    int j=0; 
    while (words[i][j] != '\0') {
      int val = (int)words[i][j];
      if ((val < 65) || (val > 90 && val < 97) || (val > 122)) {
        j++;
        continue;
      }
      if (val < 97) val += 32;
      val -= 97;
      val += 1;
      if (val > max_found) {
        max_found = val;
      }
      j++;
    }
  }
  return max_found;
}

int chain_rank(std::vector<std::string >& words) {
  int max_found = 0;
  for (int i=0; i<(int)words.size(); ++i) {
    for (int j=0; j<(int)words[i].size(); ++j) {
      int val = (int)words[i][j];
      if ((val < 65) || (val > 90 && val < 97) || (val > 122)) {
        j++;
        continue;
      }
      if (val < 97) val += 32;
      val -= 97;
      val += 1;
      if (val > max_found) {
        max_found = val;
      }
      j++;
    }
  }
  return max_found;
}


int chain_rank(std::string& word) {
  int max_found = 0;
  for (int j=0; j<(int)word.size(); ++j) {
    int val = (int)word[j];
    if ((val < 65) || (val > 90 && val < 97) || (val > 122)) {
      j++;
      continue;
    }
    if (val < 97) val += 32;
    val -= 97;
    val += 1;
    if (val > max_found) {
      max_found = val;
    }
  }
  return max_found;
}



string reverse(string s){		// returns reverse of a word
	string t,u;
	
	
	if(s.length()==0){
		return(s);
	} else {
		u=s;
		t=u.erase(0,u.length()-1) + reverse(s.erase(s.length()-1,1));
		return(t);
	};
}
