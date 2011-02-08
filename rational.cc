#include <vector>
#include <iostream>
#include <math.h>

#include "rational.h"


using namespace std;


int gcd(int a, int b) {
  int t1 = (b > a ? b : a); //max
  int t2 = (b > a ? a : b); //min
  int t3;
  while (t2 != 0) {
    t3 = t1 % t2;
    t1 = t2;
    t2 = t3;
  }
  return t1;
}


int lcm(int a, int b) {
  int g = gcd(a,b);
  return (a*b)/g;
}




/******************************************************************************/
/* member functions for the rational class                                    */
/******************************************************************************/
rational::rational() {
  mpq_init(R);
  inited = 1;
  mpq_set_si(R,0,1);
}

rational::rational(int a, int b) {
  mpq_init(R);
  inited = 1;
  mpq_set_si(R, a,b);
}


rational::rational(mpq_t q) {
  mpq_init(R);
  inited = 1;
  mpq_set(R, q);
}


rational::rational(const rational& other) {
  mpq_init(R);
  mpq_set(R, other.R);
}

rational::~rational() {
  mpq_clear(R);
}


void rational::get_mpq(mpq_t q) {
  mpq_set(q,R);
}

rational& rational::operator=(const rational& rhs) {
  if (this == &rhs) {
    return *this;
  }
  mpq_set(R, rhs.R);
  return *this;
}
    


int rational::d() {
  return mpz_get_si(mpq_denref(R));;
}

int rational::n() {
  return mpz_get_si(mpq_numref(R));;
}

double rational::get_d() {
  return mpq_get_d(R);
}



void rational::canonicalize() {
  /*
  int sign = (denom < 0 ? -1 : 1);
  num *= sign;
  denom *= sign;
  sign = (num < 0 ? -1 : 1);
  num *= sign;
  int g = gcd(num,denom);
  num /= g;
  denom /= g;
  num *= sign;
  */
  mpq_canonicalize(R);
}


rational rational::add(rational other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_add(temp, R, other.R);
  rational r = rational(temp);
  mpq_clear(temp);
  return r;
}
    
    
rational rational::div(rational other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_div(temp, R, other.R);
  rational r = rational(temp);
  mpq_clear(temp);
  return r;
}

rational operator+(rational first, rational other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_add(temp, first.R, other.R);
  rational r = rational(temp);
  mpq_clear(temp);
  return r;
}

rational rational::operator-(rational& other){
  mpq_t temp;
  mpq_init(temp);
  mpq_sub(temp, R, other.R);
  rational r = rational(temp);
  mpq_clear(temp);
  return r;
}

rational rational::operator-(){
  mpq_t temp;
  mpq_init(temp);
  mpq_neg(temp, R);
  rational r = rational(temp);
  mpq_clear(temp);
  return r;
}


rational operator/(rational first, rational other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_div(temp, first.R, other.R);
  rational r = rational(temp);
  mpq_clear(temp);
  return r;
}

rational operator/(rational first, int other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_set_si(temp, other, 1);
  mpq_div(temp, first.R, temp);
  rational r = rational(temp);
  mpq_clear(temp);
  return r;
}


rational operator*(rational first, rational other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_mul(temp, first.R, other.R);
  rational r = rational(temp);
  mpq_clear(temp);
  return r;
}

bool operator<(rational first, rational other) {
  return (mpq_cmp(first.R, other.R) == -1);
}

bool operator==(rational first, rational other) {
  return (mpq_cmp(first.R, other.R)==0);
}


ostream& operator<<(ostream& os, rational r) {
  if (mpz_cmp_si(mpq_denref(r.R), 1) == 0) {
    os << mpz_get_si(mpq_numref(r.R));
  } else { 
    os << mpz_get_si(mpq_numref(r.R)) << "/" << mpz_get_si(mpq_denref(r.R));
  }
  return os;
}


/*****************************************************************************/
/* Helper functions                                                          */
/*****************************************************************************/

rational cont_frac_value(vector<int> a) {
  int lenA = a.size();
  rational answer = rational(a[lenA-1], 1);
  int i;
  for (i=lenA-2; i>=0; i--) {
    answer = rational(a[i], 1).add( rational(1,1).div(answer));
  }
  return answer;
}

rational approxRat(double a) {
  vector<int> As(0);
  //cout << "Called with " << a << "\n"; fflush(stdout);
  As.push_back((int)floor(a));
  //cout << "Floor: " << (int)floor(a) << "\n";
  double currentRem = a - (double)As[As.size()-1];
  rational currentR = rational(As[0], 1);
  while (fabs(currentR.get_d() - a) > 0.00000000001) {
    //cout << "last conv element: " << As[As.size()-1] << " and current Rem: "<< currentRem << "\n";
    As.push_back((int)floor(1.0/currentRem));
    currentRem = (1.0/currentRem) - (double)As[As.size()-1];
    currentR =  cont_frac_value(As);
  }
  //cout << "Done!  Returning " << currentR.n() << "/" << currentR.d() << "\n";
  return currentR;
} 
