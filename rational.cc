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
Rational::Rational() {
  mpq_init(R);
  inited = 1;
  mpq_set_si(R,0,1);
}

Rational::Rational(int r) {
  mpq_init(R);
  inited = 1;
  mpq_set_si(R,r,1);
}  

Rational::Rational(int a, int b) {
  mpq_init(R);
  inited = 1;
  mpq_set_si(R, a,b);
}

Rational::Rational(mpq_t q) {
  mpq_init(R);
  inited = 1;
  mpq_set(R, q);
}

Rational::Rational(const Rational& other) {
  mpq_init(R);
  mpq_set(R, other.R);
}

Rational& Rational::operator=(const Rational& rhs) {
  if (this == &rhs) {
    return *this;
  }
  mpq_set(R, rhs.R);
  return *this;
}

Rational::~Rational() {
  mpq_clear(R);
}

void Rational::get_mpq(mpq_t q) {
  mpq_set(q,R);
}

double Rational::get_d() {
  return mpq_get_d(R);
}

void Rational::canonicalize() {
  mpq_canonicalize(R);
}

int Rational::d() {
  return mpz_get_si(mpq_denref(R));;
}

int Rational::n() {
  return mpz_get_si(mpq_numref(R));;
}

Rational Rational::add(const Rational& other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_add(temp, R, other.R);
  Rational r = Rational(temp);
  mpq_clear(temp);
  return r;
}
    
    
Rational Rational::div(const Rational& other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_div(temp, R, other.R);
  Rational r = Rational(temp);
  mpq_clear(temp);
  return r;
}

Rational Rational::operator+(const Rational& other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_add(temp, R, other.R);
  Rational r = Rational(temp);
  mpq_clear(temp);
  return r;
}

Rational Rational::operator+(int other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_set_si(temp, other, 1);
  mpq_add(temp, R, temp);
  Rational r = Rational(temp);
  mpq_clear(temp);
  return r;
}

Rational Rational::operator-(const Rational& other){
  mpq_t temp;
  mpq_init(temp);
  mpq_sub(temp, R, other.R);
  Rational r = Rational(temp);
  mpq_clear(temp);
  return r;
}

Rational Rational::operator-(int other){
  mpq_t temp;
  mpq_init(temp);
  mpq_set_si(temp, other, 1);
  mpq_sub(temp, R, temp);
  Rational r = Rational(temp);
  mpq_clear(temp);
  return r;
}

Rational Rational::operator-(){
  mpq_t temp;
  mpq_init(temp);
  mpq_neg(temp, R);
  Rational r = Rational(temp);
  mpq_clear(temp);
  return r;
}


Rational Rational::operator/(const Rational& other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_div(temp, R, other.R);
  Rational r = Rational(temp);
  mpq_clear(temp);
  return r;
}

Rational Rational::operator/(int other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_set_si(temp, other, 1);
  mpq_div(temp, R, temp);
  Rational r = Rational(temp);
  mpq_clear(temp);
  return r;
}

Rational Rational::operator*(const Rational& other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_mul(temp, R, other.R);
  Rational r = Rational(temp);
  mpq_clear(temp);
  return r;
}

Rational Rational::operator*(int other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_set_si(temp, other, 1);
  mpq_mul(temp, R, temp);
  Rational r = Rational(temp);
  mpq_clear(temp);
  return r;
}


Rational& Rational::operator+=(const Rational& other) {
  mpq_add(R, R, other.R);
  return *this;
}

Rational& Rational::operator+=(int other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_set_si(temp, other, 1);
  mpq_add(R, R, temp);
  mpq_clear(temp);
  return *this;
}

Rational& Rational::operator*=(const Rational& other) {
  mpq_mul(R, R, other.R);
  return *this;
}

Rational& Rational::operator*=(int other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_set_si(temp, other, 1);
  mpq_mul(R, R, temp);
  mpq_clear(temp);
  return *this;
}

Rational& Rational::operator/=(const Rational& other) {
  mpq_div(R, R, other.R);
  return *this;
}

Rational& Rational::operator/=(int other) {
  mpq_t temp;
  mpq_init(temp);
  mpq_set_si(temp, other, 1);
  mpq_div(R, R, temp);
  mpq_clear(temp);
  return *this;
}


bool Rational::operator<(const Rational& other) {
  return (mpq_cmp(R, other.R) < 0);
}

bool Rational::operator>(const Rational& other) {
  return (mpq_cmp(R, other.R) > 0);
}

bool Rational::operator<(int other) {
  return (mpq_cmp_si(R, other, 1) < 0);
}

bool Rational::operator>(int other) {
  return (mpq_cmp_si(R, other, 1) > 0);
}

bool Rational::operator==(const Rational& other) {
  return (mpq_cmp(R, other.R)==0);
}

bool Rational::operator==(int other) {
  return (mpq_cmp_si(R, other, 1)==0);
}


ostream& operator<<(ostream& os, Rational r) {
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

Rational cont_frac_value(vector<int> a) {
  int lenA = a.size();
  Rational answer = Rational(a[lenA-1], 1);
  int i;
  for (i=lenA-2; i>=0; i--) {
    answer = Rational(a[i], 1).add( Rational(1,1).div(answer));
  }
  return answer;
}

Rational approx_rat(double b) {
  double a;
  vector<int> As(0);
  a = b;
  //cout << "Called with " << a << "\n"; fflush(stdout);
  As.push_back((int)floor(a));
  //cout << "Floor: " << (int)floor(a) << "\n";
  double currentRem = a - (double)As[As.size()-1];
  Rational currentR = Rational(As[0], 1);
  while (fabs(currentR.get_d() - a) > 0.00000001) {
    //cout << "last conv element: " << As[As.size()-1] << " and current Rem: "<< currentRem << "\n";
    As.push_back((int)floor(1.0/currentRem));
    currentRem = (1.0/currentRem) - (double)As[As.size()-1];
    currentR =  cont_frac_value(As);
  }
  //cout << "Done!  Returning " << currentR.n() << "/" << currentR.d() << "\n";
  return currentR;
} 

Rational approx_rat_be_nice(double b) {
  double a;
  vector<int> As(0);
  a = b - 0.000000001;
  //cout << "Called with " << a << "\n"; fflush(stdout);
  As.push_back((int)floor(a));
  //cout << "Floor: " << (int)floor(a) << "\n";
  double currentRem = a - (double)As[As.size()-1];
  Rational currentR = Rational(As[0], 1);
  while (fabs(currentR.get_d() - a) > 0.00000001) {
    //cout << "last conv element: " << As[As.size()-1] << " and current Rem: "<< currentRem << "\n";
    As.push_back((int)floor(1.0/currentRem));
    currentRem = (1.0/currentRem) - (double)As[As.size()-1];
    currentR =  cont_frac_value(As);
  }
  //cout << "Done!  Returning " << currentR.n() << "/" << currentR.d() << "\n";
  return currentR;
} 

