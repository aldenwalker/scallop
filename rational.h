#ifndef rational_H
#define rational_H

#include <iostream>
#include <vector>
#include <gmp.h>

/******************************************************************************/
/** rational class  hiding gmp                                                */
/******************************************************************************/
class Rational {
  private:
    int inited;
    mpq_t R;
  public:
    Rational();
    Rational(int r);
    Rational(int a, int b);
    Rational(mpq_t q);
    Rational(const Rational& other);
    Rational& operator=(const Rational& rhs);
    ~Rational();
    void get_mpq(mpq_t q);
    double get_d();    
    void canonicalize();
    int d();
    int n();
    Rational add(const Rational& other);
    Rational div(const Rational& other);
    Rational operator+(const Rational& other);
    Rational operator+(int other);
    Rational operator-(const Rational& other);
    Rational operator-(int other);
    Rational operator-();
    Rational operator/(const Rational& other);
    Rational operator/(int other);
    Rational operator*(const Rational& other);
    Rational operator*(int other);
    Rational& operator+=(const Rational& other);
    Rational& operator+=(int other);
    Rational& operator*=(const Rational& other);
    Rational& operator*=(int other);
    Rational& operator/=(const Rational& other);
    Rational& operator/=(int other);
    bool operator<(const Rational& other);
    bool operator>(const Rational& other);
    bool operator<(int other);
    bool operator>(int other);
    bool operator==(const Rational& other);
    bool operator==(int other);
    
    
    friend std::ostream& operator<<(std::ostream& os, Rational r);

};


Rational cont_frac_value(std::vector<int> a);
Rational approx_rat(double a);
Rational approx_rat_be_nice(double a);
int lcm(int a, int b);
int gcd(int a, int b);


#endif

