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
    Rational(int a, int b);
    Rational(mpq_t q);
    Rational(const Rational& other);
    ~Rational();
    void get_mpq(mpq_t q);
    Rational& operator=(const Rational& rhs);
    double get_d();
    Rational add(Rational other);
    Rational div(Rational other);
    friend Rational operator+(Rational first,  Rational other);
    Rational operator-(Rational& other);
    Rational operator-();
    friend Rational operator/(Rational first, Rational other);
    friend Rational operator/(Rational first, int other);
    friend Rational operator*(Rational first, Rational other);
    friend bool operator==(Rational first, Rational other);
    friend bool operator<(Rational first, Rational other);
    friend std::ostream& operator<<(std::ostream& os, Rational r);
    void canonicalize();
    int d();
    int n();
};


Rational cont_frac_value(std::vector<int> a);
Rational approxRat(double a);
Rational approxRat_be_nice(double a);
int lcm(int a, int b);
int gcd(int a, int b);


#endif

