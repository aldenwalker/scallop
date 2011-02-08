#ifndef __rational__
#define __rational__

#include <iostream>
#include <vector>
#include <gmp.h>

/******************************************************************************/
/** rational class  hiding gmp                                                */
/******************************************************************************/
class rational {
  private:
    int inited;
    mpq_t R;
  public:
    rational();
    rational(int a, int b);
    rational(mpq_t q);
    rational(const rational& other);
    ~rational();
    void get_mpq(mpq_t q);
    rational& operator=(const rational& rhs);
    double get_d();
    rational add(rational other);
    rational div(rational other);
    friend rational operator+(rational first,  rational other);
    rational operator-(rational& other);
    rational operator-();
    friend rational operator/(rational first, rational other);
    friend rational operator/(rational first, int other);
    friend rational operator*(rational first, rational other);
    friend bool operator==(rational first, rational other);
    friend bool operator<(rational first, rational other);
    friend std::ostream& operator<<(std::ostream& os, rational r);
    void canonicalize();
    int d();
    int n();
};


rational cont_frac_value(std::vector<int> a);
rational approxRat(double a);
int lcm(int a, int b);
int gcd(int a, int b);


#endif

