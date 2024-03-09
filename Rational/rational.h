// rational.h
// Author: zz1918@nyu.edu

#include "rational.cpp"

class Rational;
ostream& operator<<(ostream& os, Rational* r);
Rational* make_rational(string expression, vector<char>* a);
Rational* make_rational(string expression, string s);
Rational* sum(Rational* r1, Rational* r2);
Rational* diff(Rational* r1, Rational* r2);
Rational* prod(Rational* r1, Rational* r2);
Rational* quot(Rational* r1, Rational* r2);
Rational* sum(Rational* r1, Rational* r2, Rational* r3);
Rational* prod(Rational* r1, Rational* r2, Rational* r3);