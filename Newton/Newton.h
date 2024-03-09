// Newton.h
// Author: zz1918@nyu.edu

#include "Newton.cpp"

template<int dim, int rank>
class Newton;

template<int dim, int rank>
ostream& operator<<(ostream& os, Newton<rank, dim> N);