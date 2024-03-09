// Vector.h
// Author: zz1918@nyu.edu

template<typename Scarlar>
vector<Scarlar> to_vector(Matrix<Scarlar, -1, 1> V);

template<typename Scarlar>
Matrix<Scarlar, -1, 1> to_Vector(vector<Scarlar> v);

#define tVi to_Vector<int>
#define tVd to_Vector<double>
#define tvi to_vector<int>
#define tvd to_vector<double>