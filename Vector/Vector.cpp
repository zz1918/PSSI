// Vector.cpp : This file transfers between vector<Scarlar> and Vector<Scarlar>
// Author: zz1918@nyu.edu

#include<vector>
#include<Eigen/Dense>
using namespace std;
using namespace Eigen;

template<typename Scarlar>
vector<Scarlar> to_vector(Matrix<Scarlar,-1, 1> V)
{
	vector<Scarlar> v;
	for (int i = 0; i < V.size(); ++i)
		v.push_back(V(i));
	return v;
}
template<typename Scarlar>
Matrix<Scarlar, -1, 1> to_Vector(vector<Scarlar> v)
{
	Matrix<Scarlar, -1, 1> V;
	V.resize(v.size());
	for (int i = 0; i < v.size(); ++i)
		V(i)=v[i];
	return V;
}