// Newton.cpp : This file solves non-linear equations
// by Newton method with backtracking line search.
// Author: zz1918@nyu.edu

#include<iostream>
#include<cmath>
#include<Eigen/Dense>
using namespace std;
using namespace Eigen;

template<int dim, int rank>
class Newton
{
	typedef Matrix<double, rank, 1>(*Func)(Matrix<double, dim, 1>);
	typedef Matrix<double, rank, dim>(*Jacob)(Matrix<double, dim, 1>);
	Matrix<double, dim - rank, dim> nulls;			// Newton direction requires to be perpendicular to the nulls.
	Func function;
	Jacob jacobian;
	Matrix<double, dim, 1> initial;
	Matrix<double, dim, 1> result;					// Please do not directly change result, use move_result instead.
	int step = 0;
	double epsilon = 1e-8;
	double range = 10;
	Matrix<double, rank, 1> value;
	Matrix<double, rank, dim> jacobi;
	static Matrix<double, rank, 1> ZeroF(Matrix<double, dim, 1> X)
	{
		return Matrix<double, rank, 1>::Zero();
	}
	static Matrix<double, rank, dim> ZeroJ(Matrix<double, dim, 1> X)
	{
		return Matrix<double, rank, dim>::Zero();
	}
	void move_result(Matrix<double, dim, 1> delta)
	{
		result += delta;
		value = function(result);
		jacobi = jacobian(result);
	}
public:
	Newton()
	{
		function = ZeroF;
		jacobian = ZeroJ;
	}
	Newton(Func f, Jacob j)
	{
		function = f;
		jacobian = j;
	}
	Matrix<double, dim, 1> full_value()
	{
		Matrix<double, rank, 1> V = value;
		Matrix<double, dim, 1> U;
		for (Index i = 0; i < rank; ++i)
			for (Index j = 0; j < 1; ++j)
				U(i, j) = V(i, j);
		for (Index i = 0; i < dim - rank; ++i)
			for (Index j = 0; j < 1; ++j)
				U(i + rank, j) = 0.0;
		return U;
	}
	Matrix<double, dim, dim> full_jacobi()
	{
		Matrix<double, rank, dim> J = jacobi;
		Matrix<double, dim, dim> K;
		for (Index i = 0; i < rank; ++i)
			for (Index j = 0; j < dim; ++j)
				K(i, j) = J(i, j);
		for (Index i = 0; i < dim - rank; ++i)
			for (Index j = 0; j < dim; ++j)
				K(i + rank, j) = nulls(i, j);
		return K;
	}
	bool ends()
	{
		double fnorm = value.norm();
		return isnan(fnorm) || fnorm < epsilon || (result - initial).norm() > range;
	}
	bool solved()
	{
		return value.norm() < epsilon;
	}
	int steps()
	{
		return step;
	}
	void set_ini(Matrix<double, dim, 1> X)
	{
		initial = X;
		result = X;
		value = function(result);
		jacobi = jacobian(result);
		step = 0;
	}
	void set_ini(double x)
	{
		assert(dim == 1);
		Matrix<double, dim, 1> X;
		X(0) = x;
		set_ini(X);
	}
	void set_null(Matrix<double, dim - rank, dim> N)
	{
		nulls = N;
	}
	void set_epsilon(double e)
	{
		epsilon = e;
	}
	void set_range(double r)
	{
		range = r;
	}
	void next_step(double t = 1.0, bool backtrack = true, double alpha = 0.5, double beta = 1.0)
	{												// Find direction and apply backtracking line search
		if (ends())
			return;
		Matrix<double, dim, dim> K = full_jacobi();
		Matrix<double, dim, 1> delta;
		if (K.determinant() != 0)
			delta = -K.inverse() * full_value();
		else
		{
			FullPivLU<MatrixXd> lu(K);
			MatrixXd K_null_space = lu.kernel();
			delta = -K_null_space.col(0).normalized();
		}
		if (backtrack)
		{
			double new_norm = function(result + t * delta).squaredNorm();
			double t_epsilon = pow(epsilon, dim * 1.0 / rank);
			while (t > t_epsilon && (new_norm * beta > value.squaredNorm() || isnan(new_norm)))
			{
				t *= alpha;
				new_norm = function(result + t * delta).squaredNorm();
			}
			if (t <= t_epsilon)
				cout << "No available backtracking move. Possibly no root." << endl;
		}
		move_result(t * delta);
		++step;
	}
	void solve(int stop = 20, double t = 1.0, bool backtrack = true, double alpha = 0.5, double beta = 1.0)
	{
		while (!ends() && step < stop)
			next_step(t, backtrack, alpha, beta);
	}
	void view(ostream& os = cout)
	{
		if (dim > rank)
			os << "Null space:" << nulls << endl;
		os << "Initial:" << initial.transpose() << endl;
		os << "Step: " << step << endl;
		os << "Now:" << result.transpose() << endl;
		os << "Now_value:" << value.transpose() << endl;
		if ((result - initial).norm() > range)
			os << "Out of range!" << endl;
	}
	Matrix<double, dim, 1> solution()
	{
		return result;
	}
};

template<int dim, int rank>
ostream& operator<<(ostream& os, Newton<rank, dim> N)
{
	N.view(os);
	return os;
}