// PPI.cpp

#include<iostream>
#include<iomanip>
#include<random>
#include<ctime>
#include<Eigen/Dense>
#include<interval.h>
#include<rational.h>
#include<json.h>
#include<Vector.h>
#include<Newton.h>
using namespace std;
using namespace Eigen;
#define Id intervald
#define UE 1e-12												// An universal epsilon.
#define WE 1e-3													// Epsilon for weakly monotonicity.

void Random_ini()
{
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
}
double Random(double a, double b)
{
    return a + (b - a) * (rand() / double(RAND_MAX));
}

/*{
	"initial": ["0.1111106996e-4", ".2454545146", "0.1111106996e-4", ".2454545146"],
	"target": [".9988892998", ".2460661205", ".9988892998", ".2460661205"],
	"epsilon": "0.1e-3",
	"surface1": ["0.4*v*(v^2+1.499999999)", "-0.2*u*(u+1.500000000)*(u-3.000000000)", "2.3*v*(v^2-1.956521739*v+1.304347826)"],
	"surface2": ["0.4*t*(t^2+1.499999999)", "-0.2*s*(s+1.500000000)*(s-3.000000000)", "-3.0*t+1.0+4.5*t^2-2.4*t^3+0.1*s^3*t^3"],
	"vars": ["u", "v", "s", "t"],
	"stick_size": ".1",
	"target3d": [".1531879582", "0.1000000000e-4", ".4992605938"]
}*/


Rational* solving[3];				// Use static member as the target of Newton method

Vector3d solve_f(Vector4d X)
{
	vector<double> x;
	for (int i = 0; i < 4; ++i)
		x.push_back(X(i));
	return Vector3d(solving[0]->value<double>(x), solving[1]->value<double>(x), solving[2]->value<double>(x));
}

Matrix<double, 3, 4> solve_j(Vector4d X)
{
	Matrix<double, 3, 4> J;
	vector<double> x;
	for (int i = 0; i < 4; ++i)
		x.push_back(X(i));
	J << solving[0]->deriv<double>(x, 0), solving[0]->deriv<double>(x, 1), solving[0]->deriv<double>(x, 2), solving[0]->deriv<double>(x, 3),
		solving[1]->deriv<double>(x, 0), solving[1]->deriv<double>(x, 1), solving[1]->deriv<double>(x, 2), solving[1]->deriv<double>(x, 3),
		solving[2]->deriv<double>(x, 0), solving[2]->deriv<double>(x, 1), solving[2]->deriv<double>(x, 2), solving[2]->deriv<double>(x, 3);
	return J;
}

class trial
{
	Vector4d initial;
	Vector4d target;
	double epsilon;
	vector<char>* alphabets;
	Rational* Surface[2][3];
	Rational* Diff[3];							// The W function.
	Rational* JDiff[3][4];						// Jacobian function of W function (JW).
	Rational* DT[4];							// Determinant function of T matrix given by erasing the i-th column from JW.
	string vars;
	double stick_size;
	Vector3d target3d;
public:
	trial() {}
	trial(json data)
	{
		epsilon = stod(string(data["epsilon"]));
		stick_size = stod(string(data["stick_size"]));
		initial = read_vec4(data["initial"]);
		target = read_vec4(data["target"]);
		target3d = read_vec3(data["target3d"]);
		vars = "";
		for (int i = 0; i < 4; ++i)
			vars += string(data["vars"][i])[0];
		alphabets = make_alphabets(vars);
		for (int i = 0; i < 3; ++i)
		{
			Surface[0][i] = make_rational(string(data["surface1"][i]), alphabets);
			Surface[1][i] = make_rational(string(data["surface2"][i]), alphabets);
			Diff[i] = make_rational("(" + string(data["surface1"][i]) + ")-(" + string(data["surface2"][i]) + ")", alphabets);
		}
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 4; ++j)
				JDiff[i][j] = Diff[i]->deriv(j);
		const int emits[4][3] = { {1,2,3},{0,2,3},{0,1,3},{0,1,2} };
		for (int i = 0; i < 4; ++i)
		{
			Rational* DTi = new Rational(alphabets, 0.0);
			for (int j = 0; j < 3; ++j)
			{
				Rational* DTij = new Rational(alphabets, 1.0);
				for (int k = 0; k < 3; ++k)
					DTij = prod(JDiff[k][emits[i][(j + k) % 3]], DTij);
				DTi = sum(DTi, DTij);
			}
			for (int j = 0; j < 3; ++j)
			{
				Rational* DTij = new Rational(alphabets, 1.0);
				for (int k = 0; k < 3; ++k)
					DTij = prod(JDiff[k][emits[i][(j - k + 3) % 3]], DTij);
				DTi = diff(DTi, DTij);
			}
			DT[i] = DTi;
		}
	}
	Rational* surface(int i, int j)
	{
		assert(i == 0 || i == 1);
		assert(j == 0 || j == 1 || j == 2);
		return Surface[i][j];
	}
	Rational* the_diff(int i)
	{
		assert(i == 0 || i == 1 || i == 2);
		return Diff[i];
	}
	Vector4d find_root(Vector4d ini, Vector4d null)
	{
		solving[0] = Diff[0];
		solving[1] = Diff[1];
		solving[2] = Diff[2];
		Newton<4, 3> N(solve_f, solve_j);
		N.set_ini(ini);
		N.set_null(null.transpose());
		cout << setprecision(10);
		N.solve();
		cout << "Solution: " << N.solution().transpose() << endl;
		cout << "Solution value: " << solve_f(N.solution()).transpose() << endl;
		return N.solution();
		//return Vector4d::Zero();
	}
	MatrixXd Jacob(Vector4d V)
	{
		MatrixXd J(3, 4);
		vector<double> x;
		for (int i = 0; i < 4; ++i)
			x.push_back(V(i));
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 4; ++j)
				J(i, j) = Diff[i]->deriv<double>(x, j);
		return J;
	}
	RowVector4d Deter(Vector4d V)
	{
		MatrixXd J = Jacob(V);
		RowVector4d D;
		Matrix3d Blocks;
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 3; ++j)
				for (int k = 0; k < 3; ++k)
					Blocks(j, k) = ((k < i) ? J(j, k) : J(j, k + 1));
			D(i) = Blocks.determinant() * pow(-1, i + 1);
		}
		return D;
	}
	MatrixId ValueInterval(MatrixId I)			// I is interval of Vector4d, return the interval of Value
	{
		assert(I.size() == 4);
		MatrixId V(3, 1);
		vector<Id> x;
		for (int i = 0; i < 4; ++i)
			x.push_back(I(i));
		for (int i = 0; i < 3; ++i)
			V.set(i, 1, Diff[i]->value<Id>(x));
		return V;
	}
	MatrixId FaceValueInterval(MatrixId I, Vector4d Q, int t)
	{
		assert(I.size() == 4);
		MatrixId V(3, 1);
		vector<Id> x;
		for (int i = 0; i < 4; ++i)
			x.push_back((i == t) ? Id(Q(i), Q(i)) : I(i));
		for (int i = 0; i < 3; ++i)
			V.set(i, 1, Diff[i]->value<Id>(x));
		return V;
	}
	MatrixId JacobInterval(MatrixId I)			// I is interval of Vector4d, return the interval of Jacobian
	{
		assert(I.size() == 4);
		MatrixId J(3, 4);
		vector<Id> x;
		for (int i = 0; i < 4; ++i)
			x.push_back(I(i));
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 4; ++j)
				J.set(i, j, Diff[i]->deriv<Id>(x, j));
		return J;
	}
	MatrixId DeterInterval(MatrixId I, bool JI = false)		// I is interval of Vector4d when JI=false, return the 4 intervals of determinants
	{
		MatrixId J;
		if (JI)
			J = I;
		else
			J = JacobInterval(I);
		MatrixId D(1, 4);
		J.add_row();
		for (int i = 0; i < 4; ++i)
			D.set(0, i, J.coblock(3, i).determinant() * pow(-1, i + 1));
		return D;
	}
	MatrixId FaceJacobInterval(MatrixId I, Vector4d Q, int t)	// I is I is interval of Vector4d.
	{
		assert(I.size() == 4);
		assert(t < 4);
		MatrixId J(3, 3);
		vector<Id> x;
		for (int i = 0; i < 4; ++i)
			x.push_back((i == t) ? Id(Q(i), Q(i)) : I(i));
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				J.set(i, j, (j < t) ? Diff[i]->deriv<Id>(x, j) : Diff[i]->deriv<Id>(x, j + 1));
		return J;
	}
	MatrixId FaceDeterInterval(MatrixId I, Vector4d Q)
	{
		MatrixId D(1, 4);
		for (int i = 0; i < 4; ++i)
		{
			MatrixId J = FaceJacobInterval(I, Q, i);
			D.set(0, i, J.determinant() * pow(-1, i + 1));
		}
		return D;
	}
	bool Strong_Monotonic(MatrixId I, bool DI = false)			// Strongly monotonic, the first part of Lemma 5.
	{
		MatrixId DTI;
		if (DI)													// If DTI is calculated before the function, use it directly.
			DTI = I;
		else
			DTI = DeterInterval(I);
		for (int i = 0; i < 4; ++i)
			if (DTI(i).contains(0.0))
				return false;
		return true;
	}
	bool Weak_Monotonic(Vector4d P, Vector4d Q)					// Weakly monotonic, the second part of Lemma 5.
	{
		// Step a: exists i such that 0\notin \det(T_i(PQ)).
		MatrixId DTPQ = DeterInterval(MatrixId(P, Q));
		bool exist_i = false;
		for (int i = 0; i < 4; ++i)
			if (DTPQ(i).ncontains(0.0))
				exist_i = true;
		if (!exist_i)
			return false;
		// Step b: B_{P,Q'}^{q_i} if then.
		RowVector4d DTP = Deter(P), DTQ = Deter(Q);
		Vector4d Epsilon;
		for (int i = 0; i < 4; ++i)
			if (abs(DTQ(i)) < UE)
				Epsilon(i) = 0;
			else
				Epsilon(i) = WE;
		Vector4d QQ = Q + Epsilon.cwiseProduct(Q - P);						// This is Q'.
		bool if_then = true;
		MatrixId DTPQQq = FaceDeterInterval(MatrixId(P, QQ), Q);
		for (int i = 0; i < 4; ++i)
			if (abs(DTQ(i)) >= UE && DTPQQq(i).contains(0.0))				// DTPQQq(i) needed to be implemented.
				if_then = false;
		if (!if_then)
			return false;
		// Step c: some j if then.
		if_then = true;
		for (int i = 0; i < 4; ++i)
			if (abs(DTQ(i)) <= UE)
			{
				MatrixId VIPQQq = FaceValueInterval(MatrixId(P, QQ), Q, i);
				if (VIPQQq.ncontains(Vector3d(0, 0, 0)))
					if_then = false;
			}
		if (!if_then)
			return false;
		// Step d: determination not contain 0 if then.
		MatrixId JWPQQ = JacobInterval(MatrixId(P, QQ));			// J(W)(\B_{P,Q'})
		MatrixId DTPQQ = DeterInterval(JWPQQ, true);				// \det(T_i(\B_{P,Q'}))
		if_then = true;
		JWPQQ.add_row();
		vector<Id> x;												// MatrixId of \B_{P,Q'} in variable form
		for (int i = 0; i < 4; ++i)
			x.push_back(Id(P(i), QQ(i)));
		for (int i = 0; i < 4; ++i)
			if (DTPQQ(i).contains(0.0))
			{
				MatrixId JWDTPQQ = JWPQQ;							// J(W,\det(T_i))(\B_{P,Q'})
				Rational* DTi = DT[i];								// \det(T_i);
				for (int j = 0; j < 4; ++j)
					JWDTPQQ.set(3, j, DT[i]->deriv<Id>(x, j));
				if (JWDTPQQ.determinant().contains(0.0))
					if_then = false;
			}
		if (!if_then)
			return false;
		return true;
	}
	bool Strongly_Monotonic(Vector4d P, Vector4d Q)						// Lemma 5.
	{
		if (Strong_Monotonic(MatrixId(P, Q)))
			return true;
		return Weak_Monotonic(P, Q);
	}
	vector<int> IsTerminate(Vector4d P, Vector4d Q)						// Algorithm 1.
	{
		// Step 1: initialization of S.
		vector<int> S;							// {0,1,2,3} instead of {1,2,3,4}.
		vector<int> Empty;						// This is an empty set.
		for (int i = 0; i < 4; ++i)
			S.push_back(i);
		// Step 2: quick test is used to exclude curve segments that are definitely not strongly monotonic.
		RowVector4d DTP = Deter(P), DTQ = Deter(Q);
		bool test_pass = true;
		for (int i = 0; i < 4; ++i)
			if (DTP(i) * DTQ(i) >= 0 && abs(DTQ(i)) >= UE)
				test_pass = false;
		if (test_pass)
			return S;
		// Step 3: check whether it is a strongly monotonic curve segment using Lemma 5.
		if (Strongly_Monotonic(P,Q))
			return Empty;
		// Step 4:
		MatrixId DTJ = DeterInterval(MatrixId(P,Q));
		S.clear();
		for (int j = 0; j < 4; ++j)
			if (DTJ(j).contains(0.0))					// Instead of extract, we clear S and push the unextracted instead.
				S.push_back(j);
		return S;
	}
	void view(ostream& os = cout)
	{
		os << "initial: " << initial.transpose() << endl;
		os << "target: " << target.transpose() << endl;
		os << "epsilon: " << epsilon << endl;
		for (int i = 0; i < 2; ++i)
		{
			os << "Surface" << i << ":" << endl;
			for (int j = 0; j < 3; ++j)
				os << "          " << Surface[i][j] << endl;
		}
		os << "Diff:" << endl;
		for (int j = 0; j < 3; ++j)
			os << "     " << Diff[j] << endl;
		os << "vars: " << vars << endl;
		os << "stick_size: " << stick_size << endl;
		os << "target3d: " << target3d.transpose() << endl;
	}
};

ostream& operator<<(ostream& os, trial t)
{
	t.view(os);
	return os;
}

int main(int argc, char* argv[])
{
    const char filename[1000] = "C:/Users/zhaoq/Desktop/Study/My_Geometry_Library/Curve_Tracing/Input/inputs.js";
    trial t(read_json(filename));
	//t.find_root(Vector4d(0.5, 0.5, 0.5, 0.5), Vector4d(1.0, 1.0, 1.0, 1.0));
	cout << t.JacobInterval(MatrixId(Vector4d(0.5, 0.5, 0.5, 0.5), Vector4d(1.0, 1.0, 1.0, 1.0))) << endl;
	cout << t.DeterInterval(MatrixId(Vector4d(0.5, 0.5, 0.5, 0.5), Vector4d(1.0, 1.0, 1.0, 1.0))) << endl;
	cout << t.Jacob(Vector4d(0.5, 0.5, 0.5, 0.5)) << endl;
	cout << t.Deter(Vector4d(0.5, 0.5, 0.5, 0.5)) << endl;
	//cout << "Hello world!" << endl;
	return 0;
}