// PPI.cpp

#include<iostream>
#include<iomanip>
#include<ctime>
#include<random>
#include<vector>
#include<stack>
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
#pragma warning(disable: 4244)
#pragma warning(disable: 4267)

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
	double epsilon;
	vector<char>* alphabets;
	vector<Vector4d> criticals;
	vector<Vector4d> criticals_3d;
	int critical_size;
	int critical_size_3d;
	double boundary;
	Rational* Surface[2][3];
	Rational* Diff[3];							// The W function.
	Rational* JDiff[3][4];						// Jacobian function of W function (JW).
	Rational* DT[4];							// Determinant function of T matrix given by erasing the i-th column from JW.
	Rational* DTC[3];							// Determinant function of TC matrix.
	string vars;
	Vector3d target3d;
public:
	Vector4d initial;
	Vector4d target;
	double stick_size;
	trial() {}
	trial(json data)
	{
		epsilon = stod(string(data["epsilon"]));
		stick_size = stod(string(data["stick_size"]));
		boundary = stod(string(data["boundary"]));
		initial = read_vec4(data["initial"]);
		target = read_vec4(data["target"]);
		target3d = read_vec3(data["target3d"]);
		critical_size = stoi(string(data["critical_size"]));
		for (int i = 0; i < critical_size; ++i)
			criticals.push_back(read_vec4(data["criticals"][i]));
		critical_size_3d = stoi(string(data["critical_size_3d"]));
		for (int i = 0; i < critical_size_3d; ++i)
			criticals_3d.push_back(read_vec4(data["criticals_3d"][i]));
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
		for (int i = 0; i < 3; ++i)
			DTC[i] = diff(prod(JDiff[i][0], DT[0]), prod(JDiff[i][1], DT[1]));
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
	Vector3d Value(Vector4d V)
	{
		vector<double> x;
		for (int i = 0; i < 4; ++i)
			x.push_back(V(i));
		return Vector3d(Surface[0][0]->value<double>(x), Surface[0][1]->value<double>(x), Surface[0][2]->value<double>(x));
	}
	Vector4d find_root(Vector4d ini, Vector4d null)
	{
		solving[0] = Diff[0];
		solving[1] = Diff[1];
		solving[2] = Diff[2];
		Newton<4, 3> N(solve_f, solve_j);
		N.set_ini(ini);
		N.set_null(null.transpose());
		//cout << setprecision(10);
		N.solve();
		//cout << "Solution: " << N.solution().transpose() << endl;
		//cout << "Solution value: " << solve_f(N.solution()).transpose() << endl;
		return N.solution();
		//return Vector4d::Zero();
	}
	Vector4d tan_dir(Vector4d ini)					// Tangent direction.
	{
		return Deter(ini);
	}
	Vector4d get_next(Vector4d ini,double length)
	{
		Vector4d tau = tan_dir(ini);
		Vector4d Newton_Ini = ini - tau.normalized() * length;
		return find_root(Newton_Ini, tau);
	}
	Vector4d get_last(Vector4d ini, double length)
	{
		Vector4d tau = tan_dir(ini);
		Vector4d Newton_Ini = ini + tau.normalized() * length;
		return find_root(Newton_Ini, tau);
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
	RowVector4d Deter(MatrixXd V)									// TW(V)
	{
		MatrixXd J(3, 4);
		if (V.size() == 4)
			J = Jacob(V);
		else
			J = V;
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
	RowVector3d Deter3(Vector4d V)									// TC(V)
	{
		RowVector3d D;
		MatrixXd J = Jacob(V);
		RowVector4d TW = Deter(J);
		for (int i = 0; i < 3; ++i)
			D(i) = J(i, 0) * TW(0) - J(i, 1) * TW(1);
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
		//cout << "Value finished!" << endl;
		for (int i = 0; i < 3; ++i)
			V.set(i, 0, Diff[i]->value<Id>(x));
		//cout << "Assign finished!" << endl;
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
	MatrixId Deter3Interval(MatrixId I, bool JI = false)									// TC(V)
	{
		MatrixId D(1, 3);
		MatrixId J;
		if (JI)
			J = I;
		else
			J = JacobInterval(I);
		MatrixId TW = DeterInterval(J, true);
		for (int i = 0; i < 3; ++i)
			D.set(0, i, J(i, 0) * TW(0) - J(i, 1) * TW(1));
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

	// Strongly-monotonicity for Algorithm 1.
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
		//cout << "Checking weakly monotonic (" << P.transpose() << ")-(" << Q.transpose() << ")" << endl;

		// Step a: exists i such that 0\notin \det(T_i(PQ)).
		MatrixId DTPQ = DeterInterval(MatrixId(P, Q));
		bool exist_i = false;
		for (int i = 0; i < 4; ++i)
			if (DTPQ(i).ncontains(0.0))
				exist_i = true;
		if (!exist_i)
			return false;

		//cout << "P-Q exsits a monotonic direction." << endl;

		// Step b: B_{P,Q'}^{q_i} if then.
		RowVector4d DTP = Deter(P), DTQ = Deter(Q);
		Vector4d Epsilon;
		for (int i = 0; i < 4; ++i)
			if (abs(DTQ(i)) < UE)
				Epsilon(i) = 0;
			else
				Epsilon(i) = WE;
		Vector4d QQ = Q + Epsilon.cwiseProduct(Q - P);						// This is Q'.
		//cout << "Q' is (" << QQ.transpose() << ")" << endl;
		MatrixId DTPQQq = FaceDeterInterval(MatrixId(P, QQ), Q);
		for (int i = 0; i < 4; ++i)
			if (abs(DTQ(i)) >= UE && DTPQQq(i).contains(0.0))				// DTPQQq(i) needed to be implemented.
				return false;

		//cout << "P-Q' extended is monotonic." << endl;

		// Step c: some j if then.
		for (int i = 0; i < 4; ++i)
			if (abs(DTQ(i)) < UE)
			{
				MatrixId VIPQQq = FaceValueInterval(MatrixId(P, QQ), Q, i);
				if (VIPQQq.contains(Vector3d(0, 0, 0)))
					return false;
			}

		//cout << "P-Q' has no root when Q is a critical." << endl;

		// Step d: determinant not contain 0 if then.
		MatrixId JWPQQ = JacobInterval(MatrixId(P, QQ));			// J(W)(\B_{P,Q'})
		MatrixId DTPQQ = DeterInterval(JWPQQ, true);				// \det(T_i(\B_{P,Q'}))
		JWPQQ.add_row();
		vector<Id> x;												// MatrixId of \B_{P,Q'} in variable form
		for (int i = 0; i < 4; ++i)
			x.push_back(Id(P(i), QQ(i)));
		for (int i = 0; i < 4; ++i)
			if (DTPQQ(i).contains(0.0))
			{
				MatrixId JWDTPQQ = JWPQQ;							// J(W,\det(T_i))(\B_{P,Q'})
				//Rational* DTi = DT[i];								// \det(T_i);
				for (int j = 0; j < 4; ++j)
					JWDTPQQ.set(3, j, DT[i]->deriv<Id>(x, j));
				if (JWDTPQQ.determinant().contains(0.0))
					return false;
			}
		//cout << "P-Q is monotonic with only critical Q." << endl;
		return true;
	}
	bool Strongly_Monotonic(Vector4d P, Vector4d Q)						// Lemma 5.
	{
		return Strong_Monotonic(MatrixId(P, Q)) || Weak_Monotonic(P, Q);
	}
	bool Weak_Monotonic_3d(Vector4d P, Vector4d Q)					// Second part of lemma 7.
	{
		//cout << "Checking 3d weakly monotonic for box (" << P.transpose() << ")-(" << Q.transpose() << ")" << endl;

		// Step 1, DTC(Q)=0 for some i, and DTC(P)!=0 for all i.
		bool bondary_critical = false;
		RowVector3d DTP = Deter3(P), DTQ = Deter3(Q);
		for (int i = 0; i < 3; ++i)
		{
			if (abs(DTP(i)) < UE)
				return false;
			if (abs(DTQ(i) < UE))
				bondary_critical = true;
		}
		if (!bondary_critical)
			return false;
		//cout << "Only Q is a critical point." << endl;

		// Step 2, 0 \notin DTC(PQ) for some i.
		bool exists_monotonic_direct = false;
		MatrixId DTPQ = Deter3Interval(MatrixId(P, Q));
		for (int i = 0; i < 3; ++i)
			if (!DTPQ(i).contains(0.0))
				exists_monotonic_direct = true;
		if (!exists_monotonic_direct)
			return false;
		//cout << "There is a monotonic direction." << endl;

		// Step 3, if DT(Q)!=0 for some i, then 0 \notin DTC(PQ) for the i.
		for (int i = 0; i < 3; ++i)
			if (abs(DTQ(i)) >= UE && DTPQ(i).contains(0.0))
				return false;
		//cout << "Monotonic direction has no criticals." << endl;

		// Step 4, if 0 \in DTC(PQ), then 0 \notin JWDTC(PQ).
		MatrixId JWPQ = JacobInterval(MatrixId(P, Q));					// J(W)(B)
		JWPQ.add_row();
		vector<Id> x;													// MatrixId of B in variable form
		for (int i = 0; i < 4; ++i)
			x.push_back(Id(P(i), Q(i)));
		for (int i = 0; i < 3; ++i)
			if (DTPQ(i).contains(0.0))
			{
				MatrixId JWDTPQ = JWPQ;									// J(W,\det(T_i))(B)
				for (int j = 0; j < 4; ++j)
					JWDTPQ.set(3, j, DTC[i]->deriv<Id>(x, j));
				//cout << JWDTPQ << endl;
				//cout << JWDTPQ.determinant() << endl;
				if (JWDTPQ.determinant().contains(0.0))
					return false;
			}
		//cout << "Q is uniquely critical." << endl;

		return true;
	}

	// Step 3 for Algorithm 2. This function will find the farthest K' that makes K-K' strongly monotonic and also push K' into the list S.
	Vector4d basic_trace(Vector4d K, vector<Vector4d>& S,double l)
	{
		RowVector4d DTK = Deter(K);
		Vector4d KK;									// K'
		bool needed = false;
		bool output = false;
		//output = true;
		for (int i = 0; i < 4; ++i)
			if (abs(DTK(i)) < UE)
				needed = true;
		if (needed)
		{
			do
			{
				KK = get_next(K, l);
				l /= 2;
				if (output)
					cout << "Try box (" << K.transpose() << ")-(" << KK.transpose() << ")" << endl;
			} while (!IsTerminate(KK, K).empty() && l >= UE);
			S.push_back(KK);
			return KK;
		}
		else
			return K;
	}

	// Step 5 for Algorithm 2. This function will find all critical points and insert them into S'.
	void find_critical(Vector4d P, Vector4d Q, vector<Vector4d>& SS)
	{
		MatrixId PQ(P, Q);
		for (int i = 0; i < criticals.size(); ++i)
			if (PQ.contains(criticals[i]))
				SS.push_back(criticals[i]);
	}
	void find_critical_3d(Vector4d P, Vector4d Q, vector<Vector4d>& SS)
	{
		MatrixId PQ(P, Q);
		for (int i = 0; i < criticals_3d.size(); ++i)
			if (PQ.contains(criticals_3d[i]))
				SS.push_back(criticals_3d[i]);
	}

	// Starting from now, the functions are algorithm 1,2,3,etc.
	vector<int> IsTerminate(Vector4d P, Vector4d Q)						// Algorithm 1.
	{
		// Step 1: initialization of S.
		vector<int> S;							// {0,1,2,3} instead of {1,2,3,4}.
		vector<int> Empty;						// This is an empty set.
		for (int i = 0; i < 4; ++i)
			S.push_back(i);
		//cout << "Step 1 finished!" << endl;
		
		// Step 2: quick test is used to exclude curve segments that are definitely not strongly monotonic.
		RowVector4d DTP = Deter(P), DTQ = Deter(Q);
		bool test_pass = true;
		for (int i = 0; i < 4; ++i)
			if (DTP(i) * DTQ(i) >= 0 && abs(DTQ(i)) >= UE)
				test_pass = false;
		if (test_pass)
			return S;
		//cout << "Step 2 finished!" << endl;
		
		// Step 3: check whether it is a strongly monotonic curve segment using Lemma 5.
		if (Strongly_Monotonic(P,Q))
			return Empty;
		//cout << "Step 3 finished!" << endl;

		// Step 4: check the determinant of the Jaccobian of the W with the four determinants.
		MatrixId DTJ = DeterInterval(MatrixId(P,Q));
		S.clear();
		for (int j = 0; j < 4; ++j)
			if (DTJ(j).contains(0.0))					// Instead of extraction, we clear S and push the unextracted instead.
				S.push_back(j);
		//cout << "Step 4 finished!" << endl;
		return S;
	}
	vector<Vector4d> Decompose4D(Vector4d P, Vector4d Q,double l)					// Algorithm 2.
	{
		vector<Vector4d> S;								// S
		vector<Vector4d> SS;							// S'
		Vector4d K1 = P, K2 = Q;
		int t = 0;
		S.push_back(K1);
		bool output = false;
		//output = true;
		while (l > UE && t < 20)
		{
			t++;
			if(output)
				cout << "Cycle " << t << ", now decomposing the box (" << K1.transpose() << ")-(" << K2.transpose() << ")" << endl;

			// Step 2, define S'.
			SS.clear();

			// Step 3, find the non-critical K1 to be the initial point.
			K1 = basic_trace(K1, S, l);

			// Step 4, if K1-K2 is strongly monotonic, then end the algorithm.
			if (IsTerminate(K1, K2).empty())
				S.push_back(K2);
			else
			{
				// Step 5, find all criticals between K1 and K2 and push into SS.
				find_critical(K1, K2, SS);

				// Step 5(1), half the interval to exclude the critical K2.
				if (SS.empty())
				{
					l /= 2;
					K2 = get_next(K1, l);
					if (output)
						cout << "Return to step 2 from step 5(1)" << endl;
					continue;
				}
				else
				{
					// Step 5(2), extract a point H from S'
					while (!SS.empty())
					{
						Vector4d H = SS.back();
						SS.pop_back();
						// Step 5(2)i, find K1-H box.
						if (IsTerminate(K1, H).empty())
							S.push_back(H);
						else
							continue;
						// Step 5(2)ii, check H-K2 box.
						if (IsTerminate(K2, H).empty())
							S.push_back(K2);
						break;
					}
				}
			}
			// Step 6, if S is not empty and the last box is not in P-Q, then return S.
			if (S.size() > 1 && (!MatrixId(P, Q).contains(S[S.size() - 2])) && (P != S[S.size() - 2]) && (!MatrixId(P, Q).contains(S[S.size() - 1]) && (Q != S[S.size() - 1])))
				break;

			// Step 7, trace from K2 to get a new K1,K2.
			K1 = K2;
			K2 = get_next(K1, l);
			if (output)
				cout << "Return to step 2 from step 7" << endl;
		}
		return S;
	}
	vector<pair<Vector4d, Vector4d>> Decompose3D(Vector4d P, Vector4d Q)					// Algorithm 3.
	{
		vector<pair<Vector4d, Vector4d>> S;
		stack<pair<Vector4d, Vector4d>> SS;
		SS.push(make_pair(P, Q));
		int t = 0;
		bool output = false;
		//output = true;
		while (!SS.empty())
		{
			t++;
			// Step 3, extract a box from S'.
			pair<Vector4d, Vector4d> B= SS.top();
			SS.pop();

			if ((B.first - B.second).norm() < UE)
				continue;

			if (output)
				cout << "Checking box (" << B.first.transpose() << ")-(" << B.second.transpose() << ")" << endl;

			// Step 4, if no critical points, insert B into S.
			MatrixId TC = Deter3Interval(MatrixId(B.first, B.second));
			bool is_mono = true;
			for (int i = 0; i < 3; ++i)
				if (TC(i).contains(0.0))
					is_mono = false;
			if (is_mono)
			{
				S.push_back(B);
				continue;
			}

			// Step 5, if satiesfies lemma 7 (2), then insert B into S.
			if (Weak_Monotonic_3d(B.first, B.second) || Weak_Monotonic_3d(B.second, B.first))
			{
				S.push_back(B);
				continue;
			}

			// Step 6, find all criticals.
			vector<Vector4d> SSS;
			find_critical_3d(B.first, B.second, SSS);

			// Step 7, if no critical (find_critical_3d will not include the boundary criticals).
			if (SSS.empty())
			{
				double l = (B.first - B.second).norm();
				Vector4d B_third = get_next(B.first, l / 2);
				SS.push(make_pair(B_third, B.second));
				SS.push(make_pair(B.first, B_third));
			}
			// Step 8, if there are criticals.
			else
			{
				// Bubble sort, which is easy to write. I don't think this will slow down since there are not much critical points.
				/*for (int i = 0; i < SSS.size(); ++i)
					for (int j = 0; j < i; ++j)
					{
						if (SSS[j](0) < SSS[j + 1](0))
							swap(SSS[j], SSS[j + 1]);
						if (SSS[j](0) == SSS[j + 1](0) && SSS[j](1) < SSS[j + 1](1))
							swap(SSS[j], SSS[j + 1]);
						if (SSS[j](0) == SSS[j + 1](0) && SSS[j](1) == SSS[j + 1](1) && SSS[j](2) < SSS[j + 1](2))
							swap(SSS[j], SSS[j + 1]);
						if (SSS[j](0) == SSS[j + 1](0) && SSS[j](1) == SSS[j + 1](1) && SSS[j](2) == SSS[j + 1](2) && SSS[j](3) < SSS[j + 1](3))
							swap(SSS[j], SSS[j + 1]);
					}*/
				SS.push(make_pair(SSS.back(), B.second));
				for (int i = SSS.size() - 2; i >= 0; --i)
					SS.push(make_pair(SSS[i], SSS[i + 1]));
				SS.push(make_pair(B.first, SSS.front()));
			}
		}
		return S;
	}
	vector<pair<Vector4d, Vector4d>> Curve_Trace(Vector4d ini, vector<Vector4d> term, double l,double xi)	// Algorithm 4.
	{
		//cout << term.size() << endl;

		// Step 1, initialize S and S'.
		vector<pair<Vector4d, Vector4d>> S;
		vector<pair<Vector4d, Vector4d>> SS;

		MatrixId B;
		Vector4d P = ini;

		int t = 0;
		bool output = false;
		//output = true;

		do
		{
			++t;
			// Step 2, trace a point.
			Vector4d Q = get_next(P, l);

			// Step 3, get a box sequence.
			if (output)
				cout << "4D decomposing (" << P.transpose() << ")-(" << Q.transpose() << ")" << endl;
			vector<Vector4d> SSS = Decompose4D(P, Q, l);

			// Step 4, find the cover box.
			Vector4d SSS_min = Vector4d(boundary, boundary, boundary, boundary);
			Vector4d SSS_max = Vector4d(-boundary, -boundary, -boundary, -boundary);
			for (int i = 0; i < SSS.size(); ++i)
			{
				if (SSS_min(0) > SSS[i](0))
					SSS_min(0) = SSS[i](0);
				if (SSS_min(1) > SSS[i](1))
					SSS_min(1) = SSS[i](1);
				if (SSS_min(2) > SSS[i](2))
					SSS_min(2) = SSS[i](2);
				if (SSS_min(3) > SSS[i](3))
					SSS_min(3) = SSS[i](3);
				if (SSS_max(0) < SSS[i](0))
					SSS_max(0) = SSS[i](0);
				if (SSS_max(1) < SSS[i](1))
					SSS_max(1) = SSS[i](1);
				if (SSS_max(2) < SSS[i](2))
					SSS_max(2) = SSS[i](2);
				if (SSS_max(3) < SSS[i](3))
					SSS_max(3) = SSS[i](3);
			}
			B = MatrixId(SSS_min, SSS_max);

			// Step 5, insert the boxes into S. If the last box does not contain the terminate points, repeat the steps.
			for (int i = 0; i < SSS.size() - 1; ++i)
				S.push_back(make_pair(SSS[i], SSS[i + 1]));
			if (B(0).min() >= boundary || B(1).min() >= boundary || B(2).min() >= boundary || B(3).min() >= boundary)
				break;
			if (B(0).max() <= -boundary || B(1).max() <= -boundary || B(2).max() <= -boundary || B(3).max() <= -boundary)
				break;
			for (int i = 0; i < term.size(); ++i)
				if (B.contains(term[i]))
					break;
			if(t<10)
				term.push_back(Q);;
			P = Q;
		} while (true);

		// Step 6, get S.

		// Step 7, decompose boxes in S into 3D strongly monotonic curves.
		for (int i = 0; i < S.size(); ++i)
		{
			if (output)
				cout << "3D decomposing (" << S[i].first.transpose() << ")-(" << S[i].second.transpose() << ")" << endl;
			vector<pair<Vector4d, Vector4d>> SSS = Decompose3D(S[i].first, S[i].second);
			for (int j = 0; j < SSS.size(); ++j)
				SS.push_back(SSS[j]);
		}

		// Step 8, decompose S' into smaller intervals.
		vector<pair<Vector4d, Vector4d>> SSSS;

		for (int i = 0; i < SS.size(); ++i)
		{
			double box_size = (Value(SS[i].first) - Value(SS[i].second)).norm();
			double cover_size = (SS[i].first - SS[i].second).norm();
			if (box_size < xi)
				SSSS.push_back(SS[i]);
			else
			{
				int k = 2 * int(box_size / xi);
				Vector4d Q = get_next(SS[i].first, cover_size / k);
				SSSS.push_back(make_pair(SS[i].first, Q));
				while (MatrixId(SS[i].first, SS[i].second).contains(Q))
				{
					Vector4d R = get_next(Q, cover_size / k);
					if (MatrixId(SS[i].first, SS[i].second).contains(R))
						SSSS.push_back(make_pair(Q, R));
					Q = R;
				}
				SSSS.push_back(make_pair(Q, SS[i].second));
			}
		}

		return SSSS;
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

void test0(trial t)
{
	cout << t.the_diff(0) << endl;
	cout << t.the_diff(1) << endl;
	cout << t.the_diff(2) << endl;
	//cout << t.find_root(Vector4d(0, 0, 0, 1), Vector4d(1, 1, 1, 1)).transpose() << endl;
	cout << t.get_next(t.initial, 1).transpose() << endl;
	cout << t.get_last(t.target, 1).transpose() << endl;
}

void test1(trial t)
{
	cout << t.JacobInterval(MatrixId(t.initial, t.target)) << endl;
	cout << t.DeterInterval(MatrixId(t.initial, t.target)) << endl;
	cout << t.Deter3Interval(MatrixId(t.initial, t.target)) << endl;
	cout << t.Jacob((t.initial + t.target) / 2) << endl;
	cout << t.Deter((t.initial + t.target) / 2) << endl;
	cout << t.Deter3((t.initial + t.target) / 2) << endl;
}

void test2(trial t)
{
	cout << t.FaceValueInterval(MatrixId(t.initial, t.target + Vector4d(WE, WE, WE, WE)), t.target, 0).transpose() << endl;
	cout << t.FaceValueInterval(MatrixId(t.initial, t.target + Vector4d(WE, WE, WE, WE)), t.target, 1).transpose() << endl;
	cout << t.FaceValueInterval(MatrixId(t.initial, t.target + Vector4d(WE, WE, WE, WE)), t.target, 2).transpose() << endl;
	cout << t.FaceValueInterval(MatrixId(t.initial, t.target + Vector4d(WE, WE, WE, WE)), t.target, 3).transpose() << endl;
}

void test3(trial t)
{
	vector<int> S = t.IsTerminate(t.initial, t.target);
	for (int i = 0; i < S.size(); ++i)
		cout << S[i] << " ";
	cout << endl;
}

void test4(trial t)
{
	vector<Vector4d> S = t.Decompose4D(t.initial, t.target, 0.5);
	cout << endl << "The boxes are:" << endl;
	for (int i = 0; i < S.size(); ++i)
		cout << "(" << S[i].transpose() << ")" << endl;
}

void test5(trial t)
{
	vector<pair<Vector4d, Vector4d>> S = t.Decompose3D(t.initial, t.target);
	cout << endl << "The boxes are:" << endl;
	for (int i = 0; i < S.size(); ++i)
		cout << "(" << S[i].first.transpose() << ")-(" << S[i].second.transpose() << ")" << endl;
}

void test6(trial t)
{
	vector<Vector4d> terminates;
	terminates.push_back(t.target);
	clock_t start, end;
	double duration;
	start = clock();
	vector<pair<Vector4d, Vector4d>> S = t.Curve_Trace(t.initial, terminates, t.stick_size, 1.0);
	end = clock();
	duration = (end - start) / (double)CLOCKS_PER_SEC;
	cout << endl << "The boxes are: " << endl;
	for (int i = 0; i < S.size(); ++i)
		cout << "(" << S[i].first.transpose() << ")-(" << S[i].second.transpose() << ")" << endl;
	cout << "There are " << S.size() << " boxes." << endl;
	cout << "Cost " << duration << " seconds." << endl;
}

int main(int argc, char* argv[])
{
	const char filename[1000] = "../Input/inputs.js";
	cout << std::setprecision(12) << fixed;
    trial t(read_json(filename));
	cout << "Test 0 (Newton's method):" << endl;
	test0(t);
	cout << endl << "Test 1 (DT vector):" << endl;
	test1(t);
	cout << endl << "Test 2 (Intersection value at a face):" << endl;
	test2(t);
	cout << endl << "Test 3 (Algorithm 1):" << endl;
	test3(t);
	cout << endl << "Test 4 (Algorithm 2):" << endl;
	test4(t);
	cout << endl << "Test 5 (Algorithm 3):" << endl;
	test5(t);
	cout << endl << "Test 6 (Algorithm 4):" << endl;
	test6(t);
	return 0;
}