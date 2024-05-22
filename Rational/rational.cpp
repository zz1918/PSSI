// rational.cpp : This file defines rational functions and
// gives methods of converting from an expression string
// and valuing and differentiating.
// Author: zz1918@nyu.edu

#include<iostream>
#include<vector>
#include<cassert>
#include<string>
#include<expression.h>
using namespace std;
#define opts_num 5
const char Rational_Opts[opts_num] = { 'n','+','-','*','/' };
class Rational
{
	vector<char>* alphabets;				// Storing the alphabets for variables.
	int opt = 0;							// Operator
	double num = 0.0;						// Number (if applicable)
	unsigned int pw = 1;					// Power (if applicable)
	bool root = true;						// Is the root of an expression?
	Rational* next[2] = { NULL,NULL };
	int find_opt(char o)
	{
		for (int i = 0; i < opts_num; ++i)
			if (o == Rational_Opts[i])
				return i;
		for (int i = 0; i < size(); ++i)
			if (o == (*alphabets)[i])
				return i + opts_num;
		assert(false && "Invalid function operator!");
		return -1;
	}
public:
	Rational(vector<char>* a = new vector<char>)
	{
		alphabets = a;
	}
	Rational(vector<char>* a, char o)
	{
		alphabets = a;
		opt = find_opt(o);
	}
	Rational(vector<char>* a, double n)
	{
		alphabets = a;
		num = n;
	}
	Rational(Rational* r)
	{
		alphabets = r->vars();
		opt = find_opt(r->Opt());
		num = r->Num();
		pw = r->PW();
		root = r->is_root();
		next[0] = r->left();
		next[1] = r->right();
	}
	~Rational()
	{
		delete next[0];
		delete next[1];
	}
	int size()
	{
		return alphabets->size();
	}
	bool is_primitive()
	{
		return opt == 1 || opt == 2;
	}
	bool is_root()
	{
		return root;
	}
	bool is_num()
	{
		return opt == 0;
	}
	bool is_var()
	{
		return opt >= opts_num;
	}
	vector<char>* vars()
	{
		return alphabets;
	}
	void add_var(char o)
	{
		alphabets->push_back(o);
	}
	void set_var(vector<char>* a)
	{
		alphabets = a;
	}
	void set_opt(char o)
	{
		opt = find_opt(o);
	}
	void set_num(double n)
	{
		num = n;
	}
	void set_pw(unsigned int p)
	{
		pw = p;
	}
	void set_pw(int p)
	{
		assert(p >= 0);
		set_pw((unsigned int)(p));
	}
	void set_left(Rational* r)
	{
		next[0] = r;
	}
	void set_right(Rational* r)
	{
		next[1] = r;
	}
	void set_root(bool r)
	{
		root = r;
	}
	void set_copy(Rational* r)
	{
		alphabets = r->vars();
		opt = find_opt(r->Opt());
		num = r->Num();
		pw = r->PW();
		root = r->is_root();
		next[0] = r->left();
		next[1] = r->right();
	}
	void clear()
	{
		opt = 0;
		num = 0.0;
		pw = 1;
		next[0] = next[1] = NULL;
		root = true;
	}
	char Opt()
	{
		if (opt < opts_num)
			return Rational_Opts[opt];
		else
			return (*alphabets)[opt - opts_num];
	}
	double Num()
	{
		return num;
	}
	unsigned int PW()
	{
		return pw;
	}
	Rational* left()
	{
		return next[0];
	}
	Rational* right()
	{
		return next[1];
	}
	Rational* right_most()
	{
		if (next[1] == NULL)
			return this;
		else
			return right()->right_most();
	}
	template<typename Scarlar>
	Scarlar value(vector<Scarlar> Vars)
	{
		if (opt < 0)
			return Scarlar(0);
		if (pw == 0)
			return Scarlar(1);
		Scarlar raw;
		switch (opt)
		{
		case 0:raw = Scarlar(num); break;
		case 1:raw = next[0]->value(Vars) + next[1]->value(Vars); break;
		case 2:raw = next[0]->value(Vars) - next[1]->value(Vars); break;
		case 3:raw = next[0]->value(Vars) * next[1]->value(Vars); break;
		case 4:raw = next[0]->value(Vars) / next[1]->value(Vars); break;
		default:raw = Vars[opt - opts_num];
		}
		if (pw > 1)
			return pow(raw, pw);
		else
			return raw;
	}
	template<typename Scarlar>
	Scarlar deriv(vector<Scarlar> Vars, int Opt)
	{
		if (opt < 0)
			return Scarlar(0);
		if (pw == 0)
			return Scarlar(0);
		Scarlar raw;
		switch (opt)
		{
		case 0:raw = Scarlar(0); break;
		case 1:raw = next[0]->deriv(Vars, Opt) + next[1]->deriv(Vars, Opt); break;
		case 2:raw = next[0]->deriv(Vars, Opt) - next[1]->deriv(Vars, Opt); break;
		case 3:raw = next[0]->deriv(Vars, Opt) * next[1]->value(Vars) + next[0]->value(Vars) * next[1]->deriv(Vars, Opt); break;
		case 4:raw = (next[0]->deriv(Vars, Opt) * next[1]->value(Vars) - next[0]->value(Vars) * next[1]->deriv(Vars, Opt)) / pow(next[1]->value(Vars), 2); break;
		default:raw = ((opt == Opt + opts_num) ? Scarlar(1) : Scarlar(0));
		}
		if (pw > 1)
		{
			int PW = pw;
			set_pw(PW - 1);
			Scarlar pure = Scarlar(PW) * raw * value<Scarlar>(Vars);
			set_pw(PW);
			return pure;
		}
		else
			return raw;
	}
	template<typename Scarlar>
	Scarlar deriv(vector<Scarlar> Vars, char Var)
	{
		int Opt = find_opt(Var) - opts_num;
		return deriv(Vars, Opt);
	}
	Rational* deriv(int Opt)
	{
		Rational* result = new Rational(alphabets);
		if (opt < 0)
			return result;
		if (pw == 0)
			return result;
		switch (opt)
		{
		case 0: break;
		case 1:
		{
			result->set_opt('+');
			result->set_left(next[0]->deriv(Opt));
			result->set_right(next[1]->deriv(Opt));
		}
			break;
		case 2:
		{
			result->set_opt('-');
			result->set_left(next[0]->deriv(Opt));
			result->set_right(next[1]->deriv(Opt));
		}
			break;
		case 3:
		{
			Rational* next_left = new Rational(alphabets);
			Rational* next_right = new Rational(alphabets);
			next_left->set_opt('*');
			next_left->set_left(next[0]->deriv(Opt));
			next_left->set_right(next[1]);
			next_right->set_opt('*');
			next_right->set_left(next[1]->deriv(Opt));
			next_right->set_right(next[0]);
			result->set_opt('+');
			result->set_left(next_left);
			result->set_right(next_right);
		}
			break;
		case 4:
		{
			Rational* nominator = new Rational(alphabets);
			Rational* denominator = new Rational(alphabets);
			Rational* nominator_left = new Rational(alphabets);
			Rational* nominator_right = new Rational(alphabets);
			nominator_left->set_opt('*');
			nominator_left->set_left(next[0]->deriv(Opt));
			nominator_left->set_right(next[1]);
			nominator_right->set_opt('*');
			nominator_right->set_left(next[0]);
			nominator_right->set_right(next[1]->deriv(Opt));
			nominator->set_opt('-');
			nominator->set_left(nominator_left);
			nominator->set_right(nominator_right);
			denominator->set_copy(next[1]);
			denominator->set_pw(2 * denominator->PW());
			result->set_opt('/');
			result->set_left(nominator);
			result->set_right(denominator);
		}
			break;
		default:if (opt == Opt + opts_num)
					result->set_num(1.0);
		}
		if (pw > 1)
		{
			Rational* power_part = new Rational(alphabets);
			Rational* coeff = new Rational(alphabets, double(pw));
			Rational* deriv_part = new Rational(alphabets);
			Rational* next_left = new Rational(alphabets);
			power_part->set_copy(this);
			power_part->set_pw(pw - 1);
			deriv_part->set_copy(result);
			next_left->set_opt('*');
			next_left->set_left(coeff);
			next_left->set_right(deriv_part);
			result->clear();
			result->set_opt('*');
			result->set_left(next_left);
			result->set_right(power_part);
		}
		return result;
	}
	Rational* deriv(char Var)
	{
		int Opt = find_opt(Var) - opts_num;
		return deriv(Opt);
	}
	void view(ostream& os)
	{
		if (opt < 0)
		{
			os << "Invalid function!" << endl;
			return;
		}
		switch (opt)
		{
		case 0:os << num; break;
		case 1:os << "("; next[0]->view(os); os << "+"; next[1]->view(os); os << ")"; break;
		case 2:os << "("; next[0]->view(os); os << "-"; next[1]->view(os); os << ")"; break;
		case 3:os << "("; next[0]->view(os); os << "*"; next[1]->view(os); os << ")"; break;
		case 4:os << "("; next[0]->view(os); os << "/"; next[1]->view(os); os << ")"; break;
		default:os << (*alphabets)[opt - opts_num];
		}
		if (pw > 1)
			os << "^" << pw;
	}
};
#undef opts_num

ostream& operator<<(ostream& os, Rational* r)
{
	r->view(os);
	return os;
}

Rational* make_rational(string expression, vector<char>* a)
{
	Rational* result = new Rational(a);
	for (int i = 0; i < expression.length(); ++i)
	{
		char o = expression[i];
		switch (o)
		{
		case '+':
		case '-':
		{
			Rational* new_result = new Rational(a);
			new_result->set_left(result);
			new_result->set_opt(o);
			new_result->set_root(false);
			result = new_result;
		}
		break;
		case '*':
		case '/':
		{
			if (!result->is_primitive() || result->is_root())
			{
				Rational* new_result = new Rational(a);
				new_result->set_left(result);
				new_result->set_opt(o);
				new_result->set_root(false);
				result = new_result;
			}
			else
			{
				Rational* new_right = new Rational(a, o);
				new_right->set_left(result->right());
				new_right->set_root(false);
				result->set_right(new_right);
			}
		}
		break;
		case '(':
		{
			if (result->is_num() || result->is_var())						// cover the current scalar or atom
				result = make_rational(match(expression, i), a);
			else															// there is a non-scalar and non-atom expression
			{
				Rational* inner = make_rational(match(expression, i), a);
				result->right_most()->set_right(inner);
			}
		}
		break;
		case '^':
		{
			++i;
			if (result->is_root())
				result->set_pw(catch_int(expression, i));
			else
				result->right_most()->set_pw(catch_int(expression, i));
		}
		break;
		default:
		{
			if (is_num(o) || o == '.')
			{
				if (result->is_num() || result->is_var())					// cover the current scalar or atom
					result->set_num(catch_num(expression, i));
				else                    // there is a non-scalar and non-atom expression
					result->right_most()->set_right(new Rational(a, catch_num(expression, i)));
			}
			else
				if (result->is_num() || result->is_var())					// cover the current scalar or atom
					result->set_opt(o);
				else
					result->right_most()->set_right(new Rational(a, o));
			result->set_root(false);
		}
		}
	}
	result->set_root(true);
	return result;
}

Rational* make_rational(string expression, string s)
{
	vector<char>* a = new vector<char>;
	for (int i = 0; i < s.length(); ++i)
		a->push_back(s[i]);
	return make_rational(expression, a);
}

vector<char>* make_alphabets(string s)
{
	vector<char>* a = new vector<char>;
	for (int i = 0; i < s.length(); ++i)
		a->push_back(s[i]);
	return a;
}

Rational* sum(Rational* r1, Rational* r2)
{
	assert(r1->vars() == r2->vars());
	Rational* result = new Rational(r1->vars());
	result->set_left(r1);
	result->set_right(r2);
	result->set_opt('+');
	return result;
}
Rational* diff(Rational* r1, Rational* r2)
{
	assert(r1->vars() == r2->vars());
	Rational* result = new Rational(r1->vars());
	result->set_left(r1);
	result->set_right(r2);
	result->set_opt('-');
	return result;
}
Rational* prod(Rational* r1, Rational* r2)
{
	assert(r1->vars() == r2->vars());
	Rational* result = new Rational(r1->vars());
	result->set_left(r1);
	result->set_right(r2);
	result->set_opt('*');
	return result;
}
Rational* quot(Rational* r1, Rational* r2)
{
	assert(r1->vars() == r2->vars());
	Rational* result = new Rational(r1->vars());
	result->set_left(r1);
	result->set_right(r2);
	result->set_opt('/');
	return result;
}
Rational* sum(Rational* r1, Rational* r2, Rational* r3)
{
	return sum(sum(r1, r2), r3);
}
Rational* prod(Rational* r1, Rational* r2, Rational* r3)
{
	return prod(prod(r1, r2), r3);
}