// expression.cpp : This file defines methods for parsing an expression.
// Author: zz1918@nyu.edu

#include<cassert>
#include<string>
using namespace std;

// If c is a number
bool is_num(char c)
{
	return ('0' <= c) && (c <= '9');
}
// If c is an standard operator (+,-,*,/,^)
bool is_opt(char c)
{
	return (c == '+') || (c == '-') || (c == '*') || (c == '/') || (c == '^');
}
// Remove ' ' from string.
string remove_blank(string expression)
{
	string new_expression = "";
	for (int i = 0; i < expression.length(); ++i)
		if (expression[i] != ' ')
			new_expression += expression[i];
	return new_expression;
}
// Catch the number from place i, pointer i will point to end number place j after returning
unsigned int catch_int(string expression, int& i)
{
	assert(is_num(expression[i]));
	int j = i;
	while (j < expression.length() && is_num(expression[j]))
		++j;
	int begin = i, length = j - i;
	i = j - 1;
	return (unsigned int)(stoi(expression.substr(begin, length)));
}
// Catch the number from place i, pointer i will point to end number place j after returning
double catch_num(string expression, int& i)
{
	assert(is_num(expression[i]) || expression[i] == '.');
	int j = i;
	while (j < expression.length() && (is_num(expression[j]) || expression[j] == '.' || expression[j] == 'e'))
		++j;
	int begin = i, length = j - i;
	i = j - 1;
	return stod(expression.substr(begin, length));
}
// Match the substring for the left enclosure at place i, pointer i will point to right enclosure place j after returning
string match(string expression, int& i)
{
	assert(expression[i] == '(');
	int j = i + 1;
	int deepness = 0;
	for (; j < expression.length(); ++j)
	{
		if (expression[j] == '(')
			deepness += 1;
		if (expression[j] == ')')
			deepness -= 1;
		if (deepness < 0)
			break;
	}
	int begin = i + 1, length = j - i - 1;
	i = j;
	return expression.substr(begin, length);
}
// Add "*" to string to make it readable for machines, list for the list of vars.
string add_mul(string expression, char* var_list = NULL, int list_size = 0)
{
	string new_expression = expression.substr(0, 1);
	for (int i = 1; i < expression.length(); ++i)
	{
		char last = expression[i - 1];
		char current = expression[i];
		if (is_opt(last) || is_num(current) || current == '.')
		{
			new_expression += current;
			continue;
		}
		bool unchanged = true;
		for (int j = 0; j < list_size; ++j)
			if (current == *(var_list + j))
			{
				new_expression += '*';
				new_expression += current;
				unchanged = false;
				break;
			}
		if (unchanged)
			new_expression += current;
	}
	return new_expression;
}