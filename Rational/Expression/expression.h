// expression.h
// Author: zz1918@nyu.edu

#include "expression.cpp"

bool is_num(char c);
bool is_opt(char c);
string remove_blank(string expression);
unsigned int catch_int(string expression, int& i);
double catch_num(string expression, int& i);
string match(string expression, int& i);
string add_mul(string expression, char* var_list, int list_size);