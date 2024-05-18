#pragma once

#include <glm/glm.hpp>

#include <string>
#include <iostream>
#include <vector>

#include "Vector.hpp"

#define float float

/*
Putting this here for later:

How to easily use const with these classes

#define uc(o) const_cast<decltype(o)&>(o)

*/

class SymExp;
class SymParser;
class SymExpTable;

class Product
{
public:
	float coeff = 1;
	std::vector<int> ids;
	std::vector<int> pows; //negative pows not allowed

	Product() {}
	explicit Product(float coeffFloat) : coeff(coeffFloat) {}
	Product(float coeffFloat, int id) : coeff(coeffFloat) { MultId(id); }
	Product(float coeffFloat, int id, int pow) : coeff(coeffFloat) { MultId(id, pow); }
	Product(std::vector<int> idsVec) : ids(idsVec) {}
	Product(float coeffFloat, std::vector<int> idsVec) : coeff(coeffFloat), ids(idsVec) {}
	Product(float coeffFloat, std::vector<int> idsVec, std::vector<int> powsVec) : coeff(coeffFloat), ids(idsVec), pows(powsVec) {}

	void ClearEmpty();
	Product& MultId(int id, int pow = 1);

	Product(const Product& prod) { coeff = prod.coeff; ids = prod.ids; pows = prod.pows; }
	Product& operator=(Product& prod) { coeff = prod.coeff; ids = prod.ids; pows = prod.pows; return *this; }
	Product(Product&& prodMove) { coeff = prodMove.coeff; ids = prodMove.ids; pows = prodMove.pows; }
	Product& operator=(Product&& prodMove) { coeff = prodMove.coeff; ids = prodMove.ids; pows = prodMove.pows; return *this; }
	
	SymExp operator+(Product& prod);
	SymExp operator+(Product&& prod);
	SymExp operator-(Product& prod);
	SymExp operator-(Product&& prod);
	Product operator*(Product& prod);
	SymExp operator*(Product&& prod);
	Product operator*(float co) { return Product(coeff*co, ids, pows); };
	//Product operator/(Product& prod); //negative powers are deprecated
	Product operator/(float co) { return Product(coeff/co, ids, pows); };
	Product& operator*=(Product& prod) { operator=(operator*(prod)); return *this; }
	Product& operator*=(Product&& prod) { operator=(operator*(prod)); return *this; }
	//Product& operator/=(Product& prod) { operator=(operator/(prod)); return *this; }

	bool CheckType(Product& prod);

	SymExp Eval(SymExpTable& table);
	float SclEval(SymExpTable& table);
	SymExp Eval(std::vector<int> valIds, std::vector<float> values);
	float SclEval(std::vector<int> valIds, std::vector<float> values);

	std::string ToString();
	std::string ToString(SymParser& parser, bool genDef = true);

	std::vector<int> GetIds();
	std::vector<Product> Gradient();
	void Gradient(std::vector<Product>& out);
	std::vector<Product> Gradient(std::vector<int>& ids);
	void Gradient(std::vector<int>& ids, std::vector<Product>& grad);
};

class SymExp
{
public:
	float scalar = 0;
	std::vector<Product> terms;

	SymExp() {}
	explicit SymExp(float scalarFloat) : scalar(scalarFloat) {}
	SymExp(std::vector<Product> termsVec) : terms(termsVec) {}
	SymExp(float scalarFloat, std::vector<Product> termsVec) : scalar(scalarFloat), terms(termsVec) {}
	SymExp(Product prod) { terms.push_back(prod); }

	SymExp(const SymExp& symExp) { scalar = symExp.scalar; terms = symExp.terms; }
	SymExp& operator=(SymExp& symExp) { scalar = symExp.scalar; terms = symExp.terms; return *this; }
	SymExp(SymExp&& symExpMove) { scalar = symExpMove.scalar; terms = symExpMove.terms;}
	SymExp& operator=(SymExp&& symExpMove) { scalar = symExpMove.scalar; terms = symExpMove.terms; return *this; }
	SymExp operator+(SymExp& symExp);
	SymExp operator+(SymExp&& symExp) { return operator+(symExp); };
	SymExp operator-(SymExp& symExp);
	SymExp operator-(SymExp&& symExp) { return operator-(symExp); };
	SymExp operator*(SymExp& symExp);
	SymExp operator*(SymExp&& symExp) { return operator*(symExp); };
	//SymExp operator/(Product& prod);
	//SymExp operator/(Product&& prodMove) { Product hold = prodMove; return operator/(hold); };
	SymExp& operator+=(SymExp& symExp) { operator=(operator+(symExp)); return *this; }
	SymExp& operator+=(SymExp&& symExp) { operator=(operator+(symExp)); return *this; }
	SymExp& operator-=(SymExp& symExp) { operator=(operator-(symExp)); return *this; }
	SymExp& operator-=(SymExp&& symExp) { operator=(operator-(symExp)); return *this; }
	SymExp& operator*=(SymExp& symExp) { operator=(operator*(symExp)); return *this; }
	SymExp& operator*=(SymExp&& symExp) { operator=(operator*(symExp)); return *this; }
	//SymExp& operator/=(Product& prod) { operator=(operator/(prod)); return *this; }

	void Simplify();

	std::string ToString();
	std::string ToString(SymParser& parser, bool genDef = true);
	SymExp Eval(SymExpTable& table);
	float SclEval(SymExpTable& table);
	float SclEval(SymExpTable&& table) { return SclEval(table); }
	SymExp Eval(std::vector<int> valIds, std::vector<float> values);
	float SclEval(std::vector<int> valIds, std::vector<float> values);

	std::vector<int> GetIds();
	std::vector<SymExp> Gradient();
	void Gradient(std::vector<SymExp>& grad);
	std::vector<SymExp> Gradient(std::vector<int>& gradIds);
	void Gradient(std::vector<int>& gradIds, std::vector<SymExp>& grad);
};

std::vector<int> GetIds(std::vector<SymExp>& exps);
Vector2D<SymExp> Gradient(std::vector<SymExp>& exps);
void Gradient(std::vector<SymExp>& exps, Vector2D<SymExp>& grad);
Vector2D<SymExp> Gradient(std::vector<SymExp>& exps, std::vector<int>& gradIds);
void Gradient(std::vector<SymExp>& exps, std::vector<int>& gradIds, Vector2D<SymExp>& grad);
std::vector<float> SclEvalVec(std::vector<SymExp>& exps, SymExpTable& table);
std::vector<float> SclEvalVec(std::vector<SymExp>& exps, std::vector<int>& valIds, std::vector<float>& values);
Vector2D<float> SclEvalVec2D(Vector2D<SymExp>& exps, SymExpTable& table);
Vector2D<float> SclEvalVec2D(Vector2D<SymExp>& exps, std::vector<int>& valIds, std::vector<float>& values);
std::vector<float> NMnTo1(SymExp& poly, std::vector<int> valIds, std::vector<float> initial, float threshold = 0.00001);
std::vector<float> NMnTom(std::vector<SymExp>& polys, std::vector<int> valIds, std::vector<float> initial, float threshold = 0.00001);

//takes a vec2D that represents a map that takes v1 to v2, with a column describing how a coordinate of v1 maps to a vector of v2 (like a matrix), then returns a set of linear equations which describe how to replace the og coords xi\
with new coords x'i in the form of < <x1 = x'1 + m2x2 + m3x3 + ...>, <x2 = x'2 + n3x3 + ...>, ...> . Given xi, one can calculate x'i by moving xi terms to one side. Given x'i, one can find xi by starting at xn = x'n, and working down until n = 1.\
Also writes into outNewMap the projection of the vecMap on the new basis (outNewMap has orthogonal column vectors). outNewMap can be thought of as a new map on the x'i basis, i.e. vecMap*x = outNewMap*x'.\
Note that return value shouldn't be treated like a matrix nor a vector of vectors, it represents the right side of a system of linear equations, so conversion between basis must be done following the rules mentioned above.
//Vector2D<float> CoeffFromGramSchmidt(const Vector2D<float>& vecMap, Vector2D<float>& outNewMap);
//Vector2D<float> CoeffFromGramSchmidt(const Vector2D<float>& vecMap, Vector2D<float>* outNewMap);

//resets thread_local variables used for NM and CoeffFromGramSchmidt
void ResetGlobals();

/* std::vector<SymExp> ComplexPoly(const SymExp exp); //Calculates the sytem that represents generalizing a polynomial on reals to a polynomial on complex numbers
std::vector<SymExp> ComplexPoly(const SymExp exp, SymParser& parse);
std::vector<SymExp> ComplexPoly(const SymExp exp, SymParser& parse, const std::vector<int> ids, const std::vector<int> cIds);
std::vector<SymExp> ComplexPoly(const SymExp exp, const std::vector<int> ids, const std::vector<int> cIds); */

class SymExpTable
{
public:
	std::vector<int> lookup;
	std::vector<SymExp> exps;

	void Add(int key, float val) { SymExp hold = SymExp(val); Add(key, hold); }
	void Add(int key, SymExp exp); //check if key is already in lookup, update if it is, otherwise push_back
};

std::string LettersFromNum(int num); // 0 -> a, 1 -> b, 26 -> aa

class SymParser
{
public:
	std::vector<int> lookup;
	std::vector<std::string> names;

	void Add(int key, std::string name); //check if key is already in lookup, update if it is
	void GenerateDefaults(Product& prod);
	void GenerateDefaults(SymExp& symExp);
	std::string SearchId(int id);
	int SearchName(std::string& name);
};

std::string FloatToString(float& floatRef);

int BinSearch(std::vector<int>& vec, int val);
int BinSearch(std::vector<int>& vec, int lower, int upper, int val);
#undef float