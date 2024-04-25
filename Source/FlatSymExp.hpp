#pragma once
#include "SymExp.hpp"

#define float float

#include <vector>

//contains multiple exps btw
class FlatSymExp
{
public:
	//initial destructible state
	char* data = nullptr; //raw data of SymExp
	int capacity = 0;

	//account for sizeof(float)

	//getters
	int idCount();
	int idKey(int no);
	int expCount();
	int expIndex(int exp);
	float expScalar(int exp);
	int prodCount(int exp);
	int prodIndex(int exp, int prod);
	float prodCoeff(int exp, int prod);
	int factorCount(int exp, int prod);
	int factorIndex(int exp, int prod, int factor);
	int factorId(int exp, int prod, int factor);
	int factorPow(int exp, int prod, int factor);

	//Only implement setters for the stuff needed right now

	FlatSymExp() = default;
	FlatSymExp(SymExp se);

	void GenerateSymExp(SymExp& out);

	void CopyData(char* newAddress);
	~FlatSymExp() { if (data != nullptr) delete[] data;}
};

FlatSymExp CreateFlat(const SymExp& exp);
FlatSymExp CreateFlat(const std::vector<SymExp>& exps);


#undef float