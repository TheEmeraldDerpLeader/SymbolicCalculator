#pragma once
#include "SymExp.hpp"

#define float float

#include <vector>

//plan of action:
//1. Evaluate flat (done)
//2. Gradient of flat (done)
//3. NM on flat

//contains multiple exps btw
class FlatSymExp
{
public:
	//initial destructible state
	unsigned char* data = nullptr; //raw data of SymExp
	/* Structure of data:
	(int) number of ids - i
	i*(int) association from global ids to local ids (local ids are index of global in this array)
	(int) number of symExp - n
	n*(int) symExp indices in data
	for each n:
		(float) scalar of symExp
		(int) product count - m
		for each m:
			(float) coefficient of product
			(int) product size - q
			q*(int,int) id of product, power of product
	
	*/
	int size = 0;
	int capacity = 0;
	void PushBack(size_t t) { PushBack((int)t); }
	void PushBack(int i);
	void PushBack(float f);
	int& iAt(int i);
	float& fAt(int i);
	const int& iAtC(int i) const;
	const float& fAtC(int i) const;
	void Reserve(int i);
	void Resize(int i); //calls reserve(2*i) if i is > capacity
	unsigned char* Extend(int i);

	//account for sizeof(float)

	//getters
	int idCount() const;
	int idKey(int no) const;
	int expCount() const;
	int& expIndex(int exp);
	const int& expIndexC(int exp) const;
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
	FlatSymExp(const SymExp& exp);
	FlatSymExp(const std::vector<SymExp>& exps);

	FlatSymExp(FlatSymExp&& exp);
	FlatSymExp& operator=(FlatSymExp&& exp);

	void GenerateSymExp(SymExp& out);
	void GenerateSymExps(std::vector<SymExp>& out);

	std::vector<float> SclEval(const std::vector<int>& valIds, const std::vector<float>& values) const;

	FlatSymExp Gradient() const;

	void CopyData(char* newAddress);
	~FlatSymExp() { if (data != nullptr) delete[] data; size = 0; capacity = 0; }
};

FlatSymExp CreateFlat(const SymExp& exp);
FlatSymExp CreateFlat(const std::vector<SymExp>& exps);


#undef float