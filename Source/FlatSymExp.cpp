#include "FlatSymExp.hpp"

#define float float

#define soi sizeof(int)
#define sof sizeof(float)

void pW(unsigned char*& p, int i) { *(int*)p = i; p += soi; }
void pW(unsigned char*& p, float f) { *(float*)p = f; p += sof; }
void pW(unsigned char*& p, size_t t) { pW(p, (int)t); }
int& iR(unsigned char*& p) { p += soi; return *(int*)(p-soi); }
float& fR(unsigned char*& p) { p += sof; return *(float*)(p-sof); }

FlatSymExp::FlatSymExp(const SymExp& exp)
{
	unsigned char* bufIndex;

	Reserve(1024);
	std::vector<int> ids = exp.GetIds();
	bufIndex = Extend(soi+(soi*ids.size())); //ids count, ids
	pW(bufIndex, ids.size());
	for (int i = 0; i < ids.size(); i++)
		pW(bufIndex, ids[i]);

	bufIndex = Extend(2*soi); //exp count, exp indices
	pW(bufIndex, 1);

	expIndex(0) = size;
	bufIndex = Extend(sof+soi); //scalar, number of terms
	pW(bufIndex, exp.scalar);
	pW(bufIndex, exp.terms.size());
	for (int i = 0; i < exp.terms.size(); i++)
	{
		bufIndex = Extend(sof+soi+(2*soi*exp.terms[i].ids.size())); //coeff, number of factors, ids, pows
		pW(bufIndex, exp.terms[i].coeff);
		pW(bufIndex, exp.terms[i].ids.size());
		for (int j = 0; j < exp.terms[i].ids.size(); j++)
		{
			pW(bufIndex, BinSearch(ids, exp.terms[i].ids[j]));
			pW(bufIndex, exp.terms[i].pows[j]);
		}
	}

}

FlatSymExp::FlatSymExp(const std::vector<SymExp>& exps)
{
	Reserve(1024);


	unsigned char* bufIndex;

	Reserve(1024);
	std::vector<int> ids = GetIds(exps);
	bufIndex = Extend(soi+(soi*ids.size())); //ids count, ids
	pW(bufIndex, ids.size());
	for (int i = 0; i < ids.size(); i++)
		pW(bufIndex, ids[i]);

	bufIndex = Extend(soi*(exps.size()+1)); //exp count, exp indices
	pW(bufIndex, exps.size());

	for (int e = 0; e < exps.size(); e++)
	{
		const SymExp& exp = exps[e];

		expIndex(e) = size;
		bufIndex = Extend(sof+soi); //scalar, number of terms
		pW(bufIndex, exp.scalar);
		pW(bufIndex, exp.terms.size());
		for (int i = 0; i < exp.terms.size(); i++)
		{
			bufIndex = Extend(sof+soi+(2*soi*exp.terms[i].ids.size())); //coeff, number of factors, ids, pows
			pW(bufIndex, exp.terms[i].coeff);
			pW(bufIndex, exp.terms[i].ids.size());
			for (int j = 0; j < exp.terms[i].ids.size(); j++)
			{
				pW(bufIndex, BinSearch(ids, exp.terms[i].ids[j]));
				pW(bufIndex, exp.terms[i].pows[j]);
			}
		}
	}
}



void FlatSymExp::PushBack(int i)
{
	Resize(size+sizeof(int));
	*(int*)(data+size) = i;
}

void FlatSymExp::PushBack(float f)
{
	Resize(size+sof);
	*(float*)(data+size) = f;
}

int& FlatSymExp::iAt(int i)
{
	return *(int*)(data+i);
}

float& FlatSymExp::fAt(int i)
{
	return *(float*)(data+i);
}

const int& FlatSymExp::iAtC(int i) const
{
	return *(int*)(data+i);
}

const float& FlatSymExp::fAtC(int i) const
{
	return *(float*)(data+i);
}

void FlatSymExp::Reserve(int i)
{
	if (i <= capacity)
		return;
	if (data == nullptr)
	{
		data = new unsigned char[i];
		capacity = i;
		return;
	}
	unsigned char* temp = new unsigned char[i];
	for (int j = 0; j < capacity; j++)
		temp[j] = data[j];
	capacity = i;
	delete[] data;
	data = temp;
}

void FlatSymExp::Resize(int i)
{
	if (i > capacity)
	 	Reserve(i*2);
	size = i;
}

unsigned char* FlatSymExp::Extend(int i)
{
	int s = size;
	Resize(size+i);
	return data+s;
}

int FlatSymExp::idCount() const
{
	return iAtC(0);
}

int FlatSymExp::idKey(int no) const
{
	return iAtC(soi*(1+no));
}

int FlatSymExp::expCount() const
{
	return iAtC(soi*(1+iAtC(0)));
}

int& FlatSymExp::expIndex(int exp)
{
	int index = 0;
	index += soi*(1+iAt(index)) + soi;
	return iAt(index + (soi*exp));
}

const int& FlatSymExp::expIndexC(int exp) const
{
	int index = 0;
	index += soi*(1+iAtC(index)) + soi;
	return iAtC(index + (soi*exp));
}

std::vector<float> FlatSymExp::SclEval(const std::vector<int>& valIds, const std::vector<float>& values) const
{
	std::vector<float> hold;
	std::vector<float> coords;


	int expC = expCount();
	hold.resize(expC);
	
	coords.resize(idCount());
	for (int i = 0; i < coords.size(); i++)
	{
		int index = BinSearch(valIds, idKey(i));
		if (index == valIds.size() || valIds[index] != idKey(i))
			return hold;
		else
			coords[i] = values[index];
	}

	unsigned char* bufIndex = data+expIndexC(0);

	for (int i = 0; i < expC; i++)
	{
		hold[i] = fR(bufIndex);
		int prodCount = iR(bufIndex);

		for (int j = 0; j < prodCount; j++)
		{
			float prodVal = fR(bufIndex);
			int facCount = iR(bufIndex);

			for (int k = 0; k < facCount; k++)
			{
				float facVal = coords[iR(bufIndex)];
				int pow = iR(bufIndex);
				for (int p = 0; p < pow; p++)
					prodVal *= facVal;
			}

			hold[i] += prodVal;
		}
	}

	return hold;
}

#undef float

