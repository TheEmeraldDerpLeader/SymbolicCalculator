#include "FlatSymExp.hpp"

#define float float

#define soi sizeof(int)
#define sof sizeof(float)

template<typename T>
inline T&& mov(T& t) { return reinterpret_cast<T&&>(t); }

void pW(unsigned char*& p, int i) { *(int*)p = i; p += soi; }
void pW(unsigned char*& p, float f) { *(float*)p = f; p += sof; }
void pW(unsigned char*& p, size_t t) { pW(p, (int)t); }
int& iR(unsigned char*& p) { p += soi; return *(int*)(p-soi); }
float& fR(unsigned char*& p) { p += sof; return *(float*)(p-sof); }
const int& iR(const unsigned char*& p) { p += soi; return *(int*)(p-soi); }
const float& fR(const unsigned char*& p) { p += sof; return *(float*)(p-sof); }
const int& iR(const unsigned char* const & p) { return *(int*)(p); } //very cool, 4 function defs bc we want consts

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

FlatSymExp::FlatSymExp(FlatSymExp&& exp)
{
	data = exp.data;
	size = exp.size;
	capacity = exp.capacity;
	exp.data = nullptr;
	exp.size = 0;
	exp.capacity = 0;
}

FlatSymExp& FlatSymExp::operator=(FlatSymExp&& exp)
{
	if (data != nullptr)
		delete[] data;
	data = exp.data;
	size = exp.size;
	capacity = exp.capacity;
	exp.data = nullptr;
	exp.size = 0;
	exp.capacity = 0;
	return *this;
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

FlatSymExp FlatSymExp::Gradient() const
{
	int idC = idCount();
	int expC = expCount();
	FlatSymExp f;
	f.Reserve(expIndexC(0) + idC*( size-expIndexC(0) ));

	const unsigned char* thisIndex = data;
	unsigned char* bufIndex;

	//Write ids
	bufIndex = f.Extend(soi*(1+idC)+soi*(1+(expC*idC)));
	for (int i = 0; i < soi*(1+idC); i++)
	{
		*bufIndex = *thisIndex;
		bufIndex++; thisIndex++;
	}
	
	//Write new exp count
	pW(bufIndex, expC*idC);
	bufIndex += soi*(expC*idC);
	thisIndex += soi*(1+expC);

	for (int j = 0; j < idC; j++)
	{ 
		for (int i = 0; i < expC; i++)
		{
			//set thisIndex to start of SymExp
			thisIndex = data+expIndexC(i)+sof;
			int prodC = iR(thisIndex);


			f.expIndex((i*idC)+j) = f.size;

			int coeffIndex = f.size;
			int prodCIndex = f.size+sof;
			f.Extend(sof+soi);
			pW(bufIndex, float());
			pW(bufIndex,prodC);


			for (int p = 0; p < prodC; p++)
			{
				float prodF = fR(thisIndex);
				int facC = iR(thisIndex);

				//search factors for id == j
				//if there and pow > 1, reduce pow by 1
				//if pow == 1 and facC > 1, remove factor and keep rest. Otherwise add factor to coeffIndex
				//if not there, remove product and decrease prodCIndex by one.

				//search factors
				int jIdPow = 0;
				for (int f = 0; f < facC; f++)
					if (iR(thisIndex+f*(2*soi)) == j)
					{
						jIdPow = iR(thisIndex+soi+f*(2*soi));
						break;
					}
				

				//don't add product if id isn't a factor
				if (jIdPow == 0)
				{
					thisIndex += soi*2*facC;
					(*(int*)(f.data+prodCIndex))--;
					continue;
				}

				//Extend f appropriately
				if (jIdPow == 1)
					f.Extend(sof+soi + (facC-1)*(2*soi));
				else
					f.Extend(sof+soi + facC*(2*soi));
				//Write coeff times pow (bc derivative)
				pW(bufIndex, prodF*jIdPow);

				//Write new product
				if (jIdPow == 1) //remove factor of id == j
				{
					pW(bufIndex, facC-1);
					for (int f = 0; f < facC; f++)
					{
						int id = iR(thisIndex);
						if (id == j)
						{
							thisIndex += soi;
							continue;
						}
						else
						{
							pW(bufIndex, id);
							pW(bufIndex, iR(thisIndex));
						}
					}
				}
				else //reduce factor of id == j
				{
					pW(bufIndex, facC);
					for (int f = 0; f < facC; f++)
					{
						int id = iR(thisIndex);
						if (id == j)
						{
							pW(bufIndex, id);
							pW(bufIndex, iR(thisIndex)-1);
						}
						else
						{
							pW(bufIndex, id);
							pW(bufIndex, iR(thisIndex));
						}
					}
				}

			}
		}
	}



	return mov(f);
}

#undef float

