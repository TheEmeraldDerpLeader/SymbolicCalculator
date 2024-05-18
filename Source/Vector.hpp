#pragma once
#include <vector>

template<typename T>
class Vector2D;
template<typename T>
class VectorRef;
template<typename T>
class VectorOffRef;
template<typename T>
class Vector2DRef;


template<typename T>
class Vector : public std::vector<T>
{
public:
	Vector<T>(int size)
	{
		this->resize(size);
	}

	Vector<T>() //necessary since one overriden constructor invalidates other constructors
	{

	}
	Vector<T>(std::vector<T>& vCopy) : std::vector<T>(vCopy)
	{
	
	}
	Vector<T>(std::vector<T>&& vMove) : std::vector<T>(vMove)
	{

	}
	Vector<T>& operator=(std::vector<T>& vCopy) 
	{
		*static_cast<std::vector<T>*>(this) = vCopy;
		return *this;
	}
	Vector<T>& operator=(std::vector<T>&& vMove) 
	{
		*static_cast<std::vector<T>*>(this) = vMove;
		return *this;
	}
};

template<typename T>
class Vector2D : public std::vector<T>
{
public:
	int width = 0;
	int height = 0;
	
	Vector2D<T>(int w, int h)
	{
		SetDim(w,h);
	}

	Vector2D<T>() //necessary since one overriden constructor invalidates other constructors
	{

	}
	Vector2D<T>(std::vector<T>& vCopy) : std::vector<T>(vCopy)
	{

	}
	Vector2D<T>(std::vector<T>&& vMove) : std::vector<T>(vMove)
	{

	}
	Vector2D<T>& operator=(std::vector<T>& vCopy) 
	{
		*static_cast<std::vector<T>*>(this) = vCopy;
		return *this;
	}
	Vector2D<T>& operator=(std::vector<T>&& vMove) 
	{
		*static_cast<std::vector<T>*>(this) = vMove;
		return *this;
	}
	
	void SetDim(int w, int h)
	{
		width = w;
		height = h;
		this->resize(w*h);
	}

	T& At(int x, int y) { return (*this)[(x*height)+y]; }

	VectorRef<T> GetCol(int x)
	{
		return VectorRef<T>(height, this->data()+(height*x));
	}
	VectorOffRef<T> GetRow(int y)
	{
		return VectorOffRef<T>(width, height, this->data()+y);
	}
};


template<typename T>
struct VectorRef
{
private:
public:
	int length = 0;
	T* data = nullptr;
	inline int size() {return length;}

	VectorRef(int lengthInt, T* dataPtr)
	{
		length = lengthInt;
		data = dataPtr;
	}
	VectorRef(Vector<T>& v)
	{
		length = v.size();
		data = v.data();
	}
	VectorRef(Vector2DRef<T>& v)
	{
		length = v.width*v.height;
		data = v.data;
	}

	T& operator[](int index) { return data[index]; }

	//Vector<T> operator*(const Vector2DRef<T>& v) const;
};

template<typename T>
struct VectorOffRef
{
private:
public:
	int length = 0;
	int offset = 1;
	T* data = nullptr;
	inline int size() {return length;}

	VectorOffRef(int lengthInt, int offsetInt, T* dataPtr)
	{
		length = lengthInt;
		offset = offsetInt;
		data = dataPtr;
	}

	T& operator[](int index) { return data[offset*index]; }

	//Vector<T> operator*(const Vector2DRef<T>& v) const;
};

template<typename T>
struct Vector2DRef
{
public:
	int width = 0;
	int height = 0;
	T* data = nullptr;
	inline int size() {return width*height;}

	Vector2DRef(int w, int h, T* dataPtr)
	{
		width = w;
		height = h;
		data = dataPtr;
	}
	Vector2DRef(Vector2D<T>& v)
	{
		width = v.width;
		height = v.height;
		data = v.data();
	}

	T& operator[](int index) { return data[index]; }

	T& At(int x, int y) { return data[(x*height)+y]; }

	VectorRef<T> GetCol(int x)
	{
		return VectorRef<T>(height, data+(height*x));
	}
	VectorOffRef<T> GetRow(int y)
	{
		return VectorOffRef<T>(width, height, data+y);
	}
};


template<typename T, template<typename> typename T1, template<typename> typename T2>
Vector<T> Add(T1<T>& v1, T2<T>& v2)
{
	int l = v1.size();
	Vector<T> hold; hold.resize(l);
	for (int i = 0; i < l; i++)
		hold[i] = v1[i]+v2[i];
	return hold;
}
template<typename T, template<typename> typename T1, template<typename> typename T2, template<typename> typename T3>
void Add(T1<T>& v1, T2<T>& v2, T3<T>& out)
{
	int l = v1.size();
	for (int i = 0; i < l; i++)
		out[i] = v1[i]+v2[i];
}
template<typename T, template<typename> typename T1, template<typename> typename T2>
void AddEq(T1<T>& v1, T2<T>& v2)
{
	int l = v1.size();
	for (int i = 0; i < l; i++)
		v1[i] += v2[i];
}
template<typename T, template<typename> typename T1, template<typename> typename T2>
Vector<T> Sub(T1<T>& v1, T2<T>& v2)
{
	int l = v1.size();
	Vector<T> hold; hold.resize(l);
	for (int i = 0; i < l; i++)
		hold[i] = v1[i]-v2[i];
	return hold;
}
template<typename T, template<typename> typename T1, template<typename> typename T2, template<typename> typename T3>
void Sub(T1<T>& v1, T2<T>& v2, T3<T>& out)
{
	int l = v1.size();
	for (int i = 0; i < l; i++)
		out[i] = v1[i]-v2[i];
}
template<typename T, template<typename> typename T1, template<typename> typename T2>
void SubEq(T1<T>& v1, T2<T>& v2)
{
	int l = v1.size();
	for (int i = 0; i < l; i++)
		v1[i] -= v2[i];
}
template<typename T, template<typename> typename T1>
Vector<T> Mul(T1<T>& v1, T& s)
{
	int l = v1.size();
	Vector<T> hold; hold.resize(l);
	for (int i = 0; i < l; i++)
		hold[i] = v1[i]*s;
	return hold;
}
template<typename T, template<typename> typename T1, template<typename> typename T2>
void Mul(T1<T>& v1, T& s, T2<T>& out)
{
	int l = v1.size();
	for (int i = 0; i < l; i++)
		out[i] = v1[i]*s;
}
template<typename T, template<typename> typename T1>
void MulEq(T1<T>& v1, T& s)
{
	int l = v1.size();
	for (int i = 0; i < l; i++)
		v1[i] *= s;
}

template<typename T, template<typename> typename T1>
T SumAbs(T1<T>& v1) //will error if glm::abs(T) doesn't work
{
	T sum = 0;
	int l = v1.size();
	for (int i = 0; i < l; i++)
		sum += glm::abs(v1[i]);
	return sum;
}
template<typename T, template<typename> typename T1, template<typename> typename T2>
T Dot(T1<T>& v1, T2<T>& v2)
{
	T sum = T();
	int l = v1.size();
	for (int i = 0; i < l; i++)
		sum += v1[i]*v2[i];
	return sum;
}
template<typename T, template<typename> typename T1, template<typename> typename T2>
Vector<T> CompMul(T1<T>& v1, T2<T>& v2)
{
	int l = v1.size();
	Vector<T> hold; hold.resize(l);
	for (int i = 0; i < l; i++)
		hold[i] = v1[i]*v2[i];
	return hold;
}

/*template<typename T>
inline Vector<T> VectorRef<T>::operator*(const Vector2DRef<T>& v) const
{
	Vector<T> hold; hold.resize(v.width);
	for (int i = 0; i < v.width; i++)
	{
		T sum = T();
		for (int j = 0; j < v.height; j++)
			sum += (*this)[j] * v[i][j];
		hold[i] = sum;
	}
	return hold;
}*/