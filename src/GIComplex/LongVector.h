/*
(c) 2012 Fengtao Fan
*/
#ifndef MATH_LONGVECTOR_INCLUDED // -*- C++ -*-
#define MATH_LONGVECTOR_INCLUDED

/************************************************************************
  Long Vector class.
 ************************************************************************/
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
using namespace std;
class LongVector {
	#define FEQ_EPS 1e-6
	#define FEQ_EPS2 1e-12
private:
    float *elt;
 	int dim_;

protected:
    inline void copy(const LongVector& v);

public:
    //
    // Standard constructors
    //
	LongVector() : dim_(0), elt(NULL)
	{
	}
	LongVector(const int dim)
	{
		dim_ = dim;
		//if (dim_ > 0)
			elt = new float[dim_];
		//else
		//	elt = NULL;
	}
	LongVector(const float val, const int dim)
	{
		dim_ = dim;
		elt = new float[dim_];
		for (int i = 0; i < dim_; i++)
			elt[i] = val;
	}
    LongVector(float *inArray, const int dim)
	{
		dim_ = dim;
		if (dim_ > 0)
		{
			elt = new float[dim_];
			std::memcpy(elt, inArray, sizeof(float) * dim_);
		}
		else
		{
			elt = NULL;
		}

	}

    LongVector(const LongVector& v)
	{
		dim_ = v.dim_;
		elt = new float[dim_];
		memcpy(elt, v.elt, sizeof(float) * dim_);
	}
	~LongVector()
	{
		delete [] elt;
	}
    //
    // Access methods
    //

    float& operator()(int i)       { return elt[i]; }
    float  operator()(int i) const { return elt[i]; }

    float& operator[](int i)       { return elt[i]; }
    float  operator[](int i) const { return elt[i]; }

    float *raw()             { return elt; }
    const float *raw() const { return elt; }

    //
    // Comparison operators
    //
    inline bool operator==(const LongVector& v) const;
    inline bool operator!=(const LongVector& v) const;

    //
    // Assignment and in-place arithmetic methods
    //
    inline void set(float x, float y, float z) { elt[0]=x; elt[1]=y; elt[2]=z; }
    inline LongVector& operator=(const LongVector& v);
    inline LongVector& operator+=(const LongVector& v);
    inline LongVector& operator-=(const LongVector& v);
    inline LongVector& operator*=(float s);
    inline LongVector& operator/=(float s);

    //
    // Binary arithmetic methods
    //
    inline LongVector operator+(const LongVector& v) const;
    inline LongVector operator-(const LongVector& v) const;
    inline LongVector operator-() const;

    inline LongVector operator*(float s) const;
    inline LongVector operator/(float s) const;
    inline float operator*(const LongVector& v) const;
   // inline LongVector operator^(const LongVector& v) const;

	//
	// set individual component
	//
	inline void setI(unsigned int i, float s)
	{
		elt[i] = s;
	}
	inline void set(float *inArray)
	{
		for (int i = 0; i < dim_; i++)
			elt[i] = inArray[i];
	}
	inline int dim() const
	{
		return dim_;
	}
};


////////////////////////////////////////////////////////////////////////
//
// Method definitions
//

inline void LongVector::copy(const LongVector& v)
{
	if (dim_ != v.dim_)
	{
		delete [] elt;
		dim_ = v.dim_;
		elt = new float[dim_];
		memcpy(elt, v.elt, sizeof(float) * dim_);
	}
	else
		memcpy(elt, v.elt, sizeof(float) * dim_);
}

inline bool LongVector::operator==(const LongVector& v) const
{
	bool ret = false;
	float dx_sum = 0.0f;
	//if (dim_ != v.dim_)
	if (dim_ == v.dim_)
	{
		for (int i = 0; i < dim_; i++)
			dx_sum += (elt[i] - v.elt[i]) * (elt[i] - v.elt[i]);
		if (dx_sum < FEQ_EPS2)
			ret = true;
	}
	return ret;
}

inline bool LongVector::operator!=(const LongVector& v) const
{
	bool ret = true;
	float dx_sum = 0.0f;
	//if (dim_ != v.dim_)
	if (dim_ == v.dim_)
	{
		for (int i = 0; i < dim_; i++)
			dx_sum += (elt[i] - v.elt[i]) * (elt[i] - v.elt[i]);
		if (dx_sum < FEQ_EPS2)
			ret = false;
	}
	return ret;
}

inline LongVector& LongVector::operator=(const LongVector& v)
{
    copy(v);
    return *this;
}

inline LongVector& LongVector::operator+=(const LongVector& v)
{
	if (dim_ != v.dim_)
	{
		std::cout << "DIMENSION NOT MATCHED" << std::endl;
		exit(0);
	}
	for (int i = 0; i < dim_; i++)
		elt[i] += v[i];
	//elt[0] += v[0];   elt[1] += v[1];   elt[2] += v[2];
    return *this;
}

inline LongVector& LongVector::operator-=(const LongVector& v)
{
	if (dim_ != v.dim_)
	{
		std::cout << "DIMENSION NOT MATCHED" << std::endl;
		exit(0);
	}
	for (int i = 0; i < dim_; i++)
		elt[i] -= v[i];
    //elt[0] -= v[0];   elt[1] -= v[1];   elt[2] -= v[2];
    return *this;
}

inline LongVector& LongVector::operator*=(float s)
{
	for (int i = 0; i < dim_; i++)
		elt[i] *= s;
    //elt[0] *= s;   elt[1] *= s;   elt[2] *= s;
    return *this;
}

inline LongVector& LongVector::operator/=(float s)
{
	for (int i = 0; i < dim_; i++)
		elt[i] /= s;
	//elt[0] /= s;   elt[1] /= s;   elt[2] /= s;
    return *this;
}


inline LongVector LongVector::operator+(const LongVector& v) const
{
	if (dim_ != v.dim_)
	{
		std::cout << "DIMENSION NOT MATCHED" << std::endl;
		exit(0);
	}
	LongVector ret(dim_);
 	for (int i = 0; i < dim_; i++)
		ret[i] = elt[i] + v[i];
     return ret;
}

inline LongVector LongVector::operator-(const LongVector& v) const
{
	if (dim_ != v.dim_)
	{
		std::cout << "DIMENSION NOT MATCHED" << std::endl;
		exit(0);
	}
	LongVector ret(dim_);
 	for (int i = 0; i < dim_; i++)
		ret[i] = elt[i] - v[i];
     return ret;
   // return LongVector(elt[0]-v[0], elt[1]-v[1], elt[2]-v[2]);
}

inline LongVector LongVector::operator-() const
{
	LongVector ret(dim_);
 	for (int i = 0; i < dim_; i++)
		ret[i] = -elt[i];
     return ret;
	//return LongVector(-elt[0], -elt[1], -elt[2]);
}

inline LongVector LongVector::operator*(float s) const
{
	LongVector ret(dim_);
 	for (int i = 0; i < dim_; i++)
		ret[i] = elt[i] * s;
     return ret;
    //return LongVector(elt[0]*s, elt[1]*s, elt[2]*s);
}

inline LongVector LongVector::operator/(float s) const
{
	LongVector ret(dim_);
 	for (int i = 0; i < dim_; i++)
		ret[i] = elt[i] / s;
     return ret;
 //   return LongVector(elt[0]/s, elt[1]/s, elt[2]/s);
}

inline float LongVector::operator*(const LongVector& v) const
{
	if (dim_ != v.dim_)
	{
		std::cout << "DIMENSION NOT MATCHED" << std::endl;
		exit(0);
	}
 	float sum = 0.0f;
	for (int i = 0; i < dim_; i++)
		sum += elt[i] * v[i];
	return sum;
    //return elt[0]*v[0] + elt[1]*v[1] + elt[2]*v[2];
}

//inline LongVector LongVector::operator^(const LongVector& v) const
//{
//    LongVector w( elt[1]*v[2] - v[1]*elt[2],
//	   -elt[0]*v[2] + v[0]*elt[2],
//	    elt[0]*v[1] - v[0]*elt[1] );
//    return w;
//}

// Make scalar multiplication commutative
inline LongVector operator*(float s, const LongVector& v) { return v*s; }



////////////////////////////////////////////////////////////////////////
//
// Primitive function definitions
//

inline float norm(const LongVector& v)
{
 	float sum = 0.0f;
	for (int i = 0; i < v.dim(); i++)
		sum += v[i] * v[i];
    return sqrt(sum);
}

inline float norm2(const LongVector& v)
{
 	float sum = 0.0f;
	for (int i = 0; i < v.dim(); i++)
		sum += v[i] * v[i];
	return sum;
   // return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}
inline float unitize(LongVector& v)
{
    float l=norm2(v);
    if( l!=1.0f && l!=0.0f )
    {
	l = sqrtf(l);
	v /= l;
    }
    return l;
}
inline float length(const LongVector& v) { return norm(v); }

inline ostream& operator<<(ostream& out, const LongVector& v)
{
	out << "[";
	for (int i = 0; i < v.dim(); i++)
		out << v[i] << " ";
	out << "]";
    return out;// << "[" << v[0] << " " << v[1] << " " << v[2] << "]";
}

// MATH_LONGVECTOR_INCLUDED
#endif
