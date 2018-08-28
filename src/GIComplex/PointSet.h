/*
(c) 2012 Fengtao Fan
*/
#ifndef _POINT_SET_IN_ANY_DIM_H_
#define _POINT_SET_IN_ANY_DIM_H_

#include "LongVector.h"
#include <vector>
#include <map>

class PointSet
{
	/*
	POINT CLOUD DATA FILE FORMAT
	LINE 1 : DIMENSION
	FOLLOWING LINES ARE POINTS.
	EACH POINT OCCUPIES A LINE.

	EXAMPLE: (PCD IN R^3 WHICH ONLY HAS ONE POINT)
	3
	0.0 1.0 0.0
	*/
/*
	POINT CLOUD DATA (PCD) CLASS FOR HIGH DIMENSTION POINTS:

	DATA STRUCTURES:

	LongVector : n-dimension point class

	DATA MEMBERS:
	_bdsMin		:	left bottom corner of the bounding box in R^n
	_bdsMax		:	right upper corner of the bounding box in R^n
	_size		:	the # of points in PCD;
	_dimension	:	the dimension of points = n (in R^n);
	_PointSet	:	pointer to a vector container of LongVector which contains all points
	_sqDistToBasePt : reserved for local homology

*/
public:
	LongVector _bdsMin;
	LongVector _bdsMax;
	unsigned int _size;
	unsigned int _dimension;
	vector<LongVector> *_PointSet;
	vector<double> _sqDistToBasePt;
	/*
	CONSTRUCTOR AND DECONSTRUCTOR
	*/
	/*
	DEFAULT CONSTRUCTOR
	*/
	PointSet() : _size(0), _dimension(0), _PointSet(NULL)
	{
		_bdsMin = LongVector(0);
		_bdsMax = LongVector(0);
	}
	/*
	CONSTRUCTOR BY REFERENCE
	*/
	PointSet(const PointSet &rhs);
	/*
	CONSTRUCT FROM A PCD FILE
	*/
	PointSet(char const * inFileName)
	{
		ReadPointsFromFile(inFileName);
	}
	/*
	DECONSTRUCTOR
	*/
	~PointSet()
	{
		{// DESTROY THE MEMORY HOLDINGN PCD
			std::vector<LongVector> tmp;
			_PointSet->swap(tmp);
		}
		delete _PointSet;
		_PointSet = NULL;
	}
	//
	void SampleTorus(const float c, const float a, const int sampleSize);
	void SampleSphere(float sampleSize, float radius);
	void SampleCrossLine(double radius, double interval);
	void SampleCrossPlane(double radius, double interval);
	void SubdivideSphere(const int depth, const double radius);
	void SampleCircleProductSphere(const int xgridSize,
							   const int ygridSize,
							   const int zgridSize,
							   const float radius);
	void SampleEllipticalSphere(const float a,
								const float b,
								const float c,
								const int uGridSize,
								const int vGridSize);
	void SampleCircleProductCircle(const int xgridSize,
								   const int ygridSize,
								   const float radius);
	/*
	READ POINTS FROM A PCD FILE
	*/
	void ReadPointsFromFile(char const * fileName);

	/*
	EXTRACT THE SUBSET OT PCD ACCORDING TO THE LANDMARKS
	*/
	void ResampleByLandmarks(const std::map<int, int> newToOldLandmarks,
							 PointSet& newPts);

	/*
	WRITE THE PCD TO A FILE
	*/
	void WriteBackToFile(char const* pFileName);

	/*
	WRITE COLOR MAPPING TO A FILE
	*/
	void WriteColorMappingBackToFile(char const* pFileName, const int color_num, std::vector<int> &ColorMapping);

	/*
	CLEAN DUPLICATED POINTS
	THAT IS, TO MAKE POINTS ARE UNIQUE
	*/
	void ClearDuplicatedCopy(char const *pFileName);
	//
	void CreateSubset(PointSet &subPts, std::map<int, int> & subsampling_indices);

private:


};
#endif //_POINT_SET_IN_ANY_DIM_H_
