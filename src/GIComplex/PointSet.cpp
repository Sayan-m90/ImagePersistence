/*
(c) 2012 Fengtao Fan
*/
#include "PointSet.h"
#include <fstream>
#include <iostream>
#include <sstream>

//#include <windows.h>
//#include <GL/gl.h>
#include <time.h>
#include <limits>
//
#include <set>
#include <algorithm>

PointSet::PointSet(const PointSet &rhs)
{
	_dimension = rhs._dimension;
	_size = rhs._size;
	_bdsMin = rhs._bdsMin;
	_bdsMax = rhs._bdsMax;

	_PointSet = new vector<LongVector> ;
	_PointSet->reserve(rhs._PointSet->size());
	_PointSet->assign(rhs._PointSet->begin(), rhs._PointSet->end());

	_sqDistToBasePt = rhs._sqDistToBasePt;

} 
void PointSet::SampleTorus(const float c, const float a, 
							const int sampleSize)//,
							//const int vSampleSize)
{
	/*
	x(u,v) = (c + a*cos(v))*cos(u)
	y(u,v) = (c + a*cos(v))*sin(u)
	z(u,v) = a*sin(v)
	*/
	//if (_PointSet)
	{
		delete[] _PointSet;
		_PointSet = 0;
	}
	_dimension = 3;
	float recArray[3] = {c, c, c};

	_bdsMax = LongVector(c+a, _dimension);
	_bdsMin = LongVector(-c-a, _dimension);
	_PointSet = new vector<LongVector> ;
	float u, v;

	const float local_PI = 3.141592f;
	float diffVal = (2.0f * local_PI) / sampleSize;
	u = 0.0f;
	for (int i = 0; i <  sampleSize; i++)
	{
		v = 0.0;
		for (int j = 0; j < sampleSize; j++)
		{
			recArray[0] = (c + a * cos(v)) * cos(u);
			recArray[1] = (c + a * cos(v)) * sin(u);
			recArray[2] =  a * sin(v) ;
			_PointSet->push_back( LongVector(recArray, 3) );
			v += diffVal;
		}
		u += diffVal;
	}

}
void PointSet::SampleSphere(float size, float radius)
{/*
 NOT FINISHED FUNCTION
 */
	std::cout << "NOT IMPLEMENTED " << std::endl;
	//if (_PointSet)
	{
		delete[] _PointSet;
		_PointSet = 0;
	}
	_dimension = 3;
	_bdsMax = LongVector(radius, _dimension);
	_bdsMin = LongVector(-radius, _dimension);
	_PointSet = new vector<LongVector> ;
	for (int i = 0; i <= size; i++)
		for (int j = 0; j < size; j++) {

				
			}
}
void PointSet::SampleCrossLine(double radius, double interval)
{
	// both radius and interval > 0
	//if (_PointSet)
	{
		delete[] _PointSet;
		_PointSet = 0;
	}
	_dimension = 2;
	float recArray[2] = {0};
	_bdsMax = LongVector(interval, _dimension);
	_bdsMin = LongVector(-interval, _dimension);

	int sample_num = (int)((interval) / radius);
	_PointSet = new vector<LongVector> ;
	
	_PointSet->push_back(LongVector(0.f, _dimension));
	for (int i = 0; i < sample_num; i++){
		recArray[0] = i * radius - interval; recArray[1] = 0;
		_PointSet->push_back(LongVector(recArray, 2));
		recArray[0] = interval - i * radius; recArray[1] = 0;
		_PointSet->push_back(LongVector(recArray, 2));
		recArray[0] = 0; recArray[1] = i * radius - interval;
		_PointSet->push_back(LongVector(recArray, 2));
		recArray[0] = 0; recArray[1] = interval - i * radius;
		_PointSet->push_back(LongVector(recArray, 2));
	}
	_size = _PointSet->size();
}
void PointSet::SampleCrossPlane(double radius, double interval)
{
	// both radius and interval > 0
	//if (_PointSet)
	{
		delete[] _PointSet;
		_PointSet = 0;
	}
	_dimension = 3;
	_bdsMax = LongVector(interval, _dimension);
	_bdsMin = LongVector(-interval, _dimension);

	//int sample_num = (int)((interval) / radius);
	//_PointSet = new vector<Vector3> ;
	//
	//_PointSet->push_back(Vector3(0, 0, 0));
	//for (int i = 0; i < sample_num; i++){
	//	_PointSet->push_back(Vector3(i * radius - interval, 0, 0));
	//	_PointSet->push_back(Vector3(interval - i * radius, 0, 0));
	//	_PointSet->push_back(Vector3(0, i * radius - interval, 0));
	//	_PointSet->push_back(Vector3(0, interval - i * radius, 0));
	//	_PointSet->push_back(Vector3(0, 0, i * radius - interval));
	//	_PointSet->push_back(Vector3(0, 0, interval - i * radius));

	//}
	//for (int i = 0; i < sample_num; i++){
	//	for (int j = 0; j < sample_num; j++){
	//		_PointSet->push_back(Vector3(i * radius - interval, j * radius - interval, 0));
	//		_PointSet->push_back(Vector3(i * radius - interval, interval - j * radius, 0));

	//		_PointSet->push_back(Vector3(interval - i * radius, j * radius - interval, 0));
	//		_PointSet->push_back(Vector3(interval - i * radius, interval - j * radius, 0));

	//		_PointSet->push_back(Vector3(i * radius - interval, 0, j * radius - interval));
	//		_PointSet->push_back(Vector3(i * radius - interval, 0, interval - j * radius));

	//		_PointSet->push_back(Vector3(interval - i * radius, 0, j * radius - interval));
	//		_PointSet->push_back(Vector3(interval - i * radius, 0, interval - j * radius));
	//	}
	//}
	_size = _PointSet->size();
}

void PointSet::ReadPointsFromFile(char const *pFileName)
{ 
	//
	std::cout << "Reading < " << pFileName << " > " << std::endl;
	bool firstPoint = true; 

	LongVector bdsMin, bdsMax;

	float *cVal;
	// allocate memory
	_PointSet = new vector<LongVector> ; 

	std::ifstream ifile;
	ifile.open(pFileName, std::ifstream::in);
  
	std::string sBuf; 
	//
	bool bFirstRead = true; 
	//
	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	std::stringstream locSstr(std::stringstream::in | std::stringstream::out);
	if (ifile.is_open())
	{
        /*
         * Get the size of the file
         */
        ifile.seekg(0,std::ios::end);
        long long iFileSize = ifile.tellg();
        ifile.seekg(0,std::ios::beg);
 
		//copy whole file into file buffer
		char* fBuf = new char[iFileSize + 1];
		ifile.read(fBuf, iFileSize);
		// add extra symbol
		fBuf[iFileSize] = '\n';
		//
		sBuf.assign(fBuf);
		sstr.str(sBuf);

		// close file 
		ifile.close();
		sBuf.clear();
		//deallocate memory
		delete [] fBuf; 
 

		//sstr << ifile.rdbuf();
		//sstr << "\n";
		while (sstr.good())
		{
			std::getline(sstr, sBuf);
			if (!sstr.good())
				break;
			//
			if (!sBuf.empty() && ( sBuf[0] >= '0' && sBuf[0] <= '9' || sBuf[0] == '-') )
			{

				//
				if (bFirstRead)
				{
					locSstr.str(sBuf);
					locSstr >> _dimension ;
 					if (_dimension <= 0) {
						std::cout << "Invalid dimension " << _dimension << std::endl;
						exit(0);
					}
					//
					cVal = new float[_dimension] ;
					//
					locSstr.clear();
					bFirstRead = false;
				}
				else
				{//
					LongVector tmp(_dimension);
					//
					locSstr.str(sBuf);

 					for (unsigned int i = 0; i < _dimension; i++)
						locSstr >> cVal[i];
					locSstr.clear();
					//
					if (firstPoint){
						bdsMax = LongVector(cVal, _dimension);
						bdsMin = bdsMax; 
						firstPoint = false;
					}
					else{
						for (unsigned k = 0; k < _dimension; k++){
							if (cVal[k] > bdsMax[k])
								bdsMax.setI(k, cVal[k]);
							if (cVal[k] < bdsMin[k])
								bdsMin.setI(k, cVal[k]);
						}
					} 
					tmp.set(cVal);

					_PointSet->push_back(tmp); 
				}
			}
  		}// while
		sstr.clear();
	}
	else
	{
		std::cout << "Can NOT open file " << pFileName << std::endl;
		exit(0);
	}
	_bdsMin = bdsMin;
	_bdsMax = bdsMax;
	_size = _PointSet->size();
	delete[] cVal;
	std::cout << "---- Done ---- V[" << _size << "]" << std::endl << std::endl;
	//for (int i = 0; i < (int)_PointSet->size(); i++)
	//{
	//	 for (unsigned int j = 0; j < _dimension; j++)
	//		std::cout << (*_PointSet)[i][j] << " " ;
	//	 std::cout << std::endl;
	//}
}
void PointSet::CreateSubset(PointSet &subPts, std::map<int, int> & subsampling_indices)
{
	//
	{
		delete[] subPts._PointSet;
		subPts._PointSet = 0;
	}
	// assign directly from parent point set
	subPts._bdsMax = _bdsMax;
	subPts._bdsMin = _bdsMin;
	//
	subPts._dimension = _dimension;
	subPts._size = subsampling_indices.size();
	//
	subPts._PointSet = new vector<LongVector>;
	//
	subPts._PointSet->resize(subPts._size);
	//
	for (std::map<int, int>::iterator mIter = subsampling_indices.begin();
		mIter != subsampling_indices.end();
		mIter++)
	{
		(*subPts._PointSet)[mIter->second] = (*_PointSet)[mIter->first];
	}
	return;
}
void PointSet::SubdivideSphere(const int depth, const double radius)
{
	// both radius and interval > 0
	//if (_PointSet)
	{
		delete[] _PointSet;
		_PointSet = 0;
	}
	_dimension = 3;
	float recArray[3];
	_bdsMax = LongVector(radius, _dimension);
	_bdsMin = LongVector(-radius, _dimension);

	_PointSet = new vector<LongVector> ;

	//// Initial tetrahedron
	//float init_vec3[] = {	0.f,       0.f,		1.f,		//v0
	//					    0.9428f,   0.f,		-0.3333f,   //v1
	//						-0.4714f,  0.8165f, -0.3333f,	//v2
	//						-0.4714f, -0.8165f, -0.3333f};	//v3
	////store initial points
	//for (int i = 0; i < 4; i++) {
	//	recArray[0] = init_vec3[3*i+0] * radius;
	//	recArray[1] = init_vec3[3*i+1]  * radius;
	//	recArray[2] = init_vec3[3*i+2]  * radius;
	//	_PointSet->push_back(LongVector(recArray, _dimension));
	//}
	//// edge contains 5 components
	//// 1st -- index of left endpoint L
	//// 2nd -- index of right endpoint R
	//// 3rd -- index of middle point M
	//// 4th -- index of edge [L M]
	//// 5th -- index of edge [M R]
	//int init_edges[] = {0, 1, -1, -1, -1, //e0
	//					0, 2, -1, -1, -1, //e1
	//					0, 3, -1, -1, -1, //e2
	//					1, 2, -1, -1, -1, //e3////
	//					1, 3, -1, -1, -1, //e4
	//					2, 3, -1, -1, -1}; //e5
	//					
	//int init_triangles[] = { 0, 1, 3,
	//						0, 2, 4,
	//						1, 2, 5,
	//						3, 4, 5};
	// Intial dual pyramids
	float init_vec3[] = {1.0f, 0.0f, 0.0f, //v0
						 0.0f, -1.0f, 0.0f, //v1
						 -1.0f, 0.0f, 0.0f, //v2
						 0.0f, 1.0f, 0.0f, //v3
						 0.0f, 0.0f, 1.0f, //v4
						 0.0f, 0.0f, -1.0f}; //v5
	//store initial points
	for (int i = 0; i < 6; i++) {
		recArray[0] = init_vec3[3*i+0] * radius;
		recArray[1] = init_vec3[3*i+1]  * radius;
		recArray[2] = init_vec3[3*i+2]  * radius;
		_PointSet->push_back(LongVector(recArray, _dimension));
	}
	 //edge contains 5 components
	 //1st -- index of left endpoint L
	 //2nd -- index of right endpoint R
	 //3rd -- index of middle point M
	 //4th -- index of edge [L M]
	 //5th -- index of edge [M R]
	int init_edges[] = {0, 1, -1, -1, -1, //e0
						1, 2, -1, -1, -1, //e1
						2, 3, -1, -1, -1, //e2
						3, 0, -1, -1, -1, //e3////
						0, 4, -1, -1, -1, //e4
						1, 4, -1, -1, -1, //e5
						2, 4, -1, -1, -1, //e6
						3, 4, -1, -1, -1,// e7////
						0, 5, -1, -1, -1, //e8
						1, 5, -1, -1, -1, //e9
						2, 5, -1, -1, -1, //e10
						3, 5, -1, -1, -1}; //e11
						
	int init_triangles[] = { 0, 4, 5,
							1, 5, 6,
							2, 6, 7,
							3, 7, 4,
							0, 8, 9,
							1, 9, 10,
							2, 10, 11,
							3, 11, 8};
	LongVector px, py;
	LongVector midxy;

	int leftEdge, rightEdge, commonVer, triEdgeIndx;


	int nVerNum = sizeof(init_vec3) / (3 * sizeof(float));
	int nTriNum = sizeof(init_triangles) / (3 * sizeof(int));
	int nEdgeNum = sizeof(init_edges) / (5 * sizeof(int));

	std::vector<int> triIndxVec(init_triangles, init_triangles + nTriNum * 3);
	std::vector<int> edgeIndxVec(init_edges, init_edges + 5 * nEdgeNum);
	std::vector<int> newEdgeIndx;

	int curTriNum = nTriNum;

	for (int subdDepth = 0; subdDepth < depth; subdDepth++) {
		newEdgeIndx.clear();
		curTriNum = triIndxVec.size() / 3;
		
		for ( int i = 0; i < curTriNum; i++){

			for (int jEdgeIndx = 0; jEdgeIndx < 3; jEdgeIndx++) {
				triEdgeIndx = triIndxVec[3 * i + jEdgeIndx];
				if (edgeIndxVec[5 * triEdgeIndx + 2] < 0){
					//middle point on this edge is NOT computed
					px = (*_PointSet)[edgeIndxVec[5 * triEdgeIndx]];
					py = (*_PointSet)[edgeIndxVec[5 * triEdgeIndx + 1]];

					midxy = (px + py) / 2.0;
					unitize(midxy);
					midxy = midxy * radius;
					
					_PointSet->push_back(midxy);
					edgeIndxVec[5 * triEdgeIndx + 2] = _PointSet->size() - 1;

					newEdgeIndx.push_back(edgeIndxVec[5 * triEdgeIndx]);
					newEdgeIndx.push_back(edgeIndxVec[5 * triEdgeIndx + 2]);
					newEdgeIndx.push_back(-1);
					newEdgeIndx.push_back(-1);
					newEdgeIndx.push_back(-1);
					edgeIndxVec[5 * triEdgeIndx + 3] = newEdgeIndx.size()/5 - 1;

					newEdgeIndx.push_back(edgeIndxVec[5 * triEdgeIndx + 2]);
					newEdgeIndx.push_back(edgeIndxVec[5 * triEdgeIndx + 1]);
					newEdgeIndx.push_back(-1);
					newEdgeIndx.push_back(-1);
					newEdgeIndx.push_back(-1);
					edgeIndxVec[5 * triEdgeIndx + 4] = newEdgeIndx.size()/5 - 1;
				}
			}
			for (int je = 0; je < 3; je++) {
				leftEdge = triIndxVec[3 * i + je ];
				rightEdge = triIndxVec[3 * i + (je + 1) % 3];
				if (edgeIndxVec[5 * leftEdge] == edgeIndxVec[5 * rightEdge] ||
					edgeIndxVec[5 * leftEdge] == edgeIndxVec[5 * rightEdge + 1])
					commonVer = edgeIndxVec[5 * leftEdge];
				else
					commonVer = edgeIndxVec[5 * leftEdge+1];

				newEdgeIndx.push_back(edgeIndxVec[5 * leftEdge + 2]);
				newEdgeIndx.push_back(edgeIndxVec[5 * rightEdge + 2]);
				newEdgeIndx.push_back(-1);
				newEdgeIndx.push_back(-1);
				newEdgeIndx.push_back(-1);

				triIndxVec.push_back(newEdgeIndx.size()/5 - 1);
				if (edgeIndxVec[5 * leftEdge] == commonVer) {
					triIndxVec.push_back(edgeIndxVec[5 * leftEdge + 3]);
				}
				else 
					triIndxVec.push_back(edgeIndxVec[5 * leftEdge + 4]);

				if (edgeIndxVec[5 * rightEdge] == commonVer) {
					triIndxVec.push_back(edgeIndxVec[5 * rightEdge + 3]);
				}
				else 
					triIndxVec.push_back(edgeIndxVec[5 * rightEdge + 4]);
			}
			triIndxVec[3 * i + 0] = newEdgeIndx.size() / 5 - 1;
			triIndxVec[3 * i + 1] = newEdgeIndx.size() / 5 - 2;
			triIndxVec[3 * i + 2] = newEdgeIndx.size() / 5 - 3;
		}
		edgeIndxVec.clear();
		edgeIndxVec.assign(newEdgeIndx.begin(), newEdgeIndx.end());
	}
	newEdgeIndx.clear();
	edgeIndxVec.clear();
	triIndxVec.clear();
	//
	//for (unsigned int i = 0; i < _PointSet->size(); i++)
	//{
	//	unitize((*_PointSet)[i]);
	//}
	//
	_size = _PointSet->size();
}
void PointSet::ResampleByLandmarks(std::map<int, int> newToOldLandmarks,
							 PointSet& newPts)
{
	newPts._bdsMax = _bdsMax;
	newPts._bdsMin = _bdsMin;
	newPts._size = newToOldLandmarks.size();
	newPts._dimension = _dimension;
	if (newPts._PointSet)
	{
		{// DESTROY THE MEMORY HOLDINGN PCD 
			std::vector<LongVector> tmp;
			newPts._PointSet->swap(tmp);
		}
		delete newPts._PointSet;
	}
	//
	newPts._PointSet = new vector<LongVector>;
	for (unsigned int i = 0; i < newToOldLandmarks.size(); i++)
		newPts._PointSet->push_back((*_PointSet)[newToOldLandmarks[i]]);
	//for (std::map<int, int>::const_iterator cmIter = newToOldLandmarks.begin();
	//	cmIter != newToOldLandmarks.end();
	//	cmIter++)
	//{
	//	newPts._PointSet->push_back((*_PointSet)[cmIter->second]);
	//}
}

void PointSet::SampleCircleProductSphere(const int xgridSize, 
							   const int ygridSize, 
							   const int zgridSize, 
							   const float radius)
{// sample in grid x grid x grid
	float theta = 0.0f; // angle for circle
	float alpha = 0.0f; // angle for sphere
	float beta  = 0.0f; // angle for sphere
	const float local_PI = 3.141592f;
	float dx = (2.0f * local_PI) / xgridSize;
	float dy = (2.0f * local_PI) / ygridSize;
	float dz = (2.0f * local_PI) / zgridSize;

	_dimension = 5; // first two for circle, other 3 for sphere
	_PointSet = new vector<LongVector>;
  /* initialize random seed: */
	srand ( time(NULL) );
	const int max_int = std::numeric_limits<int>::max() ;
	const int max_int_f =  (max_int) ;
	float recArray[5];
	for (int i = 0; i < xgridSize; i++)
	{
		for (int j = 0; j < ygridSize; j++)
		{
			for (int k = 0; k < zgridSize; k++)
			{
				//randomly pick one point in
				//[0 dx] x [0 dy] x [0 dz]
				float rx = 0;
				float ry = 0;
				float rz = 0;
				rx =  float(rand());
				rx = rx /  ( float (max_int));
				ry =  rand()   ;
				ry = ry /  ( max_int);
				rz =  rand();
				rz = rz /  ( max_int);
				recArray[0] = radius * cosf(theta + rx * dx);
				recArray[1] = radius * sinf(theta + rx * dx);
				//
				recArray[2] = radius * cosf(alpha + ry * dy) * cosf(beta + rz * dz);
				recArray[3] = radius * cosf(alpha + ry * dy) * sinf(beta + rz * dz);
				recArray[4] = radius * sinf(alpha + ry * dy);
				//
				_PointSet->push_back(LongVector(recArray, 5));
				//
				beta += dz;
			}
			alpha += dy;
		}
		theta += dx;
	}
	_size = _PointSet->size();
}
void PointSet::SampleCircleProductCircle(const int xgridSize, 
							   const int ygridSize, 
							   const float radius)
{// sample in grid x grid x grid
	float theta = 0.0f; // angle for circle_A
	float alpha = 0.0f; // angle for circle_B 
	const float local_PI = 3.141592f;
	float dx = (2.0f * local_PI) / xgridSize;
	float dy = (2.0f * local_PI) / ygridSize; 

	_dimension = 4; // first two for circle_A, other 2 for circle_B
	_PointSet = new vector<LongVector>;
  /* initialize random seed: */
	//srand ( time(NULL) );
	const int max_int = std::numeric_limits<int>::max() ;
	const int max_int_f =  (max_int) ;
	float recArray[4];
	for (int i = 0; i < xgridSize; i++)
	{
		for (int j = 0; j < ygridSize; j++)
		{ 
				//randomly pick one point in
				//[0 dx] x [0 dy]  
				float rx = 0;
				float ry = 0; 
				rx =  float(rand());
				rx = rx /  ( float (max_int));
				ry =  rand()   ;
				ry = ry /  ( max_int); 

				recArray[0] = radius * cosf(theta);// cosf(theta + rx * dx);
				recArray[1] = radius * sinf(theta);
				//
				recArray[2] = radius * cosf(alpha);//cosf(alpha + ry * dy) * cosf(beta + rz * dz);
				recArray[3] = radius * sinf(alpha);// + ry * dy) * sinf(beta + rz * dz); 
				//
				_PointSet->push_back(LongVector(recArray, _dimension));
				//
			 
			alpha += dy;
		}
		theta += dx;
	}
	_size = _PointSet->size();
}
void PointSet::SampleEllipticalSphere(const float a,
								const float b,
								const float c,
								const int uGridSize,
								const int vGridSize)
{
	const float local_PI = 3.141592f;
	const float du = (local_PI) / uGridSize;
	const float dv = (2.0f * local_PI) / vGridSize; 
	float u = -local_PI / 2.0f;
	float v = -local_PI;
	float recArray[3];
	//
	if (a < 0.f || b < 0.f || c < 0.f)
	{
		std::cout << "Elliptical axis ZERO" << std::endl;
		exit(0);
	}
	_dimension = 3;
	_PointSet = new vector<LongVector>;

	for (int i = 0; i <= uGridSize; i++)
	{
		for (int j = 0; j <= vGridSize; j++)
		{
			recArray[0] = a * cosf(u) * cosf(v);
			recArray[1] = b * cosf(u) * sinf(v);
			recArray[2] = c * sinf(u);
			_PointSet->push_back(LongVector(recArray, _dimension));
			//
			v += dv;
		}
		u += du;
	}
	//
	_size = _PointSet->size();
	//
	_bdsMin = LongVector(recArray, _dimension);
	_bdsMax = LongVector(recArray, _dimension);
	//
	_bdsMax[0] = a;
	if (_bdsMax[0] < b)
		_bdsMax[0] = b;
	if (_bdsMax[0] < c)
		_bdsMax[0] = c;
	_bdsMax[1] = _bdsMax[2] = _bdsMax[0];
	_bdsMin[0] = _bdsMin[1] = _bdsMin[2] = -_bdsMax[0];
	return;
}
void PointSet::WriteColorMappingBackToFile(char const* pFileName, const int color_num, std::vector<int> &ColorMapping)
{
/*
dimension
coordinates
*/
	std::cout << "wrting... " << std::endl;
 
	std::ofstream ofile ;
	ofile.open(pFileName); 

  	//
 	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		sstr << color_num << std::endl;
		for (unsigned int i = 0; i < ColorMapping.size(); i++)
		{
			sstr << ColorMapping[i] << " ";
		}
		sstr << std::endl;
		//
		//ofile << sstr.rdbuf();
		ofile.write(sstr.str().c_str(), sstr.str().size());
		//
		ofile.close();
		sstr.clear();
	}
	else
	{
		std::cout << "Can NOT open file " << pFileName << std::endl;
		exit(0);
	}
	std::cout << "Done... " << color_num << std::endl;
	// 
	return;
}
void PointSet::WriteBackToFile(char const * pFileName)
{
/*
dimension
coordinates
*/
	std::cout << "wrting... " << pFileName <<  std::endl;
 
	std::ofstream ofile ;
	ofile.open(pFileName); 

  	//
 	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		sstr << _dimension << std::endl;
		for (unsigned int i = 0; i < _PointSet->size(); i++)
		{
			for (unsigned int k = 0; k < _dimension; k++)
				sstr << (*_PointSet)[i][k] << " ";
			sstr << std::endl;
		}
		//
		//ofile << sstr.rdbuf();
		ofile.write(sstr.str().c_str(), sstr.str().size());
		//
		ofile.close();
		sstr.clear();
	}
	else
	{
		std::cout << "Can NOT open file " << pFileName << std::endl;
		exit(0);
	}
	std::cout << "Done... " << _size << std::endl;
	// 
	return;
}
struct LongVectorLessThan
{
	bool operator()(LongVector const &lhs, LongVector const &rhs) const
	{
		bool ret = false;
		for (int i = 0; i < lhs.dim(); i++)
		{
			if (lhs[i] < rhs[i])
			{
				ret = true;
				break;
			}
			else
			{
				if (lhs[i] > rhs[i])
				{
					ret = false;
					break;
				}
			}
		}
		return ret;
	}
};
void PointSet::ClearDuplicatedCopy(char const *pFileName)
{
	std::set<LongVector, LongVectorLessThan> cleanCopy;
	for (unsigned int i = 0; i < _size; i++)
	{
		cleanCopy.insert((*_PointSet)[i]);
	}
	if (cleanCopy.size() < _size)
	{
		std::cout << "DUPLICATED COPIES FOUND" << std::endl;
		{
			//clean the memory of vector of longvector
			std::vector<LongVector> tmpVec;
			_PointSet->swap(tmpVec);
		}
		_PointSet->resize(cleanCopy.size());
		std::copy(cleanCopy.begin(), cleanCopy.end(), _PointSet->begin());//std::inserter(*_PointSet,_PointSet->begin()));
		_size = _PointSet->size();
	}
	else
	{
		std::cout << "THERE IS NO COPIES" << std::endl;
	}
	WriteBackToFile(pFileName);
}