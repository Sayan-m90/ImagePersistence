/*
(c) 2012 Fengtao Fan
*/
#ifndef _ANN_SAMPLING_H_
#define _ANN_SAMPLING_H_

#include <vector>

#include <ANN/ANN.h>					// ANN declarations
#include <ANN/ANNx.h>					// more ANN declarations
#include <ANN/ANNperf.h>				// performance evaluation

#include "SimpleGraph.h"
#include "PointSet.h"  

namespace ANNSearch
{
	//DASComp::_simpGraph 
	void BuildDistanceGraph_ann( SimpleGraph& retGraph, 
								const PointSet &ptSet,
								const double sqDistance); 
/**********************************************************/ 
	// Euclidean distance subsampling
	void Subsampling_EuclideanDistance( const double sqDelta,
										const PointSet &ptSet,
										std::map<int, int> &SubsampleIndices,
										std::vector<int> &ColorMapping
										);
/**********************************************************/

/**********************************************************/
	void Subsampling_GraphDistance( const float geoDiskRadius_delta,   
									const SimpleGraph &baseGraph,
									std::map<int, int> &SubsampleIndices,
									std::vector<int> &ColorMapping);
/**********************************************************/

	// get the bi-colored graph
	void SetColorMappingAndExtractColoredGraph(	const std::vector<int>& ColorMapping,
								const SimpleGraph& RipsGraph,
 								SimpleGraph& colorGraph							
							);  
};
#endif //_ANN_SAMPLING_H_
