/*
(c) 2012 Fengtao Fan
*/
#ifndef _DIJKSTRA_ALGORITHM_H_
#define _DIJKSTRA_ALGORITHM_H_

#include <iostream>
#include <vector>
#include "FibonacciHeap.h"

namespace FIBO = FibonacciHeap;
namespace DijkstraAlgorithm
{
 
	void DijkstraShortestPath(	const double geoDiskRadius,
								const int srcIndex,
 								std::vector<FIBO::Node*> &vecNode,  
								enum FIBO::State unTouchedState,
								enum FIBO::State scannedState,
								std::vector<int> &scanned_vertex); 

}
#endif //_DIJKSTRA_ALGORITHM_H_
