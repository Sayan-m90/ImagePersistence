/*
(c) 2012 Fengtao Fan
*/
#include "DijkstraAlgorithm.h"
#include <algorithm>

namespace FIBO = FibonacciHeap;
namespace DijkstraAlgorithm
{
 
	void DijkstraShortestPath(	const double geoDiskRadius,
								const int srcIndex,
								std::vector<FIBO::Node*>& vecNode,
								enum FIBO::State unTouchedState,
								enum FIBO::State scannedState,
								std::vector<int> &visited_vertex)
	{// only for connected graph
	 
		FIBO::FibonacciHeap* heap = new FIBO::FibonacciHeap();

		vecNode[srcIndex]->state = FIBO::LABELED;
		vecNode[srcIndex]->key   = 0;

		heap->insertVertex(vecNode[srcIndex]);
			
		//
		visited_vertex.push_back(srcIndex);
		// Scan
		do
		{
			// Delete minimum path
			FIBO::Node* v = heap->deleteMin();

			if (v->key > geoDiskRadius)
				break; // because Dijkstra algorithm's output is in increasing order
						// thus no need to consider the distance greater than geoDiskRadius
				
			v->state = scannedState;//SCANNED;
				
			for(unsigned int j = 0; j < v->incomingEdges.size(); j++)
			{
				FIBO::Edge* currentEdge = v->incomingEdges[j];
				FIBO::Node* headOfCurrentEdge = currentEdge->tail;

				if(headOfCurrentEdge->state != scannedState)//SCANNED)
				{
					if(headOfCurrentEdge->state == unTouchedState)//UNLABELED)
					{
						// Insert a vertex with infinite key
						headOfCurrentEdge->state = FIBO::LABELED;
						headOfCurrentEdge->pred = v;
						headOfCurrentEdge->key = v->key + currentEdge->length;
						heap->insertVertex(headOfCurrentEdge);
						// record the visited vertex
						visited_vertex.push_back(headOfCurrentEdge->data);
					}
					else if(headOfCurrentEdge->key > v->key + currentEdge->length )
					{
						// decrease the key of a vertex with finite key
						headOfCurrentEdge->pred = v;
						heap->decreaseKey(v->key + currentEdge->length, headOfCurrentEdge);
					}
				}
			}
		} while(!heap->isEmpty());	

	while (!heap->isEmpty())
		heap->deleteMin();
	// record the distance
	//for (unsigned int col = 0; col < vecNode.size(); col++){
	//	wMat[col] = vecNode[col]->key;
	//}
	}
 

}