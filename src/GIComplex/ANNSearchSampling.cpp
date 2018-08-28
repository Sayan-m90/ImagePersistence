/*
(c) 2012 Fengtao Fan
*/ 
#include "ANNSearchSampling.h" 
#include "FibonacciHeap.h"
#include "DijkstraAlgorithm.h"

//namespace DASComp = AbstractSimplicialComplex;
namespace ANNSearch
{
	// randomly select a point from a point set
	int RandomSelect(const std::set<int>& intSet)
	{
		
		int r;
		int size = (int)intSet.size();
		r = rand();
		r = r % size;
		std::set<int>::const_iterator csIter = intSet.begin();
		for (int i = 0; i < r; i++)
			csIter++;
		return *csIter;
	} 
	//
	void Subsampling_GraphDistance( const float geoDiskRadius_delta,   
										const SimpleGraph &baseGraph,
										std::map<int, int> &SubsampleIndices,
										std::vector<int> &ColorMapping)
	{
		std::vector<int> scanned_vertex;
		std::vector<FIBO::Node*> vecNode; 
		std::vector<FIBO::Edge*> vecEdge;
		double segLength = 0.0;
		const int gNodeSize = (int)baseGraph.vecNode.size();
		// create vecNode
		for (int i = 0; i <gNodeSize; i++)
		{
			FIBO::Node* tempVer = new FIBO::Node(i, -1);
			vecNode.push_back(tempVer);
		}
		// create vecEdge
		for (int i = 0; i < gNodeSize; i++)
		{
			//FIBO::Edge* tmpEdge = new FIBO::Edge();
			for (std::map<int, float>::const_iterator mIter = baseGraph.vecNode[i].edgeWeights.begin();
				mIter != baseGraph.vecNode[i].edgeWeights.end();
				mIter++)
			{// edge (i, j) stored only once at vertex i ( i < j);
				int oppv = mIter->first;
				float w = mIter->second;
				//LongVector tempDiff;
				//tempDiff = (*ptSet._PointSet)[i] - (*ptSet._PointSet)[*sIter];
				//segLength = norm(tempDiff);
				FIBO::Edge* tempEdge = new FIBO::Edge(vecNode[i],  vecNode[oppv], w);
				tempEdge->head->addIncomingEdge(tempEdge);
				tempEdge->tail->addOutgoingEdge(tempEdge);
				vecEdge.push_back(tempEdge);

				tempEdge = new FIBO::Edge(vecNode[oppv], vecNode[i], w);
				tempEdge->head->addIncomingEdge(tempEdge);
				tempEdge->tail->addOutgoingEdge(tempEdge);
				vecEdge.push_back(tempEdge);
			}
		}
		//
		float halfGeoDiskRadius = geoDiskRadius_delta * 0.5f;
		// allocating landmarks randomly
 		std::set<int> uncolored_points; // unsampled vertex set
 
		// record current shortest path distance
		std::vector<double> shortestPathToLandmark;
		shortestPathToLandmark.reserve(vecNode.size());
		
		// Initialization
		for (unsigned int i = 0; i < vecNode.size(); i++){
			vecNode[i]->color = -1; // 0 means none assigment
			uncolored_points.insert(i);
			shortestPathToLandmark.push_back(0.0);
 		}
		
		//
 		int cur_pt_index = 0; // current landmark index
		int subsampling_index = 0; // counter for sub-sampling points
		enum FIBO::State unTouchedState, scannedState; // used for Dijkstra 
		// set up the appropriate states
		unTouchedState = FIBO::UNLABELED;
		scannedState = FIBO::SCANNED;

		////std::cout << "goes into while loop" << std::endl;
		while (!uncolored_points.empty())
		{
			cur_pt_index = RandomSelect(uncolored_points); // randomly select a vertex from these untouched vertices
 			//
			{// clean the memory of scanned_vertex
				std::vector<int> tmp;
				tmp.swap(scanned_vertex);
				scanned_vertex.clear();
			}
			DijkstraAlgorithm::DijkstraShortestPath(geoDiskRadius_delta, cur_pt_index, vecNode, 
									unTouchedState, scannedState,scanned_vertex);
			// record the shortest path distance
			//int visited_ver_cnt = 0;
			//for (unsigned int sortIter = 0; sortIter < vecNode.size(); sortIter++) 
			//{
			//	if (vecNode[sortIter]->state == FIBO::LABELED ||
			//		vecNode[sortIter]->state == FIBO::SCANNED) 
			//		visited_ver_cnt++;
			//}
			//if (visited_ver_cnt != scanned_vertex.size())
			//{
			//	std::cout << "ERroro in using scanned_ver" << std::endl;
			//	exit(0);
			//}
			//std::cout <<  "\r" << uncolored_points.size() << std::flush;
			// mark the nodes in radius
			for (unsigned int iscan = 0; iscan < scanned_vertex.size(); iscan++)
			{
				if (vecNode[scanned_vertex[iscan]]->state != FIBO::LABELED &&
					vecNode[scanned_vertex[iscan]]->state != FIBO::SCANNED)
				{
					std::cout << "ERROR in using scanned_vertex" << std::endl;
					exit(0);
				}
				if (vecNode[scanned_vertex[iscan]]->state == FIBO::SCANNED)
				{
					if (vecNode[scanned_vertex[iscan]]->color < 0)
					{// first time assignment
						vecNode[scanned_vertex[iscan]]->color = subsampling_index;
						shortestPathToLandmark[scanned_vertex[iscan]] = vecNode[scanned_vertex[iscan]]->key;
						//
						uncolored_points.erase(scanned_vertex[iscan]); 
					}
					else
					{
						if (shortestPathToLandmark[scanned_vertex[iscan]] > vecNode[scanned_vertex[iscan]]->key) 
						{
							shortestPathToLandmark[scanned_vertex[iscan]] = vecNode[scanned_vertex[iscan]]->key;
							vecNode[scanned_vertex[iscan]]->color = subsampling_index;
						}
					}
				}
				vecNode[scanned_vertex[iscan]]->state = unTouchedState;
				vecNode[scanned_vertex[iscan]]->key = 0;
			}
			 
			// record the sub sampling indices
			SubsampleIndices[cur_pt_index] = subsampling_index++; // increment the # of landmarks

 		}
		shortestPathToLandmark.clear();
		// delete data
		for (unsigned int i = 0; i < vecNode.size(); i++)
		{
			ColorMapping[i] = vecNode[i]->color;
			delete vecNode[i];
			vecNode[i] = 0;
		}
		for (unsigned int i = 0; i < vecEdge.size(); i++)
		{
			delete vecEdge[i];
			vecEdge[i] = 0;
		}
		return;
	}
	//
	void Subsampling_EuclideanDistance( const double sqDelta,
										const PointSet &ptSet,
										std::map<int, int> &SubsampleIndices,
										std::vector<int> &ColorMapping
										)
	{//
		//
		int nPts = (int)ptSet._PointSet->size();
		int dim = ptSet._dimension;
		//
		int k = 0; // parameter k in k-ann 
		//
		double eps = 0.0; // distance tolerance
		//
		ANNdist sqRadius = sqDelta;
		//
		SubsampleIndices.clear(); // initial empty subsampling
		//
		std::vector<double> curDistToSubsampling(nPts);
		//
		//initialize color mapping with negative values to indicate no mapping yet
		//
		memset(&(ColorMapping[0]), -1, sizeof(int) * nPts);
		//
		std::set<int> uncolored_points;
		//
		for (int i = 0; i < nPts; i++)
			uncolored_points.insert(i);
		// 
		int subsampling_index = 0;
		/*************/
		ANNpoint queryPt;
		ANNidxArray nnIdx = NULL;
		ANNdistArray dists = NULL;

		queryPt = annAllocPt(dim); // need to delloc
		
 		nnIdx = new ANNidx[nPts];
		dists = new ANNdist[nPts];
		/**************************/

		ANNpointArray dataPts;
		dataPts = annAllocPts(nPts, dim);
		// adding vertices
 
		for (int i = 0; i < nPts; i++){
			for (int j = 0; j < dim; j++)
			{
				dataPts[i][j] = (*ptSet._PointSet)[i][j]; 
			}
		}
		///////////////////
		ANNkd_tree* kdTree;
		kdTree = new ANNkd_tree(dataPts, 
								nPts,
								dim);
		//
		int NeighborCounter = 0;
		//
		// perform the sub-sampling here
		while (!uncolored_points.empty())
		{
			int cur_pt_idx = RandomSelect(uncolored_points);
			//
			// set up the query point as this point
			for (int j = 0; j < dim; j++)
			{
				queryPt[j] = (*ptSet._PointSet)[cur_pt_idx][j];
			}
			// get the number of points inside the fixed radius
			k = kdTree->annkFRSearch(
									queryPt,
									sqRadius,
									0, // to get the number of points within the radius
									nnIdx,
									dists,
									eps);
			NeighborCounter = kdTree->annkFRSearch(
											queryPt,
											sqRadius,
											k,
											nnIdx,
											dists,
											eps);
			// record edges
			if (NeighborCounter > k)
			{
				std::cout << "k is too small" << std::endl;
				exit(0);
			}
			if (NeighborCounter > 0)
			{// add edges
				if (nnIdx[0] != cur_pt_idx)
				{
					std::cout << "SELF is NOT contained" << std::endl;
					std::cout << "i " << cur_pt_idx << std::endl;
					std::cout << "nn " << nnIdx[0] << std::endl;
					std::cout << dists[0] << std::endl;
					exit(0);
				}
				else
				{
					ColorMapping[cur_pt_idx] = subsampling_index;
					uncolored_points.erase(cur_pt_idx);
				}
				for (int j = 1; j < NeighborCounter; j++)
				{
					if (ColorMapping[nnIdx[j]] < 0)
					{ // this point is not colored
						// color this point
						ColorMapping[nnIdx[j]] = subsampling_index;
						// remove it from uncolored vertex set
						uncolored_points.erase(nnIdx[j]);
						// record its distance to sub-sampling
						curDistToSubsampling[nnIdx[j]] = dists[j];
					}
					else
					{// this point is colored before
					// update its color based the distance
						if (dists[j] < curDistToSubsampling[nnIdx[j]] )
							// by default , if dists[j] == curDistToSubsampling[nnIdx[j]], we assign it with smaller index point
						{
							curDistToSubsampling[nnIdx[j]] = dists[j];
							ColorMapping[nnIdx[j]] = subsampling_index;
						}
					}
				}
			}
			else
			{
				std::cout << "isolated point !!!" << std::endl;
				exit(0);
			}
			//
			SubsampleIndices[cur_pt_idx] = subsampling_index++;
			//
			//std::cout << "i " << i << std::endl;
		}
 		// clean data
		//annDeallocPts(dataPts);

		delete [] nnIdx;
		delete [] dists;
		delete kdTree;
		if (dataPts)
			annDeallocPts(dataPts);
		if (queryPt)
			annDeallocPt(queryPt);
		annClose();
	}
											
	//
	void BuildDistanceGraph_ann( SimpleGraph& retGraph, 
								const PointSet &ptSet,
								const double sqDistance)
	{
		int nPts = (int)ptSet._PointSet->size();
		int dim =  ptSet._dimension;
		//
		int k = 2000;//nPts;
		if (k > nPts)
			k = nPts;
		double eps = 0.0;
		ANNdist sqRad = sqDistance;
		//
		retGraph.InitNodes(nPts);

		/*************/
		ANNpoint queryPt;
		ANNidxArray nnIdx = NULL;
		ANNdistArray dists = NULL;

		queryPt = annAllocPt(dim); // need to delloc
		
 		nnIdx = new ANNidx[nPts];
		dists = new ANNdist[nPts];
		/**************************/

		ANNpointArray dataPts;
		dataPts = annAllocPts(nPts, dim);
		// adding vertices
 
		for (int i = 0; i < nPts; i++){
			for (int j = 0; j < dim; j++)
			{
				dataPts[i][j] = (*ptSet._PointSet)[i][j]; 
			}
		}
		///////////////////
		ANNkd_tree* kdTree;
		kdTree = new ANNkd_tree(dataPts, 
								nPts,
								dim);
		//
		int NeighborCounter = 0;
		for (int i = 0; i < nPts; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				queryPt[j] = (*ptSet._PointSet)[i][j];
			}
		//
			// get the number of points inside the fixed radius
			k = kdTree->annkFRSearch(
									queryPt,
									sqRad,
									0, // to get the number of points within the radius
									nnIdx,
									dists);//,
									//eps); omitting eps para means exact search
			NeighborCounter = kdTree->annkFRSearch(
											queryPt,
											sqRad,
											k,
											nnIdx,
											dists,
											eps);
			// record edges
			if (NeighborCounter > k)
			{
				std::cout << "k is too small" << std::endl;
				exit(0);
			}
			if (NeighborCounter > 0)
			{// add edges
				if (nnIdx[0] != i)
				{
					std::cout << "SELF is NOT contained" << std::endl;
					std::cout << "i " << i << std::endl;
					std::cout << "nn " << nnIdx[0] << std::endl;
					std::cout << dists[0] << std::endl;
					exit(0);
				}
				for (int j = 1; j < NeighborCounter; j++)
				{
					//if (!retGraph.vecNode[i].is_neighbor(nnIdx[j]))
					if (!retGraph.is_edge(i, nnIdx[j]))
					{
						retGraph.AddEdge(i, nnIdx[j]);
 					}
				}
			}
			//std::cout << "i " << i << std::endl;
		}
 		// clean data
		//annDeallocPts(dataPts);

		delete [] nnIdx;
		delete [] dists;
		delete kdTree;
		if (dataPts)
			annDeallocPts(dataPts);
		if (queryPt)
			annDeallocPt(queryPt);
		annClose();
	} 

/************************************************/
 
//
	void SetColorMappingAndExtractColoredGraph(	const std::vector<int>& ColorMapping,
												const SimpleGraph& RipsGraph,
 												SimpleGraph& colorGraph							
							)
 	{
 		// iterate over the edges in RipsGraph
		//
		colorGraph.InitNodes((int)RipsGraph.vecNode.size());
		//
		for (int i = 0; i < (int) RipsGraph.vecNode.size(); i++)
		{
			// set color mapping
			colorGraph.vecNode[i].color = ColorMapping[i];
			//RipsGraph.vecNode[i].color = ColorMapping[i];
			//
			if (!RipsGraph.vecNode[i].edgeWeights.empty())
			{
				for (std::map<int, float>::const_iterator mIter = RipsGraph.vecNode[i].edgeWeights.begin();
					mIter != RipsGraph.vecNode[i].edgeWeights.end();
					mIter++)
				{
					if (ColorMapping[mIter->first] != ColorMapping[i])
					{ // bi-colored edge
						colorGraph.AddEdge(i, mIter->first);
					}
				}
			}
		} // for i size
	}

};
