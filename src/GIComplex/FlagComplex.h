/*
(c) 2012 Fengtao Fan
*/
#ifndef _FLAG_COMPLEX_H_
#define _FLAG_COMPLEX_H_
#include "SimpleGraph.h"
#include "SimplexNodeGIC.h"
#include "SimplicialComplexGIC.h"

#include <vector>
#include <set>
#include <map>
#include <algorithm>

class FlagComplex : public SimplicialTree<SimpleGraph>
{
public:
	FlagComplex()
	{
	}
	FlagComplex(const int in_dim, const int v_number, SimpleGraph* inDataPtr ) :
				SimplicialTree(in_dim, v_number, inDataPtr)
	{
	}
	bool Construction(); //virtual
};
//

//bool FlagComplex::Construction()
//{
//	// vertex has constructed in the initialization
//	//
//	// construct the edge
//	if (dim > 0)
//	{
//		int cur_dim_simplex_size = 0;
//		for (unsigned int i = 0; i < EuclideanDataPtr->vecNode.size(); i++)
//		{
//			Map_int_stnPtr_ptr curChildSetMap = vertex_array[i].get()->children_map_ptr;
//			SimpleGraphNode & gNode = EuclideanDataPtr->vecNode[i];
//			if (!gNode.adjList.empty())
//			{
//				Map_int_stnPtr_ptr pMap(new std::map<int, GICSimplicialTreeNode_ptr>);
//
//				for(std::set<int>::iterator sIter = gNode.adjList.begin(); 
//											sIter != gNode.adjList.end();
//											sIter++)
//				{
//					cur_dim_simplex_size++;// increase the # of edges
//					//
//					GICSimplicialTreeNode_ptr pStn(new GICSimplicialTreeNode(*sIter, // v_index
//																	vertex_array[i], // parent
//																	Map_int_stnPtr_ptr(),
//																	GICSimplicialTreeNode_ptr(),
//																	GICSimplicialTreeNode_ptr()
//																	)
//																	);
//					// link into the circular list of *sIter
//					InsertCircularList(head_circular_list_in_each_dim[*sIter][0], pStn);
//
//					// push it into the map
//					(*pMap.get())[*sIter] = pStn;
//				}
//				// set up the children set
//				curChildSetMap = pMap;
//			}
//		}
//		Simplicies_Cnt[1] = cur_dim_simplex_size;
//	}
//	//
//	//
//	if (dim > 1)
//	{
//		int cur_dim_simplex_size = 0;
//		int depth = 0;
//		for (int i = 1; i < dim; i++)
//		{
//			cur_dim_simplex_size = 0;
//			depth = i + 1; // want to visit all simplices at this dimension
//			if (Simplicies_Cnt[i] == 0)
//			{// no one lower dim simlex, means no current dim simplex
//				Simplicies_Cnt[i + 1] = 0;
//			}
//			else
//			{
//				for (int vid = 0; vid < Simplicies_Cnt[0]; vid++)
//				{
//					if (head_circular_list_in_each_dim[vid][depth - 1])
//					{// the circular list is not empty
//						GICSimplicialTreeNode_ptr pIter(head_circular_list_in_each_dim[vid][depth - 1]);
//						do
//						{// visit each simplex
//							// get previous (depth - 1) intersection set
//							std::set<int> prevIntersection;
//							int ver_index = pIter->last_label;
//							for (std::map<int, GICSimplicialTreeNode_ptr>::iterator 
//								mIter = pIter->parent_ptr->children_map_ptr->begin();
//								mIter != pIter->parent_ptr->children_map_ptr->end();
//								mIter++)
//							{
//								prevIntersection.insert(mIter->first);
//							}
//							//
//							unsigned int pot_size = std::min(prevIntersection.size(), EuclideanDataPtr->vecNode[ver_index].adjList.size());
//							std::vector<int> curIntersection(pot_size);
//							std::vector<int>::iterator interIter = std::set_intersection(prevIntersection.begin(),
//																			prevIntersection.end(),
//																			EuclideanDataPtr->vecNode[ver_index].adjList.begin(),
//																			EuclideanDataPtr->vecNode[ver_index].adjList.end(),
//																			curIntersection.begin());
//							if (interIter - curIntersection.begin() > 0)
//							{// create new simplicies 
//								Map_int_stnPtr_ptr pMap(new std::map<int, GICSimplicialTreeNode_ptr>);
//								for (int sid = 0; sid < interIter - curIntersection.begin(); sid++)
//								{
//									cur_dim_simplex_size++;// increase the # of edges
//									//
//									GICSimplicialTreeNode_ptr pStn(new GICSimplicialTreeNode(curIntersection[sid], // v_index
//																					pIter, // parent
//																					Map_int_stnPtr_ptr(),
//																					GICSimplicialTreeNode_ptr(),
//																					GICSimplicialTreeNode_ptr()
//																					)
//																					);
//									// link into the circular list of *sIter
//									InsertCircularList(head_circular_list_in_each_dim[curIntersection[sid]][depth - 1], pStn);
//
//									// push it into the map
//									(*pMap.get())[curIntersection[sid]] = pStn;
//								}
//								// set up the children set
//								pIter->children_map_ptr = pMap;
//							}
//							// move ahead
//							pIter = pIter->next_circular_ptr;
//						}while (pIter != head_circular_list_in_each_dim[vid][depth - 1]);
//					}
//				}
//				Simplicies_Cnt[i + 1] = cur_dim_simplex_size;
//			}//
//		}
//	}
//	return true;
//}
#endif //_FLAG_COMPLEX_H_