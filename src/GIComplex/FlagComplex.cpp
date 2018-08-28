/*
(c) 2012 Fengtao Fan
*/
#include "FlagComplex.h"

#include <set>
#include <map>
#include <algorithm>

bool FlagComplex::Construction()
{// index the simplicies by the dimension order
	int simplex_size = Simplicies_Cnt[0];
	// vertex has constructed in the initialization
	//
	// construct the edge
	if (dim > 0)
	{
		int cur_dim_simplex_size = 0;
		for (unsigned int i = 0; i < EuclideanDataPtr->vecNode.size(); i++)
		{
			Map_int_stnPtr_ptr& curChildSetMap = vertex_array[i].get()->children_map_ptr;
			SimpleGraphNode & gNode = EuclideanDataPtr->vecNode[i];
			if (!gNode.edgeWeights.empty())
			{	//
				Map_int_stnPtr_ptr pMap(new std::map<int, GICSimplicialTreeNode_ptr>);
				//
				for(std::map<int, float>::iterator mIter = gNode.edgeWeights.begin(); 
											mIter != gNode.edgeWeights.end();
											mIter++)
				{
					cur_dim_simplex_size++;// increase the # of edges
					//
					GICSimplicialTreeNode_ptr pStn(new GICSimplicialTreeNode(mIter->first, // v_index
																	simplex_size++, // order in the filtration
																	vertex_array[i], // parent
																	Map_int_stnPtr_ptr(),
																	GICSimplicialTreeNode_ptr() 
																	)
																	);
					// link into the circular list of *sIter
					InsertCircularList(head_circular_list_in_each_dim[mIter->first][0], pStn);
					//
					// push it into the map
					(*pMap.get())[mIter->first] = pStn;
				}
				// set up the children set
				curChildSetMap = pMap;
			}
		}
		Simplicies_Cnt[1] = cur_dim_simplex_size;
	}
	//
	//
	if (dim > 1)
	{
		int cur_dim_simplex_size = 0;
		int depth = 0;
		for (int i = 1; i < dim; i++)
		{
			cur_dim_simplex_size = 0;
			depth = i + 1; // want to visit all simplices at this dimension
			if (Simplicies_Cnt[i] == 0)
			{// no one lower dim simlex, means no current dim simplex
				Simplicies_Cnt[i + 1] = 0;
			}
			else
			{
				for (int vid = 0; vid < Simplicies_Cnt[0]; vid++)
				{
					if (head_circular_list_in_each_dim[vid][depth - 2])
					{// the circular list is not empty
						GICSimplicialTreeNode_ptr pIter(head_circular_list_in_each_dim[vid][depth - 2]);
						do
						{// visit each simplex
							// get previous (depth - 1) intersection set
							//std::set<int> prevIntersection;
							std::vector<int> prevIntersection;
							int ver_index = pIter->last_label;
							for (std::map<int, GICSimplicialTreeNode_ptr>::iterator 
								mIter = pIter->parent_ptr->children_map_ptr->begin();
								mIter != pIter->parent_ptr->children_map_ptr->end();
								mIter++)
							{
								if (mIter->first > pIter->last_label)//  the neighbors of pIter with index > pIter->last_label
									//prevIntersection.insert(mIter->first);
									prevIntersection.push_back(mIter->first);
							}
							//
							unsigned int pot_size = std::min(prevIntersection.size(), EuclideanDataPtr->vecNode[ver_index].edgeWeights.size());
							//
							std::set<int> curAdjList;
							for (std::map<int, float>::iterator mIter = EuclideanDataPtr->vecNode[ver_index].edgeWeights.begin();
								mIter != EuclideanDataPtr->vecNode[ver_index].edgeWeights.end(); 
								mIter++)
							{
								curAdjList.insert(mIter->first);
							}
							//
							std::vector<int> curIntersection(pot_size);
							std::vector<int>::iterator interIter = std::set_intersection(prevIntersection.begin(),
																			prevIntersection.end(),
																			curAdjList.begin(),
																			curAdjList.end(),
																			curIntersection.begin());
							if (interIter - curIntersection.begin() > 0)
							{// create new simplicies 
								Map_int_stnPtr_ptr pMap(new std::map<int, GICSimplicialTreeNode_ptr>);
								for (int sid = 0; sid < interIter - curIntersection.begin(); sid++)
								{
									//
									cur_dim_simplex_size++;// increase the # of edges
									//
									GICSimplicialTreeNode_ptr pStn(new GICSimplicialTreeNode(curIntersection[sid], // v_index
																					simplex_size++, // order in the filtration
																					pIter, // parent
																					Map_int_stnPtr_ptr(),
																					GICSimplicialTreeNode_ptr()
																					)
																					);
									// link into the circular list of *sIter
									InsertCircularList(head_circular_list_in_each_dim[curIntersection[sid]][depth - 1], pStn);
									//
									// push it into the map
									(*pMap.get())[curIntersection[sid]] = pStn;
								}
								// set up the children set
								pIter->children_map_ptr = pMap;
							}
							// move ahead
							pIter = pIter->next_circular_ptr;
						}while (pIter != head_circular_list_in_each_dim[vid][depth - 2]);
					}
				}
				Simplicies_Cnt[depth] = cur_dim_simplex_size;
			}//
		}
	}
	return true;
}