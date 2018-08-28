/*
(c) 2012 Fengtao Fan
*/
#include "GIComplex.h"
#include "FlagComplex.h"

#include <set>
#include <map>
#include <algorithm>

bool GIComplex::Construction()
{
	// get the flag complex from the graph
	FlagComplex rips_complex(dim, EuclideanDataPtr->vecNode.size(), EuclideanDataPtr);
	rips_complex.Construction();
	// index the simplicies by the dimension order
	int simplex_size = Simplicies_Cnt[0];
	// vertex has constructed in the initialization
	// here each vertex in the graph mapped to some point (its color) in the sub-sampling
	// Construct the edge
	if (dim > 0)
	{
		int cur_dim_simplex_size = 0;
		for (int d = 1; d < dim + 1; d++)
		{
			cur_dim_simplex_size = 0;
			//
			if (Simplicies_Cnt[d - 1] == 0)
			{// no one lower dim simlex, means no current dim simplex
				Simplicies_Cnt[d] = 0;
			}
			else
			{
				for (int vid = 0; vid < rips_complex.Simplicies_Cnt[0]; vid++)
				{
					if (rips_complex.head_circular_list_in_each_dim[vid][d - 1])
					{
						GICSimplicialTreeNode_ptr pIter(rips_complex.head_circular_list_in_each_dim[vid][d - 1]);
						do
						{// visit each simplex
							// get the colors of this simplex
							std::set<int> simplex_colors;
							GICSimplicialTreeNode_ptr pSimplexIter(pIter);
							int simplex_dim = 0;
							do
							{
								simplex_dim++;
								simplex_colors.insert(EuclideanDataPtr->vecNode[pSimplexIter->last_label].color);
								pSimplexIter = pSimplexIter->parent_ptr;
							}while (pSimplexIter);
							//
							if (simplex_colors.size() == simplex_dim)
							{
								
							// check if this simplex exists or not
							bool simplex_existence = true;
							//  simplex_color is not empty here and with size >= 2, otherwise error
							std::set<int>::iterator sIter = simplex_colors.begin();
							pSimplexIter = vertex_array[*simplex_colors.begin()];
							//
							for (	sIter++;
									sIter != simplex_colors.end();
									sIter++)
							{
								if (pSimplexIter->children_map_ptr)
								{
									std::map<int, GICSimplicialTreeNode_ptr>::iterator mFindIter = pSimplexIter->children_map_ptr->find(*sIter);
									if ( mFindIter == pSimplexIter->children_map_ptr->end())
									{
										simplex_existence = false;
										break;
									}
									pSimplexIter = mFindIter->second;
								}
								else
								{// empty std::map
									simplex_existence = false;
									break;
								}
							}
							if (!simplex_existence)
							{// create new one
								int ver_index = *simplex_colors.rbegin();
								cur_dim_simplex_size++;// increase the # of edges
								//
								GICSimplicialTreeNode_ptr pStn(new GICSimplicialTreeNode(ver_index, // v_index
																				simplex_size++, // order in the filtration
																				pSimplexIter, // parent
																				Map_int_stnPtr_ptr(),
																				GICSimplicialTreeNode_ptr()
																				)
																				);
								// link into the circular list of *sIter
								InsertCircularList(head_circular_list_in_each_dim[ver_index][d - 1], pStn);

								// push it into the map
								if (pSimplexIter->children_map_ptr)
								{
									(*(pSimplexIter->children_map_ptr.get()))[ver_index] = pStn;
								}
								else
								{
									Map_int_stnPtr_ptr pMap(new std::map<int, GICSimplicialTreeNode_ptr>);
									(*pMap.get())[ver_index] = pStn;
									pSimplexIter->children_map_ptr = pMap;
								}
							}
							}
							// move ahead
							pIter = pIter->next_circular_ptr;
						}while (pIter != rips_complex.head_circular_list_in_each_dim[vid][d - 1]);
					}// if 
				}// for vid
				Simplicies_Cnt[d] = cur_dim_simplex_size;
			}// else			
		}// for d
	}// if dim

	return true;
}
