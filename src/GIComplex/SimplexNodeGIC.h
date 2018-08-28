/*
(c) 2012 Fengtao Fan
*/
#ifndef _SIMPLICIAL_TREE_NODE_H_
#define _SIMPLICIAL_TREE_NODE_H_

/* boost shared pointer */
#include <boost/shared_ptr.hpp>

#include <map>

/* declaration of simpicial tree node*/

class GICSimplicialTreeNode;

/* shared pointer to simplicial tree node */

typedef boost::shared_ptr<GICSimplicialTreeNode> GICSimplicialTreeNode_ptr;

/* sibling structure for each node */

typedef std::map<int, GICSimplicialTreeNode_ptr> Map_int_stnPtr;

/* shared pointer to <int, stn_ptr> mapping*/

typedef boost::shared_ptr<Map_int_stnPtr> Map_int_stnPtr_ptr;

class GICSimplicialTreeNode
{
public:
	// constructor
	GICSimplicialTreeNode() :
		last_label(-1),
		index_in_filtration(-1),
		parent_ptr(),
		children_map_ptr(),
		next_circular_ptr() //, prev_circular_ptr()
	{
	}

	GICSimplicialTreeNode(const int v_index,
						const int idx_filtration,
		GICSimplicialTreeNode_ptr  in_parent_ptr, //& changed for ubuntu
		Map_int_stnPtr_ptr  in_children_map_ptr, //&
		GICSimplicialTreeNode_ptr in_next_circular_ptr ) : //&
				last_label(v_index),
				index_in_filtration(idx_filtration),
				parent_ptr(in_parent_ptr),
				children_map_ptr(in_children_map_ptr),
				next_circular_ptr(in_next_circular_ptr)
	{
	}

public:
	int last_label;
	int index_in_filtration;
	GICSimplicialTreeNode_ptr parent_ptr;
	Map_int_stnPtr_ptr children_map_ptr;
	GICSimplicialTreeNode_ptr next_circular_ptr;
};

#endif // _SIMPLICIAL_TREE_NODE_H_
