/*
(c) 2012 Fengtao Fan
*/
#ifndef _SIMPLE_GRAPH_H_
#define _SIMPLE_GRAPH_H_ 
/*******
Edges are recorded in the adjacent list only ONCE !
********/
#include <set>
#include <map>
#include <vector>
#include <iterator>

using std::vector;

class SimpleGraphNode 
{ 
public:
	SimpleGraphNode() : v_index(-1), color(0)
	{
	}
	SimpleGraphNode(const int index) : v_index(index), color(0)
	{
	}
	SimpleGraphNode(const int index, const int cc) : v_index(index), color(cc)
	{
	}
	// copy constructor
	SimpleGraphNode(const SimpleGraphNode& rhs)
	{
		v_index = rhs.v_index; 
		color = rhs.color;

 		if (!rhs.edgeWeights.empty())
		{
			std::copy(	rhs.edgeWeights.begin(),
						rhs.edgeWeights.end(),
						std::inserter(edgeWeights, edgeWeights.begin()));
		}
 	}
	// assignment constructor
	SimpleGraphNode& operator=(const SimpleGraphNode& rhs)
	{
		//adjList.clear();
		edgeWeights.clear();
		//
		v_index = rhs.v_index; 
		color = rhs.color;

	 
		if (!rhs.edgeWeights.empty())
		{
			std::copy(	rhs.edgeWeights.begin(),
						rhs.edgeWeights.end(),
						std::inserter(edgeWeights, edgeWeights.begin()));
		}
		return *this;
	}
	~SimpleGraphNode()
	{ 
		//
		edgeWeights.clear();
	}
/************************************/
	inline bool is_neighbor(const int u)  const 
	{ 
		std::map<int, float>::const_iterator mFindIter = edgeWeights.find(u);
		return (mFindIter != edgeWeights.end());
	}
	inline int degree()
	{
		return int(edgeWeights.size());
	}
/************************************/
public:
	int v_index; 
	int color; 
	std::map<int, float> edgeWeights;
};

class SimpleGraph
{
	public:
	// constructor
	SimpleGraph() : color_number(0)
	{
	}
	SimpleGraph(const int n) : color_number(0)
	{
		SimpleGraphNode tmpNode;
		vecNode.reserve(n);
		for (int i = 0; i < n; i++)
		{
			tmpNode.v_index = i;
			vecNode.push_back(tmpNode);
		}
	}
	SimpleGraph(const SimpleGraph& rhs)
	{
		color_number = rhs.color_number;
		vecNode.assign(rhs.vecNode.begin(),
						rhs.vecNode.end());
	}
	SimpleGraph& operator=(const SimpleGraph& rhs)
	{
		color_number = rhs.color_number;
		vecNode.clear();
		vecNode.assign(rhs.vecNode.begin(),
						rhs.vecNode.end());
		return *this;
	}
	~SimpleGraph()
	{
		vecNode.clear();
	}
	/********************************/
	inline int ColorNumber ()
	{
		return color_number;
	}
	inline void SetColorNumber(const int c_num)
	{
		color_number = c_num;
		return;
	}
	/********************************/
	void AddNode(const int v_index)
	{
		SimpleGraphNode tmpNode(v_index);
		vecNode.push_back(tmpNode);
	}
	void AddNode(const int v_index, const int color)
	{
		SimpleGraphNode tmpNode(v_index, color);
		vecNode.push_back(tmpNode);
	} 
	// 
    void AddEdge(const int i, const int j)
    {// each edge has only one copy
        // stored into the adjList of smaller vertex index
        if (i < j)
        {
            vecNode[i].edgeWeights[j] = 0.0f;
        }
        else
        {
            // i > j
            vecNode[i].edgeWeights[i] = 0.0f;
        }
    }
    
    void AddEdge(const int i, const int j, const float w)
	{// each edge has only one copy
		// stored into the adjList of smaller vertex index
		if (i < j)
		{ 
			vecNode[i].edgeWeights[j] = w;
		}
		else
		{
			// i > j  
			vecNode[i].edgeWeights[i] = w;
		}
	}
	//
	bool is_edge(const int i, const int j)
	{
		bool ret = false;
		if (i < j)
		{
			ret = vecNode[i].is_neighbor(j);
		}
		else
		{// j < i
			ret = vecNode[j].is_neighbor(i);
		}
		return ret;
	}
	float edge_weight(const int i, const int j)
	{
		float ret = 0.f;
		if (i < j)
		{
			ret = vecNode[i].edgeWeights[j];
		}
		else
		{
			ret = vecNode[j].edgeWeights[i];
		}
		return ret;
	}
	//
	void InitNodes(const int n)
	{
		if (!vecNode.empty())
		{
			vecNode.clear();
			vecNode.reserve(n);
		}
		for (int v_index = 0; v_index < n; v_index++)
		{
			SimpleGraphNode tmpNode(v_index);
			vecNode.push_back(tmpNode);
		}
	}
	//
	void ReadWeightedGaph(char const * pFileName, std::vector<std::vector<float> > &pts);
	void WriteWeightedGraph(char const * pFileName, std::vector<std::vector<float> > &pts);
	//
	void ReadColorMappingFromFile(char const * pFileName);
	void WriteColorMappingToFile(char const * pFileName);
	//
	void ReadFromFile(char const* pFileName);
	void WriteBackToFile(char const * pFileName);
	long int EdgeNum();
	//
	int CheckComponents();
    
    void SetCoords(vector<vector<float>> coords);
    
	//
	public:
		std::vector<SimpleGraphNode> vecNode;
		int color_number;
	};
#endif //_SIMPLE_GRAPH_H_