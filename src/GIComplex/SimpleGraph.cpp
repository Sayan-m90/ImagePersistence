/*
(c) 2012 Fengtao Fan
*/
#include "SimpleGraph.h"

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

#include <map>
#include <vector>
#include <queue>

#include <cstdlib>

int SimpleGraph::CheckComponents()
{// check the # of components

	//
	std::set<int> unvisited_nodes;
	std::vector<bool> colored(vecNode.size());
	for (unsigned int i = 0; i < vecNode.size(); i++)
	{
		unvisited_nodes.insert(i);
		colored[i] = false;
	}
	//
	std::vector<std::set<int> > vecAdjList(vecNode.size());
	//
	for (unsigned int i = 0; i < vecNode.size(); i++)
	{
		for (std::map<int, float>::iterator mIter  = vecNode[i].edgeWeights.begin();
				mIter != vecNode[i].edgeWeights.end();
				mIter++)
		{
			vecAdjList[i].insert(mIter->first);
			vecAdjList[mIter->first].insert(i);
		}
	}
	//
	int component_size = 0;
	//std::cout << "comp " << std::endl;
	while (!unvisited_nodes.empty())
	{
		component_size++;
		int src_node = *unvisited_nodes.begin();
		//
		// DFS visit
		std::queue<int> nodeQ;
		nodeQ.push(src_node);
		//

		while (!nodeQ.empty())
		{
			int cur_node = nodeQ.front();
			//
			colored[cur_node] = true;
			//
			nodeQ.pop();
			//
			unvisited_nodes.erase(cur_node);
			//
			for (std::set<int>::iterator sIter  = vecAdjList[cur_node].begin();
				sIter != vecAdjList[cur_node].end();
				sIter++)
			{
				if (colored[*sIter] == false)
				{
					nodeQ.push(*sIter);
					colored[*sIter] = true;
				}
			}
		}
	}
	return component_size;
}
//
void SimpleGraph::WriteWeightedGraph(char const *pFileName, std::vector<std::vector<float> > &pts)
{
	std::cout << "Writing < " << pFileName << " > " << std::endl;
	//
	std::ofstream ofile;
	ofile.open(pFileName, std::ifstream::out);
  	//
 	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		if (!pts.empty())
		{
			sstr << pts.front().size() << " ";
		}
		else
		{
			sstr << "0 " ;
		}
		sstr << vecNode.size() << std::endl;
		//
		// if necessary, write the coordinates
		if (!pts.empty())
		{
			for (unsigned int i = 0; i < pts.size(); i++)
			{
				for (unsigned int j = 0; j < pts[i].size(); j++)
				{
					sstr << pts[i][j] << " " ;
				}
				sstr << std::endl;
			}
		}
		for (unsigned int i = 0; i < vecNode.size(); i++)
		{
			if (!vecNode[i].edgeWeights.empty())
			{
				for (std::map<int, float>::iterator sIter = vecNode[i].edgeWeights.begin();
					sIter != vecNode[i].edgeWeights.end();
					sIter++)
				{
					sstr << i << " " << sIter->first << " " << sIter->second << std::endl;
				}
			}
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
	std::cout << "---- Done ---- "  << std::endl;
	//
	return;
}
//
void SimpleGraph::ReadWeightedGaph(char const * pFileName, std::vector<std::vector<float> > &pts)
{
	std::cout << "Reading < " << pFileName << " > " << std::endl;

	std::ifstream ifile;
	ifile.open(pFileName, std::ifstream::in);

	std::string sBuf;
	//
	int nPointDimension = 0;
	int vertexNum = 0;
	std::vector<float> coords;
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
			if (!sBuf.empty() &&  ( sBuf[0] >= '0' && sBuf[0] <= '9'  || sBuf[0] == '-'))
			{
				if (bFirstRead)
				{
					locSstr.str(sBuf);
					locSstr >> nPointDimension;
					locSstr >> vertexNum ;
					//resize the vector
					//vecNode.resize(vertexNum);
					if (!vecNode.empty())
					{
						vecNode.clear();
					}
					InitNodes(vertexNum);
					//
					if (nPointDimension && vertexNum > 0)
					{
						coords.resize(nPointDimension);
						pts.reserve(vertexNum);
					}
					//
					locSstr.clear();
					bFirstRead = false;
				}
				else
				{//
					if (nPointDimension && vertexNum > 0)
					{
						// read point coordinates
						locSstr.str(sBuf);
						for (int i = 0; i < nPointDimension; i++)
						{
							locSstr >> coords[i];
						}
						locSstr.clear();
						//
						pts.push_back(coords);
						//
						vertexNum--;
					}
					else
					{// read weighted edges
						int a = 0;
						int b = 0;
						float w = 0.f;
						locSstr.str(sBuf);
						locSstr >> a;
						locSstr >> b;
						locSstr >> w;
						//
						AddEdge(a, b, w);
						//
						locSstr.clear();
					}
				}
			}
			//sBuf.clear();
		}
		sstr.clear();
	}
	else
	{
		std::cout << "Can NOT open file " << pFileName << std::endl;
		exit(0);
	}
	std::cout << "---- Done ----" << std::endl;
	//
	return;
}
//
void SimpleGraph::ReadColorMappingFromFile(char const * pFileName)
{
	std::cout << " ... Reading Color Mapping ... " << std::endl;

	std::ifstream ifile;
	ifile.open(pFileName, std::ifstream::in);

	std::string sBuf;
	//
	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
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
		 /* format
		 2 // number of colors
		 0 0 1 1 1 // color mapped for 5 vertices
		 */
		//
		sstr >> color_number;
		for (unsigned int i = 0; i < vecNode.size(); i++)
			sstr >> vecNode[i].color;
		//
		sstr.clear();
	}
	else
	{
		std::cout << "Can NOT open file " << pFileName << std::endl;
		exit(0);
	}
	std::cout << "...Done... " << color_number << std::endl;
	return;
}
void SimpleGraph::WriteColorMappingToFile(char const * pFileName)
{
// N number of distinct color
// color_index_v0 color_index_v1 color_index_v2
//
 /* format
 2 // number of colors
 0 0 1 1 1 // color mapped for 5 vertices
 */
	std::cout << " ... Wrting Color Mapping... " << std::endl;
	//
	std::ofstream ofile;
	ofile.open(pFileName, std::ifstream::out);
  	//
 	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		sstr << color_number << std::endl; // first line is # of distinct colors
		//sstr << vecNode.size() << std::endl;
		for (unsigned int i = 0; i < vecNode.size(); i++)
		{
			sstr << vecNode[i].color << " ";
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
	std::cout << "Done... " << color_number << std::endl;
	//
	return;
}
void SimpleGraph::WriteBackToFile(char const* pFileName)
{
	std::cout << "wrting... " << std::endl;
	//
	std::ofstream ofile;
	ofile.open(pFileName, std::ifstream::out);
  	//
 	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		sstr << vecNode.size() << std::endl;
		for (unsigned int i = 0; i < vecNode.size(); i++)
		{
			if (!vecNode[i].edgeWeights.empty())
			{
				for (std::map<int, float>::iterator mIter = vecNode[i].edgeWeights.begin();
					mIter != vecNode[i].edgeWeights.end();
					mIter++)
				{
					sstr << 2 << " " << i << " " << mIter->first << std::endl;
				}
			}
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
	std::cout << "Done... " << vecNode.size() << std::endl;
	//
	return;
}
long int SimpleGraph::EdgeNum()
{// each edge is stored only once
	long int edgeNum = 0;
	for (unsigned  int i = 0; i < vecNode.size(); i++)
		edgeNum += vecNode[i].edgeWeights.size();//long int(vecNode[i].edgeWeights.size());
	return edgeNum;
}
void SimpleGraph::ReadFromFile(char const * pFileName)
{
	std::cout << "reading... " << std::endl;

	std::ifstream ifile;
	ifile.open(pFileName, std::ifstream::in);

	std::string sBuf;
	//
	int vertexNum = 0;
	//
	bool bFirstRead = true;
	int vertexCount = 0;
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
			if (!sBuf.empty() &&  ( sBuf[0] >= '0' && sBuf[0] <= '9' ))
			{
				if (bFirstRead)
				{
					locSstr.str(sBuf);
					locSstr >> vertexNum ;
					//resize the vector
					//vecNode.resize(vertexNum);
					vecNode.clear();
					InitNodes(vertexNum);
					//
					locSstr.clear();
					bFirstRead = false;
				}
				else
				{//

					int adjListSize = 0;
					int adjElem = 0;
					locSstr.str(sBuf);
					locSstr >> adjListSize;
					for (int k = 0; k < adjListSize; k++)
					{
						locSstr >> adjElem;
						if (adjElem > vertexCount) // make it work also on the adjacent list structure of graph
													// where (i, j) edge makes i appear in j-list but also j appear in i-list
						{
							AddEdge(vertexCount, adjElem);
						}
					}
					locSstr.clear();
					//
					//tmpGraphNode.adjList.clear();
					vertexCount++;
					if (vertexCount == vertexNum)
						break;
				}
			}
			//sBuf.clear();
		}
		sstr.clear();
	}
	else
	{
		std::cout << "Can NOT open file " << pFileName << std::endl;
		exit(0);
	}
	std::cout << "Done... " << vecNode.size() << std::endl;
	//
	return;
}

void SimpleGraph::SetCoords(vector<vector<float>> new_coords) {
    
}
