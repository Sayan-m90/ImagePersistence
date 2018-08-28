/*
(c) 2012 Fengtao Fan
*/
#ifndef _SIMPLICIAL_TREE_H_
#define _SIMPLICIAL_TREE_H_

#include "SimplexNodeGIC.h"
//#include "filtration.h"

#include <vector>
#include <set>

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <cstring>

/* declaration of simpicial tree*/
template<typename T>
class SimplicialTree
{
public:
	/*constructor */
	// default constructor
	SimplicialTree() : dim(-1), Simplicies_Cnt(NULL), EuclideanDataPtr(NULL)
	{
	}
    
	// constructor with vertex number and dim
    SimplicialTree(const int in_dim, const int v_number, T *inDataPtr ) : dim(in_dim) , EuclideanDataPtr(inDataPtr)
	{// dim >= 0
		Simplicies_Cnt = new int[dim + 1];
		Simplicies_Cnt[0] = v_number; // v_number >= 0
		/*
		intialize the zero dim (vertex) data
		*/
		vertex_array.resize(Simplicies_Cnt[0]);
		for (int i = 0; i < Simplicies_Cnt[0]; i++)
		{
			GICSimplicialTreeNode_ptr p(new GICSimplicialTreeNode(i,
															i,
															GICSimplicialTreeNode_ptr(),
															Map_int_stnPtr_ptr(),
															GICSimplicialTreeNode_ptr()
															)
															);
			vertex_array[i] = p;
		}
		/*
		set up the head_circular_list
		*/
		head_circular_list_in_each_dim.resize(Simplicies_Cnt[0]);
		//
		for (int i = 0; i < Simplicies_Cnt[0]; i++)
		{
			for (int j = 0; j < dim; j++) // dim zero is not included
				head_circular_list_in_each_dim[i].push_back(GICSimplicialTreeNode_ptr());
		}
	}

	/*deconstructor */
	~SimplicialTree()
	{
		if (Simplicies_Cnt)
		{
			delete [] Simplicies_Cnt;
		}
	}
	// return the size of the complex
	inline int ComplexSize()
	{
		int total_size = 0;
		for (int i = 0; i < dim + 1; i++)
		{
			total_size += Simplicies_Cnt[i];
		}
		return total_size;
	}
	//
	inline int dDimSimplicesSize(const int d)
	{
		if (d > dim || d < 0)
		{
			std::cout << "Querying dimension exceeding the Complex Dimension " << std::endl;
			return  -1;
		}
		return Simplicies_Cnt[d];
	}
	inline void Report()
	{
		for (int i = 0; i < dim + 1; i++)
		{
			std::cout << "DIM " << i << " : " << Simplicies_Cnt[i] << std::endl;
		}
		std::cout << "total size : " << ComplexSize() << std::endl;
		return;
	}
	//
	bool InsertCircularList(GICSimplicialTreeNode_ptr & head, const GICSimplicialTreeNode_ptr & elem) ;
//{
//	if (head) // check the shared ptr is null
//	{// not null
//		// insert at head
//		elem->next_circular_ptr = head;
//		elem->prev_circular_ptr = head->prev_circular_ptr;
//		head->prev_circular_ptr->next_circular_ptr = elem;
//		head->prev_circular_ptr = elem;
//		//
//		head = elem;
//	}
//	else
//	{// it is null
//		head = elem;
//		elem->next_circular_ptr = elem;
//		elem->prev_circular_ptr = elem;
//	}
//	return true;
//}
	// retrieve boundary
	bool Boundary(const GICSimplicialTreeNode_ptr& sigma, std::vector<GICSimplicialTreeNode_ptr> &bdries);
	// retrieve coboundary
	bool CoBoundary(const GICSimplicialTreeNode_ptr sigma, std::vector<GICSimplicialTreeNode_ptr> coBdries);
	// construct the complex
	virtual bool Construction() = 0;
	// build a filtration
	//bool BuildFiltration(filtration &out_fil);
	// make a filtration where smaller dimension simplex precedes before large dimension simplex
	// bool MakeFiltrationOrderedByDimension(const filtration &out_fil);
	void WriteBackToFile(const char* pFileName);
	void ReadFromFile(const char* pFileName);
	//
	void WriteStatisticsToFile(const char* pFileName);
	//
public:
	T* EuclideanDataPtr; // data used for constructing this complex, like the underlying graph for Rips complex
	int dim; // dimension of the simpicial complex
	int *Simplicies_Cnt; // # of simplicies in each dimension
	//
	std::vector<GICSimplicialTreeNode_ptr> vertex_array; // dimension zero simplicies
	//
	std::vector<std::vector<GICSimplicialTreeNode_ptr> > head_circular_list_in_each_dim;
};
template <class T>
bool SimplicialTree<T>::Boundary(const GICSimplicialTreeNode_ptr& sigma, std::vector<GICSimplicialTreeNode_ptr> &bdries)
{
	// traverse back to the root to get all boundary faces
	GICSimplicialTreeNode_ptr pIter(sigma);
	GICSimplicialTreeNode_ptr pBdryFace;
	// vertex simplex doesn't have boundary
	// the half vertex set
	std::set<int> tail_ver_index_set;
	if (pIter->parent_ptr)
	{// simplex of dim > 0
		do
		{// the boundary face without current vertex
			//if (pIter == sigma)
			//{// the boundary face without the last vertex
			//	bdries.push_back(pIter->parent_ptr);
			//	tail_ver_index_set.insert(pIter->last_label); // record the tail vertex for next boundary face
			//}
			//else
			// go from its parent to visit all tail vertex in tail_ver_index_set
			// at the end of visit, the boundary face pointer is reached.
			if (pIter->parent_ptr)
			{
				pBdryFace = pIter->parent_ptr;
			}
			else
			{// at the moment of removing the first vertex
				// require tail_ver_index_set NOT empty
				pBdryFace = vertex_array[*tail_ver_index_set.begin()];
				tail_ver_index_set.erase(tail_ver_index_set.begin());
			}
			//
			for (std::set<int>::iterator sIter = tail_ver_index_set.begin();
										sIter != tail_ver_index_set.end();
										sIter++)
			{//
				pBdryFace = (*pBdryFace->children_map_ptr.get())[*sIter];
			}
			//
			bdries.push_back(pBdryFace);
			//
			tail_ver_index_set.insert(pIter->last_label);
			// traverse back one more step
			pIter = pIter->parent_ptr;
		}while (pIter);
	}
	/******************************************/
	return true;
}
template <class T>
bool SimplicialTree<T>::InsertCircularList(GICSimplicialTreeNode_ptr & head, const GICSimplicialTreeNode_ptr & elem)
{
	if (head) // check the shared ptr is null
	{// not null
		//// insert at head
		//elem->next_circular_ptr = head;
		//elem->prev_circular_ptr = head->prev_circular_ptr;
		//head->prev_circular_ptr->next_circular_ptr = elem;
		//head->prev_circular_ptr = elem;
		//
		//head = elem;

		// insert just behind the head element
		elem->next_circular_ptr = head->next_circular_ptr;
		head->next_circular_ptr = elem;
	}
	else
	{// it is null
		//head = elem;
		//elem->next_circular_ptr = elem;
		//elem->prev_circular_ptr = elem;
		head = elem;
		elem->next_circular_ptr = elem;
	}
	return true;
}
//template <class T>
//bool SimplicialTree<T>::BuildFiltration(filtration &out_fil)
//{ // Preq: all simplicies are sorted in some order
//
//	std::vector<std::list<int> *> bdMat;
//	bdMat.resize(ComplexSize());
//
//	// empty boundary matrix for each simplex
//	for (unsigned int i = 0; i < bdMat.size(); i++)
//		bdMat[i] = new std::list<int>;
//
//	// vertices are trivially added as they have empty boundaries
//	// now add edges
//	for (int i = 1; i < dim + 1; i++)
//	{
//		//depth = i + 1; // want to visit all simplices at this dimension
//		if (Simplicies_Cnt[i] != 0)
//		{// visit each simplex in dimension i
//			for (int vid = 0; vid < Simplicies_Cnt[0]; vid++)
//			{
//				if (head_circular_list_in_each_dim[vid][i - 1]) // i-1 as vertex are stored in an array
//				{// the circular list is not empty
//					GICSimplicialTreeNode_ptr pIter(head_circular_list_in_each_dim[vid][i - 1]);
//					do
//					{// visit each simplex
//						std::vector<GICSimplicialTreeNode_ptr> bdries;
//						//
//						int curSimplexIndex = pIter->index_in_filtration;
//						Boundary(pIter, bdries);
//						// get the index
//						std::set<int> bdryIndices;
//						for (std::vector<GICSimplicialTreeNode_ptr>::iterator vIter = bdries.begin();
//								vIter != bdries.end();
//								vIter++)
//						{
//							bdryIndices.insert(vIter->get()->index_in_filtration);
//						}
//						// fill the boundary matrix
//						bdMat[curSimplexIndex]->resize(bdryIndices.size());
//						std::copy(bdryIndices.begin(), bdryIndices.end(), bdMat[curSimplexIndex]->begin());
//						// move ahead
//						pIter = pIter->next_circular_ptr;
//					}while (pIter != head_circular_list_in_each_dim[vid][i - 1]);
//				} // if head
//			}// for vid
//		} // if si
//	}// for i dim
//	out_fil.assign(bdMat, dim + 1);
//	// clean the memory
//	for (unsigned int i = 0; i < bdMat.size(); i++)
//	{
//		std::list<int> tmp;
//		bdMat[i]->swap(tmp);
//		delete bdMat[i];
//	}
//	//
//	return true;
//}
template <class T>
void SimplicialTree<T>::ReadFromFile(const char* pFileName)
{
// MAX_DIM VER_NUM #max dimension of simplex (0-dim vertex 1-dim edge 2-dim triangle) #vertex number
// simplex_dim ver_0 ver_1 ... ver_dim
//
 /* format
 2 5 # max_dim 2, 5 vertices
 1 0 1
 1 2 3
 .
 .
 .
 1 1 4
 */
	std::cout << " ...Reading Files ... " << std::endl;
	//
	std::ifstream ifile;
	ifile.open(pFileName, std::ifstream::in);

	//
	std::string sBuf;
	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	//
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

		// dimension and # of vertices are always there
		// read the dimension
		sstr >> dim ;
		Simplicies_Cnt = new int[dim + 1];
		memset(Simplicies_Cnt, 0, sizeof(int) * (dim + 1));
		// read the number of vertices
		sstr >> Simplicies_Cnt[0];
		//Simplicies_Cnt[0] = v_number; // v_number >= 0
		/*
		intialize the zero dim (vertex) data
		*/
		vertex_array.resize(Simplicies_Cnt[0]);
		for (int i = 0; i < Simplicies_Cnt[0]; i++)
		{
			GICSimplicialTreeNode_ptr p(new GICSimplicialTreeNode(i,
															i,
															GICSimplicialTreeNode_ptr(),
															Map_int_stnPtr_ptr(),
															GICSimplicialTreeNode_ptr()
															)
															);
			vertex_array[i] = p;
		}

		/*
		set up the head_circular_list
		*/
		head_circular_list_in_each_dim.resize(Simplicies_Cnt[0]);
		//
		for (int i = 0; i < Simplicies_Cnt[0]; i++)
		{
			for (int j = 0; j < dim; j++) // dim zero is not included
				head_circular_list_in_each_dim[i].push_back(GICSimplicialTreeNode_ptr());
		}
		// now add edges and other simplices
		int simplex_counter = Simplicies_Cnt[0];
		// read the simplex
		int simplex_dim = 1;
		std::vector<int> simplex_vertex(simplex_dim);
		while (sstr.good())
		{
			sstr >> simplex_dim ;
			if (!sstr.good())
				break;
			if (simplex_dim + 1 > simplex_vertex.size())
				simplex_vertex.resize(simplex_dim + 1);
			for (int i = 0; i < simplex_dim + 1 ; i++)
				sstr >> simplex_vertex[simplex_dim - i];// this array is stored by inverse order
			if (!sstr.good())
				break;
			// get the parent of current simplex
			GICSimplicialTreeNode_ptr curParent = vertex_array[simplex_vertex[0]];
			for (int i = 1; i < simplex_dim; i++)
			{
				curParent = (*curParent->children_map_ptr)[simplex_vertex[i]];
			}
			// construct a simplex
			Simplicies_Cnt[simplex_dim]++;
			GICSimplicialTreeNode_ptr pSimp(new GICSimplicialTreeNode(simplex_vertex[simplex_dim],
															simplex_counter++,
															curParent,
															Map_int_stnPtr_ptr(),
															GICSimplicialTreeNode_ptr()
															)
															);
			// link it into the simplicial tree
			if (curParent->children_map_ptr)
			{
				(*curParent->children_map_ptr)[simplex_vertex[simplex_dim]] = pSimp;
			}
			else
			{
				Map_int_stnPtr_ptr pMap(new std::map<int, GICSimplicialTreeNode_ptr>);
				(*pMap.get())[simplex_vertex[simplex_dim]] = pSimp;
				curParent->children_map_ptr = pMap;
			}
			//
			InsertCircularList(head_circular_list_in_each_dim[simplex_vertex[simplex_dim]][simplex_dim - 1], pSimp);
			//// insert all simplices with dim >= 1 into the simplicial complex

		}// while

		ifile.close();
		sstr.clear();
	}
	else
	{
		std::cout << "Can NOT open file " << pFileName << std::endl;
		exit(0);
	}
	std::cout << "Done... " << pFileName << std::endl;
	//
	return;
}

template <class T>
void SimplicialTree<T>::WriteStatisticsToFile(const char* pFileName)
{
	std::cout << "Writing statistics < " << pFileName << " >" << std::endl;
	std::ofstream ofile;
	ofile.open(pFileName, std::ifstream::out);
  	//
 	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		//sstr << pFileName << std::endl;

		for (int i = 0; i < dim + 1; i++)
		{
			sstr << "DIM " << i << " : " << Simplicies_Cnt[i] << std::endl;
		}
		sstr << "Total size : " << ComplexSize() << std::endl;
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
	std::cout << "---- Done ----" << std::endl << std::endl;
	return ;
}
template <class T>
void SimplicialTree<T>::WriteBackToFile(const char* pFileName)
{
// MAX_DIM VER_NUM #max dimension of simplex (0-dim vertex 1-dim edge 2-dim triangle) #vertex number
// simplex_dim ver_0 ver_1 ... ver_dim
//
 /* format
 2 5 # max_dim 2, 5 vertices
 1 0 1
 1 2 3
 .
 .
 .
 1 1 4
 */
	std::cout << " ... Wrting Color Mapping... " << std::endl;
	//
	std::ofstream ofile;
	ofile.open(pFileName, std::ifstream::out);
  	//
 	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		sstr << dim << " " << vertex_array.size() << std::endl; // first line is # of distinct colors
		//sstr << vecNode.size() << std::endl;
	// vertices are trivially added as they have empty boundaries
		// now add edges
		for (int i = 1; i < dim + 1; i++)
		{
			//depth = i + 1; // want to visit all simplices at this dimension
			if (Simplicies_Cnt[i] != 0)
			{// visit each simplex in dimension i
				for (int vid = 0; vid < Simplicies_Cnt[0]; vid++)
				{
					if (head_circular_list_in_each_dim[vid][i - 1]) // i-1 as vertex are stored in an array
					{// the circular list is not empty
						GICSimplicialTreeNode_ptr pIter(head_circular_list_in_each_dim[vid][i - 1]);
						do
						{// visit each simplex
							sstr << i << " ";
							GICSimplicialTreeNode_ptr pParentIter(pIter);
							//
							do
							{
								sstr << pParentIter->last_label << " ";
								//
								pParentIter = pParentIter->parent_ptr;
							}while (pParentIter);
							sstr << std::endl;
							// move ahead
							pIter = pIter->next_circular_ptr;
						}while (pIter != head_circular_list_in_each_dim[vid][i - 1]);
					} // if head
				}// for vid
			} // if si
		}// for i dim
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
	std::cout << "Done... " << pFileName << std::endl;
	//
	return;
}


class GIComplex : public SimplicialTree<SimpleGraph>
{
public:
    GIComplex()
    {
    };
    
    GIComplex(const int in_dim, const int v_number, SimpleGraph *inDataPtr) : SimplicialTree(in_dim, v_number, inDataPtr) {
    }
    bool Construction(); //virtual
};


#endif // _SIMPLICIAL_TREE_H_
