/*
(c) 2012 Fengtao Fan
*/
#include "ANNSearchSampling.h"
#include "SimpleGraph.h"
#include "PointSet.h"
#include "GIComplex.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <boost/program_options.hpp>

////////////////
void ConstructRipsGraph(const float eps_val, const PointSet &inPts,  SimpleGraph &outRipsGraph)
{
//	Construct Rips graph on the input points
	float sqDistance = eps_val * eps_val;
	ANNSearch::BuildDistanceGraph_ann(outRipsGraph, inPts, sqDistance);
	return;
}
void ComputeEdgeWeights(const PointSet &inPts, SimpleGraph &outRipsGraph)
{
	float segLength = 0.f;
	LongVector tempDiff;
	//
	for (unsigned int i = 0; i < outRipsGraph.vecNode.size(); i++)
	{
		for (std::map<int, float>::iterator mIter = outRipsGraph.vecNode[i].edgeWeights.begin();
				mIter != outRipsGraph.vecNode[i].edgeWeights.end();
				mIter++)
		{
			tempDiff = (*inPts._PointSet)[i] - (*inPts._PointSet)[mIter->first];
			segLength = norm(tempDiff);
			mIter->second = segLength;
		}
	}
	return;
}
void DeltaSparseSampling_EuclideanDistance(const PointSet &inPts,  const float delta_dist,
					std::map<int, int> &SubPointIndex,
					std::vector<int> &ColorMapping)
{
	//
	float sqDistance = 0.0;
	//
	// Perform the delta-sampling delta-sparse subsampling
	//
	sqDistance = delta_dist * delta_dist;
	ColorMapping.resize(inPts._PointSet->size());
	ANNSearch::Subsampling_EuclideanDistance(sqDistance,
											inPts,
											SubPointIndex,
											ColorMapping);

	return;
}
void DeltaSparseSampling_GraphDistance( const SimpleGraph &eps_rips_graph,
										const float delta_dist,
										std::map<int, int> &SubPointIndex,
										std::vector<int> &ColorMapping)
{
	//
	// Perform the delta-sampling delta-sparse subsampling
	//
	ColorMapping.resize(eps_rips_graph.vecNode.size());
	ANNSearch::Subsampling_GraphDistance(delta_dist,
										 eps_rips_graph,
										 SubPointIndex,
										 ColorMapping);
	return;
}
/////////////////

void ComputeBiColoredGraph(const float eps_sampling_dist, const float delta_dist, const bool graph_distance_flag, const PointSet &inPts,
						std::map<int, int> &SubPointIndex, SimpleGraph &biColoredGraph)
{
	SimpleGraph RipsGraph;
	std::vector<int> ColorMapping;

	ConstructRipsGraph(eps_sampling_dist, inPts, RipsGraph);

	//std::cout << "org comp " << RipsGraph.CheckComponents() << std::endl;
	// generate the batch file for creating rips using this rips graph
	// Delat-Sparse-Sampling and construct rips graph on sub-sampling points
	if (graph_distance_flag)
	{
		ComputeEdgeWeights(inPts, RipsGraph);
		DeltaSparseSampling_GraphDistance(RipsGraph, delta_dist, SubPointIndex, ColorMapping);
	}
	else
		DeltaSparseSampling_EuclideanDistance(inPts, delta_dist, SubPointIndex, ColorMapping);

	// construct the bi-colored graph from the org rips graph
	ANNSearch::SetColorMappingAndExtractColoredGraph(ColorMapping, RipsGraph, biColoredGraph);

	//std::cout << "bic comp " << biColoredGraph.CheckComponents() << std::endl;
	biColoredGraph.color_number = (int) SubPointIndex.size();
	return;
}

void ComputeBiColoredGraph(const char* pFileName, const float delta_dist, std::vector<std::vector<float> > &pts,
						std::map<int, int> &SubPointIndex, SimpleGraph &biColoredGraph)
{
	SimpleGraph orgGraph;
	std::vector<int> ColorMapping;

	orgGraph.ReadWeightedGaph(pFileName, pts);
	
	//std::cout << "org comp " << RipsGraph.CheckComponents() << std::endl;
	// generate the batch file for creating rips using this rips graph
	// Delat-Sparse-Sampling and construct rips graph on sub-sampling points

	DeltaSparseSampling_GraphDistance(orgGraph, delta_dist, SubPointIndex, ColorMapping);

	// construct the bi-colored graph from the org rips graph
	ANNSearch::SetColorMappingAndExtractColoredGraph(ColorMapping, orgGraph, biColoredGraph);

	//std::cout << "bic comp " << biColoredGraph.CheckComponents() << std::endl;
	biColoredGraph.color_number = (int)SubPointIndex.size();
	//RipsGraph.WriteBackToFile("test_rip.txt");
	return;
}

void WriteOFFformatComplex(const char* pFileName, const std::vector<std::vector<float> > &pts, std::map<int, int> &subPointIndices, GIComplex &gic)
{
	std::cout << "Writing < " << pFileName << " >" << std::endl;
	//
	std::ofstream ofile;
	ofile.open(pFileName, std::ifstream::out);
  	//
 	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		if (!pts.empty())
		{
			// write points out
			std::vector<int> vecSubPts(subPointIndices.size());
			for (std::map<int, int>::iterator mIter = subPointIndices.begin();
				mIter != subPointIndices.end(); mIter++)
			{
				vecSubPts[mIter->second] = mIter->first;
			}
			//
			sstr << pts.front().size() << " " << vecSubPts.size() << std::endl;
			for (unsigned int i = 0; i < vecSubPts.size(); i++)
			{
				for (unsigned v = 0; v < pts.front().size(); v++)
				{
					sstr << pts[vecSubPts[i]][v] << " ";
				}
				sstr << std::endl;
			}
		}
		else
		{
			sstr << "0 " << subPointIndices.size() << std::endl;
			std::vector<int> subVertexIndices(subPointIndices.size());
			for (std::map<int, int>::iterator mIter = subPointIndices.begin();
				mIter != subPointIndices.end(); mIter++)
				subVertexIndices[mIter->second] = mIter->first;
			for (unsigned int i = 0; i < subVertexIndices.size(); i++)
				sstr << subVertexIndices[i] << std::endl;
		}
		//
		// now add edges
		//
		for (int i = 1; i < gic.dim + 1; i++)
		{
			//depth = i + 1; // want to visit all simplices at this dimension
			if (gic.Simplicies_Cnt[i] != 0)
			{// visit each simplex in dimension i
				for (int vid = 0; vid < gic.Simplicies_Cnt[0]; vid++)
				{
					if (gic.head_circular_list_in_each_dim[vid][i - 1]) // i-1 as vertex are stored in an array
					{// the circular list is not empty
						GICSimplicialTreeNode_ptr pIter(gic.head_circular_list_in_each_dim[vid][i - 1]);
						do
						{// visit each simplex
							sstr << i + 1<< " ";
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
						}while (pIter != gic.head_circular_list_in_each_dim[vid][i - 1]);
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
	std::cout << "---- Done ---- " << std::endl << std::endl;
	//
}
void WriteOFFformatComplex(const char* pFileName, const PointSet &pts, std::map<int, int> &subPointIndices, GIComplex &gic)
{
	std::cout << "Writing < " << pFileName << " >" << std::endl;
	//
	std::ofstream ofile;
	ofile.open(pFileName, std::ifstream::out);
  	//
 	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		// write points out
		std::vector<int> vecSubPts(subPointIndices.size());
		for (std::map<int, int>::iterator mIter = subPointIndices.begin();
			mIter != subPointIndices.end(); mIter++)
		{
			vecSubPts[mIter->second] = mIter->first;
		}
		//
		sstr << pts._dimension << " " << vecSubPts.size() << std::endl;
		for (unsigned int i = 0; i < vecSubPts.size(); i++)
		{
			for (unsigned v = 0; v < pts._dimension; v++)
			{
				sstr << (*pts._PointSet)[vecSubPts[i]][v] << " ";
			}
			sstr << std::endl;
		}
		//
		// now add edges
		//
		for (int i = 1; i < gic.dim + 1; i++)
		{
			//depth = i + 1; // want to visit all simplices at this dimension
			if (gic.Simplicies_Cnt[i] != 0)
			{// visit each simplex in dimension i
				for (int vid = 0; vid < gic.Simplicies_Cnt[0]; vid++)
				{
					if (gic.head_circular_list_in_each_dim[vid][i - 1]) // i-1 as vertex are stored in an array
					{// the circular list is not empty
						GICSimplicialTreeNode_ptr pIter(gic.head_circular_list_in_each_dim[vid][i - 1]);
						do
						{// visit each simplex
							sstr << i + 1<< " ";
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
						}while (pIter != gic.head_circular_list_in_each_dim[vid][i - 1]);
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
	std::cout << "---- Done ---- " << std::endl << std::endl;
	//
}
/***********************************************/
char strLicenseGIC[] = "THIS SOFTWARE IS PROVIDED \"AS-IS\". THERE IS NO WARRANTY OF ANY KIND. "
"NEITHER THE AUTHORS NOR THE OHIO STATE UNIVERSITY WILL BE LIABLE FOR "
"ANY DAMAGES OF ANY KIND, EVEN IF ADVISED OF SUCH POSSIBILITY. \n"
"\n"
"This software was developed (and is copyrighted by) the Jyamiti group at "
"The Ohio State University. Please do not redistribute this software. "
"This program is for academic research use only. This software uses the "
"Boost library (www.boost.org) and Ann library "
"(www.cs.umd.edu/~mount/ANN/) which are covered under their own licenses.\n"
"\n"
"The Boost library's license "
"(which applies to the Boost library ONLY and NOT to this program itself) is "
"as follows:\n"
"\n"
"LICENSE\n"
"---------------------------------------------------------------------------\n"
"Boost Software License - Version 1.0 - August 17th, 2003\n"
"\n"
"Permission is hereby granted, free of charge, to any person or organization "
"obtaining a copy of the software and accompanying documentation covered by "
"this license (the \"Software\") to use, reproduce, display, distribute, "
"execute, and transmit the Software, and to prepare derivative works of the "
"Software, and to permit third-parties to whom the Software is furnished to "
"do so, all subject to the following: \n"
"\n"
"The copyright notices in the Software and this entire statement, including "
"the above license grant, this restriction and the following disclaimer, "
"must be included in all copies of the Software, in whole or in part, and "
"all derivative works of the Software, unless such copies or derivative "
"works are solely in the form of machine-executable object code generated by "
"a source language processor. \n"
"\n"
"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR "
"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, "
"FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT "
"SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE "
"FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, "
"ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER "
"DEALINGS IN THE SOFTWARE. \n"
"---------------------------------------------------------------------------\n"
"\n"
"ANN library's license "
"(which applies to the ANN library ONLY and NOT to this program itself) is "
"as follows: \n"
"\n"
"LICENSE\n"
"---------------------------------------------------------------------------\n"
"The ANN Library (all versions) is provided under the terms and "
"conditions of the GNU Lesser General Public Library, which is stated "
"below.  It can also be found at: \n"
"\n"
"   http:////www.gnu.org/copyleft/lesser.html \n"
"---------------------------------------------------------------------------\n";
/**********************************************************************/


bool ParseCommand(int argc, char** argv,
				int &input_type,
				std::string &InputFile,
				std::string &OutputFile,
				bool &graph_distance_flag,
				float &eps_dist,
				float &delta_dist,
				int &max_dimension)
{
	try
	{
		// Define the program options description
//		namespace po = boost::program_options;
		boost::program_options::options_description desc("GIComplex Usage");
		desc.add_options()
			(",h", "Help information;")
			(",l", "License information;")
			("type", boost::program_options::value<int>(&input_type)->required(), "Input data type: point cloud (0) or weighted graph (1);")
			(",I", boost::program_options::value<std::string>(&InputFile)->required(), "Input file name;")
			(",O", boost::program_options::value<std::string>(&OutputFile)->required(), "Output file name prefix;")
			(",g", boost::program_options::value<bool>(&graph_distance_flag)->default_value(false), "Subsample by Euclidean distance (false) or graph distance (true); for weighted graph input, it's always (true);")
			("epsilon", boost::program_options::value<float>(&eps_dist), "Parameter for building the Rips graph on input samples; Not used for weighted graph input;")
			("delta", boost::program_options::value<float>(&delta_dist)->required(), "Parameter for generating subsamples from input samples;")
			("dim", boost::program_options::value<int>(&max_dimension)->required(), "The maximum dimension of simplices in output complex;");

		// Parser map
		boost::program_options::variables_map vm;
		try
		{
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);

			//
			if (vm.count("-h"))
			{
				std::cout << desc << std::endl;
			}
			//
			if (vm.count("-l"))
			{
				std::cout << strLicenseGIC << std::endl;
			}
			//
			boost::program_options::notify(vm);
		}
		catch(boost::program_options::required_option& e)
		{
			std::cerr<< "ERROR: " << e.what() << std::endl;
			return false;
		}
		catch(boost::program_options::error& e)
		{
			std::cerr<< "ERROR: " << e.what() << std::endl;
			return false;
		}
	}
	catch(std::exception& e)
	{
		std::cerr << "Unhandled Exception reached the top of main: "
					<< e.what() << ", application will now exit" << std::endl;
		return false;

	}
	return true;
}

/*********************/
/*int main(int argc, char **argv)
{
	//
	std::string InputFileName;
	std::string OutputFileName;
	int input_type = 0;
	float eps_dist = .05f;
	float delta_dist = 0.9f;
	bool graph_distance_flag = true;
	int max_dimension = 3;

	//
	if (ParseCommand(argc, argv, input_type, InputFileName, OutputFileName, graph_distance_flag, eps_dist, delta_dist, max_dimension))
	{
		if (!input_type)
		{
			PointSet pts;
			//pts.SampleCrossLine(0.02, 1.0);
			pts.ReadPointsFromFile(InputFileName.c_str());
			std::map<int, int> SubPointIndex;
			SimpleGraph biColoredGraph;
			//
			std::cout << std::endl << "Construct graph induced compolex " << std::endl << std::endl;
			//
			ComputeBiColoredGraph(eps_dist, delta_dist, graph_distance_flag, pts, SubPointIndex, biColoredGraph);
			//
			std::string outfile_name(OutputFileName);
			//
			GIComplex gic(max_dimension, biColoredGraph.color_number, &biColoredGraph);
			gic.Construction();
			//
			outfile_name = outfile_name + "_stat.txt";
			gic.WriteStatisticsToFile(outfile_name.c_str());
			//
			outfile_name = OutputFileName + "_complex.txt";
			WriteOFFformatComplex(outfile_name.c_str(), pts, SubPointIndex, gic);
		}
		else
		{// read weighted graph
			std::vector<std::vector<float> > pts;
			std::map<int, int> SubPointIndex;
			SimpleGraph biColoredGraph;
			//
			ComputeBiColoredGraph(InputFileName.c_str(), delta_dist, pts, SubPointIndex, biColoredGraph);
			//
			//
			std::cout << std::endl << "Construct graph induced compolex " << std::endl << std::endl;
			//
			std::string outfile_name(OutputFileName);
			//
			GIComplex gic(max_dimension, biColoredGraph.color_number, &biColoredGraph);
			gic.Construction();
			//
			outfile_name = outfile_name + "_stat.txt";
			gic.WriteStatisticsToFile(outfile_name.c_str());
			//
			outfile_name = OutputFileName + "_complex.txt";
			WriteOFFformatComplex(outfile_name.c_str(), pts, SubPointIndex, gic);
		}
		//pts.WriteBackToFile("test.txt");
		//biColoredGraph.WriteBackToFile("test_graph.txt");
	}
	//
	return 0;
}*/
