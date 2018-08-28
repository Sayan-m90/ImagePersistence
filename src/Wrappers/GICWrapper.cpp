//
//  GICWrapper.cpp
//  RandomHomology
//
//  Created by Bill Varcho on 3/10/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#include "GICWrapper.hpp"

void ComputeBiColoredGraph(const char* pFileName,
                           const float delta_dist,
                           std::vector<std::vector<float> > &pts,
                           std::map<int, int> &SubPointIndex,
                           SimpleGraph &biColoredGraph);


void WriteComplex(const char* pFileName,
                  vector<vector<float> > &pts,
                  std::map<int, int> &subPointIndices,
                  GIComplex &gic);

void GICWrapper::SimpleGraphFromKNN(Graph &k_g,
                                    SimpleGraph &g,
                                    vector<vector<float>> &pts) {

    g.InitNodes(k_g.GetNumPoints());
    pts = k_g.GetPoints();
    
    for (pair<int, int> p : k_g.GetEdges()) {
        vector<float> p1 = pts[p.first];
        vector<float> p2 = pts[p.second];
        float sum = 0.0;
        for (int i = 0; i < p1.size(); i++) {
            sum += pow((p1[i] - p2[i]), 2.0);
            //cout<<"i: "<<i<<" p1[i]: "<<p1[i]<<" p2[i]: "<<p2[i]<<" sum: "<<sum<<endl;
            //getchar();
        }
        //cout<<"here"<<k_g.GetDim();
        //cout<<endl;
        //cout<<" p.first "<<p.first<<" p.second "<<p.second<<" sum: "<<sum<<endl;
        g.AddEdge(p.first, p.second, sqrt(sum));
    }
    
}

void GICWrapper::CreateBiColoredGraph(string point_file,
                                      SimpleGraph &bi_colored_graph,
                                      map<int, int> &sub_point_index,
                                      vector<int> &color_map,
                                      vector<vector<float>> &subsampled_pts,
                                      int neighbors) {
    
    SimpleGraph original_graph;
    vector<vector<float>> pts;
    
    // 1. Read in Data Points
    // 2. Create Graph on original data points
    KNN_Graph *g = new KNN_Graph(point_file, neighbors);
    
    // Normal Graph expirimentation led to unsatisfying results :(
    // NormalGraph *g = new NormalGraph(point_file, neighbors);
    
    int max_dim = (int) g->GetPoints()[0].size();
    SimpleGraphFromKNN(*g, original_graph, pts);
    //cout << "Vec Node Size:" << original_graph.vecNode.size() << endl;
    
    // 3. Subsample Data Points
    // (TODO: add sampling parameter instead of hard coding)
    MortonCode *mc = new MortonCode(g->GetPoints());
    vector<vector<float>> graph_points = g->GetPoints();
    for (int i = 0; i < graph_points.size(); i++) {
        MortonPoint n = mc->nextPoint();
        if (i%10 != 0 & i%10 != 5) { //  somewhat arbitrary, should include
            subsampled_pts.push_back(n.point);
            sub_point_index.insert(pair<int,int>(subsampled_pts.size() - 1, i));
        }
    }
    
    // 4. Creat ANN structure for Subsample DP
    ANNkd_tree *kd_tree = Utilities::ConstructKDTree(subsampled_pts, max_dim);
    
    // 5. Generate SubPointIndex + Generate Color Mapping
    for (int i = 0; i < graph_points.size(); i++) {
        // do an ann query to get index that this point is mapped to
        ANNpoint query_point;
        ANNidxArray nn_ids;
        ANNdistArray dists;
        int k = 1;
        double eps = .01;
        query_point = annAllocPt(max_dim);
        nn_ids = new ANNidx[k];
        dists = new ANNdist[k];
        
        vector<float> pt = graph_points[i];
        for (int d = 0; d < max_dim; d++) {
            query_point[d] = pt[d];
        }
        
        kd_tree->annkSearch(query_point,
                            k,
                            nn_ids,
                            dists,
                            eps);
        color_map.push_back(nn_ids[0]);
    }

    ANNSearch::SetColorMappingAndExtractColoredGraph(color_map,
                                                     original_graph,
                                                     bi_colored_graph);

    bi_colored_graph.color_number = (int)sub_point_index.size();
}

void GICWrapper::Run(string point_file, string out_file,
                     int max_dim, int neighbors) {
    map<int, int> sub_point_index;
    vector<int> color_map;
    vector<vector<float>> subsampled_pts;
    
    SimpleGraph bi_colored_graph;
    CreateBiColoredGraph(point_file,
                         bi_colored_graph,
                         sub_point_index,
                         color_map,
                         subsampled_pts,
                         neighbors);

    GIComplex gic(max_dim, bi_colored_graph.color_number, &bi_colored_graph);
    gic.Construction();
    
    // TODO(me): not writing out stats
    // string outfile_name = out_file + "_stat.txt";
    // gic.WriteStatisticsToFile(outfile_name.c_str());

    // outfile_name = out_file; // + "_complex.txt";
    
    WriteComplex(out_file.c_str(),
                          subsampled_pts,
                          sub_point_index,
                          gic);
}

// Writes the Complex to a file. Largely Fangtao's Code, with some modifications
// made specifically for this new situation
void WriteComplex(const char* pFileName,
                  vector<vector<float> > &pts,
                  std::map<int, int> &subPointIndices,
                  GIComplex &gic)
{
    std::cout << "Writing < " << pFileName << " >" << std::endl;
    std::ofstream ofile;
    ofile.open(pFileName, std::ifstream::out);
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    if (ofile.is_open())
    {
        // dim # of points
        sstr << pts.front().size() << " " << subPointIndices.size() << endl;
        
        // points
        for (vector<float> pt : pts) {
            for (int i = 0; i < pt.size(); i++) {
                sstr << pt[i];
                if (i < pt.size() - 1) {
                    sstr << " ";
                } else {
                    sstr << endl;
                }
            }
        }
        
        for (int i = 1; i < gic.dim + 1; i++)
        {
            //depth = i + 1; // want to visit all simplices at this dimension
            if (gic.Simplicies_Cnt[i] != 0)
            {// visit each simplex in dimension i
                for (int vid = 0; vid < gic.Simplicies_Cnt[0]; vid++)
                {
                    if (gic.head_circular_list_in_each_dim[vid][i - 1])
                        // i-1 as vertex are stored in an array
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
        ofile.write(sstr.str().c_str(), sstr.str().size());
        ofile.close();
        sstr.clear();
    }
    else
    {
        cout << "Can NOT open file " << pFileName << endl;
        exit(0);
    }
    cout << "---- Done GIC---- " << endl << endl;
}