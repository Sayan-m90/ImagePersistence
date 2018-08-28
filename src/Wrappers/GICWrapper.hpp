//
//  GICWrapper.hpp
//  RandomHomology
//
//  Created by Bill Varcho on 3/10/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#ifndef GICWrapper_hpp
#define GICWrapper_hpp

#include <stdio.h>
#include "../Graphs/Graph.hpp"
#include "../Graphs/KNN_Graph.hpp"
#include "../Graphs/NormalGraph.hpp"
#include "../GIC.hpp"
#include "../GIComplex/PointSet.h"
#include "../GIComplex/SimpleGraph.h"
#include "../GIComplex/GIComplex.h"
#include "../GIComplex/ANNSearchSampling.h"
#include "ANN/ANN.h"
#include "ANN/ANNperf.h"
#include "ANN/ANNx.h"


using std::map;
using std::pair;

class GICWrapper {
public:
    static void Run(Graph *g, GIC &output);
    static void Run(string point_file, string out_gic_file,
                    int dim, int neighbors);
private:
    static void CreateBiColoredGraph(string point_file,
                                     SimpleGraph &bi_colored_graph,
                                     map<int, int> &sub_point_index,
                                     vector<int> &color_map,
                                     vector<vector<float>> &subsampled_pts,
                                     int neighbors);
    static void SimpleGraphFromKNN(Graph &k_g,
                                   SimpleGraph &g,
                                   vector<vector<float>> &pts);
};

#endif /* GICWrapper_hpp */
