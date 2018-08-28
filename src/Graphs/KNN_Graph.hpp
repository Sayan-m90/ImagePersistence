//
//  KNN_Graph.hpp
//  RandomHomology
//
//  Created by Bill Varcho on 3/8/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#ifndef KNN_Graph_hpp
#define KNN_Graph_hpp

#include "ANN/ANN.h"
#include "../Utilities.hpp"
#include "Graph.hpp"
#include <set>
#include <stdio.h>
#include <vector>

using std::cout;
using std::endl;
using std::pair;
using std::set;
using std::string;
using std::vector;

/*---------------------------------------------------
 | Temp class for creating a graph on N-D points
 | based on approximate K-Nearest Neighbors. Proof
 | of concept for using this graph as input to GIC
 | as initial graph instead of epsilon-delta graph.
 *---------------------------------------------------*/

class KNN_Graph : public Graph {
private:
    int dim, K;
    ANNkd_tree* KD_tree;
    void Construct();
public:
    KNN_Graph(int K);
    KNN_Graph(string fp, int K);
    KNN_Graph(vector<vector<float>> pts, int K);
    ~KNN_Graph() {}
    void writeToFile(string filepath);
    void writeWeightedToFile(string filepath);
};


#endif /* KNN_Graph_hpp */
