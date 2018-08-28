//
//  NormalGraph.hpp
//  RandomHomology
//
//  Created by Bill Varcho on 4/8/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#ifndef NormalGraph_hpp
#define NormalGraph_hpp

#include <cmath>
#include <stdio.h>
#include "Graph.hpp"
#include "../MortonCode.hpp"
#include "../Utilities.hpp"

/*--------------------------------------------------------
 | Normal Graph is a graph consructed with the help
 | of a Morton Code ordering of the points. Edges are
 | then added by stepping throught the ordering, and for
 | each point looking at points near it in the ordering.
 | A Gaussian Function is then fit to the distances to 
 | these neighbors from the current point. Edges are then
 | added between points with distance within one standard
 | deviation of the fitted Gaussian.
 *--------------------------------------------------------*/

using std::cout;
using std::endl;

class NormalGraph : public Graph {
private:
    int K;
    void Construct();
public:
    NormalGraph(string fp, int K);
    NormalGraph(int K);
    NormalGraph(vector<vector<float>> pts, int K);
    ~NormalGraph() {}
};


#endif /* NormalGraph_hpp */
