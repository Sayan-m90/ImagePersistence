//
//  Graph.hpp
//  RandomHomology
//
//  Created by Bill Varcho on 4/8/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#ifndef Graph_hpp
#define Graph_hpp

#include <set>
#include <vector>

using std::pair;
using std::set;
using std::vector;

/*---------------------------------------------------
 | Graph (Abstract) Base Class
 *---------------------------------------------------*/

class Graph {
public:
    int dim;
    vector<vector<float>> points;
    set<pair<int, int>> edges;
    virtual void Construct() = 0;
    int GetNumPoints() {                return points.size();}
    int GetNumEdges() {                 return edges.size();}
    int GetDim() {                      return dim;}
    set<pair<int, int>> GetEdges() {    return edges;}
    vector<vector<float>> GetPoints() { return points;}
    void SetPoints(vector<vector<float>> pts_) {points = pts_;}
};

#endif /* Graph_hpp */
