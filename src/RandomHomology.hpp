//
//  RandomHomology.hpp
//  BitTests
//
//  Created by Bill Varcho on 2/24/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#ifndef RandomHomology_hpp
#define RandomHomology_hpp


#include "Utilities.hpp"

using std::set;

class RandomHomology {
public:
    int dim;
    MortonCode *mc;
    ANNkd_tree *kd_tree;
    set<int>removed_ids;
    vector<vector<int>> indices;
    vector<vector<float>> initial_points;
    RandomHomology(vector<vector<float>> &points,
                   vector<float> &min_bounds,
                   vector<float> &max_bounds, int dim);
    RandomHomology(GIC &g);
    ~RandomHomology() {};
    void run(double alpha);
    void run(double alpha, vector<Operation*> &c);
    void collapseToClosest(MortonPoint p, vector<Operation*> &collapses);
};

#endif /* RandomHomology_hpp */
