//
//  Utilities.hpp
//  BitTests
//
//  Created by Bill Varcho on 2/24/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#ifndef Utilities_hpp
#define Utilities_hpp

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <thread>
#include <vector>

#include "ANN/ANN.h"
#include "MortonCode.hpp"
#include "Collapse.hpp"
#include "Constants.hpp"
#include "GIC.hpp"
#include "MortonCode.hpp"


using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::cout;
using std::endl;
using std::min;
using std::max;
using std::sort;

class Utilities {
public:
    static void ReadInPoints(vector<vector<double>> &pts, string fp);
    static void ReadInPoints(vector<vector<float>> &pts, string fp);
    // Used mainly for testing
    static void CreateRandomPoints(vector<vector<double>> &pts,
                                   vector<double> &mins,
                                   vector<double> &maxes,
                                   int dim, int n);
    static void WriteCollapsesToFile(string fp, vector<Operation*> &collapses);
    template <typename T>
    static ANNkd_tree* ConstructKDTree(vector<vector<T>> pts, int dim) {
        int max_points = pts.size();
        ANNpointArray data_points;
        data_points = annAllocPts(max_points, dim);
        
        // pts to data points
        for (int i = 0; i < max_points; i++) {
            ANNpoint new_p = annAllocPt(dim);
            vector<T> current_pt = pts[i];
            for (int d = 0; d < dim; d++) {
                new_p[d] = current_pt[d];
            }
            data_points[i] = new_p;
        }
        // create and return the kdtree
        return new ANNkd_tree(data_points,
                              max_points,
                              dim);
    }
};

#endif /* Utilities_hpp */
