//
//  GIC.hpp
//  RandomHomology
//
//  Created by Bill Varcho on 3/1/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#ifndef GIC_hpp
#define GIC_hpp

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <boost/algorithm/string.hpp>
#include "Constants.hpp"

using std::string;
using std::vector;
using std::ifstream;
using std::stringstream;
using std::cout;
using std::endl;
using std::min;
using std::max;
using std::exception;

// As of right now this class is being strictly used for visualization,
// however I plan to use in the reading in of data and construction of
// morton ordering later
class GIC {
public:
    GIC();
    GIC(string fp);
    string fp;
    int dim;
    double bb_min, bb_max;
    vector<vector<float>> pts;
    vector<vector<int>> indices;
};

#endif /* GIC_hpp */
