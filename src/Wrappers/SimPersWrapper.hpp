//
//  SimPersWrapper.hpp
//  RandomHomology
//
//  Created by Bill Varcho on 3/10/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#ifndef SimPersWrapper_hpp
#define SimPersWrapper_hpp


#include "../Collapse.hpp"
#include "../Constants.hpp"
#include "../SimPers/SimplicialComplex.h"
#include <iostream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <sstream>
#include <unordered_set>
#include <ctime>
#include <fstream>
#include <boost/program_options.hpp>

using std::vector;

class SimpersWrapper {
public:
    static int Run(vector<Operation *> &c,
                    string fp);
};

#endif /* SimPersWrapper_hpp */
