//
//  Constants.hpp
//  RandomHomology
//
//  Created by Bill Varcho on 3/6/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#ifndef Constants_hpp
#define Constants_hpp

#include <stdio.h>
#include <string>

using std::string;

class Constants {
public:
    const static char* VERTEX_SHADER;
    const static char* FRAGMENT_SHADER;
    const static char* VERTEX_SHADER_GIC;
    const static char* FRAGMENT_SHADER_GIC;
    const static char* VERTEX_SHADER_BARCODE;
    const static char* FRAGMENT_SHADER_BARCODE;
    const static int GIC_MAX_DIM_SIZE;
    const static int MORTON_CODE_MULTIPLIER;
    const static int MOUSE_SPEED;
    const static int MAX_NUM_COLLAPSES;
    const static int MAX_NUM_BARCODES;
    const static int VIEWER_POS_INF_VALUE;
    const static string SEPARATOR_STR;
};

#endif /* Constants_hpp */
