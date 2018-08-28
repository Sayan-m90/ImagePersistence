//
//  Constants.cpp
//  RandomHomology
//
//  Created by Bill Varcho on 3/6/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#include "Constants.hpp"

const int Constants::GIC_MAX_DIM_SIZE = 5;
const int Constants::MAX_NUM_COLLAPSES = 90;
const int Constants::MAX_NUM_BARCODES = 100;
const int Constants::MORTON_CODE_MULTIPLIER = 100000;
const int Constants::MOUSE_SPEED = 100;
const int Constants::VIEWER_POS_INF_VALUE = INT32_MAX;
const string Constants::SEPARATOR_STR = "\t ";
const char *Constants::VERTEX_SHADER = "vertex.shader";
const char *Constants::FRAGMENT_SHADER =  "fragment.shader";
const char *Constants::VERTEX_SHADER_GIC = "vertex_gic.shader";
const char *Constants::FRAGMENT_SHADER_GIC = "fragment_gic.shader";
const char *Constants::VERTEX_SHADER_BARCODE = "vertex_barcode.shader";
const char *Constants::FRAGMENT_SHADER_BARCODE = "fragment_barcode.shader";