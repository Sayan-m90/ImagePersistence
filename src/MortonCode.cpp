//
//  MortonCode.cpp
//  BitTests
//
//  Created by Bill Varcho on 2/24/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#include "MortonCode.hpp"

// Custom comparator for the large bitstrings create via the encoding
bool operator<(const MortonPoint &a, const MortonPoint &b) {
    for (int i = 0; i < a.code.size() - 1; i += 1) {
        if (a.code[i] != b.code[i]) {
            return a.code[i] < b.code[i];
        }
    }
    return a.code[a.code.size() - 1] < b.code[b.code.size() - 1];
}

/*-------------------------------------*
 |  Morton Point Constructors          |
 *-------------------------------------*/
//int MortonPoint::id = 0;

MortonPoint::MortonPoint(const vector<float> pt, int min, int max, int id) {
    point = pt;
    vector<uint64_t> scaled_coords;
    for (int i = 0; i < pt.size(); i++) {
        scaled_coords.push_back(
                    (uint64_t) Constants::MORTON_CODE_MULTIPLIER * ((max - min)
                    + pt[i]) / (2.));
        code.push_back(0);
    }

    // now need to set the Morton Code of the point
    int index = 0;  // current position in 64 bit string
    int count = 0;  // which index in code we are writing to
    int coord_index = 0;
    for (int j = 0; j< 64 * pt.size(); j += 1) {
        uint64_t mask = 1;
        mask = mask << (63 - index);
        uint64_t bit = (mask & scaled_coords[coord_index]) >> (63 - index);

        // now that we have the proper bit set the code
        code[count] = (code[count] << 1) | bit;

        if (j % pt.size() == pt.size() - 1) { index += 1; }
        if (j % 64 == 63) { count += 1;}
        // decrement coord index
        coord_index = (coord_index + 1) % (pt.size());
    }
    p_id = id;
    //id++;
}

/*-----------------------------------*
 |  Morton Code Constructors         |
 *-----------------------------------*/

MortonCode::MortonCode(const vector<vector<float>> &points) {
    bbox_max = 1;
    bbox_min = -1;
    priority_queue<MortonPoint> queue;
    int i = 0;
    for (vector<float> pt : points) {
        add(pt, i);
        i++;
    }
}

MortonCode::MortonCode(const vector<vector<float>> &points,
                       const vector<float> &min_bounds,
                       const vector<float> &max_bounds) {
    bbox_min = *min_element(min_bounds.begin(), min_bounds.end()) + .1;
    bbox_max = *max_element(min_bounds.begin(), min_bounds.end()) - .1;
    priority_queue<MortonPoint> queue;
    int i = 0;
    for (vector<float> pt : points) {
        add(pt, i);
        i++;
    }
}

MortonCode::MortonCode(GIC g) {
    bbox_max = g.bb_max;
    bbox_min = g.bb_min;
    int i = 0;
    priority_queue<MortonPoint> queue;
    for (vector<float> pt : g.pts) {
        add(pt, i);
        i++;
    }
}

/*-----------------------------------*
 |  Morton Code Class Methods        |
 *-----------------------------------*/

void MortonCode::add(vector<float> pt, int id) {
    MortonPoint mp = MortonPoint(pt, bbox_min, bbox_max, id);
    queue.push(mp);
}

void MortonCode::add(MortonPoint pt) {
    queue.push(pt);
}

vector<float> MortonCode::next() {
    if (!empty()) {
        MortonPoint next_point = queue.top();
        queue.pop();
        return next_point.point;
    } else {
        vector<float> x;
        x.push_back(-1);
        return x;
    }
}

void MortonCode::setPoints(vector<MortonPoint> pts) {
    priority_queue<MortonPoint> new_queue;
    queue.swap(new_queue);
    for (MortonPoint p : pts) {
        queue.push(p);
    }
}

MortonPoint MortonCode::nextPoint() {
    MortonPoint next_point = queue.top();
    queue.pop();
    return next_point;
}

bool MortonCode::empty() {
    return !(queue.size() > 0);
}

int MortonCode::size() {
    return (int) queue.size();
}

vector<MortonPoint> MortonCode::GetOrdering() {
    vector<MortonPoint> order;
    while (!empty()) {
        order.push_back(nextPoint());
    }
    return order;
}