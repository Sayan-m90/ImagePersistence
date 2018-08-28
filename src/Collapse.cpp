//
//  Collapse.cpp
//  RandomHomology
//
//  Created by Bill Varcho on 3/5/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#include "Collapse.hpp"

/*----------------------*
 |  Insert Operation
 *----------------------*/

Insert::Insert() {
}

Insert::Insert(int id) {
    v_id.push_back(id);
}

Insert::Insert(vector<int> ids) {
    for (int id : ids) {
        v_id.push_back(id);
    }
}

void Insert::Print() {
    cout << PrintString();
}

string Insert::PrintString() {
    string s = "i ";
    for (int i = v_id.size() - 1; i > 0; i--) {
        s = s + to_string(v_id[i]) + " ";
    }
    return s + to_string(v_id[0]) + "\n";
}

OP_TYPE Insert::Type() {
    return INSERT_OP;
}

bool Insert::IsVertexInsert() {
    return v_id.size() == 1;
}


/*----------------------*
 |  Collapse Operation
 *----------------------*/

Collapse::Collapse() {
    v_start = -1;
    v_target = -1;
}

Collapse::Collapse(int start, int target) {
    v_start = start;
    v_target = target;
}

void Collapse::Print() {
    cout << PrintString();
}

string Collapse::PrintString() {
    return "c " + to_string(v_start) + " t " + to_string(v_target) + "\n";
}


OP_TYPE Collapse::Type() {
    return COLLAPSE_OP;
}


/*----------------------*
 |  Timestamp Operation
 *----------------------*/

Timestamp::Timestamp() {
    v_time = -1;
}

Timestamp::Timestamp(float time) {
    v_time = time;
}

void Timestamp::Print() {
    cout << PrintString();
}

string Timestamp::PrintString() {
    return "# " + to_string(v_time) + "\n";
}


OP_TYPE Timestamp::Type() {
    return TIME_OP;
}
