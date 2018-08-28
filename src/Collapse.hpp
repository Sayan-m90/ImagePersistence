//
//  Collapse.hpp
//  RandomHomology
//
//  Created by Bill Varcho on 3/5/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#ifndef Collapse_hpp
#define Collapse_hpp

#include <stdio.h>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::to_string;

enum OP_TYPE {INSERT_OP, COLLAPSE_OP, TIME_OP, NULL_OP};

class Operation {
public:
    virtual void Print() = 0;
    virtual string PrintString() = 0;
    virtual OP_TYPE Type() = 0;
    virtual bool IsVertexInsert() {return false;}
};

class Insert : public Operation {
public:
    vector<int> v_id;
    Insert();
    Insert(int v_id);
    Insert(vector<int> ids);
    ~Insert() {}
    virtual void Print();
    virtual bool IsVertexInsert();
    virtual string PrintString();
    virtual OP_TYPE Type();
};

class Collapse : public Operation {
public:
    int v_start, v_target;
    Collapse();
    Collapse(int start, int end);
    ~Collapse() {}
    virtual void Print();
    virtual string PrintString();
    virtual OP_TYPE Type();
};

class Timestamp : public Operation {
public:
    float v_time;
    Timestamp();
    Timestamp(float time);
    ~Timestamp() {}
    virtual void Print();
    virtual string PrintString();
    virtual OP_TYPE Type();
};

#endif /* Collapse_hpp */
