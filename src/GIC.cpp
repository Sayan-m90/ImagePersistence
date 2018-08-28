//
//  GIC.cpp
//  RandomHomology
//
//  Created by Bill Varcho on 3/1/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#include "GIC.hpp"

GIC::GIC() {}

GIC::GIC(string _fp) {
    fp = _fp;
    ifstream input(fp);     //input is the complex file generated
    
    bool done = false;
    int num_points;
    float min_n = std::numeric_limits<float>::max();
    float max_n = std::numeric_limits<float>::min();
    
    for (int i = 0; i < Constants::GIC_MAX_DIM_SIZE; i++) {
        vector<int> ndim_indicies;
        indices.push_back(ndim_indicies);
    }
    
    if (input.is_open()) {
        // get the first line (dim, num pts)
        string s_dim, s_num, tmp;
        stringstream ss(stringstream::in | stringstream::out);
        getline(input, s_dim, ' ');
        getline(input, s_num);
        ss << s_dim << s_num;
        
        dim = stoi(s_dim);
        num_points = stoi(s_num);
        
        int point_count = 0;
        while (point_count < num_points) {
            done = input.eof();
            std::string sx, sy, sz;
            //float x, y, z;
            std::stringstream ss(std::stringstream::in |
                                 std::stringstream::out);
            string s;
            getline(input, s);
            
            vector<string> strs;
            boost::split(strs, s, boost::is_any_of(Constants::SEPARATOR_STR));
            vector<float> pt;
            for (size_t i = 0; i < strs.size(); i++) {
                if (strs[i].size() > 0) {
                    float coord = stod(strs[i]);
                    pt.push_back(coord);
                    min_n = min(min_n, coord);
                    max_n = max(max_n, coord);
                }
            }
            
           // min_n = min(min_n, min(x, min(y, z)));
           // max_n = max(max_n, max(x, max(y, z)));
            
            // add point to the vector
            pts.push_back(pt);
            
            // add point to indices
            indices[0].push_back(point_count);
            
            point_count += 1;
        }
        int coord;
        
        // now read in the indices
        while (!input.eof()) {
            string s;
            vector<string> strs;
            getline(input, s);

            boost::split(strs, s, boost::is_any_of(Constants::SEPARATOR_STR));
            vector<int> ids;
            for (size_t i = 1; i < strs.size(); i++) {
                if (strs[i].size() > 0 ) {  
                    coord = stoi(strs[i]);
                    ids.push_back(coord);
                    //fflush(stdin);
                    //cout<<strs.size()<<":"<<coord<<"      ";
                    //getchar(); 
                }
            }
           //cout<<"\n";
            for (int i : ids) {
                indices[ids.size() - 1].push_back(i);
            }
        }
    }
    
    bb_max = max_n;
    bb_min = min_n;
}
