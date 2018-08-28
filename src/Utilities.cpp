//
//  Utilities.cpp
//  BitTests
//
//  Created by Bill Varcho on 2/24/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#include "Utilities.hpp"

#define EPSILON 1.0e-5

void Utilities::CreateRandomPoints(vector<vector<double>> &pts,
                                   vector<double> &mins,
                                   vector<double> &maxes,
                                   int dim, int n) {
    for (int i = 0; i < n; i++) {
        vector<double> point;
        for (int j = 0; j < dim; j+=1) {
            double coord;
            coord = rand() / ((double) INT_MAX);
            coord = 2*coord - 1;
            point.push_back(coord);
            if (mins.size() > 0) {
                if (coord < mins[0]) {
                    mins[0] = coord;
                }
            } else { mins.push_back(coord);}
            if (maxes.size() > 0) {
                if (coord > maxes[0]) {
                    maxes[0] = coord;
                }
            } else { maxes.push_back(coord);}
        }
        pts.push_back(point);
    }
}

void Utilities::ReadInPoints(vector<vector<double>> &pts, string fp) {
    ifstream input(fp);
    if (input.is_open()) {
        while (!input.eof()) {
            std::stringstream ss(std::stringstream::in |
                                 std::stringstream::out);
            string s;
            getline(input, s);
            vector<string> strs;
            boost::split(strs, s, boost::is_any_of(Constants::SEPARATOR_STR));
            vector<double> pt;
            for (size_t i = 0; i < strs.size(); i++) {
                if (strs[i].size() > 0) {
                    double coord = stod(strs[i]);
                    pt.push_back(coord);
                }
            }
            // add point to the vector
            if (pt.size() > 0) {
                pts.push_back(pt);
            }
        }
    }
    input.close();
}

void Utilities::ReadInPoints(vector<vector<float>> &pts, string fp) {
    ifstream input(fp);
    if (input.is_open()) {
        while (!input.eof()) {
            std::stringstream ss(std::stringstream::in |
                                 std::stringstream::out);
            string s;
            getline(input, s);
            vector<string> strs;
            boost::split(strs, s, boost::is_any_of(Constants::SEPARATOR_STR));
            vector<float> pt;
            for (size_t i = 0; i < strs.size(); i++) {
                if (strs[i].size() > 0) {
                    float coord = stof(strs[i]);
                    pt.push_back(coord);
                }
            }
            // add point to the vector
            if (pt.size() > 0) {
                pts.push_back(pt);
            }
        }
    }
    input.close();
}

void Utilities::WriteCollapsesToFile(string fp, vector<Operation*> &collapses) {
    cout << "Writing Collapses to File: " + fp +"_collapses"<< endl;
    ofstream out_file;
    int num_simplices = 0;
    out_file.open (fp + "_collapses");
    for (const auto t: collapses) {
        out_file << t->PrintString();
        if (t->IsVertexInsert()) {
            num_simplices += 1;
        }
    }
    out_file.close();
    
    cout << "Writing iDC File: " + fp << endl;
    ofstream out_file_iDC;
    out_file_iDC.open (fp + "_iDC");
    out_file_iDC << num_simplices << "\n";
    for (int i = 0; i < num_simplices; i++) {
        out_file_iDC << i << "\n";
    }
    out_file_iDC.close();
}

