//
//  NormalGraph.cpp
//  RandomHomology
//
//  Created by Bill Varcho on 4/8/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#include "NormalGraph.hpp"

NormalGraph::NormalGraph(int K_) {
    K = K_;
}

NormalGraph::NormalGraph(vector<vector<float>> pts_, int K_) {
    cout << pts_.size() << endl;
    points = pts_;
    K = K_;
    Construct();
}

void NormalGraph::Construct() {
    // Set Dim
    cout << points.size() << endl;
    dim = points[0].size();

    // Create Morton Code ordering of points
    MortonCode *mc = new MortonCode(points);
    vector<MortonPoint> order = mc->GetOrdering();
    
    // TODO (me): implement
    // brute force implementation for testing
    for (int i = 0; i < order.size(); i++) {
        vector<MortonPoint> neighbors;
        // get neighbors behind
        int num_behind = min(K, i);
        for (int j = 0; j < num_behind; j++) {
            neighbors.push_back(order[i-j-1]);
        }
        
        // get neighbors in front
        int num_infront = min(K, (int) order.size() - i - 1);
        for (int j = 0; j < num_infront; j++) {
            neighbors.push_back(order[i+j+1]);
        }
    
        float mu = 0.0;
        float sigma_sq = 0.0;
        
        // calculate mu
        for (MortonPoint p : neighbors) {
            mu += sqrt( pow(p.point[0] - order[i].point[0], 2.0)
                       + pow(p.point[1] - order[i].point[1], 2.0)
                       + pow(p.point[2] - order[i].point[2], 2.0));
        }
        mu = mu / neighbors.size();
        
        // calculate sigma^2
        for (MortonPoint p : neighbors) {
            float dist = sqrt( pow(p.point[0] - order[i].point[0], 2.0)
                       + pow(p.point[1] - order[i].point[1], 2.0)
                       + pow(p.point[2] - order[i].point[2], 2.0));
            
            sigma_sq += pow(mu - dist, 2);
        }
        sigma_sq = sigma_sq / neighbors.size();
        
        //cout << mu << " " << sigma_sq << endl;
        
        // add edge if within one s.d.
        float sigma = sqrt(sigma_sq);
        for (MortonPoint p : neighbors) {
            float dist = sqrt( pow(p.point[0] - order[i].point[0], 2.0)
                              + pow(p.point[1] - order[i].point[1], 2.0)
                              + pow(p.point[2] - order[i].point[2], 2.0));
            
            if (abs(mu - dist) < sigma) {
                // create edge
                edges.insert(pair<int, int>(p.p_id, order[i].p_id));
            }
        }
    }
}


NormalGraph::NormalGraph(string fp, int K_) {
    K = K_;
    vector<vector<float>> pts;
    Utilities::ReadInPoints(pts, fp);
    points = pts;
    Construct();
}

