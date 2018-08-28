//
//  KNN_Graph.cpp
//  RandomHomology
//
//  Created by Bill Varcho on 3/8/16.
//  Copyright © 2016 Bill Varcho. All rights reserved.
//

#include "KNN_Graph.hpp"

/*----------------*
 | Constructors
 *----------------*/

void KNN_Graph::Construct() {
    if (points.size() > 0) {
        dim = points[0].size();
        K = min(K, (int) points.size() - 1);
        cout << K << endl;
        
        // create ANN_KD tree from points
        KD_tree = Utilities::ConstructKDTree(points, dim);
        
        // for every point find K-NN and create an edge between them
        int index = 0;
        for (vector<float> point : points) {
            ANNpoint query_point;
            ANNidxArray nn_ids;
            ANNdistArray dists;
            double eps = .01;
            query_point = annAllocPt(KD_tree->theDim(), 0.);
            nn_ids = new ANNidx[K + 1];
            dists = new ANNdist[K + 1];
            
            for (int d = 0; d < KD_tree->theDim(); d++) {
                query_point[d] = point[d];
            }
            
            KD_tree->annkSearch(query_point,
                                K + 1,
                                nn_ids,
                                dists,
                                eps);
            
            for (int i = 0; i < K + 1; i++) {
                if (index != nn_ids[i]) {
                    int min_id, max_id;
                    min_id = min(index, nn_ids[i]);
                    max_id = max(index, nn_ids[i]);
                    edges.insert(pair<int, int>(min_id, max_id));
                }
            }
            index++;
        }
    } else {
        dim = 0;
    }
}

KNN_Graph::KNN_Graph(int K_) {
    K = K_;
    vector<vector<float>> pts;
    points = pts;
    Construct();
}

KNN_Graph::KNN_Graph(string fp, int K_) {
    K = K_;
    vector<vector<float>> pts;
    Utilities::ReadInPoints(pts, fp);
    points = pts;
    Construct();
}

KNN_Graph::KNN_Graph(vector<vector<float>> pts, int K_) {
    K = K_;
    points = pts;
    Construct();
}

void KNN_Graph::writeToFile(string filepath) {
    ofstream out_file;
    out_file.open (filepath);
    out_file << dim << " " << points.size() << endl;
    for (vector<float> pt : points) {
        for (int i = 0; i < pt.size(); i++) {
            out_file << pt[i] << " ";
        }
        out_file << endl;
    }
    for (pair<int, int> edge : edges) {
        out_file << 2 << " " << edge.first << " " << edge.second << endl;
    }
    out_file.close();
}

// Writes a wieghted graph file based on Fangtao's file format
void KNN_Graph::writeWeightedToFile(string filepath) {
    ofstream out_file;
    out_file.open (filepath);
    // dim num_vert
    out_file << dim << " " << points.size() << endl;
    // x_1 x_2 .... x_d
    for (vector<float> pt : points) {
        for (int i = 0; i < dim; i++) {
            out_file << pt[i];
            if (i < dim - 1) {
                out_file << " ";
            }
        }
        out_file << endl;
    }
    // v_1 v_2 w_1
    // v_n1 v_n2 w_n1
    bool first = true;
    for (pair<int, int> p : edges) {
        if (!first) {
            out_file << endl;
        } else {
            first = !first;
        }
        vector<float> p1 = points[p.first];
        vector<float> p2 = points[p.second];
        double sum = 0.0;
        for (int i = 0; i < dim; i++) {
            sum += pow((p1[i] - p2[i]), 2.0);
        }
        out_file << p.first << " " << p.second << " " << sqrt(sum);
    }
    out_file.close();
}
