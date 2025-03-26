#ifndef CLUSTERING_C_KMEANS2D_H
#define CLUSTERING_C_KMEANS2D_H

#include "cg_util.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <ctime>

std::tuple<double, std::vector<int>, std::vector<Point>> kmeans2D(const std::vector<Point>& points, int k, int numTrials);
std::tuple<double, std::vector<int>, std::vector<Point>> kmeans2D(const std::vector<Point>& points, int k, int numTrials, std::string &instance_name);
std::tuple<double, std::vector<int>, std::vector<Point>> evaluate_clustering(const std::vector<Point> &points, int k, std::vector<int> &assignments);
std::vector<std::pair<int, double>> point_distance_from_centroids(Point &point, std::vector<Point> &centroids);
double computePartialSSE(int cluster_j, const std::vector<Point> &points, const std::vector<int> &assignments);

struct ClusteringResult {
    double loss; // loss of the produced cluster
    std::vector<int> assignment; // assignment vector
    std::vector<Point> centroid; // vector of k centroids
    std::vector<std::vector<iPoint>> points_cluster; // data points (with ids) in cluster j for each j = 1, ..., k
};

ClusteringResult get_cluster_assignment_2D(std::vector<Point> &P, int k, int n_start);
ClusteringResult *get_cluster_assignment_2D_ptr(std::vector<Point> &P, int k, int n_start);
ClusteringResult *get_cluster_assignment_2D_ptr(std::vector<Point> &P, int k, int n_start, std::string &instance_name);
ClusteringResult evaluate_cluster_assignment_2D(std::vector<Point> &P, int k, std::vector<int> &target_assignment);
ClusteringResult *evaluate_cluster_assignment_2D_ptr(std::vector<Point> &P, int k, std::vector<int> &target_assignment);

#endif //CLUSTERING_C_KMEANS2D_H
