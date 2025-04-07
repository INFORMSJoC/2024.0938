#ifndef CLUSTERING_C_CG_UTIL_H
#define CLUSTERING_C_CG_UTIL_H

#include <gurobi_c++.h>
#include <unordered_set>
#include <iomanip>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <set>

#define TOL 1e-5

struct MustLinkMapping
{
    std::vector<std::unordered_set<int>> ml_comp;
    std::vector<int> ml_index;
};

struct Point
{
    double x;
    double y;
};

struct iPoint
{
    int id;
    Point point;
};

struct Circle
{
    Point center;
    double radius;
};

struct OrderedDistances {
    std::vector<double> dist;
    std::vector<int> index;
};

std::vector<OrderedDistances> get_ordered_distances(std::vector<Point>& P);
std::vector<std::vector<double>> get_distance_matrix(std::vector<Point> &P);
double get_cost_clustered_together(std::vector<Point> &P);
double point_distance(Point &p1, Point &p2);
int intersect_circles(Circle &circle1, Circle &circle2, std::vector<Point> &Q);

std::vector<double> disaggregate_dual_variables_uniform(MustLinkMapping &mlm, std::vector<double> &dual_c_shr);
std::vector<double> disaggregate_dual_variables_sparse(MustLinkMapping &mlm, std::vector<double> &dual_c_shr, int index);
std::vector<double> disaggregate_dual_variables_complementarity_problem(std::vector<Point> &P, MustLinkMapping &mlm,
                                                                        std::vector<std::vector<int>> &incompatible_cols,
                                                                        double aggr_dual_k, std::vector<double> &aggr_dual_c);

std::vector<std::pair<int, int>> candidate_must_link_radius(int k, double perc, std::vector<Point> &centroids, std::vector<std::vector<iPoint>> &points);

bool isCompatible(MustLinkMapping &mlm, std::vector<int> &column);
std::vector<bool> aggregate_column(MustLinkMapping &mlm, std::vector<int> &column);
std::vector<std::vector<bool>> update_columns(MustLinkMapping &mlm, std::vector<std::vector<int>> &cols);
std::vector<std::vector<bool>> get_initial_columns_aggregation(std::vector<std::vector<iPoint>> &points, MustLinkMapping &mlm, std::vector<std::vector<int>> &or_cols);
int get_number_of_incompatibilities(MustLinkMapping &mlm, std::vector<int> &column);
void update_aggregation_improved(MustLinkMapping &mlm, std::vector<int> &column);

MustLinkMapping get_ml_map(int n, std::vector<std::pair<int, int>> &ml);
MustLinkMapping get_ml_map(std::vector<int> &assignment);
MustLinkMapping get_ml_map(std::vector<int> &hierarchical_assignment, std::vector<int> &kmeans_assignment);
void display_graph(MustLinkMapping &mlm);



#endif //CLUSTERING_C_CG_UTIL_H
