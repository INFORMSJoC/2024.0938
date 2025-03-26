#ifndef CLUSTERING_C_CGSTRUCT_H
#define CLUSTERING_C_CGSTRUCT_H

#include <iostream>
#include <vector>
#include <gurobi_c++.h>
#include "cg_util.h"
#include "Kmeans2D.h"

struct CGData {

    std::vector<Point> P; // data points
    int k; // number of clusters
    std::vector<OrderedDistances> od; // sorted distances
    std::vector<std::vector<double>> D; // distance matrix
    std::vector<int> partition_assignment;
    double bigM; // cost of the cluster in which all points are clustered together

};

struct CGModel {

    MustLinkMapping mlm; // aggregating partition
    std::vector<std::vector<int>> or_cols; // inherited columns
    std::vector<std::pair<int, int>> ml_pairs; // must link constraints
    std::vector<std::pair<int, int>> cl_pairs; // cannot link constraints

};

struct CGConfig {

    // CG setting
    int n_components;
    double time_limit;
    int max_iter;
    int log_step;
    int dual_disaggregation_strategy;
    bool best_partition_update;
    double disaggregation_threshold;
    int max_comp_cols;
    int max_cols;
    bool cg_verbose;
    bool cg_heuristic_pricing_child;
    bool cg_mip_heuristic;
    int max_comp_cols_child;
    bool cg_inherit_partition_child;
    double ip_tol;

};

struct BBConfig {

    // BB setting
    std::string instance_name;
    double bb_tol;
    int bb_parallel;
    int bb_max_nodes;
    int bb_visiting_strategy;
    bool bb_verbose;

};

class BBResult {

public:

    double opt_value;
    std::vector<int> opt_assignment;
    double opt_gap;

};

class CGResult {

public:

    int n; // number of data points
    int k; // number of clusters
    int n_components; // number of components

    double start_ub; // initial loss
    double best_ub; // best loss found by the heuristic
    double best_lb; // best lower bound
    int n_covering_start; // initial number of covering constraints
    int n_covering_end; // final number of covering constraints
    int n_update; // number of times the partition has been updated
    int n_iter; // number of iterations
    double time;
    double time_pricing;
    double time_master;

    double avg_covering; // average number of covering constraints
    double avg_n_compatible; // average number of compatible columns with negative reduced cost
    double avg_n_incompatible; // average number of compatible columns with negative reduced cost
    double avg_target_degree; // best p-incompatible column;
    double avg_best_degree; // best incompatible column;
    double avg_time_pricing; // average time needed for solving the pricing problem
    double avg_time_master; // average time needed for solving the master problem

    int count_times_active; // count the number of times some bounds are active at the end of the CG

    // output structs
    GRBModel *grb_end; // gurobi restricted master problem
    CGModel *model_end; // mapping and columns at the end of CG algorithm
    ClusteringResult *cr_end; // final clustering result

    // optimal dual c
    // std::vector<double> opt_dual_c;
    // optimal dual c aggregated
    // std::vector<double> opt_dual_c_shr;

public:

    void print() const;
    std::string to_csv_line() const;
    static std::string to_csv_header();
    void save_assignment(const char *filename) const;
    ~CGResult();

};

class CGLog {

public:

    bool verbose;

    int iter; // iteration number
    int n_covering; // number of covering constraints
    int col_armp; // number of columns in ARMP
    double opt_armp; // optimal value of the ARMP
    double time_armp; // computational time in seconds for the ARMP
    int n_incompatible; // number of incompatible columns with negative reduced cost
    int n_compatible; // number of compatible columns with negative reduced cost
    double pricing_incompatible; // reduced cost of the best incompatible column
    double pricing_compatible; // reduced cost of the best compatible column
    double col_ratio; // ratio between the cost of the best compatible and incompatible
    double time_pricing; // computational time of the pricing procedure
    double lb_incomp; // valid lower bound (computed with the best incompatible column)
    double gap; // relative gap between opt_armp and lb_incomp
    double best_gap; // best relative gap
    double ip_gap; // relative gap between best 'lb_incomp' and best known ub
    double time_iter; // elapsed time

public:

    static void header();
    static void footer();
    void print() const;

};


#endif //CLUSTERING_C_CGSTRUCT_H
