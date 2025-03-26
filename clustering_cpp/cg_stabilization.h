#ifndef CLUSTERING_C_CG_STABILIZATION_H
#define CLUSTERING_C_CG_STABILIZATION_H

#include "cg_util.h"
#include "Kmeans2D.h"
#include "CGStruct.h"

struct DualBounds {
    std::vector<double> lb_box;
    std::vector<double> ub_box;
};

CGResult *solve_stabilized_master_problem_dca(CGData *data, CGConfig *config, ClusteringResult *cr_start);
// CGResult *solve_stabilized_master_problem_dca_child(CGData *data, CGModel *model, CGConfig *config, ClusteringResult *cr_start);

CGResult *solve_stabilized_master_problem_dca_heuristic(CGData *data, CGConfig *config, ClusteringResult *cr_start);

GRBModel *build_stabilized_restricted_master_problem_aggregation_ptr(std::vector<Point> &P, std::vector<std::vector<bool>> &A, int k, MustLinkMapping &mlm, std::vector<double> &lb, std::vector<double> &ub);
void update_stabilized_restricted_master_problem(GRBModel *master, std::vector<Point> &P, MustLinkMapping &mlm, std::vector<std::vector<bool>> A, std::vector<double> &lb, std::vector<double> &ub);
DualBounds compute_bounds(MustLinkMapping &mlm, std::vector<Point> &P, std::vector<int> &opt_assignment, std::vector<Point> &opt_centroids);
bool update_bounds(GRBModel *rms, std::vector<double> &lb, std::vector<double> &ub, int &active_lb, int &active_ub);

#endif //CLUSTERING_C_CG_STABILIZATION_H
