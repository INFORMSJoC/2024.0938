#ifndef CLUSTERING_C_CG_REDUCTION_H
#define CLUSTERING_C_CG_REDUCTION_H

#include <gurobi_c++.h>
#include "cg_util.h"
#include "Kmeans2D.h"
#include "CGStruct.h"

#define AVERAGE 0
#define SPARSE 1
#define COMPLEMENTARITY 2

double cost_aggregation(std::vector<Point> &P, std::vector<bool> &a, MustLinkMapping &mlm);
CGResult *solve_master_problem_dca_child(CGData *data, CGModel *model, CGConfig *cg_config, BBConfig *bb_config, ClusteringResult *cr_start);



#endif
