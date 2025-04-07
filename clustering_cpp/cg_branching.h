#ifndef CLUSTERING_C_CG_BRANCHING_H
#define CLUSTERING_C_CG_BRANCHING_H

#include "cg_util.h"
#include "CGStruct.h"
#include "cg_pricing_multi.h"
#include "Kmeans2D.h"


struct PricingMINLP {
    std::vector<double> bigM;
    std::pair<double, double> y_coord_lb;
    std::pair<double, double> y_coord_ub;
};

struct DinkelbachResult {
    PricingOutput po;
    int n_iter;
    double time;
};

DinkelbachResult dinkelbach_pricing(std::vector<bool> &isVar, double f_init, std::vector<Point> &P, std::vector<std::vector<double>> &D, std::vector<double> &dual_c, double dual_k,
                                    std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs);

DinkelbachResult dinkelbach_pricing(double f_init, std::vector<Point> &P, std::vector<std::vector<double>> &D, std::vector<double> &dual_c, double dual_k, std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs);
//std::vector<int> convex_minlp(std::vector<int> &init_column, std::vector<Point> &P, std::vector<OrderedDistances> &od, std::vector<double> &dual_c, double dual_k, std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs);
PricingOutput pricing_heuristic(GRBModel *lp_model, PricingOutput &init_po, CGData *data, std::vector<double> &price_c, double price_k);
GRBModel *build_branching_problem_decomposed_ptr(std::vector<Point> &P, std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs);

std::vector<std::vector<int>> update_columns_must_link(std::vector<std::vector<int>> &or_cols, std::pair<int, int> &branching_pair);
std::vector<std::vector<int>> update_columns_cannot_link(std::vector<std::vector<int>> &or_cols, std::pair<int, int> &branching_pair);

#endif //CLUSTERING_C_CG_BRANCHING_H
