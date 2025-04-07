#ifndef CLUSTERING_C_CG_PRICING_MULTI_H
#define CLUSTERING_C_CG_PRICING_MULTI_H

#include <queue>
#include "cg_util.h"
#include "cg_pricing.h"
#include "cg_branching.h"

// root node
std::vector<std::vector<PricingOutput>> get_compatible_columns(MustLinkMapping &mlm, std::vector<OrderedDistances> &od, std::vector<double> &price_c, double price_k, std::vector<Point> &P, int max_comp_columns, int max_columns);
// must link and cannot link
std::vector<std::vector<PricingOutput>> get_compatible_columns(double &relaxed_rc, bool &exact_pricing, GRBModel *heuristic_model,
                                                               std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs,
                                                               CGData *data, MustLinkMapping &mlm, std::vector<double> &price_c, double price_k, int max_comp_columns,
                                                               bool heuristic_pricing, bool verbose);


std::vector<std::vector<PricingOutput>> get_compatible_columns(bool &exact_pricing, MustLinkMapping &mlm, std::vector<OrderedDistances> &od, std::vector<double> &price_c, double price_k, std::vector<Point> &P, int max_compatible_columns, int max_columns, bool verbose);
#endif //CLUSTERING_C_CG_PRICING_MULTI_H
