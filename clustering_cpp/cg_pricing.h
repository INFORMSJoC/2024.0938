#ifndef CLUSTERING_C_CG_PRICING_H
#define CLUSTERING_C_CG_PRICING_H

#include <vector>
#include "cg_util.h"

struct PricingOutput
{
    double pricing_cost;
    std::vector<int> column;
    bool operator < (const PricingOutput &data) const {
        return pricing_cost < data.pricing_cost;
    }
};

struct PricingOutputComparator {
    bool operator()(const PricingOutput& lhs, const PricingOutput& rhs) const {
        return lhs.pricing_cost > rhs.pricing_cost;
    }
};


double pricing_cost(std::vector<double> &price_c, double price_k, std::vector<Point> &P, std::vector<int> &a, Point &centroid);
double change_pricing_cost(std::vector<double> &price_c, Point &new_a, int index_new_a, Point &centroid, double a_size, int enter);
void update_centroid(Point &new_a, Point &centroid, double a_size, int enter);

std::vector<PricingOutput> pricing_problem(MustLinkMapping &mlm, std::vector<OrderedDistances> &ordered_points, std::vector<double> &price_c, double price_k, std::vector<Point> &P);

PricingOutput pricing_problem(MustLinkMapping &mlm, std::vector<OrderedDistances> &ordered_points, std::vector<double> &price_c, double price_k, std::vector<Point> &P,
                              std::vector<std::vector<double>> &D, std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs);


std::vector<std::vector<PricingOutput>> pricing_problem_limited(MustLinkMapping &mlm, std::vector<OrderedDistances> &ordered_points, std::vector<double> &price_c, double price_k, std::vector<Point> &P, int max_comp_cols);
#endif //CLUSTERING_C_CG_PRICING_H
