#include "cg_pricing_multi.h"
#include "cg_pricing.h"
#include <queue>
#include <algorithm>


// RETURNS MULTIPLE NEGATIVE REDUCED COST COLUMNS (UNCONSTRAINED PRICING)
std::vector<std::vector<PricingOutput>> get_compatible_columns(MustLinkMapping &mlm, std::vector<OrderedDistances> &od, std::vector<double> &price_c, double price_k, std::vector<Point> &P, int max_compatible_columns, int max_columns) {

    std::vector<PricingOutput> columns = pricing_problem(mlm, od, price_c, price_k, P);
    // std::cout << columns.size() << "\n";

    std::sort(columns.begin(), columns.end()); // sorted negative reduced cost columns

    std::vector<PricingOutput> best_compatible;
    best_compatible.reserve(max_compatible_columns);
    std::vector<PricingOutput> best;
    best.reserve(max_columns);
    // if (!columns.empty()) best.push_back(columns[0]);

    double current_cost = 0.0;
    int count_compatible = 0;
    int count = 0;

    int cols = (int) columns.size();

    bool stop_compatible = false;
    bool stop = false;

    for (int i = 0; i < cols; i++) {

        if (std::abs(std::abs(current_cost) - std::abs(columns[i].pricing_cost)) > TOL) {
            current_cost = columns[i].pricing_cost;
            if (!stop_compatible) {
                if (isCompatible(mlm, columns[i].column)) {
                    best_compatible.push_back(columns[i]);
                    count_compatible++;
                }
            }
            if (!stop) {
                best.push_back(columns[i]);
                count++;
            }
        }

        if (count_compatible == max_compatible_columns)
            stop_compatible = true;

        if (count == max_columns)
            stop = true;

        if (stop_compatible && stop)
            break;

    }

    return {best_compatible, best};

}

std::vector<std::vector<PricingOutput>> get_compatible_columns(bool &exact_pricing, MustLinkMapping &mlm, std::vector<OrderedDistances> &od, std::vector<double> &price_c, double price_k, std::vector<Point> &P, int max_compatible_columns, int max_columns, bool verbose) {

    std::vector<std::vector<PricingOutput>> pp;

    // (compatible, incompatible)
    pp = pricing_problem_limited(mlm, od, price_c, price_k, P, max_compatible_columns);
    //std::cout << pp[0].size() << "\n";
    //std::cout << pp[1].size() << "\n";
    if (!pp[0].empty()) std::sort(pp[0].begin(), pp[0].end()); // sorted negative reduced cost columns

    // heuristic pricing executed, if there are no columns with negative reduced cost use exact pricing
    if (pp[1].empty()|| pp[1][0].pricing_cost > -TOL) {
        exact_pricing = true;
        if (verbose) std::cout << "      Running exact pricing..." << "\n";
        pp = get_compatible_columns(mlm, od, price_c, price_k, P, max_compatible_columns, max_columns);
        return pp;
    }

    return pp;

}


// RETURNS MULTIPLE NEGATIVE REDUCED COST COLUMNS (CONSTRAINED PRICING) AND PERFORMS EXACT PRICING WHEN REQUIRED
std::vector<std::vector<PricingOutput>> get_compatible_columns(double &relaxed_rc, bool &exact_pricing, GRBModel *heuristic_model,
                                                                std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs,
                                                                CGData *data, MustLinkMapping &mlm, std::vector<double> &price_c, double price_k, int max_comp_columns,
                                                                bool heuristic_pricing, bool verbose) {

    int n = (int) mlm.ml_index.size();

    if (!heuristic_pricing) {

        std::vector<PricingOutput> best_compatible;
        best_compatible.reserve(1);
        std::vector<PricingOutput> best;
        best.reserve(1);

        exact_pricing = true;
        if (verbose)
            std::cout << "      Running exact pricing..." << "\n";
        PricingOutput exact = pricing_problem(mlm, data->od, price_c, price_k, data->P, data->D, ml_pairs, cl_pairs);
        if (verbose)
            std::cout << "      Exact pricing reduced cost: " << exact.pricing_cost << "\n";
        if (exact.pricing_cost < -TOL) {
            best.push_back(exact);
            if (isCompatible(mlm, exact.column)) {
                best_compatible.push_back(exact);
            }
        }

        return {best_compatible, best};
    }

    std::vector<PricingOutput> columns = pricing_problem(mlm, data->od, price_c, price_k, data->P);

    std::sort(columns.begin(), columns.end()); // sorted negative reduced cost columns

    if (!columns.empty()) {
        // minimum reduced cost from the unconstrained (relaxed) pricing (used to compute a valid lower bound at the current CG iteration)
        relaxed_rc = columns[0].pricing_cost;
    }

    std::vector<PricingOutput> best_compatible;
    best_compatible.reserve(max_comp_columns);
    std::vector<PricingOutput> best;
    best.reserve(1);

    PricingOutput best_column;
    best_column.pricing_cost = 0.0;
    best_column.column.reserve(n);

    double current_cost = 0.0;
    int count_compatible = 0;
    int cols = (int) columns.size();

    for (int i = 0; i < cols; i++) {

        int v_count = 0;

        if (std::abs(std::abs(current_cost) - std::abs(columns[i].pricing_cost)) > TOL) {
            current_cost = columns[i].pricing_cost;

            std::vector<bool> x_or(n, false);
            for (int j : columns[i].column) {
                x_or[j] = true;
            }

            //std::cout << "CL: " << "\n";
            // loop over cl constraints
            for (auto &pair: cl_pairs) {
                //std::cout << "(" << pair.first << " " << pair.second << ") ";
                if (x_or[pair.first] && x_or[pair.second]) {
                    v_count++;
                }
            }

            //std::cout << "\nML: " << "\n";
            // lover over ml constraints
            for (auto &pair: ml_pairs) {
                //std::cout << "(" << pair.first << " " << pair.second << ") ";
                if (x_or[pair.first] != x_or[pair.second]) {
                    v_count++;
                }
            }

            // if the unconstrained column is not involved in any constraint there no need to call the heuristic
            PricingOutput po;
            if (v_count == 0) {
                po = columns[i];
            } else {
                po = pricing_heuristic(heuristic_model, columns[i], data, price_c, price_k);
            }

            // PricingOutput po = pricing_heuristic(heuristic_model, columns[i], data, price_c, price_k);

            if (po.pricing_cost < -TOL) {

                if (po.pricing_cost < best_column.pricing_cost)
                    best_column = po;

                if (isCompatible(mlm, po.column)) {
                    best_compatible.push_back(po);
                    count_compatible++;
                }

            }
        }


        if (count_compatible == max_comp_columns)
            break;

    }

    best.push_back(best_column);

    std::sort(best_compatible.begin(), best_compatible.end()); // sorted negative reduced cost columns

    auto it_best_compatible = std::unique(best_compatible.begin(), best_compatible.end(),
                                          [](PricingOutput &l, PricingOutput &r) { return std::abs(std::abs(l.pricing_cost) - std::abs(r.pricing_cost)) < TOL; });
    best_compatible.resize(distance(best_compatible.begin(), it_best_compatible));

    // heuristic pricing executed, if there are no columns with negative reduced cost use exact pricing
    if (best.empty() || best_column.pricing_cost > -TOL) {
        exact_pricing = true;
        if (verbose)
            std::cout << "      Running exact pricing..." << "\n";
        PricingOutput exact = pricing_problem(mlm, data->od, price_c, price_k, data->P, data->D, ml_pairs, cl_pairs);
        if (verbose)
            std::cout << "      Exact pricing reduced cost: " << exact.pricing_cost << "\n";
        if (exact.pricing_cost < -TOL) {
            best[0] = exact;
            if (isCompatible(mlm, exact.column)) {
                best_compatible.push_back(exact);
            }
        }
    }

    return {best_compatible, best};
}


/*
// RETURNS MULTIPLE NEGATIVE REDUCED COST COLUMNS (CONSTRAINED PRICING) AND PERFORMS EXACT PRICING WHEN NEEDED
std::vector<std::vector<PricingOutput>> get_compatible_columns(double &relaxed_rc, bool &exact_pricing, GRBModel *heuristic_model,
                                                               std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs,
                                                               CGData *data, MustLinkMapping &mlm, std::vector<double> &price_c, double price_k, int max_comp_columns, int max_columns) {

    int n = (int) mlm.ml_index.size();

    std::vector<PricingOutput> columns = pricing_problem(mlm, data->od, price_c, price_k, data->P);

    std::sort(columns.begin(), columns.end()); // sorted negative reduced cost columns

    if (!columns.empty()) {
        // minimum reduced cost from the unconstrained (relaxed) pricing (used to compute a valid lower bound at the current CG iteration)
        relaxed_rc = columns[0].pricing_cost;
    }

    std::vector<PricingOutput> best_compatible;
    best_compatible.reserve(max_comp_columns);
    std::vector<PricingOutput> best;
    best.reserve(1);

    PricingOutput best_column;
    best_column.pricing_cost = 0.0;
    best_column.column.reserve(n);

    double current_cost = 0.0;
    int count_compatible = 0;

    int cols = (int) columns.size();

    for (int i = 0; i < cols; i++) {

        if (std::abs(std::abs(current_cost) - std::abs(columns[i].pricing_cost)) > TOL) {
            current_cost = columns[i].pricing_cost;

            //std::cout << "Examining: ";
            //std::sort(columns[i].column.begin(), columns[i].column.end());
            //print_vector_int(columns[i].column);
            PricingOutput po = pricing_heuristic(heuristic_model, columns[i], data, price_c, price_k);
            //std::cout << "Returned: ";
            //print_vector_int(po.column);
            //std::cout << "Cost: " << po.pricing_cost << "\n";

            if (po.pricing_cost < -TOL) {

                if (po.pricing_cost < best_column.pricing_cost)
                    best_column = po;

                if (isCompatible(mlm, po.column)) {
                    best_compatible.push_back(po);
                    count_compatible++;
                }

            }
        }


        if (count_compatible == max_comp_columns)
            break;

    }

    best.push_back(best_column);

    std::sort(best_compatible.begin(), best_compatible.end()); // sorted negative reduced cost columns

    auto it_best_compatible = std::unique(best_compatible.begin(), best_compatible.end(),
                                          [](PricingOutput &l, PricingOutput &r) { return std::abs(std::abs(l.pricing_cost) - std::abs(r.pricing_cost)) < TOL; });
    best_compatible.resize(distance(best_compatible.begin(), it_best_compatible));

    // heuristic pricing executed, if there are no columns with negative reduced cost use exact pricing
    if (best.empty() || best_column.pricing_cost > -TOL) {
        exact_pricing = true;
        std::cout << "      Running exact pricing..." << "\n";
        // std::cout << "      Init pricing: " << best_column.pricing_cost - price_k << "\n";
        DinkelbachResult pr_exact = dinkelbach_pricing(0.0, data->P, data->D, price_c, price_k, ml_pairs, cl_pairs);
        std::cout << "      Iter: " << pr_exact.n_iter << "\t Time (s): " << pr_exact.time << "\n";
        std::cout << "      Reduced cost: " << pr_exact.po.pricing_cost << "\n";
        if (pr_exact.po.pricing_cost < -TOL) {
            best[0] = pr_exact.po;
            if (isCompatible(mlm, pr_exact.po.column)) {
                best_compatible.push_back(pr_exact.po);
            }
        }
    }

    return {best_compatible, best};
}
*/