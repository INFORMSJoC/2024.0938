#include "cg_pricing.h"
#include "cg_util.h"
#include "cg_branching.h"

#define CAPACITY 1000000000

// computes objective function from vector of point indices (pricing problem)
double pricing_cost(std::vector<double> &price_c, double price_k, std::vector<Point> &P, std::vector<int> &a, Point &centroid){

    int a_size = (int) a.size();

    double sum_x = 0.0;
    double sum_y = 0.0;
    for (int l = 0; l < a_size; l++) {
        sum_x += P[a[l]].x;
        sum_y += P[a[l]].y;
    }
    centroid.x = sum_x / a_size;
    centroid.y = sum_y / a_size;

    double cost = 0;

    for (int l = 0; l < a_size; l++){
        cost += std::pow(point_distance(P[a[l]], centroid) , 2);
    }

    for (int l = 0; l < a_size; l++){
        cost -= price_c[a[l]];
    }

    return cost + price_k;

}

// enter = 1 if new_a is inserted in 'a', otherwise enter = -1
double change_pricing_cost(std::vector<double> &price_c, Point &new_a, int index_new_a, Point &centroid, double a_size, int enter){

    double cost;
    cost = enter * (a_size / (a_size + enter) * std::pow(point_distance(new_a, centroid), 2));
    cost = cost - enter * price_c[index_new_a];

    return cost;

}

// enter = 1 if new_a is inserted in 'a' otherwise enter = -1
void update_centroid(Point &new_a, Point &centroid, double a_size, int enter) {

    centroid.x = (a_size * centroid.x + enter * new_a.x) / (a_size + enter);
    centroid.y = (a_size * centroid.y + enter * new_a.y) / (a_size + enter);

}

// plain pricing problem (definitive version)
std::vector<PricingOutput> pricing_problem(MustLinkMapping &mlm, std::vector<OrderedDistances> &ordered_points, std::vector<double> &price_c, double price_k, std::vector<Point> &P) {

    std::vector<PricingOutput> columns; // negative reduced cost columns
    columns.reserve(CAPACITY);

    int n = (int) P.size();

    double current_cost;

    std::vector<int> cluster;
    cluster.reserve(n);
    std::vector<int> profitable;
    profitable.reserve(n);

    double max_price = 0;
    for (int i = 0; i < n; i++) {
        if (price_c[i] > 0) {
            profitable.push_back(i);
            if (price_c[i] > max_price)
                max_price = price_c[i];
        }
    }

    std::vector<bool> A(n); // set of discs which intersects another one in two points
    for (int i = 0; i < n; i++) {
        A[i] = false;
    }

    std::vector<std::vector<int>> B_isContained(n); // set of discs not in A
    for (int i = 0; i < n; i++) {
        B_isContained[i].push_back(i);
    }

    Point centroid = {0.0, 0.0};
    int n_intersect;
    double r_i, r_j;

    for (int i = 0; i < (int) (profitable.size() - 1); i++) {

        for (int j = 0; j < n - 1; j++) {

            // skips if the dual values of the variable is negative or if we have already considered this pair
            if (price_c[ordered_points[profitable[i]].index[j]] <= 0 ||
                profitable[i] > ordered_points[profitable[i]].index[j])
                continue;

            if (ordered_points[profitable[i]].dist[j] <=
                std::sqrt(price_c[profitable[i]]) + std::sqrt(price_c[ordered_points[profitable[i]].index[j]])) {

                Circle ci = {{P[profitable[i]].x, P[profitable[i]].y}, std::sqrt(price_c[profitable[i]])};
                Circle cj = {{P[ordered_points[profitable[i]].index[j]].x, P[ordered_points[profitable[i]].index[j]].y},
                             std::sqrt(price_c[ordered_points[profitable[i]].index[j]])};

                std::vector<Point> Q(2); // returns the intersection points
                n_intersect = intersect_circles(ci, cj, Q);

                if (n_intersect == 2) {

                    // put points i and j in A
                    A[profitable[i]] = true;
                    A[ordered_points[profitable[i]].index[j]] = true;

                    for (int q = 0; q < 2; q++) { // for each intersection point

                        cluster.clear();

                        for (int k = 0; k < n - 1; k++) {
                            double dist = ordered_points[profitable[i]].dist[k];
                            int exact_k = ordered_points[profitable[i]].index[k];

                            if (std::sqrt(price_c[profitable[i]]) + std::sqrt(max_price) < dist) break;
                            if (exact_k != ordered_points[profitable[i]].index[j] &&
                                point_distance(P[exact_k], Q[q]) <= std::sqrt(price_c[exact_k]))
                                cluster.push_back(exact_k);
                        }

                        // case where the only area to be considered is that with i and j
                        if (cluster.empty()) {
                            cluster.push_back(profitable[i]);
                            cluster.push_back(ordered_points[profitable[i]].index[j]);
                            current_cost = pricing_cost(price_c, price_k, P, cluster, centroid);
                            if (current_cost < -TOL)
                                columns.push_back({current_cost, cluster});
                        } else {  // i and j intersects other areas
                            current_cost = pricing_cost(price_c, price_k, P, cluster, centroid);
                            if (current_cost < -TOL)
                                columns.push_back({current_cost, cluster});
                            // evaluate only the difference of inserting i
                            current_cost += change_pricing_cost(price_c, P[profitable[i]], profitable[i], centroid, (double) cluster.size(), 1);
                            update_centroid(P[profitable[i]], centroid, (double) cluster.size(), 1);
                            cluster.push_back(profitable[i]);
                            if (current_cost < -TOL)
                                columns.push_back({current_cost, cluster});
                            // evaluate only the difference of inserting j
                            current_cost += change_pricing_cost(price_c, P[ordered_points[profitable[i]].index[j]], ordered_points[profitable[i]].index[j], centroid, (double) cluster.size(), 1);
                            update_centroid(P[ordered_points[profitable[i]].index[j]], centroid, (double) cluster.size(), 1);
                            cluster.push_back(ordered_points[profitable[i]].index[j]);
                            if (current_cost < -TOL)
                                columns.push_back({current_cost, cluster});
                            // evaluate only the difference of removing i
                            current_cost += change_pricing_cost(price_c, P[profitable[i]], profitable[i], centroid, (double) cluster.size(), -1);
                            // removes i and j and reinserts j
                            cluster.pop_back();
                            cluster.pop_back();
                            cluster.push_back(ordered_points[profitable[i]].index[j]);
                            if (current_cost < -TOL)
                                columns.push_back({current_cost, cluster});
                        }
                    }
                } else {

                    if (n_intersect == 0) { // one circle is contained into another

                        r_i = std::sqrt(price_c[profitable[i]]);
                        r_j = std::sqrt(price_c[ordered_points[profitable[i]].index[j]]);
                        if (r_i < r_j) {
                            B_isContained[profitable[i]].push_back(ordered_points[profitable[i]].index[j]);
                        } else {
                            B_isContained[ordered_points[profitable[i]].index[j]].push_back(profitable[i]);
                        }

                    }
                }
            } else {
                if (std::sqrt(price_c[profitable[i]]) + std::sqrt(max_price) <
                    ordered_points[profitable[i]].dist[j])
                    break;
            }

        }
    }

    // treat discs in B
    double cost_yolk;
    for (int i = 0; i < n; i++) {
        if (!A[i]) {
            cost_yolk = pricing_cost(price_c, price_k, P, B_isContained[i], centroid);
            if (cost_yolk < -TOL)
                columns.push_back({cost_yolk, B_isContained[i]});
        }
    }

    return columns;

}

// leave the pricing as soon as a maximum number of compatible columns with negative reduced costs are found
std::vector<std::vector<PricingOutput>> pricing_problem_limited(MustLinkMapping &mlm, std::vector<OrderedDistances> &ordered_points, std::vector<double> &price_c, double price_k, std::vector<Point> &P, int max_comp_cols) {

    std::vector<PricingOutput> columns(1); // negative reduced cost columns

    std::vector<PricingOutput> compatible_columns; // negative reduced cost compatible columns
    compatible_columns.reserve(max_comp_cols);

    double best_column_pricing_cost = 0.0;
    double best_compatible_column_pricing_cost = 0.0;
    int count_compatible = 0;

    int n = (int) P.size();

    double current_cost;

    std::vector<int> cluster;
    cluster.reserve(n);
    std::vector<int> profitable;
    profitable.reserve(n);

    double max_price = 0;
    for (int i = 0; i < n; i++) {
        if (price_c[i] > 0) {
            profitable.push_back(i);
            if (price_c[i] > max_price)
                max_price = price_c[i];
        }
    }

    std::vector<bool> A(n); // set of discs which intersects another one in two points
    for (int i = 0; i < n; i++) {
        A[i] = false;
    }

    std::vector<std::vector<int>> B_isContained(n); // set of discs not in A
    for (int i = 0; i < n; i++) {
        B_isContained[i].push_back(i);
    }

    Point centroid = {0.0, 0.0};
    int n_intersect;
    double r_i, r_j;

    bool exit = false;

    for (int i = 0; i < (int) (profitable.size() - 1); i++) {

        if (exit) break;

        for (int j = 0; j < n - 1; j++) {

            if (exit) break;

            // skips if the dual values of the variable is negative or if we have already considered this pair
            if (price_c[ordered_points[profitable[i]].index[j]] <= 0 ||
                profitable[i] > ordered_points[profitable[i]].index[j])
                continue;

            if (ordered_points[profitable[i]].dist[j] <=
                std::sqrt(price_c[profitable[i]]) + std::sqrt(price_c[ordered_points[profitable[i]].index[j]])) {

                Circle ci = {{P[profitable[i]].x, P[profitable[i]].y}, std::sqrt(price_c[profitable[i]])};
                Circle cj = {
                        {P[ordered_points[profitable[i]].index[j]].x, P[ordered_points[profitable[i]].index[j]].y},
                        std::sqrt(price_c[ordered_points[profitable[i]].index[j]])};

                std::vector<Point> Q(2); // returns the intersection points
                n_intersect = intersect_circles(ci, cj, Q);

                if (n_intersect == 2) {

                    // put points i and j in A
                    A[profitable[i]] = true;
                    A[ordered_points[profitable[i]].index[j]] = true;

                    for (int q = 0; q < 2; q++) { // for each intersection point

                        cluster.clear();

                        for (int k = 0; k < n - 1; k++) {
                            double dist = ordered_points[profitable[i]].dist[k];
                            int exact_k = ordered_points[profitable[i]].index[k];

                            if (std::sqrt(price_c[profitable[i]]) + std::sqrt(max_price) < dist) break;
                            if (exact_k != ordered_points[profitable[i]].index[j] &&
                                point_distance(P[exact_k], Q[q]) <= std::sqrt(price_c[exact_k]))
                                cluster.push_back(exact_k);
                        }

                        // case where the only area to be considered is that with i and j
                        if (cluster.empty()) {
                            cluster.push_back(profitable[i]);
                            cluster.push_back(ordered_points[profitable[i]].index[j]);
                            current_cost = pricing_cost(price_c, price_k, P, cluster, centroid);
                            if (current_cost < -TOL && current_cost < best_compatible_column_pricing_cost &&
                                count_compatible < max_comp_cols) {
                                if (isCompatible(mlm, cluster)) {
                                    compatible_columns.push_back({current_cost, cluster});
                                    best_compatible_column_pricing_cost = current_cost;
                                    count_compatible++;
                                }
                                if (count_compatible >= max_comp_cols) exit = true;
                            }
                            if (current_cost < -TOL && current_cost < best_column_pricing_cost) {
                                columns[0].pricing_cost = current_cost;
                                columns[0].column = cluster;
                                best_column_pricing_cost = current_cost;
                            }
                        } else {  // i and j intersects other areas
                            current_cost = pricing_cost(price_c, price_k, P, cluster, centroid);
                            if (current_cost < -TOL && current_cost < best_compatible_column_pricing_cost &&
                                count_compatible < max_comp_cols) {
                                if (isCompatible(mlm, cluster)) {
                                    compatible_columns.push_back({current_cost, cluster});
                                    best_compatible_column_pricing_cost = current_cost;
                                    count_compatible++;
                                }
                                if (count_compatible >= max_comp_cols) exit = true;
                            }
                            if (current_cost < -TOL && current_cost < best_column_pricing_cost) {
                                columns[0].pricing_cost = current_cost;
                                columns[0].column = cluster;
                                best_column_pricing_cost = current_cost;
                            }
                            // evaluate only the difference of inserting i
                            current_cost += change_pricing_cost(price_c, P[profitable[i]], profitable[i], centroid,
                                                                (double) cluster.size(), 1);
                            update_centroid(P[profitable[i]], centroid, (double) cluster.size(), 1);
                            cluster.push_back(profitable[i]);
                            if (current_cost < -TOL && current_cost < best_compatible_column_pricing_cost &&
                                count_compatible < max_comp_cols) {
                                if (isCompatible(mlm, cluster)) {
                                    compatible_columns.push_back({current_cost, cluster});
                                    best_compatible_column_pricing_cost = current_cost;
                                    count_compatible++;
                                }
                                if (count_compatible >= max_comp_cols) exit = true;
                            }
                            if (current_cost < -TOL && current_cost < best_column_pricing_cost) {
                                columns[0].pricing_cost = current_cost;
                                columns[0].column = cluster;
                                best_column_pricing_cost = current_cost;
                            }
                            // evaluate only the difference of inserting j
                            current_cost += change_pricing_cost(price_c, P[ordered_points[profitable[i]].index[j]],
                                                                ordered_points[profitable[i]].index[j], centroid,
                                                                (double) cluster.size(), 1);
                            update_centroid(P[ordered_points[profitable[i]].index[j]], centroid,
                                            (double) cluster.size(), 1);
                            cluster.push_back(ordered_points[profitable[i]].index[j]);
                            if (current_cost < -TOL && current_cost < best_compatible_column_pricing_cost &&
                                count_compatible < max_comp_cols) {
                                if (isCompatible(mlm, cluster)) {
                                    compatible_columns.push_back({current_cost, cluster});
                                    best_compatible_column_pricing_cost = current_cost;
                                    count_compatible++;
                                }
                                if (count_compatible >= max_comp_cols) exit = true;
                            }
                            if (current_cost < -TOL && current_cost < best_column_pricing_cost) {
                                columns[0].pricing_cost = current_cost;
                                columns[0].column = cluster;
                                best_column_pricing_cost = current_cost;
                            }
                            // evaluate only the difference of removing i
                            current_cost += change_pricing_cost(price_c, P[profitable[i]], profitable[i], centroid,
                                                                (double) cluster.size(), -1);
                            // removes i and j and reinserts j
                            cluster.pop_back();
                            cluster.pop_back();
                            cluster.push_back(ordered_points[profitable[i]].index[j]);
                            if (current_cost < -TOL && current_cost < best_compatible_column_pricing_cost &&
                                count_compatible < max_comp_cols) {
                                if (isCompatible(mlm, cluster)) {
                                    compatible_columns.push_back({current_cost, cluster});
                                    best_compatible_column_pricing_cost = current_cost;
                                    count_compatible++;
                                }
                                if (count_compatible >= max_comp_cols) exit = true;
                            }
                            if (current_cost < -TOL && current_cost < best_column_pricing_cost) {
                                columns[0].pricing_cost = current_cost;
                                columns[0].column = cluster;
                                best_column_pricing_cost = current_cost;
                            }
                        }
                    }
                } else {

                    if (n_intersect == 0) { // one circle is contained into another

                        r_i = std::sqrt(price_c[profitable[i]]);
                        r_j = std::sqrt(price_c[ordered_points[profitable[i]].index[j]]);
                        if (r_i < r_j) {
                            B_isContained[profitable[i]].push_back(ordered_points[profitable[i]].index[j]);
                        } else {
                            B_isContained[ordered_points[profitable[i]].index[j]].push_back(profitable[i]);
                        }

                    }
                }
            } else {
                if (std::sqrt(price_c[profitable[i]]) + std::sqrt(max_price) <
                    ordered_points[profitable[i]].dist[j])
                    break;
            }

        }
    }

    if (count_compatible < max_comp_cols) {
        // treat discs in B
        double cost_yolk;
        for (int i = 0; i < n; i++) {
            if (!A[i]) {
                cost_yolk = pricing_cost(price_c, price_k, P, B_isContained[i], centroid);
                if (cost_yolk < -TOL && cost_yolk < best_compatible_column_pricing_cost &&
                    count_compatible < max_comp_cols) {
                    if (isCompatible(mlm, B_isContained[i])) {
                        compatible_columns.push_back({cost_yolk, B_isContained[i]});
                        best_compatible_column_pricing_cost = cost_yolk;
                        count_compatible++;
                    }
                }
                if (cost_yolk < -TOL && cost_yolk < best_column_pricing_cost) {
                    columns[0].pricing_cost = cost_yolk;
                    columns[0].column = B_isContained[i];
                    best_column_pricing_cost = cost_yolk;
                }
            }
        }
    }


    return {compatible_columns, columns};

}

void updateBestCluster(std::vector<double> &price_c, double price_k, std::vector<Point> &P, std::vector<std::vector<double>> &D,
                       std::vector<int> &best_cluster, double &best_cost, std::vector<int> &a,
                       std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs) {

    Point centroid = {0.0, 0.0};

    double current_cost = pricing_cost(price_c, price_k, P, a, centroid);
    if (current_cost > best_cost)
        return;

    int n = (int) P.size();
    int a_size = (int) a.size();
    int v_count = 0;

    if (ml_pairs.empty() && cl_pairs.empty()) {
        current_cost = pricing_cost(price_c, price_k, P, a, centroid);
        if (current_cost < best_cost) {
            best_cost = current_cost;
            best_cluster.clear();
            for (int l = 0; l < a_size; l++) {
                best_cluster.push_back(a[l]);
            }
        }
    } else {

        std::vector<bool> isV(n, false);

        std::vector<bool> x_or(n, false);
        for (int i = 0; i < a_size; i++) {
            x_or[a[i]] = true;
        }

        // loop over cl constraints
        for (auto &pair: cl_pairs) {
            if (x_or[pair.first] && x_or[pair.second]) {
                if (!isV[pair.first]) {
                    isV[pair.first] = true;
                    v_count++;
                }
                if (!isV[pair.second]) {
                    isV[pair.second] = true;
                    v_count++;
                }
            }
        }

        // lover over ml constraints
        for (auto &pair: ml_pairs) {
            if (x_or[pair.first] != x_or[pair.second]) {
                if (!isV[pair.first]) {
                    isV[pair.first] = true;
                    v_count++;
                }
                if (!isV[pair.second]) {
                    isV[pair.second] = true;
                    v_count++;
                }
            }
        }

        if (v_count == 0) {
            // std::cout << "      No need to run exact pricing" << "\n";
            current_cost = pricing_cost(price_c, price_k, P, a, centroid);
            if (current_cost < best_cost) {
                best_cost = current_cost;
                best_cluster.clear();
                for (int l = 0; l < a_size; l++) {
                    best_cluster.push_back(a[l]);
                }
            }
        } else {

            std::vector<bool> isVar(n, false);
            int n_var = 0;
            for (int i = 0; i < n; i++) {
                if (x_or[i] || isV[i]) {
                    isVar[i] = true;
                    n_var++;
                }
            }

            // std::cout << "      Running exact pricing on a subset with " << n_var << " variables" << "\n";
            DinkelbachResult pr_exact = dinkelbach_pricing(isVar, 0.0, P, D, price_c, price_k, ml_pairs, cl_pairs);

            //std::cout << "      Iter: " << pr_exact.n_iter << "\t Time (s): " << pr_exact.time << "\n";
            //std::cout << "      Column: ";
            //print_vector_int(pr_exact.po.column);
            //std::cout << "      Reduced cost: " << pr_exact.po.pricing_cost << "\n";

            if (pr_exact.po.pricing_cost < best_cost) {
                best_cost = pr_exact.po.pricing_cost;
                best_cluster.clear();
                for (auto &point: pr_exact.po.column) {
                    best_cluster.push_back(point);
                }
            }
        }
    }
}

PricingOutput pricing_problem(MustLinkMapping &mlm, std::vector<OrderedDistances> &ordered_points, std::vector<double> &price_c, double price_k,
                              std::vector<Point> &P, std::vector<std::vector<double>> &D,
                              std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs) {

    int n = (int) P.size();

    std::vector<int> best_cluster;
    best_cluster.reserve(n);
    double best_cost = 0.0;

    std::vector<int> cluster;
    cluster.reserve(n);
    std::vector<int> profitable;
    profitable.reserve(n);

    double max_price = 0;
    for (int i = 0; i < n; i++) {
        if (price_c[i] > 0) {
            profitable.push_back(i);
            if (price_c[i] > max_price)
                max_price = price_c[i];
        }
    }

    std::vector<bool> A(n); // set of discs which intersects another one in two points
    for (int i = 0; i < n; i++) {
        A[i] = false;
    }

    std::vector<std::vector<int>> B_isContained(n); // set of discs not in A
    for (int i = 0; i < n; i++) {
        B_isContained[i].push_back(i);
    }

    int n_intersect;
    double r_i, r_j;

    for (int i = 0; i < (int) (profitable.size() - 1); i++) {

        for (int j = 0; j < n - 1; j++) {

            // skips if the dual values of the variable is negative or if we have already considered this pair
            if (price_c[ordered_points[profitable[i]].index[j]] <= 0 ||
                profitable[i] > ordered_points[profitable[i]].index[j])
                continue;

            if (ordered_points[profitable[i]].dist[j] <=
                std::sqrt(price_c[profitable[i]]) + std::sqrt(price_c[ordered_points[profitable[i]].index[j]])) {

                Circle ci = {{P[profitable[i]].x, P[profitable[i]].y}, std::sqrt(price_c[profitable[i]])};
                Circle cj = {{P[ordered_points[profitable[i]].index[j]].x, P[ordered_points[profitable[i]].index[j]].y},
                             std::sqrt(price_c[ordered_points[profitable[i]].index[j]])};

                std::vector<Point> Q(2); // returns the intersection points
                n_intersect = intersect_circles(ci, cj, Q);

                if (n_intersect == 2) {

                    // put points i and j in A
                    A[profitable[i]] = true;
                    A[ordered_points[profitable[i]].index[j]] = true;

                    for (int q = 0; q < 2; q++) { // for each intersection point

                        cluster.clear();

                        for (int k = 0; k < n - 1; k++) {
                            double dist = ordered_points[profitable[i]].dist[k];
                            int exact_k = ordered_points[profitable[i]].index[k];

                            if (std::sqrt(price_c[profitable[i]]) + std::sqrt(max_price) < dist) break;
                            if (exact_k != ordered_points[profitable[i]].index[j] && point_distance(P[exact_k], Q[q]) <= std::sqrt(price_c[exact_k]))
                                cluster.push_back(exact_k);
                        }

                        // case where the only area to be considered is that with i and j
                        if (cluster.empty()) {
                            cluster.push_back(profitable[i]);
                            cluster.push_back(ordered_points[profitable[i]].index[j]);
                            updateBestCluster(price_c, price_k, P, D, best_cluster, best_cost, cluster, ml_pairs, cl_pairs);
                        } else {  // i and j intersects other areas
                            updateBestCluster(price_c, price_k, P, D, best_cluster, best_cost, cluster, ml_pairs, cl_pairs);
                            // evaluate only the difference of inserting i
                            cluster.push_back(profitable[i]);
                            updateBestCluster(price_c, price_k, P, D, best_cluster, best_cost, cluster, ml_pairs, cl_pairs);
                            // evaluate only the difference of inserting j
                            cluster.push_back(ordered_points[profitable[i]].index[j]);
                            updateBestCluster(price_c, price_k, P, D, best_cluster, best_cost, cluster, ml_pairs, cl_pairs);
                            // removes i and j and reinserts j
                            cluster.pop_back();
                            cluster.pop_back();
                            cluster.push_back(ordered_points[profitable[i]].index[j]);
                            updateBestCluster(price_c, price_k, P, D, best_cluster, best_cost, cluster, ml_pairs, cl_pairs);
                        }
                    }
                } else {

                    if (n_intersect == 0) { // one circle is contained into another

                        r_i = std::sqrt(price_c[profitable[i]]);
                        r_j = std::sqrt(price_c[ordered_points[profitable[i]].index[j]]);
                        if (r_i < r_j) {
                            B_isContained[profitable[i]].push_back(ordered_points[profitable[i]].index[j]);
                        } else {
                            B_isContained[ordered_points[profitable[i]].index[j]].push_back(profitable[i]);
                        }

                    }
                }
            } else {
                if (std::sqrt(price_c[profitable[i]]) + std::sqrt(max_price) <
                    ordered_points[profitable[i]].dist[j])
                    break;
            }

        }
    }

    // treat discs in B
    for (int i = 0; i < n; i++) {
        if (!A[i]) {
            updateBestCluster(price_c, price_k, P, D, best_cluster, best_cost, B_isContained[i], ml_pairs, cl_pairs);
        }
    }

    return {best_cost, best_cluster};

}