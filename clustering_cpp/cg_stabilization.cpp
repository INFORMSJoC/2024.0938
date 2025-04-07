#include <unordered_set>
#include <chrono>
#include <fstream>
#include <random>
#include "cg_stabilization.h"
#include "cg_reduction.h"
#include "cg_util.h"
#include "cg_pricing.h"
#include "cg_pricing_multi.h"
#include "CGStruct.h"
#include "Kmeans2D.h"


GRBModel *build_stabilized_restricted_master_problem_aggregation_ptr(std::vector<Point> &P, std::vector<std::vector<bool>> &A, int k, MustLinkMapping &mlm, std::vector<double> &lb, std::vector<double> &ub) {

    // P: n x d matrix of data points
    // A: n_aggregated x n_col matrix of initial columns
    // k: number of clusters
    // ml_map: mapping between points ad aggregated points
    // lb: lower bound for dual variables
    // ub: upper bound for dual variables

    int n_cols = (int) A.size();
    // number of super-points
    int n_aggregated = (int) mlm.ml_comp.size();

    auto env = GRBEnv(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    auto master = new GRBModel(env);

    // Variables of master problem
    std::vector<GRBVar> VarMaster(n_cols);
    for (int i = 0; i < n_cols; i++){
        VarMaster[i] = master->addVar(0.0, GRB_INFINITY, cost_aggregation(P, A[i], mlm), GRB_CONTINUOUS, "C" + std::to_string(i));
    }
    std::vector<GRBVar> xi(n_aggregated);
    std::vector<GRBVar> eta(n_aggregated);
    for (int i = 0; i < n_aggregated; i++) {
        xi[i] = master->addVar(0.0, GRB_INFINITY, -lb[i], GRB_CONTINUOUS, "XI" + std::to_string(i));
        eta[i] = master->addVar(0.0, GRB_INFINITY, ub[i], GRB_CONTINUOUS, "ETA" + std::to_string(i));
    }

    // Constraints of the master problem
    GRBLinExpr sum_k = 0;
    for (int i = 0; i < n_cols; i++)
        sum_k += VarMaster[i];
    master->addConstr(-sum_k >= -k); // we want k clusters

    GRBLinExpr lhs_sum = 0;
    for (int i = 0; i < n_aggregated; i++) {
        lhs_sum = 0;
        for (int j = 0; j < n_cols; j++) {
            if (A[j][i])
                lhs_sum += VarMaster[j];
        }
        master->addConstr(- xi[i] + eta[i] + lhs_sum >= 1); // each data point must be assigned to a single cluster
    }

    master->update();

    master->set("OutputFlag", "0");
    master->set("Method", "1"); // -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier

    return master;

}

void update_stabilized_restricted_master_problem(GRBModel *master, std::vector<Point> &P, MustLinkMapping &mlm, std::vector<std::vector<bool>> A, std::vector<double> &lb, std::vector<double> &ub) {

    int old_n_covering = (int) master->get(GRB_IntAttr_NumConstrs) - 1;
    int new_n_covering = (int) mlm.ml_comp.size();
    int old_n_cols = (int) master->get(GRB_IntAttr_NumVars) - 2 * old_n_covering;
    int new_n_cols = (int) A.size();

    // update old xi and eta
    for (int i = 0; i < old_n_covering; i++) {
        master->getVarByName("XI" + std::to_string(i)).set(GRB_DoubleAttr_Obj, -lb[i]);
        master->getVarByName("ETA" + std::to_string(i)).set(GRB_DoubleAttr_Obj, ub[i]);
    }

    // add covering constraints
    GRBLinExpr lhs_sum = 0;
    for (int i = old_n_covering; i < new_n_covering; i++) {
        lhs_sum = 0;
        for (int j = 0; j < old_n_cols; j++) {
            if (A[j][i]) {
                lhs_sum += master->getVarByName("C" + std::to_string(j));
            }
        }
        GRBVar xi = master->addVar(0.0, GRB_INFINITY, -lb[i], GRB_CONTINUOUS, "XI" + std::to_string(i));
        GRBVar eta = master->addVar(0.0, GRB_INFINITY, ub[i], GRB_CONTINUOUS, "ETA" + std::to_string(i));
        master->addConstr(- xi + eta + lhs_sum >= 1);
    }

    master->update();

    for (int i = old_n_cols; i < new_n_cols; i++) {
        GRBConstr *tempCons = master->getConstrs();
        GRBColumn col;
        col.addTerm(-1.0, tempCons[0]);
        for (int j = 1; j < new_n_covering + 1; j++)
            col.addTerm(int(A[i][j - 1]), tempCons[j]);
        master->addVar(0.0, GRB_INFINITY, cost_aggregation(P, A[i], mlm), GRB_CONTINUOUS, col, "C" + std::to_string(i));
    }

}


DualBounds compute_bounds(MustLinkMapping &mlm, std::vector<Point> &P, std::vector<int> &opt_assignment, std::vector<Point> &opt_centroids) {

    int n_comp = (int) mlm.ml_comp.size();

    std::vector<double> lb(n_comp);
    std::vector<double> ub(n_comp);

    int j_first, j_second;
    double c_bar_i, c_tilde_i;
    std::vector<int> tmp_assignment;
    std::vector<std::pair<int, double>> pair;

    for (int i = 0; i < n_comp; i++) {

        // find the cluster where all the points of the component are assigned (use their centroid to get the cluster index)
        Point super_point = {0.0, 0.0};
        for (auto &p : mlm.ml_comp[i]) {
            super_point.x = super_point.x + P[p].x;
            super_point.y = super_point.y + P[p].y;
        }
        super_point.x = super_point.x / (int) mlm.ml_comp[i].size();
        super_point.y = super_point.y / (int) mlm.ml_comp[i].size();
        pair = point_distance_from_centroids(super_point, opt_centroids); // sorted vector or pairs (cluster_id, distance_from_cluster_id)

        /*********************************** LOWER BOUND ESTIMATION ********************************/

        // the index of the closest cluster to the super point
        j_first = pair[0].first;
        // compute the cost of the clustering 'j' where the super point of the component 'i' is assigned
        c_bar_i = computePartialSSE(j_first, P, opt_assignment);
        tmp_assignment = opt_assignment;
        for (auto &p : mlm.ml_comp[i])
            tmp_assignment[p] = -1;
        // compute the cost of the clustering 'j' where we omit all the points in component 'i'
        c_tilde_i = computePartialSSE(j_first, P, tmp_assignment);
        lb[i] = c_bar_i - c_tilde_i;
        lb[i] = lb[i] - lb[i] * 0.0001;

        /********************************** UPPER BOUND ESTIMATION *********************************/

        // the index of the second-best cluster for the super point
        j_second = pair[1].first;
        // compute the cost of the clustering 'j_second'
        c_bar_i = computePartialSSE(j_second, P, opt_assignment);
        tmp_assignment = opt_assignment;
        for (auto &p : mlm.ml_comp[i])
            tmp_assignment[p] = j_second;
        // compute the cost of the clustering 'j' where we moved all the points of component 'i'
        c_tilde_i = computePartialSSE(j_second, P, tmp_assignment);
        ub[i] = c_tilde_i - c_bar_i;
        ub[i] = ub[i] + ub[i] * 0.0001;

        if (ub[i] < lb[i]) lb[i] = 0.0; // avoid infeasibility when the quality of the initial clustering is not good

    }

    return {lb, ub};
}

// update estimate lb and ub on lambda variables
bool update_bounds(GRBModel *rms, std::vector<double> &lb, std::vector<double> &ub, int &active_lb, int &active_ub) {

    bool updated = false;
    int n_comp = (int) lb.size();
    double alpha;

    for (int i = 0; i < n_comp; i++) {
        alpha = ub[i] - lb[i];
        if (std::abs(rms->getConstr(i+1).get(GRB_DoubleAttr_Pi) - lb[i]) < TOL) {
            // update the lower bound
            if (lb[i] > 0) {
                active_lb++;
                lb[i] = std::max(0.0, lb[i] - alpha);
                rms->getVarByName("XI" + std::to_string(i)).set(GRB_DoubleAttr_Obj, -lb[i]);
                updated = true;
            }
        }
        if (std::abs(rms->getConstr(i+1).get(GRB_DoubleAttr_Pi) - ub[i]) < TOL) {
            // update the upper bound
            if (ub[i] > 0) {
                // std::cout << "OLD UB[i]: " << ub[i] << "\t ALPHA[i]: " << alpha;
                ub[i] = ub[i] + alpha;
                // std::cout << "\t NEW UB[i]: " << ub[i] << "\n";
                active_ub++;
                rms->getVarByName("ETA" + std::to_string(i)).set(GRB_DoubleAttr_Obj, ub[i]);
                updated = true;
            }
        }
    }

    return updated;
}

std::vector<int> cluster_recovery_stabilized(int n, std::vector<std::vector<int>> &or_cols, GRBModel *rms, int k) {

    std::vector<int> assignment(n);
    int n_covering = (int) rms->get(GRB_IntAttr_NumConstrs) - 1;
    int n_cols = (int) rms->get(GRB_IntAttr_NumVars) - 2 * n_covering;
    int cluster_id = 0;
    for (int i = 0; i < n_cols; i++) {
        double v = rms->getVarByName("C" + std::to_string(i)).get(GRB_DoubleAttr_X);
        if (v <= 1.0 + 0.01 && v >= 1.0 - 0.01) {
            for (auto &elem: or_cols[i])
                assignment[elem] = cluster_id;
            cluster_id++;
        }
    }

    if (cluster_id == k) {
        return assignment;
    } else {
        return {};
    }
}

GRBModel *mip_post_processing_stabilized(GRBModel *rms, int &bin_vars) {

    auto *rms_discrete = new GRBModel(*rms);
    rms_discrete->optimize();

    int n_covering = (int) rms->get(GRB_IntAttr_NumConstrs) - 1;
    int n_cols = (int) rms->get(GRB_IntAttr_NumVars) - 2 * n_covering;
    bool is_integer = true;
    for (int i = 0; i < n_cols; i++) {
        double v = rms_discrete->getVarByName("C" + std::to_string(i)).get(GRB_DoubleAttr_X);
        if (v > 0.01 && v < 1 - 0.01) {
            bin_vars++;
            rms_discrete->getVarByName("C" + std::to_string(i)).set(GRB_CharAttr_VType, GRB_BINARY);
            is_integer = false;
        }
    }
    if (!is_integer) {
        rms_discrete->set("OutputFlag", "0");
        rms_discrete->optimize();
        // std::cout << "MIP Objective: " << rms_discrete.get(GRB_DoubleAttr_ObjVal) << "\n";
        return rms_discrete;
    }

    delete (rms_discrete);
    return nullptr;

}

// MAIN CG
CGResult *solve_stabilized_master_problem_dca(CGData *data, CGConfig *config, ClusteringResult *cr_start) {

    // CG result
    auto *result = new CGResult();
    result->n = (int) data->P.size();
    result->k = data->k;
    result->n_components = config->n_components;

    result->start_ub = cr_start->loss;
    result->best_ub = cr_start->loss;

    auto start = std::chrono::high_resolution_clock::now();

    MustLinkMapping mlm = get_ml_map(data->partition_assignment, cr_start->assignment);
    // display_graph(mlm);

    DualBounds bounds = compute_bounds(mlm, data->P, cr_start->assignment, cr_start->centroid);

    // CG log
    CGLog log{};
    log.verbose = config->cg_verbose;

    // Initial number of covering constraints
    log.n_covering = (int) mlm.ml_comp.size();
    result->n_covering_start = log.n_covering;

    std::vector<std::vector<bool>> init_A;
    // keep track of the added columns at each iteration
    std::vector<std::vector<int>> or_cols;
    or_cols.reserve(config->max_iter);
    // add the initial set of columns
    init_A = get_initial_columns_aggregation(cr_start->points_cluster, mlm, or_cols);
    // std::cout << "Initial columns: " << or_cols.size() << "\n";

    GRBModel *rms = build_stabilized_restricted_master_problem_aggregation_ptr(data->P, init_A, data->k, mlm, bounds.lb_box, bounds.ub_box);

    if (config->cg_verbose) CGLog::header();

    log.n_compatible = 0;
    log.n_incompatible = 0;
    log.col_ratio = 0.0;
    log.opt_armp = 0.0;
    log.col_armp = 0;
    log.iter = 0;
    log.gap = std::numeric_limits<double>::infinity();
    log.best_gap = std::numeric_limits<double>::infinity();
    log.ip_gap = std::numeric_limits<double>::infinity();

    result->time_master = 0.0;
    result->time_pricing = 0.0;
    result->n_update = 0;
    result->count_times_active = 0;
    result->best_lb = -std::numeric_limits<double>::infinity();

    std::vector<std::vector<int>> incompatible_cols;
    incompatible_cols.reserve(config->max_iter);

    int sum_best_degree = 0; // best incompatible column;
    int sum_target_degree = 0; // best p-incompatible column;
    int sum_covering = 0;
    int count_n_incompatible = 0;
    int count_n_compatible = 0;

    double partial_time = 0.0;

    while (true) {

        auto start_iter = std::chrono::high_resolution_clock::now();

        sum_covering += log.n_covering;
        log.col_armp = (int) or_cols.size();
        auto start_master = std::chrono::high_resolution_clock::now();
        rms->optimize();

        auto stop_master = std::chrono::high_resolution_clock::now();
        auto elapsed_master = (std::chrono::duration_cast<std::chrono::microseconds>(stop_master - start_master).count());
        log.time_armp = (double) elapsed_master / 1000000;
        log.opt_armp = rms->get(GRB_DoubleAttr_ObjVal);
        result->time_master += log.time_armp;

        //Dual variables associated to point constraints
        double dual_k = rms->getConstr(0).get(GRB_DoubleAttr_Pi);
        std::vector<double> dual_c_shr(log.n_covering);
        for (int j = 1; j < log.n_covering + 1; j++)
            dual_c_shr[j - 1] = rms->getConstr(j).get(GRB_DoubleAttr_Pi);

        // disaggregate dual variables before calling the pricing problem
        std::vector<double> dual_c;
        switch (config->dual_disaggregation_strategy) {
            case AVERAGE:
                dual_c = disaggregate_dual_variables_uniform(mlm, dual_c_shr);
                break;
            case SPARSE:
                dual_c = disaggregate_dual_variables_sparse(mlm, dual_c_shr, 0);
                break;
            case COMPLEMENTARITY:
                if (incompatible_cols.empty()) {// if incompatible columns have not been generated then disaggregate with the mean
                    dual_c = disaggregate_dual_variables_uniform(mlm, dual_c_shr);
                } else { // solve the complementary problem
                    dual_c = disaggregate_dual_variables_complementarity_problem(data->P, mlm, incompatible_cols,
                                                                                 dual_k, dual_c_shr);
                    break;
                    default:
                        dual_c = disaggregate_dual_variables_uniform(mlm, dual_c_shr);
                    break;
                }
        }

        // (compatible, incompatible)
        auto start_pricing = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<PricingOutput>> p_result;
        // unconstrained pricing
        p_result = get_compatible_columns(mlm, data->od, dual_c, dual_k, data->P, config->max_comp_cols, config->max_cols);
        auto stop_pricing = std::chrono::high_resolution_clock::now();
        auto elapsed_pricing = (std::chrono::duration_cast<std::chrono::microseconds>(stop_pricing - start_pricing).count());
        log.time_pricing = (double) elapsed_pricing / 1000000;

        if (!p_result[1].empty()) {
            log.n_incompatible = (int) p_result[1].size();
            count_n_compatible += log.n_compatible;
            log.pricing_incompatible = p_result[1][0].pricing_cost;
            if (!p_result[0].empty()) {
                log.n_compatible = (int) p_result[0].size();
                count_n_incompatible += log.n_incompatible;
                log.pricing_compatible = p_result[0][0].pricing_cost;
                log.col_ratio = p_result[0][0].pricing_cost / p_result[1][0].pricing_cost;
            } else {
                log.n_compatible = 0;
                log.pricing_compatible = 0.0;
                log.col_ratio = 0.0;
            }
        } else {
            log.n_incompatible = 0;
            log.n_compatible = 0;
            log.pricing_incompatible = 0.0;
            log.pricing_compatible = 0.0;
            log.col_ratio = 0.0;
        }
        result->time_pricing += log.time_pricing;

        if (config->dual_disaggregation_strategy == COMPLEMENTARITY) {
            for (auto &elem : p_result[1]) {
                if (!isCompatible(mlm, elem.column) && elem.pricing_cost < 0) {
                    incompatible_cols.push_back(elem.column);
                }
            }
        }

        if (p_result[0].empty() || std::abs(p_result[0][0].pricing_cost) <= TOL || log.col_ratio <= config->disaggregation_threshold) {

            log.lb_incomp = log.opt_armp + data->k*log.pricing_incompatible;
            if (log.lb_incomp > result->best_lb)
                result->best_lb = log.lb_incomp;
            log.gap = std::abs(log.opt_armp - log.lb_incomp) / std::abs(log.opt_armp);
            if (log.gap < log.best_gap) {
                log.best_gap = log.gap;
            }
            log.ip_gap = std::abs(result->best_ub - result->best_lb) / std::abs(result->best_ub);
            if (log.ip_gap <= config->ip_tol || log.gap <= config->ip_tol * 0.1) {
                log.print();
                int count_lb_update = 0;
                int count_ub_update = 0;
                if (update_bounds(rms, bounds.lb_box, bounds.ub_box, count_lb_update, count_ub_update)) {
                    if (config->cg_verbose)
                        std::cout << "      Active LB: " << count_lb_update << "\tActive UB: " << count_ub_update << "\n";
                    log.best_gap = std::numeric_limits<double>::infinity();
                    result->best_lb = -std::numeric_limits<double>::infinity();
                    result->count_times_active++;
                    log.iter++;
                    continue;
                } else {
                    if (config->cg_verbose)
                        std::cout << "      Stop: problem solved!" << "\n";
                    break;
                }
            }

            if (log.iter % config->log_step == 0) {
                log.print();
            }

            if (config->best_partition_update) {

                int inc = get_number_of_incompatibilities(mlm, p_result[1][0].column);
                // use the best incompatible column (regardless of the degree) to update the aggregating partition
                sum_target_degree += inc;
                sum_best_degree += inc;
                update_aggregation_improved(mlm, p_result[1][0].column);
                // add the previously generated incompatible column that now is compatible
                or_cols.push_back(p_result[1][0].column);

            } else {

                // find the most negative reduced cost incompatible column with the smallest degree
                PricingOutput min_p_column = p_result[1][0];
                int min_p_inc = get_number_of_incompatibilities(mlm, p_result[1][0].column);
                for (auto &elem : p_result[1]) {
                    int inc = get_number_of_incompatibilities(mlm, elem.column);
                    if (inc < min_p_inc && inc != 0) {
                        min_p_column = elem;
                        min_p_inc = inc;
                    }
                }

                sum_target_degree += min_p_inc;
                sum_best_degree += get_number_of_incompatibilities(mlm, p_result[1][0].column);
                update_aggregation_improved(mlm, min_p_column.column);
                // add the previously generated incompatible column that now is compatible
                or_cols.push_back(min_p_column.column);

            }

            // display_graph(mlm);
            log.n_covering = (int) mlm.ml_comp.size();
            init_A = update_columns(mlm, or_cols);
            // compute bounds
            bounds = compute_bounds(mlm, data->P, cr_start->assignment, cr_start->centroid);
            // update ARMP
            update_stabilized_restricted_master_problem(rms, data->P, mlm, init_A, bounds.lb_box, bounds.ub_box);
            result->n_update++;
            log.iter++;
            if (config->dual_disaggregation_strategy == COMPLEMENTARITY) {
                incompatible_cols.clear();
            }
            auto stop_iter = std::chrono::high_resolution_clock::now();
            auto elapsed_iter = (std::chrono::duration_cast<std::chrono::microseconds>(stop_iter - start_iter).count());
            log.time_iter = partial_time += ((double) elapsed_iter / 1000000);
            continue;

        }

        log.lb_incomp = log.opt_armp + data->k*log.pricing_incompatible;
        if (log.lb_incomp > result->best_lb)
            result->best_lb = log.lb_incomp;
        log.gap = std::abs(log.opt_armp - log.lb_incomp) / std::abs(log.opt_armp);
        if (log.gap < log.best_gap)
            log.best_gap = log.gap;
        log.ip_gap = std::abs(result->best_ub - result->best_lb) / std::abs(result->best_ub);
        if (log.ip_gap <= config->ip_tol || log.gap <= config->ip_tol * 0.1) {
            log.print();
            int count_lb_update = 0;
            int count_ub_update = 0;
            if (update_bounds(rms, bounds.lb_box, bounds.ub_box, count_lb_update, count_ub_update)) {
                if (config->cg_verbose)
                    std::cout << "      Active LB: " << count_lb_update << "\tActive UB: " << count_ub_update << "\n";
                log.best_gap = std::numeric_limits<double>::infinity();
                result->best_lb = -std::numeric_limits<double>::infinity();
                result->count_times_active++;
                log.iter++;
                auto stop_iter = std::chrono::high_resolution_clock::now();
                auto elapsed_iter = (std::chrono::duration_cast<std::chrono::microseconds>(stop_iter - start_iter).count());
                log.time_iter = partial_time += ((double) elapsed_iter / 1000000);
                continue;
            } else {
                if (config->cg_verbose)
                    std::cout << "      Stop: problem solved!" << "\n";
                break;
            }
        }

        if (log.iter % config->log_step == 0) {
            log.print();
        }

        int n_covering = (int) rms->get(GRB_IntAttr_NumConstrs) - 1;
        int n_cols = (int) rms->get(GRB_IntAttr_NumVars) - 2 * n_covering;

        for (auto &elem : p_result[0]) {
            // keep track for the newly added column (with original indices)
            or_cols.push_back(elem.column);
            // aggregate and add compatible column
            std::vector<bool> a = aggregate_column(mlm, elem.column);
            GRBConstr *tempCons = rms->getConstrs();
            GRBColumn col;
            col.addTerm(-1.0, tempCons[0]);
            for (int j = 1; j < log.n_covering + 1; j++)
                col.addTerm(int(a[j - 1]), tempCons[j]);
            rms->addVar(0.0, GRB_INFINITY, cost_aggregation(data->P, a, mlm), GRB_CONTINUOUS, col, "C" + std::to_string(n_cols));
            n_cols++;
        }

        log.iter++;
        auto stop_iter = std::chrono::high_resolution_clock::now();
        auto elapsed_iter = (std::chrono::duration_cast<std::chrono::microseconds>(stop_iter - start_iter).count());
        log.time_iter = partial_time += ((double) elapsed_iter / 1000000);

        // stop CG algorithm with additional iterations in order to get a valid dual bound if the dual variables are active
        if (log.iter >= config->max_iter || log.time_iter >= config->time_limit) {
            int count_lb_update = 0;
            int count_ub_update = 0;
            if (update_bounds(rms, bounds.lb_box, bounds.ub_box, count_lb_update, count_ub_update)) {
                if (config->cg_verbose)
                    std::cout << "      Active LB: " << count_lb_update << "\tActive UB: " << count_ub_update << "\n";
                log.best_gap = std::numeric_limits<double>::infinity();
                result->best_lb = -std::numeric_limits<double>::infinity();
                result->count_times_active++;
                continue;
            } else {
                if (config->cg_verbose) {
                    if (log.iter >= config->max_iter)
                        std::cout << "      Stop: maximum number of iterations!" << "\n";
                    if (log.time_iter >= config->time_limit)
                        std::cout << "      Stop: time limit reached!" << "\n";
                }
                break;
            }
        }

    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto elapsed = (std::chrono::duration_cast<std::chrono::microseconds> (stop - start).count());

    // find a feasible assignment from the optimal solution
    std::vector<int> assignment = cluster_recovery_stabilized((int) mlm.ml_index.size(), or_cols, rms, data->k);
    if (!assignment.empty()) {
        ClusteringResult *cr = evaluate_cluster_assignment_2D_ptr(data->P, data->k, assignment);
        if (cr->loss < result->start_ub) {
            result->best_ub = cr->loss;
            result->cr_end = cr;
        } else {
            result->cr_end = new ClusteringResult(*cr_start);
        }
        if (config->cg_verbose)
            std::cout << "      Optimal solution is integral!" << "\n";
    } else {
        result->cr_end = new ClusteringResult(*cr_start);
        if (config->cg_verbose)
            std::cout << "      Optimal solution is not integral" << "\n";
    }
    if (log.ip_gap > config->ip_tol) {
        if (config->cg_mip_heuristic) {
            int bin_vars = 0;
            GRBModel *rms_discrete = mip_post_processing_stabilized(rms, bin_vars);
            if (config->cg_verbose)
                std::cout << "      Binary variables in MIP: " << bin_vars << "\n";
            if (rms_discrete != nullptr) {
                assignment = cluster_recovery_stabilized((int) mlm.ml_index.size(), or_cols, rms_discrete, data->k);
                if (!assignment.empty()) {
                    ClusteringResult *cr = evaluate_cluster_assignment_2D_ptr(data->P, data->k, assignment);
                    std::cout << "      Feasible solution by the MIP heuristic: " << cr->loss << "\n";
                    if (cr->loss < result->start_ub) {
                        if (config->cg_verbose)
                            std::cout << "      Better solution found by the MIP heuristic!" << "\n";
                        result->best_ub = cr->loss;
                        delete (result->cr_end);
                        result->cr_end = cr;
                    }
                }
            }
        }

    }

    // save CG statistics
    result->n_covering_end = log.n_covering;
    result->n_iter = log.iter + 1;
    result->time = (double) elapsed / 1000000;

    result->avg_best_degree = (double) sum_best_degree / result->n_update;
    result->avg_target_degree = (double) sum_target_degree / result->n_update;
    result->avg_covering = (double) sum_covering / result->n_iter;
    result->avg_n_compatible = (double) count_n_compatible / result->n_iter;
    result->avg_n_incompatible = (double) count_n_incompatible / result->n_iter;
    result->avg_time_master = (double) result->time_master / result->n_iter;
    result->avg_time_pricing = (double) result->time_pricing / result->n_iter;

    // save output data
    auto *model = new CGModel();
    model->mlm = mlm;
    model->or_cols = or_cols;
    model->ml_pairs = {};
    model->cl_pairs = {};
    result->model_end = model;
    result->grb_end = rms;

    if (config->cg_verbose) CGLog::footer();

    return result;

}