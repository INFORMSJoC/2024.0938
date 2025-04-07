#include <unordered_set>
#include <chrono>
#include <fstream>
#include "cg_reduction.h"
#include "cg_util.h"
#include "cg_pricing.h"
#include "cg_pricing_multi.h"
#include "CGStruct.h"
#include "Kmeans2D.h"

// compute clustering cost (master problem)
double cost_aggregation(std::vector<Point> &P, std::vector<bool> &a, MustLinkMapping &mlm) {

    int new_n = (int) mlm.ml_comp.size();
    int card = 0;
    Point centroid = {0.0, 0.0};
    for (int i = 0; i < new_n; i++) {
        if (a[i]) {
            for (auto &idx : mlm.ml_comp[i]) {
                centroid.x += P[idx].x;
                centroid.y += P[idx].y;
                card++;
            }
        }
    }
    centroid.x = centroid.x / card;
    centroid.y = centroid.y / card;

    double cost = 0;
    for (int i = 0; i < new_n; i++)
        if (a[i]) {
            for (auto &idx : mlm.ml_comp[i]) {
                cost += std::pow(point_distance(P[idx], centroid), 2);
            }
        }

    return cost;
}

std::vector<int> cluster_recovery(int n, std::vector<std::vector<int>> &or_cols, GRBModel *rms, int k) {

    GRBVar *vars = rms->getVars();
    int n_vars = rms->get(GRB_IntAttr_NumVars);
    std::vector<int> assignment(n);
    int cluster_id = 0;
    for (int i = 0; i < n_vars; i++) {
        double v = vars[i].get(GRB_DoubleAttr_X);
        if (v <= 1.0 + 0.01 && v >= 1.0 - 0.01) {
            for (auto &elem : or_cols[i])
                assignment[elem] = cluster_id;
            cluster_id++;
        }
    }

    if (cluster_id == k) {
        return assignment;
    }
    else {
        return {};
    }

}


GRBModel *build_restricted_master_problem_aggregation_ptr(std::vector<Point> &P, std::vector<std::vector<bool>> &A, int k, MustLinkMapping &mlm) {

    // P: n x d matrix of data points
    // A: n_aggregated x n_col matrix of initial columns
    // k: number of clusters
    // ml_map: mapping between points ad aggregated points

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
        master->addConstr(lhs_sum >= 1); // each data point must be assigned to a single cluster
    }

    master->update();

    master->set("OutputFlag", "0");
    master->set("Method", "1"); // -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier

    return master;

}

// new warm start (do not build a new model)
void update_restricted_master_problem(GRBModel *master, std::vector<Point> &P, MustLinkMapping &mlm, std::vector<std::vector<bool>> A) {

    int old_n_covering = (int) master->get(GRB_IntAttr_NumConstrs);
    int new_n_covering = (int) mlm.ml_comp.size();
    int old_n_cols = (int) master->get(GRB_IntAttr_NumVars);
    int new_n_cols = (int) A.size();

    GRBLinExpr lhs_sum = 0;
    for (int i = old_n_covering - 1; i < new_n_covering; i++) {
        lhs_sum = 0;
        for (int j = 0; j < old_n_cols; j++) {
            if (A[j][i])
                lhs_sum += master->getVar(j);
        }
        master->addConstr(lhs_sum >= 1);
    }

    master->update();

    for (int i = old_n_cols; i < new_n_cols; i++) {
        GRBConstr *tempCons = master->getConstrs();
        GRBColumn col;
        col.addTerm(-1.0, tempCons[0]);
        for (int j = 1; j < new_n_covering + 1; j++)
            col.addTerm(int(A[i][j - 1]), tempCons[j]);
        master->addVar(0.0, GRB_INFINITY, cost_aggregation(P, A[i], mlm), GRB_CONTINUOUS, col);
    }

}

GRBModel *mip_post_processing(GRBModel *rms, int &bin_vars) {

    auto *rms_discrete = new GRBModel(*rms);
    rms_discrete->optimize();

    int n_vars = (int) rms->get(GRB_IntAttr_NumVars);
    bool is_integer = true;
    for (int i = 0; i < n_vars; i++) {
        double v = rms_discrete->getVar(i).get(GRB_DoubleAttr_X);
        if (v > 0.01 && v < 1 - 0.01) {
            bin_vars++;
            rms_discrete->getVar(i).set(GRB_CharAttr_VType, GRB_BINARY);
            is_integer = false;
        }
    }
    if (!is_integer) {
        rms_discrete->set("OutputFlag", "0");
        rms_discrete->set("Threads", "4");
        rms_discrete->optimize();
        // std::cout << "MIP Objective: " << rms_discrete.get(GRB_DoubleAttr_ObjVal) << "\n";
        return rms_discrete;
    }

    delete (rms_discrete);
    return nullptr;

}

CGResult *solve_master_problem_dca_child(CGData *data, CGModel *model, CGConfig *config, BBConfig *bb_config, ClusteringResult *cr_start) {

    // CG result
    auto *result = new CGResult();
    result->n = (int) data->P.size();
    result->k = data->k;
    result->n_components = config->n_components;

    result->start_ub = cr_start->loss;
    result->best_ub = cr_start->loss;

    if (!config->cg_inherit_partition_child) {
        // do not inherit the aggregating partition: aggregate from scratch
        model->mlm = get_ml_map(data->partition_assignment, cr_start->assignment);
        std::vector<std::vector<int>> new_or_cols;
        new_or_cols.reserve(model->or_cols.size());
        for (auto &col : model->or_cols) {
            if (isCompatible(model->mlm, col))
                new_or_cols.push_back(col);
        }
        model->or_cols = new_or_cols;
    }

    auto start = std::chrono::high_resolution_clock::now();

    // CG log
    CGLog log{};
    log.verbose = config->cg_verbose;

    // Initial number of covering constraints
    log.n_covering = (int) model->mlm.ml_comp.size();
    result->n_covering_start = log.n_covering;

    // add artificial column to prevent unbounded dual
    int or_n = (int) data->P.size();
    std::vector<int> artificial_column(or_n);
    for (int i = 0; i < or_n; i++) artificial_column.push_back(i);
    model->or_cols.push_back(artificial_column);

    std::vector<std::vector<bool>> init_A;
    model->or_cols.reserve(config->max_iter);
    // use provided columns and aggregate them
    init_A = update_columns(model->mlm, model->or_cols);
    // std::cout << "Initial columns: " << model->or_cols.size() << "\n";

    GRBModel *constr_model = build_branching_problem_decomposed_ptr(data->P, model->ml_pairs, model->cl_pairs);
    GRBModel *rms = build_restricted_master_problem_aggregation_ptr(data->P, init_A, data->k, model->mlm);

    if (config->cg_verbose) CGLog::header();

    log.n_compatible = 0;
    log.n_incompatible = 0;
    log.col_ratio = 0.0;
    log.opt_armp = 0.0;
    log.col_armp = 0;
    log.iter = 0;
    log.gap = std::numeric_limits<double>::infinity();
    log.best_gap = std::numeric_limits<double>::infinity();

    result->time_master = 0.0;
    result->time_pricing = 0.0;
    result->n_update = 0;
    result->count_times_active = 0;
    result->best_lb = -std::numeric_limits<double>::infinity();

    int sum_covering = 0;
    int count_n_incompatible = 0;
    int count_n_compatible = 0;

    double partial_time = 0.0;

    while (true) {

        auto start_iter = std::chrono::high_resolution_clock::now();

        sum_covering += log.n_covering;
        log.col_armp = (int) model->or_cols.size();
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
        std::vector<double> dual_c = disaggregate_dual_variables_uniform(model->mlm, dual_c_shr);

        // (compatible, incompatible)
        auto start_pricing = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<PricingOutput>> p_result;
        // constrained pricing
        double relaxed_rc = 0.0;
        bool exact_pricing = false;
        p_result = get_compatible_columns(relaxed_rc, exact_pricing, constr_model, model->ml_pairs, model->cl_pairs, data, model->mlm,
                                          dual_c, dual_k, config->max_comp_cols_child, config->cg_heuristic_pricing_child, bb_config->bb_verbose);
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

        if (p_result[0].empty() || std::abs(p_result[0][0].pricing_cost) <= TOL || log.col_ratio <= config->disaggregation_threshold) {

            if (exact_pricing)
                log.lb_incomp = log.opt_armp + data->k*log.pricing_incompatible;
            else
                log.lb_incomp = log.opt_armp + data->k*relaxed_rc;

            if (log.lb_incomp > result->best_lb)
                result->best_lb = log.lb_incomp;
            log.gap = (std::abs(log.opt_armp - log.lb_incomp) / std::abs(log.opt_armp));
            if (log.gap < log.best_gap) {
                log.best_gap = log.gap;
            }
            if (std::abs(log.gap) <= TOL) {
                log.print();
                if (config->cg_verbose)
                    std::cout << "      Stop: problem solved!" << "\n";
                break;
            }
            if ((result->best_ub - result->best_lb) / result->best_ub <= bb_config->bb_tol) {
                log.print();
                if (config->cg_verbose)
                    std::cout << "      Stop: problem pruned!" << "\n";
                break;
            }

            if (log.iter % config->log_step == 0) {
                log.print();
            }

            if (config->best_partition_update) {

                // use the best incompatible column (regardless of the degree) to update the aggregating partition
                update_aggregation_improved(model->mlm, p_result[1][0].column);
                // add the previously generated incompatible column that now is compatible
                model->or_cols.push_back(p_result[1][0].column);

            } else {

                // find the most negative reduced cost incompatible column with the smallest degree
                PricingOutput min_p_column = p_result[1][0];
                int min_p_inc = get_number_of_incompatibilities(model->mlm, p_result[1][0].column);
                for (auto &elem : p_result[1]) {
                    int inc = get_number_of_incompatibilities(model->mlm, elem.column);
                    if (inc < min_p_inc && inc != 0) {
                        min_p_column = elem;
                        min_p_inc = inc;
                    }
                }

                update_aggregation_improved(model->mlm, min_p_column.column);
                // add the previously generated incompatible column that now is compatible
                model->or_cols.push_back(min_p_column.column);

            }

            // display_graph(mlm);
            log.n_covering = (int) model->mlm.ml_comp.size();
            init_A = update_columns(model->mlm, model->or_cols);
            // update ARMP
            update_restricted_master_problem(rms, data->P, model->mlm, init_A);
            result->n_update++;
            log.iter++;
            auto stop_iter = std::chrono::high_resolution_clock::now();
            auto elapsed_iter = (std::chrono::duration_cast<std::chrono::microseconds>(stop_iter - start_iter).count());
            log.time_iter = partial_time += ((double) elapsed_iter / 1000000);
            continue;

        }

        if (exact_pricing)
            log.lb_incomp = log.opt_armp + data->k*log.pricing_incompatible;
        else
            log.lb_incomp = log.opt_armp + data->k*relaxed_rc;

        if (log.lb_incomp > result->best_lb)
            result->best_lb = log.lb_incomp;
        log.gap = (std::abs(log.opt_armp - log.lb_incomp) / std::abs(log.opt_armp));
        if (log.gap < log.best_gap) {
            log.best_gap = log.gap;
        }
        if ((result->best_ub - result->best_lb) / result->best_ub <= bb_config->bb_tol) {
            log.print();
            if (config->cg_verbose)
                std::cout << "      Stop: problem pruned!" << "\n";
            break;
        }


        if (log.iter % config->log_step == 0) {
            log.print();
        }

        for (auto &elem : p_result[0]) {

            // keep track for the newly added column (with original indices)
            model->or_cols.push_back(elem.column);
            // aggregate and add compatible column
            std::vector<bool> a = aggregate_column(model->mlm, elem.column);
            GRBConstr *tempCons = rms->getConstrs();
            GRBColumn col;
            col.addTerm(-1.0, tempCons[0]);
            for (int j = 1; j < log.n_covering + 1; j++)
                col.addTerm(int(a[j - 1]), tempCons[j]);
            rms->addVar(0.0, GRB_INFINITY, cost_aggregation(data->P, a, model->mlm), GRB_CONTINUOUS, col);

        }

        log.iter++;
        auto stop_iter = std::chrono::high_resolution_clock::now();
        auto elapsed_iter = (std::chrono::duration_cast<std::chrono::microseconds>(stop_iter - start_iter).count());
        log.time_iter = partial_time += ((double) elapsed_iter / 1000000);

        if (log.iter >= config->max_iter) {
            log.print();
            rms->optimize();
            if (config->cg_verbose)
                std::cout << "      Stop: maximum number of iterations!" << "\n";
            break;
        }

    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto elapsed = (std::chrono::duration_cast<std::chrono::microseconds> (stop - start).count());

    // find a feasible assignment
    std::vector<int> assignment = cluster_recovery(result->n, model->or_cols, rms, data->k);
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
    }
    else {
        result->cr_end = new ClusteringResult(*cr_start);
        if (config->cg_verbose)
            std::cout << "      Optimal solution is not integral" << "\n";
        if (config->cg_mip_heuristic) {
            int bin_vars = 0;
            GRBModel *rms_discrete = mip_post_processing(rms, bin_vars);
            if (config->cg_verbose)
                std::cout << "      Binary variables in MIP: " << bin_vars << "\n";
            if (rms_discrete != nullptr) {
                assignment = cluster_recovery(result->n, model->or_cols, rms_discrete, data->k);
                if (!assignment.empty()) {
                    ClusteringResult *cr = evaluate_cluster_assignment_2D_ptr(data->P, data->k, assignment);
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

    result->avg_best_degree = NAN;
    result->avg_target_degree = NAN;
    result->avg_covering = (double) sum_covering / result->n_iter;
    result->avg_n_compatible = (double) count_n_compatible / result->n_iter;
    result->avg_n_incompatible = (double) count_n_incompatible / result->n_iter;
    result->avg_time_master = (double) result->time_master / result->n_iter;
    result->avg_time_pricing = (double) result->time_pricing / result->n_iter;

    // save output data
    auto *model_end = new CGModel();
    model_end->mlm = model->mlm;
    model_end->or_cols = model->or_cols;
    model_end->cl_pairs = model->cl_pairs;
    model_end->ml_pairs = model->ml_pairs;
    result->model_end = model_end;
    result->grb_end = rms;

    if (config->cg_verbose) CGLog::footer();

    return result;

}


/*
GRBModel *build_restricted_master_problem_aggregation_ptr(std::vector<Point> &P, std::vector<std::vector<bool>> &A, int k, MustLinkMapping &mlm, double bigM) {

    // P: n x d matrix of data points
    // A: n_aggregated x n_col matrix of initial columns
    // k: number of clusters
    // ml_map: mapping between points ad aggregated points

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
    GRBVar s = master->addVar(0.0, GRB_INFINITY, bigM, GRB_CONTINUOUS, "S");

    // Constraints of the master problem
    GRBLinExpr sum_k = 0;
    for (int i = 0; i < n_cols; i++)
        sum_k += VarMaster[i];
    master->addConstr(-sum_k - s >= -k); // we want k clusters

    GRBLinExpr lhs_sum = 0;
    for (int i = 0; i < n_aggregated; i++) {
        lhs_sum = 0;
        for (int j = 0; j < n_cols; j++) {
            if (A[j][i])
                lhs_sum += VarMaster[j];
        }
        lhs_sum += s;
        master->addConstr(lhs_sum >= 1); // each data point must be assigned to a single cluster
    }

    master->set("OutputFlag", "0");
    master->set("Method", "1"); // -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier

    return master;

}
*/


/*
CGResult *solve_master_problem_dca(CGData *data, CGConfig *config, ClusteringResult *cr_start) {

    // CG result
    auto *result = new CGResult();
    result->n = (int) data->P.size();
    result->k = data->k;
    result->r_perc = config->perc;

    result->start_ub = cr_start->loss;
    result->best_ub = cr_start->loss;

    std::vector<std::pair<int, int>> candidate_ml_pairs = candidate_must_link_radius(data->k, config->perc, cr_start->centroid, cr_start->points_cluster);
    MustLinkMapping mlm = get_ml_map(result->n, candidate_ml_pairs);

    auto start = std::chrono::high_resolution_clock::now();

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
    // get_initial_columns_2_best(cr_start->assignment, cr_start->centroid, data->P, or_cols);
    init_A = get_initial_columns_aggregation(cr_start->points_cluster, mlm, or_cols);
    // std::cout << "Initial columns: " << or_cols.size() << "\n";

    GRBModel *rms = build_restricted_master_problem_aggregation_ptr(data->P, init_A, data->k, mlm);

    if (config->cg_verbose) CGLog::header();

    log.col_ratio = 0.0;
    log.opt_armp = 0.0;
    log.col_armp = 0;
    log.iter = 0;
    log.gap = std::numeric_limits<double>::infinity();
    log.best_gap = std::numeric_limits<double>::infinity();

    result->time_master = 0.0;
    result->time_pricing = 0.0;
    result->n_update = 0;
    result->count_times_active = 0;
    result->best_lb = -std::numeric_limits<double>::infinity();

    std::vector<int> sum_best_degree; // best incompatible column;
    sum_best_degree.reserve(config->max_iter);

    std::vector<int> sum_target_degree; // best p-incompatible column;
    sum_target_degree.reserve(config->max_iter);

    std::vector<int> sum_covering;
    sum_covering.reserve(config->max_iter);

    int count_n_incompatible = 0;
    int count_n_compatible = 0;

    while (true) {

        sum_covering.push_back(log.n_covering);
        log.col_armp = (int) or_cols.size();
        auto start_master = std::chrono::high_resolution_clock::now();
        rms->optimize();

        auto stop_master = std::chrono::high_resolution_clock::now();
        auto elapsed_master = (std::chrono::duration_cast<std::chrono::microseconds>(
                stop_master - start_master).count());
        log.time_armp = (double) elapsed_master / 1000000;
        log.opt_armp = rms->get(GRB_DoubleAttr_ObjVal);
        result->time_master += log.time_armp;

        //Dual variables associated to point constraints
        double dual_k = rms->getConstr(0).get(GRB_DoubleAttr_Pi);
        std::vector<double> dual_c_shr(log.n_covering);
        for (int j = 1; j < log.n_covering + 1; j++)
            dual_c_shr[j - 1] = rms->getConstr(j).get(GRB_DoubleAttr_Pi);

        // disaggregate dual variables before calling the pricing problem
        std::vector<double> dual_c = disaggregate_dual_variables_uniform(mlm, dual_c_shr);

        // (compatible, incompatible)
        auto start_pricing = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<PricingOutput>> p_result = get_compatible_columns(mlm, data->od, dual_c, dual_k, data->P, config->max_comp_cols, config->max_cols);
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

        if (p_result[0].empty() || std::abs(p_result[0][0].pricing_cost) <= TOL || log.col_ratio <= config->disaggregation_threshold) {

            log.lb_incomp = log.opt_armp + data->k*log.pricing_incompatible;
            if (log.lb_incomp > result->best_lb)
                result->best_lb = log.lb_incomp;
            log.gap = (std::abs(log.opt_armp - log.lb_incomp) / std::abs(log.opt_armp));
            if (log.gap < log.best_gap) {
                log.best_gap = log.gap;
            }
            if (std::abs(log.gap) <= TOL) {
                log.print();
                if (config->cg_verbose)
                    std::cout << "      Stop: problem solved!" << "\n";
                break;
            }

            if (log.iter % config->log_step == 0) {
                log.print();
            }

            if (config->best_partition_update) {

                int inc = get_number_of_incompatibilities(mlm, p_result[1][0].column);
                // use the best incompatible column (regardless of the degree) to update the aggregating partition
                sum_target_degree.push_back(inc);
                sum_best_degree.push_back(inc);
                update_aggregation_improved(mlm, p_result[1][0].column);
                or_cols.push_back(p_result[1][0].column);

            } else {

                // find the minimum p-incompatible column
                PricingOutput min_p_column = p_result[1][0];
                int min_p_inc = get_number_of_incompatibilities(mlm, p_result[1][0].column);
                for (auto &elem : p_result[1]) {
                    int inc = get_number_of_incompatibilities(mlm, elem.column);
                    if (inc < min_p_inc) {
                        min_p_column = elem;
                        min_p_inc = inc;
                    }
                }

                sum_target_degree.push_back(min_p_inc);
                sum_best_degree.push_back(get_number_of_incompatibilities(mlm, p_result[1][0].column));
                update_aggregation_improved(mlm, min_p_column.column);
                // add the previously generated incompatible column that now is compatible
                or_cols.push_back(min_p_column.column);

            }

            // display_graph(mlm);
            log.n_covering = (int) mlm.ml_comp.size();
            init_A = update_columns(mlm, or_cols);
            // update ARMP
            update_restricted_master_problem(rms, data->P, mlm, init_A);
            result->n_update++;
            log.iter++;
            continue;

        }

        log.lb_incomp = log.opt_armp + data->k*log.pricing_incompatible;
        if (log.lb_incomp > result->best_lb)
            result->best_lb = log.lb_incomp;
        log.gap = (std::abs(log.opt_armp - log.lb_incomp) / std::abs(log.opt_armp));
        if (log.gap < log.best_gap) {
            log.best_gap = log.gap;
        }
        if (std::abs(log.gap) <= TOL) {
            log.print();
            if (config->cg_verbose)
                std::cout << "      Stop: problem solved!" << "\n\n";
            break;
        }

        if (log.iter % config->log_step == 0) {
            log.print();
        }

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
            rms->addVar(0.0, GRB_INFINITY, cost_aggregation(data->P, a, mlm), GRB_CONTINUOUS, col);

        }

        log.iter++;

        if (log.iter >= config->max_iter) {
            log.print();
            rms->optimize();
            if (config->cg_verbose)
                std::cout << "      Stop: maximum number of iterations!" << "\n";
            break;
        }

    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto elapsed = (std::chrono::duration_cast<std::chrono::microseconds> (stop - start).count());

    // find a feasible assignment
    std::vector<int> assignment = cluster_recovery(result->n, or_cols, rms, data->k);
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
    }
    else {
        result->cr_end = new ClusteringResult(*cr_start);
        if (config->cg_verbose)
            std::cout << "      Optimal solution is not integral" << "\n";
    }

    // save CG statistics
    result->n_covering_end = log.n_covering;
    result->n_iter = log.iter;
    result->time = (double) elapsed / 1000000;

    result->avg_best_degree = NAN;
    if (!sum_best_degree.empty()) {
        double tmp_sum = 0.0;
        for (auto &elem: sum_best_degree) tmp_sum += elem;
        result->avg_best_degree = tmp_sum / (double) sum_best_degree.size();
    }

    result->avg_target_degree = NAN;
    if (!sum_target_degree.empty()) {
        double tmp_sum = 0.0;
        for (auto &elem: sum_target_degree) tmp_sum += elem;
        result->avg_target_degree = tmp_sum / (double) sum_target_degree.size();
    }

    double tmp_sum = 0.0;
    for (auto &elem : sum_covering) tmp_sum += elem;
    result->avg_covering = tmp_sum / (double) sum_covering.size();

    result->avg_n_compatible = (double) count_n_compatible / result->n_iter;
    result->avg_n_incompatible = (double) count_n_incompatible / result->n_iter;

    result->avg_time_master = (double) result->time_master / result->n_iter;
    result->avg_time_pricing = (double) result->time_pricing / result->n_iter;

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
*/
