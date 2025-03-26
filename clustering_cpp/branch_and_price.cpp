#include <fstream>
#include "branch_and_price.h"
#include "cg_stabilization.h"
#include "cg_pricing.h"
#include "CGStruct.h"
#include "ThreadPool.h"
#include "cg_branching.h"
#include "cg_reduction.h"

std::pair<int, int> find_branching_pair(int n, std::vector<int> &t1, std::vector<int> &t2) {

    std::vector<bool> a1(n, false);
    std::vector<bool> a2(n, false);

    for (auto &elem : t1)
        a1[elem] = true;
    for (auto &elem : t2)
        a2[elem] = true;

    int i_index = -1;
    int j_index = -1;
    bool found = false;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (a1[i] && a1[j]) {
                if ((a2[i] && !a2[j]) || (!a2[i] && a2[j])) {
                    i_index = i;
                    j_index = j;
                    found = true;
                }
            }
            if (found) break;
        }
        if (found) break;
    }


    return {i_index, j_index};
}

std::vector<std::vector<int>> get_basic_variables(int node_id, GRBModel *grbModel, std::vector<std::vector<int>> &or_cols) {

    std::vector<std::vector<int>> basic_col;

    grbModel->optimize(); // optimize again to get an optimal basis

    if (node_id == 0) {

        int n_covering = (int) grbModel->get(GRB_IntAttr_NumConstrs) - 1;
        int n_cols = (int) grbModel->get(GRB_IntAttr_NumVars) - 2 * n_covering;

        basic_col.reserve(n_cols);
        for (int i = 0; i < n_cols; i++) {
            if (grbModel->getVarByName("C" + std::to_string(i)).get(GRB_IntAttr_VBasis) == GRB_BASIC) {
                basic_col.push_back(or_cols[i]);
            }
        }

    } else {

        GRBVar *vars = grbModel->getVars();
        int n_cols = (int) grbModel->get(GRB_IntAttr_NumVars);
        basic_col.reserve(n_cols);
        for (int i = 0; i < n_cols; i++) {
            if (vars[i].get(GRB_IntAttr_VBasis) == GRB_BASIC) {
                basic_col.push_back(or_cols[i]);
            }
        }

    }

    return basic_col;
}

std::pair<int, int> get_branching_decision(int node_id, int n, GRBModel *model, std::vector<std::vector<int>> &or_cols) {

    std::vector<std::pair<int, double>> pairs;

    if (node_id == 0) { // root node
        int n_covering = (int) model->get(GRB_IntAttr_NumConstrs) - 1;
        int n_cols = (int) model->get(GRB_IntAttr_NumVars) - 2 * n_covering;
        pairs.reserve(n_cols);
        for (int i = 0; i < n_cols; i++) {
            double v = model->getVarByName("C" + std::to_string(i)).get(GRB_DoubleAttr_X);
            if (v > 0.01 && v < 1 - 0.01) {
                pairs.emplace_back(i, std::min(v, 1 - v));
            }
        }
    } else { // ML or CL child
        GRBVar *vars = model->getVars();
        int n_cols = (int) model->get(GRB_IntAttr_NumVars);
        pairs.reserve(n_cols);
        for (int i = 0; i < n_cols; i++) {
            double v = vars[i].get(GRB_DoubleAttr_X);
            if (v > 0.01 && v < 1 - 0.01) {
                pairs.emplace_back(i, std::min(v, 1 - v));
            }
        }
    }

    if (pairs.size() < 2)
        return std::make_pair(-1, -1);

    // sort in descend order
    std::sort(pairs.begin(), pairs.end(), [](auto &left, auto &right) {
        return left.second > right.second;
    });


    int n_pairs = (int) pairs.size();
    int i_index = -1;
    int j_index = -1;
    bool found = false;
    for (int i = 0; i < n_pairs; i++) {
        for (int j = i + 1; j < n_pairs; j++) {
            std::pair<int, int> pair = find_branching_pair(n, or_cols[pairs[i].first], or_cols[pairs[j].first]);
            if (pair.first != -1 && pair.second != -1) {
                i_index = pair.first;
                j_index = pair.second;
                found = true;
            }
            if (found) break;
        }
        if (found) break;
    }

    return std::make_pair(i_index, j_index);

}

std::vector<JobData *> create_children_jobs(Node *node, SharedData *shared_data, BBConfig *bb_config) {

    int n = (int) node->cg_result->model_end->mlm.ml_index.size();

    std::pair<int, int> branching_pair = get_branching_decision(node->id, n, node->cg_result->grb_end, node->cg_result->model_end->or_cols);

    if (branching_pair.first == -1 && branching_pair.second == -1) {

        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        // std::cout << "      PRUNING BY OPTIMALITY " << node->id << "\n";
        if (node->cg_result->cr_end->loss - shared_data->global_ub <= -bb_config->bb_tol) {
            // update global upper bound
            shared_data->global_ub = node->cg_result->cr_end->loss;
            shared_data->global_assignment = node->cg_result->cr_end->assignment;
        }

        delete (node->cg_result);
        delete (node);

        return {};

        // mutex is automatically released when lock goes out of scope
    }

    std::vector<JobData *> branch;

    Node *child_node_ml = new Node(*node); // copy data for ML child
    child_node_ml->cg_result = new CGResult(*(node->cg_result));
    child_node_ml->cg_result->grb_end = new GRBModel(*(node->cg_result->grb_end));
    child_node_ml->cg_result->model_end = new CGModel(*(node->cg_result->model_end));
    child_node_ml->cg_result->cr_end = new ClusteringResult(*(node->cg_result->cr_end));

    Node *child_node_cl = new Node(*node); // copy data for CL child
    child_node_cl->cg_result = new CGResult(*(node->cg_result));
    child_node_cl->cg_result->grb_end = new GRBModel(*(node->cg_result->grb_end));
    child_node_cl->cg_result->model_end = new CGModel(*(node->cg_result->model_end));
    child_node_cl->cg_result->cr_end = new ClusteringResult(*(node->cg_result->cr_end));

    auto *ml_job_data = new JobData();
    ml_job_data->type = MUST_LINK;
    ml_job_data->branching_pair = branching_pair;
    ml_job_data->node = child_node_ml;

    auto *cl_job_data = new JobData();
    cl_job_data->type = CANNOT_LINK;
    cl_job_data->branching_pair = branching_pair;
    cl_job_data->node = child_node_cl;

    branch.emplace_back(ml_job_data);
    branch.emplace_back(cl_job_data);

    delete (node->cg_result);
    delete (node);

    return branch;

}

std::vector<JobData *> solve_child_problem(JobData *job, CGData *cg_data, SharedData *shared_data, CGConfig *cg_config, BBConfig *bb_config) {

    /*
    std::cout << "Inherited columns: " << "\n";
    for (auto &col : job->node->cg_result->model_end->or_cols) {
        sort(col.begin(), col.end());
        print_vector_int(col);
    }
    std::cout << "\n";
    */

    // choose branching subproblem
    switch (job->type) {
        case MUST_LINK:
            job->node->cg_result->model_end->ml_pairs.push_back(job->branching_pair);
            // job->node->cg_result->model_end->or_cols = get_basic_variables(job->node->id, job->node->cg_result->grb_end, job->node->cg_result->model_end->or_cols);
            job->node->cg_result->model_end->or_cols = update_columns_must_link(job->node->cg_result->model_end->or_cols, job->branching_pair);
            break;
        case CANNOT_LINK:
            job->node->cg_result->model_end->cl_pairs.push_back(job->branching_pair);
            // job->node->cg_result->model_end->or_cols = get_basic_variables(job->node->id, job->node->cg_result->grb_end, job->node->cg_result->model_end->or_cols);
            job->node->cg_result->model_end->or_cols = update_columns_cannot_link(job->node->cg_result->model_end->or_cols, job->branching_pair);
            break;
        default:
            std::cerr << "Invalid branching pair" << "\n";
            exit(EXIT_FAILURE);
    }

    CGResult *result = solve_master_problem_dca_child(cg_data, job->node->cg_result->model_end, cg_config, bb_config, job->node->cg_result->cr_end);

    // init child data
    Node *child = new Node();
    child->lb = std::max(result->best_lb, job->node->lb);
    child->cg_result = result;

    double time = result->time;
    int n_cols = (int) result->model_end->or_cols.size();
    int n_iter = result->n_iter;

    double node_gap;

    {
        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        if (child->cg_result->cr_end->loss - shared_data->global_ub <= -bb_config->bb_tol) {
            // update global upper bound
            shared_data->global_ub = child->cg_result->cr_end->loss;
            shared_data->global_assignment = child->cg_result->cr_end->assignment;
        }

        child->id = shared_data->n_nodes;
        shared_data->n_nodes++;
        int open = shared_data->queue->getSize();

        node_gap = (shared_data->global_ub - child->lb) / shared_data->global_ub;

        Node *min_lb_node = shared_data->queue->getMinLB();
        if (min_lb_node != nullptr)
            shared_data->gap = (shared_data->global_ub - min_lb_node->lb) / std::abs(shared_data->global_ub);

        BBLog log{child->id, job->type, job->branching_pair.first, job->branching_pair.second, child->lb, n_iter, n_cols,
                  result->n_covering_start, result->n_covering_end,time, result->best_ub, shared_data->global_ub, node_gap, shared_data->gap, open};
        log.print();

    }

    if (node_gap < bb_config->bb_tol) {
        // std::cout << "      PRUNING " << child->id << "\n";
        delete (child->cg_result);
        delete (child);
        return {};
    }

    return create_children_jobs(child, shared_data, bb_config);
}

std::vector<JobData *> solve_root_problem(CGData *cg_data, SharedData *shared_data, CGConfig *cg_config, BBConfig *bb_config, ClusteringResult *init_cr) {

    CGResult *result = solve_stabilized_master_problem_dca(cg_data, cg_config, init_cr);

    // init root data
    Node *node = new Node();
    node->id = 0;
    node->lb = result->best_lb;
    node->cg_result = result;

    double time = result->time;
    int n_cols = (int) result->model_end->or_cols.size();
    int n_iter = (int) result->n_iter;

    // update shared data
    shared_data->global_ub = result->best_ub;
    shared_data->global_assignment = result->cr_end->assignment;
    shared_data->n_nodes++;

    int open = shared_data->queue->getSize();

    double node_gap = (shared_data->global_ub - result->best_lb) / shared_data->global_ub;
    shared_data->gap = node_gap;

    BBLog log{node->id, ROOT, -1, -1, result->best_lb, n_iter, n_cols, result->n_covering_start, result->n_covering_end,
              time, result->best_ub, shared_data->global_ub, node_gap, shared_data->gap, open};
    log.print();

    if (node_gap <= bb_config->bb_tol) {
        delete (node->cg_result);
        delete (node);
        return {};
    }

    return create_children_jobs(node, shared_data, bb_config);
}

bool is_thread_pool_working(std::vector<bool> &thread_state) {
    int count = 0;
    for (auto && i : thread_state) {
        if (i) count++;
    }
    if (count == 0)
        return false;
    return true;
}

BBResult *branch_and_bound(CGData *cg_data, ClusteringResult *init_cr, CGConfig *cg_config, BBConfig *bb_config) {

    int n_thread = bb_config->bb_parallel;

    JobAbstractQueue *queue;
    switch (bb_config->bb_visiting_strategy) {
        case DEPTH_FIRST:
            queue = new JobStack();
            break;
        case BEST_FIRST:
            queue = new JobPriorityQueue();
            break;
        case BREADTH_FIRST:
            queue = new JobQueue();
            break;
        default:
            queue = new JobPriorityQueue();
    }

    auto *shared_data = new SharedData();
    shared_data->global_ub = init_cr->loss;
    shared_data->global_assignment = init_cr->assignment;
    shared_data->n_nodes = 0;
    shared_data->queue = queue;

    shared_data->threadStates.reserve(n_thread);
    for (int i = 0; i < n_thread; i++) {
        shared_data->threadStates.push_back(false);
    }

    ThreadPool pool(cg_config, bb_config, shared_data, cg_data, n_thread);

    BBLog::header();

    auto start_all = std::chrono::high_resolution_clock::now();

    std::vector<JobData *> jobs = solve_root_problem(cg_data, shared_data, cg_config, bb_config, init_cr);
    delete (init_cr);

    // submit jobs to the thread pool
    for (auto & job : jobs) {
        pool.addJob(job);
    }

    while (true) {

        {
            std::unique_lock<std::mutex> l(shared_data->queueMutex);
            while (is_thread_pool_working(shared_data->threadStates) && shared_data->n_nodes < bb_config->bb_max_nodes) {
                shared_data->mainConditionVariable.wait(l);
            }

            if (shared_data->queue->empty() || shared_data->n_nodes >= bb_config->bb_max_nodes)
                break;
        }

    }

    auto end_all = std::chrono::high_resolution_clock::now();
    auto elapsed = (std::chrono::duration_cast<std::chrono::microseconds> (end_all - start_all).count());

    pool.quitPool();

    BBLog::footer();

    if (queue->empty()) {
        shared_data->gap = 0.0;
    }

    std::cout << "Final MSSC solution: " << shared_data->global_ub << "\n";
    std::cout << "Gap [%]: " << shared_data->gap * 100 << "\n";
    std::cout << "Nodes: " << shared_data->n_nodes << "\n";
    std::cout << "Time [s]: " << (double) elapsed / 1000000 << "\n\n";

    auto *bb_result = new BBResult();
    bb_result->opt_assignment = shared_data->global_assignment;
    bb_result->opt_value = shared_data->global_ub;
    bb_result->opt_gap = shared_data->gap;

    // free memory
    delete (queue);
    delete (shared_data);

    return bb_result;

}

void BBLog::header() {

    std::cout << "\n" << "      ------------------------------------------------------"
                         "  BRANCH AND BOUND WITH COLUMN GENERATION AND DCA FOR THE MSSC PROBLEM  "
                         "--------------------------------------------------\n" <<
              std::setw(8) << "ID" << " " <<
              std::setw(8) << "TYPE" << " " <<
              std::setw(8) << "I" << " " <<
              std::setw(8) << "J" << " " <<
              std::setw(15) << "LB" << " " <<
              std::setw(10) << "ITER" << " " <<
              std::setw(10) << "COLS" << " " <<
              std::setw(10) << "N START" << " " <<
              std::setw(10) << "N END" << " " <<
              std::setw(10) << "TIME" << " " <<
              std::setw(15) << "UB" << " " <<
              std::setw(15) << "GUB" << " " <<
              std::setw(15) << "GAP [%]" << " " <<
              std::setw(15) << "GGAP [%]" <<
              std::setw(10) << "QUEUE" << std::endl;

}

void BBLog::footer() {

    std::cout << "      ---------------------------------------------------------------------"
                 "------------------------------------------------------"
                 "-----------------------------------------------------\n\n";

}

void BBLog::print() const {

    std::string s_type;
    switch (this->type) {
        case MUST_LINK:
            s_type = "ML";
            break;
        case CANNOT_LINK:
            s_type = "CL";
            break;
        default:
            s_type = "RT";
            break;
    }

    std::cout <<
              std::setw(8) << this->id << " " <<
              std::setw(8) << s_type << " " <<
              std::setw(8) << this->i << " " <<
              std::setw(8) << this->j << " " <<
              std::setw(15) << this->lb << " " <<
              std::setw(10) << this->iter << " " <<
              std::setw(10) << this->columns << " " <<
              std::setw(10) << this->n_covering_start << " " <<
              std::setw(10) << this->n_covering_end << " " <<
              std::setw(10) << std::fixed << std::setprecision(2) << this->time << " " <<
              std::setw(15) << this->ub << " " <<
              std::setw(15) << this->gub << " " <<
              std::setw(15) << std::fixed << std::setprecision(3) << this->gap * 100 << " " <<
              std::setw(15) << std::fixed << std::setprecision(3) << this->ggap * 100 << " " <<
              std::setw(10) << this->open;
    std::cout << "\n";

}
