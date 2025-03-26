#include "cg_branching.h"
#include <unordered_set>
#include <chrono>
#include "cg_util.h"

GRBModel *build_bqp_problem_ptr(std::vector<bool> &isVar,
                                double f_old, std::vector<std::vector<double>> &D, std::vector<double> &dual_c,
                                std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs) {

    // isVar: bool vector indicating the variables to consider in the problem
    // D: n x n distance matrix
    // dual_c: vector of dual variables (covering constraints)
    // ml_pairs, cl_pairs: must-link and cannot-link constraints

    int n = (int) D.size();

    auto env = GRBEnv(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    auto model = new GRBModel(env);

    std::vector<GRBVar> x(n);
    for (int i = 0; i < n; i++) {
        if (isVar[i])
            x[i] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, "x" + std::to_string(i));
    }

    GRBQuadExpr q_temp = 0.0;
    GRBLinExpr l_temp = 0.0;
    for (int i = 0; i < n ; i++) {
        for (int j = i + 1; j < n; j++) {
            if (isVar[i] && isVar[j])
                q_temp += (D[i][j] - dual_c[i] - dual_c[j]) * x[i] * x[j];
        }
        if (isVar[i])
            l_temp += (dual_c[i] + f_old) * x[i];
    }

    if (!cl_pairs.empty()) {
        // cannot-link constraints
        for (auto &pair: cl_pairs) {
            if (isVar[pair.first] && isVar[pair.second])
                model->addConstr(x[pair.first] + x[pair.second] <= 1);
        }
    }

    if (!ml_pairs.empty()) {
        // must-link constraints
        for (auto &pair: ml_pairs) {
            if (isVar[pair.first] && isVar[pair.second])
                model->addConstr(x[pair.first] - x[pair.second] == 0);
        }
    }

    model->setObjective(q_temp - l_temp);

    model->set("OutputFlag", "0");
    model->set("Threads", "4");

    return model;

}

GRBModel *build_bqp_problem_ptr(double f_old, std::vector<std::vector<double>> &D, std::vector<double> &dual_c,
                                std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs) {

    // P: n x d matrix of data points
    // D: n x n distance matrix
    // dual_c: vector of dual variables (covering constraints)
    // ml_pairs, cl_pairs: must-link and cannot-link constraints

    int n = (int) D.size();

    auto env = GRBEnv(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    auto model = new GRBModel(env);

    std::vector<GRBVar> x(n);
    for (int i = 0; i < n; i++) {
        x[i] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, "x" + std::to_string(i));
    }

    GRBQuadExpr q_temp = 0.0;
    GRBLinExpr l_temp = 0.0;
    for (int i = 0; i < n ; i++) {
        for (int j = i + 1; j < n; j++) {
            q_temp += (D[i][j] - dual_c[i] - dual_c[j]) * x[i] * x[j];
        }
        l_temp += (dual_c[i] + f_old) * x[i];
    }

    if (!cl_pairs.empty()) {
        // cannot-link constraints
        for (auto &pair: cl_pairs) {
            model->addConstr(x[pair.first] + x[pair.second] <= 1);
        }
    }

    if (!ml_pairs.empty()) {
        // must-link constraints
        for (auto &pair: ml_pairs) {
            model->addConstr(x[pair.first] - x[pair.second] == 0);
        }
    }

    model->setObjective(q_temp - l_temp);

    model->set("OutputFlag", "0");

    return model;

}

double compute_f(std::vector<bool> &isVar, std::vector<double> &x_opt, std::vector<std::vector<double>> &D, std::vector<double> &dual_c) {

    // x_opt: optimal solution
    // D: distance matrix
    // dual_c: vector of dual variables

    int n = (int) dual_c.size();
    double q_temp = 0.0;
    double l_temp = 0.0;
    double den = 0.0;
    for (int i = 0; i < n ; i++) {
        for (int j = i + 1; j < n; j++) {
            if (isVar[i] && isVar[j])
                q_temp += (D[i][j] - dual_c[i] - dual_c[j]) * x_opt[i] * x_opt[j];
        }
        if (isVar[i]) {
            l_temp += dual_c[i] * x_opt[i];
            den += x_opt[i];
        }
    }

    return (q_temp - l_temp) / den;
}

double compute_f(std::vector<double> &x_opt, std::vector<std::vector<double>> &D, std::vector<double> &dual_c) {

    // x_opt: optimal solution
    // D: distance matrix
    // dual_c: vector of dual variables

    int n = (int) dual_c.size();
    double q_temp = 0.0;
    double l_temp = 0.0;
    double den = 0.0;
    for (int i = 0; i < n ; i++) {
        for (int j = i + 1; j < n; j++) {
            q_temp += (D[i][j] - dual_c[i] - dual_c[j]) * x_opt[i] * x_opt[j];
        }
        l_temp += dual_c[i] * x_opt[i];
        den += x_opt[i];
    }

    return (q_temp - l_temp) / den;
}

// pricing problem considering only a subset of data points
DinkelbachResult dinkelbach_pricing(std::vector<bool> &isVar,
                                    double f_init, std::vector<Point> &P, std::vector<std::vector<double>> &D, std::vector<double> &dual_c, double dual_k,
                                    std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs) {

    int n = (int) P.size();

    std::vector<double> x_opt(n);
    std::vector<int> column;
    column.reserve(n);

    int iter = 0;
    double f = f_init;
    bool stop = false;

    auto start_pricing = std::chrono::high_resolution_clock::now();

    while (!stop) {

        GRBModel *m = build_bqp_problem_ptr(isVar, f, D, dual_c, ml_pairs, cl_pairs);
        m->optimize();
        iter++;
        double f_new = m->get(GRB_DoubleAttr_ObjVal);
        // std::cout << "f_" << iter << ": " << f_new << "\n";
        if (f_new < -TOL) {
            x_opt.clear();
            column.clear();
            for (int i = 0; i < n; i++) {
                if (isVar[i])
                    x_opt[i] = m->getVarByName("x" + std::to_string(i)).get(GRB_DoubleAttr_X);
            }
            f = compute_f(isVar, x_opt, D, dual_c);
            for (int i = 0; i < n; i++) {
                if (x_opt[i] == 1.0)
                    column.push_back(i);
            }
        } else {
            stop = true;
            //std::cout << "\n STOPPED WITH RC: " << dual_k + f << "\n";
            //std::cout << "ITER: " << iter << "\n\n";
        }

        delete (m);

    }

    auto stop_pricing = std::chrono::high_resolution_clock::now();
    auto elapsed_pricing = (std::chrono::duration_cast<std::chrono::microseconds>(stop_pricing - start_pricing).count());
    double elapsed_time = (double) elapsed_pricing / 1000000;

    return {{dual_k + f, column}, iter, elapsed_time};

}

// pricing problem considering all the data points
DinkelbachResult dinkelbach_pricing(double f_init, std::vector<Point> &P, std::vector<std::vector<double>> &D, std::vector<double> &dual_c, double dual_k,
                                    std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs) {

    int n = (int) P.size();

    std::vector<double> x_opt(n);
    std::vector<int> column;
    column.reserve(n);

    int iter = 0;
    double f = f_init;
    bool stop = false;

    auto start_pricing = std::chrono::high_resolution_clock::now();

    while (!stop) {

        GRBModel *m = build_bqp_problem_ptr(f, D, dual_c, ml_pairs, cl_pairs);
        m->optimize();
        iter++;
        double f_new = m->get(GRB_DoubleAttr_ObjVal);
        // std::cout << "f_" << iter << ": " << f_new << "\n";
        if (f_new < -TOL) {
            x_opt.clear();
            column.clear();
            for (int i = 0; i < n; i++) {
                x_opt[i] = m->getVarByName("x" + std::to_string(i)).get(GRB_DoubleAttr_X);
            }
            f = compute_f(x_opt, D, dual_c);
            for (int i = 0; i < n; i++) {
                if (x_opt[i] == 1.0)
                    column.push_back(i);
            }
        } else {
            stop = true;
            //std::cout << "\n STOPPED WITH RC: " << dual_k + f << "\n";
            //std::cout << "ITER: " << iter << "\n\n";
        }

        delete (m);

    }

    auto stop_pricing = std::chrono::high_resolution_clock::now();
    auto elapsed_pricing = (std::chrono::duration_cast<std::chrono::microseconds>(stop_pricing - start_pricing).count());
    double elapsed_time = (double) elapsed_pricing / 1000000;

    return {{dual_k + f, column}, iter, elapsed_time};

}

GRBModel *build_branching_problem_decomposed_ptr(std::vector<Point> &P, std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs) {

    // P: n x d matrix of data points
    // ml_pairs, cl_pairs: must-link and cannot-link constraints

    int n = (int) P.size();

    auto env = GRBEnv(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    auto model = new GRBModel(env);

    // add constraints of the binary problem
    std::vector<GRBVar> v(n);
    GRBLinExpr obj = 0;
    for (int i = 0; i < n; i++) {
        v[i] = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
    }

    if (!cl_pairs.empty()) {
        // cannot-link constraints
        for (auto &pair: cl_pairs) {
            model->addConstr(v[pair.first] + v[pair.second] <= 1);
        }
    }

    if (!ml_pairs.empty()) {
        // must-link constraints
        for (auto &pair: ml_pairs) {
            model->addConstr(v[pair.first] - v[pair.second] == 0);
        }
    }


    model->set("OutputFlag", "0");
    model->update();

    return model;

}


PricingOutput pricing_heuristic(GRBModel *lp_model, PricingOutput &init_po, CGData *data, std::vector<double> &price_c, double price_k) {

    int n = (int) data->P.size();
    PricingOutput po = init_po;
    po.column.reserve(n);

    // int iter = 0;
    GRBVar *vars = lp_model->getVars();
    double old_rc = INFINITY;
    while (std::abs(old_rc - po.pricing_cost) > TOL) {
        // iter++;

        Point y = {};
        int col_size = (int) po.column.size();
        for (auto &elem: po.column) {
            y.x += data->P[elem].x;
            y.y += data->P[elem].y;
        }
        y.x = y.x / col_size;
        y.y = y.y / col_size;

        double dist;
        for (int ii = 0; ii < n; ii++) {
            dist = point_distance(data->P[ii], y);
            vars[ii].set(GRB_DoubleAttr_Obj, (dist * dist) - price_c[ii]);
        }

        lp_model->optimize();
        old_rc = po.pricing_cost;
        po.pricing_cost = price_k + lp_model->get(GRB_DoubleAttr_ObjVal);
        /*
        std::cout << "Iter: " << iter << "    Pricing: " << po.pricing_cost << "\n";
        */
        po.column.clear();
        for (int j = 0; j < n; j++) {
            if (std::abs(vars[j].get(GRB_DoubleAttr_X) - 1.0) < TOL)
                po.column.push_back(j);
        }

        if (po.column.empty())
            break;
    }

    return po;
}

std::vector<std::vector<int>> update_columns_must_link(std::vector<std::vector<int>> &or_cols, std::pair<int, int> &branching_pair) {

    std::vector<std::vector<int>> upd_or_cols;
    upd_or_cols.reserve(or_cols.capacity());

    for (auto & or_col : or_cols) {
        bool first_found = false;
        bool second_found = false;
        for (auto &elem : or_col) {
            if (elem == branching_pair.first)
                first_found = true;
            if (elem == branching_pair.second)
                second_found = true;
        }
        if ((first_found && second_found) || (!first_found && !second_found))
            upd_or_cols.push_back(or_col);
    }

    return upd_or_cols;
}

std::vector<std::vector<int>> update_columns_cannot_link(std::vector<std::vector<int>> &or_cols, std::pair<int, int> &branching_pair) {

    std::vector<std::vector<int>> upd_or_cols;
    upd_or_cols.reserve(or_cols.capacity());

    for (auto & or_col : or_cols) {
        bool first_found = false;
        bool second_found = false;
        for (auto &elem : or_col) {
            if (elem == branching_pair.first)
                first_found = true;
            if (elem == branching_pair.second)
                second_found = true;
        }
        if (!(first_found && second_found))
            upd_or_cols.push_back(or_col);
    }
    return upd_or_cols;
}

/*
// PRICING PROBLEM WITH CONVEX MINLP (BIG-M FORMULATION)
GRBModel *build_branching_problem_ptr(std::vector<Point> &P, PricingMINLP &epc, double dual_k, std::vector<double> &dual_c, std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs) {

    // P: n x d matrix of data points

    int n = (int) P.size();

    auto env = GRBEnv(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    auto model = new GRBModel(env);

    // Centroid coordinates
    GRBVar y_coord_x = model->addVar(epc.y_coord_lb.first, epc.y_coord_ub.first, 0.0, GRB_CONTINUOUS, "y0");
    GRBVar y_coord_y = model->addVar(epc.y_coord_lb.second, epc.y_coord_ub.second, 0.0, GRB_CONTINUOUS, "y1");

    // add constraints of the convex MINLP reformulation
    std::vector<GRBVar> v(n);
    std::vector<GRBVar> w(n);
    GRBQuadExpr expr = 0.0;
    GRBLinExpr obj_w = 0.0;
    GRBLinExpr obj_v = 0.0;
    for (int i = 0; i < n; i++) {
        v[i] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, "v" + std::to_string(i));
        w[i] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "w" + std::to_string(i));
        expr = (P[i].x - y_coord_x) * (P[i].x - y_coord_x) + (P[i].y - y_coord_y) * (P[i].y - y_coord_y);
        model->addConstr(w[i] <= epc.bigM[i] * v[i], "bigM_le" + std::to_string(i));
        model->addQConstr(w[i] >= expr - epc.bigM[i] * (1.0 - v[i]), "bigM_ge" + std::to_string(i));
        obj_w += w[i];
        obj_v += dual_c[i] * v[i];
    }

    if (!cl_pairs.empty()) {
        // cannot-link constraints
        for (auto &pair: cl_pairs) {
            model->addConstr(v[pair.first] + v[pair.second] <= 1);
        }
    }

    if (!ml_pairs.empty()) {
        // must-link constraints
        for (auto &pair: ml_pairs) {
            model->addConstr(v[pair.first] - v[pair.second] == 0);
        }
    }

    model->setObjective(dual_k + obj_w - obj_v);

    model->set("OutputFlag", "1");

    return model;

}

PricingMINLP compute_pricing_constants(std::vector<Point> &P, std::vector<OrderedDistances> &od) {

    int n = (int) P.size();
    PricingMINLP epc;
    epc.bigM.reserve(n);
    for (int i = 0; i < n; i++) {
        epc.bigM[i] = od[i].dist[n - 2] * od[i].dist[n - 2];
    }

    epc.y_coord_lb.first = std::numeric_limits<double>::infinity();
    epc.y_coord_lb.second = std::numeric_limits<double>::infinity();
    epc.y_coord_ub.first = -std::numeric_limits<double>::infinity();
    epc.y_coord_ub.second = -std::numeric_limits<double>::infinity();

    for (int i = 0; i < n; i++) {
        if (P[i].x < epc.y_coord_lb.first)
            epc.y_coord_lb.first = P[i].x;
        if (P[i].y < epc.y_coord_lb.second)
            epc.y_coord_lb.second = P[i].y;
        if (P[i].x > epc.y_coord_ub.first)
            epc.y_coord_ub.first = P[i].x;
        if (P[i].y > epc.y_coord_ub.second)
            epc.y_coord_ub.second = P[i].y;
    }

    return epc;
}

std::vector<int> convex_minlp(std::vector<int> &init_column, std::vector<Point> &P, std::vector<OrderedDistances> &od, std::vector<double> &dual_c, double dual_k, std::vector<std::pair<int, int>> &ml_pairs, std::vector<std::pair<int, int>> &cl_pairs) {
    int n = (int) P.size();
    PricingMINLP epc = compute_pricing_constants(P, od);
    std::vector<int> column;
    GRBModel *m = build_branching_problem_ptr(P, epc, dual_k, dual_c, ml_pairs, cl_pairs);
    m->update();

    Point y = {};
    int col_size = (int) init_column.size();
    std::vector<bool> a(n, false);
    for (auto &elem: init_column) {
        y.x += P[elem].x;
        y.y += P[elem].y;
        a[elem] = true;
    }
    y.x = y.x / col_size;
    y.y = y.y / col_size;

    m->getVarByName("y0").set(GRB_DoubleAttr_Start, y.x);
    m->getVarByName("y1").set(GRB_DoubleAttr_Start, y.y);

    for (int i = 0; i < n; i++) {
        if (a[i]) m->getVarByName("v" + std::to_string(i)).set(GRB_DoubleAttr_Start, 1.0);
        else m->getVarByName("v" + std::to_string(i)).set(GRB_DoubleAttr_Start, 0.0);
    }

    m->optimize();
    double rc = m->get(GRB_DoubleAttr_ObjVal);
    std::cout << "  \nCONVEX MINLP RC: " << rc << "\n\n";
    for (int i = 0; i < (int) P.size(); i++) {
        double v = m->getVarByName("v" + std::to_string(i)).get(GRB_DoubleAttr_X);
        if (v == 1.0)
            column.push_back(i);
    }
    delete (m);
    return column;
}
*/