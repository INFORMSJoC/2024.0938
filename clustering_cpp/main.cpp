#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include "cg_stabilization.h"
#include <cstring>
#include <random>
#include "Kmeans2D.h"
#include "branch_and_price.h"

std::map<std::string, std::string> read_params(std::string &config_file) {

    std::map<std::string, std::string> config_map = {};

    std::ifstream cFile (config_file);
    if (cFile.is_open()) {
        std::string line;
        while (getline(cFile, line)){
            line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
            if(line[0] == '#' || line.empty())
                continue;
            auto delimiterPos = line.find('=');
            auto key = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);
            config_map.insert(std::pair<std::string, std::string>(key, value));
        }

    }
    else {
        std::cerr << "Couldn't open config file for reading.\n";
    }

    return config_map;
}

std::vector<Point> read_data(const char *filename, int &n, int &d) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << strerror(errno) << "\n";
        exit(EXIT_FAILURE);
    }
    // read the header n (number of points), d (dimension)
    file >> n >> d;
    if (d != 2) {
        std::cerr << "Error: this algorithm works only for data points in the plane" << "\n";
        exit(EXIT_FAILURE);
    }
    std::vector<Point> P(n);
    for (int i = 0; i < n; i++) {
        Point p = {0, 0};
        file >> p.x;
        file >> p.y;
        P[i] = p;
    }
    return P;
}

std::vector<Point> scale_data(std::vector<Point> &P) {
    int n = (int) P.size();
    std::vector<Point> P_scaled(n);
    double x_min = std::numeric_limits<double>::infinity();
    double x_max = -std::numeric_limits<double>::infinity();
    double y_min = std::numeric_limits<double>::infinity();
    double y_max = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; i++) {
        if (P[i].x < x_min)
            x_min = P[i].x;
        if (P[i].x > x_max)
            x_max = P[i].x;
        if (P[i].y < y_min)
            y_min = P[i].y;
        if (P[i].y > y_max)
            y_max = P[i].y;
    }

    for (int i = 0; i< n; i++) {
        P_scaled[i].x = (P[i].x - x_min)/(x_max - x_min);
        P_scaled[i].y = (P[i].y - y_min)/(y_max - y_min);
    }

    return P_scaled;

}

std::vector<int> read_assignment(const char *filename, int &n) {
    std::ifstream file(filename);
    if (!file) {
        //std::cerr << strerror(errno) << "\n";
        //exit(EXIT_FAILURE);
        return {};
    }
    std::vector<int> v(n);
    for (int i = 0; i < n; i++) {
        file >> v[i];
    }
    return v;
}

void run_cg(int argc, char **argv) {

    if (argc < 5) {
        std::cerr << "Usage: <DATASET> <K> <AGGREGATION_LEVEL> <CONFIG_FILE>" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string dname(argv[1]);
    std::string p = "./instances/" + dname + ".txt";

    const char *data_path = p.c_str();

    int n, d;
    std::vector<Point> P = read_data(data_path, n, d);

    int k = std::stoi(argv[2]);
    int n_comp = std::stoi(argv[3]);

    std::string cname = argv[4];
    std::string c = cname + ".txt";
    std::map<std::string, std::string> config_map = read_params(c);

    /*
    std::string csv_file = argv[5];
    csv_file = csv_file + ".csv";
    std::ofstream outfile;
    outfile.open(csv_file, std::ios_base::app); // append

    std::ifstream f(csv_file);
    if (f.peek() == std::ifstream::traits_type::eof()) {
        outfile << CGResult::to_csv_header() << "\n";
        outfile.flush();
    }
    f.close();
    */

    int n_init = std::stoi(config_map["KMEANS_INIT"]);
    std::cout << "\nK-MEANS INIT: " << n_init << "\n";


    auto *cg_config = new CGConfig();
    cg_config->n_components = n_comp;
    cg_config->time_limit = std::stod(config_map["CG_TIME_LIMIT"]);
    cg_config->max_iter = std::stoi(config_map["CG_MAX_ITER"]);
    cg_config->log_step = std::stoi(config_map["CG_LOG_STEP"]);
    cg_config->dual_disaggregation_strategy = std::stoi(config_map["CG_DUAL_DISAGGREGATION"]);
    cg_config->best_partition_update = std::stoi(config_map["CG_PARTITION_UPDATE"]) != 0;
    cg_config->disaggregation_threshold = std::stod(config_map["CG_DISAGGREGATION_THRESHOLD"]);
    cg_config->max_comp_cols = std::stoi(config_map["CG_MAX_COMP_COLS"]);
    cg_config->max_cols = std::stoi(config_map["CG_MAX_COLS"]);
    cg_config->cg_verbose = std::stoi(config_map["CG_VERBOSE"]) != 0;
    cg_config->cg_inherit_partition_child = std::stoi(config_map["CG_INHERIT_PARTITION_CHILD"]) != 0;
    cg_config->cg_heuristic_pricing_child = std::stoi(config_map["CG_HEURISTIC_PRICING_CHILD"]) != 0;
    cg_config->cg_mip_heuristic = std::stoi(config_map["CG_MIP_HEURISTIC"]) != 0;
    cg_config->max_comp_cols_child = std::stoi(config_map["CG_MAX_COMP_COLS_CHILD"]);

    auto *bb_config = new BBConfig();
    bb_config->instance_name = dname;
    bb_config->bb_tol = std::stod(config_map["BB_TOL"]);
    bb_config->bb_parallel = std::stoi(config_map["BB_THREADS"]);
    bb_config->bb_max_nodes = std::stoi(config_map["BB_MAX_NODES"]);
    bb_config->bb_visiting_strategy = std::stoi(config_map["BB_VISITING_STRATEGY"]);
    bb_config->bb_verbose = std::stoi(config_map["BB_VERBOSE"]) != 0;

    cg_config->ip_tol = bb_config->bb_tol;

    std::cout << "\n";
    std::cout << "DATA: " << dname << "\n";
    std::cout << "K: " << k << "\n";
    std::cout << "N_COMPONENTS: " << n_comp << "\n";
    std::cout << "\n";
    std::cout << "CG_TIME_LIMIT: " << cg_config->time_limit << "\n";
    std::cout << "CG_MAX_ITER: " << cg_config->max_iter << "\n";
    std::cout << "CG_LOG_STEP: " << cg_config->log_step << "\n";
    std::cout << "CG_DUAL_DISAGGREGATION: " << cg_config->dual_disaggregation_strategy << "\n";
    std::cout << "CG_PARTITION_UPDATE: " << cg_config->best_partition_update << "\n";
    std::cout << "CG_DISAGGREGATION_THRESHOLD: " << cg_config->disaggregation_threshold << "\n";
    std::cout << "CG_MAX_COMP_COLS: " << cg_config->max_comp_cols << "\n";
    std::cout << "CG_MAX_COLS: " << cg_config->max_cols << "\n";
    std::cout << "CG_MIP_HEURISTIC: " << cg_config->cg_mip_heuristic << "\n";
    std::cout << "CG_VERBOSE: " << cg_config->cg_verbose << "\n";
    std::cout << "CG_INHERIT_PARTITION_CHILD: " << cg_config->cg_inherit_partition_child << "\n";
    std::cout << "CG_HEURISTIC_PRICING_CHILD: " << cg_config->cg_heuristic_pricing_child << "\n";
    std::cout << "CG_MAX_COMP_COLS_CHILD: " << cg_config->max_comp_cols_child << "\n";
    std::cout << "\n";
    std::cout << "BB_TOL: " << bb_config->bb_tol << "\n";
    std::cout << "BB_THREADS: " << bb_config->bb_parallel << "\n";
    std::cout << "BB_MAX_NODES: " << bb_config->bb_max_nodes << "\n";
    std::cout << "BB_VISITING_STRATEGY: " << bb_config->bb_visiting_strategy << "\n";
    std::cout << "BB_VERBOSE: " << bb_config->bb_verbose << "\n";
    std::cout << "\n";



    std::vector<int> init_assignment;
    ClusteringResult *init_cr;

    if (n_init > 0) {

        std::string command = "python run_kmeans.py ";
        std::string out_assignment = "./initial_assignment/" + dname + "_" + std::to_string(k) + ".txt";
        std::string args = p + " " + std::to_string(k) + " " + std::to_string(n_init) + " " + out_assignment;
        command += args;
        std::cout << command << "\n";
        int system_val = system(command.c_str());
        if (system_val == -1) {
            // The system method failed
            std::cout << "Failed to call Python script" << "\n";
            exit(EXIT_FAILURE);
        }


        /*
        std::string command = "python run_kmeans_multistart.py ";
        std::string out_assignment = "./kmeans_large_c++/" + dname + "_" + std::to_string(k) + ".txt";
        std::string args = p + " " + std::to_string(k) + " " + std::to_string(n_init) + " " + out_assignment;
        command += args;
        std::cout << command << "\n";
        int system_val = system(command.c_str());
        if (system_val == -1) {
            // The system method failed
            std::cout << "Failed to call Python script" << "\n";
            exit(EXIT_FAILURE);
        }

        exit(EXIT_SUCCESS);
        */

        std::string s = "./initial_assignment/" + dname + "_" + std::to_string(k) + ".txt";
        init_assignment = read_assignment(s.c_str(), n);
        init_cr = evaluate_cluster_assignment_2D_ptr(P, k, init_assignment);
        std::cout << "Using kmeans (python) solution with value: " << init_cr->loss << "\n";

    } else {

        // use pre-computed assignment provided by the user
        std::string s = "./initial_assignment/" + dname + "_" + std::to_string(k) + ".txt";
        init_assignment = read_assignment(s.c_str(), n);
        init_cr = evaluate_cluster_assignment_2D_ptr(P, k, init_assignment);
        std::cout << "Using solution provided by the user with value: " << init_cr->loss << "\n";

    }

    std::vector<OrderedDistances> od = get_ordered_distances(P);
    std::vector<std::vector<double>> D = get_distance_matrix(P);

    std::vector<int> Q_assignment;

    if (n_comp == k) {

        Q_assignment = init_assignment;

    } else if (n_comp == n) {

        for (int i = 0; i < n; i++) {
            Q_assignment.push_back(i);
        }

    } else {

        std::string out_assignment = "./initial_assignment/" + dname + "_" + std::to_string(n_comp) + ".txt";

        if (n_init > 0) {
            std::string command = "python run_kmeans.py ";
            std::string args = p + " " + std::to_string(n_comp) + " " + std::to_string(n_init) + " " + out_assignment;
            command += args;
            std::cout << command << "\n";
            int system_val = system(command.c_str());
            if (system_val == -1) {
                // The system method failed
                std::cout << "Failed to call Python script" << "\n";
                exit(EXIT_FAILURE);
            }
        }

        Q_assignment = read_assignment(out_assignment.c_str(), n);

    }

    auto *cg_data = new CGData({P, k, od, D, Q_assignment});

    if (bb_config->bb_max_nodes == 1) {

        // run only the CG algorithm
        CGResult *r = solve_stabilized_master_problem_dca(cg_data, cg_config, init_cr);
        r->print();


        //outfile << r->to_csv_line() << "\n";
        //outfile.flush();

        double gap = std::abs(r->best_ub - r->best_lb) / std::abs(r->best_ub);
        std::string out_file = "./optimal_assignment/" + bb_config->instance_name + "_" + std::to_string(k) + ".txt";
        if (!r->cr_end->assignment.empty() && gap <= bb_config->bb_tol) {
            std::ofstream outFile(out_file);
            for (const auto &e : r->cr_end->assignment) outFile << e << "\n";
            outFile.close();
        }

        delete (r);
        delete (init_cr);

    } else {

        BBResult *b = branch_and_bound(cg_data, init_cr, cg_config, bb_config);

        std::string out_file = "./optimal_assignment/" + bb_config->instance_name + "_" + std::to_string(k) + ".txt";
        if (!b->opt_assignment.empty() && b->opt_gap <= bb_config->bb_tol) {
            std::ofstream outFile(out_file);
            for (const auto &e : b->opt_assignment) outFile << e << "\n";
            outFile.close();
        }

        delete (b);

    }

    delete (cg_data);
    delete (cg_config);
    delete (bb_config);

}

int main(int argc, char **argv) {

    run_cg(argc, argv);

    return EXIT_SUCCESS;
}