#ifndef CLUSTERING_C_BRANCH_AND_PRICE_H
#define CLUSTERING_C_BRANCH_AND_PRICE_H

#include <condition_variable>
#include "JobQueue.h"
#include "cg_util.h"
#include "CGStruct.h"
#include "Kmeans2D.h"

#define BEST_FIRST 0
#define DEPTH_FIRST 1
#define BREADTH_FIRST 2

#define ROOT 0
#define MUST_LINK 1
#define CANNOT_LINK -1


struct SharedData {

    // Between workers and main
    std::condition_variable mainConditionVariable;
    std::vector<bool> threadStates;

    // Queue of requests waiting to be processed
    JobAbstractQueue *queue;
    // This condition variable is used for the threads to wait until there is work to do
    std::condition_variable queueConditionVariable;
    // Mutex to protect queue
    std::mutex queueMutex;

    double global_ub;
    std::vector<int> global_assignment;
    double gap;
    int n_nodes;

};

class BBLog {

public:

    int id; // node id
    int type; // node type
    int i; // first branching index
    int j; // second branching index
    double lb; // node lower bound
    int iter; // number of CG iterations
    int columns; // number of columns
    int n_covering_start; // initial size of the problem
    int n_covering_end; // final size of the problem
    double time; // elapsed time in seconds
    double ub; // upper bound at the node
    double gub; // global upper bound
    double gap; // gap at the node
    double ggap; // global gap
    int open; // number of open nodes

public:

    static void header();
    static void footer();
    void print() const;

};

std::vector<JobData *> solve_child_problem(JobData *job, CGData *cg_data, SharedData *shared_data, CGConfig *cg_config, BBConfig *bb_config);
BBResult *branch_and_bound(CGData *cg_data, ClusteringResult *init_cr, CGConfig *cg_config, BBConfig *bb_config);
#endif //CLUSTERING_C_BRANCH_AND_PRICE_H
