#ifndef CLUSTERING_C_THREADPOOL_H
#define CLUSTERING_C_THREADPOOL_H

#include <thread>
#include "branch_and_price.h"
#include "CGStruct.h"

class ThreadPool {

private:

    CGConfig *cg_config;
    BBConfig *bb_config;
    CGData *input_data;
    SharedData  *shared_data;

    // We store the threads in a vector, so we can later stop them gracefully
    std::vector<std::thread> threads;

    // This will be set to true when the thread pool is shutting down. This tells
    // the threads to stop looping and finish
    bool done;

    void doWork(int id);


public:

    ThreadPool(CGConfig *cg_config, BBConfig *bb_config, SharedData *shared_data, CGData *input_data, int n_thread);
    void quitPool();
    void addJob(JobData *job_data);

};


#endif //CLUSTERING_C_THREADPOOL_H
