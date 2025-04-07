#ifndef CLUSTERING_C_NODE_H
#define CLUSTERING_C_NODE_H

#include <vector>
#include "CGStruct.h"


struct Node {

    // lower bound
    double lb;
    // node id
    int id;
    // result
    CGResult *cg_result;

};

struct JobData {

    int type;
    std::pair<int, int> branching_pair;
    Node *node;

};

#endif //CLUSTERING_C_NODE_H
