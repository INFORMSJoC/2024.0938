#include <algorithm>
#include <iostream>
#include "JobQueue.h"

bool JobAbstractQueue::JobDataComparator(JobData *a, JobData *b) {
    return a->node->lb < b->node->lb;
}

JobData *JobAbstractQueue::pop() {
    JobData *node = queue.front();
    queue.pop_front();
    return node;
}

bool JobAbstractQueue::empty() {
    return queue.empty();
}

void JobAbstractQueue::sort() {
    std::sort(queue.begin(), queue.end(), JobDataComparator);
}

Node *JobAbstractQueue::getMinLB() {
    double min_lb = std::numeric_limits<double>::infinity();
    Node *min_node = nullptr;
    for (auto current : queue) {
        if (current->node->lb < min_lb) {
            min_lb = current->node->lb;
            min_node = current->node;
        }
    }
    return min_node;
}

Node *JobAbstractQueue::getMaxLB() {
    double max_lb = -std::numeric_limits<double>::infinity();
    Node *max_node = nullptr;
    for (auto current : queue) {
        if (current->node->lb > max_lb) {
            max_lb = current->node->lb;
            max_node = current->node;
        }
    }
    return max_node;
}

void JobAbstractQueue::print() {
    for (auto &elem : queue) {
        std::cout << elem->node->lb << " ";
    }
    std::cout << "\n";
}

int JobAbstractQueue::getSize() {
    return queue.size();
}

void JobQueue::push(JobData *node) {
    queue.push_back(node);
}

void JobStack::push(JobData *node) {
    queue.push_front(node);
}

void JobPriorityQueue::push(JobData *node) {
    queue.push_back(node);
    this->sort();
}
