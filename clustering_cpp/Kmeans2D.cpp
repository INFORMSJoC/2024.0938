#include <random>
#include <fstream>
#include "Kmeans2D.h"

#define KMEANS_SEED 12345


// Calculate the Euclidean distance between two points
double distance(const Point& p1, const Point& p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    return std::sqrt(dx * dx + dy * dy);
}

// Compute the total sum of squared errors (SSE) for a given set of points and centroids
double computeSSE(const std::vector<Point> &points, const std::vector<Point> &centroids, const std::vector<int>& assignments) {
    int n = (int) points.size();
    double sse = 0.0;
    for (int i = 0; i < n; i++) {
        double dist = distance(points[i], centroids[assignments[i]]);
        sse += dist * dist;
    }
    return sse;
}

// Compute the partial sum of squared errors (SSE) for a given cluster 'j',
double computePartialSSE(int j, const std::vector<Point> &points, const std::vector<int> &assignments) {

    int n = (int) points.size();

    // centroid of cluster j
    Point centroid = {0.0, 0.0};
    int count = 0;

    // data points in cluster j
    std::vector<Point> points_cluster;
    points_cluster.reserve(n);

    for (int i = 0; i < n; i++) {
        if (assignments[i] == j) {
            points_cluster.push_back(points[i]);
            centroid.x += points[i].x;
            centroid.y += points[i].y;
            count++;
        }
    }
    if (count > 0) {
        centroid.x = centroid.x / count;
        centroid.y = centroid.y / count;
    }

    // compute the cost of cluster j
    double sse = 0.0;
    for (auto &point : points_cluster) {
        double dist = distance(point, centroid);
        sse += dist * dist;
    }

    return sse;
}

// Initialize the centroids using k-means++ initialization
std::vector<Point> kmeansPP(const std::vector<Point>& points, int k, bool fixed_seed) {
    int n = (int) points.size();
    std::vector<Point> centroids;
    centroids.reserve(k);
    // Select the first centroid randomly from the points
    std::random_device rd;
    std::mt19937 gen(rd());
    if (fixed_seed)
        gen.seed(KMEANS_SEED);
    std::uniform_int_distribution<> distribution(0, n - 1);
    int idx = distribution(gen);
    centroids.push_back(points[idx]);
    // Select the remaining k-1 centroids using k-means++ initialization
    for (int i = 1; i < k; i++) {
        // Compute the distance to the closest centroid for each point
        std::vector<double> dists(n, std::numeric_limits<double>::max());
        for (int j = 0; j < i; j++) {
            for (int l = 0; l < n; l++) {
                double dist = distance(points[l], centroids[j]);
                dists[l] = std::min(dists[l], dist);
            }
        }
        // Select the next centroid with probability proportional to the distance squared
        double totalDist = 0.0;
        for (int j = 0; j < n; j++) {
            totalDist += dists[j] * dists[j];
        }
        std::uniform_real_distribution<double> dist(0.0, totalDist);
        double targetDist = dist(gen);
        double sumDist = 0.0;
        for (int j = 0; j < n; j++) {
            sumDist += dists[j] * dists[j];
            if (sumDist >= targetDist) {
                centroids.push_back(points[j]);
                break;
            }
        }
    }
    return centroids;
}

// Compute k-means clustering for a given set of points and centroids
std::vector<int> kmeans(const std::vector<Point>& points, std::vector<Point>& centroids) {
    int n = (int) points.size();
    int k = (int) centroids.size();
    std::vector<int> assignments(n, -1);
    bool changed = true;
    while (changed) {
        changed = false;
        for (int i = 0; i < n; i++) {
            double minDist = std::numeric_limits<double>::max();
            int minIndex = -1;
            for (int j = 0; j < k; j++) {
                double dist = distance(points[i], centroids[j]);
                if (dist < minDist) {
                    minDist = dist;
                    minIndex = j;
                }
            }
            if (assignments[i] != minIndex) {
                assignments[i] = minIndex;
                changed = true;
            }
        }
        for (int i = 0; i < k; i++) {
            int count = 0;
            double sumX = 0.0;
            double sumY = 0.0;
            for (int j = 0; j < n; j++) {
                if (assignments[j] == i) {
                    count++;
                    sumX += points[j].x;
                    sumY += points[j].y;
                }
            }
            if (count > 0) {
                centroids[i].x = sumX / count;
                centroids[i].y = sumY / count;
            }
        }
    }
    return assignments;
}


// Run k-means multiple times with different initial centroids and return the best clustering
std::tuple<double, std::vector<int>, std::vector<Point>> kmeans2D(const std::vector<Point>& points, int k, int numTrials) {
    bool fixed_seed = false;
    if (numTrials == 1)
        fixed_seed = true;
    std::vector<int> bestAssignments(points.size(), -1);
    std::vector<Point> bestCentroids;
    double bestSSE = std::numeric_limits<double>::max();
    for (int i = 0; i < numTrials; i++) {
        std::vector<Point> centroids = kmeansPP(points, k, fixed_seed);
        std::vector<int> assignments = kmeans(points, centroids);
        double sse = computeSSE(points, centroids, assignments);
        if (sse < bestSSE) {
            std::cout << "Trial: " << i << "\t MSSC objective: " << sse << "\n";
            bestAssignments = assignments;
            bestCentroids = centroids;
            bestSSE = sse;
        }
    }
    return std::make_tuple(bestSSE, bestAssignments, bestCentroids);
}

// Run k-means multiple times with different initial centroids and return the best clustering
std::tuple<double, std::vector<int>, std::vector<Point>> kmeans2D(const std::vector<Point>& points, int k, int numTrials, std::string &instance_name) {
    bool fixed_seed = false;
    if (numTrials == 1)
        fixed_seed = true;
    std::vector<int> bestAssignments(points.size(), -1);
    std::vector<Point> bestCentroids;
    double bestSSE = std::numeric_limits<double>::max();
    for (int i = 0; i < numTrials; i++) {
        std::vector<Point> centroids = kmeansPP(points, k, fixed_seed);
        std::vector<int> assignments = kmeans(points, centroids);
        double sse = computeSSE(points, centroids, assignments);
        if (sse < bestSSE) {
            std::cout << "Trial: " << i << "\t MSSC objective: " << sse << "\n";
            std::string out_file = "./kmeans_large_c++/" + instance_name + "_" + std::to_string(k) + ".txt";
            std::ofstream outFile(out_file);
            for (const auto &e : assignments) outFile << e << "\n";
            outFile.close();
            bestAssignments = assignments;
            bestCentroids = centroids;
            bestSSE = sse;
        }
    }
    return std::make_tuple(bestSSE, bestAssignments, bestCentroids);
}

std::tuple<double, std::vector<int>, std::vector<Point>> evaluate_clustering(const std::vector<Point> &points, int k, std::vector<int> &assignments) {

    int n = (int) points.size();
    std::vector<Point> centroids(k);

    // compute the cluster centers from the initial assignment vector
    for (int i = 0; i < k; i++) {
        int count = 0;
        double sumX = 0.0;
        double sumY = 0.0;
        for (int j = 0; j < n; j++) {
            if (assignments[j] == i) {
                count++;
                sumX += points[j].x;
                sumY += points[j].y;
            }
        }
        if (count > 0) {
            centroids[i].x = sumX / count;
            centroids[i].y = sumY / count;
        }
    }

    double sse = computeSSE(points, centroids, assignments);
    return std::make_tuple(sse, assignments, centroids);


}

std::vector<std::pair<int, double>> point_distance_from_centroids(Point &point, std::vector<Point> &centroids) {
    int k = (int) centroids.size();
    std::vector<std::pair<int, double>> closest(k); // vector or pairs (cluster_id, distance_from_cluster_id)
    for (int j = 0; j < k; j++) {
        double d = distance(point, centroids[j]);
        closest[j] = std::make_pair(j, d);
    }
    std::sort(closest.begin(), closest.end(), [](auto &left, auto &right) { return left.second < right.second; });
    return closest;
}



ClusteringResult get_cluster_assignment_2D(std::vector<Point> &P, int k, int n_start) {
    ClusteringResult cr;
    int n = (int) P.size();
    std::tuple<double, std::vector<int>, std::vector<Point>> t = kmeans2D(P, k, n_start);
    cr.loss = std::get<0>(t);
    std::vector<std::vector<iPoint>> assigned_points(k);
    for (int j = 0; j < k; j++) {
        std::vector<iPoint> points_cluster_j;
        for (int i = 0; i < n; i++) {
            if (std::get<1>(t)[i] == j) {
                points_cluster_j.push_back({i, P[i]});
            }
        }
        assigned_points[j] = points_cluster_j;
    }
    cr.assignment = std::get<1>(t);
    cr.centroid = std::get<2>(t);
    cr.points_cluster = assigned_points;
    return cr;
}

ClusteringResult *get_cluster_assignment_2D_ptr(std::vector<Point> &P, int k, int n_start) {
    ClusteringResult *cr = new ClusteringResult();
    int n = (int) P.size();
    std::tuple<double, std::vector<int>, std::vector<Point>> t = kmeans2D(P, k, n_start);
    cr->loss = std::get<0>(t);
    std::vector<std::vector<iPoint>> assigned_points(k);
    for (int j = 0; j < k; j++) {
        std::vector<iPoint> points_cluster_j;
        for (int i = 0; i < n; i++) {
            if (std::get<1>(t)[i] == j) {
                points_cluster_j.push_back({i, P[i]});
            }
        }
        assigned_points[j] = points_cluster_j;
    }
    cr->assignment = std::get<1>(t);
    cr->centroid = std::get<2>(t);
    cr->points_cluster = assigned_points;
    return cr;
}

ClusteringResult *get_cluster_assignment_2D_ptr(std::vector<Point> &P, int k, int n_start, std::string &instance_name) {
    ClusteringResult *cr = new ClusteringResult();
    int n = (int) P.size();
    std::tuple<double, std::vector<int>, std::vector<Point>> t = kmeans2D(P, k, n_start, instance_name);
    cr->loss = std::get<0>(t);
    std::vector<std::vector<iPoint>> assigned_points(k);
    for (int j = 0; j < k; j++) {
        std::vector<iPoint> points_cluster_j;
        for (int i = 0; i < n; i++) {
            if (std::get<1>(t)[i] == j) {
                points_cluster_j.push_back({i, P[i]});
            }
        }
        assigned_points[j] = points_cluster_j;
    }
    cr->assignment = std::get<1>(t);
    cr->centroid = std::get<2>(t);
    cr->points_cluster = assigned_points;
    return cr;
}

ClusteringResult evaluate_cluster_assignment_2D(std::vector<Point> &P, int k, std::vector<int> &target_assignment) {
    ClusteringResult cr;
    int n = (int) P.size();
    std::tuple<double, std::vector<int>, std::vector<Point>> t = evaluate_clustering(P, k, target_assignment);
    cr.loss = std::get<0>(t);
    std::vector<std::vector<iPoint>> assigned_points(k);
    for (int j = 0; j < k; j++) {
        std::vector<iPoint> points_cluster_j;
        for (int i = 0; i < n; i++) {
            if (std::get<1>(t)[i] == j) {
                points_cluster_j.push_back({i, P[i]});
            }
        }
        assigned_points[j] = points_cluster_j;
    }
    cr.assignment = target_assignment;
    cr.centroid = std::get<2>(t);
    cr.points_cluster = assigned_points;
    return cr;
}

ClusteringResult *evaluate_cluster_assignment_2D_ptr(std::vector<Point> &P, int k, std::vector<int> &target_assignment) {
    ClusteringResult *cr = new ClusteringResult();
    int n = (int) P.size();
    std::tuple<double, std::vector<int>, std::vector<Point>> t = evaluate_clustering(P, k, target_assignment);
    cr->loss = std::get<0>(t);
    std::vector<std::vector<iPoint>> assigned_points(k);
    for (int j = 0; j < k; j++) {
        std::vector<iPoint> points_cluster_j;
        for (int i = 0; i < n; i++) {
            if (std::get<1>(t)[i] == j) {
                points_cluster_j.push_back({i, P[i]});
            }
        }
        assigned_points[j] = points_cluster_j;
    }
    cr->assignment = target_assignment;
    cr->centroid = std::get<2>(t);
    cr->points_cluster = assigned_points;
    return cr;
}

/*
std::pair<std::vector<Point>, std::vector<std::vector<iPoint>>> get_cluster_assignment_2D(std::vector<Point> &P, int k, double &loss) {
    int n = (int) P.size();
    std::tuple<double, std::vector<int>, std::vector<Point>> t = kmeans2D(P, k, 1);
    loss = std::get<0>(t);
    std::vector<std::vector<iPoint>> assigned_points(k);
    for (int j = 0; j < k; j++) {
        std::vector<iPoint> points_cluster_j;
        for (int i = 0; i < n; i++) {
            if (std::get<1>(t)[i] == j) {
                points_cluster_j.push_back({i, P[i]});
            }
        }
        assigned_points[j] = points_cluster_j;
    }

    return std::make_pair(std::get<2>(t), assigned_points);

}

std::pair<std::vector<Point>, std::vector<std::vector<iPoint>>> evaluate_cluster_assignment_2D(std::vector<Point> &P, int k, double &loss, std::vector<int> &target_assignment) {
    int n = (int) P.size();
    std::tuple<double, std::vector<int>, std::vector<Point>> t = evaluate_clustering(P, k, target_assignment);
    loss = std::get<0>(t);
    std::vector<std::vector<iPoint>> assigned_points(k);
    for (int j = 0; j < k; j++) {
        std::vector<iPoint> points_cluster_j;
        for (int i = 0; i < n; i++) {
            if (std::get<1>(t)[i] == j) {
                points_cluster_j.push_back({i, P[i]});
            }
        }
        assigned_points[j] = points_cluster_j;
    }

    return std::make_pair(std::get<2>(t), assigned_points);

}
*/