#include <unordered_set>
#include <algorithm>
#include "cg_util.h"

void add_both(std::map<int, std::set<int>> &graph, int i, int j) {
    graph[i].insert(j);
    graph[j].insert(i);
}

void dfs(int i, std::map<int, std::set<int>> &graph, std::vector<bool> &visited, std::vector<int> &component) {
    visited[i] = true;
    for (auto &j : graph[i]) {
        if (!visited[j]) {
            dfs(j, graph, visited, component);
        }
    }
    component.push_back(i);
}

std::unordered_set<int> find_item(std::vector<int> &v, int id) {
    std::unordered_set<int> comp;
    for (int i = 0; i < (int) v.size(); i++) {
        if (v[i] == id)
            comp.insert(i);
    }
    return comp;
}

MustLinkMapping get_ml_map(std::vector<int> &assignment) {

    int n = (int) assignment.size();
    int max_id = *max_element(assignment.begin(), assignment.end());

    std::vector<std::unordered_set<int>> ml_comp; // mapping
    std::vector<int> ml_index(n); // index

    for (int i = 0; i < max_id + 1; i++) {
        std::unordered_set<int> comp = find_item(assignment, i);
        ml_comp.push_back(comp);
    }

    int n_components = (int) ml_comp.size();
    for (int i = 0; i < n_components; i++) {
        for (const auto &elem : ml_comp[i]) {
            ml_index[elem] = i;
        }
    }


    return {ml_comp, ml_index};

}

void split_components(std::unordered_set<int> &comp, std::vector<int> &kmeans_assignment, std::vector<std::unordered_set<int>> &comp_list) {
    int max_id = *max_element(kmeans_assignment.begin(), kmeans_assignment.end());
    for (int i = 0; i < max_id + 1; i++) {
        std::unordered_set<int> comp_local = find_item(kmeans_assignment, i);
        std::unordered_set<int> inter;
        for (int element : comp) {
            if (comp_local.count(element) > 0) {
                inter.insert(element);
            }
        }
        if (!inter.empty())
            comp_list.push_back(inter);
    }
}

MustLinkMapping get_ml_map(std::vector<int> &hierarchical_assignment, std::vector<int> &kmeans_assignment) {

    int n = (int) hierarchical_assignment.size();
    int max_id = *max_element(hierarchical_assignment.begin(), hierarchical_assignment.end());

    std::vector<std::unordered_set<int>> ml_comp; // mapping
    std::vector<int> ml_index(n); // index

    for (int i = 0; i < max_id + 1; i++) {
        std::unordered_set<int> comp = find_item(hierarchical_assignment, i);
        split_components(comp, kmeans_assignment, ml_comp);
    }

    int n_components = (int) ml_comp.size();
    for (int i = 0; i < n_components; i++) {
        for (const auto &elem : ml_comp[i]) {
            ml_index[elem] = i;
        }
    }

    return {ml_comp, ml_index};
}

MustLinkMapping get_ml_map(int n, std::vector<std::pair<int, int>> &ml) {

    std::map<int, std::set<int>> ml_graph; // the graph induced by must-link constraints
    std::vector<std::unordered_set<int>> ml_comp; // mapping
    std::vector<int> ml_index(n); // index

    for (int i = 0; i < n; i++) {
        ml_graph.insert(std::pair<int, std::set<int>> (i, {}));
    }

    for (auto &pair_ml : ml) {
        add_both(ml_graph, pair_ml.first, pair_ml.second);
    }

    std::vector<bool> visited(n, false);
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            std::vector<int> component;
            dfs(i, ml_graph, visited, component);
            std::unordered_set<int> component_set(component.begin(), component.end());
            ml_comp.push_back(component_set);
        }
    }

    int n_components = (int) ml_comp.size();
    for (int i = 0; i < n_components; i++) {
        for (const auto &elem : ml_comp[i]) {
            ml_index[elem] = i;
        }
    }

    return {ml_comp, ml_index};
}

void display_graph(MustLinkMapping &mlm) {

    // display ml_map
    int n_comp = (int) mlm.ml_comp.size();
    for (int i = 0; i < n_comp; i++) {
        std::unordered_set<int> comp = mlm.ml_comp[i];
        if (comp.empty())
            continue;
        std::printf("%d: ", i);
        std::printf("{");
        for (auto &set_elem : comp) {
            std::printf(" %d ", set_elem);
        }
        std::printf("}\n");
    }

    // display ml_index
    for (auto &elem : mlm.ml_index) {
        std::cout << elem << " ";
    }
    std::cout << "\n";

}

// Compare two distances based on their values
bool compareDistances(const std::pair<int, double>& a, const std::pair<int, double>& b) {
    return a.second < b.second;
}

// bigM
double get_cost_clustered_together(std::vector<Point> &P) {
    Point centroid = {0.0, 0.0};
    int n = (int) P.size();
    for (int i = 0; i < n; i++) {
        centroid.x += P[i].x;
        centroid.y += P[i].y;
    }
    centroid.x = centroid.x / n;
    centroid.y = centroid.y / n;
    double cost = 0.0;
    for (int i = 0; i < n; i++) {
        cost += std::pow(point_distance(P[i], centroid), 2);
    }
    return cost;
}

std::vector<std::vector<double>> get_distance_matrix(std::vector<Point> &P) {

    int n = (int) P.size();
    std::vector<std::vector<double>> D(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double dist = point_distance(P[i], P[j]);
            D[i].push_back(dist * dist);
        }
    }

    return D;
}

// Compute sorted distances with the remaining points for each point
std::vector<OrderedDistances> get_ordered_distances(std::vector<Point>& points) {
    int n = (int) points.size();
    std::vector<OrderedDistances> orderedDistances(n);
    for (int i = 0; i < n; i++) {
        std::vector<std::pair<int, double>> distances;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                double dist = point_distance(points[i], points[j]);
                distances.emplace_back(j, dist);
            }
        }
        std::sort(distances.begin(), distances.end(), compareDistances);
        orderedDistances[i].dist.resize(n-1);
        orderedDistances[i].index.resize(n-1);
        for (int j = 0; j < n-1; j++) {
            orderedDistances[i].dist[j] = distances[j].second;
            orderedDistances[i].index[j] = distances[j].first;
        }
    }
    return orderedDistances;
}


double point_distance(Point &p1, Point &p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return std::sqrt((dx * dx) + (dy * dy));
}

int circle_circle_intersection(double x0, double y0, double r0,
                               double x1, double y1, double r1,
                               double &x1_out, double &y1_out,
                               double &x2_out, double &y2_out)
{
    double a, dx, dy, d, h, rx, ry;
    double x2, y2;

    /* dx and dy are the vertical and horizontal distances between the circle centers. */
    dx = x1 - x0;
    dy = y1 - y0;

    /* Determine the straight-line distance between the centers. */
    d = hypot(dx, dy); // Suggested by Keith Briggs

    if (d <= fabs(r0 - r1)) {
        /* no solution. one circle is contained in the other */
        return 0;
    }

    /* Determine the distance from point 0 to point 2. */
    a = ((r0*r0) - (r1*r1) + (d*d)) / (2.0 * d) ;

    /* Determine the coordinates of point 2. */
    x2 = x0 + (dx * a/d);
    y2 = y0 + (dy * a/d);

    /* Determine the distance from point 2 to either of the intersection points. */
    h = sqrt((r0*r0) - (a*a));

    /* Now determine the offsets of the intersection points from * point 2. */
    rx = -dy * (h/d);
    ry = dx * (h/d);

    /* Determine the absolute intersection points. */
    x1_out = x2 + rx;
    x2_out = x2 - rx;
    y1_out = y2 + ry;
    y2_out = y2 - ry;

    if(x1_out == x2_out && y1_out == y2_out)  {
        /* there is only one intersection point */
        return 1;
    }

    return 2;
}

int intersect_circles(Circle &circle1, Circle &circle2, std::vector<Point> &Q) {

    double x1, y1, x2, y2;
    int n_intersect = circle_circle_intersection(circle1.center.x, circle1.center.y, circle1.radius,
                                                 circle2.center.x, circle2.center.y, circle2.radius,
                                                 x1, y1, x2, y2);

    Q[0].x = x1;
    Q[0].y = y1;
    Q[1].x = x2;
    Q[1].y = y2;

    return n_intersect;

}

// The radius of a cluster is the maximum distance between all the points and the centroid.
double compute_cluster_radius(Point &centroid, std::vector<iPoint> &points) {

    double max_dist = 0.0;
    for (auto &p : points) {
        double dist = point_distance(p.point, centroid);
        if (dist > max_dist) {
            max_dist = dist;
        }
    }
    return max_dist;

}

std::vector<std::pair<int, int>> candidate_must_link_radius(int k, double perc, std::vector<Point> &centroids, std::vector<std::vector<iPoint>> &points) {

    std::vector<std::pair<int, int>> ml_pairs;
    std::vector<std::vector<iPoint>> ml_points(k); // ml points for each cluster
    for (int j = 0; j < k; j++) {
        double radius_j = compute_cluster_radius(centroids[j], points[j]);
        // std::cout << "Radius cluster " << j << ": " << radius_j << "\n";
        for (auto &p : points[j]) {
            if (point_distance(p.point, centroids[j]) <= perc * radius_j) {
                ml_points[j].push_back(p);
            }
        }
    }
    // generate ml pairs by looking at the selected points within the radius of each cluster
    for (int j = 0; j < k ; j++) {
        std::vector<iPoint> points_j = ml_points[j];
        int n_points_j = (int) points_j.size();
        // std::cout << "Selected points in cluster " << j << ": " << n_points_j << "\n";
        for (int i = 0; i < n_points_j; i++) {
            for (int h = i + 1; h < n_points_j; h++) {
                ml_pairs.emplace_back(points_j[i].id, points_j[h].id);
            }
        }
    }
    return ml_pairs;

}


/********************************************** START - DUAL VARIABLES DISAGGREGATION *****************************************/

std::vector<double> disaggregate_dual_variables_uniform(MustLinkMapping &mlm, std::vector<double> &dual_c_shr) {

    int n = (int) mlm.ml_index.size();
    std::vector<double> dual_c(n, 0); // disaggregated dual variables

    for (int i = 0; i < n; i++) {
        int j = mlm.ml_index[i];
        if (dual_c_shr[j] == 0)
            dual_c[i] = 0;
        else
            dual_c[i] = dual_c_shr[j] / (double) mlm.ml_comp[j].size(); // get the size of the component
    }

    return dual_c;
}

// choose index within the set of the component to set equal to the dual variable
std::vector<double> disaggregate_dual_variables_sparse(MustLinkMapping &mlm, std::vector<double> &dual_c_shr, int index) {

    int n = (int) mlm.ml_index.size();
    std::vector<double> dual_c(n, 0.0); // disaggregated dual variables
    int n_comp = (int) dual_c_shr.size();

    int idx;
    for (int i = 0; i < n_comp; i++) {
        int comp_size = (int) mlm.ml_comp[i].size();
        if (index < comp_size)
            idx = *std::next(mlm.ml_comp[i].begin(), index);
        else
            idx = *(mlm.ml_comp[i].begin());
        for (auto &j : mlm.ml_comp[i]) {
            if (dual_c_shr[i] != 0) {
                if (j == idx)
                    dual_c[j] = dual_c_shr[i] - (comp_size - 1) * TOL;
                else
                    dual_c[j] = TOL;
            }
        }
    }
    return dual_c;

}

double cost(std::vector<Point> &P, std::vector<bool> &a) {

    int n = (int) P.size();

    int card = 0;
    Point centroid = {0.0, 0.0};
    for (int i = 0; i < n; i++) {
        if (a[i]) {
            centroid.x += P[i].x;
            centroid.y += P[i].y;
            card++;
        }
    }
    centroid.x = centroid.x / card;
    centroid.y = centroid.y / card;

    double cost = 0;
    for (int i = 0; i < n; i++)
        if (a[i]) {
            cost += std::pow(point_distance(P[i], centroid), 2);
        }

    return cost;
}

std::vector<double> disaggregate_dual_variables_complementarity_problem(std::vector<Point> &P, MustLinkMapping &mlm,
                                                                        std::vector<std::vector<int>> &incompatible_cols,
                                                                        double aggr_dual_k, std::vector<double> &aggr_dual_c) {

    // P: n x d matrix of data points
    // ml_map: mapping between points ad aggregated points

    int n = (int) mlm.ml_index.size();
    int n_aggregated = (int) mlm.ml_comp.size();
    int n_incompatible = (int) incompatible_cols.size();

    std::vector<double> disaggregated_dual(n);

    auto env = GRBEnv(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    auto problem = new GRBModel(env);

    // Variables of the complementary problem
    std::vector<GRBVar> lambda(n);
    for (int i = 0; i < n; i++) {
        lambda[i] = problem->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "l" + std::to_string(i));
    }

    GRBVar z = problem->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "z");

    for (int t = 0; t < n_incompatible; t++) {
        std::vector<bool> incompatible_A(n, false);
        for (auto &elem : incompatible_cols[t])
            incompatible_A[elem] = true;
        GRBLinExpr rhs_sum = 0.0;
        for (int i = 0; i < n; i++) {
            if (incompatible_A[i])
                rhs_sum += lambda[i];
        }
        problem->addConstr(z <= cost(P, incompatible_A) + aggr_dual_k - rhs_sum);
    }


    for (int i = 0; i < n_aggregated; i++) {
        GRBLinExpr lhs = 0.0;
        for (auto &elem : mlm.ml_comp[i]) {
            lhs += lambda[elem];
        }
        problem->addConstr(lhs == aggr_dual_c[i]);
    }

    GRBLinExpr obj = z;
    problem->setObjective(obj, GRB_MAXIMIZE);

    problem->set("OutputFlag", "0");
    problem->set("Method", "1"); // -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier

    problem->optimize();
    GRBVar *vars = problem->getVars();
    for (int i = 0; i < n; i++) {
        disaggregated_dual[i] = vars[i].get(GRB_DoubleAttr_X);
    }

    delete (problem);

    return disaggregated_dual;
}


/********************************************** END - DUAL VARIABLES DISAGGREGATION *****************************************/


/*
 KEY FUNCTION: IT IS CALLED WHEN:
 i) A NEW (DISAGGREGATED) COLUMN IS GENERATED AT EACH ITERATION
 ii) OLD COLUMNS MUST BE CONVERTED WHEN THE PROBLEM CHANGES IN THE DCA PHASE
*/
std::vector<bool> aggregate_column(MustLinkMapping &mlm, std::vector<int> &column) {

    int col_size = (int) column.size();
    std::vector<bool> a(mlm.ml_comp.size(), false);

    for (int i = 0; i < col_size; i++) {
        // select the component where column[i] is present
        int j = mlm.ml_index[column[i]];
        a[j] = true;
    }

    return a;

}


// KEY FUNCTION: IT IS CALLED SEVERAL TIMES IN THE PRICING PROBLEM
bool isCompatible(MustLinkMapping &mlm, std::vector<int> &column) {

    int col_size = (int) column.size();
    std::vector<bool> examined(mlm.ml_comp.size(), false);

    for (int i = 0; i < col_size; i++) {
        // select the component where column[i] is present
        int j = mlm.ml_index[column[i]];
        int comp_size = (int) mlm.ml_comp[j].size();
        if (comp_size != 1 && !examined[j]) {
            examined[j] = true;
            // Go through all items of the column, checking that they are present in the component
            int n_inter = (int) count_if(column.begin(), column.end(),
                                         [&](int k) {return mlm.ml_comp[j].find(k) != mlm.ml_comp[j].end();});
            if (n_inter != comp_size) {
                // column is not compatible
                return false;
            }
        }
    }
    return true;
}


std::vector<std::vector<bool>> update_columns(MustLinkMapping &mlm, std::vector<std::vector<int>> &or_cols) {

    // Columns that were compatible with the old aggregation are still compatible when the old aggregation is decomposed
    // thus, no need to check compatibility. Update indices in old columns accordingly

    int n_cols = (int) or_cols.size();
    std::vector<std::vector<bool>> new_aggr_cols_bool(n_cols);
    for (int i = 0; i < n_cols; i++)
        new_aggr_cols_bool[i] = aggregate_column(mlm, or_cols[i]);

    return new_aggr_cols_bool;
}

std::vector<std::vector<bool>> get_initial_columns_aggregation(std::vector<std::vector<iPoint>> &points, MustLinkMapping &mlm, std::vector<std::vector<int>> &or_cols) {

    int k = (int) points.size();

    for (int j = 0; j < k; j++) {
        // the j-th cluster represents a single column
        std::vector<int> temp_col;
        for (auto &p : points[j]) {
            temp_col.push_back(p.id);
        }
        or_cols.push_back(temp_col);
    }

    // Update indices inside the column according to the current aggregation
    std::vector<std::vector<bool>> A = update_columns(mlm, or_cols);

    return A;
}

int get_number_of_incompatibilities(MustLinkMapping &mlm, std::vector<int> &column) {

    int count = 0;
    int n_comp = (int) mlm.ml_comp.size();

    for (int i = 0; i < n_comp; i++) {
        int comp_size = (int) mlm.ml_comp[i].size();
        if (comp_size != 1) {
            // Go through all items of the column, checking that they are present in the component
            int n_inter = (int) count_if(column.begin(), column.end(),
                                         [&](int k) {return mlm.ml_comp[i].find(k) != mlm.ml_comp[i].end();});
            if (n_inter != comp_size && n_inter != 0) {
                // column is not compatible
                count++;
            }
        }
    }
    return count;
}


void update_aggregation_improved(MustLinkMapping &mlm, std::vector<int> &column) {

    /*
     * for each i in incompatible column
     - find the component j where i belongs to and mark component j as visited
     - if component has size == 1 or component j has been visited then skip that element i
     - else find the intersection between component j and the whole column
      - if the number of common elements is equal to the size of the component then i is not a source of incompatibility
      - else i is the source of incompatibility and needs to be removed from the component j
      - actually consider the set diff e put it as a new component.
     */

    int col_size = (int) column.size();
    for (int i = 0; i < col_size; i++) {
        // select the component where column[i] is present
        int j = mlm.ml_index[column[i]];
        int comp_size = (int) mlm.ml_comp[j].size();
        if (comp_size > 1) {
            // Go through all items of the column, checking that they are present in the component
            int n_inter = (int) count_if(column.begin(), column.end(),
                                         [&](int k) {return mlm.ml_comp[j].find(k) != mlm.ml_comp[j].end();});
            if (n_inter != comp_size && n_inter != 0) {
                // column[i] is a source of incompatibility: add the set difference between the column and the component as a new component
                std::unordered_set<int> new_comp;
                for (const auto &elem : column) {
                    auto erase = mlm.ml_comp[j].erase(elem);
                    if (erase == 1) {
                        new_comp.insert(elem);
                    }
                }
                mlm.ml_comp.push_back(new_comp);
            }
        }
    }

    int comp_count = (int) mlm.ml_comp.size();
    for (int i = 0; i < comp_count; i++) {
        for (const auto &elem : mlm.ml_comp[i]) {
            mlm.ml_index[elem] = i;
        }
    }
}