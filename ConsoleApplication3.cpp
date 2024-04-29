//#include <vector>
//#include <utility>
//#include <iostream>
//
//#define UNDEFINED -1
//
//class Point;
//typedef std::vector<Point> PointsMatrix;
//
//
//class Point
//{
//public:
//    std::vector<double> coordinates;
//    double reachable_dist;
//    PointsMatrix neighbours;
//
//
//    Point(const std::vector<double>& v) {
//        coordinates = v;
//        reachable_dist = UNDEFINED;
//    }
//};
//
//double euclideanDistance(std::vector<double> point1, std::vector<double> point2) {
//    double distance = 0.0;
//    for (int i = 0; i < point1.size(); ++i) {
//        distance += abs(point1[i] - point2[i]);
//    }
//
//    //return sqrt(distance);
//    return distance;
//}
//
//void getNeighbours(PointsMatrix& db, double eps) {
//    int size = db.size();
//
//    for (int i = 0; i < size; i++) {
//        Point& p1 = db[i];
//
//        for (int j = i + 1; j < size; j++) {
//            Point& p2 = db[j];
//
//            double distance = euclideanDistance(p1.coordinates, p2.coordinates);
//
//            if (distance <= eps) {
//                std::cout << "p" << i << ' ' << "p" << j << ' ' << distance << std::endl;
//
//                p1.neighbours.push_back(p2);
//                p2.neighbours.push_back(p1);
//            }
//        }
//    }
//}
//
//
//void optics(PointsMatrix db, double eps, int min_pts) {
//
//}
////#Initializing the reachability distance of the selected point
////pt.reachable_dist = UNDEFINED                                 T
////for each unprocessed point pt of DB                           T
////
////#Getting the neighbours of the selected point                 -
////#according to the definitions of epsilon and
////#minPts in DBSCAN
////Nbrs = getNbrs(pt, eps)
////
////mark pt as processed
////output pt to the ordered list
////
////#Checking if the selected point is not noise
////if (core_dist(pt, eps, Minpts) != UNDEFINED)
////
////#Initializing a priority queue to get the closest data point
////#in terms of Reachability distance
////Seeds = empty priority queue
////
////#Calling the update function
////update(Nbrs, pt, Seeds, eps, Minpts)
////
////#Repeating the process for the next closest point
////for each next q in Seeds
////Nbrs' = getNbrs(q, eps)
////mark q as processed
////output q to the ordered list
////if (core_dist(q, eps, Minpts) != UNDEFINED)
////update(Nbrs', q, Seeds, eps, Minpts)
////    The pseudo - code for the update function is given below :
////
////update(Nbrs, pt, Seeds, eps, MinPts)
////
////#Calculating the core distance for the given point
////coredist = core_dist(pt, eps, MinPts)
////
////#Updating the Reachability distance for each neighbour of p
////for each obj in Nbrs
////if (obj is not processed)
////new_reach_distance = max(coredist, dist(pt, obj))
////
////#Checking if the neighbour point is in seeds
////if (obj.reachable_dist == UNDEFINED)
////
////#Updation step
////obj.reachabled_dist = new_reach_distance
////Seeds.insert(obj, new_reach_distance)
////else
////if (new_reach_distance < obj.reachable_dist)
////
////    #Updation step
////    o.reachable_dist = new_reach_distance
////    Seeds.move - up(obj, new_reach_distance)
//int main() {
//    //std::cout << "hello we are in!!";
//
//
//    // Example usage
//    std::vector<std::vector<double>> v2d = {
//        {1, 1},
//        {2, 1},
//        {1, 2},
//        {2, 2},
//        {3, 5},
//        {3, 9},
//        {3, 10},
//        {4, 10},
//        {4, 11},
//        {5, 10},
//        {7, 10},
//        {10, 9},
//        {10, 6},
//        {9, 5},
//        {10, 5},
//        {11, 5},
//        {9, 4},
//        {10, 4},
//        {11, 4},
//        {10, 3}
//    };
//
//    double eps = 5.0;
//    int MinPts = 2;
//
//    PointsMatrix points;
//
//    for (const std::vector<double>& v : v2d) {
//        points.push_back(v);
//    }
//
//    getNeighbours(points, eps);
//
//    for (int i = 0; i < points.size(); i++) {
//        Point p = points[i];
//
//        std::cout << "p" << i << ": ";
//
//        for (int j = 0; j < p.neighbours.size(); j++)
//        {
//            if (j != 0) std::cout << ", ";
//            std::cout << "p" << j;
//        }
//
//        std::cout << std::endl;
//    }
//
//    //optics(DB, eps, MinPts);
//    return 0;
//}
//
//

//#include <iostream>
//#include <vector>
//#include <queue>
//#include <cmath>
//
//// Define a point structure
//struct Point {
//    double x, y; // assuming 2D points for simplicity
//    double reachability_distance;
//    bool processed;
//
//    Point(double _x, double _y) : x(_x), y(_y), reachability_distance(-1.0), processed(false) {}
//};
//
//// Function to calculate distance between two points
//double distance(Point p1, Point p2) {
//    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
//}
//
//// Function to get neighbors of a point within epsilon distance
//std::vector<Point*> getNeighbors(Point* p, std::vector<Point>& DB, double epsilon) {
//    std::vector<Point*> neighbors;
//    for (Point& q : DB) {
//        if (distance(*p, q) <= epsilon) {
//            neighbors.push_back(&q);
//        }
//    }
//    return neighbors;
//}
//
//// Function to calculate core distance of a point
//double core_distance(Point* p, std::vector<Point>& DB, double epsilon, int MinPts) {
//    std::vector<Point*> neighbors = getNeighbors(p, DB, epsilon);
//    if (neighbors.size() >= MinPts) {
//        // Sort distances
//        std::vector<double> distances;
//        for (Point* q : neighbors) {
//            distances.push_back(distance(*p, *q));
//        }
//        std::sort(distances.begin(), distances.end());
//        return distances[MinPts - 1];
//    }
//    return -1.0; // UNDEFINED
//}
//
//// Update function
//void update(std::vector<Point*>& N, Point* p, std::priority_queue<Point*, std::vector<Point*>,
//    bool(*)(const Point*, const Point*)>& Seeds, std::vector<Point>& DB, double epsilon, int MinPts) {
//    for (Point* q : N) {
//        if (!q->processed) {
//            double new_reachability_distance = std::max(core_distance(p, DB, epsilon, MinPts), distance(*p, *q));
//            if (q->reachability_distance == -1.0) {
//                q->reachability_distance = new_reachability_distance;
//                Seeds.push(q);
//            }
//            else {
//                if (new_reachability_distance < q->reachability_distance) {
//                    q->reachability_distance = new_reachability_distance;
//                    Seeds.push(q);
//                }
//            }
//        }
//    }
//}
//
//// OPTICS function
//void OPTICS(std::vector<Point>& DB, double epsilon, int MinPts) {
//    for (Point& p : DB) {
//        p.reachability_distance = -1.0; // UNDEFINED
//    }
//
//    std::priority_queue<Point*, std::vector<Point*>, bool(*)(const Point*, const Point*)> Seeds(
//        [](const Point* p1, const Point* p2) {
//            return p1->reachability_distance < p2->reachability_distance;
//        }
//    );
//
//    std::vector<Point*> ordered_list;
//    for (Point& p : DB) {
//        if (!p.processed) {
//            std::vector<Point*> N = getNeighbors(&p, DB, epsilon);
//            p.processed = true;
//            ordered_list.push_back(&p);
//            double core_dist = core_distance(&p, DB, epsilon, MinPts);
//            if (core_dist != -1.0) {
//                update(N, &p, Seeds, DB, epsilon, MinPts);
//                while (!Seeds.empty()) {
//                    Point* q = Seeds.top();
//                    Seeds.pop();
//                    std::vector<Point*> N_prime = getNeighbors(q, DB, epsilon);
//                    q->processed = true;
//                    ordered_list.push_back(q);
//                    if (core_distance(q, DB, epsilon, MinPts) != -1.0) {
//                        update(N_prime, q, Seeds, DB, epsilon, MinPts);
//                    }
//                }
//            }
//        }
//    }
//
//    // Output ordered list
//    std::cout << "Ordered List:\n";
//    for (Point* p : ordered_list) {
//        std::cout << "(" << p->x << ", " << p->y << ") - Reachability Distance: " << p->reachability_distance << "\n";
//    }
//}
//
//int main() {
//    // Example usage
//    std::vector<Point> DB = { 
//    
//                {1, 1},
//            {2, 1},
//            {1, 2},
//            {2, 2},
//            {3, 5},
//            {3, 9},
//            {3, 10},
//            {4, 10},
//            {4, 11},
//            {5, 10},
//            {7, 10},
//            {10, 9},
//            {10, 6},
//            {9, 5},
//            {10, 5},
//            {11, 5},
//            {9, 4},
//            {10, 4},
//            {11, 4},
//            {10, 3}
//    
//    }; // Example dataset
//    double epsilon = 5.0; // Example epsilon
//    int MinPts = 2; // Example MinPts
//    OPTICS(DB, epsilon, MinPts);
//    return 0;
//}

//#include <iostream>
//#include <vector>
//#include <queue>
//#include <cmath>
//
//struct Point {
//    double x, y;
//    double reachability_distance;
//    bool processed;
//
//    Point(double _x, double _y) : x(_x), y(_y), reachability_distance(std::numeric_limits<double>::infinity()), processed(false) {}
//};
//
//double distance(const Point& p1, const Point& p2) {
//    return std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2));
//}
//
//std::vector<Point*> getNeighbors(const Point& p, const std::vector<Point>& DB, double epsilon) {
//    std::vector<Point*> neighbors;
//    for (const Point& q : DB) {
//        if (distance(p, q) <= epsilon) {
//            neighbors.push_back(const_cast<Point*>(&q));
//        }
//    }
//    return neighbors;
//}
//
//double core_distance(const Point& p, const std::vector<Point>& DB, double epsilon, int MinPts) {
//    std::vector<double> distances;
//    for (const Point& q : DB) {
//        distances.push_back(distance(p, q));
//    }
//    std::sort(distances.begin(), distances.end());
//    return distances[MinPts - 1];
//}
//
//void update(std::vector<Point*>& N, Point* p, std::priority_queue<Point*, std::vector<Point*>,
//    bool(*)(const Point*, const Point*)>& Seeds, std::vector<Point>& DB, double epsilon, int MinPts) {
//    for (Point* q : N) {
//        if (!q->processed) {
//            double new_reachability_distance = std::max(core_distance(*p, DB, epsilon, MinPts), distance(*p, *q));
//            if (new_reachability_distance < q->reachability_distance) {
//                q->reachability_distance = new_reachability_distance;
//                Seeds.push(q);
//            }
//        }
//    }
//}
//
//void OPTICS(std::vector<Point>& DB, double epsilon, int MinPts) {
//    std::priority_queue<Point*, std::vector<Point*>, bool(*)(const Point*, const Point*)> Seeds(
//        [](const Point* p1, const Point* p2) {
//            return p1->reachability_distance < p2->reachability_distance;
//        }
//    );
//
//    std::vector<Point*> ordered_list;
//    for (Point& p : DB) {
//        if (!p.processed) {
//            std::vector<Point*> N = getNeighbors(p, DB, epsilon);
//            p.processed = true;
//            ordered_list.push_back(&p);
//            double core_dist = core_distance(p, DB, epsilon, MinPts);
//            if (core_dist != -1.0) {
//                update(N, &p, Seeds, DB, epsilon, MinPts);
//                while (!Seeds.empty()) {
//                    Point* q = Seeds.top();
//                    Seeds.pop();
//                    std::vector<Point*> N_prime = getNeighbors(*q, DB, epsilon);
//                    q->processed = true;
//                    ordered_list.push_back(q);
//                    if (core_distance(*q, DB, epsilon, MinPts) != -1.0) {
//                        update(N_prime, q, Seeds, DB, epsilon, MinPts);
//                    }
//                }
//            }
//        }
//    }
//
//    // Output ordered list
//    std::cout << "Ordered List:\n";
//    for (Point* p : ordered_list) {
//        std::cout << "(" << p->x << ", " << p->y << ") - Reachability Distance: " << p->reachability_distance << "\n";
//    }
//}
//
//int main() {
//    // Example usage
//    std::vector<Point> DB = { 
//    
//                {1, 1},
//            {2, 1},
//            {1, 2},
//            {2, 2},
//   /*         {3, 5},
//            {3, 9},
//            {3, 10},
//            {4, 10},
//            {4, 11},
//            {5, 10},
//            {7, 10},
//            {10, 9},
//            {10, 6},
//            {9, 5},
//            {10, 5},
//            {11, 5},
//            {9, 4},
//            {10, 4},
//            {11, 4},
//            {10, 3}*/
//    
//    }; // Example dataset
//    double epsilon = 5.0; // Example epsilon
//    int MinPts = 2; // Example MinPts
//    OPTICS(DB, epsilon, MinPts);
//    return 0;
//}

#include <iostream>
#include <vector>
#include <queue>
#include <cmath>

struct Point {
    std::vector<double> coordinates;
    double reachability_distance;
    bool processed;

    Point(std::vector<double> _coordinates) : coordinates(_coordinates), reachability_distance(std::numeric_limits<double>::infinity()), processed(false) {}
};

double distance(const Point& p1, const Point& p2) {
    double dist = 0;
    for (size_t i = 0; i < p1.coordinates.size(); ++i) {
        dist += std::abs(p2.coordinates[i] - p1.coordinates[i]);
    }
    return dist;
}

std::vector<Point*> getNeighbors(const Point& p, const std::vector<Point>& DB, double epsilon) {
    std::vector<Point*> neighbors;
    for (const Point& q : DB) {
        if (distance(p, q) <= epsilon) {
            neighbors.push_back(const_cast<Point*>(&q));
        }
    }
    return neighbors;
}

double core_distance(const Point& p, const std::vector<Point>& DB, double epsilon, int MinPts) {
    std::vector<double> distances;
    
    for (const Point& q : DB) {
        distances.push_back(distance(p, q));
    }
    std::sort(distances.begin(), distances.end()); // []

    return distances[MinPts - 1];
}

void update(std::vector<Point*>& N, Point* p, std::priority_queue<Point*, std::vector<Point*>,
    bool(*)(const Point*, const Point*)>& Seeds, std::vector<Point>& DB, double epsilon, int MinPts) {

    for (Point* q : N) {
        if (!q->processed) {
            double new_reachability_distance = std::max(core_distance(*p, DB, epsilon, MinPts), distance(*p, *q));
            if (new_reachability_distance < q->reachability_distance) {
                q->reachability_distance = new_reachability_distance;
                Seeds.push(q);
            }
        }
    }
}

void OPTICS(std::vector<Point>& DB, double epsilon, int MinPts) {
    std::priority_queue<Point*, std::vector<Point*>, bool(*)(const Point*, const Point*)> Seeds(
        [](const Point* p1, const Point* p2) {
            return p1->reachability_distance < p2->reachability_distance;
        }
    );

    std::vector<Point*> ordered_list;

    for (Point& p : DB) {
        if (!p.processed) {
            std::vector<Point*> N = getNeighbors(p, DB, epsilon);
            p.processed = true;
            ordered_list.push_back(&p);
            double core_dist = core_distance(p, DB, epsilon, MinPts);
           
            if (core_dist != -1.0) {
                update(N, &p, Seeds, DB, epsilon, MinPts);
                while (!Seeds.empty()) { 
                    Point* q = Seeds.top(); // B
                    Seeds.pop();
                    std::vector<Point*> N_prime = getNeighbors(*q, DB, epsilon);
                    q->processed = true;
                    ordered_list.push_back(q); // A D
                    core_dist = core_distance(*q, DB, epsilon, MinPts);
                    if (core_dist != -1.0) {
                        update(N_prime, q, Seeds, DB, epsilon, MinPts);
                    }
                }
            }
        }
    }

    // Output ordered list
    std::cout << "Ordered List:\n";
    for (Point* p : ordered_list) {
        std::cout << "Coordinates: ";
        for (double coord : p->coordinates) {
            std::cout << coord << " ";
        }
        std::cout << "- Reachability Distance: " << p->reachability_distance << "\n";
    }
}

int main() {
    // Example usage
    std::vector<Point> DB;
    DB.push_back({ { 1, 1 } });
    DB.push_back({ { 2, 1 } });
    DB.push_back({ { 1, 2 } });
    DB.push_back({ { 2, 2 } });
    DB.push_back({ { 3, 5 } });
    DB.push_back({ { 3, 9 } });
    DB.push_back({ { 3, 10 } });
    DB.push_back({ { 4, 10 } });
    DB.push_back({ { 4, 11 } });
    DB.push_back({ { 5, 10 } });
    DB.push_back({ { 7, 10 } });
    DB.push_back({ { 10, 9 } });
    DB.push_back({ { 10, 6 } });
    DB.push_back({ { 9, 5 } });
    DB.push_back({ { 10, 5 } });
    DB.push_back({ { 11, 5 } });
    DB.push_back({ { 9, 4 } });
    DB.push_back({ { 10, 4 } });
    DB.push_back({ { 11, 4 } });
    DB.push_back({{ 10, 3 }});
    DB.push_back({ { 1, 1 } });
    DB.push_back({ { 2, 1 } });

    double epsilon = 5.0; // Example epsilon
    int MinPts = 2; // Example MinPts
    OPTICS(DB, epsilon, MinPts);
    return 0;
}
