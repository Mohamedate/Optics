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
