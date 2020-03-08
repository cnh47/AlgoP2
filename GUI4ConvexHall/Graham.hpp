#ifndef GRAHAM_HPP_INCLUDED
#define GRAHAM_HPP_INCLUDED
// A C++ program to find convex hull of a set of points. Refer
// https://www.geeksforgeeks.org/orientation-3-ordered-points/
// for explanation of orientation()

#include "convexHull435.hpp"

class Graham{
private:
    // A global point needed for sorting points with reference
    // to the first point Used in compare function of qsort()
    Point p0;
public:

    Point GetPoint(){
        return p0;
    }
    // A utility function to find next to top in a stack
    Point nextToTop(std::stack<Point> &S);

    // A utility function to swap two points
    int swap(Point &p1, Point &p2);

    // A utility function to return square of distance
    // between p1 and p2
    int distSq(Point p1, Point p2);

    // To find orientation of ordered triplet (p, q, r).
    // The function returns following values
    // 0 --> p, q and r are colinear
    // 1 --> Clockwise
    // 2 --> Counterclockwise
    int orientation(Point p, Point q, Point r);

    // Prints convex hull of a set of n points.
    void convexHull(std::vector<Point> points, int n, std::ofstream &output);
};
#endif
