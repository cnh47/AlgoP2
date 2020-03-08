// A C++ program to find convex hull of a set of points. Refer
// https://www.geeksforgeeks.org/orientation-3-ordered-points/
// for explanation of orientation()

#include <stack>
#include <iostream>
#include <fstream>

struct Point
{
    int x, y;
};

class Graham{
public:

    // A global point needed for sorting points with reference
    // to the first point Used in compare function of qsort()
    Point p0;

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

    // A function used by library function qsort() to sort an array of
    // points with respect to the first point
    int compare(const void *vp1, const void *vp2);

    // Prints convex hull of a set of n points.
    void convexHull(Point points[], int n);
};
