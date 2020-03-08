// A C++ program to find convex hull of a set of points. Refer
// https://www.geeksforgeeks.org/orientation-3-ordered-points/
// for explanation of orientation()
//#include <bits/stdc++.h>
using namespace std;
class Jarvis{
public:
    Jarvis();

    struct Point
    {
        int x, y;
    };

    // To find orientation of ordered triplet (p, q, r).
    // The function returns following values
    // 0 --> p, q and r are colinear
    // 1 --> Clockwise
    // 2 --> Counterclockwise
    int orientation(Point p, Point q, Point r);

    // Prints convex hull of a set of n points.
    void convexHull(Point points[], int n);
};
