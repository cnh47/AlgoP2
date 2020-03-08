// A C++ program to find convex hull of a set of points. Refer
// https://www.geeksforgeeks.org/orientation-3-ordered-points/
// for explanation of orientation()

//This implementation is taken from
//https://www.geeksforgeeks.org/convex-hull-set-1-jarviss-algorithm-or-wrapping/
//with permission from Prof. Duan to use in this project.

#ifndef JARVISMARCH_HPP_INCLUDED
#define JARVISMARCH_HPP_INCLUDED

#include "convexHull435.hpp"

class Jarvis
{
    public:
        //find orientation of ordered triplet
        int orientation(Point p, Point q, Point r);

        //prints convex hull of a set of n points
        void convexHull(std::vector<Point> points, int n, std::ofstream &output);

        //point structure


};

#endif
