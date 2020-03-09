// You need to complete this program for your second project.
// Standard libraries
#include <stack>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>

// #include "JarvisMarch.hpp"
//#include "Graham.hpp"
// #include "convexHull435.hpp"

struct Point
{
	int x, y;
};

int JarvisOrientation(Point p, Point q, Point r);
void JarvisConvexHull(std::vector<Point> points, int n, std::ofstream &output);

Point GrahamNextToTop(std::stack<Point> &S);
void GrahamSwap(Point &p1, Point &p2);
int GrahamDistSq(Point p1, Point p2);
int GrahamOrientation(Point p, Point q, Point r);
int GrahamCompare(const void *vp1, const void *vp2);
void GrahamConvexHull(Point points[], int n, std::ofstream &output);

// iPair is integer pairs
#define iPair std::pair<int, int>

// Stores the result (points of convex hull)
std::set<iPair> hull;

int QuickFindSide(iPair p1, iPair p2, iPair p);
int QuickLineDist(iPair p1, iPair p2, iPair p);
void QuickHull(iPair a[], int n, iPair p1, iPair p2, int side);
void QuickPrintHull(iPair a[], int n, std::ofstream &output);

int main(int argc, char *argv[])
{
    if (argc < 3)
        std::cout << "wrong format! should be \"a.exe algType dataFile\"";
    else
    {
        std::string algType = argv[1];
        std::string dataFilename = argv[2];
        std::string outputFile = "";
        std::ifstream dataFile;
        dataFile.open(dataFilename);
        //read your data points from dataFile (see class example for the format)

        switch(algType[0])
        {
            case 'G':
            {
                //call your Graham Scan algorithm to solve the problem
                Point tmp;
                std::vector<Point> vPoints;
                int n=0, x=0, y=0;

                while(dataFile >> tmp.x >> tmp.y)
                {
                    vPoints.push_back(tmp);
                    ++n;
                }

                int vSize = vPoints.size();
                Point *aPoints = new Point[vSize];
                for(int i = 0; i < vSize; ++i){
                    aPoints[vSize - i] = vPoints[vSize - i];
                }

                outputFile = "hull_G_" + dataFilename;
                std::ofstream hullFile(outputFile);
                GrahamConvexHull(aPoints, n, hullFile);
                hullFile.close();
				delete [] aPoints;
				aPoints = NULL;
                break;
            }
            case 'J':
                {
                    //call your Javis March algorithm to solve the problem
                    Point tmp;
                    std::vector<Point> points;
                    int n=0, x=0, y=0;

                    while(dataFile >> tmp.x >> tmp.y)
                    {
                        points.push_back(tmp);
                        ++n;
                    }

                    outputFile = "hull_J_" + dataFilename;
                    std::ofstream hullFile(outputFile);
                    JarvisConvexHull(points, n, hullFile);
                    hullFile.close();
                    break;
                }
            case 'Q':
            {
                //call your Quickhull algorithm to solve the problem

                iPair tmp;
                std::vector<iPair> vPair;
                int n=0, x=0, y=0;

                while(dataFile >> tmp.first >> tmp.second)
                {
                    vPair.push_back(tmp);
                    ++n;
                }

                int vSize = vPair.size();
                iPair *aPoints = new iPair [vSize];
                for(int i = 0; i < vSize; ++i){
                    aPoints[vSize - i] = vPair[vSize - i];
                }

                outputFile = "hull_Q_" + dataFilename;
                std::ofstream hullFile(outputFile);
                QuickPrintHull(aPoints, n, hullFile);
                hullFile.close();
                break;

                outputFile = "hull_Q.txt";
                break;
            }
            default:
                //any other parameter is called
                std::cout << "You done messed up AA Ron!" << std::endl;
                break;
        }

        //write your convex hull to the outputFile (see class example for the format)
        //you should be able to visulize your convex hull using the "ConvexHull_GUI" program.
        dataFile.close();
    }

    return 0;
}

////////////////////////////////////////////////////////////////
//                      JARVIS MARCH                          //
////////////////////////////////////////////////////////////////
// A C++ program to find convex hull of a set of points. Refer
// https://www.geeksforgeeks.org/orientation-3-ordered-points/
// for explanation of orientation()

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int JarvisOrientation(Point p, Point q, Point r)
{
	int val = (q.y - p.y) * (r.x - q.x) -
			(q.x - p.x) * (r.y - q.y);

	if (val == 0) return 0; // colinear
	return (val > 0)? 1: 2; // clock or counterclock wise
}

// Prints convex hull of a set of n points.
void JarvisConvexHull(std::vector<Point> points, int n, std::ofstream &output)
{
	// There must be at least 3 points
	if (n < 3)
        return;

	// Initialize Result
	std::vector<Point> hull;

	// Find the leftmost point
	int l = 0;
	for (int i = 1; i < n; i++)
		if (points[i].x < points[l].x)
			l = i;

	// Start from leftmost point, keep moving counterclockwise
	// until reach the start point again. This loop runs O(h)
	// times where h is number of points in result or output.
	int p = l, q;
	do
	{
		// Add current point to result
		hull.push_back(points[p]);

		// Search for a point 'q' such that orientation(p, x,
		// q) is counterclockwise for all points 'x'. The idea
		// is to keep track of last visited most counterclock-
		// wise point in q. If any point 'i' is more counterclock-
		// wise than q, then update q.
		q = (p+1)%n;
		for (int i = 0; i < n; i++)
		{
		// If i is more counterclockwise than current q, then
		// update q
		if (JarvisOrientation(points[p], points[i], points[q]) == 2)
			q = i;
		}

		// Now q is the most counterclockwise with respect to p
		// Set p as q for next iteration, so that q is added to
		// result 'hull'
		p = q;

	} while (p != l); // While we don't come to first point

	// Print Result
	for (int i = 0; i < hull.size(); i++)
		output << hull[i].x << '\t' << hull[i].y << '\n';
}

////////////////////////////////////////////////////////////////
//                      End of JARVIS MARCH                   //
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
//                      Graham Scan                           //
////////////////////////////////////////////////////////////////
// A C++ program to find convex hull of a set of points. Refer
// https://www.geeksforgeeks.org/orientation-3-ordered-points/
// for explanation of orientation()

// A global point needed for sorting points with reference
// to the first point Used in compare function of qsort()
Point p0;

// A utility function to find next to top in a stack
Point GrahamNextToTop(std::stack<Point> &S)
{
	Point p = S.top();
	S.pop();
	Point res = S.top();
	S.push(p);
	return res;
}

// A utility function to swap two points
void GrahamSwap(Point &p1, Point &p2)
{
	Point temp = p1;
	p1 = p2;
	p2 = temp;
}

// A utility function to return square of distance
// between p1 and p2
int GrahamDistSq(Point p1, Point p2)
{
	return (p1.x - p2.x)*(p1.x - p2.x) +
		(p1.y - p2.y)*(p1.y - p2.y);
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int GrahamOrientation(Point p, Point q, Point r)
{
	int val = (q.y - p.y) * (r.x - q.x) -
			(q.x - p.x) * (r.y - q.y);

	if (val == 0) return 0; // colinear
	return (val > 0)? 1: 2; // clock or counterclock wise
}

// A function used by library function qsort() to sort an array of
// points with respect to the first point
int GrahamCompare(const void *vp1, const void *vp2)
{
    Point *p1 = (Point *)vp1;
    Point *p2 = (Point *)vp2;

    // Find orientation
    int o = GrahamOrientation(p0, *p1, *p2);
    if (o == 0)
    	return (GrahamDistSq(p0, *p2) >= GrahamDistSq(p0, *p1))? -1 : 1;

    return (o == 2)? -1: 1;
}

// Prints convex hull of a set of n points.
void GrahamConvexHull(Point points[], int n, std::ofstream &output)
{
    // Find the bottommost point
    int ymin = points[0].y, min = 0;
    for (int i = 1; i < n; i++)
    {
    	int y = points[i].y;

    	// Pick the bottom-most or chose the left
    	// most point in case of tie
    	if ((y < ymin) || (ymin == y &&
    		points[i].x < points[min].x))
    		ymin = points[i].y, min = i;
    }

    // Place the bottom-most point at first position
    GrahamSwap(points[0], points[min]);

    // Sort n-1 points with respect to the first point.
    // A point p1 comes before p2 in sorted output if p2
    // has larger polar angle (in counterclockwise
    // direction) than p1
    p0 = points[0];
    qsort(&points[1], n-1, sizeof(Point), GrahamCompare);

    // If two or more points make same angle with p0,
    // Remove all but the one that is farthest from p0
    // Remember that, in above sorting, our criteria was
    // to keep the farthest point at the end when more than
    // one points have same angle.
    int m = 1; // Initialize size of modified array
    for (int i=1; i<n; i++)
    {
    	// Keep removing i while angle of i and i+1 is same
    	// with respect to p0
    	while (i < n-1 && GrahamOrientation(p0, points[i],
    									points[i+1]) == 0)
    		i++;


    	points[m] = points[i];
    	m++; // Update size of modified array
    }

    // If modified array of points has less than 3 points,
    // convex hull is not possible
    if (m < 3) return;

    // Create an empty stack and push first three points
    // to it.
    std::stack<Point> S; ;
    S.push(points[0]);
    S.push(points[1]);
    S.push(points[2]);

    // Process remaining n-3 points
    for (int i = 3; i < m; i++)
    {
    	// Keep removing top while the angle formed by
    	// points next-to-top, top, and points[i] makes
    	// a non-left turn
        while (GrahamOrientation(GrahamNextToTop(S), S.top(), points[i]) != 2)
         S.pop();
      S.push(points[i]);
    }

    // Now stack has the output points, print contents of stack
    while (!S.empty())
    {
    	Point p = S.top();
    	output << p.x << '\t' << p.y << std::endl;
    	S.pop();
    }
}

////////////////////////////////////////////////////////////////
//                    End of Graham Scan                      //
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
//                       Quickhull                            //
////////////////////////////////////////////////////////////////
// C++ program to implement Quick Hull algorithm
// to find convex hull.
//using namespace std;

// Returns the side of point p with respect to line
// joining points p1 and p2.
int QuickFindSide(iPair p1, iPair p2, iPair p)
{
	int val = (p.second - p1.second) * (p2.first - p1.first) -
			(p2.second - p1.second) * (p.first - p1.first);

	if (val > 0)
		return 1;
	if (val < 0)
		return -1;
	return 0;
}

// returns a value proportional to the distance
// between the point p and the line joining the
// points p1 and p2
int QuickLineDist(iPair p1, iPair p2, iPair p)
{
	return abs ((p.second - p1.second) * (p2.first - p1.first) -
			(p2.second - p1.second) * (p.first - p1.first));
}

// End points of line L are p1 and p2. side can have value
// 1 or -1 specifying each of the parts made by the line L
void QuickHull(iPair a[], int n, iPair p1, iPair p2, int side)
{
	int ind = -1;
	int max_dist = 0;

	// finding the point with maximum distance
	// from L and also on the specified side of L.
	for (int i=0; i<n; i++)
	{
		int temp = QuickLineDist(p1, p2, a[i]);
		if (QuickFindSide(p1, p2, a[i]) == side && temp > max_dist)
		{
			ind = i;
			max_dist = temp;
		}
	}

	// If no point is found, add the end points
	// of L to the convex hull.
	if (ind == -1)
	{
		hull.insert(p1);
		hull.insert(p2);
		return;
	}

	// Recur for the two parts divided by a[ind]
	QuickHull(a, n, a[ind], p1, -QuickFindSide(a[ind], p1, p2));
	QuickHull(a, n, a[ind], p2, -QuickFindSide(a[ind], p2, p1));
}

void QuickPrintHull(iPair a[], int n, std::ofstream &output)
{
	// a[i].second -> y-coordinate of the ith point
	if (n < 3)
	{
		std::cout << "Convex hull not possible\n";
		return;
	}

	// Finding the point with minimum and
	// maximum x-coordinate
	int min_x = 0, max_x = 0;
	for (int i=1; i<n; i++)
	{
		if (a[i].first < a[min_x].first)
			min_x = i;
		if (a[i].first > a[max_x].first)
			max_x = i;
	}

	// Recursively find convex hull points on
	// one side of line joining a[min_x] and
	// a[max_x]
	QuickHull(a, n, a[min_x], a[max_x], 1);

	// Recursively find convex hull points on
	// other side of line joining a[min_x] and
	// a[max_x]
	QuickHull(a, n, a[min_x], a[max_x], -1);
	
	while (!hull.empty())
	{
		output << (*hull.begin()).first << "\t" << (*hull.begin()).second << std::endl;
		hull.erase(hull.begin());
	}
}

////////////////////////////////////////////////////////////////
//                    End of Quickhull                        //
////////////////////////////////////////////////////////////////
