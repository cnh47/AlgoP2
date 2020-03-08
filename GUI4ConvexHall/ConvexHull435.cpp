// You need to complete this program for your second project.

// Standard libraries
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "JarvisMarch.hpp"
#include "Graham.hpp"
#include "QuickHull.hpp"

void ReadInPoint(std::string);
int GetArrySize(std::string FileName);

int main(int argc, char *argv[])
{
    Jarvis Jarvis;
    if (argc < 3)
      std::cout << "wrong format! should be \"a.exe algType dataFile\"";
    else {
      std::string algType = argv[1];
      std::string dataFilename = argv[2];

      std::string outputFile = "";
      //read your data points from dataFile (see class example for the format)

        if (algType[0]=='G') {
           //call your Graham Scan algorithm to solve the problem

           outputFile = "hull_G.txt";
        }
        else if (algType[0]=='J') {
           //call your Javis March algorithm to solve the problem
           outputFile = "hull_J.txt";
        }
        else { //default
           //call your Quickhull algorithm to solve the problem
           outputFile = "hull_Q.txt";
        }

      //write your convex hull to the outputFile (see class example for the format)
      //you should be able to visulize your convex hull using the "ConvexHull_GUI" program.
    }

    return 0;
}

int GetArrySize(std::string FileName){
    ifstream InFile;
    int size;
    std::string x;
    InFile.open(FileName);
    while(InFile >> x){
        size++;
    }
    return size;
}

void ReadInPoints(std::string FileName){
    ifstream InFile;
    int x, y;
    InFile.open(FileName);

    while(!InFile >> x >> y){

    }
}
