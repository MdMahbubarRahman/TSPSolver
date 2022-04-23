#include <iostream>
#include <iterator>
#include <cmath>
#include <vector>
#include <list>
#include <map>
#include <fstream>
#include <string>
#include <sstream>

#include "BranchAndBound.h"
#include "HungarianAlg.h"
#include "OperatorTheory.h"

int main() {

    std::fstream myFile;
    //open raod distance data file
    myFile.open("DistanceMatrix.txt", std::ios::in);//read
    std::vector<std::vector<double>> roadDistance;
    if (myFile.is_open()) {
        std::string line, word;
        std::istringstream iss;
        int rowNum = 0;
        while (std::getline(myFile, line)) {
            if (rowNum > 0) {
                std::vector<double> dist;
                iss.clear();
                iss.str(line);
                int colNum = 0;
                while (iss.good()) {
                    iss >> word;
                    //&& colNum <= NUMBER_OF_NODES
                    if (colNum > 0) {
                        double abc = std::stod(word);//abc == road distance
                        //std::cout << abc << std::endl;
                        dist.push_back(abc);
                    }
                    colNum += 1;
                }
                roadDistance.push_back(dist);
            }
            rowNum += 1;
        }
        myFile.close();
    }

    std::vector<int> initialtour;
    initialtour.push_back(4);
    initialtour.push_back(10);
    initialtour.push_back(13);
    initialtour.push_back(40);
    initialtour.push_back(50);
    initialtour.push_back(43);
    initialtour.push_back(14);
    initialtour.push_back(16);
    initialtour.push_back(17);

    BranchAndBoundSolver bbSolver(initialtour, roadDistance);
    bbSolver.runBranchAndBoundSolver();


    return 0;
}

