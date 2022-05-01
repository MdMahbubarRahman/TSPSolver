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
    initialtour.push_back(0);
    initialtour.push_back(10);
    initialtour.push_back(9);
    initialtour.push_back(15);
    initialtour.push_back(8);
    initialtour.push_back(7);
    initialtour.push_back(6);

    std::cout << "\nSize of the initial tour : " << initialtour.size() << std::endl;

    BranchAndBoundSolver bbSolver(initialtour, roadDistance);
    bbSolver.runBranchAndBoundSolver();
    TSPSolution tspSol = bbSolver.getTSPSolution();
    std::cout << "\nInitial tsp tour : " << std::endl;
    for (auto it : initialtour) {
        std::cout << it << " ";
    }
    std::cout << " " << std::endl;
    double cost = 0;
    for (int i = 0; i < initialtour.size() - 1; i++) {
        cost += roadDistance[initialtour.at(i)][initialtour.at(i + 1)];
    }
    cost += roadDistance[initialtour.at(initialtour.size() - 1)][initialtour.at(0)];
    std::cout << "\nInitial tsp tour cost : " << cost << std::endl;
    std::cout << "\nFinal tsp tour : " << std::endl;
    for (auto& it : tspSol.tour) {
        std::cout << it << " ";
    }
    std::cout << " " << std::endl;
    std::cout << "\nFinal tsp solution cost : " << tspSol.cost << std::endl;

    return 0;
}

