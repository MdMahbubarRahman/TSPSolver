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
    /*
	TODOs:
	1. Framework for variable number of nodes
	2. Framework for importing distance matrix
	3. Modular Hungarian Algorithm to solve assignment problem.
	4. Framework for generating basic solution from assignment solution.
	5. Framework for cost matrix.
	6. Framework for operator theory for parametric programming
	7. Framework for scanning routine
	8. Framework for branch and bound tree	
	*/

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
    /*
    std::cout << "Show road distance matrix" << std::endl;
    std::cout << "size of the roadDistance matrix : " << roadDistance.size() << std::endl;
    for (int i = 0; i < roadDistance.size(); i++) {
        for (int j = 0; j < roadDistance.size(); j++) {
            std::cout << roadDistance[i][j] << " ";
        }
        std::cout << ";" << std::endl;
    }
    */

    //open aerial distance data file
    myFile.open("GCDistanceMatrix.txt", std::ios::in);//read
    std::vector<std::vector<double>> aerialDistance;
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
                    if (colNum > 0) {
                        double abc = std::stod(word);//abc == road distance
                        //std::cout << abc << std::endl;
                        dist.push_back(abc);
                    }
                    colNum += 1;
                }
                aerialDistance.push_back(dist);
            }
            rowNum += 1;
        }
        myFile.close();
    }
    /*
    std::cout << "Show aerial distance matrix" << std::endl;
    std::cout << "size of the aerialDistance matrix : " << aerialDistance.size() << std::endl;
    for (int i = 0; i < aerialDistance.size(); i++) {
        for (int j = 0; j < aerialDistance.size(); j++) {
            std::cout << aerialDistance[i][j] << " ";
        }
        std::cout << ";" << std::endl;
    }
    */


    std::vector<int> tour;
 
    tour.push_back(15);
    tour.push_back(9);
    tour.push_back(8);
    tour.push_back(6);
    tour.push_back(0);
    tour.push_back(10);
    tour.push_back(56);
    tour.push_back(45);
    tour.push_back(40);
    tour.push_back(52);
    tour.push_back(30);
    tour.push_back(31);
    tour.push_back(32);
    tour.push_back(33);
    tour.push_back(34);
    /*
    tour.push_back(35);
    tour.push_back(36);
    tour.push_back(37);
    tour.push_back(38);
    tour.push_back(39);
    
    HungarianAlg abc(tour, roadDistance);
    abc.runHungarianAlg();
    abc.showAssignmentSolution();
    */

    OperatorTheory oper;
    oper.generateCycleAndUpdateBasicSolution();

	return 0;
}