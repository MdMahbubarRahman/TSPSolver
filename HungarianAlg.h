#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <Queue>
#include <map>

#ifndef HUNGARIANALG_H
#define HUNGARIANALG_H

//The Hungarian algorithm solves assignment problem which is a relaxation of the traveling salesman problem.
//Hungarian algorithm generates either a TSP tour or a number of subtours.

class HungarianAlg {
private:
	double reducedCost;
	double baseValue;
	double initSolCost;
	double finalSolCost;
	bool tspOptimal;
	int reducedCostMatrixSize;
	struct IntersectionPoint {
		int row;
		int column;
		IntersectionPoint(int a, int b) {
			row = a;
			column = b;
		}
	};
	std::map<int, int> boxPoints;
	bool assignmentOptimal;
	std::list<std::vector<int>> listOfRoutes;
	std::vector<int> initialTsp;//tsp presentation node/node/.../node/depotnode
	std::vector<int> incumbentTsp;
	std::vector<int> optimalTsp;
	std::vector<std::vector<double>> costMatrix;
	std::vector<std::vector<double>> reducedCostMatrix;
public:
	HungarianAlg();
	HungarianAlg(const HungarianAlg & hunAlg);
	HungarianAlg(std::vector<int> initialTsp, std::vector<std::vector<double>> costMatrix);
	void performRowColumnReduction();
	void rowScanning(std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints);
	void columnScanning(std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints);
	void diagonalSelection(std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints);
	void diagonalSelection1(std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints);
	void diagonalSelection2(std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints);
	void diagonalSelectionAdv(std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints);
	int uncoveredZeroScanning(std::set<int>& coveredRows, std::set<int>& coveredColumns);
	void coverMinimumValueAssignmentsAndCheckOptimality(bool& assignmentOptimal, std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints, std::list<IntersectionPoint>& twiceCoveredPoints);
	void performRowColumnReductionForUncoveredCells(bool& assignmentOptimal, std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints, std::list<IntersectionPoint>& twiceCoveredPoints);
	//void generateUpperBoundTsp(std::list<std::vector<int>> listOfRoutes);
	void populateListOfRoutes(std::map<int, int> boxPoints);
	void checkListOfRoutes(std::map<int, int> boxPoints, std::list<std::vector<int>>& routeList);
	void showAssignmentSolution();
	void runHungarianAlg();
};

#endif