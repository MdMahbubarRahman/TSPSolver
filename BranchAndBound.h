#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <Queue>
#include <map>

#include "HungarianAlg.h"
#include "OperatorTheory.h"

#ifndef BRANCHANDBOUND_H
#define BRANCHANDBOUND_H

enum BBSolverStatus {
	Optimal, SubOptimal, LocalOptimal, TerminatedWithErrors
};

enum CostOperatorStatus {
	CostOperatorTerminatedSuccessfully, CostOperatorTerminatedWithErrors
};

enum HungarianAlgStatus {
	HungarianAlgTerminatedSuccessfully, HungarianAlgTerminatedWithErrors
};

struct RootNode {
	int nodeIndex;
	bool isFinalSolutionATour;
	double initialWeakerLowerBound;
	double finalWeakerLowerBound;
	HungarianAlgStatus hungarianAlgStatus;
	std::vector<int> listOfCities;
	std::vector<std::vector<double>> initialCostTableau;
	std::vector<std::vector<double>> finalCostCostTableau;
	std::list<BasicCell> initialSolution;
	std::list<BasicCell> finalSolution;
};

struct Node {
	int nodeIndex;
	int parentNodeIndex;
	bool isFinalSolutionATour;
	double initialWeakerLowerBound;
	double finalWeakerLowerBound;
	Cell branchOnCell;
	CostOperatorStatus costOperatorStatus;
	std::vector<int> listOfCities;
	std::vector<std::vector<double>> initialCostTableau;
	std::vector<std::vector<double>> finalCostCostTableau;
	std::list<BasicCell> initialSolution;
	std::list<BasicCell> finalSolution;
	std::vector<double> initialRowDualSolution;
	std::vector<double> initialColumnDualSolution;
	std::vector<double> finalRowDualSolution;
	std::vector<double> finalColumnDualSolution;
	std::list<Node> childNodes;
};

struct TourSolution {
	std::list<BasicCell> tourSolution;
	double objValue;
};

struct Incumbent {
	std::list<BasicCell> incumbentSolution;
	std::vector<int> tour;
	double objValue;
};

struct TSPSolution {
	double cost;
	std::vector<int> tour;
	BBSolverStatus status;
};

struct BBTree {
	double totalTime;
	double maxTimeLimit;
	BBSolverStatus status;
	Incumbent incumbent;
	int numOfTourSolution;
	std::list<TourSolution> tourSolutions;
	int numOfNodePrunedByBound;
	int numOfNodePrunedByIntegrality;
	int numOfNodePrunedWithErrorsOrInfeasibility;
	double lowerBound;
	double weakerLowerBound;
	double upperBound;
	std::list<Node> branchNodes;
	int currentNumOfNodes;
};

//Branch and Bound Solver for solvinig traveling salesman problem
class BranchAndBoundSolver {
private:
	std::vector<std::vector<double>> costTableau;
	std::vector<int> oldListOfCities;
	std::vector<int> newListOfCities;
	std::map<int, int> newCityToOldCityMap;
	std::vector<BasicCell> assignmentSolution;
	std::vector<BasicCell> transportationSolution;
	std::list<std::vector<int>> routeLists;
	TSPSolution tspSolution;
	BBTree bbTree;
public:
	BranchAndBoundSolver();
	BranchAndBoundSolver(std::vector<int> initialTourOfCities, std::vector<std::vector<double>> wholeCostMatrix);
	BranchAndBoundSolver(const BranchAndBoundSolver& depthFirstBBSolver);
	void solveAssignmentProblem();
	void generateBasicTransportationSolution();
	void solveTransportationProblem();
	TSPSolution getTSPSolution();
	void runBranchAndBoundSolver();
};

#endif




