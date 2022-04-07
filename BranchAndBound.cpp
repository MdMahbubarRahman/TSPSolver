#include "BranchAndBound.h"


//default constructor
BranchAndBoundSolver::BranchAndBoundSolver() {
	bbTree.totalTime = 0;
	bbTree.maxTimeLimit = 1000;//1000sec
	bbTree.numOfTourSolution = 0;
	bbTree.numOfNodePrunedByBound = 0;
	bbTree.numOfNodePrunedByIntegrality = 0;
	bbTree.numOfNodePrunedWithErrorsOrInfeasibility = 0;
	bbTree.lowerBound = -INFINITY;
	bbTree.weakerLowerBound = -INFINITY;
	bbTree.upperBound = INFINITY;
	bbTree.currentNumOfNodes = 0;
}

//constructor
BranchAndBoundSolver::BranchAndBoundSolver(std::vector<int> initialTourOfCities, std::vector<std::vector<double>> wholeCostMatrix) {
	std::vector<double> costVec;
	oldListOfCities = initialTourOfCities;
	std::cout << "\nPopulate new city to old city map" << std::endl;
	for (int i = 0; i < oldListOfCities.size(); i++) {
		newCityToOldCityMap.insert(std::pair<int, int>(i, oldListOfCities.at(i)));
		newListOfCities.push_back(i);
	}
	std::cout << "\nPopulate the cost matrix for the current cities." << std::endl;
	double val = 0;
	for (int i = 0; i < oldListOfCities.size(); i++) {
		for (int j = 0; j < oldListOfCities.size(); j++) {
			val = wholeCostMatrix[oldListOfCities.at(i)][oldListOfCities.at(j)];
			costVec.push_back(val);
		}
		costTableau.push_back(costVec);
		costVec.clear();
	}
	std::cout << "\nPopulate bbTree with initial values." << std::endl;
	bbTree.totalTime = 0;
	bbTree.maxTimeLimit = 1000;//1000sec
	bbTree.numOfTourSolution = 0;
	bbTree.numOfNodePrunedByBound = 0;
	bbTree.numOfNodePrunedByIntegrality = 0;
	bbTree.numOfNodePrunedWithErrorsOrInfeasibility = 0;
	bbTree.lowerBound = -INFINITY;
	bbTree.weakerLowerBound = -INFINITY;
	bbTree.upperBound = INFINITY;
	bbTree.currentNumOfNodes = 0;
}

//copy constructor
BranchAndBoundSolver::BranchAndBoundSolver(const BranchAndBoundSolver& depthFirstBBSolver) {
	costTableau = depthFirstBBSolver.costTableau;
	oldListOfCities = depthFirstBBSolver.oldListOfCities;
	newListOfCities = depthFirstBBSolver.newListOfCities;
	newCityToOldCityMap = depthFirstBBSolver.newCityToOldCityMap;
	tspSolution = depthFirstBBSolver.tspSolution;
	bbTree = depthFirstBBSolver.bbTree;
	assignmentSolution = depthFirstBBSolver.assignmentSolution;
	transportationSolution = depthFirstBBSolver.transportationSolution;
	routeLists = depthFirstBBSolver.routeLists;
}

//solves assignment problem using Hungarian algorithm
void BranchAndBoundSolver::solveAssignmentProblem() {
	HungarianAlg hungAlg(newListOfCities, costTableau);
	hungAlg.runHungarianAlg();
	std::map<int, int> assignmentSol;
	routeLists = hungAlg.getListOfRoutes();
	assignmentSol = hungAlg.getAssignmentSolution();
	for (auto& it : assignmentSol) {
		BasicCell basicCell = BasicCell();
		basicCell.rowID = it.first;
		basicCell.colID = it.second;
		basicCell.value = 1.0;
		assignmentSolution.push_back(basicCell);
	}
}

//generate basic solution for equivalent transportation problem
void BranchAndBoundSolver::generateBasicTransportationSolution() {
	int size = costTableau.size();
	std::multimap<int, int> assignmentSolution;
	std::multimap<int, int> transportationSolution;
	transportationSolution = assignmentSolution;
	std::cout << "Row scanning for the acyclic connected graph." << std::endl;
	std::vector<int> columnCells;
	std::vector<int> rowCells;
	for (int i = 0; i < size; i++) {
		columnCells.push_back(1);
	}
	for (int i = 0; i < size; i++) {
		auto col = assignmentSolution.find(i);
		double val1 = INFINITY;
		double val2 = 0;
		int indx = 0;
		for (int j = 1; j < size; j++) {
			auto it = assignmentSolution.find(i);
			j != i ? val2 = costTableau[j][(*it).second] : val2 = INFINITY;
			if (val2 < val1) {
				val1 = val2;
				indx = j;
			}
		}
		auto itLow = transportationSolution.lower_bound(indx);
		auto itUp = transportationSolution.upper_bound(indx);
		int count = 0;
		for (auto& it = itLow; it != itUp; it++) {
			count += 1;
		}
		if (count <= 1 && indx != i) {
			//exclude any assignment sol?
			transportationSolution.insert(std::pair<int, int>(indx, (*col).second));
			int numColCells = columnCells.at(indx);
			columnCells[(*col).second] = numColCells + 1;
		}
	}
	for (auto& it : transportationSolution) {
		std::cout << "\n" << it.first << " " << it.second << std::endl;
	}
	//total number of cell assigned to the basic solution
	int numOfBasicCells = 0;
	for (int i = 0; i < size; i++) {
		numOfBasicCells += columnCells.at(i);
	}
	//populate rowcell vector
	for (int i = 0; i < size; i++) {
		auto itLow = transportationSolution.lower_bound(i);
		auto itUp = transportationSolution.upper_bound(i);
		int val = 0;
		for (auto& it = itLow; it != itUp; it++) {
			val += 1;
		}
		rowCells.push_back(val);
	}
	//columnMap
	std::multimap<int, int> columnMap;
	for (auto& it : transportationSolution) {
		columnMap.insert(std::pair<int, int>(it.second, it.first));
	}
	//column scanning for the acyclic connected graph
	if (numOfBasicCells < (2 * size - 1)) {
		std::cout << "\nStart of column scanning" << std::endl;
		for (int i = 0; i < size; i++) {
			if (columnCells.at(i) == 1) {
				//find row ID, You still has column ID
				auto it = columnMap.find(i);
				int rowID = (*it).second;
				if (rowCells.at(rowID) == 1) {
					int col = 0;
					double val3 = INFINITY;
					double val4 = 0;
					for (int k = 0; k < size; k++) {
						k != i ? val4 = costTableau[rowID][k] : val4 = INFINITY;
						if (val4 < val3) {
							val3 = val4;
							col = k;
						}
					}
					transportationSolution.insert(std::pair<int, int>(rowID, col));
					columnMap.insert(std::pair<int, int>(col, rowID));
					columnCells[i] = columnCells.at(i) + 1;
					rowCells[rowID] = rowCells.at(rowID) + 1;
					numOfBasicCells += 1;
					if (numOfBasicCells == (2 * size - 1)) {
						break;
					}
				}
			}
		}
	}
	std::cout << "\nThe final transportation basic solution is : " << std::endl;
	for (auto& it : transportationSolution) {
		std::cout << "\n" << it.first << " " << it.second << std::endl;
	}
	/*
	//print rows with number of contained basic cells
	for (int i = 0; i < size; i++) {
		std::cout << "\nRow number : " << i << " Number of contained basic cells : " << rowCells.at(i);
	}
	//print columns with number of contained basic cells
	for (int i = 0; i < size; i++) {
		std::cout << "\nColumn number : " << i << " Number of contained basic cells : " << columnCells.at(i);
	}
	*/
}

//solve equivalent transportation problem using Cost operator
void BranchAndBoundSolver::solveTransportationProblem() {


}

//returns tsp solution if any exists
TSPSolution BranchAndBoundSolver::getTSPSolution() {
	return tspSolution;
}

//run B&B algorithm
void BranchAndBoundSolver::runBranchAndBoundSolver() {


}

