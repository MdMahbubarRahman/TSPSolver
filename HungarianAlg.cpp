#include "HungarianAlg.h"


//default constructor
HungarianAlg::HungarianAlg() {
	reducedCost = 0.0;
	baseValue = 0.0;
	initSolCost = 0.0;
	finalSolCost = 0.0;
	tspOptimal = false;
	assignmentOptimal = false;
	reducedCostMatrixSize = 0;
}


//copy constructor
HungarianAlg::HungarianAlg(const HungarianAlg& hunAlg) {
	reducedCost = hunAlg.reducedCost;
	baseValue = hunAlg.baseValue;
	initSolCost = hunAlg.initSolCost;
	finalSolCost = hunAlg.finalSolCost;
	tspOptimal = hunAlg.tspOptimal;
	reducedCostMatrixSize = hunAlg.reducedCostMatrixSize;
	struct IntersectionPoint {
		int row;
		int column;
		IntersectionPoint(int a, int b) {
			row = a;
			column = b;
		}
	};
	boxPoints = hunAlg.boxPoints;
	assignmentOptimal = hunAlg.assignmentOptimal;
	listOfRoutes = hunAlg.listOfRoutes;
	initialTsp = hunAlg.initialTsp;//tsp presentation node/node/.../node/depotnode
	incumbentTsp = hunAlg.incumbentTsp;
	optimalTsp = hunAlg.optimalTsp;
	costMatrix = hunAlg.costMatrix;
	reducedCostMatrix = hunAlg.reducedCostMatrix;
}


//constructor
HungarianAlg::HungarianAlg(std::vector<int> initTsp, std::vector<std::vector<double>> cMatrix) {
	reducedCost = 0.0;
	baseValue = 0.0;
	initSolCost = 0.0;
	finalSolCost = 0.0;
	tspOptimal = false;
	assignmentOptimal = false;
	initialTsp = initTsp;
	costMatrix = cMatrix;
	//populate reduced const matrix
	reducedCostMatrixSize = initialTsp.size();
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		std::vector<double> costVec;
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			double val = 0;
			i == j ? val = INFINITY : val = costMatrix[initialTsp.at(i)][initialTsp.at(j)];
			costVec.push_back(val);
		}
		reducedCostMatrix.push_back(costVec);
	}
	//initial tsp cost
	for (int i = 0; i < reducedCostMatrixSize - 1; i++) {
		initSolCost += costMatrix[initialTsp.at(i)][initialTsp.at(i + 1)];
	}
	initSolCost += costMatrix[initialTsp.at(reducedCostMatrixSize - 1)][initialTsp.at(0)];
}


//performs row reduction
void HungarianAlg::performRowColumnReduction() {
	std::cout << "\nShow reduced cost matrix ." << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			std::cout << reducedCostMatrix[i][j] << " ";
		}
		std::cout << ";" << std::endl;
	}
	//perform row reduction
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		double minVal = INFINITY;
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			if (minVal > reducedCostMatrix[i][j]) {
				minVal = reducedCostMatrix[i][j];
			}
		}
		//subtract min val from each the corresponding row
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			reducedCostMatrix[i][j] = reducedCostMatrix[i][j] - minVal;
		}
		reducedCost += minVal;
	}
	std::cout << "\nShow reduced cost matrix after row reduction." << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			std::cout << reducedCostMatrix[i][j] << " ";
		}
		std::cout << ";" << std::endl;
	}


	//perform column reduction
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		double minVal = INFINITY;
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			if (minVal > reducedCostMatrix[j][i]) {
				minVal = reducedCostMatrix[j][i];
			}
		}
		//subtract min val from each the corresponding column
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			reducedCostMatrix[j][i] = reducedCostMatrix[j][i] - minVal;
		}
		reducedCost += minVal;
	}
	reducedCost == initSolCost ? tspOptimal = true : tspOptimal = false;
	std::cout << "\nShow reduced cost matrix after row column reduction." << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			std::cout << reducedCostMatrix[i][j] << " ";
		}
		std::cout << ";" << std::endl;
	}
}

//perform row scanning
void HungarianAlg::rowScanning(std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints) {
	//row scaning
	std::cout << "\nRow scanning in progress!" << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		auto it = coveredRows.find(i);
		if (it == coveredRows.end()) {
			int numOfZeros = 0;
			int colPos = 0;
			for (int j = 0; j < reducedCostMatrixSize; j++) {
				auto iter = coveredColumns.find(j);
				if (iter == coveredColumns.end() && reducedCostMatrix[i][j] == 0) {
					numOfZeros += 1;
					colPos = j;
				}
			}
			if (numOfZeros == 1) {
				coveredColumns.insert(colPos);
				boxPoints.insert(std::pair<int, int>(i, colPos));
				std::cout << "The new box point is : " << i << " , " << colPos << std::endl;
			}
		}
	}
}

//perform column scanning
void HungarianAlg::columnScanning(std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints) {
	//column scaning
	std::cout << "\nColumn scanning in progress!" << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		auto it = coveredColumns.find(i);
		if (it == coveredColumns.end()) {
			int numOfZeros = 0;
			int rowPos = 0;
			for (int j = 0; j < reducedCostMatrixSize; j++) {
				auto iter = coveredRows.find(j);
				if (iter == coveredRows.end() && reducedCostMatrix[j][i] == 0) {
					numOfZeros += 1;
					rowPos = j;
				}
			}
			if (numOfZeros == 1) {
				coveredRows.insert(rowPos);
				boxPoints.insert(std::pair<int, int>(rowPos, i));
				std::cout << "The new box point is : " << rowPos << " , " << i << std::endl;
			}
		}
	}
}

//scan for the number of uncovered zeros
int HungarianAlg::uncoveredZeroScanning(std::set<int>& coveredRows, std::set<int>& coveredColumns) {
	//check if there is any uncovered zero
	std::cout << "\nNumber of uncovered zero scanning in progress!" << std::endl;
	std::map<int, int> numZeroInRow;
	std::map<int, int> numZeroInCol;
	std::cout << "\nUpdate number of uncovered zeros in each row." << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		auto itr = coveredRows.find(i);
		if (itr == coveredRows.end()) {
			int numZero = 0;
			for (int j = 0; j < reducedCostMatrixSize; j++) {
				auto it = coveredColumns.find(j);
				if (it == coveredColumns.end() && reducedCostMatrix[i][j] == 0) {
					numZero += 1;
				}
			}
			numZeroInRow.insert(std::pair<int, int>(i, numZero));
		}
	}
	std::cout << "\nUpdate number of uncovered zeros in each column." << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		auto itr = coveredColumns.find(i);
		if (itr == coveredColumns.end()) {
			int numZero = 0;
			for (int j = 0; j < reducedCostMatrixSize; j++) {
				auto it = coveredRows.find(j);
				if (it == coveredRows.end() && reducedCostMatrix[j][i] == 0) {
					numZero += 1;
				}
			}
			numZeroInCol.insert(std::pair<int, int>(i, numZero));
		}
	}
	//Number of uncovered zeros
	int numOfZeros = 0;
	for (auto& it : numZeroInCol) {
		numOfZeros += it.second;
	}
	std::cout << "\nNumber of zeros : " << numOfZeros << std::endl;

	return numOfZeros;
}


//perform diagonal selection in case more than one zero in a row
void HungarianAlg::diagonalSelectionAdv(std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints) {
	//1. first choose a zero in a row who has multiple zeros
	//2. //then if the sum of the number of covered rows and the number of covered columns is less than the size of the reduced matrix,
	//2. look for a second zero in a row which is positioned in diagonally opposite location
	//3. if the second zero containing row has two or more zeros, go to step 2

	//map of uncovered rows and columns
	std::map<int, int> uncoveredRowMapWithZeros;
	std::map<int, int> uncoveredColumnMapWithZeros;
	int iter1 = 0;
	int iter2 = 0;
	//uncovered rows with multiple zeros
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		auto it = coveredRows.find(i);
		if (it == coveredRows.end()) {
			int numZero = 0;
			for (int j = 0; j < reducedCostMatrixSize; j++) {
				auto itt = coveredColumns.find(j);
				if (itt == coveredColumns.end() && reducedCostMatrix[i][j] == 0) {
					numZero += 1;
				}
			}
			if (numZero > 1) {
				uncoveredRowMapWithZeros.insert(std::pair<int, int>(iter1, i));
				iter1 += 1;
			}
		}
	}
	//uncovered columns with multiple zeros
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		auto it = coveredColumns.find(i);
		if (it == coveredColumns.end()) {
			int numZero = 0;
			for (int j = 0; j < reducedCostMatrixSize; j++) {
				auto itt = coveredRows.find(j);
				if (itt == coveredRows.end() && reducedCostMatrix[j][i] == 0) {
					numZero += 1;
				}
			}
			if (numZero > 1) {
				uncoveredColumnMapWithZeros.insert(std::pair<int, int>(iter2, i));
				iter2 += 1;
			}
		}
	}
	std::cout << "\nDiagonal selection in progress!" << std::endl;
	while (!uncoveredRowMapWithZeros.empty()) {
		bool flag1 = true;
		bool flag2 = true;
		int rowID = 10000;
		int colID = 10000;
		int difVal = 0;
		//perform first zero selection
		for (auto& it : uncoveredRowMapWithZeros) {
			for (auto& itt : uncoveredColumnMapWithZeros) {
				if (reducedCostMatrix[it.second][itt.second] == 0) {
					boxPoints.insert(std::pair<int, int>(it.second, itt.second));
					//coveredRows.insert(it.second);
					coveredColumns.insert(itt.second);
					std::cout << "The new box point is : " << it.second << " , " << itt.second << std::endl;
					rowID = it.first;
					colID = itt.first;
					difVal = itt.first - it.first;
					flag1 = false;
					break;
				}
			}
			if (flag1 == false) {
				break;
			}
		}
		rowID != 10000 ? uncoveredRowMapWithZeros.erase(rowID) : rowID = 10000;
		//perform diagonal opposite zero selection
		if (!uncoveredRowMapWithZeros.empty()) {
			for (auto& it : uncoveredRowMapWithZeros) {
				if (it.first > rowID) {
					for (auto& itt : uncoveredColumnMapWithZeros) {
						int val = itt.first - it.first;
						if (itt.first > colID && val == difVal && reducedCostMatrix[it.second][itt.second] == 0) {
							boxPoints.insert(std::pair<int, int>(it.second, itt.second));
							//coveredRows.insert(it.second);
							coveredColumns.insert(itt.second);
							rowID = it.first;
							std::cout << "The new box point is : " << it.second << " , " << itt.second << std::endl;
							flag2 = false;
							break;
						}
					}
				}
				if (flag2 == false) {
					break;
				}
			}
		}
		rowID != 10000 ? uncoveredRowMapWithZeros.erase(rowID) : rowID = 10000;
	}
}



//cover zero values and check for optimal assignment problem solution
void HungarianAlg::coverMinimumValueAssignmentsAndCheckOptimality(bool& assignmentOptimal, std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints, std::list<IntersectionPoint>& twiceCoveredPoints) {
	int numZero1 = 1000;
	int numZero2 = 0;
	int iter = 0;
	while (numZero1 > 0) {
		//row scaning
		rowScanning(coveredRows, coveredColumns, boxPoints);
		//column scaning
		columnScanning(coveredRows, coveredColumns, boxPoints);
		//check if there is any uncovered zero
		numZero2 = uncoveredZeroScanning(coveredRows, coveredColumns);
		if (numZero1 == numZero2) {
			iter += 1;
		}
		else {
			numZero1 = numZero2;
			iter = 0;
		}
		if (numZero1 > 1 && iter == 1) {
			//diagonal selection
			diagonalSelectionAdv(coveredRows, coveredColumns, boxPoints);
			iter = 0;
		}
	}
	
	//show intersection points
	std::cout << "\nShow the intersection points" << std::endl;
	for (auto i : coveredRows) {
		for (auto j : coveredColumns) {
			IntersectionPoint point = IntersectionPoint(i, j);
			std::cout << "Row : " << i << " Column : " << j << std::endl;
			twiceCoveredPoints.push_back(point);
		}
	}
	
	//check assignment optimality
	std::cout << "\nNumber of rows covered : " << coveredRows.size() << ", Number of columns covered : " << coveredColumns.size() << std::endl;
	//std::cout << "\nReduced cost matrix size : " << reducedCostMatrixSize << std::endl;
	if ((coveredColumns.size() + coveredRows.size()) == reducedCostMatrixSize) {
		assignmentOptimal = true;
		std::cout << "\nAssignment optimal : " << assignmentOptimal << std::endl;
	}
}

//performs row column reduction for uncovered nonzero cells
void HungarianAlg::performRowColumnReductionForUncoveredCells(bool& assignmentOptimal, std::set<int>& coveredRows, std::set<int>& coveredColumns, std::map<int, int>& boxPoints, std::list<IntersectionPoint>& twiceCoveredPoints) {
	//perform row reduction
	std::set<int> uncoveredRows, uncoveredColumns;
	//update uncovered rows
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		auto it = coveredRows.find(i);
		if (it == coveredRows.end()) {
			uncoveredRows.insert(i);
		}
	}
	//update uncovered columns
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		auto it = coveredColumns.find(i);
		if (it == coveredColumns.end()) {
			uncoveredColumns.insert(i);
		}
	}
	//find minimum value from the uncovered cells
	double minVal = INFINITY;
	for (auto it : uncoveredRows) {
		for (auto itt : uncoveredColumns) {
			if (minVal > reducedCostMatrix[it][itt]) {
				minVal = reducedCostMatrix[it][itt];
			}
		}
	}
	//subtract min value from all uncovered cells
	for (auto it : uncoveredRows) {
		for (auto itt : uncoveredColumns) {
			reducedCostMatrix[it][itt] = reducedCostMatrix[it][itt] - minVal;
		}
		reducedCost += minVal;
	}
	//add min value to all intersection points
	for (auto& point : twiceCoveredPoints) {
		reducedCostMatrix[point.row][point.column] = reducedCostMatrix[point.row][point.column] + minVal;
	}
	/*
	std::cout << "\nShow reduced cost matrix after updating uncovered cells." << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			std::cout << reducedCostMatrix[i][j] << " ";
		}
		std::cout << ";" << std::endl;
	}
	*/
}

//generates routes from assignment solution
void HungarianAlg::populateListOfRoutes(std::map<int, int> boxPoints) {
	std::cout << "\nPrint the contents of the newAssignment" << std::endl;
	for (auto& it : boxPoints) {
		std::cout << it.first << " " << it.second << std::endl;
	}
	std::cout << "\nThe routes can be formed from the above optimal assignments!" << std::endl;
	//start route from depot
	/*
	std::vector<int> route;
	int sourceNode = 0;
	int destinationNode = 0;
	int routeStartNode = 0;
	std::map<int, int> newAssignments = boxPoints;
	int i = 0;
	while (!newAssignments.empty()) {
		if (i == 0) {
			auto it = newAssignments.begin();
			routeStartNode = (*it).first;
			destinationNode = (*it).second;
			route.push_back(destinationNode);
			newAssignments.erase(routeStartNode);
			i += 1;
		}
		else {
			sourceNode = destinationNode;
			destinationNode = newAssignments[sourceNode];
			route.push_back(destinationNode);
			newAssignments.erase(sourceNode);
			if (destinationNode == routeStartNode) {
				listOfRoutes.push_back(route);
				route.clear();
				i = 0;
			}
			else {
				i += 1;
			}
		}
	}
	*/
}

//finds the number of routes generated from assignment solution
void HungarianAlg::checkListOfRoutes(std::map<int, int> boxPoints, std::list<std::vector<int>>& routeList) {
	//start route from depot
	routeList.clear();
	std::vector<int> route;
	int sourceNode = 0;
	int destinationNode = 0;
	int routeStartNode = 0;
	std::map<int, int> newAssignments = boxPoints;
	int i = 0;
	std::cout << "\nStart routes construction." << std::endl;
	while (!newAssignments.empty()) {
		if (i == 0) {
			auto it = newAssignments.begin();
			routeStartNode = (*it).first;
			destinationNode = (*it).second;
			route.push_back(destinationNode);
			newAssignments.erase(routeStartNode);
			i += 1;
		}
		else {
			sourceNode = destinationNode;
			destinationNode = newAssignments[sourceNode];
			route.push_back(destinationNode);
			newAssignments.erase(sourceNode);
			if (destinationNode == routeStartNode) {
				routeList.push_back(route);
				route.clear();
				i = 0;
			}
			else {
				i += 1;
			}
		}
	}
	std::cout << "\nEnd routes construction." << std::endl;
}


//solves relaxed tsp problem which is an assignment problem and generates lower bound of TSP
void HungarianAlg::runHungarianAlg() {
	std::set<int> coveredRows;
	std::set<int> coveredColumns;
	std::list<IntersectionPoint> twiceCoveredPoints;
	performRowColumnReduction();
	if (tspOptimal == true) {
		std::cout << "\nThe initial tsp solution is optimal! No need to look for optimal tsp solution." << std::endl;
	}
	else {
		while (assignmentOptimal != true) {
			coverMinimumValueAssignmentsAndCheckOptimality(assignmentOptimal, coveredRows, coveredColumns, boxPoints, twiceCoveredPoints);
			if (assignmentOptimal != true) {
				performRowColumnReductionForUncoveredCells(assignmentOptimal, coveredRows, coveredColumns, boxPoints, twiceCoveredPoints);
				coveredColumns.clear();
				coveredRows.clear();
				boxPoints.clear();
				twiceCoveredPoints.clear();
			}
			else {
				std::cout << "\nOptimum assignment is found!" << std::endl;
				populateListOfRoutes(boxPoints);
			}
		}
	}
}

//shows results
void HungarianAlg::showAssignmentSolution() {
	std::cout << "\nShow the reduced cost matrix : " << std::endl;
	for (int i = 0; i < reducedCostMatrixSize; i++) {
		for (int j = 0; j < reducedCostMatrixSize; j++) {
			std::cout << reducedCostMatrix[i][j] << " ";
		}
		std::cout << ";" << std::endl;
	}
	std::cout << "\nShow assignments : " << std::endl;
	for (auto& it : boxPoints) {
		std::cout << it.first << " --> " << it.second << std::endl;
	}
	std::cout << "\nShow the routes : " << std::endl;
	for (auto& it : listOfRoutes) {
		for (auto itt : it) {
			std::cout << itt << " ";
		}
		std::cout << ";" << std::endl;
	}
}

