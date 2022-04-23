#include "BranchAndBound.h"


//default constructor
BranchAndBoundSolver::BranchAndBoundSolver() {
	totalTime = 0;
	maxTimeLimit = 1000;//1000sec
	lowerBound = -INFINITY;
	weakerLowerBound = -INFINITY;
	upperBound = INFINITY;
}

//constructor
BranchAndBoundSolver::BranchAndBoundSolver(std::vector<int> initialTourOfCities, std::vector<std::vector<double>> wholeCostMatrix) {
	oldListOfCities = initialTourOfCities;
	std::cout << "\nPopulate new city to old city map" << std::endl;
	for (int i = 0; i < oldListOfCities.size(); i++) {
		newCityToOldCityMap.insert(std::pair<int, int>(i, oldListOfCities.at(i)));
		newListOfCities.push_back(i);
	}
	std::cout << "\nPopulate the cost matrix for the current cities." << std::endl;
	double val = 0;
	for (int i = 0; i < oldListOfCities.size(); i++) {
		std::vector<double> costVec;
		for (int j = 0; j < oldListOfCities.size(); j++) {
			i == j ? val = INFINITY : val = wholeCostMatrix[oldListOfCities.at(i)][oldListOfCities.at(j)];
			costVec.push_back(val);
		}
		costTableau.push_back(costVec);
	}
	//populate variables with initial values
	totalTime = 0;
	maxTimeLimit = 1000;//1000sec
	lowerBound = -INFINITY;
	weakerLowerBound = -INFINITY;
	upperBound = INFINITY;
}

//copy constructor
BranchAndBoundSolver::BranchAndBoundSolver(const BranchAndBoundSolver& depthFirstBBSolver) {
	lowerBound = depthFirstBBSolver.lowerBound;
	weakerLowerBound = depthFirstBBSolver.weakerLowerBound;
	upperBound = depthFirstBBSolver.upperBound;
	totalTime = depthFirstBBSolver.totalTime;
	maxTimeLimit = depthFirstBBSolver.maxTimeLimit;
	bbTree = depthFirstBBSolver.bbTree;
	bStatus = depthFirstBBSolver.bStatus;
	incumbent = depthFirstBBSolver.incumbent;
	tourSolutions = depthFirstBBSolver.tourSolutions;
	tspSolution = depthFirstBBSolver.tspSolution;
	costTableau = depthFirstBBSolver.costTableau;
	oldListOfCities = depthFirstBBSolver.oldListOfCities;
	newListOfCities = depthFirstBBSolver.newListOfCities;
	newCityToOldCityMap = depthFirstBBSolver.newCityToOldCityMap;
	assignmentSolution = depthFirstBBSolver.assignmentSolution;
	transportationBasicSolution = depthFirstBBSolver.transportationBasicSolution;
	routeLists = depthFirstBBSolver.routeLists;
}

//solves assignment problem using Hungarian algorithm to get lower bound of tsp
void BranchAndBoundSolver::solveAssignmentProblem() {
	HungarianAlg hungAlg(newListOfCities, costTableau);
	hungAlg.runHungarianAlg();
	routeLists = hungAlg.getListOfRoutes();
	assignmentSolution = hungAlg.getAssignmentSolution();
}

//generate basic solution for equivalent transportation problem
void BranchAndBoundSolver::generateBasicTransportationSolution() {
	double M = INFINITY;
	int size = costTableau.size();
	std::multimap<int, int> transportationSolution;
	for (auto& it : assignmentSolution) {
		transportationSolution.insert(std::pair<int, int>(it.first, it.second));
	}
	std::vector<int> columnCells;
	std::vector<int> rowCells;
	for (int i = 0; i < size; i++) {
		columnCells.push_back(1);
	}
	std::cout << "Row scanning for the acyclic connected graph." << std::endl;
	for (int i = 0; i < size; i++) {
		auto col = assignmentSolution.find(i);
		double val2 = 0;
		for (int j = 1; j < size; j++) {
			if (j != i) {
				val2 = costTableau[j][(*col).second];
				if (val2 < M) {
					auto itLow = transportationSolution.lower_bound(j);
					auto itUp = transportationSolution.upper_bound(j);
					int count = 0;
					for (auto& it = itLow; it != itUp; it++) {
						count += 1;
					}
					if (count <= 1) {
						transportationSolution.insert(std::pair<int, int>(j, (*col).second));
						columnCells[(*col).second] = 2;
						break;
					}
				}
			}
		}
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
	int counter = 0;
	std::multimap<int, int> columnMap;
	for (auto& it : transportationSolution) {
		columnMap.insert(std::pair<int, int>(it.second, it.first));
		counter += 1;
	}
	//column scanning for the acyclic connected graph
	if (counter < (2 * size - 1)) {
		std::cout << "\nStart of column scanning" << std::endl;
		for (int i = 0; i < size; i++) {
			if (columnCells.at(i) == 1) {
				//find row ID, You still has column ID
				auto it = columnMap.find(i);
				int rowID = (*it).second;
				if (rowCells.at(rowID) == 1) {
					int col = 0;
					double val3 = M;
					double val4 = 0;
					for (int k = 0; k < size; k++) {
						k != i ? val4 = costTableau[rowID][k] : val4 = M;
						if (val4 <= val3) {
							val3 = val4;
							col = k;
						}
					}
					transportationSolution.insert(std::pair<int, int>(rowID, col));
					columnMap.insert(std::pair<int, int>(col, rowID));
					columnCells[i] = columnCells.at(i) + 1;
					rowCells[rowID] = rowCells.at(rowID) + 1;
					counter += 1;
					if (counter == (2 * size - 1)) {
						break;
					}
				}
			}
		}
	}
	//keep the assigned value in a map for later use
	std::map<int, std::map<int, double>> assignedValueMap;
	for (auto& it : assignmentSolution) {
		assignedValueMap[it.first][it.second] = 1.0;
	}
	//populate the transportation basic solution
	for (auto& it : transportationSolution) {
		BasicCell base = BasicCell();
		base.rowID = it.first;
		base.colID = it.second;
		assignedValueMap[it.first][it.second] == 1.0 ? base.value = 1.0 : base.value = 0.0;
		transportationBasicSolution.push_back(base);
	}
}

//initialize the branch and bound tree
void BranchAndBoundSolver::initBranchAndBoundTree() {
	solveAssignmentProblem();
	generateBasicTransportationSolution();
	OperatorTheory operThry = OperatorTheory(transportationBasicSolution, costTableau);
	operThry.runCostOperatorForGeneratingRootNodes();
	//operThry.showChildNodes();
	std::list<Node> rootNodes = operThry.getChildNodes();
	//populate branch and bound tree
	int counter = 0;
	bbTree.numOfNodePrunedByBound = 0;
	bbTree.numOfNodePrunedByIntegrality = 0;
	bbTree.numOfNodePrunedWithErrorsOrInfeasibility = 0;
	for (auto& it : rootNodes) {
		bbTree.branchNodes.push_back(it);
		counter += 1;
	}
	bbTree.currentNumOfNodes = counter;
}

//returns tsp solution if any exists
TSPSolution BranchAndBoundSolver::getTSPSolution() {
	return tspSolution;
}

//prune node consisting of a single tour
void BranchAndBoundSolver::pruneNodeByIntegrality() {
	bool singleTourExist = false;
	std::list<Node>::iterator nodeId;
	for (auto it = bbTree.branchNodes.begin(); it != bbTree.branchNodes.end(); it++) {
		if ((*it).isSolutionATour == true) {
			nodeId = it;
			singleTourExist = true;
		}
	}
	if (singleTourExist == true) {
		std::map<int, int> tour;
		std::vector<int> route;
		int sourceNode = 0;
		int destinationNode = 0;
		int routeStartNode = 0;
		std::map<int, int> newAssignments;
		for (auto& it : (*nodeId).basicSolution) {
			if (it.value == 1.0) {
				newAssignments.insert(std::pair<int, int>(it.rowID, it.colID));
			}
		}
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
					break;
				}
				else {
					i += 1;
				}
			}
		}
		TourSolution sol = TourSolution();
		sol.objValue = (*nodeId).weakerLowerBound;
		sol.tourSolution = route;
		tourSolutions.push_back(sol);
		bbTree.branchNodes.erase(nodeId);
	}
}

//node selection based on best weaker lower bound
Node BranchAndBoundSolver::selectNodeBasedOnBestWeakerLowerBound() {
	std::list<Node>::iterator nodeId;
	double val = INFINITY;
	for (auto it = bbTree.branchNodes.begin(); it != bbTree.branchNodes.end(); it++) {
		if ((*it).weakerLowerBound < val) {
			val = (*it).weakerLowerBound;
			nodeId = it;
		}
	}
	Node node = (*nodeId);
	bbTree.branchNodes.erase(nodeId);
	return node;
}

//solve a single node using cost operator
void BranchAndBoundSolver::solveNodeByCostOperator(Node node) {
	OperatorTheory opThry = OperatorTheory(node);
	opThry.runCostOperatorForSolvingANode();
	std::list<Node> nodes = opThry.getChildNodes();
	for (auto& it : nodes) {
		bbTree.branchNodes.push_back(it);
	}
}

//run B&B algorithm
void BranchAndBoundSolver::runBranchAndBoundSolver() {
	initBranchAndBoundTree();
	while (!bbTree.branchNodes.empty()) {
		pruneNodeByIntegrality();
		Node node = selectNodeBasedOnBestWeakerLowerBound();
		solveNodeByCostOperator(node);
		//need works
	}
}

