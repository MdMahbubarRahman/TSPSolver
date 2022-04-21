#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <set>
#include <Queue>
#include <map>

#ifndef OPERATORTHEORY_H
#define OPERATORTHEORY_H

enum CostOperatorStatus {
	cTerminatedSuccessfully, cTerminatedWithErrors
};

struct BasicCell {
	int rowID;
	int colID;
	double value;
};

struct Cell {
	int rowID;
	int colID;
};

enum Direction {
	East, West, North, South, Origin
};

enum CellType {
	Getter, Giver
};

struct AllocatedCell {
	BasicCell cellProperty;
	CellType cellType;
	Cell prevCell;
	Cell postCell;
};

struct Node {
	int nodeIndex;
	int parentNodeIndex;
	bool isSolutionATour;
	double weakerLowerBound;
	double lowerBound;
	double delta;
	std::set<int> Ip;
	std::set<int> Iq;
	std::set<int> Jp;
	std::set<int> Jq;
	Cell branchOnCell;
	Cell enteringCell;
	std::list<BasicCell> basicSolution;
	std::vector<double> rowDualSolution;
	std::vector<double> columnDualSolution;
	std::vector<std::vector<double>> costTableau;
};


class OperatorTheory {
private:
	int nodeIndex;
	int parentNodeIndex;
	bool isSolutionATour;
	double weakerLowerBound;
	double lowerBound;
	double delta;
	std::set<int> Ip;
	std::set<int> Iq;
	std::set<int> Jp;
	std::set<int> Jq;
	Cell branchOnCell;
	Cell enteringCell;
	std::list<BasicCell> basicSolution;
	std::vector<double> rowWiseDualSolution;
	std::vector<double> columnWiseDualSolution;
	std::vector<std::vector<double>> costTableau;
	std::list<Node> childNodes;
	CostOperatorStatus cStatus;
	std::list<std::vector<int>> clistOfRoutes;
	std::list<AllocatedCell> cycleOfCells;

public:
	OperatorTheory();
	OperatorTheory(Node node);
	OperatorTheory(std::list<BasicCell> basicSolution, std::vector<std::vector<double>> costTableau);
	OperatorTheory(std::list<BasicCell> basicSolution, std::vector<double> rowWiseDualSolution, std::vector<double> columnWiseDualSolution, std::vector<std::vector<double>> costTableau);
	OperatorTheory(const OperatorTheory& opThr);
	void generateInitialDualSolution();
	void scanningRoutine(int p, int q);
	void updateDualSolution(double val);
	void findMaxDeltaAndEnteringCell(int p, int q);
	void generateCycleAndUpdateBasicSolution(int enteringCellRowId, int enteringCellColID);
	void generateListOfRoutes(std::map<int, int> boxPoints);
	void updateBounds(std::map<int, int> boxPoints);
	void runCostOperatorForGeneratingRootNodes();
	void runCostOperatorForSolvingANode();
	std::list<Node> getChildNodes();
	void showChildNodes();
};

#endif






