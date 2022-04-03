#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <set>
#include <Queue>
#include <map>

#ifndef OPERATORTHEORY_H
#define OPERATORTHEORY_H

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

class OperatorTheory {
private:
	std::set<int> Ip;
	std::set<int> Iq;
	std::set<int> Jp;
	std::set<int> Jq;
	std::list<BasicCell> basicSolution;
	std::vector<double> rowWiseDualSolution;
	std::vector<double> columnWiseDualSolution;
	std::vector<std::vector<double>> costTableau;

public:
	OperatorTheory();
	OperatorTheory(std::list<BasicCell> basicSolution, std::vector<std::vector<double>> costTableau);
	OperatorTheory(std::list<BasicCell> basicSolution, std::vector<double> rowWiseDualSolution, std::vector<double> columnWiseDualSolution, std::vector<std::vector<double>> costTableau);
	OperatorTheory(const OperatorTheory& opThr);
	void generateInitialDualSolution();
	void scanningRoutine(int p, int q);
	void updateDualSolution(double val);
	double findMaxDeltaAndEnteringCell(int p, int q);
	void generateCycleAndUpdateBasicSolution(int enteringCellRowId, int enteringCellColID);
	void generateAcyclicConnectedGraphSolution();
	void updateWeakerLowerBound();
	void runCostOperator();
};

#endif






