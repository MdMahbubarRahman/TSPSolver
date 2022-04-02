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
	int rowId;
	int colID;
	double value;
};

enum Direction {
	East, West, North, South, Origin
};

enum CellType {
	Getter, Giver
};

struct AllocatedCell {
	Direction prevDirection;
	Direction postDirection;
	BasicCell cellProperty;
	CellType cellType;
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
	void updateDualSolution();
	void findMaxDeltaAndEnteringCell(int p, int q);
	void generateCycleAndUpdateBasicSolution();
	void generateAcyclicConnectedGraphSolution();
	void updateWeakerLowerBound();
	void runCostOperator();
};

#endif






