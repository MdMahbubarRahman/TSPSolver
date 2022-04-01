#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <Queue>
#include <map>


#ifndef OPERATORTHEORY_H
#define OPERATORTHEORY_H


class OperatorTheory{
private:
	int tableauSize = 0;
	std::set<int> Ip;
	std::set<int> Iq;
	std::set<int> Jp;
	std::set<int> Jq;
	std::vector<double> currentRowDualVariables;
	std::vector<double> currentColumnDualVariables;
	std::vector<double> transformedRowDualVariables;
	std::vector<double> transformedColumnDualVariables;
	std::multimap<int, int> basicSolnCells; //map columns to rows for the basic solution
	std::vector<std::vector<double>> tpTableau;

public:
	OperatorTheory();
	OperatorTheory(std::multimap<int, int> basicSolnCells, std::vector<std::vector<double>> tpTableau);
	OperatorTheory(const OperatorTheory & opThr);
	void scanningRoutine();
	void updateDualVariables();
	void findMaxDeltaAndEnteringCell();
	void generateCycleAndUpdateBasicSolution();
	void generateAcyclicConnectedGraphSolution();
	void updateWeakerLowerBound();
	void runOperatorTheory();
};

#endif




