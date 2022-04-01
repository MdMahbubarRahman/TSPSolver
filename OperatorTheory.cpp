#include "OperatorTheory.h"

OperatorTheory::OperatorTheory() {

}

OperatorTheory::OperatorTheory(std::multimap<int, int> basicCells, std::vector<std::vector<double>> tableau) {

}

OperatorTheory::OperatorTheory(const OperatorTheory& opThr) {

}

void OperatorTheory::scanningRoutine() {
	std::vector<int> rowLabel;
	std::vector<int> columnLabel;
	int cellRowID = 0;
	int cellColumnID = 0;
	int size = tableauSize;
	bool stopFlag = false;
	int iter = 0;
	//Populate rowlabel and columnlabel with  default values;
	for (int i = 0; i < size; i++) {
		rowLabel.push_back(0);
		columnLabel.push_back(0);
	}
	//Run scanning routine
	while (!stopFlag) {
		iter += 1;
		bool newColumnLabel = false;
		bool newRowLabel = false;
		if (iter == 1) {
			//0. Label with 1 the columns of basis cells in row p cexept for basis cell (p,q); label row p with 2
			for (auto it = basicSolnCells.begin(); it != basicSolnCells.end(); it++) {
				if ((*it).first == cellRowID && (*it).second != cellColumnID) {
					columnLabel[(*it).second] = 1;
					newColumnLabel = true;
				}
			}
			rowLabel[cellRowID] = 2;
		}
		//1. If no new columns were labelled 1 on the preceding step go to 5. otherwise go to 2.
		if (newColumnLabel != true) {
			break;
		}
		//2. For each column labelled with 1, label with 1 each row containing a basis cell in that column unless the row is 
		//labelled with 2; then change the label of the column to 2.
		for (int i = 0; i < size; i++) {
			if (columnLabel[i] == 1) {
				for (auto it = basicSolnCells.begin(); it != basicSolnCells.end(); ++it) {
					if ((*it).second == i && rowLabel[(*it).first] != 2) {
						rowLabel[(*it).first] = 1;
						newRowLabel = true;
					}
				}
				columnLabel[i] = 2;
			}
		}
		//3. If no new rows were labelled with 1 on the preceeding step go to 5. Otherwise go to 4.
		if (newRowLabel != true) {
			break;
		}
		//4. For each row labelled with 1, label with 1 each column containing a basis cell in that row unless the column
		//is labelled with 2; then change the label of the row to 2. Go to 1.
		for (int i = 0; i < size; i++) {
			if (rowLabel[i] == 1) {
				for (auto it = basicSolnCells.begin(); it != basicSolnCells.end(); ++it) {
					if ((*it).first == i && columnLabel[(*it).second] != 2) {
						columnLabel[(*it).second] = 1;
						newColumnLabel = true;
					}
				}
				rowLabel[i] = 2;
			}
		}
	}
	//5. Stop. The set Ip(Jp) consists of the rows(columns) labelled with 2; the set Iq(Jq) consists of the 
	//unlabelled rows(columns).
	for (int i = 0; i < size; i++) {
		rowLabel[i] == 2 ? Ip.insert(i) : Iq.insert(i);
		columnLabel[i] == 2 ? Jp.insert(i) : Jq.insert(i);
	}
}


void OperatorTheory::updateDualVariables() {
	//Dual solutions for basis preserving cost operators \sigma Cpq+
	double delta = 0;
	//update row dual variables
	for (int i = 0; i < tableauSize; i++) {
		auto it = Ip.find(i);
		if (it != Ip.end()) {
			transformedRowDualVariables[i] = currentRowDualVariables[i] + delta;
		}
		else {
			transformedRowDualVariables[i] = currentRowDualVariables[i];
		}
	}
	//update column dual variables
	for (int i = 0; i < tableauSize; i++) {
		auto it = Jp.find(i);
		if (it != Jp.end()) {
			transformedColumnDualVariables[i] = currentColumnDualVariables[i] - delta;
		}
		else {
			transformedColumnDualVariables[i] = currentColumnDualVariables[i];
		}
	}
}


void OperatorTheory::findMaxDeltaAndEnteringCell() {
	double maxDelta = INFINITY;
	int cellRowID = 0;
	int cellColumnID = 0;
	int enteringCellRowID = 0;
	int enteringCellColumnID = 0;
	double val = 0;
	for (auto it: Ip) {
		for (auto itt: Jq) {
			if (it == cellRowID && itt == cellColumnID) {
				//do nothing
			}
			else {
				val = tpTableau[it][itt] - currentRowDualVariables[it] - currentColumnDualVariables[itt];
				if (val < maxDelta) {
					maxDelta = val;
					enteringCellRowID = it;
					enteringCellColumnID = itt;
				}
			}
		}
	}
}

//generates cycle or loop containing entering and leaving cells
void OperatorTheory::generateCycleAndUpdateBasicSolution() {
	std::multimap<int, int> basicSolution;
	int enteringCellRowID = 0;
	int enteringCellColumnID = 1;
	int leavingCellRowID = 3;
	int leavingCellColumnID = 1;
	std::map<int, std::map<int, double>> solutionToSupplyMap;
	std::map<int, std::map<int, double>> getterBasicCells;
	std::map<int, std::map<int, double>> giverBasicCells;
	solutionToSupplyMap[0][3] = 1;
	solutionToSupplyMap[1][0] = 1;
	solutionToSupplyMap[1][3] = 0;
	solutionToSupplyMap[2][4] = 1;
	solutionToSupplyMap[2][5] = 0;
	solutionToSupplyMap[3][1] = 1;
	solutionToSupplyMap[3][4] = 0;
	solutionToSupplyMap[4][0] = 0;
	solutionToSupplyMap[4][5] = 1;
	solutionToSupplyMap[5][2] = 1;
	solutionToSupplyMap[5][4] = 0;
	//This is how we can loop over the two dimensional map
	std::cout << "\nThe basic solution is : " << std::endl;
	for (auto it = solutionToSupplyMap.begin(); it != solutionToSupplyMap.end(); it++) {
		for (auto ptr = it->second.begin(); ptr != it->second.end(); ptr++) {
			std::cout << (*it).first << " " << (*ptr).first << " " << (*ptr).second << std::endl;
		}
	}
	//now lets form a cycle taking entering cell and leaving cell of the basic solution.
	std::cout <<"\nThe value of cell : (1, 0) is : " << solutionToSupplyMap[1][0] << std::endl;
	//solutionToSupplyMap[1][0] = 0;
	//now lets form a cycle taking entering cell and leaving cell of the basic solution.
	//std::cout << "\nThe value of cell : (1, 0) is : " << solutionToSupplyMap[1][0] << std::endl;
	//










}

//generates basic solution equivalent to the basic solution of transportation problem
void OperatorTheory::generateAcyclicConnectedGraphSolution() {
	std::vector<std::vector<double>> costTableau;
	std::multimap<int, int> assignmentSolution;
	std::multimap<int, int> transportationSolution; 
	//Need to generate acyclic connected graph which would represent the basic solution 
	//of the equivalent transportation problem.
	//populate the cost tableau
	std::vector<double> rowVec;
	rowVec.push_back(100000);
	rowVec.push_back(27); 
	rowVec.push_back(43); 
	rowVec.push_back(16); 
	rowVec.push_back(30); 
	rowVec.push_back(26);
	costTableau.push_back(rowVec);
	rowVec.clear();
	rowVec.push_back(7);
	rowVec.push_back(100000);
	rowVec.push_back(16);
	rowVec.push_back(1);
	rowVec.push_back(30);
	rowVec.push_back(25);
	costTableau.push_back(rowVec);
	rowVec.clear();
	rowVec.push_back(20);
	rowVec.push_back(13);
	rowVec.push_back(100000);
	rowVec.push_back(35);
	rowVec.push_back(5);
	rowVec.push_back(0);
	costTableau.push_back(rowVec);
	rowVec.clear();
	rowVec.push_back(21);
	rowVec.push_back(16);
	rowVec.push_back(25);
	rowVec.push_back(100000);
	rowVec.push_back(18);
	rowVec.push_back(18);
	costTableau.push_back(rowVec);
	rowVec.clear();
	rowVec.push_back(12);
	rowVec.push_back(46);
	rowVec.push_back(27);
	rowVec.push_back(48);
	rowVec.push_back(100000);
	rowVec.push_back(5);
	costTableau.push_back(rowVec);
	rowVec.clear();
	rowVec.push_back(23);
	rowVec.push_back(5);
	rowVec.push_back(5);
	rowVec.push_back(9);
	rowVec.push_back(5);
	rowVec.push_back(100000);
	costTableau.push_back(rowVec);
	rowVec.clear();
	//show the cost tableau
	for (auto &it: costTableau) {
		for (auto iter: it) {
			std::cout << " " << iter << " ";
		}
		std::cout << "\n";
	}

	//populate the assignment solution
	assignmentSolution.insert(std::pair<int, int>(0, 3));
	assignmentSolution.insert(std::pair<int, int>(1, 0));
	assignmentSolution.insert(std::pair<int, int>(2, 4));
	assignmentSolution.insert(std::pair<int, int>(3, 1));
	assignmentSolution.insert(std::pair<int, int>(4, 5));
	assignmentSolution.insert(std::pair<int, int>(5, 2));
	//for (auto & it: assignmentSolution) {
	//	std::cout << "\n" << it.first << " " << it.second << std::endl;
	//}
	//populate transportation solution with assignment solutions
	transportationSolution = assignmentSolution;
	//row scanning for the acyclic connected graph
	std::vector<int> columnCells;
	std::vector<int> rowCells;
	for (int i = 0; i < 6; i++) {
		columnCells.push_back(1);
	}
	for (int i = 0; i < 6; i++) {
		auto col = assignmentSolution.find(i);
		double val1 = INFINITY;
		double val2 = 0;
		int indx = 0;
		for (int j = 1; j < 6; j++) {
			auto it = assignmentSolution.find(i);
			j != i?	val2 = costTableau[j][(*it).second]: val2 = INFINITY;
			if (val2 < val1) {
				val1 = val2;
				indx = j;
			}
		}
		auto itLow = transportationSolution.lower_bound(indx);
		auto itUp = transportationSolution.upper_bound(indx);
		int count = 0;
		for (auto &it = itLow; it != itUp; it++) {
			count += 1;
		}
		if (count <= 1 && indx != i) {
			//exclude any assignment sol?
			transportationSolution.insert(std::pair<int,int>(indx,(*col).second));
			int numColCells = columnCells.at(indx);
			columnCells[(*col).second] = numColCells + 1;
		}
	}
	for (auto& it : transportationSolution) {
		std::cout << "\n" << it.first << " " << it.second << std::endl;
	}
	//total number of cell assigned to the basic solution
	int numOfBasicCells = 0;
	for (int i = 0; i < 6; i++) {
		numOfBasicCells += columnCells.at(i);
	}
	//populate rowcell vector
	for (int i = 0; i < 6; i++) {
		auto itLow = transportationSolution.lower_bound(i);
		auto itUp = transportationSolution.upper_bound(i);
		int val = 0;
		for (auto &it = itLow; it != itUp; it++) {
			val += 1;
		}
		rowCells.push_back(val);
	}	
	//columnMap
	std::multimap<int, int> columnMap;
	for (auto & it: transportationSolution) {
		columnMap.insert(std::pair<int, int>(it.second, it.first));
	}
	//column scanning for the acyclic connected graph
	if (numOfBasicCells < 11) {
		std::cout << "\nStart of column scanning" << std::endl;
		for (int i = 0; i < 6; i++) {
			if (columnCells.at(i) == 1) {
				//find row ID, You still has column ID
				auto it = columnMap.find(i);
				int rowID = (*it).second;
				if (rowCells.at(rowID) == 1) {
					int col = 0;
					double val3 = INFINITY;
					double val4 = 0;
					for (int k = 0; k < 6; k++) {
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
					if (numOfBasicCells == 11) {
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
	for (int i = 0; i < 6; i++) {
		std::cout << "\nRow number : " << i << " Number of contained basic cells : " << rowCells.at(i);
	}
	//print columns with number of contained basic cells
	for (int i = 0; i < 6; i++) {
		std::cout << "\nColumn number : " << i << " Number of contained basic cells : " << columnCells.at(i);
	}
	*/	
}


void OperatorTheory::updateWeakerLowerBound() {



}

void OperatorTheory::runCostOperator() {

}

