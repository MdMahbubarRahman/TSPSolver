#include "OperatorTheory.h"

//default constructor
OperatorTheory::OperatorTheory() {
	std::cout << "\nThe default constructor has been called!" << std::endl;
}

//default operator to perform cost operator from tsp solution
OperatorTheory::OperatorTheory(std::list<BasicCell> basicSol, std::vector<std::vector<double>> cTableau) {
	//populate all other values using these two input values
	basicSolution = basicSol;
	costTableau = cTableau;
	for (int i = 0; i < costTableau.size(); i++) {
		rowWiseDualSolution.push_back(0);
		columnWiseDualSolution.push_back(0);
	}
}

//constructor
OperatorTheory::OperatorTheory(std::list<BasicCell> basicSol, std::vector<double> rowWiseDualSol, std::vector<double> columnWiseDualSol, std::vector<std::vector<double>> cTableau) {
	std::cout << "\nThe operator theory constructor has been called!" << std::endl;
	basicSolution = basicSol;
	rowWiseDualSolution = rowWiseDualSol;
	columnWiseDualSolution = columnWiseDualSol;
	costTableau = cTableau;
}

//copy constructor
OperatorTheory::OperatorTheory(const OperatorTheory& opThr) {
	std::cout << "\nThe copy constructor has been called!" << std::endl;
	basicSolution = opThr.basicSolution;
	rowWiseDualSolution = opThr.rowWiseDualSolution;
	columnWiseDualSolution = opThr.columnWiseDualSolution;
	costTableau = opThr.costTableau;
}

//generate initial dual solution
void OperatorTheory::generateInitialDualSolution() {
	std::multimap<int, int> rowColMapForBasicSolution;
	std::multimap<int, int> colRowMapForBasicSolution;
	for (auto& it : basicSolution) {
		rowColMapForBasicSolution.insert(std::pair<int, int>(it.rowId, it.colID));
		colRowMapForBasicSolution.insert(std::pair<int, int>(it.colID, it.rowId));
	}
	//generate dual solution
	bool dualflag = false;
	std::map<int, int> mapForRowDualSolution;
	std::map<int, int> mapForColumnDualSolution;
	for (int i = 0; i < costTableau.size(); i++) {
		mapForRowDualSolution.insert(std::pair<int, int>(i, 0));
		mapForColumnDualSolution.insert(std::pair<int, int>(i, 0));
	}
	mapForRowDualSolution.erase(0);
	mapForRowDualSolution.insert(std::pair<int, int>(0, 1));
	rowWiseDualSolution[0] = 0;
	while (!dualflag) {
		//update column dual solution
		for (auto& it : mapForRowDualSolution) {
			if (it.second == 1) {
				for (auto itt = rowColMapForBasicSolution.lower_bound(it.first); itt != rowColMapForBasicSolution.upper_bound(it.first); itt++) {
					int col = (*itt).second;
					if (mapForColumnDualSolution[col] == 0) {
						mapForColumnDualSolution.erase(col);
						mapForColumnDualSolution.insert(std::pair<int, int>(col, 1));
						double rowDualVal = rowWiseDualSolution.at(it.first);
						double cellCost = costTableau[(*itt).first][col];
						double colDualVal = cellCost - rowDualVal;
						columnWiseDualSolution[col] = colDualVal;
					}
				}
			}
		}
		//update row dual solution
		for (auto& it : mapForColumnDualSolution) {
			if (it.second == 1) {
				for (auto itt = colRowMapForBasicSolution.lower_bound(it.first); itt != colRowMapForBasicSolution.upper_bound(it.first); itt++) {
					int row = (*itt).second;
					if (mapForRowDualSolution[row] == 0) {
						mapForRowDualSolution.erase(row);
						mapForRowDualSolution.insert(std::pair<int, int>(row, 1));
						double colDualVal = columnWiseDualSolution.at(it.first);
						double cellCost = costTableau[row][(*itt).first];
						double rowDualVal = cellCost - colDualVal;
						rowWiseDualSolution[row] = rowDualVal;
					}
				}
			}
		}
		//check whether dual solution generation is complete
		int counter = 0;
		for (auto& it : mapForColumnDualSolution) {
			if (it.second == 1) {
				counter += 1;
			}
		}
		for (auto& it : mapForRowDualSolution) {
			if (it.second == 1) {
				counter += 1;
			}
		}
		std::cout << "\nThe number of dual solutions found : " << counter << std::endl;
		if (counter == (2 * costTableau.size())) {
			dualflag = true;
		}
	}
	std::cout << "\nShow the dual solutions." << std::endl;
	std::cout << "\nThe row dual solution are : " << std::endl;
	for (auto& it : rowWiseDualSolution) {
		std::cout << it << std::endl;
	}
	std::cout << "\nThe column dual solution are : " << std::endl;
	for (auto& it : columnWiseDualSolution) {
		std::cout << it << std::endl;
	}
}

//scanning routine to populate Ip, Iq, Jp, Jq sets
void OperatorTheory::scanningRoutine(int p, int q) {
	std::vector<int> rowLabel;
	std::vector<int> columnLabel;
	int cellRowID = p;
	int cellColumnID = q;
	int size = costTableau.size();
	bool stopFlag = false;
	int iter = 0;
	//Populate rowlabel and columnlabel with  default values;
	for (int i = 0; i < size; i++) {
		rowLabel.push_back(0);
		columnLabel.push_back(0);
	}
	//Run scanning routine
	bool newColumnLabel = false;
	bool newRowLabel = false;
	while (!stopFlag) {
		iter += 1;
		if (iter == 1) {
			std::cout << "\nLabel with 1 the columns of basis cells in row p cexept for basis cell (p,q); label row p with 2" << std::endl;
			for (auto it = basicSolution.begin(); it != basicSolution.end(); it++) {
				if ((*it).rowId == cellRowID && (*it).colID != cellColumnID) {
					columnLabel[(*it).colID] = 1;
					newColumnLabel = true;
				}
			}
			rowLabel[cellRowID] = 2;
		}
		std::cout << "\nIf no new columns were labelled 1 on the preceding step go to 5. otherwise go to 2." << std::endl;
		if (newColumnLabel != true) {
			break;
		}
		std::cout << "\nFor each column labelled with 1, label with 1 each row containing a basis cell in that column unless the row is labelled with 2; then change the label of the column to 2." << std::endl;
		for (int i = 0; i < size; i++) {
			if (columnLabel[i] == 1) {
				for (auto it = basicSolution.begin(); it != basicSolution.end(); ++it) {
					if ((*it).colID == i && rowLabel[(*it).rowId] != 2) {
						rowLabel[(*it).rowId] = 1;
						newRowLabel = true;
					}
				}
				columnLabel[i] = 2;
			}
		}
		newColumnLabel = false;
		std::cout << "\nIf no new rows were labelled with 1 on the preceeding step go to 5. Otherwise go to 4." << std::endl;
		if (newRowLabel != true) {
			break;
		}
		std::cout << "\nFor each row labelled with 1, label with 1 each column containing a basis cell in that row unless the column is labelled with 2; then change the label of the row to 2. Go to 1." << std::endl;
		for (int i = 0; i < size; i++) {
			if (rowLabel[i] == 1) {
				for (auto it = basicSolution.begin(); it != basicSolution.end(); ++it) {
					if ((*it).rowId == i && columnLabel[(*it).colID] != 2) {
						columnLabel[(*it).colID] = 1;
						newColumnLabel = true;
					}
				}
				rowLabel[i] = 2;
			}
		}
		newRowLabel = false;
	}
	std::cout << "\nStop. The set Ip(Jp) consists of the rows(columns) labelled with 2; the set Iq(Jq) consists of the unlabelled rows(columns)." << std::endl;
	for (int i = 0; i < size; i++) {
		rowLabel[i] == 2 ? Ip.insert(i) : Iq.insert(i);
		columnLabel[i] == 2 ? Jp.insert(i) : Jq.insert(i);
	}
}

//update dual variables
void OperatorTheory::updateDualSolution() {
	//Dual solutions for basis preserving cost operators \sigma Cpq+
	std::vector<double> transformedRowWiseDualSolution;
	std::vector<double> transformedColumnWiseDualSolution;
	int size = costTableau.size();
	double delta = 0;
	//update row dual variables
	for (int i = 0; i < size; i++) {
		auto it = Ip.find(i);
		if (it != Ip.end()) {
			transformedRowWiseDualSolution[i] = rowWiseDualSolution[i] + delta;
		}
		else {
			transformedRowWiseDualSolution[i] = rowWiseDualSolution[i];
		}
	}
	//update column dual variables
	for (int i = 0; i < size; i++) {
		auto it = Jp.find(i);
		if (it != Jp.end()) {
			transformedColumnWiseDualSolution[i] = columnWiseDualSolution[i] - delta;
		}
		else {
			transformedColumnWiseDualSolution[i] = columnWiseDualSolution[i];
		}
	}
}

//find max delta and find the potential entering cells
void OperatorTheory::findMaxDeltaAndEnteringCell(int p, int q) {
	double maxDelta = INFINITY;
	int cellRowID = p;
	int cellColumnID = q;
	int enteringCellRowID = 0;
	int enteringCellColumnID = 0;
	double val = 0;
	for (auto it : Ip) {
		for (auto itt : Jq) {
			if (it == cellRowID && itt == cellColumnID) {
				//do nothing
			}
			else {
				val = costTableau[it][itt] - rowWiseDualSolution[it] - columnWiseDualSolution[itt];
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
	std::map<BasicCell, double> cellToAllocatedValueMap;
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
	std::cout << "\nThe value of cell : (1, 0) is : " << solutionToSupplyMap[1][0] << std::endl;
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
	for (auto& it : costTableau) {
		for (auto iter : it) {
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
	for (int i = 0; i < 6; i++) {
		numOfBasicCells += columnCells.at(i);
	}
	//populate rowcell vector
	for (int i = 0; i < 6; i++) {
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

//update weaker lower bound
void OperatorTheory::updateWeakerLowerBound() {


}

//run cost operator
void OperatorTheory::runCostOperator() {

}

