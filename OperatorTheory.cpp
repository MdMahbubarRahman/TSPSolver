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
		rowColMapForBasicSolution.insert(std::pair<int, int>(it.rowID, it.colID));
		colRowMapForBasicSolution.insert(std::pair<int, int>(it.colID, it.rowID));
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
			//std::cout << "\nLabel with 1 the columns of basis cells in row p cexept for basis cell (p,q); label row p with 2" << std::endl;
			for (auto it = basicSolution.begin(); it != basicSolution.end(); it++) {
				if ((*it).rowID == cellRowID && (*it).colID != cellColumnID) {
					columnLabel[(*it).colID] = 1;
					newColumnLabel = true;
				}
			}
			rowLabel[cellRowID] = 2;
		}
		//std::cout << "\nIf no new columns were labelled 1 on the preceding step go to 5. otherwise go to 2." << std::endl;
		if (newColumnLabel != true) {
			break;
		}
		//std::cout << "\nFor each column labelled with 1, label with 1 each row containing a basis cell in that column unless the row is labelled with 2; then change the label of the column to 2." << std::endl;
		for (int i = 0; i < size; i++) {
			if (columnLabel[i] == 1) {
				for (auto it = basicSolution.begin(); it != basicSolution.end(); ++it) {
					if ((*it).colID == i && rowLabel[(*it).rowID] != 2) {
						rowLabel[(*it).rowID] = 1;
						newRowLabel = true;
					}
				}
				columnLabel[i] = 2;
			}
		}
		newColumnLabel = false;
		//std::cout << "\nIf no new rows were labelled with 1 on the preceeding step go to 5. Otherwise go to 4." << std::endl;
		if (newRowLabel != true) {
			break;
		}
		//std::cout << "\nFor each row labelled with 1, label with 1 each column containing a basis cell in that row unless the column is labelled with 2; then change the label of the row to 2. Go to 1." << std::endl;
		for (int i = 0; i < size; i++) {
			if (rowLabel[i] == 1) {
				for (auto it = basicSolution.begin(); it != basicSolution.end(); ++it) {
					if ((*it).rowID == i && columnLabel[(*it).colID] != 2) {
						columnLabel[(*it).colID] = 1;
						newColumnLabel = true;
					}
				}
				rowLabel[i] = 2;
			}
		}
		newRowLabel = false;
	}
	//std::cout << "\nStop. The set Ip(Jp) consists of the rows(columns) labelled with 2; the set Iq(Jq) consists of the unlabelled rows(columns)." << std::endl;
	for (int i = 0; i < size; i++) {
		rowLabel[i] == 2 ? Ip.insert(i) : Iq.insert(i);
		columnLabel[i] == 2 ? Jp.insert(i) : Jq.insert(i);
	}
	//show the contents of the different sets
	std::cout << "\nIp set elements : " << std::endl;
	for (auto& it : Ip) {
		std::cout << it << std::endl;
	}
	std::cout << "\nIq set elements : " << std::endl;
	for (auto& it : Iq) {
		std::cout << it << std::endl;
	}
	std::cout << "\nJp set elements : " << std::endl;
	for (auto& it : Jp) {
		std::cout << it << std::endl;
	}
	std::cout << "\nJq set elements : " << std::endl;
	for (auto& it : Jq) {
		std::cout << it << std::endl;
	}
}

//update dual variables
void OperatorTheory::updateDualSolution(double val) {
	//Dual solutions for basis preserving cost operators \sigma Cpq+
	std::vector<double> transformedRowWiseDualSolution;
	std::vector<double> transformedColumnWiseDualSolution;
	for (int i = 0; i < costTableau.size(); i++) {
		transformedRowWiseDualSolution.push_back(0);
		transformedColumnWiseDualSolution.push_back(0);
	}
	int size = costTableau.size();
	double delta = val;
	//update row dual variables
	for (int i = 0; i < size; i++) {
		auto it = Ip.find(i);
		if (it != Ip.end()) {
			transformedRowWiseDualSolution[i] = rowWiseDualSolution.at(i) + delta;
		}
		else {
			transformedRowWiseDualSolution[i] = rowWiseDualSolution.at(i);
		}
	}
	//update column dual variables
	for (int i = 0; i < size; i++) {
		auto it = Jp.find(i);
		if (it != Jp.end()) {
			transformedColumnWiseDualSolution[i] = columnWiseDualSolution.at(i) - delta;
		}
		else {
			transformedColumnWiseDualSolution[i] = columnWiseDualSolution.at(i);
		}
	}
	std::cout << "\nShow updated dual solutions." << std::endl;
	std::cout << "\nShow row dual solution : " << std::endl;
	for (auto it : transformedRowWiseDualSolution) {
		std::cout << it << std::endl;
	}
	std::cout << "\nShow column dual solution : " << std::endl;
	for (auto it : transformedColumnWiseDualSolution) {
		std::cout << it << std::endl;
	}
}

//find max delta and find the potential entering cells
double OperatorTheory::findMaxDeltaAndEnteringCell(int p, int q) {
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
	std::cout << "\nThe entering cell row id : " << enteringCellRowID << " " << "column id : " << enteringCellColumnID << std::endl;
	return maxDelta;
}

//generates cycle or loop containing entering and leaving cells
void OperatorTheory::generateCycleAndUpdateBasicSolution(int enCellRowId, int enCellColID) {
	std::multimap<int, int> rowColMapForBasicSolution;
	std::multimap<int, int> colRowMapForBasicSolution;
	std::map<int, std::map<int, double>> cellToValueMap;
	for (auto& it : basicSolution) {
		rowColMapForBasicSolution.insert(std::pair<int, int>(it.rowID, it.colID));
		colRowMapForBasicSolution.insert(std::pair<int, int>(it.colID, it.rowID));
		cellToValueMap[it.rowID][it.colID] = it.value;
	}
	std::list<AllocatedCell> cycleOfAllocatedCells;
	int enteringCellRowID = enCellRowId;
	int enteringCellColumnID = enCellColID;
	cellToValueMap[enteringCellRowID][enteringCellColumnID] = 0.0;
	rowColMapForBasicSolution.insert(std::pair<int, int>(enteringCellRowID, enteringCellColumnID));
	colRowMapForBasicSolution.insert(std::pair<int, int>(enteringCellColumnID, enteringCellRowID));
	int leavingCellRowID = 0;
	int leavingCellColumnID = 0;
	//Entering cell cofig as allocated cell
	AllocatedCell enteringCell = AllocatedCell();
	enteringCell.cellProperty.rowID = enteringCellRowID;
	enteringCell.cellProperty.colID = enteringCellColumnID;
	enteringCell.cellProperty.value = 0.0;
	enteringCell.cellType = Getter;
	enteringCell.prevCell.rowID = enteringCellRowID;
	enteringCell.prevCell.colID = enteringCellColumnID;
	enteringCell.postCell.rowID = NULL;
	enteringCell.postCell.colID = NULL;
	AllocatedCell currentCell = enteringCell;
	bool cycleComplete = false;
	//int counter = 0;
	//std::cout << "\nRow id : " << currentCell.cellProperty.rowID << " Column id : " << currentCell.cellProperty.colID << " Cell type : " << currentCell.cellType << std::endl;
	while (!cycleComplete) {
		if ((currentCell.prevCell.rowID == currentCell.cellProperty.rowID) && (currentCell.prevCell.colID == currentCell.cellProperty.colID)) {
			//searach row 
			std::list<int> listOfColumns;
			std::list<int> listOfRows;
			for (auto it = rowColMapForBasicSolution.lower_bound(currentCell.cellProperty.rowID); it != rowColMapForBasicSolution.upper_bound(currentCell.cellProperty.rowID); it++) {
				if ((*it).second != currentCell.cellProperty.colID) {
					listOfColumns.push_back((*it).second);
				}
			}
			// search column
			for (auto it = colRowMapForBasicSolution.lower_bound(currentCell.cellProperty.colID); it != colRowMapForBasicSolution.upper_bound(currentCell.cellProperty.colID); it++) {
				if ((*it).second != currentCell.cellProperty.rowID) {
					listOfRows.push_back((*it).second);
				}
			}
			//choose post cell and update current cell
			if (listOfRows.size() <= listOfColumns.size()) {
				double val = NULL;
				for (auto it : listOfRows) {
					val = cellToValueMap[it][currentCell.cellProperty.colID];
					if (val == 1.0) {
						currentCell.postCell.rowID = it;
						currentCell.postCell.colID = currentCell.cellProperty.colID;
						break;
					}
				}
			}
			else {
				double val = NULL;
				for (auto it : listOfColumns) {
					val = cellToValueMap[currentCell.cellProperty.rowID][it];
					if (val == 1.0) {
						currentCell.postCell.rowID = currentCell.cellProperty.rowID;
						currentCell.postCell.colID = it;
						break;
					}
				}
			}
			//populate post cell as an allocated cell
			AllocatedCell newAllocatedCell = AllocatedCell();
			newAllocatedCell.cellProperty.rowID = currentCell.postCell.rowID;
			newAllocatedCell.cellProperty.colID = currentCell.postCell.colID;
			newAllocatedCell.cellProperty.value = 1.0;
			newAllocatedCell.cellType = Giver;
			newAllocatedCell.prevCell.rowID = currentCell.cellProperty.rowID;
			newAllocatedCell.prevCell.colID = currentCell.cellProperty.colID;
			newAllocatedCell.postCell.rowID = NULL;
			newAllocatedCell.postCell.colID = NULL;
			cycleOfAllocatedCells.push_back(currentCell);
			currentCell = newAllocatedCell;
			/*
			counter += 1;
			std::cout << "\nFirst counter no: " << counter << std::endl;
			std::cout << "\nRow id : " << currentCell.cellProperty.rowID << " Column id : " << currentCell.cellProperty.colID << " Cell type : " << currentCell.cellType << std::endl;
			*/
		}
		else if ((currentCell.prevCell.rowID == currentCell.cellProperty.rowID) && (currentCell.prevCell.colID != currentCell.cellProperty.colID)) {
			int rowId = NULL;
			std::vector<int> rowList;
			for (auto it = colRowMapForBasicSolution.lower_bound(currentCell.cellProperty.colID); it != colRowMapForBasicSolution.upper_bound(currentCell.cellProperty.colID); it++) {
				if (((*it).second != currentCell.cellProperty.rowID) && ((*it).second == enteringCellRowID)) {
					rowId = enteringCellRowID;
				}
			}
			if (rowId == NULL) {
				for (auto it = colRowMapForBasicSolution.lower_bound(currentCell.cellProperty.colID); it != colRowMapForBasicSolution.upper_bound(currentCell.cellProperty.colID); it++) {
					if ((*it).second != currentCell.cellProperty.rowID) {
						rowList.push_back((*it).second);
					}
				}
				if (rowList.size() == 0) {
					//needs to implement later
				}
				else if (rowList.size() == 1) {
					rowId = rowList.at(0);
				}
				else {
					if (currentCell.cellType == Giver) {
						double val = NULL;
						for (auto it : rowList) {
							val = cellToValueMap[it][currentCell.cellProperty.colID];
							if (val == 0.0) {
								rowId = it;
								break;
							}
						}

					}
					else if (currentCell.cellType == Getter) {
						double val = NULL;
						for (auto it : rowList) {
							val = cellToValueMap[it][currentCell.cellProperty.colID];
							if (val == 1.0) {
								rowId = it;
								break;
							}
						}
					}
				}
			}
			currentCell.postCell.colID = currentCell.cellProperty.colID;
			currentCell.postCell.rowID = rowId;
			if ((currentCell.postCell.rowID == enteringCellRowID) && (currentCell.postCell.colID == enteringCellColumnID)) {
				cycleComplete = true;
				cycleOfAllocatedCells.push_back(currentCell);
			}
			else {
				//populate new allocated cell
				AllocatedCell newAllocatedCell = AllocatedCell();
				newAllocatedCell.cellProperty.rowID = currentCell.postCell.rowID;
				newAllocatedCell.cellProperty.colID = currentCell.postCell.colID;
				newAllocatedCell.cellProperty.value = cellToValueMap[newAllocatedCell.cellProperty.rowID][newAllocatedCell.cellProperty.colID];
				currentCell.cellType == Getter ? newAllocatedCell.cellType = Giver : newAllocatedCell.cellType = Getter;
				newAllocatedCell.prevCell.rowID = currentCell.cellProperty.rowID;
				newAllocatedCell.prevCell.colID = currentCell.cellProperty.colID;
				newAllocatedCell.postCell.rowID = NULL;
				newAllocatedCell.postCell.colID = NULL;
				cycleOfAllocatedCells.push_back(currentCell);
				currentCell = newAllocatedCell;
			}
			/*
			counter += 1;
			std::cout << "\nMiddle counter no: " << counter << std::endl;
			std::cout << "\nRow id : " << currentCell.cellProperty.rowID << " Column id : " << currentCell.cellProperty.colID << " Cell type : " << currentCell.cellType << std::endl;
			*/
		}
		else if ((currentCell.prevCell.colID == currentCell.cellProperty.colID) && (currentCell.prevCell.rowID != currentCell.cellProperty.rowID)) {
			int colId = NULL;
			std::vector<int> colList;
			for (auto it = rowColMapForBasicSolution.lower_bound(currentCell.cellProperty.rowID); it != rowColMapForBasicSolution.upper_bound(currentCell.cellProperty.rowID); it++) {
				if (((*it).second != currentCell.cellProperty.colID) && ((*it).second == enteringCellColumnID)) {
					colId = enteringCellColumnID;
				}
			}
			if (colId == NULL) {
				for (auto it = rowColMapForBasicSolution.lower_bound(currentCell.cellProperty.rowID); it != rowColMapForBasicSolution.upper_bound(currentCell.cellProperty.rowID); it++) {
					if ((*it).second != currentCell.cellProperty.colID) {
						colList.push_back((*it).second);
					}
				}
				if (colList.size() == 0) {
					//needs to implement later
				}
				else if (colList.size() == 1) {
					colId = colList.at(0);
				}
				else {
					if (currentCell.cellType == Giver) {
						double val = NULL;
						for (auto it : colList) {
							val = cellToValueMap[currentCell.cellProperty.rowID][it];
							if (val == 0.0) {
								colId = it;
								break;
							}
						}

					}
					else if (currentCell.cellType == Getter) {
						double val = NULL;
						for (auto it : colList) {
							val = cellToValueMap[currentCell.cellProperty.rowID][it];
							if (val == 1.0) {
								colId = it;
								break;
							}
						}
					}
				}
			}
			currentCell.postCell.colID = colId;
			currentCell.postCell.rowID = currentCell.cellProperty.rowID;
			if ((currentCell.postCell.rowID == enteringCellRowID) && (currentCell.postCell.colID == enteringCellColumnID)) {
				cycleComplete = true;
				cycleOfAllocatedCells.push_back(currentCell);
			}
			else {
				//populate new allocated cell
				AllocatedCell newAllocatedCell = AllocatedCell();
				newAllocatedCell.cellProperty.rowID = currentCell.postCell.rowID;
				newAllocatedCell.cellProperty.colID = currentCell.postCell.colID;
				newAllocatedCell.cellProperty.value = cellToValueMap[newAllocatedCell.cellProperty.rowID][newAllocatedCell.cellProperty.colID];
				currentCell.cellType == Getter ? newAllocatedCell.cellType = Giver : newAllocatedCell.cellType = Getter;
				newAllocatedCell.prevCell.rowID = currentCell.cellProperty.rowID;
				newAllocatedCell.prevCell.colID = currentCell.cellProperty.colID;
				newAllocatedCell.postCell.rowID = NULL;
				newAllocatedCell.postCell.colID = NULL;
				cycleOfAllocatedCells.push_back(currentCell);
				currentCell = newAllocatedCell;
			}
			/*
			counter += 1;
			std::cout << "\nLast counter no: " << counter << std::endl;
			std::cout << "\nRow id : " << currentCell.cellProperty.rowID << " Column id : " << currentCell.cellProperty.colID << " Cell type : " << currentCell.cellType << std::endl;
			*/
		}
	}
	std::cout << "\nShow the cycle : " << std::endl;
	for (auto& it : cycleOfAllocatedCells) {
		std::cout << "\nrow id : " << it.cellProperty.rowID << " column id : " << it.cellProperty.colID << " cell type : " << it.cellType << " allocated value : " << it.cellProperty.value << std::endl;
	}
}

//generates basic solution equivalent to the basic solution of transportation problem
void OperatorTheory::generateAcyclicConnectedGraphSolution() {
	int size = costTableau.size();
	std::multimap<int, int> assignmentSolution;
	std::multimap<int, int> transportationSolution;
	std::cout << "\nNeed to generate acyclic connected graph which would represent the basic solution of the equivalent transportation problem." << std::endl;
	std::cout << "\nPopulate the assignment solution." << std::endl;
	assignmentSolution.insert(std::pair<int, int>(0, 3));
	assignmentSolution.insert(std::pair<int, int>(1, 0));
	assignmentSolution.insert(std::pair<int, int>(2, 4));
	assignmentSolution.insert(std::pair<int, int>(3, 1));
	assignmentSolution.insert(std::pair<int, int>(4, 5));
	assignmentSolution.insert(std::pair<int, int>(5, 2));
	std::cout << "\nPopulate transportation solution with assignment solutions." << std::endl;
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

//update weaker lower bound
void OperatorTheory::updateWeakerLowerBound() {


}

//run cost operator
void OperatorTheory::runCostOperator() {

}

