#include "OperatorTheory.h"

//default constructor
OperatorTheory::OperatorTheory() {
	std::cout << "\nThe default constructor has been called!" << std::endl;
	nodeIndex = 0;
	parentNodeIndex = 0;
	isSolutionATour = false;
	weakerLowerBound = -INFINITY;
	lowerBound = -INFINITY;
	delta = 0.0;
}

OperatorTheory::OperatorTheory(Node node) {
	nodeIndex = node.nodeIndex;
	parentNodeIndex = node.parentNodeIndex;
	isSolutionATour = node.isSolutionATour;
	weakerLowerBound = node.weakerLowerBound;
	lowerBound = node.lowerBound;
	delta = node.delta;
	Ip = node.Ip;
	Iq = node.Iq;
	Jp = node.Jp;
	Jq = node.Jq;
	branchOnCell = node.branchOnCell;
	basicSolution = node.basicSolution;
	rowWiseDualSolution = node.rowDualSolution;
	columnWiseDualSolution = node.columnDualSolution;
	costTableau = node.costTableau;
}

//default operator to perform cost operator from tsp solution
OperatorTheory::OperatorTheory(std::list<BasicCell> basicSol, std::vector<std::vector<double>> cTableau) {
	//populate all other values using these two input values
	nodeIndex = 0;
	parentNodeIndex = 0;
	isSolutionATour = false;
	weakerLowerBound = -INFINITY;
	lowerBound = -INFINITY;
	delta = 0.0;
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
	nodeIndex = 0;
	parentNodeIndex = 0;
	isSolutionATour = false;
	weakerLowerBound = -INFINITY;
	lowerBound = -INFINITY;
	delta = 0.0;
	basicSolution = basicSol;
	rowWiseDualSolution = rowWiseDualSol;
	columnWiseDualSolution = columnWiseDualSol;
	costTableau = cTableau;
}

//copy constructor
OperatorTheory::OperatorTheory(const OperatorTheory& opThr) {
	std::cout << "\nThe copy constructor has been called!" << std::endl;
	nodeIndex = opThr.nodeIndex;
	parentNodeIndex = opThr.parentNodeIndex;
	isSolutionATour = opThr.isSolutionATour;
	weakerLowerBound = opThr.weakerLowerBound;
	lowerBound = opThr.lowerBound;
	delta = opThr.delta;
	Ip = opThr.Ip;
	Iq = opThr.Iq;
	Jp = opThr.Jp;
	Jq = opThr.Jq;
	branchOnCell = opThr.branchOnCell;
	basicSolution = opThr.basicSolution;
	rowWiseDualSolution = opThr.rowWiseDualSolution;
	columnWiseDualSolution = opThr.columnWiseDualSolution;
	costTableau = opThr.costTableau;
	childNodes = opThr.childNodes;
	cStatus = opThr.cStatus;
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
	std::cout << "\nThe entering cell row id : " << enteringCellRowID << " " << "column id : " << enteringCellColumnID << std::endl;
	delta = maxDelta;
	Cell cell = Cell();
	cell.rowID = enteringCellRowID;
	cell.colID = enteringCellColumnID;
	enteringCell = cell;
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
	cycleOfCells = cycleOfAllocatedCells;
}

void OperatorTheory::generateListOfRoutes(std::map<int, int> boxPoints) {
	std::vector<int> route;
	int sourceNode = 0;
	int destinationNode = 0;
	int routeStartNode = 0;
	std::map<int, int> newAssignments;
	newAssignments = boxPoints;
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
				clistOfRoutes.push_back(route);
				route.clear();
				i = 0;
			}
			else {
				i += 1;
			}
		}
	}
}

//update weaker lower bound
void OperatorTheory::runCostOperatorForGeneratingRootNodes() {
	generateInitialDualSolution();
	std::map<int, int> boxPoints;
	for (auto& it : basicSolution) {
		if (it.value == 1.0) {
			boxPoints.insert(std::pair<int, int>(it.rowID, it.colID));
		}
	}
	generateListOfRoutes(boxPoints);
	updateBounds(boxPoints);
	//generate nodes to start branch and bound search
	std::vector<int> route;
	if (clistOfRoutes.size() > 1) {
		for (auto& it : clistOfRoutes) {
			route = it;
			break;
		}
		isSolutionATour = false;
		std::map<int, int> cellMaps;
		for (int i = 0; i < route.size() - 1; i++) {
			cellMaps.insert(std::pair<int, int>(route.at(i), route.at(i + 1)));
		}
		//create nodes
		double nodeNum = parentNodeIndex;
		for (auto& it : cellMaps) {
			nodeNum += 1;
			Cell cell = Cell();
			cell.rowID = it.first;
			cell.colID = it.second;
			branchOnCell = cell;
			scanningRoutine(branchOnCell.rowID, branchOnCell.colID);
			findMaxDeltaAndEnteringCell(branchOnCell.rowID, branchOnCell.colID);
			double weakerLB = 0;
			double lb = 0;
			weakerLB = lowerBound + delta;
			lb = lowerBound;
			Node node = Node();
			node.nodeIndex = nodeNum;
			node.parentNodeIndex = parentNodeIndex;
			node.isSolutionATour = false;
			node.weakerLowerBound = weakerLB;
			node.lowerBound = lb;
			node.delta = delta;
			node.Ip = Ip;
			node.Iq = Iq;
			node.Jp = Jp;
			node.Jq = Jq;
			node.branchOnCell = branchOnCell;
			node.enteringCell = enteringCell;
			node.basicSolution = basicSolution;
			node.rowDualSolution = rowWiseDualSolution;
			node.columnDualSolution = columnWiseDualSolution;
			node.costTableau = costTableau;
			childNodes.push_back(node);
		}
	}
	else if (clistOfRoutes.size() == 1) {
		Node node = Node();
		node.nodeIndex = nodeIndex;
		node.parentNodeIndex = parentNodeIndex;
		node.isSolutionATour = true;
		node.weakerLowerBound = weakerLowerBound;
		node.lowerBound = lowerBound;
		node.delta = 0;
		node.basicSolution = basicSolution;
		node.rowDualSolution = rowWiseDualSolution;
		node.columnDualSolution = columnWiseDualSolution;
		node.costTableau = costTableau;
		childNodes.push_back(node);
	}
}

//update bounds
void OperatorTheory::updateBounds(std::map<int, int> bPoints) {
	double cost = 0;
	std::map<int, int> boxPoints = bPoints;
	for (auto& it : boxPoints) {
		cost += costTableau[it.first][it.second];
	}
	weakerLowerBound = cost;
	lowerBound = cost;
}

//run cost operator
void OperatorTheory::runCostOperatorForSolvingANode() {
	generateCycleAndUpdateBasicSolution(enteringCell.rowID, enteringCell.colID);
	//need works
}

std::list<Node> OperatorTheory::getChildNodes() {
	return childNodes;
}

void OperatorTheory::showChildNodes() {
	std::cout << "\nShow the child nodes : " << std::endl;
	for (auto& it : childNodes) {
		std::cout << "Node index : " << it.nodeIndex << std::endl;
		std::cout << "Parent node index : " << it.parentNodeIndex << std::endl;
		std::cout << "Is solution a tour : " << it.isSolutionATour << std::endl;
		std::cout << "Weaker lower bound : " << it.weakerLowerBound << std::endl;
		std::cout << "Lower bound : " << it.lowerBound << std::endl;
		std::cout << "Delta : " << it.delta << std::endl;
		std::cout << "Ip : " << std::endl;
		for (auto ip : it.Ip) {
			std::cout << ip << " ";
		}
		std::cout << "\nIq : " << std::endl;
		for (auto iq : it.Iq) {
			std::cout << iq << " ";
		}
		std::cout << "\nJp : " << std::endl;
		for (auto jp : it.Jp) {
			std::cout << jp << " ";
		}
		std::cout << "\nJq : " << std::endl;
		for (auto jq : it.Jq) {
			std::cout << jq << " ";
		}
		std::cout << "\nBranch on cell : " << "Row id : " << it.branchOnCell.rowID << ", Column id : " << it.branchOnCell.colID << std::endl;
		std::cout << "Entering cell : " << "Row id : " << it.enteringCell.rowID << ", Column id : " << it.enteringCell.colID << std::endl;
		std::cout << "Show basic solution : " << std::endl;
		for (auto& bSol : it.basicSolution) {
			std::cout << "Row id : " << bSol.rowID << " Column id : " << bSol.colID << " Value : " << bSol.value << std::endl;
		}
		std::cout << "Show rowwise dual solution : " << std::endl;
		for (auto& rd : it.rowDualSolution) {
			std::cout << rd << " ";
		}
		std::cout << "\nShow column wise dual solution : " << std::endl;
		for (auto& cd : it.columnDualSolution) {
			std::cout << cd << " ";
		}
		std::cout << "\nShow the cost tableau : " << std::endl;
		for (int i = 0; i < it.costTableau.size(); i++) {
			for (int j = 0; j < it.costTableau.size(); j++) {
				std::cout << it.costTableau[i][j] << " ";
			}
			std::cout << " " << std::endl;
		}
	}
}



