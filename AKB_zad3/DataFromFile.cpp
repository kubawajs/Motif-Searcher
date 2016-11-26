#include "stdafx.h"
#include "DataFromFile.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;

vector<Sequence> DataFromFile::getSeqData()
{
	return DataFromFile::seqData;
}

string DataFromFile::getDataName()
{
	return DataFromFile::dataName;
}

Matrix DataFromFile::getMatrix()
{
	return DataFromFile::matrix;
}

int DataFromFile::getReliability()
{
	return DataFromFile::reliability;
}

void DataFromFile::loadFromFile(string dataName, vector <Sequence> seqData) {
	fstream file;
	string line;
	string lastSeq = "";
	string path = dataName + ".txt";

	file.open(path);
	if (!file.good()) {
		cout << "Error! Can't open the file." << endl;
		exit(0);
	}
	else {
		Sequence s1;
		while (getline(file, line)) {
			if (line[0] == '>') {
				s1.setName(line);
			}
			else {
				if (lastSeq == "") {
					lastSeq = line;
				}
				else {
					lastSeq += line;
					s1.setSequence(lastSeq);
					seqData.push_back(s1);
					lastSeq = "";
				}
			}
		}
	}
	DataFromFile::minConnections = seqData.size() * 0.7;//seventy percents of sequences must have connections between each other
	DataFromFile::seqData = seqData;
}

void DataFromFile::loadQualFile(string dataName, vector<Sequence> seqData)
{
	fstream file;
	string line;
	vector <int> qualities;
	string number = "";
	int seqNo = 0;
	string path = dataName + "-qual.txt";

	file.open(path);
	if (!file.good()) {
		cout << "Error! Can't open the file." << endl;
		exit(0);
	}
	else {
		while (getline(file, line)) {
			if (line[0] == '>') {
				if (line != seqData[seqNo].getName()) {
					cout << "Error! Fasta file and qual file not equal." << endl;
					exit(0);
				}
				else {
					if (!qualities.empty()) {
						seqData[seqNo - 1].setQual(qualities);
						qualities.clear();
					}
				}
				seqNo++;
			}
			else {
				for (int i = 0; i < line.size(); i++) {
					if (line[i] == ' ' && number != "") {
						qualities.push_back(atoi(number.c_str()));
						number = "";
					}
					else if (line[i] != ' ' && number == "") {
						number = line[i];
					}
					else if (line[i] != ' ' && number != "") {
						number += line[i];
						qualities.push_back(atoi(number.c_str()));
						number = "";
					}
				}
			}
		}
		if (seqNo == seqData.size()) {
			seqData[seqNo - 1].setQual(qualities);
			qualities.clear();
		}
	}
	DataFromFile::seqData = seqData;
}

void DataFromFile::createGraph(Matrix matrix, vector <Sequence> data)
{
	for (int i = 0; i < data.size(); i++) {
		data[i].createSubstrings(data[i].getSequence(), data[i].getQual(), substrLength, reliability);//create all substrings of sequence
	}
	matrix.countMatrixSizeFromSeq(data, 0);//sum all sequences
	matrix.initializeMatrix(matrix.getSize());
	DataFromFile::seqData = data;
	DataFromFile::matrix = matrix;
}

void DataFromFile::filterLowSubstrs(Matrix matrix, vector<Sequence> data, int reliability)
{
	vector <Vertex> checking;
	vector <int> checkQuals, infoTable;
	int countGoodNucs, edge, limit;
	infoTable = matrix.getInfoTable();
	//loop over sequences
	for (int i = 0; i < data.size(); i++) {
		//loop over substrings in sequence
		for (int j = 0; j < data[i].getSubstrSize(); j++) {
			countGoodNucs = 0;
			checking = data[i].getSubstrings();
			checkQuals = checking[j].getQual();
			limit = (checkQuals.size() / 2) + 1; //min number of good nucs in substring
			//loop over elements of substring
			for (int k = 0; k < checkQuals.size(); k++) {
				if (checkQuals[k] >= reliability) {
					countGoodNucs++;
					if (countGoodNucs >= limit) {
						break;
					}
				}
				if (k == (limit - 1) && countGoodNucs == 0) {
					break;
				}
			}
			if (countGoodNucs < limit) {
				if (i == 0) {
					edge = j;
				}
				else {
					edge = j + infoTable[i - 1];
				}
				matrix.delVertex(edge);
			}
		}
	}
	DataFromFile::matrix = matrix;
}

void DataFromFile::setSeqData(vector<Sequence> seqData)
{
	DataFromFile::seqData = seqData;
}

void DataFromFile::createEdges(Matrix matrix, vector <Sequence> data, vector <int> infoTable, int reliability) {
	int actual = 0, createdEdges = 0;
	vector <vector <int>> graph = matrix.getMatrix();
	string waiting = "";

	cout << "Tworze polaczenia, prosze czekac...";

	for (int i = 0; i < matrix.getSize(); i++) {
		if (i == infoTable[actual]) {
			actual++;
		}
		if (graph[i][0] != -1) {
			for (int j = infoTable[actual]; j < matrix.getSize(); j++) {
				if (graph[0][j] != -1) {
					int noSeq1, noSubstr1, noSeq2, noSubstr2;
					noSeq1 = actual;//check 1st sequence
					if (actual == 0) {//check which substring of actual sequence is being analyzed
						noSubstr1 = i;
					}
					else {
						noSubstr1 = i - infoTable[actual - 1];
					}
					for (int k = 0; k < infoTable.size(); k++) {//check 2nd sequence
						if (j < infoTable[k]) {
							noSeq2 = k;
							if (k == 0) {//check which substring of actual sequence is being analyzed
								noSubstr2 = j;
							}
							else {
								noSubstr2 = j - infoTable[k - 1];
							}
							break;
						}
					}
					Vertex v1, v2;
					v1 = data[noSeq1].getSubstrById(noSubstr1);//get first vertex
					v2 = data[noSeq2].getSubstrById(noSubstr2);//get 2nd vertex

					if (data[noSeq1].compareSubstrs(v1, v2, v1.getSubstrLength(), reliability)) {
						matrix.createEdge(i, j);
						data[noSeq1].vertexLvlUp(noSubstr1);//increment level of 1st vertex
						data[noSeq2].vertexLvlUp(noSubstr2);//increment level of 2nd vertex
						createdEdges++;
					}
				}
			}
		}
	}
	graph = matrix.getMatrix();
	DataFromFile::setSeqData(data);
	DataFromFile::matrix.setMatrix(graph);
	cout << endl << "Utworzono polaczen: " << createdEdges << endl;
}

void DataFromFile::checkIfHasMinConnections(Matrix matrix)
{
	vector <int> infoTable = matrix.getInfoTable();
	vector <vector <int>> m = matrix.getMatrix();
	int actuali = 0;
	int actualj = 0;
	int noSubstr;
	bool test = false;

	for (int i = 0; i < matrix.getSize(); i++) {
		if (i == infoTable[actuali]) {
			actuali++;
		}
		int connections = 0;
		if (m[i][0] != -1) {
			actualj = 0;
			for (int j = 0; j < matrix.getSize(); j++) {
				if (j >= infoTable[actualj]) {
					actualj++;
				}
				if (m[i][j] == 1) {
					connections++;
					j = infoTable[actualj];
				}
			}
			if (connections >= DataFromFile::minConnections) {
				if (actuali == 0) {//check which substring of actual sequence is being analyzed
					noSubstr = i;
				}
				else {
					noSubstr = i - infoTable[actuali - 1];
				}
				DataFromFile::seqData[actuali].setVertexHasMinConnections(noSubstr);
			}
		}
	}
}

void DataFromFile::printSequences()
{
	for (int i = 0; i < DataFromFile::seqData.size(); i++) { //for each sequence
		int seqSize = seqData[i].getSequence().size(); //get length of sequence
		string sequence(seqSize, '-'); //create string with '-' chars with length of sequence
		vector <Vertex> substrs = DataFromFile::seqData[i].getSubstrings();

		cout << "\n" << DataFromFile::seqData[i].getName() << endl; //print name of sequence
		cout << "O: " << DataFromFile::seqData[i].getSequence() << endl;

		for (int j = 0; j < seqData[i].getSubstrSize(); j++) { //for each substring in sequence
			int x = j;
			if (substrs[j].getHasMinConnections()) {
				for (int k = 0; k < substrs[j].getSubstrLength(); k++) { //for each char in substring
					if (substrs[j].getQual()[k] >= DataFromFile::reliability) {
						sequence[x + k] = substrs[j].getSubstring()[k];
					}
					else {
						sequence[x + k] = '*';
					}
				}
			}
		}

		cout << "D: " << sequence << "\n" << endl;
	}
}

void DataFromFile::createListOfVerticesSorted()
{
	vector <vector <int>> graph = DataFromFile::matrix.getMatrix();
	vector <int> infoTable = DataFromFile::getInfoTable(DataFromFile::getMatrix());
	int seqId, noSubstr;
	Vertex v1;

	//create list
	for (int i = 0; i < graph.size(); i++) {
		if (graph[i][0] != -1) {
			seqId = DataFromFile::matrix.getSequenceIdFromMatrix(i); //get sequence number
			
			if (seqId == 0) {//check which substring of actual sequence is being analyzed
				noSubstr = i;
			}
			else {
				noSubstr = i - infoTable[seqId - 1];
			}
			
			v1 = DataFromFile::seqData[seqId].getSubstrById(noSubstr);
			
			if (v1.getHasMinConnections()) {
				VertexInList vertexInList(v1, i, seqId);
				DataFromFile::vertexByLevel.push_back(vertexInList);
			}	
		}
	}

	//sort list
	DataFromFile::sortByVertexLvl(DataFromFile::vertexByLevel, 0, DataFromFile::vertexByLevel.size()-1);
}

void DataFromFile::buildMaxClique() {
	//TODO: build max clique
	//dopóki przyrost wielkoœci kliki > 4 + 1 pêtla dodatkowo
	//buduj klike dla danego zbioru wierzcholkow na podstawie listy wierzcholkow posortowanej wg vertex level
	//		> dany zbior wierzcholkow => 1. iteracja normalnie na vertexByLevel, 2. i kolejne z budowanego zbioru wierzcholkow do sprawdzenia
	//			> zbior wierzcholkow do sprawdzenia - dla kazdego wierzcholka w rozwiazaniu +- 1, 2 indeksy, posortuj wedlug poziomu wierzcholkow
	//dodaj klike do aktualnego rozwiazania
	vector <VertexInList> result;
	vector <VertexInList> vertexToCheck;
	vector <VertexInList> temporaryResult;
	vector <VertexInList> verticesToAdd;

	int increase = 0;
	bool increaseUnderMin = false;

	result = DataFromFile::buildClique(DataFromFile::vertexByLevel);

	while (!increaseUnderMin) {
		int lastIncrease = increase;

		vertexToCheck.clear();
		vertexToCheck = DataFromFile::prepareVertexSet(result, 2);// !!! important - sensitivity
		temporaryResult = DataFromFile::buildClique(vertexToCheck);
		
		verticesToAdd.clear();
		//check if not in actual result
		for (int i = 0; i < vertexToCheck.size(); i++) {
			for (int j = 0; j < result.size(); j++) {
				if (vertexToCheck[i].getIndex() == result[j].getIndex()) {
					break;
				}
				else if (j == result.size() - 1) {
					verticesToAdd.push_back(vertexToCheck[i]);
				}
			}
		}
		
		//add to result
		if (!verticesToAdd.empty()) {
			for (int i = 0; i < verticesToAdd.size(); i++) {
				result.push_back(verticesToAdd[i]);
			}
		}

		//TODO: change sensitivity of building result
		//check if increase is more than 4
		increase = result.size() - lastIncrease;

		cout << "Result series of cliques size: " << result.size() << endl;
	}

	cout << "Result status: Ready" << endl;
}

vector <VertexInList> DataFromFile::prepareVertexSet(vector <VertexInList> actualResult, int sensitivity) {
	//TODO: zmien t¹ funkcjê ¿eby nie dublowa³a dodawanych wierzcho³ków
	//tzn. jak sprawdza dla wierzcho³ka i dodawanie i +- j to dla i+1 niech ju¿ nie sprawdza/nie dodaje wierzcho³ka!!!

	vector <VertexInList> vertexToCheck;
	vector <int> infoTable = DataFromFile::matrix.getInfoTable();

	for (int i = 0; i < actualResult.size(); i++) {
		int index = actualResult[i].getIndex();
		int noSeq = actualResult[i].getSeqIndex();
		for (int j = -sensitivity; j <= sensitivity; j++) {
			if (j != 0) {
				if (actualResult[i].getIndex() + j >= 0 && actualResult[i].getIndex() + j < infoTable[infoTable.size() - 1]) {
					if (noSeq - 1 >= 0) { //check if is in the range of infoTable
						if (index + j < infoTable[noSeq] && index + j > infoTable[noSeq - 1]) {
							int noSubstr = index - infoTable[noSeq - 1] + j;
							Vertex v1 = DataFromFile::seqData[noSeq].getSubstrById(noSubstr);
							VertexInList vertexInList(v1, index + j, noSeq);
							vertexToCheck.push_back(vertexInList);
						}
					}
					else {
						if (index + j < infoTable[noSeq] && index + j >= 0) {
							int noSubstr = index + j;
							Vertex v1 = DataFromFile::seqData[noSeq].getSubstrById(noSubstr);
							VertexInList vertexInList(v1, index + j, noSeq);
							vertexToCheck.push_back(vertexInList);
						}
					}
				}
			}
		}
	}
	if (vertexToCheck.size() > 0) {
		DataFromFile::sortByVertexLvl(vertexToCheck, 0, vertexToCheck.size() - 1);
	}
	cout << "Vertex set size: " << vertexToCheck.size() << " Ready" << endl;
	return vertexToCheck;
}

vector <VertexInList> DataFromFile::buildClique(vector<VertexInList> vertexByLevel) {
	vector <VertexInList> clique;

	for (int i = 0; i < vertexByLevel.size(); i++) {
		VertexInList analyzedVertex = vertexByLevel[i];
		if (clique.size() == 0) {
			clique.push_back(analyzedVertex);
		}
		else if (checkConnectionsInClique(clique, analyzedVertex, DataFromFile::getMatrix())) {
			clique.push_back(analyzedVertex);
		}
	}

	cout << "Clique built." << endl;//TODO: del this line
	return clique;
}

bool DataFromFile::checkConnectionsInClique(vector <VertexInList> result, VertexInList analyzedVertex, Matrix matrix) {
	vector <vector <int>> graph = matrix.getMatrix();
	for (int i = 0; i < result.size(); i++) {
		if (graph[analyzedVertex.getIndex()][result[i].getIndex()] < 1) {//TODO: przemyslec - && analyzedVertex.getSeqIndex() != result[i].getSeqIndex()) { 
			return false;
		}
	}
	return true;
}

vector <int> DataFromFile::getInfoTable(Matrix matrix) {
	return matrix.getInfoTable();
}

void DataFromFile::sortByVertexLvl(vector <VertexInList> &vertexInLvlList, int left, int right) {
	int i = left;
	int j = right;
	int x = vertexInLvlList[(left + right) / 2].getVertex().getVertexLvl();
	do {
		while (vertexInLvlList[i].getVertex().getVertexLvl() > x)
			i++;
		while (vertexInLvlList[j].getVertex().getVertexLvl() < x)
			j--;
		if (i <= j) {
			swap(vertexInLvlList[i], vertexInLvlList[j]);
			i++;
			j--;
		}
	} while (i <= j);

	if (left < j) sortByVertexLvl(vertexInLvlList, left, j);
	if (right > i) sortByVertexLvl(vertexInLvlList, i, right);
}

DataFromFile::DataFromFile()
{
}

DataFromFile::DataFromFile(string dataName, int substrLength, int reliability)
{
	DataFromFile::dataName = dataName;
	DataFromFile::substrLength = substrLength;
	DataFromFile::reliability = reliability;
}

DataFromFile::~DataFromFile()
{
}