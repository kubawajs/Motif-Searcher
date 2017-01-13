#include "stdafx.h"
#include "DataFromFile.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <map>

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

void DataFromFile::printSequences(vector <Sequence> seqData)
{
	for (int i = 0; i < seqData.size(); i++) { //for each sequence
		int seqSize = seqData[i].getSequence().size(); //get length of sequence
		string sequence(seqSize, '-'); //create string with '-' chars with length of sequence
		vector <Vertex> substrs = seqData[i].getSubstrings();

		cout << "\n" << seqData[i].getName() << endl; //print name of sequence
		cout << "O: " << DataFromFile::seqData[i].getSequence() << endl;

		for (int j = 0; j < seqData[i].getSubstrSize(); j++) { //for each substring in sequence
			if (substrs[j].getHasMinConnections()) {
				for (int k = 0; k < substrs[j].getSubstrLength(); k++)
				{ //for each char in substring
					if (substrs[j].getQual()[k] >= DataFromFile::reliability) {
						sequence[j + k] = substrs[j].getSubstring()[k];
					}
					else {
						sequence[j + k] = '*';
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
				v1.setIndex(i);
				v1.setSeqIndex(seqId);
				DataFromFile::vertexByLevel.push_back(v1);
			}
		}
	}

	//sort list
	DataFromFile::sortByVertexLvl(DataFromFile::vertexByLevel, 0, DataFromFile::vertexByLevel.size() - 1);
}

void DataFromFile::buildMaxClique() {
	vector <Vertex> result;
	vector <Vertex> vertexToCheckLeft, vertexToCheckRight;
	vector <Vertex> temporaryResult;
	vector <Vertex> verticesToAdd;
	string motif;

	result = DataFromFile::buildClique(DataFromFile::vertexByLevel);
	motif = DataFromFile::buildMotif(result, DataFromFile::reliability);

	//rozbuduj w lewo
	do
	{
		//TODO: zapetla sie na ostatnim podciagu jesli jest spelniony warunek (moze dac cos zeby sprawdzalo czy juz nalezy do rozwiazania?
		//dodawaj do result, albo na input do szukania kliki dawaj output z poprzedniego
		vertexToCheckLeft = DataFromFile::prepareVertexSetLeft(result, SENSITIVITY);//TODO: moze lepiej sprawdzi sie przesuniecie o max dlugosc podciagu - 1?
		temporaryResult = DataFromFile::buildClique(vertexToCheckLeft);
		string tempMotif = DataFromFile::buildMotif(temporaryResult, DataFromFile::reliability);
		motif = DataFromFile::parseMotifLeft(motif, tempMotif);
		//add tempRes to result but check if not exist yet
	} while (temporaryResult.size() >= (0.6 * DataFromFile::seqData.size()));//TODO: przemysl czy nie zamienic na 0,55 -> da to min 4/7 sekwencji a nie 5/7

	//rozbuduj w prawo
	do
	{
		vertexToCheckRight = DataFromFile::prepareVertexSetRight(result, SENSITIVITY);
		temporaryResult = DataFromFile::buildClique(vertexToCheckRight);
		string tempMotif = DataFromFile::buildMotif(temporaryResult, DataFromFile::reliability);
		motif = DataFromFile::parseMotifRight(motif, tempMotif);
		//add tempRes to result
	} while (temporaryResult.size() >= (0.6 * DataFromFile::seqData.size()));

	cout << "Result status: Ready" << endl;

	//Creating results
	Result readyResult(result, DataFromFile::seqData.size());
	readyResult.parseSequences(DataFromFile::reliability);

	cout << "Printed";
}

string DataFromFile::buildMotif(vector <Vertex> verticesToAlign, int reliability)//build motif for clique
{
	int actSeq, it = 0, maxValue = 0;
	int seqId = verticesToAlign[0].getSeqIndex();
	Vertex v1;
	map <string, int> submotifs = {};
	vector <char> preMotif;
	string consensusMotif;

	while (it < verticesToAlign.size())
	{
		v1 = verticesToAlign[it];
		seqId = v1.getSeqIndex();
		actSeq = seqId;
		preMotif.clear();


		for (int i = 0; i < v1.getQual().size(); i++)
		{
			if (v1.getQual()[i] >= reliability)
			{
				preMotif.push_back(v1.getSubstring()[i]);
			}
		}
		string motif(preMotif.begin(), preMotif.end());
		submotifs[motif]++;//add to map

		while (actSeq == seqId && it < verticesToAlign.size())
		{
			it++;
			if (it < verticesToAlign.size())
			{
				actSeq = verticesToAlign[it].getSeqIndex();
			}
		}
	}

	for (auto const& motif : submotifs)//foreach motif in motifs
	{
		if (motif.second > maxValue)
		{
			consensusMotif = motif.first;
			maxValue = motif.second;
		}
		else if (motif.second == maxValue && motif.first.length() > consensusMotif.length())
		{
			consensusMotif = motif.first;
			maxValue = motif.second;
		}
	}

	return consensusMotif; //temporary
}

string DataFromFile::parseMotifLeft(string existingMotif, string motifToAdd)
{
	string toAdd = "";

	for (int i = 0; i <= SENSITIVITY; i++)
	{
		if (i != 0) {
			toAdd += motifToAdd[0];
		}

		motifToAdd.erase(0, i);
		if (existingMotif.find(motifToAdd) != string::npos) {
			existingMotif = toAdd + existingMotif;
			i = SENSITIVITY;
			break;
		}
	}

	//cout << existingMotif << endl;
	return existingMotif;
}

string DataFromFile::parseMotifRight(string existingMotif, string motifToAdd)
{
	string toAdd = "";

	for (int i = 0; i <= SENSITIVITY; i++)
	{
		if (i != 0) {
			toAdd += motifToAdd[motifToAdd.size() - 1];
		}

		motifToAdd = motifToAdd.substr(0, motifToAdd.size() - 1);
		if (existingMotif.find(motifToAdd) != string::npos) {
			existingMotif = toAdd + existingMotif;
			i = SENSITIVITY;
			break;
		}
	}

	//cout << existingMotif << endl;
	return existingMotif;
}

vector <Vertex> DataFromFile::prepareVertexSetLeft(vector <Vertex> actualResult, int sensitivity) {

	vector <Vertex> resultSortedByIndex = actualResult;
	vector <Vertex> vertexSet;
	vector <int> infoTable = DataFromFile::getInfoTable(DataFromFile::matrix);

	DataFromFile::sortByIndex(resultSortedByIndex, 0, actualResult.size() - 1);

	int seqId, actInd, actSeq;
	vector <Vertex> verticesToAdd;
	Vertex v1;

	int it = 0;

	while (it < resultSortedByIndex.size())
	{
		seqId = resultSortedByIndex[it].getSeqIndex();
		actSeq = seqId;
		v1 = resultSortedByIndex[it];
		verticesToAdd.push_back(v1);

		while (actSeq == seqId && it < resultSortedByIndex.size())
		{
			it++;
			if (it < resultSortedByIndex.size())
			{
				actSeq = resultSortedByIndex[it].getSeqIndex();
			}
		}
	}

	for (int i = 0; i < verticesToAdd.size(); i++)
	{
		//TODO: why index -1
		int index = verticesToAdd[i].getIndex();
		int noSubstr, j = 0;
		seqId = verticesToAdd[i].getSeqIndex();

		while (j < sensitivity)
		{
			index -= 1;
			if (seqId > 0)
			{
				if (index >= infoTable[seqId - 1])
				{
					noSubstr = index - infoTable[seqId - 1];
					v1 = DataFromFile::seqData[seqId].getSubstrById(noSubstr);

					if (v1.getHasMinConnections())
					{
						v1.setIndex(index);
						v1.setSeqIndex(seqId);
						vertexSet.push_back(v1);
						j = sensitivity;
					}
				}
			}
			else
			{
				if (index >= 0)
				{
					noSubstr = index;
					v1 = DataFromFile::seqData[seqId].getSubstrById(noSubstr);

					if (v1.getHasMinConnections())
					{
						v1.setIndex(index);
						v1.setSeqIndex(seqId);
						vertexSet.push_back(v1);
						j = sensitivity;
					}
				}
			}
			j++;
		}
	}

	if (!vertexSet.empty()) {
		DataFromFile::sortByVertexLvl(vertexSet, 0, vertexSet.size() - 1);
	}

	return vertexSet;
}

vector <Vertex> DataFromFile::prepareVertexSetRight(vector <Vertex> actualResult, int sensitivity) {//TODO: sprawdz czy zamiast sensitivity lepiej nie daæ szerokoœæ okna

	vector <Vertex> resultSortedByIndex = actualResult;
	vector <Vertex> vertexSet;
	vector <int> infoTable = DataFromFile::getInfoTable(DataFromFile::matrix);

	DataFromFile::sortByIndex(resultSortedByIndex, 0, actualResult.size() - 1);

	int seqId, actInd, actSeq;
	vector <Vertex> verticesToAdd;
	Vertex v1;

	int it = 0;

	while (it < resultSortedByIndex.size())
	{
		seqId = resultSortedByIndex[it].getSeqIndex();
		actSeq = seqId;
		v1 = resultSortedByIndex[it];
		verticesToAdd.push_back(v1);

		while (actSeq == seqId && it < resultSortedByIndex.size())
		{
			it++;
			if (it < resultSortedByIndex.size())
			{
				actSeq = resultSortedByIndex[it].getSeqIndex();
			}
		}
	}

	for (int i = 0; i < verticesToAdd.size(); i++)
	{
		int index = verticesToAdd[i].getIndex();
		int noSubstr, j = 0;
		seqId = verticesToAdd[i].getSeqIndex();

		while (j < sensitivity)
		{
			index += 1;
			if (seqId < DataFromFile::seqData.size() - 1)
			{
				if (index < infoTable[seqId])
				{
					if (seqId == 0)
					{
						noSubstr = index;
					}
					else
					{
						noSubstr = index - infoTable[seqId - 1];
					}

					v1 = DataFromFile::seqData[seqId].getSubstrById(noSubstr);

					if (v1.getHasMinConnections())
					{
						v1.setIndex(index);
						v1.setSeqIndex(seqId);
						vertexSet.push_back(v1);
						j = sensitivity;
					}
				}
			}
			else
			{
				if (index < DataFromFile::matrix.getMatrix()[0].size())
				{
					noSubstr = index - infoTable[seqId - 1];
					v1 = DataFromFile::seqData[seqId].getSubstrById(noSubstr);

					if (v1.getHasMinConnections())
					{
						v1.setIndex(index);
						v1.setSeqIndex(seqId);
						vertexSet.push_back(v1);
						j = sensitivity;
					}
				}
			}
			j++;
		}
	}

	if (!vertexSet.empty()) {
		DataFromFile::sortByVertexLvl(vertexSet, 0, vertexSet.size() - 1);
	}

	return vertexSet;
}

vector <Vertex> DataFromFile::buildClique(vector<Vertex> vertexByLevel) {
	vector <Vertex> clique;

	for (int i = 0; i < vertexByLevel.size(); i++) {
		Vertex analyzedVertex = vertexByLevel[i];
		if (clique.size() == 0) {
			clique.push_back(analyzedVertex);
		}
		else if (checkConnectionsInClique(clique, analyzedVertex, DataFromFile::getMatrix())) {
			clique.push_back(analyzedVertex);
		}
	}

	return clique;
}

bool DataFromFile::checkConnectionsInClique(vector <Vertex> result, Vertex analyzedVertex, Matrix matrix) {
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

void DataFromFile::sortByVertexLvl(vector <Vertex> &vertexInLvlList, int left, int right) {
	int i = left;
	int j = right;
	int x = vertexInLvlList[(left + right) / 2].getVertexLvl();
	do {
		while (vertexInLvlList[i].getVertexLvl() > x)
			i++;
		while (vertexInLvlList[j].getVertexLvl() < x)
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

void DataFromFile::sortByIndex(vector<Vertex> &vertexInLvlList, int left, int right)
{
	int i = left;
	int j = right;
	int x = vertexInLvlList[(left + right) / 2].getIndex();
	do {
		while (vertexInLvlList[i].getIndex() < x)
			i++;
		while (vertexInLvlList[j].getIndex() > x)
			j--;
		if (i <= j) {
			swap(vertexInLvlList[i], vertexInLvlList[j]);
			i++;
			j--;
		}
	} while (i <= j);

	if (left < j) sortByIndex(vertexInLvlList, left, j);
	if (right > i) sortByIndex(vertexInLvlList, i, right);
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