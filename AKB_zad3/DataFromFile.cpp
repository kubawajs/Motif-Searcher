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
		int seqId = 0;
		while (getline(file, line)) {
			if (line[0] == '>') {
				s1.setName(line);
				s1.setSeqId(seqId);
				seqId++;
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
			limit = checking[j].getQual().size() - PERMITTED_DELETIONS; //min number of good nucs in substring
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

		int connections = 0, linkedSeq = 0;
		vector<bool> sequencesToMark;
		for (int i = 0; i < DataFromFile::getSeqData().size(); i++)
		{
			sequencesToMark.push_back(false);
		}

		if (m[i][0] != -1) {
			actualj = 0;
			for (int j = 0; j < matrix.getSize(); j++) {
				if (j >= infoTable[actualj]) {
					actualj++;
				}
				if (m[i][j] == 1) {
					connections++;
					sequencesToMark[actualj] = true;
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
				for (int k = 0; k < sequencesToMark.size(); k++)
				{
					if (sequencesToMark[k])
					{
						linkedSeq++;
					}
				}
				DataFromFile::seqData[actuali].setVertexNumOfConSeq(noSubstr, linkedSeq);
				DataFromFile::seqData[actuali].setVertexHasMinConnections(noSubstr);
			}
		}
	}
}

void DataFromFile::printSequences(vector <Sequence> seqData, vector <Sequence> resultSeq)
{
	ResultMotif resultMotif;
	for (int i = 0; i < seqData.size(); i++) { //for each sequence
		int seqSize = seqData[i].getSequence().size(); //get length of sequence
		string sequence(seqSize, '-'); //create string with '-' chars with length of sequence
		vector <Vertex> substrs = seqData[i].getSubstrings();

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

		cout << "\n" << seqData[i].getName() << endl; //print name of sequence
		cout << "O: " << DataFromFile::seqData[i].getSequence() << endl;
		cout << "D: " << sequence << endl;

		//print motifs
		for (int j = 0; j < resultSeq.size(); j++)
		{
			if (seqData[i].getSeqId() == resultSeq[j].getSeqId())
			{
				vector <Vertex> resultVertices = resultSeq[j].getSubstrings();
				resultMotif.printMotifOnSeq(resultVertices, seqSize, reliability);
				j = resultSeq.size();
				break;
			}
			else if (j == resultSeq.size() - 1)
			{
				string seqMotifEmpty(seqSize, '-');
				cout << "M: " << seqMotifEmpty << '\n' << endl;
			}
		}
	}

	//print vertices
	cout << "Wierzcholki biorace udzial w budowie rozwiazania (id_sekwencji:id_podciagu): " << endl;
	for (int i = 0; i < resultSeq.size(); i++)
	{
		resultMotif.printVerticesInMotif(resultSeq[i].getSubstrings(), reliability);
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
	std::sort(vertexByLevel.begin(), vertexByLevel.end());
}

void DataFromFile::buildResults()
{
	vector <Vertex> startingClique;
	ResultMotif resultMotif;
	//jesli results puste
	startingClique = DataFromFile::buildClique(vertexByLevel);
	DataFromFile::addResult(DataFromFile::buildResult(startingClique));

	//jesli ju¿ jest jakieœ rozwi¹zanie usun pierwotna klike z wektora startowego

	int i = 1;
	while (i < NUMBER_OF_RESULTS)
	{
		vertexByLevel = DataFromFile::filterVector(vertexByLevel, startingClique);
		if (!vertexByLevel.empty())
		{
			startingClique = DataFromFile::buildClique(vertexByLevel);
			resultMotif = DataFromFile::buildResult(startingClique);
			if (std::find(DataFromFile::results.begin(), DataFromFile::results.end(), resultMotif) == DataFromFile::results.end())
			{
				DataFromFile::addResult(resultMotif);
				i++;
			}
		}
		else
		{
			i = NUMBER_OF_RESULTS;
		}
	}

	DataFromFile::sortResults(DataFromFile::results);
}

void DataFromFile::addResult(ResultMotif result)
{
	DataFromFile::results.push_back(result);
}

ResultMotif DataFromFile::buildResult(vector <Vertex> startingClique) {
	vector <Vertex> result, vertexToCheckLeft, vertexToCheckRight, temporaryResult, verticesToAdd;
	string motif;

	result = startingClique;
	motif = DataFromFile::buildMotif(result, DataFromFile::reliability);

	bool isInBuild = true;
	//extend on left
	do
	{
		vertexToCheckLeft = DataFromFile::prepareVertexSetLeft(result, SENSITIVITY);
		temporaryResult = DataFromFile::buildClique(vertexToCheckLeft);
		if (temporaryResult.size() < 0.55 * DataFromFile::seqData.size())
		{
			isInBuild = false;
			break;
		}
		string tempMotif = DataFromFile::buildMotif(temporaryResult, DataFromFile::reliability);
		motif = DataFromFile::parseMotifLeft(motif, tempMotif);
		result = DataFromFile::sumResult(result, temporaryResult);
	} while (isInBuild);

	//extend on right
	isInBuild = true;
	do
	{
		vertexToCheckRight = DataFromFile::prepareVertexSetRight(result, SENSITIVITY);
		temporaryResult = DataFromFile::buildClique(vertexToCheckRight);
		if (temporaryResult.size() < 0.55 * DataFromFile::seqData.size())
		{
			isInBuild = false;
			break;
		}
		string tempMotif = DataFromFile::buildMotif(temporaryResult, DataFromFile::reliability);
		motif = DataFromFile::parseMotifRight(motif, tempMotif);
		result = DataFromFile::sumResult(result, temporaryResult);
	} while (isInBuild);

	//Creating results
	ResultMotif readyResult;
	readyResult.setResult(result, DataFromFile::seqData.size());
	readyResult.setMotif(motif);
	readyResult.parseSequences(DataFromFile::reliability);
	readyResult.setStartingClique(startingClique);

	return readyResult;
}

vector<ResultMotif> DataFromFile::getResults()
{
	return DataFromFile::results;
}

string DataFromFile::buildMotif(vector <Vertex> verticesToAlign, int reliability)//build motif for clique
{
	int actSeq, it = 0, maxValue = 0;
	int seqId = verticesToAlign[0].getSeqIndex();
	Vertex v1;
	map <string, int> submotifs = {};
	vector <char> preMotif;
	string consensusMotif;
	DataFromFile::sortByIndex(verticesToAlign, 0, verticesToAlign.size() - 1);

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
	string toCompare = "";

	if (existingMotif.size() >= motifToAdd.size())
	{
		toCompare = existingMotif.substr(0, motifToAdd.size());
	}

	for (int i = 1; i <= SENSITIVITY; i++)
	{
		if (toCompare.find(motifToAdd) != string::npos) {
			existingMotif = toAdd + existingMotif;
			i = SENSITIVITY;
			break;
		}

		toCompare = toCompare.substr(0, toCompare.size() - 1);
		toAdd += motifToAdd[0];
		motifToAdd.erase(0, i);
	}

	return existingMotif;
}

string DataFromFile::parseMotifRight(string existingMotif, string motifToAdd)
{
	string toAdd = "";
	string toCompare = "";

	if (existingMotif.size() >= motifToAdd.size())
	{
		int lengthDiff = existingMotif.size() - motifToAdd.size();
		toCompare = existingMotif.substr(lengthDiff, motifToAdd.size());
	}

	for (int i = 1; i <= SENSITIVITY; i++)
	{
		if (toCompare.find(motifToAdd) != string::npos) {
			existingMotif = existingMotif + toAdd;
			i = SENSITIVITY;
			break;
		}
		else if (toCompare != "")
		{
			toCompare = toCompare.substr(1, toCompare.size());
			toAdd += motifToAdd[motifToAdd.size() - 1];
			motifToAdd.pop_back();
		}
	}

	return existingMotif;
}

void DataFromFile::printResult(vector<ResultMotif> result)
{
	int longestMotif = 0;
	for (int i = 0; i < result.size(); i++)
	{
		if (i == 0)
		{
			longestMotif = result[0].getMotif().size();
		}

		if (longestMotif == result[i].getMotif().size())
		{
			cout << "Znalezione rozwiazanie id: " << i << endl;
			DataFromFile::printSequences(DataFromFile::seqData, DataFromFile::results[i].getSequences());
		}
		else
		{
			i = result.size();
			break;
		}
	}
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
	std:sort(vertexSet.begin(), vertexSet.end());
	}

	return vertexSet;
}

vector <Vertex> DataFromFile::prepareVertexSetRight(vector <Vertex> actualResult, int sensitivity) {//TODO: sprawdz czy zamiast sensitivity lepiej nie daæ szerokoœæ okna

	vector <Vertex> resultSortedByIndex = actualResult;
	vector <Vertex> vertexSet;
	vector <int> infoTable = DataFromFile::getInfoTable(DataFromFile::matrix);


	DataFromFile::sortByIndex(resultSortedByIndex, 0, actualResult.size() - 1);
	reverse(resultSortedByIndex.begin(), resultSortedByIndex.end());

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
	std:sort(vertexSet.begin(), vertexSet.end());
	}

	return vertexSet;
}

vector <Vertex> DataFromFile::buildClique(vector<Vertex> vertexByLevel) {
	vector <Vertex> clique;
	vector <bool> usedSequences;

	if (vertexByLevel.empty())
	{
		return clique;
	}
std:sort(vertexByLevel.begin(), vertexByLevel.end());

	for (int i = 0; i < DataFromFile::seqData.size(); i++)
	{
		usedSequences.push_back(false);
	}

	for (int i = 0; i < vertexByLevel.size(); i++) {
		Vertex analyzedVertex = vertexByLevel[i];
		if (clique.size() == 0) {
			clique.push_back(analyzedVertex);
			usedSequences[analyzedVertex.getSeqIndex()] = true;
		}
		else if (!usedSequences[analyzedVertex.getSeqIndex()] && checkConnectionsInClique(clique, analyzedVertex, DataFromFile::getMatrix())) {
			clique.push_back(analyzedVertex);
			usedSequences[analyzedVertex.getSeqIndex()] = true;
		}
	}

	return clique;
}

vector<Vertex> DataFromFile::sumResult(vector<Vertex> actualResult, vector<Vertex> tempResult)
{
	for (int i = 0; i < tempResult.size(); i++)
	{
		bool found = (std::find(actualResult.begin(), actualResult.end(), tempResult[i]) != actualResult.end());

		if (!found)
		{
			actualResult.push_back(tempResult[i]);
		}
	}

	return actualResult;
}

vector<Vertex> DataFromFile::filterVector(vector<Vertex> toFilter, vector<Vertex> filtering)
{
	vector<Vertex> result;
	for (auto &elementToFilter : toFilter) // access by reference to avoid copying
	{
		if (std::find(filtering.begin(), filtering.end(), elementToFilter) == filtering.end()) {
			result.push_back(elementToFilter);
		}
	}

	return result;
}

bool DataFromFile::checkConnectionsInClique(vector <Vertex> result, Vertex analyzedVertex, Matrix matrix) {
	vector <vector <int>> graph = matrix.getMatrix();
	for (int i = 0; i < result.size(); i++)
	{
		if (graph[analyzedVertex.getIndex()][result[i].getIndex()] < 1)
		{
			return false;
		}
	}
	return true;
}

vector <int> DataFromFile::getInfoTable(Matrix matrix) {
	return matrix.getInfoTable();
}

void DataFromFile::sortByIndex(vector<Vertex> &vertexInLvlList, int left, int right)
{
	int i = left;
	int j = right;
	int x = vertexInLvlList[(left + right) / 2].getIndex();
	do
	{
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

void DataFromFile::sortResults(vector<ResultMotif>& results)
{
	std::sort(results.begin(), results.end());
	std::reverse(results.begin(), results.end());
}

void DataFromFile::printBestMotifs(vector <ResultMotif> results)
{
	int i = 1;
	if (!results.empty())
	{
		cout << endl << "Znaleziono motyw:" << endl;
		cout << results[0].getMotif() << endl << endl;
		while (i < results.size())
		{
			if (results[i - 1].getMotif().size() == results[i].getMotif().size())
			{
				cout << endl << "Znaleziono motyw:" << endl;
				cout << results[i].getMotif() << endl << endl;
				i++;
			}
			else
			{
				i = results.size();
				break;
			}
		}
	}
	else
	{
		cout << endl << "Nie znaleziono motywu w podanych sekwencjach";
	}
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