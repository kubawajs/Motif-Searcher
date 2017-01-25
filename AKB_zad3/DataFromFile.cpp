#include "stdafx.h"
#include "DataFromFile.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <map>

using namespace std;

vector<Sequence> DataFromFile::getSeqData() const
{
	return seqData;
}

string DataFromFile::getDataName() const
{
	return dataName;
}

Matrix DataFromFile::getMatrix() const
{
	return matrix;
}

int DataFromFile::getReliability() const
{
	return reliability;
}

void DataFromFile::loadFromFile(string dataName, vector<Sequence> seqData)
{
	fstream file;
	string line;
	string lastSeq = "";
	string path = dataName + ".txt";

	file.open(path);
	if (!file.good())
	{
		cout << "Error! Can't open the file." << endl;
		exit(0);
	}
	else
	{
		Sequence s1;
		int seqId = 0;
		while (getline(file, line))
		{
			if (line[0] == '>')
			{
				if (lastSeq != "")
				{
					s1.setSequence(lastSeq);
					seqData.push_back(s1);
					lastSeq = "";
				}
				s1.setName(line);
				s1.setSeqId(seqId);
				seqId++;
			}
			else
			{
				if (lastSeq == "")
				{
					lastSeq = line;
				}
				else
				{
					lastSeq += line;
				}
			}
		}
		s1.setSequence(lastSeq);
		seqData.push_back(s1);
		lastSeq = "";
	}
	minConnections = seqData.size() * 0.55;//55 percents of sequences must have connections between each other
	DataFromFile::seqData = seqData;
}

void DataFromFile::loadQualFile(string dataName, vector<Sequence> seqData)
{
	fstream file;
	string line;
	vector<int> qualities;
	string number = "";
	int seqNo = 0;
	string path = dataName + "-qual.txt";

	file.open(path);
	if (!file.good())
	{
		cout << "Error! Can't open the file." << endl;
		exit(0);
	}
	else
	{
		while (getline(file, line))
		{
			if (line[0] == '>')
			{
				if (line != seqData[seqNo].getName())
				{
					cout << "Error! Fasta file and qual file not equal." << endl;
					exit(0);
				}
				else
				{
					if (!qualities.empty())
					{
						seqData[seqNo - 1].setQual(qualities);
						qualities.clear();
					}
				}
				seqNo++;
			}
			else
			{
				for (int i = 0; i < line.size(); i++)
				{
					if (line[i] == ' ' && number != "")
					{
						qualities.push_back(atoi(number.c_str()));
						number = "";
					}
					else if (line[i] != ' ' && number == "")
					{
						number = line[i];
					}
					else if (line[i] != ' ' && number != "")
					{
						number += line[i];
						qualities.push_back(atoi(number.c_str()));
						number = "";
					}
				}
			}
		}
		if (seqNo == seqData.size())
		{
			seqData[seqNo - 1].setQual(qualities);
			qualities.clear();
		}
	}
	DataFromFile::seqData = seqData;
}

void DataFromFile::createGraph(Matrix matrix, vector<Sequence> data)
{
	for (int i = 0; i < data.size(); i++)
	{
		data[i].createSubstrings(data[i].getSequence(), data[i].getQual(), substrLength, reliability);//create all substrings of sequence
	}
	matrix.countMatrixSizeFromSeq(data, 0);//sum all sequences
	matrix.initializeMatrix(matrix.getSize());
	seqData = data;
	DataFromFile::matrix = matrix;
}

void DataFromFile::filterLowSubstrs(Matrix matrix, vector<Sequence> data, int reliability)
{
	vector<Vertex> checking;
	vector<int> checkQuals, infoTable;
	int countGoodNucs, edge, limit;
	infoTable = matrix.getInfoTable();
	//loop over sequences
	for (int i = 0; i < data.size(); i++)
	{
		//loop over substrings in sequence
		for (int j = 0; j < data[i].getSubstrSize(); j++)
		{
			countGoodNucs = 0;
			checking = data[i].getSubstrings();
			checkQuals = checking[j].getQual();
			limit = checking[j].getQual().size() - PERMITTED_DELETIONS; //min number of good nucs in substring
			//loop over elements of substring
			for (int k = 0; k < checkQuals.size(); k++)
			{
				if (checkQuals[k] >= reliability)
				{
					countGoodNucs++;
					if (countGoodNucs >= limit)
					{
						break;
					}
				}
				if (k == (limit - 1) && countGoodNucs == 0)
				{
					break;
				}
			}
			if (countGoodNucs < limit)
			{
				if (i == 0)
				{
					edge = j;
				}
				else
				{
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
	if (!seqData.empty())
	{
		DataFromFile::seqData = seqData;
	}
}

void DataFromFile::createEdges(Matrix matrix, vector<Sequence> data, vector<int> infoTable, int reliability)
{
	int actual = 0, createdEdges = 0;
	auto graph = matrix.getMatrix();
	string waiting = "";

	cout << "Tworze polaczenia, prosze czekac...";

	for (int i = 0; i < matrix.getSize(); i++)
	{
		if (i == infoTable[actual])
		{
			actual++;
		}
		if (graph[i][0] != -1)
		{
			for (int j = infoTable[actual]; j < matrix.getSize(); j++)
			{
				if (graph[0][j] != -1)
				{
					int noSeq1, noSubstr1, noSeq2 = -1, noSubstr2 = -1;
					noSeq1 = actual;//check 1st sequence
					if (actual == 0)
					{//check which substring of actual sequence is being analyzed
						noSubstr1 = i;
					}
					else
					{
						noSubstr1 = i - infoTable[actual - 1];
					}
					for (int k = 0; k < infoTable.size(); k++)
					{//check 2nd sequence
						if (j < infoTable[k])
						{
							noSeq2 = k;
							if (k == 0)
							{//check which substring of actual sequence is being analyzed
								noSubstr2 = j;
							}
							else
							{
								noSubstr2 = j - infoTable[k - 1];
							}
							break;
						}
					}

					if (noSeq2 != -1)
					{
						auto v1 = data[noSeq1].getSubstrById(noSubstr1);//get first vertex
						auto v2 = data[noSeq2].getSubstrById(noSubstr2);//get 2nd vertex

						if (data[noSeq1].compareSubstrs(v1, v2, v1.getSubstrLength(), reliability))
						{
							matrix.createEdge(i, j);
							data[noSeq1].vertexLvlUp(noSubstr1);//increment level of 1st vertex
							data[noSeq2].vertexLvlUp(noSubstr2);//increment level of 2nd vertex
							createdEdges++;
						}
					}
				}
			}
		}
	}
	graph = matrix.getMatrix();
	setSeqData(data);
	DataFromFile::matrix.setMatrix(graph);
	cout << endl << "Utworzono polaczen: " << createdEdges << endl;
}

void DataFromFile::checkIfHasMinConnections(Matrix matrix)
{
	auto infoTable = matrix.getInfoTable();
	auto m = matrix.getMatrix();
	int actuali = 0, actualj, noSubstr;

	for (int i = 0; i < matrix.getSize(); i++)
	{
		if (i == infoTable[actuali])
		{
			actuali++;
		}

		int connections = 0, linkedSeq = 0;
		vector<bool> sequencesToMark;
		for (int x = 0; x < getSeqData().size(); x++)
		{
			sequencesToMark.push_back(false);
		}

		if (m[i][0] != -1)
		{
			actualj = 0;
			for (int j = 0; j < matrix.getSize(); j++)
			{
				if (j >= infoTable[actualj])
				{
					actualj++;
				}
				if (m[i][j] == 1)
				{
					connections++;
					sequencesToMark[actualj] = true;
					j = infoTable[actualj];
				}
			}
			if (connections >= minConnections)
			{
				if (actuali == 0)
				{//check which substring of actual sequence is being analyzed
					noSubstr = i;
				}
				else
				{
					noSubstr = i - infoTable[actuali - 1];
				}
				for (int k = 0; k < sequencesToMark.size(); k++)
				{
					if (sequencesToMark[k])
					{
						linkedSeq++;
					}
				}
				seqData[actuali].setVertexNumOfConSeq(noSubstr, linkedSeq);
				seqData[actuali].setVertexHasMinConnections(noSubstr);
			}
		}
	}
}

void DataFromFile::printSequences(vector<Sequence> seqData, vector<Sequence> resultSeq)
{
	ResultMotif resultMotif;
	for (int i = 0; i < seqData.size(); i++)
	{ //for each sequence
		int seqSize = seqData[i].getSequence().size(); //get length of sequence
		string sequence(seqSize, '-'); //create string with '-' chars with length of sequence
		auto substrs = seqData[i].getSubstrings();

		for (int j = 0; j < seqData[i].getSubstrSize(); j++)
		{ //for each substring in sequence
			if (substrs[j].getHasMinConnections())
			{
				for (int k = 0; k < substrs[j].getSubstrLength(); k++)
				{ //for each char in substring
					if (substrs[j].getQual()[k] >= reliability)
					{
						sequence[j + k] = substrs[j].getSubstring()[k];
					}
					else
					{
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
				auto resultVertices = resultSeq[j].getSubstrings();
				resultMotif.printMotifOnSeq(resultVertices, seqSize, reliability);
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
	auto graph = matrix.getMatrix();
	auto infoTable = getInfoTable(getMatrix());
	int seqId, noSubstr;
	Vertex v1;

	//create list
	for (int i = 0; i < graph.size(); i++)
	{
		if (graph[i][0] != -1)
		{
			seqId = matrix.getSequenceIdFromMatrix(i); //get sequence number

			if (seqId == 0)
			{//check which substring of actual sequence is being analyzed
				noSubstr = i;
			}
			else
			{
				noSubstr = i - infoTable[seqId - 1];
			}

			v1 = seqData[seqId].getSubstrById(noSubstr);

			if (v1.getHasMinConnections())
			{
				v1.setIndex(i);
				v1.setSeqIndex(seqId);
				vertexByLevel.push_back(v1);
			}
		}
	}

	//sort list
	sort(vertexByLevel.begin(), vertexByLevel.end());
}

void DataFromFile::buildResults()
{
	vector<Vertex> startingClique;
	ResultMotif resultMotif;
	//jesli results puste
	startingClique = buildClique(vertexByLevel);
	addResult(buildResult(startingClique));

	//jesli ju¿ jest jakieœ rozwi¹zanie usun pierwotna klike z wektora startowego

	int i = 1;
	while (i < NUMBER_OF_RESULTS)
	{
		vertexByLevel = filterVector(vertexByLevel, startingClique);
		if (!vertexByLevel.empty())
		{
			startingClique = buildClique(vertexByLevel);
			resultMotif = buildResult(startingClique);
			if (std::find(results.begin(), results.end(), resultMotif) == results.end())
			{
				addResult(resultMotif);
				i++;
			}
		}
		else
		{
			i = NUMBER_OF_RESULTS;
		}
	}

	sortResults(results);
}

void DataFromFile::addResult(ResultMotif result)
{
	results.push_back(result);
}

ResultMotif DataFromFile::buildResult(vector<Vertex> startingClique)
{
	vector<Vertex> result, vertexToCheckLeft, vertexToCheckRight, temporaryResult, verticesToAdd;
	string motif;

	result = startingClique;
	motif = buildMotif(result, reliability);

	bool isInBuild = true;
	//extend on left
	do
	{
		vertexToCheckLeft = prepareVertexSetLeft(result, SENSITIVITY);
		temporaryResult = buildClique(vertexToCheckLeft);
		if (temporaryResult.size() < 0.55 * seqData.size())
		{
			break;
		}
		string tempMotif = buildMotif(temporaryResult, reliability);
		motif = parseMotifLeft(motif, tempMotif);
		result = sumResult(result, temporaryResult);
	}
	while (isInBuild);

	//extend on right
	isInBuild = true;
	do
	{
		vertexToCheckRight = prepareVertexSetRight(result, SENSITIVITY);
		temporaryResult = buildClique(vertexToCheckRight);
		if (temporaryResult.size() < 0.55 * seqData.size())
		{
			break;
		}
		string tempMotif = buildMotif(temporaryResult, reliability);
		motif = parseMotifRight(motif, tempMotif);
		result = sumResult(result, temporaryResult);
	}
	while (isInBuild);

	//Creating results
	ResultMotif readyResult;
	readyResult.setResult(result, seqData.size());
	readyResult.setMotif(motif);
	readyResult.parseSequences(reliability);
	readyResult.setStartingClique(startingClique);

	return readyResult;
}

vector<ResultMotif> DataFromFile::getResults() const
{
	return results;
}

string DataFromFile::buildMotif(vector<Vertex> verticesToAlign, int reliability)
{
	int actSeq, it = 0, maxValue = 0, seqId;
	Vertex v1;
	map<string, int> submotifs = {};
	vector<char> preMotif;
	string consensusMotif;

	//std::sort(verticesToAlign.begin(), verticesToAlign.end(), [](const auto& A, const auto& B) {
	//return A.getIndex() < B.getIndex(); });

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
		if (toCompare.find(motifToAdd) != string::npos)
		{
			existingMotif = toAdd + existingMotif;
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
		if (toCompare.find(motifToAdd) != string::npos)
		{
			existingMotif = existingMotif + toAdd;
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
			printSequences(seqData, results[i].getSequences());
		}
		else
		{
			break;
		}
	}
}

vector<Vertex> DataFromFile::prepareVertexSetLeft(vector<Vertex> actualResult, int sensitivity)
{
	auto resultSortedByIndex = actualResult;
	auto infoTable = getInfoTable(matrix);
	vector<Vertex> vertexSet;

	std::sort(resultSortedByIndex.begin(), resultSortedByIndex.end(), [](const auto& A, const auto& B)
	          {
		          return A.getIndex() < B.getIndex();
	          });

	int seqId, actSeq;
	vector<Vertex> verticesToAdd;
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
					v1 = seqData[seqId].getSubstrById(noSubstr);

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
					v1 = seqData[seqId].getSubstrById(noSubstr);

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

	/*if (!vertexSet.empty()) {
		std::sort(vertexSet.begin(), vertexSet.end());
	}*/

	return vertexSet;
}

vector<Vertex> DataFromFile::prepareVertexSetRight(vector<Vertex> actualResult, int sensitivity)
{//TODO: sprawdz czy zamiast sensitivity lepiej nie daæ szerokoœæ okna

	auto resultSortedByIndex = actualResult;
	auto infoTable = getInfoTable(matrix);
	vector<Vertex> vertexSet;

	std::sort(resultSortedByIndex.begin(), resultSortedByIndex.end(), [](const auto& A, const auto& B)
	          {
		          return A.getIndex() < B.getIndex();
	          });
	reverse(resultSortedByIndex.begin(), resultSortedByIndex.end());

	int seqId, actSeq;
	vector<Vertex> verticesToAdd;
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
			if (seqId < seqData.size() - 1)
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

					v1 = seqData[seqId].getSubstrById(noSubstr);

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
				if (index < matrix.getMatrix()[0].size())
				{
					noSubstr = index - infoTable[seqId - 1];
					v1 = seqData[seqId].getSubstrById(noSubstr);

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

	if (!vertexSet.empty())
	{
		sort(vertexSet.begin(), vertexSet.end());
	}

	return vertexSet;
}

vector<Vertex> DataFromFile::buildClique(vector<Vertex> vertexByLevel) const
{
	vector<Vertex> clique;
	vector<bool> usedSequences;

	if (vertexByLevel.empty())
	{
		return clique;
	}
	sort(vertexByLevel.begin(), vertexByLevel.end());

	for (int i = 0; i < seqData.size(); i++)
	{
		usedSequences.push_back(false);
	}

	for (int i = 0; i < vertexByLevel.size(); i++)
	{
		Vertex analyzedVertex = vertexByLevel[i];
		if (clique.size() == 0)
		{
			clique.push_back(analyzedVertex);
			usedSequences[analyzedVertex.getSeqIndex()] = true;
		}
		else if (!usedSequences[analyzedVertex.getSeqIndex()] && checkConnectionsInClique(clique, analyzedVertex, DataFromFile::getMatrix()))
		{
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
		bool found = (find(actualResult.begin(), actualResult.end(), tempResult[i]) != actualResult.end());

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
	for (auto& elementToFilter : toFilter) // access by reference to avoid copying
	{
		if (find(filtering.begin(), filtering.end(), elementToFilter) == filtering.end())
		{
			result.push_back(elementToFilter);
		}
	}

	return result;
}

bool DataFromFile::checkConnectionsInClique(vector<Vertex> result, Vertex analyzedVertex, Matrix matrix)
{
	auto graph = matrix.getMatrix();

	for (int i = 0; i < result.size(); i++)
	{
		if (graph[analyzedVertex.getIndex()][result[i].getIndex()] < 1)
		{
			return false;
		}
	}
	return true;
}

vector<int> DataFromFile::getInfoTable(Matrix matrix)
{
	return matrix.getInfoTable();
}

void DataFromFile::sortResults(vector<ResultMotif>& results) const
{
	sort(results.begin(), results.end());
	reverse(results.begin(), results.end());
}

void DataFromFile::printBestMotifs(vector<ResultMotif> results)
{
	int i = 1;
	if (!results.empty() && results.begin()->getMotif().size() != 0)
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
	substrLength = 4;
	reliability = 10;
	minConnections = 1;
}

DataFromFile::DataFromFile(string _dataName, int _substrLength, int _reliability)
{
	dataName = _dataName;
	substrLength = _substrLength;
	reliability = _reliability;
	minConnections = 0;
}

DataFromFile::~DataFromFile()
{
}
