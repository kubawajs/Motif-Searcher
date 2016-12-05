#pragma once

#include "Sequence.h"
#include "Matrix.h"
#include "VertexInList.h"

using namespace std;

class DataFromFile
{
	vector <Sequence> seqData;
	vector <VertexInList> vertexByLevel;
	Matrix matrix;
	string dataName;
	int substrLength;
	int reliability;
	int minConnections;

public:
	vector <Sequence> getSeqData();
	string getDataName();
	Matrix getMatrix();
	int getReliability();
	void loadFromFile(string dataName, vector <Sequence> seqData);
	void loadQualFile(string dataName, vector <Sequence> seqData);
	void createGraph(Matrix matrix, vector <Sequence> data);
	void filterLowSubstrs(Matrix matrix, vector<Sequence> data, int reliability);
	//filter substrings with lower qual and notice in matrix
	void setSeqData(vector <Sequence> seqData);
	void createEdges(Matrix matrix, vector<Sequence> data, vector<int> infoTable, int reliability);
	void checkIfHasMinConnections(Matrix matrix);//function for checking if vertex has connection with min 5 other sequences
	void printSequences();
	void createListOfVerticesSorted();
	void buildMaxClique();
	vector<VertexInList> prepareVertexSet(vector<VertexInList> actualResult, int sensitivity);
	vector<VertexInList> buildClique(vector <VertexInList> vertexByLevel);
	vector<int> getInfoTable(Matrix matrix);
	void sortByVertexLvl(vector<VertexInList> &vertexInLvlList, int left, int right);
	void sortByIndex(vector<VertexInList> &vertexInLvlList, int left, int right);
	bool checkConnectionsInClique(vector<VertexInList> result, VertexInList analyzedVertex, Matrix matrix);
	DataFromFile();
	DataFromFile(string dataName, int substrLength, int reliability);
	~DataFromFile();
};

