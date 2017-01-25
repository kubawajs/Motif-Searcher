#pragma once

#include "Matrix.h"
#include "Sequence.h"
#include "ResultMotif.h"

const int SENSITIVITY = 5;
const int NUMBER_OF_RESULTS = 10;
const int PERMITTED_DELETIONS = 2;

class DataFromFile
{
	vector <Sequence> seqData;
	vector <Vertex> vertexByLevel;
	vector <ResultMotif> results;
	Matrix matrix;
	string dataName;
	int substrLength;
	int reliability;
	int minConnections;

public:
	vector<Sequence> getSeqData() const;
	string getDataName() const;
	Matrix getMatrix() const;
	int getReliability() const;
	void loadFromFile(string dataName, vector<Sequence> seqData);
	void loadQualFile(string dataName, vector<Sequence> seqData);
	void createGraph(Matrix matrix, vector<Sequence> data);
	void filterLowSubstrs(Matrix matrix, vector<Sequence> data, int reliability);//filter substrings with lower qual and notice in matrix
	void setSeqData(vector<Sequence> seqData);
	void createEdges(Matrix matrix, vector<Sequence> data, vector<int> infoTable, int reliability);
	void checkIfHasMinConnections(Matrix matrix);//function for checking if vertex has connection with min 5 other sequences
	void printSequences(vector<Sequence> seqData, vector<Sequence> resultSeq);
	void createListOfVerticesSorted();
	void buildResults();
	void addResult(ResultMotif result);
	ResultMotif buildResult(vector<Vertex> startingClique);
	vector<ResultMotif> getResults() const;
	vector<Vertex> prepareVertexSetLeft(vector<Vertex> actualResult, int sensitivity);
	vector<Vertex> prepareVertexSetRight(vector<Vertex> actualResult, int sensitivity);
	vector<Vertex> buildClique(vector<Vertex> vertexByLevel) const;
	static vector<Vertex> sumResult(vector<Vertex> actualResult, vector<Vertex> tempResult);
	static vector<Vertex> filterVector(vector<Vertex> toFilter, vector<Vertex> filtering);
	static vector<int> getInfoTable(Matrix matrix);
	void sortResults(vector<ResultMotif>& results) const;
	static void printBestMotifs(vector<ResultMotif> results);
	static bool checkConnectionsInClique(vector<Vertex> result, Vertex analyzedVertex, Matrix matrix);
	string buildMotif(vector<Vertex> verticesToAlign, int reliability);
	static string parseMotifLeft(string existingMotif, string motifToAdd);
	static string parseMotifRight(string existingMotif, string motifToAdd);
	void printResult(vector<ResultMotif> result);
	DataFromFile();
	DataFromFile(string dataName, int substrLength, int reliability);
	~DataFromFile();
};
