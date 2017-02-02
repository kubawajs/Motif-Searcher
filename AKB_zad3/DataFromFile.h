#pragma once

#include "Sequence.h"
#include "ResultMotif.h"

const int SENSITIVITY = 3;
const int ITERATIONS = 5;
const int PERMITTED_DELETIONS = 2;
const double SEQUENCES_PERCENT = 0.60;

class DataFromFile
{
	string dataName;
	int substrLength;
	int reliability;
	int minConnections;
	vector<Sequence> seqData;
	vector<Vertex> vertexByLevel;
	vector<ResultMotif> results;

public:
	//GETTERS
	string getDataName() const;
	int getReliability() const;
	Vertex getVertexByID(int ID);
	vector<Sequence> getSeqData() const;
	vector<ResultMotif> getResults() const;

	//SETTERS
	void setSeqData(vector<Sequence> _seqData);
	void setVertexConWithOtherSeq(Vertex& vertex);

	//METHODS
	void addResult(ResultMotif result);
	vector<Vertex> buildClique(vector<Vertex> candidates) const;
	string buildMotif(vector<Vertex> verticesToAlign);
	ResultMotif buildResult();
	void buildResults();
	static bool checkConnectionsInClique(vector<Vertex> clique, Vertex toCheck);
	void createEdges();
	void createGraph();
	void createListOfVerticesSorted();
	static vector<Vertex> filterVector(vector<Vertex> toFilter, vector<Vertex> filtering);
	void loadFromFile(string dataName, vector<Sequence> seqData);
	void loadQualFile(string dataName, vector<Sequence> seqData);
	static string parseMotifLeft(string genMotif, string toParse);
	static string parseMotifRight(string genMotif, string toParse);
	vector<Vertex> prepareVertexSetLeft(vector<Vertex> actualResult);
	vector<Vertex> prepareVertexSetRight(vector<Vertex> actualResult);
	void printSequences(vector<Sequence> seqData, vector<Sequence> resultSeq);
	void printResult();
	vector<Sequence> sequencesIntoSubstrings(vector<Sequence> _data) const;
	void sortResults(vector<ResultMotif>& results) const;
	static vector<Vertex> sumResult(vector<Vertex> actualResult, vector<Vertex> tempResult);

	//CONSTRUCTORS, DESTRUCTORS
	DataFromFile(string _dataName, int _substrLength, int _reliability);
	~DataFromFile();
};
