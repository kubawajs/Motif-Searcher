#pragma once

#include "Sequence.h"
#include "Matrix.h"

using namespace std;

class DataFromFile
{
	vector <Sequence> seqData;
	Matrix matrix;
	string dataName;
	int substrLength;
	int reliability;

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
	void printSequences();
	vector<int> getInfoTable(Matrix matrix);
	DataFromFile();
	DataFromFile(string dataName, int substrLength, int reliability);
	~DataFromFile();
};

