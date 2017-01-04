#pragma once
#include "DataFromFile.h"

class Result
{
	vector <Sequence> result;
	vector <string> alignedSequences;
	string motif;

public:
	Result();
	void sortByIndex(vector<Vertex>& vertexInLvlList, int left, int right);
	void parseSequences(int reliability);
	void alignSequences(vector <Sequence> &result);
	vector <Sequence> getSequences();
	Result(vector<Vertex> result, int numOfSeqs);
	~Result();
};

