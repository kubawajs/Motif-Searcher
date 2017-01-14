#pragma once
#include "DataFromFile.h"

class Result
{
	vector <Sequence> result;
	string motif;

public:
	void sortByIndex(vector<Vertex>& vertexInLvlList, int left, int right);
	void parseSequences(int reliability);
	vector <Sequence> getSequences();
	Result();
	Result(vector<Vertex> result, int numOfSeqs, string motif);
	~Result();
};

