#pragma once
#include "DataFromFile.h"

class Result
{
	vector <Sequence> result;
	vector <Vertex> startingClique;
	string motif;

public:
	void sortByIndex(vector<Vertex>& vertexInLvlList, int left, int right);
	void parseSequences(int reliability);
	void setStartingClique(vector <Vertex> clique);
	vector <Sequence> getSequences();
	Result();
	Result(vector<Vertex> result, int numOfSeqs, string motif);
	~Result();
};

