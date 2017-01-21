#pragma once

#include "Sequence.h"

class ResultMotif
{
	vector <Sequence> result;
	vector <Vertex> startingClique;
	string motif;

public:
	void sortByIndex(vector<Vertex>& vertexInLvlList, int left, int right);
	void parseSequences(int reliability);
	void setStartingClique(vector <Vertex> clique);
	void setResult(vector <Vertex> result, int numOfSeqs);
	void setMotif(string motif);
	void printMotifOnSeq(vector <Vertex> result, int seqSize, int reliability);
	string getMotif();
	vector <Vertex> getStartingClique();
	vector <Sequence> getSequences();
	bool operator==(const ResultMotif &r);
	bool operator<(const ResultMotif &r);
	ResultMotif();
	~ResultMotif();
};
