#pragma once

#include "Sequence.h"

class ResultMotif
{
	vector <Sequence> result;
	vector <Vertex> startingClique;
	string motif;

public:
	void sortByIndex(vector<Vertex>& vertexInLvlList, int left, int right) const;
	void parseSequences(int reliability);
	void setStartingClique(vector <Vertex> clique);
	void setResult(vector <Vertex> result, int numOfSeqs);
	void setMotif(string motif);
	static void printMotifOnSeq(vector <Vertex> result, int seqSize, int reliability);
	static void printVerticesInMotif(vector<Vertex> result, int reliability);
	string getMotif() const;
	vector <Vertex> getStartingClique() const;
	vector <Sequence> getSequences() const;
	bool operator==(const ResultMotif &r) const;
	bool operator<(const ResultMotif &r) const;
	ResultMotif();
	~ResultMotif();
};
