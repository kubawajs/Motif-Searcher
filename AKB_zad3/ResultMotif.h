#pragma once

#include "Sequence.h"

class ResultMotif
{
	vector<Sequence> result;
	vector<Vertex> startingClique;
	vector<bool> usedSequences;
	string motif;

public:
	//GETTERS
	string getMotif() const;
	vector<Sequence> getResult() const;
	vector<Vertex> getStartingClique() const;
	bool getUsedSeqByIndex(int index);

	//SETTERS
	void setMotif(string motif);
	void setResult(vector<Vertex> result, int numOfSeqs);
	void setStartingClique(vector<Vertex> clique);

	//METHODS
	void markSequence(int index);
	static void printMotifOnSeq(vector<Vertex> result, int seqSize, int reliability);
	static void printVerticesInMotif(vector<Vertex> result, int reliability);
	void resetUsedSequences(int size);

	//OPERATORS
	bool operator==(const ResultMotif& r) const;
	bool operator<(const ResultMotif& r) const;

	//CONSTRUCTORS, DESTRUCTORS
	ResultMotif();
	~ResultMotif();
};
