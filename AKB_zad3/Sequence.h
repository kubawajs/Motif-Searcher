#pragma once

#include "Vertex.h"

#include <iostream>
#include <vector>

using namespace std;

class Sequence
{
	int seqId;
	string name; //sequence name
	string sequence; //whole sequence as one string
	vector<int> qual; //quality
	vector<Vertex> substrings; //vector of all substrings from sequence

public:
	//GETTERS
	int getSeqId() const;
	string getName() const;
	string getSequence() const;
	vector<int> getQual() const;
	vector<Vertex> getSubstrings() const;
	Vertex* getSubstrById(int substrID);

	//SETTERS
	void setSeqId(int _id);
	void setName(string _name);
	void setSequence(string _seq);
	void setQual(vector<int> _qual);
	void setSubstrings(vector<Vertex> _substrings);

	//METHODS
	void createSubstrings(string sequence, vector<int> qual, int seqId,
	                      int startingIndex, int substrLength, int reliability);

	//CONSTRUCTORS, DESTRUCTORS
	Sequence();
	~Sequence();

	//BEFORE
	void addSubstr(Vertex substring);
	void setVertexNumOfConSeq(int noSubstr, int conSeq);
};
