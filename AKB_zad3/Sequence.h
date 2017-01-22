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
	vector <int> qual; //quality
	vector <Vertex> substrings; //vector of all substrings from sequence

public:
	void setName(string name);
	string getName();
	void setSequence(string seq);
	string getSequence();
	void setQual(vector <int> qual);
	vector <int> getQual();
	void setSubstr(vector <Vertex> substrings);
	void setSeqId(int seqId);
	int getSeqId();
	void addSubstr(Vertex substring);
	void createSubstrings(string sequence, vector <int> qual, int substrLength, int reliability);
	int getSubstrSize();
	vector <Vertex> getSubstrings();
	Vertex getSubstrById(int noSubstr);
	void setVertexHasMinConnections(int noSubstr);
	void setVertexNumOfConSeq(int noSubstr, int conSeq);
	bool compareSubstrs(Vertex substr1, Vertex substr2, int usersSubstrLength, int reliability);//compare two substrings including deletions
	void vertexLvlUp(int noSubstr);//increment level of vertex
	Sequence();
	~Sequence();
};

