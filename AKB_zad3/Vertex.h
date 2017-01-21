#pragma once

#include <vector>
#include <iostream>

using namespace std;

class Vertex
{
	vector <char> substring;
	vector <int> qual;
	int level;
	int index;
	int indexInSeq;
	int seqIndex;
	bool hasMinConnections;

public:
	void setSubstr(vector <char> substr);
	vector<char> getSubstring();
	void setQual(vector <int> qual);
	vector <int> getQual();
	int getIndex();
	void setIndex(int _index);
	int getSeqIndex();
	void setSeqIndex(int _seqIndex);
	int getIndexInSeq();
	void setIndexInSeq(int indexInSeq);
	bool getHasMinConnections();
	int getSubstrLength();
	int getVertexLvl();
	void printSubstr();
	void printQual();
	void lvlUp(int n);
	void setHasMinConnections(bool set);
	bool operator==(const Vertex &v);
	Vertex();
	~Vertex();
};