#pragma once

#include <vector>

using namespace std;

class Vertex
{
	int index;
	int indexInSeq;
	int seqIndex;
	int level;
	int conWithOtherSeq;
	bool hasMinConnections;
	vector <char> substring;
	vector <int> qual;

public:
	void setSubstr(vector <char> substr);
	vector<char> getSubstring() const;
	void setQual(vector <int> qual);
	vector <int> getQual() const;
	int getIndex() const;
	void setIndex(int _index);
	int getSeqIndex() const;
	void setSeqIndex(int _seqIndex);
	int getIndexInSeq() const;
	void setIndexInSeq(int indexInSeq);
	bool getHasMinConnections() const;
	int getSubstrLength() const;
	int getVertexLvl() const;
	void setConWithOtherSeq(int conWithOtherSeq);
	int getConWithOtherSeq() const;
	void printSubstr();
	void printQual();
	void lvlUp(int n);
	void setHasMinConnections(bool set);
	bool operator==(const Vertex &v) const;
	bool operator<(const Vertex &v) const;
	Vertex();
	~Vertex();
};