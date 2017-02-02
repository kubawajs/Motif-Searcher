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
	string substring;
	string substringWithDeletions;
	vector<int> qual;
	vector<int> neighboursList;

public:
	//GETTERS
	int getIndex() const;
	int getIndexInSeq() const;
	int getSeqIndex() const;
	int getVertexLvl() const;
	int getConWithOtherSeq() const;
	bool getHasMinConnections() const;
	string getSubstring() const;
	string getSubstringWithDeletions() const;
	vector<int> getQual() const;
	vector<int> getNeighboursList() const;

	//SETTERS
	void setIndex(int _index);
	void setIndexInSeq(int indexInSeq);
	void setSeqIndex(int _seqIndex);
	void setConWithOtherSeq(int _conWithOtherSeq);
	void setHasMinConnections(bool _set);
	void setSubstring(string _substring);
	void setQual(vector<int> qual);
	void setSubstringWithDeletions(string _substring);
	void setNeighboursList(vector<int> _neighboursList);

	//METHODS
	void addToNeighboursList(int _neighbourToAddID);
	void lvlUp();

	//CONSTRUCTORS, DESTRUCTORS
	Vertex();
	Vertex(int _index, int _indexInSeq, int _seqIndex, string _substring, string _substringWithDeletions, vector<int> qual);
	~Vertex();

	//OPERATORS
	bool operator==(const Vertex& v) const;
	bool operator<(const Vertex& v) const;
};
