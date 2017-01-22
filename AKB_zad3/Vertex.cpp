#include "stdafx.h"
#include "Vertex.h"

#include <vector>

using namespace std;

void Vertex::setSubstr(vector <char> substr)
{
	Vertex::substring = substr;
}

vector<char> Vertex::getSubstring()
{
	return vector<char>(Vertex::substring);
}

void Vertex::setQual(vector<int> qual)
{
	Vertex::qual = qual;
}

vector<int> Vertex::getQual()
{
	return Vertex::qual;
}

int Vertex::getIndex()
{
	return Vertex::index;
}

void Vertex::setIndex(int _index)
{
	Vertex::index = _index;
}

int Vertex::getSeqIndex()
{
	return Vertex::seqIndex;
}

void Vertex::setSeqIndex(int _seqIndex)
{
	Vertex::seqIndex = _seqIndex;
}

int Vertex::getIndexInSeq()
{
	return Vertex::indexInSeq;
}

void Vertex::setIndexInSeq(int indexInSeq)
{
	Vertex::indexInSeq = indexInSeq;
}

bool Vertex::getHasMinConnections()
{
	return Vertex::hasMinConnections;
}

int Vertex::getSubstrLength()
{
	return Vertex::substring.size();
}

int Vertex::getVertexLvl()
{
	return Vertex::level;
}

void Vertex::setConWithOtherSeq(int conWithOtherSeq)
{
	Vertex::conWithOtherSeq = conWithOtherSeq;
}

int Vertex::getConWithOtherSeq()
{
	return Vertex::conWithOtherSeq;
}

void Vertex::printSubstr()
{
	cout << "S: ";
	for (int i = 0; i < Vertex::getSubstrLength(); i++) {
		cout << Vertex::substring[i];
	}
	cout << " ";
}

void Vertex::printQual()
{
	cout << "Q: ";
	for (int i = 0; i < Vertex::getSubstrLength(); i++) {
		cout <<  Vertex::qual[i] << " ";
	}
	cout << endl;
}

void Vertex::lvlUp(int n)
{
	Vertex::level += n;
}

void Vertex::setHasMinConnections(bool set)
{
	Vertex::hasMinConnections = set;
}

bool Vertex::operator==(const Vertex &v)
{
	if (this->index == v.index)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool Vertex::operator<(const Vertex & v)
{
	if (this->conWithOtherSeq > v.conWithOtherSeq)
	{
		return true;
	}
	else if (this->conWithOtherSeq == v.conWithOtherSeq)
	{
		if (this->level > v.level)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

Vertex::Vertex()
{
	Vertex::level = 0;
	Vertex::index = -1;
	Vertex::hasMinConnections = false;
}


Vertex::~Vertex()
{
}
