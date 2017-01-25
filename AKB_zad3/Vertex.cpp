#include "stdafx.h"
#include "Vertex.h"

#include <vector>
#include <iostream>

using namespace std;

void Vertex::setSubstr(vector <char> substr)
{
	substring = substr;
}

vector<char> Vertex::getSubstring() const
{
	return substring;
}

void Vertex::setQual(vector<int> _qual)
{
	qual = _qual;
}

vector<int> Vertex::getQual() const
{
	return qual;
}

int Vertex::getIndex() const
{
	return index;
}

void Vertex::setIndex(int _index)
{
	Vertex::index = _index;
}

int Vertex::getSeqIndex() const
{
	return seqIndex;
}

void Vertex::setSeqIndex(int _seqIndex)
{
	seqIndex = _seqIndex;
}

int Vertex::getIndexInSeq() const
{
	return indexInSeq;
}

void Vertex::setIndexInSeq(int _indexInSeq)
{
	indexInSeq = _indexInSeq;
}

bool Vertex::getHasMinConnections() const
{
	return hasMinConnections;
}

int Vertex::getSubstrLength() const
{
	return substring.size();
}

int Vertex::getVertexLvl() const
{
	return level;
}

void Vertex::setConWithOtherSeq(int _conWithOtherSeq)
{
	conWithOtherSeq = _conWithOtherSeq;
}

int Vertex::getConWithOtherSeq() const
{
	return conWithOtherSeq;
}

void Vertex::printSubstr()
{
	cout << "S: ";
	for (int i = 0; i < getSubstrLength(); i++) {
		cout << substring[i];
	}
	cout << " ";
}

void Vertex::printQual()
{
	cout << "Q: ";
	for (int i = 0; i <getSubstrLength(); i++) {
		cout <<  qual[i] << " ";
	}
	cout << endl;
}

void Vertex::lvlUp(int n)
{
	level += n;
}

void Vertex::setHasMinConnections(bool set)
{
	hasMinConnections = set;
}

bool Vertex::operator==(const Vertex &v) const
{
	if (this->index == v.index)
	{
		return true;
	}
	return false;
}

bool Vertex::operator<(const Vertex & v) const
{
	if (this->conWithOtherSeq > v.conWithOtherSeq)
	{
		return true;
	}
	if (this->conWithOtherSeq == v.conWithOtherSeq)
	{
		if (this->level > v.level)
		{
			return true;
		}
		return false;
	}
	return false;
}

Vertex::Vertex()
{
	level = 0;
	hasMinConnections = false;
	conWithOtherSeq = 0;
}


Vertex::~Vertex()
{
}
