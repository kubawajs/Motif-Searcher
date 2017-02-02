#include "stdafx.h"
#include "Vertex.h"

#include <vector>
#include <iostream>

using namespace std;

/**
 * \brief Return index of vertex in matrix
 * \return 
 */
int Vertex::getIndex() const
{
	return index;
}

/**
 * \brief Return index of vertex in sequence
 * \return 
 */
int Vertex::getIndexInSeq() const
{
	return indexInSeq;
}

/**
 * \brief Return index of vertex's sequence
 * \return 
 */
int Vertex::getSeqIndex() const
{
	return seqIndex;
}

/**
 * \brief Return level of vertex
 * \return 
 */
int Vertex::getVertexLvl() const
{
	return level;
}

/**
 * \brief Return number of connections with other sequences
 * \return 
 */
int Vertex::getConWithOtherSeq() const
{
	return conWithOtherSeq;
}

/**
 * \brief Return if vertex has minimum connections with other sequences
 * \return 
 */
bool Vertex::getHasMinConnections() const
{
	return hasMinConnections;
}

/**
 * \brief Return substring
 * \return 
 */
string Vertex::getSubstring() const
{
	return substring;
}

/**
 * \brief Return substring without nucleotides under reliability level
 * \return 
 */
string Vertex::getSubstringWithDeletions() const
{
	return substringWithDeletions;
}

/**
* \brief Return qual vector
* \return
*/
vector<int> Vertex::getQual() const
{
	return qual;
}

/**
 * \brief Return neighbours list
 * \return 
 */
vector<int> Vertex::getNeighboursList() const
{
	return neighboursList;
}

/**
 * \brief Set index of vertex
 * \param _index 
 */
void Vertex::setIndex(int _index)
{
	if(_index >= 0)
	{
		index = _index;
	}
	else
	{
		cout << "Can't assign index of vertex. Wrong value (lower than zero)." << endl;
		exit(0);
	}
}

/**
 * \brief Set index of vertex in sequence
 * \param _indexInSeq 
 */
void Vertex::setIndexInSeq(int _indexInSeq)
{
	if(_indexInSeq >= 0)
	{
		indexInSeq = _indexInSeq;
	}
	else
	{
		cout << "Can't assign index of vertex in sequence. Wrong value (lower than zero)." << endl;
		exit(0);
	}
}

/**
 * \brief Set index of sequence for vertex
 * \param _seqIndex 
 */
void Vertex::setSeqIndex(int _seqIndex)
{
	if(_seqIndex >= 0)
	{
		seqIndex = _seqIndex;
	}
	else
	{
		cout << "Can't assign index of sequence for vertex. Wrong value (lower than zero)." << endl;
		exit(0);
	}
}

/**
* \brief Set number of connections with other sequences
* \param _conWithOtherSeq
*/
void Vertex::setConWithOtherSeq(int _conWithOtherSeq)
{
	if(_conWithOtherSeq >= 0)
	{
		conWithOtherSeq = _conWithOtherSeq;
	}
	else
	{
		cout << "Can't assign number of connections with other sequences. Wrong value (number of connections lower than zero)" << endl;
		exit(0);
	}
}

/**
 * \brief Set if vertex has minimum connections w/ other sequences
 * \param _set 
 */
void Vertex::setHasMinConnections(bool _set)
{
	hasMinConnections = _set;
}

/**
 * \brief Set substring of vertex
 * \param _substring 
 */
void Vertex::setSubstring(string _substring)
{
	if(_substring.length() > 0)
	{
		substring = _substring;
	}
	else
	{
		cout << "Can't assign substring. Wrong value (substring is empty)." << endl;
		exit(0);
	}
}

/**
 * \brief Set qual vector for vertex
 * \param _qual 
 */
void Vertex::setQual(vector<int> _qual)
{
	if(!_qual.empty())
	{
		qual = _qual;
	}
	else
	{
		cout << "Can't assign quals. Wrong value (qual vector is empty)." << endl;
		exit(0);
	}
}

/**
 * \brief Set substring without deleted nucleotides
 * \param _substring 
 */
void Vertex::setSubstringWithDeletions(string _substring)
{
	substringWithDeletions = _substring;
}

/**
 * \brief Set list of neighbours
 * \param _neighboursList 
 */
void Vertex::setNeighboursList(vector<int> _neighboursList)
{
	neighboursList = _neighboursList;
}

/**
 * \brief Add neighbour vertex to neighbours list
 * \param _neighbourToAddID 
 */
void Vertex::addToNeighboursList(int _neighbourToAddID)
{
	neighboursList.push_back(_neighbourToAddID);
}

/**
* \brief Level up vertex level
*/
void Vertex::lvlUp()
{
	level += 1;
}

/**
 * \brief Overload operator ==
 * \param v 
 * \return 
 */
bool Vertex::operator==(const Vertex& v) const
{
	if(this->index == v.index)
	{
		return true;
	}
	return false;
}

/**
 * \brief Overload operator <
 * \param v 
 * \return 
 */
bool Vertex::operator<(const Vertex& v) const
{
	if(this->conWithOtherSeq > v.conWithOtherSeq)
	{
		return true;
	}
	if(this->conWithOtherSeq == v.conWithOtherSeq)
	{
		if(this->level > v.level)
		{
			return true;
		}
	}
	return false;
}

/**
 * \brief Constructor
 */
Vertex::Vertex() {}

/**
 * \brief Main parameter constructor
 * \param _index 
 * \param _indexInSeq 
 * \param _seqIndex 
 * \param _substring 
 * \param _qual 
 */
Vertex::Vertex(int _index, int _indexInSeq, int _seqIndex, string _substring, string _substringWithDeletions, vector<int> _qual)
{
	setIndex(_index);
	setIndexInSeq(_indexInSeq);
	setSeqIndex(_seqIndex);
	setSubstring(_substring);
	setSubstringWithDeletions(_substringWithDeletions);
	setQual(_qual);
	level = 0;
	hasMinConnections = false;
}

/**
 * \brief Destructor for object Vertex
 */
Vertex::~Vertex()
{
	substring.clear();
	substringWithDeletions.clear();
	qual.clear();
}
