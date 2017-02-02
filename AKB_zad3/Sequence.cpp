#include "stdafx.h"
#include "Sequence.h"

using namespace std;

/**
 * \brief Set name of Sequence object
 * \param _name 
 */
void Sequence::setName(string _name)
{
	if(!_name.empty())
		name = _name;
	else
		name = "No name";
}

/**
 * \brief Set sequence of Sequence object
 * \param _sequence 
 */
void Sequence::setSequence(string _sequence)
{
	if(!_sequence.empty())
		sequence = _sequence;
	else
	{
		cout << "Can't assign sequence to object Sequence.";
		exit(0);
	}
}

/**
 * \brief Set qual vector of Sequence object
 * \param _qual 
 */
void Sequence::setQual(vector<int> _qual)
{
	if(!_qual.empty())
		qual = _qual;
	else
	{
		cout << "Can't assign qual. Qual vector is empty.";
		exit(0);
	}
}

/**
 * \brief Set number of connections with other sequence of input vertex
 * \param noSubstr 
 * \param conSeq 
 */
void Sequence::setVertexNumOfConSeq(int noSubstr, int conSeq)
{
	substrings[noSubstr].setConWithOtherSeq(conSeq);
}

/**
 * \brief Set sequence ID
 * \param _id 
 */
void Sequence::setSeqId(int _id)
{
	if(_id >= 0)
		seqId = _id;
	else
	{
		cout << "Invalid seqId." << endl;
		exit(0);
	}
}

/**
 * \brief Set substrings of Sequence object
 * \param _substrings 
 */
void Sequence::setSubstrings(vector<Vertex> _substrings)
{
	if(!_substrings.empty())
		substrings = _substrings;
	else
	{
		cout << "Can't assign substrings. Substring vector is empty.";
		exit(0);
	}
}

/**
 * \brief Return sequence id attribute
 * \return 
 */
int Sequence::getSeqId() const
{
	return seqId;
}

/**
 * \brief Return name of Sequence object
 * \return 
 */
string Sequence::getName() const
{
	return name;
}

/**
 * \brief Return sequence of Sequence object
 * \return 
 */
string Sequence::getSequence() const
{
	return sequence;
}

/**
 * \brief Return qual vector of Sequence object
 * \return 
 */
vector<int> Sequence::getQual() const
{
	return qual;
}

/**
 * \brief Return substrings of Sequence object
 * \return 
 */
vector<Vertex> Sequence::getSubstrings() const
{
	return substrings;
}

/**
 * \brief Return reference to substring in sequence
 * \param substrID 
 * \return 
 */
Vertex* Sequence::getSubstrById(int substrID)
{
	return &substrings[substrID];
}

/**
 * \brief 
 * \param substring 
 */
void Sequence::addSubstr(Vertex substring)
{
	substrings.push_back(substring);
}

/**
* \brief Create substrings from sequence and save into Sequence.substrings
* \param _sequence
* \param _qual
* \param _seqId
* \param _startingIndex
* \param _substrLength
* \param _reliability
*/
void Sequence::createSubstrings(string _sequence, vector<int> _qual, int _seqId, int _startingIndex, int _substrLength, int _reliability)
{
	string tempSub = "", tempSubWDels = "";
	vector<int> tempQual;

	for(auto x = 0; x <= (_sequence.length() - _substrLength); x++)
	{
		for(auto i = 0; i < _substrLength; i++)
		{
			tempSub += _sequence[x + i];
			if(_qual[x + i] >= _reliability)
			{
				tempSubWDels += _sequence[x + i];
			}
			tempQual.push_back(_qual[x + i]);
		}

		Vertex v1(x + _startingIndex, x, _seqId, tempSub, tempSubWDels, tempQual);
		substrings.push_back(v1);

		tempSub = "";
		tempSubWDels = "";
		tempQual.clear();
	}
}

/**
 * \brief Without params constructor
 */
Sequence::Sequence()
{
	seqId = 0;
	name = "No name";
}

/**
 * \brief Sequence element destructor
 */
Sequence::~Sequence()
{
	name.clear();
	sequence.clear();
	qual.clear();
	substrings.clear();
}
