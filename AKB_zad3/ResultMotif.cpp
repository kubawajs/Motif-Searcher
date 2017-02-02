#include "stdafx.h"
#include "ResultMotif.h"

#include <iostream>
#include <string>
#include <algorithm>

/**
* \brief Return motif
* \return
*/
string ResultMotif::getMotif() const
{
	return motif;
}

/**
 * \brief Return result sequences
 * \return 
 */
vector<Sequence> ResultMotif::getResult() const
{
	return result;
}

/**
* \brief Return starting clique
* \return
*/
vector<Vertex> ResultMotif::getStartingClique() const
{
	return startingClique;
}

/**
* \brief Return value of index from usedSequences
* \param index
* \return
*/
bool ResultMotif::getUsedSeqByIndex(int index)
{
	if(index < usedSequences.size() && index >= 0)
	{
		return usedSequences[index];
	}
	cout << "Can't return usedSequences[index]. Wrong index value." << endl;
	exit(0);
}

/**
* \brief Set motif
* \param motif
*/
void ResultMotif::setMotif(string motif)
{
	ResultMotif::motif = motif;
}

/**
 * \brief Set result sequences
 * \param _result 
 * \param numOfSeqs 
 */
void ResultMotif::setResult(vector<Vertex> _result, int numOfSeqs)
{
	if(!_result.empty())
	{
		Sequence sequence;

		sort(_result.begin(), _result.end(), [](const auto& A, const auto& B)
		     {
			     return A.getIndex() < B.getIndex();
		     });

		for(auto i = 0; i < numOfSeqs; i++)
		{
			sequence.setSeqId(i);
			result.push_back(sequence);
		}

		for(auto i = 0; i < _result.size(); i++)
		{
			auto seqIndex = _result[i].getSeqIndex();
			result[seqIndex].addSubstr(_result[i]);
		}
	}
}

/**
* \brief Set starting clique
* \param clique
*/
void ResultMotif::setStartingClique(vector<Vertex> clique)
{
	startingClique = clique;
}

/**
* \brief Mark sequence as used
* \param index
*/
void ResultMotif::markSequence(int index)
{
	usedSequences[index] = true;
}

/**
* \brief Print motif based on original sequences
* \param result
* \param seqSize
* \param reliability
*/
void ResultMotif::printMotifOnSeq(vector<Vertex> result, int seqSize, int reliability)
{
	string sequence(seqSize, '-');

	for(auto j = 0; j < result.size(); j++)
	{
		for(auto k = 0; k < result[j].getSubstring().length(); k++)
		{
			if(result[j].getQual()[k] >= reliability)
			{
				sequence[result[j].getIndexInSeq() + k] = result[j].getSubstring()[k];
			}
			else
			{
				sequence[result[j].getIndexInSeq() + k] = '*';
			}
		}
	}

	cout << "M: " << sequence << "\n" << endl;
}

/**
 * \brief Print vertices used for build result
 * \param result 
 * \param reliability 
 */
void ResultMotif::printVerticesInMotif(vector<Vertex> result, int reliability)
{
	for(auto j = 0; j < result.size(); j++)
	{
		cout << result[j].getSeqIndex() << ":" << result[j].getIndexInSeq() << " ";
		for(auto k = 0; k < result[j].getSubstring().length(); k++)
		{
			cout << result[j].getSubstring()[k];
		}
		cout << " ";
	}
	cout << endl;
}

/**
* \brief Reset used sequences
* \param size
*/
void ResultMotif::resetUsedSequences(int size)
{
	usedSequences.clear();
	for(auto i = 0; i < size; i++)
	{
		usedSequences.push_back(false);
	}
}

/**
 * \brief Overload operator ==
 * \param r 
 * \return 
 */
bool ResultMotif::operator==(const ResultMotif& r) const
{
	if(this->motif == r.motif)
	{
		return true;
	}
	return false;
}

/**
 * \brief Overload operator <
 * \param r 
 * \return 
 */
bool ResultMotif::operator<(const ResultMotif& r) const
{
	if(this->motif.size() < r.motif.size())
	{
		return true;
	}
	return false;
}

ResultMotif::ResultMotif() {}

/**
 * \brief Destructor
 */
ResultMotif::~ResultMotif()
{
	startingClique.clear();
	usedSequences.clear();
	motif.clear();
}
