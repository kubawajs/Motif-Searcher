#include "stdafx.h"
#include "ResultMotif.h"

#include <iostream>
#include <string>
#include <algorithm>

ResultMotif::ResultMotif()
{
}

ResultMotif::~ResultMotif()
{
}

void ResultMotif::sortByIndex(vector<Vertex> &vertexInLvlList, int left, int right) const
{
	int i = left;
	int j = right;
	int x = vertexInLvlList[(left + right) / 2].getIndex();
	do {
		while (vertexInLvlList[i].getIndex() < x)
			i++;
		while (vertexInLvlList[j].getIndex() > x)
			j--;
		if (i <= j) {
			swap(vertexInLvlList[i], vertexInLvlList[j]);
			i++;
			j--;
		}
	} while (i <= j);

	if (left < j) sortByIndex(vertexInLvlList, left, j);
	if (right > i) sortByIndex(vertexInLvlList, i, right);
}


void ResultMotif::parseSequences(int reliability)
{
	for (int i = 0; i < result.size(); i++) //for each sequence
	{
		string sequence = "";
		vector<Vertex> substrs = result[i].getSubstrings();

		for (int j = 0; j < substrs.size(); j++)//for each substring of sequence
		{
			int indexDiff;

			if (j == 0)
			{
				indexDiff = substrs[j].getSubstrLength();
			}
			else {
				indexDiff = substrs[j].getIndex() - substrs[j - 1].getIndex();
				if (indexDiff > substrs[j].getSubstrLength())
				{
					indexDiff = substrs[j].getSubstrLength();
				}
			}

			int startCopying = substrs[j].getSubstrLength() - indexDiff;

			for (int k = startCopying; k < substrs[j].getSubstrLength(); k++) { //for each char in substring
				if (substrs[j].getQual()[k] >= reliability)
				{
					sequence += substrs[j].getSubstring()[k];
				}
			}
			//cut start of string according to indexDiff
		}
		result[i].setSequence(sequence);
		//cout << sequence << endl;
	}
}

void ResultMotif::setStartingClique(vector<Vertex> clique)
{
	startingClique = clique;
}

void ResultMotif::setResult(vector<Vertex> result, int numOfSeqs)
{
	if(!result.empty())
	{
		Sequence sequence;

		std::sort(result.begin(), result.end(), [](const auto& A, const auto& B) {
			return A.getIndex() < B.getIndex(); });

		for (int i = 0; i < numOfSeqs; i++)
		{
			sequence.setSeqId(i);
			ResultMotif::result.push_back(sequence);
		}

		for (int i = 0; i < result.size(); i++)
		{
			int seqIndex = result[i].getSeqIndex();
			ResultMotif::result[seqIndex].addSubstr(result[i]);
		}
	}
}

void ResultMotif::setMotif(string motif)
{
	ResultMotif::motif = motif;
}

void ResultMotif::resetUsedSequences(int size)
{
	for(int i=0; i<size; i++)
	{
		usedSequences.push_back(false);
	}
}

void ResultMotif::markSequence(int index)
{
	usedSequences[index] = true;
}

bool ResultMotif::getUsedSeqByIndex(int index)
{
	return usedSequences[index];
}

void ResultMotif::printMotifOnSeq(vector<Vertex> result, int seqSize, int reliability)
{
	string sequence(seqSize, '-');

	for (int j = 0; j < result.size(); j++)
	{
		for (int k = 0; k < result[j].getSubstrLength(); k++)
		{ //for each char in substring
			if (result[j].getQual()[k] >= reliability) {
				sequence[result[j].getIndexInSeq() + k] = result[j].getSubstring()[k];
			}
			else {
				sequence[result[j].getIndexInSeq() + k] = '*';
			}
		}
	}
	
	cout << "M: " << sequence << "\n" << endl;
}

void ResultMotif::printVerticesInMotif(vector<Vertex> result, int reliability)
{
	for (int j = 0; j < result.size(); j++)
	{
		cout << result[j].getSeqIndex() << ":" << result[j].getIndexInSeq() << " ";
		for (int k = 0; k < result[j].getSubstrLength(); k++)
		{ //for each char in substring
			cout << result[j].getSubstring()[k];
		}
		cout << " ";
	}
	cout << endl;
}

string ResultMotif::getMotif() const
{
	return motif;
}

vector<Vertex> ResultMotif::getStartingClique() const
{
	return startingClique;
}

vector<Sequence> ResultMotif::getSequences() const
{
	return result;
}

bool ResultMotif::operator==(const ResultMotif&r) const
{
	if (this->motif == r.motif)
	{
		return true;
	}
	return false;
}

bool ResultMotif::operator<(const ResultMotif & r) const
{
	if (this->motif.size() < r.motif.size())
	{
		return true;
	}
	return false;
}
