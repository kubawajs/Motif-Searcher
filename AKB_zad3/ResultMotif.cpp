#include "stdafx.h"
#include "ResultMotif.h"

#include <iostream>
#include <string>

ResultMotif::ResultMotif()
{
}

ResultMotif::~ResultMotif()
{
}

void ResultMotif::sortByIndex(vector<Vertex> &vertexInLvlList, int left, int right)
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
	for (int i = 0; i < ResultMotif::result.size(); i++) //for each sequence
	{
		string sequence = "";
		vector<Vertex> substrs = ResultMotif::result[i].getSubstrings();

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
		ResultMotif::result[i].setSequence(sequence);
		//cout << sequence << endl;
	}
}

void ResultMotif::setStartingClique(vector<Vertex> clique)
{
	ResultMotif::startingClique = clique;
}

void ResultMotif::setResult(vector<Vertex> result, int numOfSeqs)
{
	Sequence sequence;
	ResultMotif::sortByIndex(result, 0, result.size() - 1);

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

void ResultMotif::setMotif(string motif)
{
	ResultMotif::motif = motif;
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

string ResultMotif::getMotif()
{
	return ResultMotif::motif;
}

vector<Vertex> ResultMotif::getStartingClique()
{
	return ResultMotif::startingClique;
}

vector<Sequence> ResultMotif::getSequences()
{
	return ResultMotif::result;
}

bool ResultMotif::operator==(const ResultMotif&r)
{
	if (this->motif == r.motif)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool ResultMotif::operator<(const ResultMotif & r)
{
	if (this->motif.size() < r.motif.size())
	{
		return true;
	}
	else
	{
		return false;
	}
}
