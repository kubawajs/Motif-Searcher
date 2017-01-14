#include "stdafx.h"
#include "Result.h"


Result::Result()
{
}

void Result::sortByIndex(vector<Vertex> &vertexInLvlList, int left, int right)
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


void Result::parseSequences(int reliability)
{
	for (int i = 0; i < Result::result.size(); i++) //for each sequence
	{
		string sequence = "";
		vector<Vertex> substrs = Result::result[i].getSubstrings();

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
		Result::result[i].setSequence(sequence);
		//cout << sequence << endl;
	}
}

vector<Sequence> Result::getSequences()
{
	return Result::result;
}

Result::Result(vector<Vertex> result, int numOfSeqs, string motif)
{
	Sequence sequence;
	Result::sortByIndex(result, 0, result.size() - 1);

	for (int i = 0; i < numOfSeqs; i++)
	{
		Result::result.push_back(sequence);
	}

	for (int i = 0; i < result.size(); i++)
	{
		int seqIndex = result[i].getSeqIndex();
		Result::result[seqIndex].addSubstr(result[i]);
	}

	Result::motif = motif;
}


Result::~Result()
{
}
