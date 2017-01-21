#include "stdafx.h"
#include "Sequence.h"

using namespace std;

void Sequence::setName(string name)
{
	Sequence::name = name;
}

string Sequence::getName()
{
	return string(Sequence::name);
}

void Sequence::setSequence(string seq)
{
	Sequence::sequence = seq;
}

string Sequence::getSequence()
{
	return string(Sequence::sequence);
}

void Sequence::setQual(vector<int> qual)
{
	Sequence::qual = qual;
}

vector<int> Sequence::getQual()
{
	return Sequence::qual;
}

void Sequence::setSubstr(vector<Vertex> substrings)
{
	Sequence::substrings = substrings;
}

void Sequence::setSeqId(int seqId)
{
	Sequence::seqId = seqId;
}

int Sequence::getSeqId()
{
	return Sequence::seqId;
}

void Sequence::addSubstr(Vertex substring)
{
	Sequence::substrings.push_back(substring);
}

void Sequence::createSubstrings(string sequence, vector<int> qual, int substrLength, int reliability)
{
	Vertex v1; // temp object for data collection
	vector <char> tempSub;
	vector <int> tempQual;

	for (int x = 0; x <= (sequence.size() - substrLength); x++) {
		for (int i = 0; i < substrLength; i++) {
			tempSub.push_back(sequence[x + i]);
			tempQual.push_back(qual[x + i]);
		}
		v1.setSubstr(tempSub); // set substring in object
		v1.setQual(tempQual); // set qualities of all letters in substring
		v1.setIndexInSeq(x);
		Sequence::substrings.push_back(v1); // add substring to vector in sequence data
		//clear temp vectors for next substring
		tempSub.clear();
		tempQual.clear();
	}
}

int Sequence::getSubstrSize()
{
	return int(substrings.size());
}

vector<Vertex> Sequence::getSubstrings()
{
	return Sequence::substrings;
}

Vertex Sequence::getSubstrById(int noSubstr)
{
	return Sequence::substrings[noSubstr];
}

void Sequence::setVertexHasMinConnections(int noSubstr)
{
	Sequence::substrings[noSubstr].setHasMinConnections(true);
}

bool Sequence::compareSubstrs(Vertex substr1, Vertex substr2, int usersSubstrLength, int reliability)
{
	vector <char> substrSh, substrLn;
	int lengthSh, lengthLn, diffInLength, similarities = 0;

	vector <int> tempQ;
	vector <char> tempS;

	//set data at 1st sequence
	vector <int> qual = substr1.getQual();
	vector <char> substr = substr1.getSubstring();

	for (int i = 0; i < substr1.getSubstrLength(); i++) {
		if (qual[i] >= reliability) {
			tempQ.push_back(qual[i]);
			tempS.push_back(substr[i]);
		}
	}
	//save all data in object and clear temporary vectors
	substr1.setQual(tempQ);
	substr1.setSubstr(tempS);
	tempQ.clear();
	tempS.clear();
	//set data at 2nd sequence
	qual = substr2.getQual();
	substr = substr2.getSubstring();
	
	for (int i = 0; i < substr2.getSubstrLength(); i++) {
		if (qual[i] >= reliability) {
			tempQ.push_back(qual[i]);
			tempS.push_back(substr[i]);
		}
	}

	substr2.setQual(tempQ);
	substr2.setSubstr(tempS);

	//check which sequence is shorter
	if (substr1.getSubstrLength() >= substr2.getSubstrLength()) {
		substrSh = substr2.getSubstring();
		substrLn = substr1.getSubstring();
		lengthSh = substr2.getSubstrLength();
		lengthLn = substr1.getSubstrLength();
	}
	else {
		substrSh = substr1.getSubstring();
		substrLn = substr2.getSubstring();
		lengthSh = substr1.getSubstrLength();
		lengthLn = substr2.getSubstrLength();
	}
	diffInLength = lengthLn - lengthSh;
	for (int i = 0; i <= diffInLength; i++) {
		similarities = 0;
		for (int j = 0; j < lengthSh; j++) {
			if (substrSh[j] == substrLn[j + i]) {
				similarities++;
			}
		}
		if (similarities == lengthSh) {
			return true;
		}
	}
	return false;
}

void Sequence::vertexLvlUp(int noSubstr)
{
	Sequence::substrings[noSubstr].lvlUp(1);
}

Sequence::Sequence()
{
}


Sequence::~Sequence()
{
}