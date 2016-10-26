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

int Vertex::getSubstrLength()
{
	return Vertex::substring.size();
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

Vertex::Vertex()
{
}


Vertex::~Vertex()
{
}
