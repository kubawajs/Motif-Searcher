#pragma once

#include <vector>
#include <iostream>

using namespace std;

class Vertex
{
	vector <char> substring;
	vector <int> qual;
	int level;
	bool hasMinConnections;

public:
	void setSubstr(vector <char> substr);
	vector<char> getSubstring();
	void setQual(vector <int> qual);
	vector <int> getQual();
	bool getHasMinConnections();
	int getSubstrLength();
	int getVertexLvl();
	void printSubstr();//debugging time
	void printQual();
	void lvlUp(int n);
	void setHasMinConnections(bool set);
	Vertex();
	~Vertex();
};