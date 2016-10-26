#pragma once

#include <vector>
#include <iostream>

using namespace std;

class Vertex
{
	vector <char> substring;
	vector <int> qual;
	int level;

public:
	void setSubstr(vector <char> substr);
	vector<char> getSubstring();
	void setQual(vector <int> qual);
	vector <int> getQual();
	int getSubstrLength();
	void printSubstr();//debugging time
	void printQual();
	void lvlUp(int n);
	Vertex();
	~Vertex();
};