#pragma once

#include "Sequence.h"
#include <vector>
#include <iostream>

using namespace std;

class Matrix
{
	vector <vector <int>> matrix;
	vector <int> infoTable;//table for index of start of each sequence in matrix
	int size;

public:
	int getSize();
	vector <int> getInfoTable();
	vector <vector <int>> getMatrix();
	void setMatrix(vector <vector <int>> matrix);
	void countMatrixSizeFromSeq(vector<Sequence> data, int size);
	void initializeMatrix(int size);
	void createEdge(vector<vector<int>> matrix, int edgeIn, int edgeOut);
	void delVertex(int edge);
	void createEdges(vector<vector<int>> matrix, vector<Sequence> data, vector<int> infoTable, int reliability);
	//TODO: create edges
	Matrix();
	~Matrix();
};

