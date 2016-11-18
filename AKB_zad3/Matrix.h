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
	void createEdge(int edgeIn, int edgeOut);
	void delVertex(int edge);
	int getSequenceIdFromMatrix(int index);
	Matrix();
	~Matrix();
};

