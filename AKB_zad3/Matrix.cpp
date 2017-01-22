#include "stdafx.h"
#include "Matrix.h"

using namespace std;

int Matrix::getSize() const
{
	return Matrix::size;
}

vector<int> Matrix::getInfoTable() const
{
	return Matrix::infoTable;
}

vector<vector<int>> Matrix::getMatrix() const
{
	return Matrix::matrix;
}

void Matrix::setMatrix(vector<vector<int>> matrix)
{
	Matrix::matrix = matrix;
}

void Matrix::countMatrixSizeFromSeq(vector<Sequence> data, int size)
{
	for (int i = 0; i < data.size(); i++) {
		size += data[i].getSubstrSize();
		Matrix::infoTable.push_back(size);
	}
	Matrix::size = size;
}

void Matrix::initializeMatrix(int size)
{
	vector <vector <int> > tempMatrix(size, vector <int>(size));
	Matrix::matrix = tempMatrix;
}

void Matrix::createEdge(int edgeIn, int edgeOut)
{
	Matrix::matrix[edgeIn][edgeOut] = 1;
	Matrix::matrix[edgeOut][edgeIn] = 1;
}

void Matrix::delVertex(int edge)
{
	int limit = Matrix::size;
	for (int i = 0; i < limit; i++) {
		Matrix::matrix[i][edge] = -1;
		Matrix::matrix[edge][i] = -1;
	}
}

int Matrix::getSequenceIdFromMatrix(int index)
{
	int seqId = 0;
	while (index >= infoTable[seqId]) {
		seqId++;
	}
	return seqId;
}

Matrix::Matrix()
{
}


Matrix::~Matrix()
{
}
