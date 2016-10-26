#include "stdafx.h"
#include "Matrix.h"

using namespace std;

int Matrix::getSize()
{
	return Matrix::size;
}

vector<int> Matrix::getInfoTable()
{
	return Matrix::infoTable;
}

vector<vector<int>> Matrix::getMatrix()
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

void Matrix::createEdge(vector <vector <int>> matrix, int edgeIn, int edgeOut)
{
	matrix[edgeIn][edgeOut] = 1;
	matrix[edgeOut][edgeIn] = 1;
	Matrix::setMatrix(matrix);
}

void Matrix::delVertex(int edge)
{
	int limit = Matrix::size;
	for (int i = 0; i < limit; i++) {
		Matrix::matrix[i][edge] = -1;
		Matrix::matrix[edge][i] = -1;
	}
}

void Matrix::createEdges(vector <vector <int>> matrix, vector <Sequence> data, vector <int> infoTable, int reliability) {
	int actual = 0, createdEdges = 0;
	string waiting = "";

	cout << "Tworze polaczenia, prosze czekac...";

	for (int i = 0; i < matrix.size(); i++) {
		if (i == infoTable[actual]) {
			actual++;
		}
		if (matrix[i][0] != -1) {
			for (int j = infoTable[actual]; j < matrix.size(); j++) {
				if (matrix[0][j] != -1) {
					int noSeq1, noSubstr1, noSeq2, noSubstr2;
					noSeq1 = actual;//check 1st sequence
					if (actual == 0) {//check which substring of actual sequence is being analyzed
						noSubstr1 = i;
					}
					else {
						noSubstr1 = i - infoTable[actual - 1];
					}
					for (int k = 0; k < infoTable.size(); k++) {//check 2nd sequence
						if (j < infoTable[k]) {
							noSeq2 = k;
							if (k == 0) {//check which substring of actual sequence is being analyzed
								noSubstr2 = j;
							}
							else {
								noSubstr2 = j - infoTable[k - 1];
							}
							break;
						}
					}
					Vertex v1, v2;
					v1 = data[noSeq1].getSubstrById(noSubstr1);//get first vertex
					v2 = data[noSeq2].getSubstrById(noSubstr2);//get 2nd vertex

					if (data[noSeq1].compareSubstrs(v1, v2, v1.getSubstrLength(), reliability)) {
						Matrix::createEdge(matrix, i, j);
						createdEdges++;
					}
				}
			}
		}
	}
	cout << endl << "Utworzono polaczen: " << createdEdges << endl;
}

Matrix::Matrix()
{
}


Matrix::~Matrix()
{
}
