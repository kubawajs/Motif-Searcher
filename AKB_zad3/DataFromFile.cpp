#include "stdafx.h"
#include "DataFromFile.h"


#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

vector<Sequence> DataFromFile::getSeqData()
{
	return DataFromFile::seqData;
}

string DataFromFile::getDataName()
{
	return DataFromFile::dataName;
}

Matrix DataFromFile::getMatrix()
{
	return DataFromFile::matrix;
}

int DataFromFile::getReliability()
{
	return DataFromFile::reliability;
}

void DataFromFile::loadFromFile(string dataName, vector <Sequence> seqData) {
	fstream file;
	string line;
	string lastSeq = "";
	string path = dataName + ".txt";

	file.open(path);
	if (!file.good()) {
		cout << "Error! Can't open the file." << endl;
		exit(0);
	}
	else {
		Sequence s1;
		while (getline(file, line)) {
			if (line[0] == '>') {
				s1.setName(line);
			}
			else {
				if (lastSeq == "") {
					lastSeq = line;
				}
				else {
					lastSeq += line;
					s1.setSequence(lastSeq);
					seqData.push_back(s1);
					lastSeq = "";
				}
			}
		}
	}
	DataFromFile::seqData = seqData;
}

void DataFromFile::loadQualFile(string dataName, vector<Sequence> seqData)
{
	fstream file;
	string line;
	vector <int> qualities;
	string number = "";
	int seqNo = 0;
	string path = dataName + "-qual.txt";

	file.open(path);
	if (!file.good()) {
		cout << "Error! Can't open the file." << endl;
		exit(0);
	}
	else {
		while (getline(file, line)) {
			if (line[0] == '>') {
				if (line != seqData[seqNo].getName()) {
					cout << "Error! Fasta file and qual file not equal." << endl;
					exit(0);
				}
				else {
					if (!qualities.empty()) {
						seqData[seqNo - 1].setQual(qualities);
						qualities.clear();
					}
				}
				seqNo++;
			}
			else {
				for (int i = 0; i < line.size(); i++) {
					if (line[i] == ' ' && number != "") {
						qualities.push_back(atoi(number.c_str()));
						number = "";
					}
					else if (line[i] != ' ' && number == "") {
						number = line[i];
					}
					else if (line[i] != ' ' && number != "") {
						number += line[i];
						qualities.push_back(atoi(number.c_str()));
						number = "";
					}
				}
			}
		}
		if (seqNo == seqData.size()) {
			seqData[seqNo - 1].setQual(qualities);
			qualities.clear();
		}
	}
	DataFromFile::seqData = seqData;
}

void DataFromFile::createGraph(Matrix matrix, vector <Sequence> data)
{
	for (int i = 0; i < data.size(); i++) {
		data[i].createSubstrings(data[i].getSequence(), data[i].getQual(), substrLength, reliability);//create all substrings of sequence
	}
	matrix.countMatrixSizeFromSeq(data, 0);//sum all sequences
	matrix.initializeMatrix(matrix.getSize());
	DataFromFile::seqData = data;
	DataFromFile::matrix = matrix;
}

void DataFromFile::filterLowSubstrs(Matrix matrix, vector<Sequence> data, int reliability)
{
	vector <Vertex> checking;
	vector <int> checkQuals, infoTable;
	int countGoodNucs, edge, limit;
	infoTable = matrix.getInfoTable();
	//loop over sequences
	for (int i = 0; i < data.size(); i++) {
		//loop over substrings in sequence
		for (int j = 0; j < data[i].getSubstrSize(); j++) {
			countGoodNucs = 0;
			checking = data[i].getSubstrings();
			checkQuals = checking[j].getQual();
			limit = (checkQuals.size() / 2) + 1; //min number of good nucs in substring
			//loop over elements of substring
			for (int k = 0; k < checkQuals.size(); k++) {
				if (checkQuals[k] >= reliability) {
					countGoodNucs++;
					if (countGoodNucs >= limit) {
						break;
					}
				}
				if (k == (limit - 1) && countGoodNucs == 0) {
					break;
				}
			}
			if (countGoodNucs < limit) {
				if (i == 0) {
					edge = j;
				}
				else {
					edge = j + infoTable[i - 1];
				}
				matrix.delVertex(edge);
			}
		}
	}
	DataFromFile::matrix = matrix;
}

DataFromFile::DataFromFile()
{
}

DataFromFile::DataFromFile(string dataName, int substrLength, int reliability)
{
	DataFromFile::dataName = dataName;
	DataFromFile::substrLength = substrLength;
	DataFromFile::reliability = reliability;
}

DataFromFile::~DataFromFile()
{
}
