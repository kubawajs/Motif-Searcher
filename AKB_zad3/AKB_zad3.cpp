// akb_zad3.cpp : Defines the entry point for the console application.
//
#pragma once

#include "stdafx.h"
#include "DataFromFile.h"

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

using namespace std;

int main()
{
	//TODO: input file name
	string fileName = "fasta"; // name of file in fasta (without .txt)
	int substrLength;
	int reliability;

	//USERS PARAMETERS
	do {
		cout << "Podaj dlugosc generowanych podciagow (4-7): " << endl;
		cin >> substrLength;
		if (substrLength < 4 || substrLength > 7) {
			cout << "Podano bledna dlugosc podciagu. Podaj wartosc z zakresu 4-7: " << endl;
			cin >> substrLength;
		}
	} while (substrLength < 4 || substrLength > 7);

	do {
		cout << "Podaj minimalny poziom wiarygodnosci (10-40): " << endl;
		cin >> reliability;
		if (reliability < 10 || reliability > 40) {
			cout << "Podano bledny poziom wiarygodnosci. Podaj minimalny poziom wiarygodnosci(10 - 40): " << endl;
			cin >> reliability;
		}
	} while (reliability < 10 || reliability > 40);

	//LOADING DATA	
	DataFromFile data(fileName, substrLength, reliability);//create object data
	data.loadFromFile(data.getDataName(), data.getSeqData());//load data from fasta file
	data.loadQualFile(data.getDataName(), data.getSeqData());//load data from qual file

	//CREATING GRAPH
	data.createGraph(data.getMatrix(), data.getSeqData());

	//MARKING TOO SHORT SUBSTRINGS IN MATRIX
	data.filterLowSubstrs(data.getMatrix(), data.getSeqData(), data.getReliability());

	//CREATING EDGES FOR NOT TO SHORT VERTICES
	data.createEdges(data.getMatrix(), data.getSeqData(), data.getInfoTable(data.getMatrix()), reliability);

	//CHECK IF EVERY VERTEX HAS CONNECTION WITH MIN SEQUENCES
	data.checkIfHasMinConnections(data.getMatrix());

	//TEST - PRINT SEQUENCES
	data.printSequences(data.getSeqData());

	//CREATE LIST OF VERTICES SORTED BY LEVEL
	data.createListOfVerticesSorted();

	//BUILD A CLIQUE BASED ON VERTICES LEVEL
	//TODO: musi budowac wiêcej ni¿ jedno rozwiazanie i wg rankingu wybierac najlepsze
	data.buildResults();

	//TODO: print result (max clique)
	//TODO: searching for clique/series of cliques
	//TODO: code refactoring
	//TODO: add conditions on all functions having access to private attributes

	cout << "dpa";

	return 0;
}
