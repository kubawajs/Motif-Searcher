// akb_zad3.cpp : Defines the entry point for the console application.
//
#pragma once

#include "stdafx.h"
#include "iostream"
#include "string"

#include "DataFromFile.h"

using namespace std;

int main()
{
	//DEFAULT PARAMETERS
	string FILE_NAME = "fasta"; // name of file in fasta (without .txt)
	int SUBSTRING_LENGTH, RELIABILITY;

	//USER PARAMETERS
	cout << "Podaj nazwe pliku fasta (bez rozszerzenia): " << endl;
	cin >> FILE_NAME;

	do
	{
		cout << "Podaj dlugosc generowanych podciagow (4-7): " << endl;
		cin >> SUBSTRING_LENGTH;
		if(SUBSTRING_LENGTH < 4 || SUBSTRING_LENGTH > 7)
		{
			cout << "Podano bledna dlugosc podciagu. Podaj wartosc z zakresu 4-7: " << endl;
			cin >> SUBSTRING_LENGTH;
		}
	} while(SUBSTRING_LENGTH < 4 || SUBSTRING_LENGTH > 7);

	do
	{
		cout << "Podaj minimalny poziom wiarygodnosci (10-40): " << endl;
		cin >> RELIABILITY;
		if(RELIABILITY < 10 || RELIABILITY > 40)
		{
			cout << "Podano bledny poziom wiarygodnosci. Podaj minimalny poziom wiarygodnosci(10 - 40): " << endl;
			cin >> RELIABILITY;
		}
	} while(RELIABILITY < 10 || RELIABILITY > 40);

	//LOADING DATA	
	DataFromFile data(FILE_NAME, SUBSTRING_LENGTH, RELIABILITY);

	//CREATING GRAPH
	data.createGraph();

	//CREATING EDGES FOR NOT TO SHORT VERTICES
	data.createEdges();

	//CREATE LIST OF VERTICES SORTED BY LEVEL
	//najdluzej trwa
	data.createListOfVerticesSorted();

	//BUILD A CLIQUE BASED ON VERTICES LEVEL
	data.buildResults();

	////PRINT RESULT
	data.printResult();

	string wait;
	cout << "Press any letter to end";
	cin >> wait;

	return 0;
}
