#include "stdafx.h"
#include "DataFromFile.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <map>
#include <iterator>
#include <sstream>
#include <iomanip>

using namespace std;

/**
 * \brief Return dataName
 * \return
 */
string DataFromFile::getDataName() const
{
	return dataName;
}

/**
 * \brief Return reliability
 * \return
 */
int DataFromFile::getReliability() const
{
	return reliability;
}

/**
 * \brief Return vertex object by index of vertex
 * \param ID 
 * \return 
 */
Vertex DataFromFile::getVertexByID(int ID)
{
	auto it = 0, _lastSum = 0;
	auto _sum = seqData[it].getSubstrings().size();
	Vertex _vertex;

	while(ID >= _sum)
	{
		_lastSum = _sum;
		it++;
		_sum += seqData[it].getSubstrings().size();
	}

	_vertex = *seqData[it].getSubstrById(ID - _lastSum);

	return _vertex;
}

/**
 * \brief Return sequence data
 * \return
 */
vector<Sequence> DataFromFile::getSeqData() const
{
	return seqData;
}

/**
 * \brief Return results
 * \return
 */
vector<ResultMotif> DataFromFile::getResults() const
{
	return results;
}

/**
 * \brief Set sequence data from vector
 * \param _seqData
 */
void DataFromFile::setSeqData(vector<Sequence> _seqData)
{
	if(!_seqData.empty())
	{
		seqData = _seqData;
	}
	else
	{
		cout << "Can't assign. _seqData vector is empty." << endl;
	}
}

/**
 * \brief Set number of connected sequences in vertex
 * \param vertex 
 */
void DataFromFile::setVertexConWithOtherSeq(Vertex& vertex)
{
	auto neighbours = vertex.getNeighboursList();
	vector<int> conSequences;

	for(auto& neighbourID : neighbours)
	{
		auto neighbour = getVertexByID(neighbourID);
		if(find(conSequences.begin(), conSequences.end(), neighbour.getSeqIndex()) == conSequences.end())
		{
			conSequences.push_back(neighbour.getSeqIndex());
		}
	}

	auto original = seqData[vertex.getSeqIndex()].getSubstrById(vertex.getIndexInSeq());
	original->setConWithOtherSeq(conSequences.size());
	if(conSequences.size() >= minConnections)
	{
		original->setHasMinConnections(true);
	}
}

/**
 * \brief Add new result to vector results
 * \param result 
 */
void DataFromFile::addResult(ResultMotif result)
{
	/*auto counter = 0;
	for (auto i = 0; i < seqData.size(); i++)
	{
		if (result.getUsedSeqByIndex(i))
		{
			counter++;
		}
	}*/
	if(result.getMotif().size() >= minConnections)// && counter >= minConnections)
	{
		results.push_back(result);
	}
}

/**
 * \brief Build clique based on input vertex vector
 * \param candidates 
 * \return 
 */
vector<Vertex> DataFromFile::buildClique(vector<Vertex> candidates) const
{
	vector<Vertex> clique;
	vector<bool> usedSequences;

	sort(candidates.begin(), candidates.end());

	if(candidates.empty())
	{
		return clique;
	}
	for(auto i = 0; i < seqData.size(); i++)
	{
		usedSequences.push_back(false);
	}

	for(auto& candidate : candidates)
	{
		if(candidate.getVertexLvl() < clique.size())
		{
			break;
		}
		if(clique.empty())
		{
			clique.push_back(candidate);
			usedSequences[candidate.getSeqIndex()] = true;
		}
		else if(!usedSequences[candidate.getSeqIndex()] && checkConnectionsInClique(clique, candidate))
		{
			clique.push_back(candidate);
			usedSequences[candidate.getSeqIndex()];
		}
	}

	if(clique.size() >= minConnections)
	{
		return clique;
	}

	clique.clear();
	return clique;
}

/**
* \brief Create consensus motif for set of vertices
* \param verticesToAlign
* \return
*/
string DataFromFile::buildMotif(vector<Vertex> verticesToAlign)
{
	auto maxValue = 0;
	map<string, int> submotifs = {};
	string consensusMotif;

	sort(verticesToAlign.begin(), verticesToAlign.end(), [](const auto& A, const auto& B)
	     {
		     return A.getIndex() < B.getIndex();
	     });

	for(auto& vertex : verticesToAlign)
	{
		submotifs[vertex.getSubstringWithDeletions()]++;
	}

	for(auto const& motif : submotifs)//foreach motif in motifs
	{
		if(motif.second > maxValue)
		{
			consensusMotif = motif.first;
			maxValue = motif.second;
		}
		else if(motif.second == maxValue && motif.first.length() > consensusMotif.length())
		{
			consensusMotif = motif.first;
			maxValue = motif.second;
		}
	}

	return consensusMotif;
}

/**
 * \brief Build result
 */
ResultMotif DataFromFile::buildResult()
{
	ResultMotif result;
	string generalMotif, motif;
	vector<Vertex> vertexSet, resultSet;
	auto clique = buildClique(vertexByLevel);
	generalMotif = buildMotif(clique);
	result.setStartingClique(clique);
	resultSet = clique;

	auto isInBuild = true;
	do
	{
		vertexSet = prepareVertexSetLeft(clique);
		clique = buildClique(vertexSet);
		if(clique.size() < minConnections)
		{
			break;
		}
		motif = buildMotif(clique);
		generalMotif = parseMotifLeft(generalMotif, motif);
		resultSet = sumResult(resultSet, clique);
	} while(isInBuild);

	clique = result.getStartingClique();
	do
	{
		vertexSet = prepareVertexSetRight(clique);
		clique = buildClique(vertexSet);
		if(clique.size() < minConnections)
		{
			break;
		}
		motif = buildMotif(clique);
		generalMotif = parseMotifRight(generalMotif, motif);
		resultSet = sumResult(resultSet, clique);
	} while(isInBuild);

	result.setMotif(generalMotif);
	result.setResult(resultSet, seqData.size());

	return result;
}

/**
 * \brief Build all results
 */
void DataFromFile::buildResults()
{
	cout << "Generuje rozwiazania, prosze czekac..." << endl;

	vector<Vertex> startingClique;
	vector<bool> usedSequences;

	ResultMotif resultMotif;
	resultMotif.resetUsedSequences(seqData.size());

	//jesli results puste
	for(auto i = 0; i < seqData.size(); i++)
	{
		usedSequences.push_back(false);
	}

	startingClique = buildClique(vertexByLevel);
	addResult(buildResult());

	//jesli ju¿ jest jakieœ rozwi¹zanie usun pierwotna klike z wektora startowego
	auto i = 1;
	while(i < ITERATIONS)
	{
		vertexByLevel = filterVector(vertexByLevel, startingClique);
		if(!vertexByLevel.empty())
		{
			startingClique = buildClique(vertexByLevel);
			resultMotif = buildResult();
			if(find(results.begin(), results.end(), resultMotif) == results.end())
			{
				addResult(resultMotif);
				i++;
			}
		}
		else
		{
			i = ITERATIONS;
		}
	}

	sortResults(results);
}

/**
 * \brief Check if vertex can be add to clique
 * \param clique 
 * \param toCheck 
 * \return 
 */
bool DataFromFile::checkConnectionsInClique(vector<Vertex> clique, Vertex toCheck)
{
	for(auto& vertex : clique)
	{
		auto neighbours = vertex.getNeighboursList();
		if(find(neighbours.begin(), neighbours.end(), toCheck.getIndex()) == neighbours.end())
		{
			return false;
		}
	}
	return true;
}

/**
* \brief Create neighbours lists for all substrings
*/
void DataFromFile::createEdges()
{
	auto createdEdges = 0;

	cout << "Tworze polaczenia, prosze czekac..." << endl;

	for(auto i = 0; i < seqData.size() - 1; i++)
	{
		for(auto j = i + 1; j < seqData.size(); j++)
		{
			for(auto& substring1 : seqData[i].getSubstrings())
			{
				auto substr1 = substring1.getSubstringWithDeletions();

				if(substr1.size() < substrLength - PERMITTED_DELETIONS)
				{
					continue;
				}
				for(auto& substring2 : seqData[j].getSubstrings())
				{
					auto substr2 = substring2.getSubstringWithDeletions();

					if(substr2.size() >= substrLength - PERMITTED_DELETIONS)
					{
						auto vertex1 = seqData[i].getSubstrById(substring1.getIndexInSeq());
						auto vertex2 = seqData[j].getSubstrById(substring2.getIndexInSeq());

						if(substr1.size() >= substr2.size() && substr1.find(substr2) != string::npos)
						{
							vertex1->addToNeighboursList(vertex2->getIndex());
							vertex2->addToNeighboursList(vertex1->getIndex());
							seqData[i].getSubstrById(substring1.getIndexInSeq())->lvlUp();
							seqData[j].getSubstrById(substring2.getIndexInSeq())->lvlUp();
							createdEdges++;
						}
						else if(substr2.find(substr1) != string::npos)
						{
							vertex1->addToNeighboursList(vertex2->getIndex());
							vertex2->addToNeighboursList(vertex1->getIndex());
							seqData[i].getSubstrById(substring1.getIndexInSeq())->lvlUp();
							seqData[j].getSubstrById(substring2.getIndexInSeq())->lvlUp();
							createdEdges++;
						}
					}
				}
			}
		}
	}

	if(createdEdges < 1)
	{
		cout << "Nie utworzono zadnych polaczen miedzy podciagami w sekwencjach." << endl;
		exit(0);
	}
	cout << "Utworzono polaczen: " << createdEdges << endl;
}

/**
 * \brief Create graph of substrings from loaded data
 */
void DataFromFile::createGraph()
{
	seqData = sequencesIntoSubstrings(seqData);
}

/**
 * \brief Create sorted list of vertices
 */
void DataFromFile::createListOfVerticesSorted()
{
	cout << "Generuje liste wierzcholkow do startu budowy rozwiazania, prosze czekac..." << endl;
	for(auto i = 0; i < seqData.size(); i++)
	{
		for(auto& substring : seqData[i].getSubstrings())
		{
			setVertexConWithOtherSeq(substring);

			if(substring.getSubstringWithDeletions().size() >= substrLength - PERMITTED_DELETIONS)
			{
				vertexByLevel.push_back(substring);
			}
		}
	}

	//sort list
	sort(vertexByLevel.begin(), vertexByLevel.end());
	cout << "Pomyslnie wygenerowano liste wierzcholkow." << endl;
}

/**
 * \brief Filter vector toFilter from elements from vector filtering
 * \param toFilter 
 * \param filtering 
 * \return 
 */
vector<Vertex> DataFromFile::filterVector(vector<Vertex> toFilter, vector<Vertex> filtering)
{
	vector<Vertex> result;
	for(auto& elementToFilter : toFilter)
	{
		if(find(filtering.begin(), filtering.end(), elementToFilter) == filtering.end())
		{
			result.push_back(elementToFilter);
		}
	}

	return result;
}

/**
* \brief Load data from fasta file
* \param _dataName
* \param _seqData
*/
void DataFromFile::loadFromFile(string _dataName, vector<Sequence> _seqData)
{
	fstream file;
	string line;
	string lastSeq = "";
	auto path = _dataName + ".txt";

	file.open(path);
	if(!file.good())
	{
		cout << "Error! Can't open the file." << endl;
		exit(0);
	}

	Sequence s1;
	auto seqId = 0;

	while(getline(file, line))
	{
		if(line[0] == '>')
		{
			if(lastSeq != "")
			{
				s1.setSequence(lastSeq);
				_seqData.push_back(s1);
				lastSeq = "";
			}
			s1.setName(line);
			s1.setSeqId(seqId);
			seqId++;
		}
		else
		{
			if(lastSeq == "")
			{
				lastSeq = line;
			}
			else
			{
				lastSeq += line;
			}
		}
	}

	s1.setSequence(lastSeq);
	_seqData.push_back(s1);
	lastSeq = "";

	if(!_seqData.empty())
	{
		seqData = _seqData;
		cout << "Pomyslnie wczytano plik fasta" << endl;
	}
	else
	{
		cout << "Can't load fasta file." << endl;
	}
}

/**
* \brief Load data from qual file
* \param _dataName
* \param _seqData
*/
void DataFromFile::loadQualFile(string _dataName, vector<Sequence> _seqData)
{
	fstream file;
	auto seqNo = 0;
	string line;
	vector<int> qualities;

	auto path = _dataName + "-qual.txt";

	file.open(path);
	if(!file.good())
	{
		cout << "Error! Can't open the file." << endl;
		exit(0);
	}

	while(getline(file, line))
	{
		if(line[0] == '>')
		{
			if(line != _seqData[seqNo].getName())
			{
				cout << "Error! Fasta file and qual file not equal." << endl;
				exit(0);
			}
			if(!qualities.empty())
			{
				_seqData[seqNo - 1].setQual(qualities);
				qualities.clear();
			}
			seqNo++;
		}
		else
		{
			stringstream iss(line);
			int number;

			while(iss >> number)
			{
				qualities.push_back(number);
			}
		}
	}

	_seqData[seqNo - 1].setQual(qualities);

	if(!_seqData.empty())
	{
		seqData = _seqData;
		cout << "Pomyslnie wczytano plik qual" << endl;
	}
	else
	{
		cout << "Can't load qual file." << endl;
	}
}

/**
 * \brief Parse general motif left with extension
 * \param genMotif 
 * \param toParse 
 * \return 
 */
string DataFromFile::parseMotifLeft(string genMotif, string toParse)
{
	string toAdd = "", motif;
	motif = genMotif;
	if(toParse.size() > genMotif.size())
	{
		if(toParse.find(genMotif) != string::npos)
		{
			return toParse;
		}
		return (toParse + genMotif);
	}

	auto diff = motif.size() - toParse.size();
	motif = motif.substr(0, motif.size() - diff);

	while(motif != toParse && toParse != "")
	{
		toAdd += toParse[0];
		toParse = toParse.substr(1, toParse.size());
		motif.pop_back();
	}
	if(toParse.empty())
	{
		return (toParse + genMotif);
	}
	return (toAdd + genMotif);
}

/**
* \brief Parse general motif with right extension
* \param genMotif
* \param toParse
* \return
*/
string DataFromFile::parseMotifRight(string genMotif, string toParse)
{
	string toAdd = "", motif;
	motif = genMotif;
	if(toParse.size() > genMotif.size())
	{
		if(toParse.find(genMotif) != string::npos)
		{
			return toParse;
		}
		return (genMotif + toParse);
	}

	auto diff = motif.size() - toParse.size();
	motif = motif.substr(diff, motif.size());

	while(motif != toParse && toParse != "")
	{
		toAdd = toParse[toParse.size() - 1] + toAdd;
		toParse.pop_back();
		motif = motif.substr(1, motif.size());
	}
	return (genMotif + toAdd);
}

/**
 * \brief Create set of vertices to extend motif on left
 * \param actualResult 
 * \return 
 */
vector<Vertex> DataFromFile::prepareVertexSetLeft(vector<Vertex> actualResult)
{
	auto resultSortedByIndex = actualResult;
	vector<Vertex> vertexSet;

	sort(resultSortedByIndex.begin(), resultSortedByIndex.end(), [](const auto& a, const auto& b)
	     {
		     return a.getIndex() < b.getIndex();
	     });

	for(auto& vertex : actualResult)
	{
		for(auto i = 1; i <= SENSITIVITY; i++)
		{
			auto index = vertex.getIndex() - i;
			if(index >= 0
			   && getVertexByID(index).getHasMinConnections()
			   && getVertexByID(index).getSeqIndex() == vertex.getSeqIndex())
			{
				vertexSet.push_back(getVertexByID(index));
			}
		}
	}

	return vertexSet;
}

/**
* \brief Create set of vertices to extend motif on right
* \param actualResult
* \return
*/
vector<Vertex> DataFromFile::prepareVertexSetRight(vector<Vertex> actualResult)
{
	auto resultSortedByIndex = actualResult;
	vector<Vertex> vertexSet;

	sort(resultSortedByIndex.begin(), resultSortedByIndex.end(), [](const auto& a, const auto& b)
	     {
		     return a.getIndex() > b.getIndex();
	     });

	auto size = 0;
	for(auto j = 0; j < seqData.size(); j++)
	{
		size += seqData[j].getSubstrings().size();
	}

	for(auto& vertex : actualResult)
	{
		for(auto i = 1; i <= SENSITIVITY; i++)
		{
			auto index = vertex.getIndex() + i;
			if(index < size
			   && getVertexByID(index).getHasMinConnections()
			   && getVertexByID(index).getSeqIndex() == vertex.getSeqIndex())
			{
				vertexSet.push_back(getVertexByID(index));
			}
		}
	}

	return vertexSet;
}

/**
 * \brief Print sequences
 * \param seqData 
 * \param resultSeq 
 */
void DataFromFile::printSequences(vector<Sequence> seqData, vector<Sequence> resultSeq)
{
	ResultMotif resultMotif;
	for(auto i = 0; i < seqData.size(); i++)
	{ //for each sequence
		int seqSize = seqData[i].getSequence().size(); //get length of sequence
		string sequence(seqSize, '-'); //create string with '-' chars with length of sequence
		auto substrs = seqData[i].getSubstrings();

		for(auto j = 0; j < seqData[i].getSubstrings().size(); j++)
		{ //for each substring in sequence
			if(substrs[j].getHasMinConnections())
			{
				for(auto k = 0; k < substrs[j].getSubstring().length(); k++)
				{ //for each char in substring
					if(substrs[j].getQual()[k] >= reliability)
					{
						sequence[j + k] = substrs[j].getSubstring()[k];
					}
					else
					{
						sequence[j + k] = '*';
					}
				}
			}
		}

		cout << "\n" << seqData[i].getName() << endl; //print name of sequence
		cout << "O: " << DataFromFile::seqData[i].getSequence() << endl;
		cout << "D: " << sequence << endl;

		//print motifs
		for(auto j = 0; j < resultSeq.size(); j++)
		{
			if(seqData[i].getSeqId() == resultSeq[j].getSeqId())
			{
				auto resultVertices = resultSeq[j].getSubstrings();
				resultMotif.printMotifOnSeq(resultVertices, seqSize, reliability);
				break;
			}
			if(j == resultSeq.size() - 1)
			{
				string seqMotifEmpty(seqSize, '-');
				cout << "M: " << seqMotifEmpty << '\n' << endl;
			}
		}
	}

	//print vertices
	cout << "Wierzcholki biorace udzial w budowie rozwiazania (id_sekwencji:id_podciagu): " << endl;
	for(auto i = 0; i < resultSeq.size(); i++)
	{
		resultMotif.printVerticesInMotif(resultSeq[i].getSubstrings(), reliability);
	}
}

/**
 * \brief Print all results
 */
void DataFromFile::printResult()
{
	auto longestMotif = 0;
	if(!results.empty())
	{
		for(auto i = 0; i < results.size(); i++)
		{
			if(i == 0)
			{
				longestMotif = results[0].getMotif().size();
			}

			if(longestMotif == results[i].getMotif().size())
			{
				cout << "Znaleziono motyw: " << results[i].getMotif() << endl;
				printSequences(seqData, results[i].getResult());
			}
			else
			{
				break;
			}
		}
	}
	else
	{
		cout << "\nNie znaleziono motywu w podanych sekwencjach." << endl;
	}
}

/**
* \brief Transform sequences data into substrings for each sequence
* \param _data
* \return
*/
vector<Sequence> DataFromFile::sequencesIntoSubstrings(vector<Sequence> _data) const
{
	auto startingIndex = 0;
	for(auto i = 0; i < _data.size(); i++)
	{
		_data[i].createSubstrings(_data[i].getSequence(), _data[i].getQual(),
		                          i, startingIndex, substrLength, reliability);
		startingIndex += _data[i].getSubstrings().size();
	}

	return _data;
}

/**
 * \brief Sort results
 * \param results 
 */
void DataFromFile::sortResults(vector<ResultMotif>& results) const
{
	sort(results.begin(), results.end());
	reverse(results.begin(), results.end());
}

/**
* \brief Sum two result vectors
* \param actualResult
* \param tempResult
* \return
*/
vector<Vertex> DataFromFile::sumResult(vector<Vertex> actualResult, vector<Vertex> tempResult)
{
	for(auto i = 0; i < tempResult.size(); i++)
	{
		auto found = (find(actualResult.begin(), actualResult.end(), tempResult[i]) != actualResult.end());

		if(!found)
		{
			actualResult.push_back(tempResult[i]);
		}
	}
	return actualResult;
}

/**
 * \brief Public constructor for object with data loaded from file
 * \param _dataName
 * \param _substrLength
 * \param _reliability
 */
DataFromFile::DataFromFile(string _dataName, int _substrLength, int _reliability)
{
	//load users parameters
	dataName = _dataName;
	substrLength = _substrLength;
	reliability = _reliability;

	//load data from fasta and qual file
	loadFromFile(dataName, seqData);
	loadQualFile(dataName, seqData);

	//set minConnections parameter
	minConnections = seqData.size() * SEQUENCES_PERCENT;
}

/**
 * \brief Destructor for object with data loaded from file
 */
DataFromFile::~DataFromFile()
{
	dataName.clear();
	seqData.clear();
	vertexByLevel.clear();
	results.clear();
}
