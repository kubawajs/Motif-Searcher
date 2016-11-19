#pragma once

#include "VertexInList.h"
#include "Vertex.h"

using namespace std;

class VertexInList
{
	Vertex vertex;
	int index;
	int seqIndex;
	//TODO: add flag if Vertex is in result??

public:
	Vertex getVertex();
	int getIndex();
	int getSeqIndex();
	VertexInList(Vertex _vertex, int _index, int _seqIndex);
	~VertexInList();
};
