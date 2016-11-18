#pragma once

#include "VertexInList.h"
#include "Vertex.h"

using namespace std;

class VertexInList
{
	Vertex vertex;
	int index;
	int seqIndex;

public:
	Vertex getVertex();
	VertexInList(Vertex _vertex, int _index, int _seqIndex);
	~VertexInList();
};
