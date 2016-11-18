#include "stdafx.h"
#include "VertexInList.h"


Vertex VertexInList::getVertex()
{
	return VertexInList::vertex;
}


VertexInList::VertexInList(Vertex _vertex, int _index, int _seqIndex)
{
	VertexInList::vertex = _vertex;
	VertexInList::index = _index;
	VertexInList::seqIndex = _seqIndex;
}

VertexInList::~VertexInList()
{
}
