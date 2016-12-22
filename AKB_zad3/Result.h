#pragma once
#include "DataFromFile.h"

class Result :
	public DataFromFile
{
	vector <Sequence> result;

public:
	Result();
	~Result();
};

