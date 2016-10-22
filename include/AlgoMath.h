#pragma once

#include <vector>

using namespace std;

struct Point
{
	int id;
	QVector3D coord;
	vector<int> sides;
	Point(int idi, QVector3D pt, vector<int> s)
	{
		id = idi;
		coord = pt;
		sides = s;
	}
	Point(QVector3D pt)
	{
		coord = pt;
	}
};

struct Side
{
	int id;
	int pLow;
	int pHigh;
	int fLeft;
	int fRight;
};

struct Face
{
	int id;
	vector<int> sides;
};