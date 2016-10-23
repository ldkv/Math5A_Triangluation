#pragma once

#include <vector>

using namespace std;

static int globalId = 0;

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
		id = globalId++;
	}
};

struct Side
{
	int id;
	int pLow;
	int pHigh;
	int fLeft;
	int fRight;

	Side(int h, int l)
	{
		pLow = l;
		pHigh = h;
	}
};

struct Face
{
	int id;
	vector<int> sides;
};

Point *getPointfromID(vector<Point> pts, int id);
vector<QVector3D> EnvelopeJarvis(vector<Point> pts);
vector<Side> TriangulationSimple(vector<Point> pts);
int getPointIndex(vector<Point> pts, int id);