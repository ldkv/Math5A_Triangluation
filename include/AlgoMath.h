#pragma once

#include <vector>

using namespace std;

static int globalId_Points = 0;
static int globalId_Sides = 0;
static int globalId_Faces = 0;
struct Point
{
	int id;
	QVector3D coord;
	vector<int> sides;
	Point() 
	{
		id = globalId_Points++;
	};
	Point(int idi, QVector3D pt, vector<int> s)
	{
		id = idi;
		coord = pt;
		sides = s;
	}
	Point(QVector3D pt)
	{
		coord = pt;
		id = globalId_Points++;
	}
};

struct Side
{
	int id;
	int pLow;
	int pHigh;
	int fLeft;
	int fRight;
	Side()
	{
		id = globalId_Sides++;
	};
	Side(int l, int h)
	{
		pLow = l;
		pHigh = h;
		id = globalId_Sides++;
	}
};

struct Face
{
	int id;
	vector<int> sides;
};

Point *getPointfromID(vector<Point> pts, int id);
vector<Point> EnveloppeJarvis(vector<Point> pts);
vector<Side> TriangulationSimple(vector<Point> pts);
int getPointIndex(vector<Point> pts, int id);
//vector<Side> GrahamScan(vector<Point> pts);
vector<Point> GrahamScan(vector<Point> pts);
void Delaunay_addPoint(vector<Point> &pts, vector<Side> &sides, vector<Face> &faces, QVector3D P);