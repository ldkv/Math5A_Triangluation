#pragma once

#include <vector>

using namespace std;

static int globalId = 0;
static int globalSideId = 0;
static int globalFaceId = 0;

struct Point
{
	int id;
	QVector3D coord;
	vector<int> sides;
	Point() 
	{
		id = globalId++;
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

	vector<Point> points;

	Side(int h, int l)
	{
		pLow = l;
		pHigh = h;
		id = globalSideId++;
	}
	Side(Point p1, Point p2) {
		points.push_back(p1);
		points.push_back(p2);
		id = globalSideId++;
	}
	bool operator==(const Side& e)
	{
		if (points[0].coord == e.points[0].coord && points[1].coord == e.points[1].coord)
			return true;

		return false;
	}
	bool operator!=(const Side& e) { return !(this->operator==(e)); }
};

struct Face
{
	int id;
	vector<Side> sides;
	vector<Point> points;

	/*Face(int sideId) {
		sides.push_back(sideId);
		id = globalFaceId++;
	}
	Face(int sideId1, int sideId2, int sideiD3) {
		sides.push_back(sideId1);
		sides.push_back(sideId2);
		sides.push_back(sideiD3);
		id = globalFaceId++;
	}*/
	Face(Point p1, Point p2, Point p3) {
		points.push_back(p1);
		points.push_back(p2);
		points.push_back(p3);
		sides.push_back(Side(p1, p2));
		sides.push_back(Side(p2, p3));
		sides.push_back(Side(p3, p1));
		id = globalFaceId++;
	}
};

Point *getPointfromID(vector<Point> pts, int id);
int getPointIndex(vector<Point> pts, int id);
int getSideIDFromPoints(vector<Side> s, Point x, Point y);
int getSideFromID(vector<Side> s, int id);
vector<Side> FindExternSides();
//bool checkVisibilitySide(Side side, Point p);
vector<QVector3D> EnvelopeJarvis(vector<Point> pts);
vector<Face> TriangulationSimple(vector<Point> pts);
int getPointIndex(vector<Point> pts, int id);
//vector<Side> GrahamScan(vector<Point> pts);
vector<Point> GrahamScan(vector<Point> pts);
vector<Side> getViewedEdge(int nextIdVert, vector<Point> pts, list<Side> &convexHull);
bool isEdgeViewed(QVector3D P, QVector3D A, QVector3D B, QVector3D n);
QVector3D crossProductNormalized(QVector3D p, QVector3D op);
