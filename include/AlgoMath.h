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
		id = -1;
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
	int idFace1 = -1;
	int idFace2 = -1;

	Side()
	{
		pLow = -1;
		pHigh = -1;
		fLeft = -1;
		fRight = -1;
		id = -1;
	}

	Side(int l, int h)
	{
		pLow = l;
		pHigh = h;
		fLeft = -1;
		fRight = -1;
		id = globalSideId++;
	}
	Side(Point p1, Point p2) {
		points.clear();
		points.push_back(p1);
		points.push_back(p2);
		pLow = p1.id;
		pHigh = p2.id;
		fLeft = -1;
		fRight = -1;
		id = globalSideId++;
	}
	Side(Point p1, Point p2, int idFace) {
		points.clear();
		points.push_back(p1);
		points.push_back(p2);
		idFace1 = idFace;
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
	vector<int> sidesID;
	vector<Point> points;

	Face()
	{
		id = -1;
	}
	Face(int sideId) {
		sidesID.push_back(sideId);
		id = globalFaceId++;
	}
	Face(int sideId1, int sideId2, int sideiD3) {
		sidesID.push_back(sideId1);
		sidesID.push_back(sideId2);
		sidesID.push_back(sideiD3);
		id = globalFaceId++;
	}
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

void resetGlobalID();
Point getPointfromID(vector<Point> pts, int id);
Side getSidefromID(vector<Side> sds, int id);
Face getFacefromID(vector<Face> faces, int id);
int getSideIDfromPoints(vector<Side> sides, int p1, int p2);
void deleteSidefromID(int id, vector<Side> &sds, vector<Point> &pts);
void deleteSidefromPoint(vector<Point> &pts, int sideID, int pointID);
void deleteFacefromID(int id, vector<Face> &faces, vector<Side> &sides);
void deleteFacefromSide(vector<Side> &sides, int sideID, int faceID);
vector<Point> getVertexesfromFace(Face F, vector<Point> pts, vector<Side> sides);
int addSide(int p1, int p2, vector<Side> &sides, vector<Point> &pts);
int addFace(int id1, int id2, int id3, vector<Face> &faces, vector<Side> &sides);
void removeSidefromPoint(vector<Point> &pts, int ptID, int sideID);
void addSidetoPoint(vector<Point> &pts, int ptID, int sideID);
bool coordsSort(Point i, Point j);
vector<Side> getIncidentEdgesOriented(Point P, vector<Side> sides, vector<Point> pts);
vector<Face> getIncidentFacesOriented(vector<Side> La1, vector<Face> faces, vector<Side> sides, vector<Side> &La2);
int mutualPoint(Side S1, Side S2);
int findPointConvexDelaunay(bool closedPoly, vector<Side> La2, vector<Point> pts, QVector3D P, Point &S1, Point &S2);


vector<Point> EnvelopeJarvis(vector<Point> pts);
bool Collinear(QVector3D v1, QVector3D v2);
bool insideTriangle(QVector3D pt, QVector3D v1, QVector3D v2, QVector3D v3);
bool insideCircumCircle(QVector3D pt, QVector3D v1, QVector3D v2, QVector3D v3);
vector<Side> getViewedEdges(vector<Side> sides, vector<Point> pts, Point P); // methode utilisant les id
void Delaunay_addPoint(vector<Point> &pts, vector<Side> &sides, vector<Face> &faces, QVector3D P);
void Delaunay_deletePoint(vector<Point> &pts, vector<Side> &sides, vector<Face> &faces, int deleteID);

//===========================================================================================================
int getPointIndex(vector<Point> pts, int id);
int getSideIDFromPoints(vector<Side> s, Point x, Point y);
//int getSideFromID(vector<Side> s, int id);
//vector<Side> FindExternSides();
//bool checkVisibilitySide(Side side, Point p);
vector<Face> TriangulationSimple(vector<Point> pts);

int getPointIndex(vector<Point> pts, int id);
//vector<Side> GrahamScan(vector<Point> pts);
vector<Point> GrahamScan(vector<Point> pts);
vector<Side> getViewedEdge(int nextIdVert, vector<Point> pts, list<Side> &convexHull);

bool isEdgeViewed(QVector3D P, QVector3D A, QVector3D B, QVector3D n);
QVector3D crossProductNormalized(QVector3D p, QVector3D op);
vector<Side> Fliping(vector<Face> faces);
vector<Face> Fliping2(vector<Face> faces);
bool inCircumCircle(Face f, QVector3D v);
vector<Point> Voronoi(vector<Point> pts);
QVector3D CircumCircleCenter(Face f);
