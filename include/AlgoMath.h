#pragma once

#include <vector>

using namespace std;

static int globalId = 0;
static int globalSideId = 0;
static int globalFaceId = 0;


#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>

#define MAXN 1010

typedef long long vtype;

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

//static list<Side>  _convexHull;

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

// VORONOI
vector<Point> diagramVoronoi(vector<Point> pts, vector<Side> sides, vector<Face> faces);

// DELAUNAY
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
QVector3D circumCircleCenter(QVector3D v1, QVector3D v2, QVector3D v3);
vector<Side> getViewedEdges(vector<Side> sides, vector<Point> pts, Point P); // methode utilisant les id
void Delaunay_addPoint(vector<Point> &pts, vector<Side> &sides, vector<Face> &faces, QVector3D P);
void Delaunay_deletePoint(vector<Point> &pts, vector<Side> &sides, vector<Face> &faces, int deleteID);

//===========================================================================================================
int getPointIndex(vector<Point> pts, int id);
int getSideIDFromPoints(vector<Side> s, Point x, Point y);
//int getSideFromID(vector<Side> s, int id);
//vector<Side> FindExternSides();
//bool checkVisibilitySide(Side side, Point p);
vector<Face> TriangulationSimple(vector<Point> pts, list<Side> &_convexHull);

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


/* Basic 3D vector implementation */
struct vec3 {
	vec3() { X[0] = X[1] = X[2] = 0; }
	vec3(vtype x, vtype y, vtype z) { X[0] = x; X[1] = y; X[2] = z; }

	/* 3D cross product */
	vec3 operator*(const vec3& v) const {
		return vec3(X[1] * v.X[2] - X[2] * v.X[1],
			X[2] * v.X[0] - X[0] * v.X[2],
			X[0] * v.X[1] - X[1] * v.X[0]);
	}

	vec3 operator-(const vec3& v) const {
		return vec3(X[0] - v.X[0], X[1] - v.X[1], X[2] - v.X[2]);
	}

	vec3 operator-() const {
		return vec3(-X[0], -X[1], -X[2]);
	}

	vtype dot(const vec3& v) const {
		return X[0] * v.X[0] + X[1] * v.X[1] + X[2] * v.X[2];
	}

	vtype X[3];
};


struct face {
	vec3 norm;
	vtype disc;
	int I[3];
};

vector<face> convexHull3D(vector<Point> pts);


static vector<QVector3D> fpoints;

struct Face2
{
	int i0, i1, i2;
	float a, b, c, d;

	Face2(int i0o, int i1o, int i2o)
	{

		i0 = i0o;
		i1 = i1o;
		i2 = i2o;

		computePlane();

	}

	void computePlane()
	{
		QVector3D v1 = fpoints[i0];
		QVector3D v2 = fpoints[i1];
		QVector3D v3 = fpoints[i2];

		a = v1.y() * (v2.z() - v3.z()) + v2.y() * (v3.z() - v1.z()) + v3.y() * (v1.z() - v2.z());
		b = v1.z() * (v2.x() - v3.x()) + v2.z() * (v3.x() - v1.x()) + v3.z() * (v1.x() - v2.x());
		c = v1.x() * (v2.y() - v3.y()) + v2.x() * (v3.y() - v1.y()) + v3.x() * (v1.y() - v2.y());
		d = -(v1.x() * (v2.y() * v3.z() - v3.y() * v2.z()) + v2.x() * (v3.y() * v1.z() - v1.y() * v3.z()) + v3.x() * (v1.y() * v2.z() - v2.y() * v1.z()));
	}

	bool isVisible(QVector3D p)
	{

		return (a * p.x() + b * p.y() + c * p.z() + d) > 0;

	}

	QVector3D centroid()
	{

		QVector3D p0 = fpoints[i0];
		QVector3D p1 = fpoints[i1];
		QVector3D p2 = fpoints[i2];
		return QVector3D((p0.x() + p1.x() + p2.x()) / 3, (p0.y() + p1.y() + p2.y()) / 3, (p0.z() + p1.z() + p2.z()) / 3);

	}

	void flip()
	{
		int t = i0;
		i0 = i1;
		i1 = t;
		computePlane();
	}

};

struct ConvexHull
{
	vector<Face2> validFaces;
	vector<Face2>  visibleFaces;
	vector<Face2>  tmpFaces;

	QVector3D centroid(vector<QVector3D> points, int index, Face2 face)
	{

		QVector3D p = points[index];
		QVector3D p0 = points[face.i0];
		QVector3D p1 = points[face.i1];
		QVector3D p2 = points[face.i2];
		return QVector3D((p.x() + p0.x() + p1.x() + p2.x()) / 4, (p.y() + p0.y() + p1.y() + p2.y()) / 4, (p.z() + p0.z() + p1.z() + p2.z()) / 4);

	}

	vector<int> process(vector<QVector3D> points)
	{
		vector<int> result;

		if (points.size() < 4)
		{
			return result;
		}

		//local copy of the point set
		//vertices = .concat();
		//Face2.points = points;
		fpoints = points;


		//calculates the first convex tetrahedron

		//creates a face with the first 3 vertices
		Face2 face(0, 1, 2);

		//this is the center of the tetrahedron, all face should point outwards:
		//they should not be visible to the centroid
		QVector3D v = centroid(points, 3, face);

		if (face.isVisible(v)) face.flip();

		Face2 face0 = Face2(3, face.i0, face.i1);
		if (face0.isVisible(v)) face0.flip();

		Face2 face1 = Face2(3, face.i1, face.i2);
		if (face1.isVisible(v)) face1.flip();

		Face2 face2 = Face2(3, face.i2, face.i0);
		if (face2.isVisible(v)) face2.flip();


		//store the tetrahedron faces in the valid faces list
		vector<Face2> validFaces;
		validFaces.push_back(face);
		validFaces.push_back(face0);
		validFaces.push_back(face1);
		validFaces.push_back(face2);

		visibleFaces.clear();
		tmpFaces.clear();



		//so as we have a convex tetrahedron, we can skip the first 4 points
		for (int i = 4; i < points.size(); i++)
		{
			//for each avaiable vertices
			v = points[i];

			//checks the point's visibility from all faces
			visibleFaces.clear();
			for each(face in validFaces)
			{
				if (face.isVisible(v))
				{
					visibleFaces.push_back(face);
				}
			}

			//the vertex is not visible : it is inside the convex hull, keep on
			if (visibleFaces.size() == 0)
			{
				continue;
			}

			//the vertex is outside the convex hull
			//delete all visible faces from the valid List
			for each (face in visibleFaces)
			{
				//int pos = find(validFaces.begin(), validFaces.end(), face) - validFaces.begin();
				//vector<int> index;
				int index;
				for (int y = 0; y < validFaces.size(); y++)
				{
					if (validFaces[y].i0 == face.i0 && validFaces[y].i1 == face.i1 && validFaces[y].i2 == face.i2) {
						//index.push_back(y);
						index = y;
						break;
					}
				}
				/*int offeset = 0;
				for (int y = 0; y < index.size(); y++)
				{
					validFaces.erase(validFaces.begin() + y- offeset);
					offeset++;
				}*/
				validFaces.erase(validFaces.begin() + index);
				//validFaces.splice(validFaces.indexOf(face), 1);
			}

			//special case : only one face is visible
			//it's ok to create 3 faces directly for they won't enclose any other point
			if (visibleFaces.size() == 1)
			{
				face = visibleFaces[0];
				validFaces.push_back(Face2(i, face.i0, face.i1));
				validFaces.push_back(Face2(i, face.i1, face.i2));
				validFaces.push_back(Face2(i, face.i2, face.i0));
				continue;
			}

			//creates all possible new faces from the visibleFaces
			tmpFaces.clear();
			for each(face in visibleFaces)
			{
				tmpFaces.push_back(Face2(i, face.i0, face.i1));
				tmpFaces.push_back(Face2(i, face.i1, face.i2));
				tmpFaces.push_back(Face2(i, face.i2, face.i0));
			}

			//Face2 other(0, 0, 0);
			for each(face in tmpFaces)
			{
				//search if there is a point in front of the face : 
				//this means the face doesn't belong to the convex hull
				bool search = false;
				for each(Face2 other in tmpFaces)
				{
					if (face.i0 != other.i0 || face.i1 != other.i1 || face.i2 != other.i2)
					{
						if (face.isVisible(other.centroid()))
						{
							search = true;
							break;
						}
					}
				}
				//the face has no point in front of it
				if (!search) validFaces.push_back(face);
			}
		}

		for each(face in validFaces)
		{
			result.push_back(face.i0);
			result.push_back(face.i1);
			result.push_back(face.i2);
		}
		return result;
	}
};