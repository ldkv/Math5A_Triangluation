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
vector<Side> getViewedEdges(vector<Side> sides, vector<Point> pts, Point P); // methode utilisant les id
void Delaunay_addPoint(vector<Point> &pts, vector<Side> &sides, vector<Face> &faces, QVector3D P);
void Delaunay_deletePoint(vector<Point> &pts, vector<Side> &sides, vector<Face> &faces, int deleteID);

vector<Point> EnvelopeJarvis(vector<Point> pts);
bool Collinear(QVector3D v1, QVector3D v2);
bool insideTriangle(QVector3D pt, QVector3D v1, QVector3D v2, QVector3D v3);
bool insideCircumCircle(QVector3D pt, QVector3D v1, QVector3D v2, QVector3D v3);
QVector3D circumCircleCenter(QVector3D v1, QVector3D v2, QVector3D v3);

//--------------------- TriangulationSimple ------------------------------------------------------------
int getPointIndex(vector<Point> pts, int id);
int getSideIDFromPoints(vector<Side> s, Point x, Point y);
vector<Face> TriangulationSimple(vector<Point> pts, list<Side> &_convexHull);

//--------------------- Graham Scan --------------------------------------------------------
int getPointIndex(vector<Point> pts, int id);
vector<Point> GrahamScan(vector<Point> pts);
vector<Side> getViewedEdge(int nextIdVert, vector<Point> pts, list<Side> &convexHull);

//--------------------- Flipping et Voronoi --------------------------------------------------------
bool isEdgeViewed(QVector3D P, QVector3D A, QVector3D B, QVector3D n);
QVector3D crossProductNormalized(QVector3D p, QVector3D op);
vector<Side> Flipping(vector<Face> faces);
vector<Face> Flipping2(vector<Face> faces);
bool inCircumCircle(Face f, QVector3D v);
vector<Point> Voronoi(vector<Point> pts);
QVector3D CircumCircleCenter(Face f);

//---------------------------------------- Enveloppe convexe 3D -------------------------------------------------------
static vector<QVector3D> cloudPoints;

struct Face3D
{
	int index0, index1, index2;
	float a, b, c, d;

	Face3D(int id0, int id1, int id2)
	{
		index0 = id0;
		index1 = id1;
		index2 = id2;

		recalculFace();
	}

	void recalculFace()
	{
		QVector3D v1 = cloudPoints[index0];
		QVector3D v2 = cloudPoints[index1];
		QVector3D v3 = cloudPoints[index2];

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
		QVector3D p0 = cloudPoints[index0];
		QVector3D p1 = cloudPoints[index1];
		QVector3D p2 = cloudPoints[index2];
		return QVector3D((p0.x() + p1.x() + p2.x()) / 3, (p0.y() + p1.y() + p2.y()) / 3, (p0.z() + p1.z() + p2.z()) / 3);
	}

	void flip()
	{
		int t = index0;
		index0 = index1;
		index1 = t;
		recalculFace();
	}

};

struct ConvexHull
{
	// Structures de base
	vector<Face3D> validFaces;
	vector<Face3D>  viewFaces;
	vector<Face3D>  tmpFaces;

	// Calcul du centre des points
	QVector3D centroid(vector<QVector3D> points, int index, Face3D face)
	{
		QVector3D p = points[index];
		QVector3D p0 = points[face.index0];
		QVector3D p1 = points[face.index1];
		QVector3D p2 = points[face.index2];
		return QVector3D((p.x() + p0.x() + p1.x() + p2.x()) / 4, (p.y() + p0.y() + p1.y() + p2.y()) / 4, (p.z() + p0.z() + p1.z() + p2.z()) / 4);
	}

	vector<int> calculate(vector<QVector3D> points)
	{
		vector<int> result;
		if (points.size() < 4)
		{
			return result;
		}

		// Référence pour les faces
		cloudPoints = points;

		// Calcul du tetrahedron
		Face3D face(0, 1, 2);
		QVector3D v = centroid(points, 3, face);
		if (face.isVisible(v)) face.flip();
		Face3D face0 = Face3D(3, face.index0, face.index1);
		if (face0.isVisible(v)) face0.flip();
		Face3D face1 = Face3D(3, face.index1, face.index2);
		if (face1.isVisible(v)) face1.flip();
		Face3D face2 = Face3D(3, face.index2, face.index0);
		if (face2.isVisible(v)) face2.flip();


		// Enregistrement des faces du tetrahedron comme faces valides
		vector<Face3D> trueFaces;
		trueFaces.push_back(face);
		trueFaces.push_back(face0);
		trueFaces.push_back(face1);
		trueFaces.push_back(face2);

		viewFaces.clear();
		tmpFaces.clear();

		// Parcourt des points
		for (int i = 4; i < points.size(); i++)
		{
			v = points[i];
			viewFaces.clear();
			for each(face in trueFaces)
			{
				if (face.isVisible(v))
				{
					viewFaces.push_back(face);
				}
			}

			if (viewFaces.size() == 0)
			{
				continue;
			}

			// Suppression des faces visibles des faces valides
			for each (face in viewFaces)
			{
				int index;
				for (int y = 0; y < trueFaces.size(); y++)
				{
					if (trueFaces[y].index0 == face.index0 && trueFaces[y].index1 == face.index1 && trueFaces[y].index2 == face.index2) {
						index = y;
						break;
					}
				}
				trueFaces.erase(trueFaces.begin() + index);
			}

			if (viewFaces.size() == 1)
			{
				face = viewFaces[0];
				trueFaces.push_back(Face3D(i, face.index0, face.index1));
				trueFaces.push_back(Face3D(i, face.index1, face.index2));
				trueFaces.push_back(Face3D(i, face.index2, face.index0));
				continue;
			}

			tmpFaces.clear();
			// Calcul de toutes les faces possibles
			for each(face in viewFaces)
			{
				tmpFaces.push_back(Face3D(i, face.index0, face.index1));
				tmpFaces.push_back(Face3D(i, face.index1, face.index2));
				tmpFaces.push_back(Face3D(i, face.index2, face.index0));
			}

			for each(face in tmpFaces)
			{
				// Recherche des faces qui ne croiseraient pas l'enveloppe convexe
				bool search = false;
				for each(Face3D other in tmpFaces)
				{
					if (face.index0 != other.index0 || face.index1 != other.index1 || face.index2 != other.index2)
					{
						if (face.isVisible(other.centroid()))
						{
							search = true;
							break;
						}
					}
				}
				if (!search) trueFaces.push_back(face);
			}
		}

		// Ajout des faces correctes
		for each(face in trueFaces)
		{
			result.push_back(face.index0);
			result.push_back(face.index1);
			result.push_back(face.index2);
		}
		return result;
	}
};