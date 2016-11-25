#include "stdafx.h"
#include "AlgoMath.h"

#include <algorithm>  
#include <QDebug>

void resetGlobalID()
{
	globalId = 0;
	globalSideId = 0;
	globalFaceId = 0;
}

Point getPointfromID(vector<Point> pts, int id)
{
	for (int i = 0; i < pts.size(); i++)
	{
		if (pts[i].id == id)
			return pts[i];
	}
	return Point();	// id = -1 => condition pour vérifier si on a bien trouvé id désiré
}

Side getSidefromID(vector<Side> sds, int id)
{
	for (int i = 0; i < sds.size(); i++)
	{
		if (sds[i].id == id)
			return sds[i];
	}
	return Side(); // id = -1 => condition pour vérifier si on a bien trouvé id désiré
}

Face getFacefromID(vector<Face> faces, int id)
{
	for (int i = 0; i < faces.size(); i++)
	{
		if (faces[i].id == id)
			return faces[i];
	}
	return Face();	// id = -1 => condition pour vérifier si on a bien trouvé id désiré
}

void deleteSidefromID(int id, vector<Side> &sds, vector<Point> &pts)
{
	for (int i = 0; i < sds.size(); i++)
	{
		if (sds[i].id == id)
		{
			deleteSidefromPoint(pts, sds[i].id, sds[i].pLow);
			deleteSidefromPoint(pts, sds[i].id, sds[i].pHigh);
			sds.erase(sds.begin() + i);
			break;
		}
	}
}

void deleteSidefromPoint(vector<Point> &pts, int sideID, int pointID)
{
	for (int i = 0; i < pts.size(); i++)
	{
		if (pts[i].id == pointID)
		{
			for (int j = 0; j < pts[i].sides.size(); j++)
			{
				if (pts[i].sides[j] == sideID)
				{
					pts[i].sides.erase(pts[i].sides.begin() + j);
					break;
				}
			}
			break;
		}
	}
}

void deleteFacefromSide(vector<Side> &sides, int sideID, int faceID)
{
	for (int i = 0; i < sides.size(); i++)
	{
		if (sides[i].id == sideID)
		{
			if (sides[i].fLeft == faceID)
				sides[i].fLeft = -1;
			else
				sides[i].fRight = -1;
			break;
		}
	}
}

void deleteFacefromID(int faceID, vector<Face> &faces, vector<Side> &sides)
{
	for (int i = 0; i < faces.size(); i++)
	{
		if (faces[i].id == faceID)
		{
			for each (int sideID in faces[i].sides)
				deleteFacefromSide(sides, sideID, faceID);
			faces.erase(faces.begin() + i);
			break;
		}
	}
}

void Delaunay_addPoint(vector<Point> &pts, vector<Side> &sides, vector<Face> &faces, QVector3D P)
{
	// A) T ne contient pas de triangle
	if (faces.size() <= 0)
	{
		// A.1) T est vide
		if (pts.size() < 1)
		{
			pts.push_back(Point(P));
			return;
		}
		// A.2) T contient 1 seul point
		if (pts.size() == 1)
		{
			pts.push_back(Point(P));
			sides.push_back(Side(pts[0].id, pts[1].id));
			pts[0].sides.push_back(sides[0].id);
			pts[1].sides.push_back(sides[0].id);
			return;
		}
		// A.3) T contient plusieurs points et ils sont tous colinéaires
		else
		{
			// A.3.1) P est colinéaires aux sommets de T
			if (Collinear(pts[0].coord - P, pts[0].coord - pts[1].coord))
			{
				//Point bary = barycenter(pts);
			}
			// A.3.2) P n'est pas colinéaire
			else
			{
				Point Snew = Point(P);
				vector<Side> La;
				for (int i = 0; i < pts.size(); i++)
				{
					La.push_back(Side(pts[i].id, Snew.id));
					pts[i].sides.push_back(La[La.size() - 1].id);
					Snew.sides.push_back(La[La.size() - 1].id);
					//if (i < N - 1)
						//faces.push_back(Face(pts[i], pts[i + 1], Snew));
				}

				for (int i = 0; i < sides.size(); i++)
				{
					int iL,iH;
					for (int j = 0; j < La.size(); j++)
					{
						if (La[j].pLow == sides[i].pLow)
							iL = j;
						if (La[j].pLow == sides[i].pHigh)
							iH = j;
					}
					faces.push_back(Face(sides[i].id, La[iL].id, La[iH].id));
					sides[i].fLeft = faces[faces.size() - 1].id;
					La[iL].fLeft = faces[faces.size() - 1].id;
					La[iH].fLeft = faces[faces.size() - 1].id;
				}
				for each (Side s in La)
					sides.push_back(s);
				pts.push_back(Snew);
			}
		}
	}
	// B) T contient des triangles
	else
	{
		// B.1) Déterminer la liste des arêtes La
		Point Snew = Point(P);
		vector<Side> La;
		int idInside;
		for (idInside = 0; idInside < faces.size(); idInside++)
		{
			vector<Point> v = getVertexesfromFace(faces[idInside], pts, sides);
			if (v.size() == 3 && insideTriangle(P, v[0].coord, v[1].coord, v[2].coord))
				break;
		}
		// B.1.1) P inside a triangle
		if (idInside < faces.size())
		{
			vector <int> sideIds = faces[idInside].sides;
			deleteFacefromID(faces[idInside].id, faces, sides);
			for each(int id in sideIds)
				La.push_back(getSidefromID(sides, id));
		}
		// B.1.2) P hors de la triangulation
		else
		{
			La = getViewedEdges(sides, pts, Snew);
			qDebug() << "La size = " + La.size() << " && nbPoints = " + pts.size();
		}

		// B.2) Traitement des arêtes La
		while (La.size() > 0)
		{
			Side a = La[La.size() - 1];
			La.pop_back();
			Face f = getFacefromID(faces, a.fLeft);
			vector<Point> v = getVertexesfromFace(f, pts, sides);
			// B.2.1) P est intérieur du triangle d'incidient de l'arête a
			if (v.size() == 3 && insideCircumCircle(P, v[0].coord, v[1].coord, v[2].coord))
			{
				vector <int> sideIds = f.sides;
				deleteFacefromID(f.id, faces, sides);// IL FAUT AUSSI DEAFFECTER LES ARETES !!!
				deleteSidefromID(a.id, sides, pts);
				// Ajouter a1 et a2 dans La (sans a)
				for each (int id in sideIds)
					if (id != a.id)
						La.push_back(getSidefromID(sides, id));
			}
			else
			{
				addSide(a.pLow, Snew.id, sides, pts);
				addSide(a.pHigh, Snew.id, sides, pts);
				int s1 = getSideIDfromPoints(sides, a.pLow, Snew.id);
				int s2 = getSideIDfromPoints(sides, a.pHigh, Snew.id);
				int k = 0;
				if (s1 != -1 && s2 != -1)
				{
					faces.push_back(Face(s1, s2, a.id));
					int newFaceID = faces[faces.size() - 1].id;
					for (int i = 0; i < sides.size(); i++)
					{
						if (sides[i].id == a.id || sides[i].id == s1 || sides[i].id == s2)
						{
							if (sides[i].fLeft == -1)
								sides[i].fLeft = newFaceID;
							else
								sides[i].fRight = newFaceID;
						}
					}
				}
			}
		}
		pts.push_back(Snew);
	}
}

int getSideIDfromPoints(vector<Side> sides, int p1, int p2)
{
	for (int i = 0; i < sides.size(); i++)
		if ((sides[i].pLow == p1 && sides[i].pHigh == p2) || (sides[i].pLow == p2 && sides[i].pHigh == p1))
			return sides[i].id;
	return -1;
}

int addSide(int p1, int p2, vector<Side> &sides, vector<Point> &pts)
{
	if (getSideIDfromPoints(sides, p1, p2) >= 0)
		return -1;
	sides.push_back(Side(p1, p2));
	int id = sides[sides.size() - 1].id;
	for (int i = 0; i < pts.size(); i++)
		if (pts[i].id == p1 || pts[i].id == p2)
			pts[i].sides.push_back(id);
	return id;
}

vector<Point> getVertexesfromFace(Face F, vector<Point> pts, vector<Side> sides)
{
	vector<Point> pts_list;
	if (F.sides.size() == 3)
	{
		int id1, id2, id3;
		id1 = getSidefromID(sides, F.sides[0]).pHigh;
		id2 = getSidefromID(sides, F.sides[0]).pLow;
		id3 = getSidefromID(sides, F.sides[1]).pHigh;
		pts_list.push_back(getPointfromID(pts, id1));
		pts_list.push_back(getPointfromID(pts, id2));
		if (id3 != id1 && id3 != id2)
			pts_list.push_back(getPointfromID(pts, id3));
		else
			pts_list.push_back(getPointfromID(pts, getSidefromID(sides, F.sides[1]).pLow));
	}
	return pts_list;
}

float sign(QVector3D p1, QVector3D p2, QVector3D p3)
{
	return (p1.x() - p3.x()) * (p2.y() - p3.y()) - (p2.x() - p3.x()) * (p1.y() - p3.y());
}

bool insideTriangle(QVector3D pt, QVector3D v1, QVector3D v2, QVector3D v3)
{
	bool b1, b2, b3;
	b1 = sign(pt, v1, v2) < 0.0f;
	b2 = sign(pt, v2, v3) < 0.0f;
	b3 = sign(pt, v3, v1) < 0.0f;

	return ((b1 == b2) && (b2 == b3));
}

bool insideCircumCircle(QVector3D pt, QVector3D v1, QVector3D v2, QVector3D v3)
{
	float alpha = (v2 - v3).lengthSquared() * QVector3D::dotProduct(v1 - v2, v1 - v3) / (2 * QVector3D::crossProduct(v1 - v2, v2 - v3).lengthSquared());
	float beta = (v1 - v3).lengthSquared() * QVector3D::dotProduct(v2 - v1, v2 - v3) / (2 * QVector3D::crossProduct(v1 - v2, v2 - v3).lengthSquared());
	float theta = (v1 - v2).lengthSquared() * QVector3D::dotProduct(v3 - v1, v3 - v2) / (2 * QVector3D::crossProduct(v1 - v2, v2 - v3).lengthSquared());

	QVector3D center = alpha*v1 + beta*v2 + theta*v3;
	float radius = (v1 - v2).length()*(v2 - v3).length()*(v3 - v1).length() / (2 * QVector3D::crossProduct(v1 - v2, v2 - v3).length());

	return center.distanceToPoint(pt) <= radius;
}

bool Collinear(QVector3D v1, QVector3D v2)
{
	return (QVector3D::crossProduct(v1, v2) == QVector3D(0, 0, 0));
}

// BUG: infinite loop if 2 points confondus
vector<Point> EnvelopeJarvis(vector<Point> pts)
{
	vector<Point> poly;
	int N = pts.size();
	if (N <= 0)
		return poly;
	
	QVector3D minPoint = pts[0].coord;
	int i0 = 0;
	for (int i = 1; i < N; i++)
	{
		if (pts[i].coord.x() < minPoint.x() || (pts[i].coord.x() == minPoint.x() && pts[i].coord.y() < minPoint.y()))
		{
			minPoint = pts[i].coord;
			i0 = i;
		}
	}

	QLineF v(0, 0, 0, -1);
	int i = i0;
	do
	{
		Point Pi = pts[i];
		poly.push_back(Pi);
		if (N <= 1)
			break;
		// recherche du point suivant
		// initialisation de angle_min et lmax avec le premier point d'indice différent de i
		int j = (i == 0) ? 1 : 0;
		QLineF PiPj(Pi.coord.x(), Pi.coord.y(), pts[j].coord.x(), pts[j].coord.y());
		qreal angle_min = v.angleTo(PiPj);
		float lmax = Pi.coord.distanceToPoint(pts[j].coord);
		int inew = j;
		// recherche du point le plus proche
		for (j = inew + 1; j < N; j++)
		{
			if (j != i)
			{
				PiPj = QLineF(Pi.coord.x(), Pi.coord.y(), pts[j].coord.x(), pts[j].coord.y());
				qreal angle = v.angleTo(PiPj);
				float l = Pi.coord.distanceToPoint(pts[j].coord);
				if (angle_min > angle || (angle_min == angle && lmax < l))
				{
					angle_min = angle;
					lmax = l;
					inew = j;
				}
			}
		}
		// mise à jour du pivot et du vecteur directeur
		v = QLineF(pts[i].coord.x(), pts[i].coord.y(), pts[inew].coord.x(), pts[inew].coord.y());
		i = inew;
	} while (i != i0);
	return poly;	
}

bool coordsSort(Point i, Point j) { return (i.coord.x()==j.coord.x() ? i.coord.y() < j.coord.y() : i.coord.x() < j.coord.x()); }

/*QVector3D calculVector(Point p1, Point p2)
{
	int x = p2.coord.x() - p1.coord.x();
	int y = p2.coord.y() - p1.coord.y();
	int z = p2.coord.z() - p1.coord.z();
	return QVector3D(x, y, z);
}*/

vector<Face> TriangulationSimple(vector<Point> pts) {
	// 1
	vector<Side> sides;
	list<Side> convexHull;
	vector<Face> faces;
	vector<Point> colinearPoints;
	std::sort(pts.begin(), pts.end(), coordsSort);
	// 2
	int i = 1;
	if (pts.size() > 2)
	{
		// 2a
		colinearPoints.push_back(pts[0]);
		for (i; i < pts.size() - 1; i++)
		{
			//QVector3D vec1 = calculVector(pts[i], pts[i + 1]);
			QVector3D vec1 = pts[0].coord - pts[i].coord;
			QVector3D vec2 = pts[0].coord - pts[i + 1].coord;

			if (vec1.x() * vec2.y() == vec1.y() * vec2.x())
			{
				colinearPoints.push_back(pts[i]);
			}
			else
			{
				break;
			}
		}


		// 2.b) Construction of the root Triangulation
		if (i < pts.size())
		{
			if (colinearPoints.size() > 1)
			{
				unsigned int j;
				for (j = 1; j < colinearPoints.size(); j++)
				{
					convexHull.push_back(Side(pts[j], pts[j - 1]));
					convexHull.push_back(Side(pts[j - 1], pts[i]));
					faces.push_back(Face(pts[j], pts[j - 1], pts[i]));
				}
				convexHull.push_back(Side(pts[i], pts[j - 1]));
			}
			else {
				faces.push_back(Face(pts[0], pts[1], pts[2]));
				convexHull.push_back(Side(pts[0], pts[1]));
				convexHull.push_back(Side(pts[1], pts[2]));
				convexHull.push_back(Side(pts[2], pts[0]));
				i = 2;
			}
		}

		// 3) Search for all edges viewed by Pi to build the next Triangulation
		while (++i < pts.size())
		{
			Point nextVert = pts[i];

			vector<Side> viewedEdges;
			viewedEdges = getViewedEdge(i, pts, convexHull);

			// Add a new face with the each viewed edge
			vector<Side>::iterator edgeView;
			for (edgeView = viewedEdges.begin(); edgeView != viewedEdges.end(); edgeView++)
			{
				//int idVertA = edgeView->points[0];
				//int idVertB = edgeView->points[1];
				Point vertA = edgeView->points[0];
				Point vertB = edgeView->points[1];

				faces.push_back(Face(vertA, nextVert, vertB));
			}
		}
	}

	return faces;
}

// Return an array of all the edge "viewed" by the given Vertex
vector<Side> getViewedEdges(vector<Side> sides, vector<Point> pts, Point P)
{
	vector<Side> viewedEdges;
	// For each Edge in our list of edges
	for (int i = 0; i < sides.size(); i++)
	{
		if (sides[i].fLeft >= 0 && sides[i].fRight >= 0)
			continue;
		// Get components the current edge
		Point A = getPointfromID(pts, sides[i].pLow);
		Point B = getPointfromID(pts, sides[i].pHigh);

		Point C = pts[0];
		int o = 0;
		while ((C.coord == A.coord || C.coord == B.coord) && o + 1 < pts.size()) {
			C = pts[++o];
		}

		// Compute face normal formed by those 3 Vertex
		QVector3D u = B.coord - A.coord;
		QVector3D v = C.coord - A.coord;
		QVector3D normal = crossProductNormalized(u, v);

		// check if AB edge is viewed by the vertex C
		if (isEdgeViewed(P.coord, A.coord, B.coord, normal))
			viewedEdges.push_back(sides[i]);
	}
	return viewedEdges;
}

// Return an array of all the edge "viewed" by the given Vertex
vector<Side> getViewedEdge(int nextIdVert, vector<Point> pts, list<Side> &convexHull)
{
	vector<Side> viewedEdge;

	list<Side>::iterator itE;
	list<Side>::iterator fromEdge = convexHull.end();
	list<Side>::iterator toEdge = convexHull.begin();

	Point nextVert = pts[nextIdVert];

	// For each Edge in our convex hull
	for (itE = convexHull.begin(); itE != convexHull.end(); ++itE)
	{
		// Get components the current edge
		Point A = itE->points[0];
		Point B = itE->points[1];

		// Get any Vertex in our convex Hull different from A or B
		// Cstart = convexHull.begin();
		// std::advance(Cstart, ++o);
		Point C = pts[0];
		int o = 0;
		while ((C.coord == A.coord || C.coord == B.coord) && o+1 < pts.size()) {
			C = pts[++o];
		}

		// Compute face normal formed by those 3 Vertex
		QVector3D u = B.coord - A.coord;
		QVector3D v = C.coord - A.coord;
		QVector3D normal = crossProductNormalized(u, v);
		//QVector3D normal ( 0.f, 0.f, 1.f);

		// Get the real Position for the 3 vertex
		QVector3D a = A.coord;
		QVector3D b = B.coord;
		QVector3D c = C.coord;

		// check if AB edge is viewed by the vertex C
		if (isEdgeViewed(nextVert.coord, a, b, normal))
		{
			viewedEdge.push_back(Side(A, B));

			// save position of the first viewed edge.
			if (fromEdge == convexHull.end())
				fromEdge = itE;
			// and store the last one as well
			toEdge = itE;
		}
	}


	vector<Side>::iterator itviewEdg;
	Point firstVert = viewedEdge.front().points[0];
	Point lastVert = viewedEdge.back().points[1];

	convexHull.insert(toEdge, Side(firstVert, nextVert));
	convexHull.insert(toEdge, Side(nextVert, lastVert));
	
	for (itviewEdg = viewedEdge.begin(); itviewEdg != viewedEdge.end(); itviewEdg++)
		convexHull.remove(*itviewEdg);
	
	return viewedEdge;
}


float dot(QVector3D p, QVector3D op) { return p.x()*op.x() + p.y()*op.y() + p.z()*op.z(); } // Produit scalaire de 2 points
float norme(QVector3D p) { return float(sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z())); } // Norme du point (longueur du vecteur.)
//T distance(const Point<T>& p) { return sqrt(pow(p.x - x, 2) + pow(p.y - y, 2) + pow(p.z - z, 2)); } // Distance avec le vecteur en paramètre
//T angle(Point<T> p) { return acos(static_cast<T>(dot(p) / (norme() * p.norme()))); } // Return l'angle non orienté formé par les 2 vecteurs.
QVector3D crossProduct(QVector3D p, QVector3D op) { QVector3D t(p.y()*op.z() - p.z()*op.y(), p.z()*p.x() - p.x()*op.z(), p.x()*op.y() - p.y()*op.x()); return t; } //Produit vectoriel de 2 points.
QVector3D normalize(QVector3D p) { float n = norme(p); return QVector3D(p.x() / n, p.y() / n, p.z() / n);} //Normalisation du point.
QVector3D crossProductNormalized(QVector3D p, QVector3D op) 
{ 
	QVector3D final;
	final = crossProduct(p, op);
	final = normalize(final);
	return final;
} // Produit vectoriel normalisé de 2 points.
bool isEdgeViewed(QVector3D P, QVector3D A, QVector3D B, QVector3D n)
{
	QVector3D u = B - A; // Build a vector interior to AB
	QVector3D v = P - A;

	// Compute the normal interior to the Edge AB
	QVector3D normalAB = crossProductNormalized(u, n);

	return dot(normalAB, v) > 0 ? true : false;
}
/*
void circumCenter2D(Point<T> A, Point<T> B, Point<T> C) // Compute the circumCenter of the triangle ABC
{
	// explanation : http://en.wikipedia.org/wiki/Circumscribed_circle
	T ax = A.x, ay = A.y, bx = B.x, by = B.y, cx = C.x, cy = C.y;
	T ax2 = A.x*A.x, ay2 = A.y*A.y, bx2 = B.x*B.x, by2 = B.y*B.y, cx2 = C.x*C.x, cy2 = C.y*C.y;

	T d = 2 * (ay*cx + by*ax - by*cx - ay*bx - cy*ax + cy*bx);

	x = by*ax2 - cy*ax2 - by2*ay + cy2*ay + bx2*cy + ay2*by + cx2*ay - cy2*by - cx2*by - bx2*ay + by2*cy - ay2*cy;
	y = ax2*cx + ay2*cx + bx2*ax - bx2*cx + by2*ax - by2*cx - ax2*bx - ay2*bx - cx2*ax + cx2*bx - cy2*ax + cy2*bx;

	x /= d;
	y /= d;
	z = 0;
}
*/


float* CrossProduct(float *a, float *b)
{
	float Product[3];

	//Cross product formula 
	Product[0] = (a[1] * b[2]) - (a[2] * b[1]);
	Product[1] = (a[2] * b[0]) - (a[0] * b[2]);
	Product[2] = (a[0] * b[1]) - (a[1] * b[0]);

	return Product;
}

/*
QVector3D normalize(QVector3D v) {

	double length = sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());

	v.x() = v.x() / length;
	v.y() = v.y() / length;
	v.z() = v.z() / length;
}
*/
/*double dotProduct(glm::vec2 vecA, glm::vec2 vecB)
{
	//return normVector(vecA) * normVector(vecB) * glm::cos (vecA, vecB);
	//return (p1.x * p2.x) + (p1.y * p2.y);
	return glm::dot(vecA, vecB);
}*/
/*
bool checkVisibilitySide(Side side, Point px, vector<Point> pts, Side side2) {
	QVector3D p1 = QVector3D(getPointfromID(pts, side.pLow)->coord);
	QVector3D p2 = QVector3D(getPointfromID(pts, side.pHigh)->coord);
	QVector3D p3 = QVector3D(getPointfromID(pts, side.pHigh)->coord);
	float dx = p2.x() - p1.x();
	float dy = p2.y() - p2.y();
	//QVector3D dz = z2 - z1;
	QVector3D normal1 = normalize(QVector3D (-dy, dx, 0));
	QVector3D normal2 = normalize(QVector3D(dy, -dx, 0));

	QVector3D finalNormal;
	if (QVector3D::dotProduct(p3 - p1, normal1) > 0) {
		finalNormal = normal1;
	}
	else {
		finalNormal = normal2;
	}

	float dp = QVector3D::dotProduct(finalNormal
		, QVector3D(getPointfromID(pts, side.pLow)->coord - px.coord));
	return dp < 180 ? true : false;
}
*/
/*bool checkVisibilityEdge(Side &edge, Point &point)
{
	QVector3D a; a.dotProduct();
	int value = dotProduct(edge, makeVector(edge., point));
	if (value < 0)
		return true;
	return false;
}*/



int getSideIDFromPoints(vector<Side> s, Point x, Point y) {
	for (int i = 0; i < s.size(); i++)
	{
		if ((s[i].pLow == x.id && s[i].pHigh == y.id) || (s[i].pHigh == x.id && s[i].pLow == y.id)) {
			return s[i].id;
		}
	}
}

int getPointIndex(vector<Point> pts, int id) {
	for (int i = 0; i < pts.size(); i++)
	{
		if (pts[i].id == id) {
			return i;
		}
	}
	return -1;
}

Point barycenter(vector<Point> points)
{
	Point bary(QVector3D(0,0,0));
	int s = points.size();
	int sX=0, sY=0, sZ=0;

	for (int i = 0; i < s; i++)
	{
		sX += points[i].coord.x();
		sY += points[i].coord.y();
		sZ += points[i].coord.z();
	}

	bary = Point(QVector3D(sX/s, sY/s, sZ/s));

	return bary;
}
struct angleSort {
	angleSort(Point pointsBarycenter) { this->pointsBarycenter = pointsBarycenter; }
	bool operator () (Point i, Point j) {
		QLineF lx(pointsBarycenter.coord.x(), pointsBarycenter.coord.y(), 1, 0);
		QLineF li(pointsBarycenter.coord.x(), pointsBarycenter.coord.y(), i.coord.x(), i.coord.y());
		QLineF lj(pointsBarycenter.coord.x(), pointsBarycenter.coord.y(), j.coord.x(), j.coord.y());
		float anglei = lx.angleTo(li);
		float anglej = lx.angleTo(lj);
		return anglei ==  anglej ? i.coord.distanceToPoint(pointsBarycenter.coord) < j.coord.distanceToPoint(pointsBarycenter.coord) : anglei < anglej;
	}
	Point pointsBarycenter;
};

/*bool angleSort(Point i, Point j) {
	/*float angle = atan2(vecX.x(), i.coord.x()) - atan2(vecX.y(), i.coord.y());
	if (angle >= M_PI) angle -= 2 * M_PI;
	if (angle <= -M_PI) angle += 2 * M_PI;
	float angle2 = atan2(vecX.x(), j.coord.x()) - atan2(vecX.y(), j.coord.y());
	if (angle2 >= M_PI) angle2 -= 2 * M_PI;
	if (angle2 <= -M_PI) angle2 += 2 * M_PI;
	QLineF lx(0, 0, 1, 0);
	QLineF li(0, 0, i.coord.x(), i.coord.y());
	QLineF lj(0, 0, j.coord.x(), j.coord.y());
	float anglei = lx.angleTo(li);
	float anglej = lx.angleTo(lj);
	return anglei < anglej; 
}*/

int magnitudeVector(QVector3D v1, QVector3D v2)
{
	return sqrt(pow(v2.x() - v1.x(), 2)+ pow(v2.y() - v1.y(), 2));
}
/*
double angleSigned(QVector3D v1, QVector3D v2)
{
	std::inner_product()
	dotProduct(v1, v2);
	//det(u,v)=ux*vy-uy*vx
	dot / (magnitudeVector(v1) * magnitudeVector(v2));

	float angle = acos(result);

	return angle * 180 / M_PI;
}*/

int nextIndexPoint(int currentIndex, int vectorSize) {
	return (currentIndex+1 >= vectorSize ?  currentIndex+1 - vectorSize : currentIndex+1);
}
int previousIndexPoint(int currentIndex, int vectorSize) {
	return (currentIndex-1 < 0 ? currentIndex-1 + vectorSize : currentIndex-1);
}

vector<Point> GrahamScan(vector<Point> pts) {
	if (pts.size() > 2) {
		Point bary = barycenter(pts);
		std::sort(pts.begin(), pts.end(), angleSort(bary));

		vector<Point> L = pts;
		//Point Sinit = L.at(0);
		//Point pivot = Sinit;
		int Sinit = 0;
		int pivot = Sinit;
		bool pass = false;
		do
		{
			QLineF PiPj(L[pivot].coord.x(), L[pivot].coord.y(), L[nextIndexPoint(pivot, L.size())].coord.x(), L[nextIndexPoint(pivot, L.size())].coord.y());
			QLineF PiPk(L[pivot].coord.x(), L[pivot].coord.y(), L[previousIndexPoint(pivot, L.size())].coord.x(), L[previousIndexPoint(pivot, L.size())].coord.y());
			if (PiPk.angleTo(PiPj) > 180 )
			{
				pivot = nextIndexPoint(pivot, L.size());
				pass = true;
			}
			else
			{
				L.erase(L.begin() + pivot);
				Sinit = previousIndexPoint(pivot, L.size());
				pivot = Sinit;
				pass = false;
			}
		} while (pivot != Sinit || pass == false);
		return L;
	}
	return pts;
}
