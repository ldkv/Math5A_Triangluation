#include "stdafx.h"
#include "AlgoMath.h"

#include <algorithm>  

Point *getPointfromID(vector<Point> pts, int id)
{
	for (int i = 0; i < pts.size(); i++)
	{
		if (pts[i].id == id)
			return &pts[i];
	}
	return nullptr;
}

void Delaunay_addPoint(vector<Point> &pts, vector<Side> &sides, vector<Face> &faces, QVector3D P)
{
	// Si T ne contient pas de triangle
	if (faces.size() <= 0)
	{
		// Si T est vide
		if (pts.size() < 1)
		{
			pts.push_back(Point(P));
			return;
		}
		// Si T contient 1 seul point
		if (pts.size() == 1)
		{
			pts.push_back(Point(P));
			sides.push_back(Side(pts[0].id, pts[1].id));
			return;
		}
		// Si T contient plusieurs points et ils sont tous colinéaires
		else
		{
			// Si P est colinéaires aux sommets de T
			if (Collinear(pts[0].coord - P, pts[0].coord - pts[1].coord))
			{
				//Point bary = barycenter(pts);
			}
			// P n'est pas colinéaire
			else
			{
				pts.push_back(Point(P));
				int N = pts.size();
				int S = sides.size();	// Nombre original des côtés
				for (int i = 0; i < N - 1; i++)
				{
					sides.push_back(Side(pts[i].id, pts[N - 1].id));
					pts[i].sides.push_back(sides[sides.size() - 1].id);
					pts[N - 1].sides.push_back(sides[sides.size() - 1].id);
				}

				for (int i = 0; i < S; i++)
				{
					
				}
			}
		}
	}
}

bool Collinear(QVector3D v1, QVector3D v2)
{
	//return (v1.y()*v2.z() == v1.z()*v2.y() && v1.z()*v2.x() == v1.x()*v2.z() && v1.x()*v2.y() == v1.y()*v2.x());
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


vector<Side> getViewedEdge(int nextIdVert, std::vector<Point> pts, std::list<Side> &convexHull)
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

	std::vector<Side>::iterator itviewEdg;
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

vector<Face> Fliping(vector<Face> faces)
{
	std::vector<Side> & sides = faces[0].sides
	if (sides.size() > 2)
	{
		std::vector<Side>::iterator itE;
		for (itE = edges.begin(); itE != edges.end(); ++itE)
		{
			int directFace = itE->getFaceA();
			int adjacentFace = itE->getFaceB();

			if (directFace != -1 && adjacentFace != -1) // Edge is between 2 faces
			{
				Face & faceA = m_mesh.getFace(directFace);
				Face & faceB = m_mesh.getFace(adjacentFace);

				// Get the oposite vertex to the edge
				int idVertDirect = faceA.getThirdVertex(*itE);
				int idVertAdjac = faceB.getThirdVertex(*itE);
				assert(idVertDirect != -1 && idVertAdjac != -1);

				Vertex vertEdgeA = m_mesh.getVertice(itE->getVertexA());
				Vertex vertEdgeB = m_mesh.getVertice(itE->getVertexB());
				Vertex vertDirect = m_mesh.getVertice(idVertDirect);
				Vertex vertAdja = m_mesh.getVertice(idVertAdjac);

				// check if vertex respect Delaunay critera 
				float angleDirect = (vertEdgeA - vertDirect).angle(vertEdgeB - vertDirect);
				float angleAdja = (vertEdgeA - vertAdja).angle(vertEdgeB - vertAdja);

				// if nott respect critera
				if (angleDirect + angleAdja > M_PI)
				{
					// Update the faces
					faceA.setVertex(idVertDirect, itE->getVertexA(), idVertAdjac);
					faceB.setVertex(idVertDirect, idVertAdjac, itE->getVertexB());

					// update Edge
					// Edge[vertA, Adj] => Face A (direct)
					m_mesh.setEdge(itE->getVertexA(), idVertAdjac, directFace);
					// Edge[Dir, vertA] => Face A
					m_mesh.setEdge(idVertDirect, itE->getVertexA(), directFace);

					// Edge[vertB, Dir] => Face B (Adjacent)
					m_mesh.setEdge(itE->getVertexB(), idVertDirect, adjacentFace);
					// Edge[Adj, vertB] => Face B (Adjacent)
					m_mesh.setEdge(idVertAdjac, itE->getVertexB(), adjacentFace);

					//	flip the Edge
					itE->setVertexA(idVertAdjac);
					itE->setVertexB(idVertDirect);
				}
			}
		}
	}
}

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

Point barycenter(std::vector<Point> points)
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

int magnitudeVector(QVector3D v1, QVector3D v2)
{
	return sqrt(pow(v2.x() - v1.x(), 2)+ pow(v2.y() - v1.y(), 2));
}


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


/*
double angleSigned(QVector3D v1, QVector3D v2)
{
	std::inner_product()
		dotProduct(v1, v2);
	//det(u,v)=ux*vy-uy*vx
	dot / (magnitudeVector(v1) * magnitudeVector(v2));

	float angle = acos(result);

	return angle * 180 / M_PI;
}


bool checkVisibilitySide(Side side, Point px, vector<Point> pts, Side side2) {
	QVector3D p1 = QVector3D(getPointfromID(pts, side.pLow)->coord);
	QVector3D p2 = QVector3D(getPointfromID(pts, side.pHigh)->coord);
	QVector3D p3 = QVector3D(getPointfromID(pts, side.pHigh)->coord);
	float dx = p2.x() - p1.x();
	float dy = p2.y() - p2.y();
	//QVector3D dz = z2 - z1;
	QVector3D normal1 = normalize(QVector3D(-dy, dx, 0));
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

bool checkVisibilityEdge(Side &edge, Point &point)
{
	QVector3D a; a.dotProduct();
	int value = dotProduct(edge, makeVector(edge., point));
	if (value < 0)
		return true;
	return false;
]
*/
