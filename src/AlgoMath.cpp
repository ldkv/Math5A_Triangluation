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

vector<QVector3D> EnvelopeJarvis(vector<Point> pts)
{
	vector<QVector3D> poly;
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
		poly.push_back(Pi.coord);
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

bool myfunction(Point i, Point j) { return (i.coord.x()==j.coord.x() ? i.coord.y() < j.coord.y() : i.coord.x() < j.coord.x()); }

QVector3D calculVector(Point p1, Point p2)
{
	int x = p2.coord.x() - p1.coord.x();
	int y = p2.coord.y() - p1.coord.y();
	int z = p2.coord.z() - p1.coord.z();
	return QVector3D(x, y, z);
}

vector<Side> TriangulationSimple(vector<Point> pts) {
	vector<Side> sortedPts;
	std::sort(pts.begin(), pts.end(), myfunction);
	int i=0;
	if (pts.size() >= 2)
	{
		for (i = 0; i < pts.size() - 2; i++)
		{
			QVector3D vec1 = calculVector(pts[i], pts[i + 1]);
			QVector3D vec2 = calculVector(pts[i], pts[i + 2]);

			sortedPts.push_back(Side(pts[i].id, pts[i+1].id));

			if (vec1.x() * vec2.y() != vec1.y() * vec2.x())
			{
				i+=2;
				break;
			}
			else
			{
				sortedPts.push_back(Side(pts[i+1].id, pts[i+2].id));
			}
		}
		if (i >= 1 && i < pts.size())
			for (int y = i-1; y >= 0; y--)
			{
				sortedPts.push_back(Side(pts[y].id, pts[i].id));
			}
	}

	return sortedPts;
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