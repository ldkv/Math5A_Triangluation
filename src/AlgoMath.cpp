#include "stdafx.h"
#include "AlgoMath.h"

#include <algorithm>  

Point* getPointfromID(vector<Point> pts, int id)
{
	for (int i = 1; i < pts.size(); i++)
	{
		if (pts[i].id == id)
			return &pts[i];
	}
	return nullptr;
}

vector<Point> EnvelopeJarvis(vector<Point> pts)
{
	vector<Point> poly;

	if (pts.size() <= 0)
		return poly;
	
	float minX = pts[0].coord.x(), minY = pts[0].coord.y();
	int minID = pts[0].id;
	for (int i = 1; i < pts.size(); i++)
	{
		if (pts[i].coord.x() < minX || (pts[i].coord.x() == minX && pts[i].coord.y() < minY))
		{
			minX = pts[i].coord.x();
			minY = pts[i].coord.y();
			minID = pts[i].id;
		}
	}

	QLineF d(0, 0, 0, -1);
	int i = minID;
	do
	{
		Point Pi = *getPointfromID(pts, i);
		poly.push_back(Pi);
		int j = (i == 0) ? 1 : 0;
		QLineF PiPj(Pi.coord.x(), Pi.coord.y(), pts[j].coord.x(), pts[j].coord.y());
		qreal angle_min = d.angleTo(PiPj);


	} while (true);
	
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