#include "stdafx.h"
#include "AlgoMath.h"

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