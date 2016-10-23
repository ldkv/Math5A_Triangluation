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
bool coordsSort(Point i, Point j) { return (i.coord.x()==j.coord.x() ? i.coord.y() < j.coord.y() : i.coord.x() < j.coord.x()); }

/*QVector3D calculVector(Point p1, Point p2)
{
	int x = p2.coord.x() - p1.coord.x();
	int y = p2.coord.y() - p1.coord.y();
	int z = p2.coord.z() - p1.coord.z();
	return QVector3D(x, y, z);
}*/

vector<Side> TriangulationSimple(vector<Point> pts) {
	vector<Side> sortedPts;
	std::sort(pts.begin(), pts.end(), coordsSort);
	int i=0;
	if (pts.size() >= 2)
	{
		for (i = 0; i < pts.size() - 2; i++)
		{
			//QVector3D vec1 = calculVector(pts[i], pts[i + 1]);
			QVector3D vec1 = pts[i].coord - pts[i + 1].coord;
			QVector3D vec2 = pts[i].coord - pts[i + 2].coord;

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
		return anglei < anglej;
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

vector<Point> GrahamScan(vector<Point> pts) {
	if (pts.size() > 0) {
		Point bary = barycenter(pts);
		std::sort(pts.begin(), pts.end(), angleSort(bary));
	}
	return pts;
}
