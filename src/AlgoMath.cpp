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
	return Point();	// id = -1 => condition pour v�rifier si on a bien trouv� id d�sir�
}

Side getSidefromID(vector<Side> sds, int id)
{
	for (int i = 0; i < sds.size(); i++)
	{
		if (sds[i].id == id)
			return sds[i];
	}
	return Side(); // id = -1 => condition pour v�rifier si on a bien trouv� id d�sir�
}

Face getFacefromID(vector<Face> faces, int id)
{
	for (int i = 0; i < faces.size(); i++)
	{
		if (faces[i].id == id)
			return faces[i];
	}
	return Face();	// id = -1 => condition pour v�rifier si on a bien trouv� id d�sir�
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
			for each (int sideID in faces[i].sidesID)
				deleteFacefromSide(sides, sideID, faceID);
			faces.erase(faces.begin() + i);
			break;
		}
	}
}

vector<Point> diagramVoronoi(vector<Point> pts, vector<Side> sides, vector<Face> faces)
{
	vector<Point> results;
	float scale = 1000;
	if (pts.size() <= 1)
		return results;
	if (faces.size() == 0)
	{
		std::sort(pts.begin(), pts.end(), coordsSort);
		for (int i = 0; i < pts.size() - 1; i++)
		{
			QVector3D AB = pts[i + 1].coord - pts[i].coord;
			QVector3D M = (pts[i + 1].coord + pts[i].coord) / 2;
			float y1 = scale, y2 = -y1;
			float x1 = AB.y() / AB.x()*(M.y() - y1) + M.x();
			float x2 = AB.y() / AB.x()*(M.y() - y2) + M.x();
			results.push_back(Point(QVector3D(x1, y1, 0)));
			results.push_back(Point(QVector3D(x2, y2, 0)));
		}
		return results;
	}

	map<int, QVector3D> Ct;
	for (int i = 0; i < faces.size(); i++)
	{
		vector<Point> v = getVertexesfromFace(faces[i], pts, sides);
		Ct[i] = circumCircleCenter(v[0].coord, v[1].coord, v[2].coord);
	}
	
	vector<Point> A;
	for (int i = 0; i < sides.size(); i++)
	{
		QVector3D A1, A2;
		if (sides[i].fLeft != -1)
			A1 = Ct[sides[i].fLeft];
		if (sides[i].fRight != -1)
			A2 = Ct[sides[i].fRight];
		else
		{
			QVector3D S1 = getPointfromID(pts, sides[i].pLow).coord;
			QVector3D S2 = getPointfromID(pts, sides[i].pHigh).coord;
			QVector3D n = (S1 + S2) / 2 - A1;
			float x = A1.x() + n.x() * scale;
			float y = A1.y() + n.y() * scale;
			A2 = QVector3D(x, y, 0);
		}
		results.push_back(Point(A1));
		results.push_back(Point(A2));
	}

	for (int i = 0; i < pts.size(); i++)
	{
		 
	}

	return results;
}

void Delaunay_addPoint(vector<Point> &pts, vector<Side> &sides, vector<Face> &faces, QVector3D P)
{
	Point Snew = Point(P);
	// A) T ne contient pas de triangle
	if (faces.size() <= 0)
	{
		// A.1) T est vide
		if (pts.size() < 1)
		{
			pts.push_back(Snew);
			return;
		}
		// A.2) T contient 1 seul point
		if (pts.size() == 1)
		{
			pts.push_back(Snew);
			addSide(pts[0].id, Snew.id, sides, pts);
			return;
		}
		// A.3) T contient plusieurs points et ils sont tous colin�aires
		else
		{
			// A.3.1) P est colin�aires aux sommets de T
			if (Collinear(pts[0].coord - P, pts[0].coord - pts[1].coord))
			{
				std::sort(pts.begin(), pts.end(), coordsSort);
				// A.3.1.1)
				if (P.x() < pts[0].coord.x())
				{
					pts.push_back(Snew);
					addSide(Snew.id, pts[0].id, sides, pts);
					return;
				}
				// A.3.1.2)
				if (P.x() > pts[pts.size() - 1].coord.x())
				{
					pts.push_back(Snew);
					addSide(pts[pts.size() - 1].id, Snew.id, sides, pts);
					return;
				}
				// A.3.1.3)

			}
			// A.3.2) P n'est pas colin�aire
			else
			{
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
					int iL, iH;
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
		// B.1) D�terminer la liste des ar�tes La
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
			vector <int> sideIds = faces[idInside].sidesID;
			deleteFacefromID(faces[idInside].id, faces, sides);
			for each(int id in sideIds)
				La.push_back(getSidefromID(sides, id));
		}
		// B.1.2) P hors de la triangulation
		else
			La = getViewedEdges(sides, pts, Snew);

		pts.push_back(Snew);
		// B.2) Traitement des ar�tes La
		while (La.size() > 0)
		{
			Side a = La[La.size() - 1];
			La.pop_back();
			Face f = getFacefromID(faces, a.fLeft);
			vector<Point> v = getVertexesfromFace(f, pts, sides);
			// B.2.1) P est int�rieur du triangle d'incidient de l'ar�te a
			if (v.size() == 3 && insideCircumCircle(P, v[0].coord, v[1].coord, v[2].coord))
			{
				vector <int> sideIds = f.sidesID;
				deleteFacefromID(f.id, faces, sides);// IL FAUT AUSSI DEAFFECTER LES ARETES !!!
				deleteSidefromID(a.id, sides, pts);
				// Ajouter a1 et a2 dans La (sans a)
				for each (int id in sideIds)
					if (id != a.id)
						La.push_back(getSidefromID(sides, id));
			}
			// B.2.2) P est � l'ext�rieur du triangle
			else
			{
				addSide(a.pLow, Snew.id, sides, pts);
				addSide(a.pHigh, Snew.id, sides, pts);
				int s1 = getSideIDfromPoints(sides, a.pLow, Snew.id);
				int s2 = getSideIDfromPoints(sides, a.pHigh, Snew.id);
				int k = 0;
				if (s1 != -1 && s2 != -1)
					addFace(s1, s2, a.id, faces, sides);
			}
		}
	}
}

int addFace(int id1, int id2, int id3, vector<Face> &faces, vector<Side> &sides)
{
	faces.push_back(Face(id1, id2, id3));
	int newFaceID = faces[faces.size() - 1].id;
	for (int i = 0; i < sides.size(); i++)
	{
		if (sides[i].id == id1 || sides[i].id == id2 || sides[i].id == id3)
		{
			if (sides[i].fLeft == -1)
				sides[i].fLeft = newFaceID;
			else
				sides[i].fRight = newFaceID;
		}
	}
	return newFaceID;
}

void addSidetoPoint(vector<Point> &pts, int ptID, int sideID)
{
	for (int i = 0; i < pts.size(); i++)
		if (pts[i].id == ptID)
		{
			pts[i].sides.push_back(sideID);
			break;
		}
}

void removeSidefromPoint(vector<Point> &pts, int ptID, int sideID)
{
	for (int i = 0; i < pts.size(); i++)
		if (pts[i].id == ptID)
		{
			int k = 0;
			while (k++ < pts[i].sides.size() && pts[i].sides[k] != sideID);
			if (k < pts[i].sides.size())
				pts[i].sides.erase(pts[i].sides.begin() + k);
			break;
		}
}

void Delaunay_deletePoint(vector<Point> &pts, vector<Side> &sides, vector<Face> &faces, int deleteID)
{
	// A) T ne contient pas de triangle
	if (faces.size() <= 0)
	{
		// A.1) T contient 1 seul point
		if (pts.size() == 1)
		{
			pts.clear();
			return;
		}
		// A.2) T contient plusieurs points et ils sont tous colin�aires
		// A.2.1) P appartient � une seule ar�te
		if (pts[deleteID].sides.size() < 2)
		{
			if (pts[deleteID].sides.size() == 1)
				deleteSidefromID(pts[deleteID].sides[0], sides, pts);
			pts.erase(pts.begin() + deleteID);
			return;
		}
		// A.2.2) P appartient � 2 ar�tes
		Side a1 = getSidefromID(sides, pts[deleteID].sides[0]);
		int s1 = (a1.pHigh == pts[deleteID].id) ? a1.pLow : a1.pHigh;
		for (int i = 0; i < sides.size(); i++)
		{
			if (sides[i].id == pts[deleteID].sides[1])
			{
				if (sides[i].pHigh == pts[deleteID].id)
					sides[i].pHigh = s1;
				else
					sides[i].pLow = s1;
				addSidetoPoint(pts, s1, sides[i].id);
				break;
			}
		}
		deleteSidefromID(pts[deleteID].sides[0], sides, pts);
		return;
	}
	// B) T contient des triangles
	// B.1) D�terminer les ar�tes et triangles contenant P et les supprimer
	vector<Side> La2, La1 = getIncidentEdgesOriented(pts[deleteID], sides, pts);
	vector<Face> Lt = getIncidentFacesOriented(La1, faces, sides, La2);
	for each (Face f in Lt)
		deleteFacefromID(f.id, faces, sides);
	for each (Side s in La1)
		deleteSidefromID(s.id, sides, pts);
	
	QVector3D deletedP = pts[deleteID].coord;
	pts.erase(pts.begin() + deleteID);
	int N = La2.size();
	Point S1, S2;
	int foundVertex;
	// B.2.1) La2 forme un polygone ferm�
	if (N > 2 && mutualPoint(La2[0], La2[N-1]) >= 0)
	{
		while (N > 3)
		{
			foundVertex = findPointConvexDelaunay(true, La2, pts, deletedP, S1, S2);
			if (foundVertex != -1)
			{
				Side a1 = La2[foundVertex];
				Side a2 = La2[foundVertex + 1 < N ? foundVertex + 1 : 0];
				Side a3 = getSidefromID(sides, addSide(S1.id, S2.id, sides, pts));
				addFace(a1.id, a2.id, a3.id, faces, sides);

				for (int i = 0; i < La2.size(); i++)
					if (La2[i].id == a1.id)
					{
						La2.erase(La2.begin() + i);
						break;
					}
				for (int i = 0; i < La2.size(); i++)
					if (La2[i].id == a2.id)
					{
						La2.erase(La2.begin() + i);
						break;
					}
				La2.insert(La2.begin() + foundVertex, a3);
				N = La2.size();
			}
		}
		addFace(La2[0].id, La2[1].id, La2[2].id, faces, sides);
	}
	// B.2.2) La2 ne forme pas un polygone ferm�
	else
	{
		foundVertex = findPointConvexDelaunay(false, La2, pts, deletedP, S1, S2);
		while (foundVertex != -1)
		{
			Side a1 = La2[foundVertex];
			Side a2 = La2[foundVertex + 1];
			Side a3 = getSidefromID(sides, addSide(S1.id, S2.id, sides, pts));
			addFace(a1.id, a2.id, a3.id, faces, sides);

			for (int i = 0; i < La2.size(); i++)
				if (La2[i].id == a1.id)
				{
					La2.erase(La2.begin() + i);
					break;
				}
			for (int i = 0; i < La2.size(); i++)
				if (La2[i].id == a2.id)
				{
					La2.erase(La2.begin() + i);
					break;
				}
			La2.insert(La2.begin() + foundVertex, a3);
			N = La2.size();
			foundVertex = findPointConvexDelaunay(false, La2, pts, deletedP, S1, S2);
		}
	}
}

// Chercher un sommet convexe tel que le circonscrit du triangle form� de ses 2 ar�tes incidentes
// ne contient aucun sommet de La2
int findPointConvexDelaunay(bool closedPoly, vector<Side> La2, vector<Point> pts, QVector3D P, Point &S1, Point &S2)
{
	int N = La2.size();
	if (!closedPoly)
		N--;
	int i;
	for (i = 0; i < N; i++)
	{
		Side a1 = La2[i];
		Side a2 = closedPoly ? La2[i + 1 < N ? i + 1 : 0] : La2[i + 1];
		Point S = getPointfromID(pts, mutualPoint(a1, a2));
		S1 = getPointfromID(pts, a1.pHigh == S.id ? a1.pLow : a1.pHigh);
		S2 = getPointfromID(pts, a2.pHigh == S.id ? a2.pLow : a2.pHigh);
		QLineF SS1(S.coord.x(), S.coord.y(), S1.coord.x(), S1.coord.y());
		QLineF SS2(S.coord.x(), S.coord.y(), S2.coord.x(), S2.coord.y());
		QLineF SP(S.coord.x(), S.coord.y(), P.x(), P.y());
		qreal angle1 = SP.angleTo(SS1) < 180 ? SP.angleTo(SS1) : SS1.angleTo(SP);
		qreal angle2 = SP.angleTo(SS2) < 180 ? SP.angleTo(SS2) : SS2.angleTo(SP);
		//qreal angle = SS1.angleTo(SS2);
		qDebug() << angle1 + angle2;
		if (angle1 + angle2 >= 180)
			continue;
		bool notInside = true;
		for each (Side a in La2)
		{
			if (a.pHigh != S.id && a.pHigh != S1.id && a.pHigh != S2.id
				&& insideCircumCircle(getPointfromID(pts, a.pHigh).coord, S.coord, S1.coord, S2.coord))
			{
				notInside = false;
				break;
			}
			if (a.pLow != S.id && a.pLow != S1.id && a.pLow != S2.id
				&& insideCircumCircle(getPointfromID(pts, a.pLow).coord, S.coord, S1.coord, S2.coord))
			{
				notInside = false;
				break;
			}
		}
		if (notInside)
			break;
	}
	return i < N ? i : -1;
}

// Retourner l'ID du point mutuel entre 2 ar�tes, sinon -1
int mutualPoint(Side a1, Side a2)
{
	if (a1.pHigh == a2.pHigh || a1.pHigh == a2.pLow)
		return a1.pHigh;
	if (a1.pLow == a2.pHigh || a1.pLow == a2.pLow)
		return a1.pLow;
	return -1;
}
// Chercher les triangles orient�s incidents des ar�tes La1
// Retourner �galement la liste orient�e La2 des ar�tes incidentes aux triangles mais non � S
vector<Face> getIncidentFacesOriented(vector<Side> La1, vector<Face> faces, vector<Side> sides, vector<Side> &La2)
{
	vector<Face> fs;
	La1.push_back(La1[0]);
	for (int i = 0; i < La1.size() - 1; i++)
	{
		for (int j = 0; j < faces.size(); j++)
		{
			int foundS1 = -1;
			int foundS2 = -1;
			for (int k = 0; k < 3; k++)
			{
				if (faces[j].sidesID[k] == La1[i].id)
					foundS1 = k;
				if (faces[j].sidesID[k] == La1[i+1].id)
					foundS2 = k;
			}
			if (foundS1 != -1 && foundS2 != -1)
			{
				bool faceAlreadyExisted = false;
				for (int i = 0; i < fs.size(); i++)
				{
					if (faces[j].id == fs[i].id)
						faceAlreadyExisted = true;
				}
				if (!faceAlreadyExisted)
				{
					fs.push_back(faces[j]);
					for (int k = 0; k < 3; k++)
						if (k != foundS1 && k != foundS2)
						{
							La2.push_back(getSidefromID(sides, faces[j].sidesID[k]));
							break;
						}
				}
				break;
			}
		}
	}
	return fs;
}

// Chercher les ar�tes incidentes depuis le point P
vector<Side> getIncidentEdgesOriented(Point P, vector<Side> sides, vector<Point> pts)
{
	vector<Side> edges;
	// Chercher une ar�te incidente externe depuis P
	for (int i = 0; i < P.sides.size(); i++)
	{
		if (getSidefromID(sides, P.sides[i]).fRight == -1)
		{
			// Echanger la place de cette ar�te au d�but
			int temp = P.sides[0];
			P.sides[0] = P.sides[i];
			P.sides[i] = temp;
			break;
		}
	}
	// La premi�re ar�te sera choisie pour la r�f�rence de calcul des angles
	Side a = getSidefromID(sides, P.sides[0]);
	Point Q = getPointfromID(pts, a.pHigh == P.id ? a.pLow : a.pHigh);
	QLineF v(P.coord.x(), P.coord.y(), Q.coord.x(), Q.coord.y());
	vector<qreal> angles;
	angles.push_back(0); 
	edges.push_back(a);
	
	for (int i = 1; i < P.sides.size(); i++)
	{
		a = getSidefromID(sides, P.sides[i]);
		Q = getPointfromID(pts, a.pHigh == P.id ? a.pLow : a.pHigh);
		QLineF PQ(P.coord.x(), P.coord.y(), Q.coord.x(), Q.coord.y());
		qreal angle = v.angleTo(PQ);
		int k = 0;
		for (k = 0; k < angles.size(); k++)
			if (angle <= angles[k])
				break;
		angles.insert(angles.begin() + k, angle);
		edges.insert(edges.begin() + k, a);
	}
	return edges;
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
	if (F.sidesID.size() == 3)
	{
		int id1, id2, id3;
		id1 = getSidefromID(sides, F.sidesID[0]).pHigh;
		id2 = getSidefromID(sides, F.sidesID[0]).pLow;
		id3 = getSidefromID(sides, F.sidesID[1]).pHigh;
		pts_list.push_back(getPointfromID(pts, id1));
		pts_list.push_back(getPointfromID(pts, id2));
		if (id3 != id1 && id3 != id2)
			pts_list.push_back(getPointfromID(pts, id3));
		else
			pts_list.push_back(getPointfromID(pts, getSidefromID(sides, F.sidesID[1]).pLow));
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

QVector3D circumCircleCenter(QVector3D v1, QVector3D v2, QVector3D v3)
{
	float alpha = (v2 - v3).lengthSquared() * QVector3D::dotProduct(v1 - v2, v1 - v3) / (2 * QVector3D::crossProduct(v1 - v2, v2 - v3).lengthSquared());
	float beta = (v1 - v3).lengthSquared() * QVector3D::dotProduct(v2 - v1, v2 - v3) / (2 * QVector3D::crossProduct(v1 - v2, v2 - v3).lengthSquared());
	float theta = (v1 - v2).lengthSquared() * QVector3D::dotProduct(v3 - v1, v3 - v2) / (2 * QVector3D::crossProduct(v1 - v2, v2 - v3).lengthSquared());

	return alpha*v1 + beta*v2 + theta*v3;
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

bool isVisible(Side side, Point point)
{
	bool result;
	QVector3D temp;
	QVector3D normal = crossProductNormalized(side.points[0].coord, side.points[1].coord);
	QVector3D newVector(point.coord.x() - side.points[0].coord.x(), point.coord.y() - side.points[0].coord.y(), point.coord.z() - side.points[0].coord.z());
	float value = temp.dotProduct(normal, newVector);
	if (value < 0) {
		result = true;
	}
	else
	{
		result = false;
	}
	return result;
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
		// initialisation de angle_min et lmax avec le premier point d'indice diff�rent de i
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
		// mise � jour du pivot et du vecteur directeur
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

///-----------------------------------------------------------------------------------
vector<Face> TriangulationSimple(vector<Point> pts, list<Side> &_convexHull) {
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
				i++;
				break;
			}
		}

		// 2.b) Construction of the root Triangulation
		if (i < pts.size())
		{
			if (colinearPoints.size() >= 2)
			{
				unsigned int j = 0;
				/*for (j = 1; j < colinearPoints.size()+1; j++)
				{
					convexHull.push_back(Side(pts[j], pts[j - 1]));
					convexHull.push_back(Side(pts[j - 1], pts[i]));
					faces.push_back(Face(pts[j], pts[j - 1], pts[i]));
				}
				convexHull.push_back(Side(pts[i], pts[j - 1]));*/
				convexHull.push_back(Side(pts[0], pts[i]));
				for (j = 0; j < colinearPoints.size(); j++)
				{
					convexHull.push_back(Side(pts[j], pts[j+1]));
					//convexHull.push_back(Side(pts[j+1], pts[j+2]));
					faces.push_back(Face(pts[j], pts[j+1], pts[i]));
				}
				convexHull.push_back(Side(pts[i], pts[j])); // Pas (j,i) !
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

	list<Side>::iterator itE;
	list<Side>::iterator fromEdge = convexHull.end();
	list<Side>::iterator toEdge = convexHull.begin();

	_convexHull.clear();
	for (itE = convexHull.begin(); itE != convexHull.end(); ++itE)
	{
		_convexHull.push_back(*itE);
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

	if (viewedEdge.size() > 0) {
		std::vector<Side>::iterator itviewEdg;
		Point firstVert = viewedEdge.front().points[0];
		Point lastVert = viewedEdge.back().points[1];

		convexHull.insert(toEdge, Side(firstVert, nextVert));
		convexHull.insert(toEdge, Side(nextVert, lastVert));

		for (itviewEdg = viewedEdge.begin(); itviewEdg != viewedEdge.end(); itviewEdg++)
			convexHull.remove(*itviewEdg);
	}
	
	return viewedEdge;
}


float dot(QVector3D p, QVector3D op) { return p.x()*op.x() + p.y()*op.y() + p.z()*op.z(); } // Produit scalaire de 2 points
float norme(QVector3D p) { return float(sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z())); } // Norme du point (longueur du vecteur.)
//T distance(const Point<T>& p) { return sqrt(pow(p.x - x, 2) + pow(p.y - y, 2) + pow(p.z - z, 2)); } // Distance avec le vecteur en param�tre
//T angle(Point<T> p) { return acos(static_cast<T>(dot(p) / (norme() * p.norme()))); } // Return l'angle non orient� form� par les 2 vecteurs.
QVector3D crossProduct(QVector3D p, QVector3D op) { QVector3D t(p.y()*op.z() - p.z()*op.y(), p.z()*p.x() - p.x()*op.z(), p.x()*op.y() - p.y()*op.x()); return t; } //Produit vectoriel de 2 points.
QVector3D normalize(QVector3D p) { float n = norme(p); return QVector3D(p.x() / n, p.y() / n, p.z() / n);} //Normalisation du point.
QVector3D crossProductNormalized(QVector3D p, QVector3D op) 
{ 
	QVector3D final;
	final = crossProduct(p, op);
	final = normalize(final);
	return final;
} // Produit vectoriel normalis� de 2 points.

bool isEdgeViewed(QVector3D P, QVector3D A, QVector3D B, QVector3D n)
{
	QVector3D u = B - A; // Build a vector interior to AB
	QVector3D v = P - A;

	// Compute the normal interior to the Edge AB
	QVector3D normalAB = crossProductNormalized(u, n);

	return dot(normalAB, v) > 0 ? true : false;
}

///-----------------------------------------------------------------------------------
vector<Side> Fliping(vector<Face> faces)
{
	vector<Side> result;
	vector<Side> AcTemp;
	vector<Side> Ac;
	for (int i = 0; i < faces.size(); i++)
	{
		AcTemp.push_back(Side(faces[i].points[0], faces[i].points[1], i));
		AcTemp.push_back(Side(faces[i].points[1], faces[i].points[2], i));
		AcTemp.push_back(Side(faces[i].points[2], faces[i].points[0], i));
	}
	vector<int> sideToDestroy;
	for (int i = 0; i < AcTemp.size(); i++) {
		for (int j = 0; j < AcTemp.size(); j++) {
			if ((i != j)) {
				bool pass = false;
				for (int h = 0; h < sideToDestroy.size(); h++) {
					if (sideToDestroy[h] == i) {
						pass = true;
					}
				}
				if (!pass && ((AcTemp[i].points[0].coord == AcTemp[j].points[0].coord && AcTemp[i].points[1].coord == AcTemp[j].points[1].coord) 
					|| (AcTemp[i].points[1].coord == AcTemp[j].points[0].coord && AcTemp[i].points[0].coord == AcTemp[j].points[1].coord))) {
					AcTemp[i].idFace2 = AcTemp[j].idFace1;
					sideToDestroy.push_back(j);
				}
			}
		}
	}
	for (int i = 0; i < AcTemp.size(); i++) {
		bool pass = false;
		for (int h = 0; h < sideToDestroy.size(); h++) {
			if (sideToDestroy[h] == i) {
				pass = true;
			}
		}
		if(!pass)
			Ac.push_back(AcTemp[i]);
	}
	while (Ac.size()!=0)
	{
		Side A = Ac.back();
		Ac.pop_back();
		if (A.idFace2 != -1) {
			QVector3D p1;
			QVector3D p2;
			for (int i = 0; i <faces[A.idFace1].points.size(); i++) {
				if (faces[A.idFace1].points[i].coord != A.points[0].coord && faces[A.idFace1].points[i].coord != A.points[1].coord) {
					p1 = faces[A.idFace1].points[i].coord;
					break;
				}
			}
			for (int i = 0; i <faces[A.idFace2].points.size(); i++) {
				if (faces[A.idFace2].points[i].coord != A.points[0].coord && faces[A.idFace2].points[i].coord != A.points[1].coord) {
					p2 = faces[A.idFace2].points[i].coord;
					break;
				}
			}
			if (inCircumCircle(faces[A.idFace1], p2) || inCircumCircle(faces[A.idFace2], p1)) {
				QVector3D p3;
				QVector3D p4;
				for (int i = 0; i <faces[A.idFace1].points.size(); i++) {
					if (faces[A.idFace1].points[i].coord != p1 && faces[A.idFace1].points[i].coord != p2) {
						p3 = faces[A.idFace1].points[i].coord;
						break;
					}
				}
				for (int i = 0; i <faces[A.idFace2].points.size(); i++) {
					if (faces[A.idFace2].points[i].coord != p1 && faces[A.idFace2].points[i].coord != p2) {
						p4 = faces[A.idFace2].points[i].coord;
						break;
					}
				}
				faces[A.idFace1].points[0] = p1;
				faces[A.idFace1].points[1] = p2;
				faces[A.idFace1].points[2] = p3;
				faces[A.idFace2].points[0] = p1;
				faces[A.idFace2].points[1] = p2;
				faces[A.idFace2].points[2] = p4;
				A.points[0] = p1;
				A.points[1] = p2;
			}
		}
		result.push_back(A);
	}
	return result;
}

bool inCircumCircle(Face f, QVector3D v)
{
	float ab = (f.points[0].coord.x() * f.points[0].coord.x()) + (f.points[0].coord.y() * f.points[0].coord.y());
	float cd = (f.points[1].coord.x() * f.points[1].coord.x()) + (f.points[1].coord.y() * f.points[1].coord.y());
	float ef = (f.points[2].coord.x() * f.points[2].coord.x()) + (f.points[2].coord.y() * f.points[2].coord.y());

	float circum_x = (ab * (f.points[2].coord.y() - f.points[1].coord.y()) + cd * (f.points[0].coord.y() - f.points[2].coord.y()) + ef * (f.points[1].coord.y() - f.points[0].coord.y())) / (f.points[0].coord.x() * (f.points[2].coord.y() - f.points[1].coord.y()) + f.points[1].coord.x() * (f.points[0].coord.y() - f.points[2].coord.y()) + f.points[2].coord.x() * (f.points[1].coord.y() - f.points[0].coord.y())) / 2.f;
	float circum_y = (ab * (f.points[2].coord.x() - f.points[1].coord.x()) + cd * (f.points[0].coord.x() - f.points[2].coord.x()) + ef * (f.points[1].coord.x() - f.points[0].coord.x())) / (f.points[0].coord.y() * (f.points[2].coord.x() - f.points[1].coord.x()) + f.points[1].coord.y() * (f.points[0].coord.x() - f.points[2].coord.x()) + f.points[2].coord.y() * (f.points[1].coord.x() - f.points[0].coord.x())) / 2.f;
	float circum_radius = sqrtf(((f.points[0].coord.x() - circum_x) * (f.points[0].coord.x() - circum_x)) + ((f.points[0].coord.y() - circum_y) * (f.points[0].coord.y() - circum_y)));

	float dist = sqrtf(((v.x() - circum_x) * (v.x() - circum_x)) + ((v.y() - circum_y) * (v.y() - circum_y)));
	return dist <= circum_radius;
}

int getSideIDFromPoints(vector<Side> s, Point x, Point y) {
	for (int i = 0; i < s.size(); i++)
	{
		if ((s[i].pLow == x.id && s[i].pHigh == y.id) || (s[i].pHigh == x.id && s[i].pLow == y.id)) {
			return s[i].id;
		}
	}
	return -1;
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

///-----------------------------------------------------------------------------------
vector<Point> Voronoi(vector<Point> pts)
{
	vector<Face> tgs;
	list<Side> _convexHull;
	tgs = TriangulationSimple(pts, _convexHull);
	vector<Face> triangles = Fliping2(tgs);
	vector<Side> sides = Fliping(tgs);

	vector<QVector3D> centers;
	vector<Point> result;
	vector<QVector2D> alreadyChecked;

	for (int i = 0; i < triangles.size(); i++)
	{
		//centers.push_back(CircumCircleCenter(fp2[i]));
		QVector3D ccc = CircumCircleCenter(triangles[i]);
		//result.push_back(Point(ccc));

		for (int j = 0; j < sides.size(); j++)
		{
			int corr = 0;
			// Pour chaque point du triangle courant
			for (int h = 0; h < triangles[i].points.size(); h++)
			{
				if (sides[j].points[0].coord == triangles[i].points[h].coord || sides[j].points[1].coord == triangles[i].points[h].coord) {
					corr++;
				}
			}
			// Si on est bien sur un cot� du triangle
			if (corr == 2) {
				// Si on est sur un cot� exterieur
				if (sides[j].idFace2 == -1) {
					Point A = sides[j].points[0];
					Point B = sides[j].points[1];

					Point C;
					for (int y = 0; y < triangles[i].points.size(); y++)
					{
						if (A.coord != triangles[i].points[y].coord && B.coord != triangles[i].points[y].coord) {
							C = triangles[i].points[y];
						}
					}

					/*QVector3D u = B.coord - A.coord;
					QVector3D v = C.coord - A.coord;
					QVector3D normal = crossProductNormalized(u, v);
					normal = crossProduct(u, v);*/

					QVector3D AB = B.coord - A.coord;
					QVector3D AC = C.coord - A.coord;
					QVector3D N = QVector3D(AB.y(), -AB.x(), 0);
					float D = N.x()*AC.x() + N.y()*AC.y();
					if (D > 0)
						N = -N;

					result.push_back(Point(ccc));
					//QVector3D cc = (A.coord + B.coord) / 2;
					//QVector3D x = QVector3D(cc.x() - ccc.x(), cc.y() - ccc.y(), cc.z() - ccc.z());
					//result.push_back(Point(ccc+x*200));
					result.push_back(Point(ccc + N*10));
					//result.push_back(Point(ccc + normal*200));
					//result.push_back(Point(ccc));
				}
				// Si on est sur un cot� adjacent
				else {
					QVector2D link; link.setX(i); link.setY(sides[j].idFace2);
					if (std::find(alreadyChecked.begin(), alreadyChecked.end(), link) != alreadyChecked.end()) {
					}
					else {
						QVector3D cccBis = CircumCircleCenter(triangles[sides[j].idFace2]);
						result.push_back(Point(ccc));
						result.push_back(Point(cccBis));
						QVector2D _link; _link.setY(i); _link.setX(sides[j].idFace2);
						alreadyChecked.push_back(_link);
					}
				}
			}
		}

		//result.push_back(Point(ccc));
	}


	/*std::list<Edge>::iterator it = aretes.begin();
	for (; it != aretes.end(); ++it)
	{
		//Cas d'une ar�te interne
		if (it->T1() != NULL && it->T2() != NULL)
		{
			Point2D center1 = it->T1()->getCircumCircleCenter();
			voronoi.push_back(center1);
			Point2D center2 = it->T2()->getCircumCircleCenter();
			voronoi.push_back(center2);
		}
		//Cas d'une ar�te externe
		else if (it->T1() != NULL && it->T2() == NULL || it->T1() == NULL && it->T2() != NULL)
		{
			Triangle *triangle;
			if (it->T1() != NULL)
				triangle = it->T1();
			else
				triangle = it->T2();

			//on push le premier point (centre triangle)
			Point2D center1;
			center1 = triangle->getCircumCircleCenter();
			voronoi.push_back(center1);

			//calcul normal
			glm::vec2 normal = triangle->getNormal(&(*it));

			Point2D center2;
			glm::vec2 normalized = glm::normalize(normal);
			normalized *= -3000;
			normalized = it->GetCenter() + normalized;
			center2 = Point2D(normalized.x, normalized.y);

			voronoi.push_back(center2);
		}
	}*/

	return result;
}

vector<Face> Fliping2(vector<Face> faces)
{
	vector<Side> AcTemp;
	vector<Side> Ac;
	for (int i = 0; i < faces.size(); i++)
	{
		AcTemp.push_back(Side(faces[i].points[0], faces[i].points[1], i));
		AcTemp.push_back(Side(faces[i].points[1], faces[i].points[2], i));
		AcTemp.push_back(Side(faces[i].points[2], faces[i].points[0], i));
	}
	vector<int> sideToDestroy;
	for (int i = 0; i < AcTemp.size(); i++) {
		for (int j = 0; j < AcTemp.size(); j++) {
			if ((i != j)) {
				bool pass = false;
				for (int h = 0; h < sideToDestroy.size(); h++) {
					if (sideToDestroy[h] == i) {
						pass = true;
					}
				}
				if (!pass && ((AcTemp[i].points[0].coord == AcTemp[j].points[0].coord && AcTemp[i].points[1].coord == AcTemp[j].points[1].coord)
					|| (AcTemp[i].points[1].coord == AcTemp[j].points[0].coord && AcTemp[i].points[0].coord == AcTemp[j].points[1].coord))) {
					AcTemp[i].idFace2 = AcTemp[j].idFace1;
					sideToDestroy.push_back(j);
				}
			}
		}
	}
	for (int i = 0; i < AcTemp.size(); i++) {
		bool pass = false;
		for (int h = 0; h < sideToDestroy.size(); h++) {
			if (sideToDestroy[h] == i) {
				pass = true;
			}
		}
		if (!pass)
			Ac.push_back(AcTemp[i]);
	}
	while (Ac.size() != 0)
	{
		Side A = Ac.back();
		Ac.pop_back();
		if (A.idFace2 != -1) {
			QVector3D p1;
			int idP1;
			QVector3D p2;
			int idP2;
			for (int i = 0; i <faces[A.idFace1].points.size(); i++) {
				if (faces[A.idFace1].points[i].coord != A.points[0].coord && faces[A.idFace1].points[i].coord != A.points[1].coord) {
					p1 = faces[A.idFace1].points[i].coord;
					idP1 = i;
					break;
				}
			}
			for (int i = 0; i <faces[A.idFace2].points.size(); i++) {
				if (faces[A.idFace2].points[i].coord != A.points[0].coord && faces[A.idFace2].points[i].coord != A.points[1].coord) {
					p2 = faces[A.idFace2].points[i].coord;
					idP2 = i;
					break;
				}
			}
			if (inCircumCircle(faces[A.idFace1], p2) && inCircumCircle(faces[A.idFace2], p1)) {
				QVector3D p3;
				int idP3;
				QVector3D p4;
				int idP4;
				for (int i = 0; i <faces[A.idFace1].points.size(); i++) {
					if (faces[A.idFace1].points[i].coord != p1 && faces[A.idFace1].points[i].coord != p2) {
						p3 = faces[A.idFace1].points[i].coord;
						idP3 = i;
						break;
					}
				}
				for (int i = 0; i <faces[A.idFace2].points.size(); i++) {
					if (faces[A.idFace2].points[i].coord != p1 && faces[A.idFace2].points[i].coord != p2 && faces[A.idFace2].points[i].coord != p3) {
						p4 = faces[A.idFace2].points[i].coord;
						idP4 = i;
						break;
					}
				}
				faces[A.idFace1].points[0] = p1;
				faces[A.idFace1].points[1] = p2;
				faces[A.idFace1].points[2] = p3;
				faces[A.idFace2].points[0] = p1;
				faces[A.idFace2].points[1] = p2;
				faces[A.idFace2].points[2] = p4;
				A.points[0] = p1;
				A.points[1] = p2;
			}
		}
	}
	return faces;
}

QVector3D CircumCircleCenter(Face f)
{
	float ab = (f.points[0].coord.x() * f.points[0].coord.x()) + (f.points[0].coord.y() * f.points[0].coord.y());
	float cd = (f.points[1].coord.x() * f.points[1].coord.x()) + (f.points[1].coord.y() * f.points[1].coord.y());
	float ef = (f.points[2].coord.x() * f.points[2].coord.x()) + (f.points[2].coord.y() * f.points[2].coord.y());

	float circum_x = (ab * (f.points[2].coord.y() - f.points[1].coord.y()) + cd * (f.points[0].coord.y() - f.points[2].coord.y()) + ef * (f.points[1].coord.y() - f.points[0].coord.y())) / (f.points[0].coord.x() * (f.points[2].coord.y() - f.points[1].coord.y()) + f.points[1].coord.x() * (f.points[0].coord.y() - f.points[2].coord.y()) + f.points[2].coord.x() * (f.points[1].coord.y() - f.points[0].coord.y())) / 2.f;
	float circum_y = (ab * (f.points[2].coord.x() - f.points[1].coord.x()) + cd * (f.points[0].coord.x() - f.points[2].coord.x()) + ef * (f.points[1].coord.x() - f.points[0].coord.x())) / (f.points[0].coord.y() * (f.points[2].coord.x() - f.points[1].coord.x()) + f.points[1].coord.y() * (f.points[0].coord.x() - f.points[2].coord.x()) + f.points[2].coord.y() * (f.points[1].coord.x() - f.points[0].coord.x())) / 2.f;
	return QVector3D(circum_x, circum_y, 0);
}





/* Original points in the input. */
vec3 A[MAXN];

/* E[i][j] indicates which (up to two) other points combine with the edge i and
* j to make a face in the hull.  Only defined when i < j.
*/
struct twoset {
	void insert(int x) { (a == -1 ? a : b) = x; }
	bool contains(int x) { return a == x || b == x; }
	void erase(int x) { (a == x ? a : b) = -1; }
	int size() { return (a != -1) + (b != -1); }
	int a, b;
} E[MAXN][MAXN];

/* Compute the half plane {x : c^T norm < disc}
* defined by the three points A[i], A[j], A[k] where
* A[inside_i] is considered to be on the 'interior' side of the face. */
face make_face(int i, int j, int k, int inside_i) {
	E[i][j].insert(k); E[i][k].insert(j); E[j][k].insert(i);

	face f;
	f.I[0] = i; f.I[1] = j; f.I[2] = k;
	f.norm = (A[j] - A[i]) * (A[k] - A[i]);
	f.disc = f.norm.dot(A[i]);
	if (f.norm.dot(A[inside_i]) > f.disc) {
		f.norm = -f.norm;
		f.disc = -f.disc;
	}
	return f;
}

vector<face> convexHull3D(vector<Point> pts) {
	/*for (int i = 0; i < N; i++) {
		cin >> A[i].X[0] >> A[i].X[1] >> A[i].X[2];
	}*/
	for (int i = 0; i < pts.size(); i++) {
		A[i].X[0] = pts[i].coord.x();
		A[i].X[1] = pts[i].coord.y();
		A[i].X[2] = pts[i].coord.z();
	}

	/* Initially construct the hull as containing only the first four points. */
	face f;
	vector<face> faces;
	memset(E, -1, sizeof(E));
	for (int i = 0; i < 16; i+=4)
		for (int j = i + 4; j < 16; j+=4)
			for (int k = j + 4; k < 16; k+=4) {
				faces.push_back(make_face(i+10, j+5, k+7, 20 - i - j - k+10));
			}

	/* Now add a point into the hull one at a time. */
	for (int i = 4; i < pts.size(); i++) {
		/* Find and delete all faces with their outside 'illuminated' by this
		* point. */
		for (int j = 0; j < faces.size(); j++) {
			f = faces[j];
			if (f.norm.dot(A[i]) > f.disc) {
				E[f.I[0]][f.I[1]].erase(f.I[2]);
				E[f.I[0]][f.I[2]].erase(f.I[1]);
				E[f.I[1]][f.I[2]].erase(f.I[0]);
				faces[j--] = faces.back();
				faces.resize(faces.size() - 1);
			}
		}
		/* Now for any edge still in the hull that is only part of one face
		* add another face containing the new point and that edge to the hull. */
		int nfaces = faces.size();
		for (int j = 0; j < nfaces; j++) {
			f = faces[j];
			for (int a = 0; a < 3; a++) for (int b = a + 1; b < 3; b++) {
				int c = 3 - a - b;
				if (E[f.I[a]][f.I[b]].size() == 2) continue;
				faces.push_back(make_face(f.I[a], f.I[b], i, f.I[c]));
			}
		}
	}
	return faces;
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
*/
