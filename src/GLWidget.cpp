#include "stdafx.h"
#include "GLWidget.h"

#include <QDebug>
#include <QKeyEvent>
#include <math.h>
#include <time.h>
#include "AlgoMath.h"

// Les données pour l'envelope convexe 3D
vector<int> indices;
vector<QVector3D> point3D;
list<Side>  _convexHull;
ConvexHull ch;

// Initialisation de la scène OpenGL
GLWidget::GLWidget(QWidget *parent) :
	QOpenGLWidget(parent), 
	m_theta(180.0f),
	m_phi(0.0f),
	m_aspectRatio(1.0)
{
	int seconde = 1000; // 1 seconde = 1000 ms
	int timerInterval = seconde / 60;
	t_Timer = new QTimer(this);
	connect(t_Timer, SIGNAL(timeout()), this, SLOT(timeOutSlot()));
	t_Timer->start(timerInterval);
	setMouseTracking(true);
	setFocusPolicy(Qt::StrongFocus);
}

GLWidget::~GLWidget()
{
	resetData(); 
	delete[] t_Timer;
}

// Pour garder les FPS à 60
void GLWidget::timeOutSlot()
{
	update();
}

// Initialisation du module OpenGL
void GLWidget::initializeGL()
{
	glClearColor(0, 0, 0, 0.0f);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	range = 100.0;
	resetPoints3D();
}

// Redimensionner de la scène pour adapter à la fenêtre principale
void GLWidget::resizeGL(int width, int height)
{
	screenH = height;
	screenW = width;

	if (height == 0)
		height = 1;

	glViewport(0, 0, width, height);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	m_aspectRatio = double(width) / double(height);
	//gluOrtho2D(0, width, height, 0);
	if (width <= height)
		glOrtho(-range, range, -range / m_aspectRatio, range / m_aspectRatio, range, -range);
	else
		glOrtho(-range * m_aspectRatio, range * m_aspectRatio, -range, range, range, -range);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

// Fonction mettre à jour de la scène OpenGL
void GLWidget::paintGL()
{
	glClearColor(0, 0, 0, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glColor3f(1.0f, 0.0f, 0.0f);

	// Shape Marking
	glPushMatrix();
	glRotatef(m_theta, 1.0f, 0.0f, 0.0f);
	glRotatef(m_phi, 0.0f, 0.0f, 1.0f);

	// Afficher la grille et les Axes dans la scène
	if (showGrid)
		drawGridandAxes();
	
	// Affichage de la Triangulation
	QElapsedTimer timer;
	int time;
	switch (modeTriangulation)
	{
	case 1:	// Triangulation simple
		if (needUpdate)
		{
			faces.clear();
			timer.start();
			faces = TriangulationSimple(points, _convexHull);
			time = timer.nsecsElapsed() / 1000;
			labelTimer[2] = QString::number(time) + " us";
			emit labelChanged(2);
			needUpdate = false;
		}
		drawFaces(faces);
		break;
	case 2:	// Triangulation simple avec Flipping
		if (needUpdate)
		{
			sides.clear();
			faces.clear();
			timer.start();
			faces = TriangulationSimple(points, _convexHull);
			sides = Flipping(faces);
			time = timer.nsecsElapsed() / 1000;
			labelTimer[3] = QString::number(time) + " us";
			emit labelChanged(3);
			needUpdate = false;
		}
		drawLinesFromSides(sides);
		break;
	case 3:	// Triangulation Delaunay
		drawFacesWithID(faces, false);
		break;
	default:
		break;
	}

	// Affichage du diagramme de Voronoï
	if (showVoronoi)
	{
		vector<Point> voronoi;
		if (modeTriangulation != 3)
			recalculateDelaunay(points);
		voronoi = diagramVoronoi(points, sides, faces);
		drawLinesFromPoints(voronoi);
		drawPoints(voronoi, QVector3D(150.0f, 0, 0));
	}

	// Affichage de l'envelope convexe
	vector<Point> envelop;
	switch (modeEnvelop)
	{
	case 1:	// Enveloppe par méthode marche de Jarvis
		timer.start();
		envelop = EnvelopeJarvis(points);
		time = timer.nsecsElapsed() / 1000;
		labelTimer[0] = QString::number(time) + " us";
		emit labelChanged(0);
		drawPoly(envelop, QVector3D(150.0f, 150.0f, 150.0f), 3);
		break;
	case 2:	// Enveloppe par méthode Graham-Scan
		timer.start();
		envelop = GrahamScan(points);
		time = timer.nsecsElapsed() / 1000;
		labelTimer[1] = QString::number(time) + " us";
		emit labelChanged(1);
		drawPoly(envelop, QVector3D(150.0f, 0, 150.0f), 2);
		break;
	default:
		break;
	}

	// Animation automatique des points
	if (movePoint)
	{
		movePoints(points);
		if (modeTriangulation == 3)
			recalculateDelaunay(points);
		needUpdate = true;
	}
	drawPoints(points, QVector3D(255.0f, 255.0f, 255.0f));

	// Affichage de l'envelope convexe en 3D
	if (showEnvelop3D) 
	{
		if (modeTriangulation == 0) 
		{
			if (movePoint)
			{
				move3DPoints(point3D);
			}
			drawConvexHull(ch.calculate(point3D), point3D);
			drawPointsch(point3D);
		}
		else
		{
			vector<QVector3D > t;
			for (size_t i = 0; i < points.size(); i++)
			{
				points[i].coord.setZ((i + 1) * 2);
				t.push_back(points[i].coord);
			}
			drawConvexHull(ch.calculate(t), points);
		}
	}
	glPopMatrix();
}

// Conversion de coordonnées d'écran à coordonnées de la scène OPENGL
QVector3D GLWidget::convertXY(int X, int Y)
{
	return QVector3D((int)((float)X * 2.0 * range * m_aspectRatio / screenW - range * m_aspectRatio), (int)((float)Y * 2.0 * range / screenH - range), 0);
}

// Callback pour les click de la souris
void GLWidget::mousePressEvent(QMouseEvent *event)
{
	pointSelected = findNearestPoint(event->pos());
	QElapsedTimer timer;
	// Ajout d'un point
	if (event->buttons() & Qt::LeftButton)
	{
		if (pointSelected == -1) 
		{  
			// Mode Triangulation Delaunay
			if (modeTriangulation == 3) 
			{
				timer.start();
				Delaunay_addPoint(points, sides, faces, convertXY(event->pos().x(), event->pos().y()));
				int time = timer.nsecsElapsed() / 1000;
				labelTimer[4] = QString::number(time) + " us";
				emit labelChanged(4);
			}
			else
				points.push_back(Point(convertXY(event->pos().x(), event->pos().y())));
			needUpdate = true;
		}
	}
	// Suppression d'un point
	else if (event->buttons() & Qt::RightButton)
	{
		if (pointSelected != -1) 
		{
			// Mode Triangulation Delaunay
			if (modeTriangulation == 3)
			{
				timer.start();
				Delaunay_deletePoint(points, sides, faces, pointSelected);
				int time = timer.nsecsElapsed() / 1000;
				labelTimer[4] = QString::number(time) + " us";
				emit labelChanged(4);
			}
			else
				points.erase(points.begin() + pointSelected);
			needUpdate = true;
		}
	}
}

// Callback pour le mouvement de la souris
void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
	mousePos = event->pos();
	if (event->buttons() & Qt::LeftButton)
	{
		if (pointSelected >= 0)
		{
			// Mode Triangulation Delaunay
			if (modeTriangulation == 3)
			{
				QElapsedTimer timer;
				timer.start();
				Delaunay_deletePoint(points, sides, faces, pointSelected);
				Delaunay_addPoint(points, sides, faces, convertXY(event->pos().x(), event->pos().y()));
				pointSelected = points.size() - 1;
				int time = timer.nsecsElapsed() / 1000;
				labelTimer[4] = QString::number(time) + " us";
				emit labelChanged(4);
			}
			else
				points[pointSelected] = Point(convertXY(event->pos().x(), event->pos().y()));
			needUpdate = true;
		}
	}
}

// Callback pour la relâche de la souris
void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
	if (event->button() == Qt::LeftButton && pointSelected >= 0)
	{
		//update();
	}
}

// Chercher du point (dans la nuage existante) la plus proche de la souris
int GLWidget::findNearestPoint(QPoint p)
{
	int nbPoints = points.size();
	for (int i = 0; i < nbPoints; i++)
	{
		QVector3D d = convertXY(p.x(), p.y()) - QVector3D(points[i].coord.x(), points[i].coord.y(), 0);
		if (sqrt(pow(d.x(), 2) + pow(d.y(), 2) <= POINT_SIZE*4))
			return i;
	}
	return -1;
}

// Dessiner la grille et les axes dans la scène
void GLWidget::drawGridandAxes()
{
	// Grid
	glLineWidth(1);
	glBegin(GL_LINES);
	glColor3f(0.1, 0.1, 0.1);
	for (float x = -100; x < 100; x += 5)
	{
		glVertex3d(x, -100, 0);
		glVertex3d(x, 100, 0);
	}
	glEnd();
	glBegin(GL_LINES);
	glColor3f(0.1, 0.1, 0.1);
	for (float z = -100; z < 100; z += 5)
	{
		glVertex3d(-100, z, 0);
		glVertex3d(100, z, 0);
	}
	glEnd();

	// Axis
	glColor3f(0, 1, 0);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 50, 0);
	glEnd();

	glColor3f(1, 0, 0);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(50, 0, 0);
	glEnd();

	glColor3f(0, 0, 1);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 0, 50);
	glEnd();
}

// Dessiner des côtés à partir des couples de points
void  GLWidget::drawLinesFromPoints(vector<Point> pts)
{
	int nbSides = pts.size();
	if (nbSides == 0)
		return;
	glColor3f(150.0f, 0.0f, 0.0f);
	//glPushAttrib is done to return everything to normal after drawing
	glPushAttrib(GL_ENABLE_BIT);
	glLineStipple(10, 0xAAAA);
	//glEnable(GL_LINE_STIPPLE);
	glBegin(GL_LINES);
	for (int i = 0; i < nbSides; i+=2)
	{
		glVertex3f(pts[i].coord.x(), pts[i].coord.y(), pts[i].coord.z());
		glVertex3f(pts[i + 1].coord.x(), pts[i + 1].coord.y(), pts[i + 1].coord.z());
	}
	glEnd();
	glPopAttrib();
}

// Dessiner des côtés à partir de la structure Side
void GLWidget::drawLinesFromSides(vector<Side> sides)
{
	int nbSides = sides.size();
	if (nbSides == 0)
		return;
	glColor3f(255.0f, 255.0f, 0.0f);
	for (int i = 0; i < nbSides; i++)
	{
		glBegin(GL_LINES);
		glVertex3f(sides[i].points[0].coord.x(), sides[i].points[0].coord.y(), sides[i].points[0].coord.z());
		glVertex3f(sides[i].points[1].coord.x(), sides[i].points[1].coord.y(), sides[i].points[1].coord.z());
		glEnd();
	}
}

// Dessiner des triangles
void GLWidget::drawFaces(vector<Face> faces)
{
	int nbPoints = faces.size();
	if (nbPoints == 0)
		return;
	glColor3f(150.0f, 150.0f, 150.0f);
	for (int i = 0; i < faces.size(); i++)
	{
		if ((faces[i].points.size() == 3)) {
			glBegin(GL_LINES);
			glVertex3f(faces[i].points[0].coord.x(), faces[i].points[0].coord.y(), faces[i].points[0].coord.z());
			glVertex3f(faces[i].points[1].coord.x(), faces[i].points[1].coord.y(), faces[i].points[1].coord.z());
			glEnd();
			glBegin(GL_LINES);
			glVertex3f(faces[i].points[1].coord.x(), faces[i].points[1].coord.y(), faces[i].points[1].coord.z());
			glVertex3f(faces[i].points[2].coord.x(), faces[i].points[2].coord.y(), faces[i].points[2].coord.z());
			glEnd();
			glBegin(GL_LINES);
			glVertex3f(faces[i].points[2].coord.x(), faces[i].points[2].coord.y(), faces[i].points[2].coord.z());
			glVertex3f(faces[i].points[0].coord.x(), faces[i].points[0].coord.y(), faces[i].points[0].coord.z());
			glEnd();
		}
	}
}

// Dessiner des triangles en utilisant les IDs
void GLWidget::drawFacesWithID(vector<Face> faces, bool stipple)
{
	int nbFaces = faces.size();
	if (nbFaces == 0)
		return;
	glColor3f(0, 255.0f, 0);
	//glPushAttrib is done to return everything to normal after drawing
	glPushAttrib(GL_ENABLE_BIT);
	glLineStipple(2, 0x00FF);
	if (stipple)
		glEnable(GL_LINE_STIPPLE);
	for (int k = 0; k < faces.size(); k++)
	{
		Point pts[6];
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 3; i++)
		{
			pts[i * 2] = getPointfromID(points, getSidefromID(sides, faces[k].sidesID[i]).pLow);
			pts[i * 2 + 1] = getPointfromID(points, getSidefromID(sides, faces[k].sidesID[i]).pHigh);
			glVertex3f(pts[i * 2].coord.x(), pts[i * 2].coord.y(), pts[i * 2].coord.z());
			glVertex3f(pts[i * 2 + 1].coord.x(), pts[i * 2 + 1].coord.y(), pts[i * 2 + 1].coord.z());
		}
		glEnd();
	}
	glPopAttrib();
}

// Dessiner d'un polygone à partir d'un ensemble des points
void GLWidget::drawPoly(vector<Point> pts, QVector3D color, float width)
{
	int nbPoints = pts.size();
	if (nbPoints == 0)
		return;
	glColor3f(color.x(), color.y(), color.z());
	glLineWidth(width);
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < nbPoints; i++)
	{
		glVertex3f(pts[i].coord.x(), pts[i].coord.y(), pts[i].coord.z());
	}
	glEnd();
}

// Dessiner des points
void GLWidget::drawPoints(vector<Point> points, QVector3D color)
{
	int nbPoints = points.size();
	if (nbPoints == 0)
		return;
	glColor3f(color.x(), color.y(), color.z());
	glPointSize(POINT_SIZE);
	glBegin(GL_POINTS);
	for (int i = 0; i < nbPoints; i++)
		glVertex3f(points[i].coord.x(), points[i].coord.y(), points[i].coord.z());
	glEnd();
}


void GLWidget::keyPressEvent(QKeyEvent* e)
{
	switch (e->key())
	{
	case Qt::Key_Escape:
		QCoreApplication::instance()->quit();
		break;

	case Qt::Key_Left:
		m_phi += 2.0f;
		//update();
		break;

	case Qt::Key_Right:
		m_phi -= 2.0f;
		//update();
		break;

	case Qt::Key_Up:
		m_theta += 2.0f;
		//update();
		break;

	case Qt::Key_Down:
		m_theta -= 2.0f;
		//update();
		break;
	}
}

// Mettre à jour le mode de triangulation grâce à l'UI
void GLWidget::changeModeTriangulation(int mode)
{
	modeTriangulation = mode;
	if (mode == 3)
		recalculateDelaunay(points);
	needUpdate = true;
}

// Réinitialiser la structure des données de Delaunay à partir d'une nuage de points
void GLWidget::recalculateDelaunay(vector<Point> pts)
{
	resetData();
	QElapsedTimer timer;
	timer.start();
	for (int i = 0; i < pts.size(); i++)
		Delaunay_addPoint(points, sides, faces, pts[i].coord);
	int time = timer.nsecsElapsed() / 1000;
	labelTimer[4] = QString::number(time) + " us";
	emit labelChanged(4);
	needUpdate = true;
}

// Réinitialiser les données
void GLWidget::resetData() 
{
	points.clear();
	sides.clear();
	faces.clear();
	resetGlobalID();
	needUpdate = true;
}

// Réinitialiser le caméra au paramètres par défaut
void GLWidget::resetCamera() 
{
	m_theta = 180.0f;
	m_phi = 0.0f;
	QApplication::setOverrideCursor(Qt::PointingHandCursor);
}


float a = 0.0025f;
float b = -0.0010f;
float c = 0.0025f;
// Animer automatiquement les points en 2D
void GLWidget::movePoints(vector<Point> &pts) {
	for (int i = 0; i < pts.size(); i++)
	{
		float x_old = pts[i].coord.x(); float y_old = pts[i].coord.y(); float z_old = pts[i].coord.z();
		if (i % 2 == 0) {
			pts[i].coord.setX(x_old * cos(a) - y_old * sin(a));
			pts[i].coord.setY(x_old * sin(a) + y_old * cos(a));
		}
		else {
			pts[i].coord.setX(x_old * cos(b) - y_old * sin(b));
			pts[i].coord.setY(x_old * sin(b) + y_old * cos(b));
		}
	}
}

// Animer automatiquement les points en 3D
void GLWidget::move3DPoints(vector<QVector3D> &pts) {
	for (int i = 0; i < pts.size(); i++)
	{
		float x_old = pts[i].x(); float y_old = pts[i].y(); float z_old = pts[i].z();
		if (i % 2 == 0) {
			pts[i].setX(x_old * cos(a) - y_old * sin(a));
			pts[i].setY(x_old * sin(a) + y_old * cos(a));
			//pts[i].setZ(z_old * sin(a) + x_old * cos(a));
		}
		else {
			pts[i].setX(x_old * cos(b) - y_old * sin(b));
			pts[i].setY(x_old * sin(b) + y_old * cos(b));
			//pts[i].setZ(z_old * sin(b) + x_old * cos(b));
		}
	}
}

/****** LES FONCTIONS POUR DESSINER L'ENVELOPE CONVEXE 3D ******/
// A AJOUTER
QVector3D GLWidget::randomVector(float offset)
{
	float radius = (0.5f + rand() * 0.5) * offset;

	float theta = rand() * 2 * 3.14159;
	float phi = rand() * 3.14159;

	return QVector3D(cos(theta) * sin(phi) * radius*1.5,
		cos(phi) * radius * 0.1,
		sin(theta) * sin(phi) * radius * 15);
}

// A AJOUTER
void GLWidget::resetPoints3D()
{
	srand(time(NULL));
	point3D.clear();
	float size = 0.03;
	int i = 100;
	while (i--)
		point3D.push_back(randomVector(size));
	indices = ch.calculate(point3D);
}

// A AJOUTER
void GLWidget::drawConvexHullFromSides()
{
	int nbSides = _convexHull.size();
	if (nbSides == 0)
		return;
	glColor3f(0.0f, 150.0f, 0.0f);

	list<Side>::iterator itE;
	list<Side>::iterator fromEdge = _convexHull.end();
	list<Side>::iterator toEdge = _convexHull.begin();

	for (itE = _convexHull.begin(); itE != _convexHull.end(); itE++)
	{
		glBegin(GL_LINES);
		glVertex3f(itE->points[0].coord.x(), itE->points[0].coord.y(), itE->points[0].coord.z());
		glVertex3f(itE->points[1].coord.x(), itE->points[1].coord.y(), itE->points[1].coord.z());
		glEnd();
	}
}

// A AJOUTER
void GLWidget::drawConvexHull(vector<int> ids, vector<Point> pts) {
	int nbPoints = ids.size();
	if (nbPoints == 0)
		return;
	glColor3f(0.0f, 150.0f, 0.0f);
	for (int i = 0; i < ids.size(); i += 3)
	{
		glBegin(GL_LINES);
		glVertex3f(pts[ids[i]].coord[0], pts[ids[i]].coord[1], pts[ids[i]].coord[2]);
		glVertex3f(pts[ids[i + 1]].coord[0], pts[ids[i + 1]].coord[1], pts[ids[i + 1]].coord[2]);
		glEnd();
		glBegin(GL_LINES);
		glVertex3f(pts[ids[i + 1]].coord[0], pts[ids[i + 1]].coord[1], pts[ids[i + 1]].coord[2]);
		glVertex3f(pts[ids[i + 2]].coord[0], pts[ids[i + 2]].coord[1], pts[ids[i + 2]].coord[2]);
		glEnd();
		glBegin(GL_LINES);
		glVertex3f(pts[ids[i + 2]].coord[0], pts[ids[i + 2]].coord[1], pts[ids[i + 2]].coord[2]);
		glVertex3f(pts[ids[i]].coord[0], pts[ids[i]].coord[1], pts[ids[i]].coord[2]);
		glEnd();
	}
}

// A AJOUTER
void GLWidget::drawConvexHull(vector<int> ids, vector<QVector3D> pts) {
	int nbPoints = ids.size();
	if (nbPoints == 0)
		return;
	glColor3f(0.0f, 150.0f, 0.0f);
	for (int i = 0; i < ids.size(); i += 3)
	{
		glBegin(GL_LINES);
		glVertex3f(pts[ids[i]].x(), pts[ids[i]].y(), pts[ids[i]].z());
		glVertex3f(pts[ids[i + 1]].x(), pts[ids[i + 1]].y(), pts[ids[i + 1]].z());
		glEnd();
		glBegin(GL_LINES);
		glVertex3f(pts[ids[i + 1]].x(), pts[ids[i + 1]].y(), pts[ids[i + 1]].z());
		glVertex3f(pts[ids[i + 2]].x(), pts[ids[i + 2]].y(), pts[ids[i + 2]].z());
		glEnd();
		glBegin(GL_LINES);
		glVertex3f(pts[ids[i + 2]].x(), pts[ids[i + 2]].y(), pts[ids[i + 2]].z());
		glVertex3f(pts[ids[i]].x(), pts[ids[i]].y(), pts[ids[i]].z());
		glEnd();
	}
}

// A AJOUTER
void GLWidget::drawPointsch(vector<QVector3D> points)
{
	int nbPoints = points.size();
	if (nbPoints == 0)
		return;
	glColor3f(1, 0, 1);
	glPointSize(POINT_SIZE);
	glBegin(GL_POINTS);
	for (int i = 0; i < nbPoints; i++)
		glVertex3f(points[i].x(), points[i].y(), points[i].z());
	glEnd();
}