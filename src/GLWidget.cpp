#include "stdafx.h"
#include "GLWidget.h"

#include <QDebug>
#include <QKeyEvent>
#include <math.h>
#include "AlgoMath.h"

vector<Point> points;

list<Side>  _convexHull;

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
	glClearColor(bgColor.red() / 255.0f, bgColor.green() / 255.0f, bgColor.blue() / 255.0f, 0.0f);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	range = 100.0;
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

	//qDebug() << width << " " << height;
	//qDebug() << -range * m_aspectRatio << " " << range * m_aspectRatio;
}

void drawConvexHullFromSides();

void drawConvexHull3D(vector<face> faces);

// Fonction mettre à jour de la scène OpenGL
void GLWidget::paintGL()
{
	glClearColor(bgColor.red() / 255.0f, bgColor.green() / 255.0f, bgColor.blue() / 255.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glColor3f(1.0f, 0.0f, 0.0f);

	// Shape Marking
	glPushMatrix();
	glRotatef(m_theta, 1.0f, 0.0f, 0.0f);
	glRotatef(m_phi, 0.0f, 0.0f, 1.0f);

	drawGridandAxes();

	// Triangulation
	vector<Face> tgs;
	switch (modeTriangulation)
	{
	case 1:	// Triangulation simple (avec flipping ou non)
		tgs = TriangulationSimple(points, _convexHull);
		//drawFaces(Fliping2(tgs));
		drawFaces(tgs);
		if (flipping)
		{
			drawLinesFromSides(Fliping(tgs));
		}
		//drawLinesFromPoints(Voronoi(points));
		//drawPoints(Voronoi(points));
		//drawConvexHullFromSides();
		break;
	case 2:	// Triangulation Delaunay
		drawFaces(faces);
		break;
	case 3:	// Voronoi
		drawLinesFromPoints(Voronoi(points));
		drawPoints(Voronoi(points));
		break;
	default:
		break;
	}

	// Enveloppe
	switch (modeEnvelop)
	{
	case 1:	// Enveloppe par méthode marche de Jarvis
		drawPoly(EnvelopeJarvis(points), QVector3D(150.0f, 150.0f, 150.0f), 6);
		break;
	case 2:	// Enveloppe par méthode Graham-Scan
		drawPoly(GrahamScan(points), QVector3D(150.0f, 0, 150.0f), 2);
		break;
	default:
		break;
	}

	//movePoints(points);

	drawPoints(points);


	drawConvexHull3D(convexHull3D(points));

	glPopMatrix();
}

QVector3D GLWidget::convertXY(int X, int Y)
{
	return QVector3D((int)((float)X * 2.0 * range * m_aspectRatio / screenW - range * m_aspectRatio), (int)((float)Y * 2.0 * range / screenH - range), 0);
}

// Callback pour les click de la souris
void GLWidget::mousePressEvent(QMouseEvent *event)
{
	pointSelected = findNearestPoint(event->pos());
	if (event->buttons() & Qt::LeftButton)
	{
		if (pointSelected == -1) {
			points.push_back(Point(convertXY(event->pos().x(), event->pos().y())));
			//Delaunay_addPoint(points, sides, faces, convertXY(event->pos().x(), event->pos().y()));
			//qDebug() <<;
			update();
		}
	}
	else if (event->buttons() & Qt::RightButton)
	{
		if (pointSelected != -1) {
			points.erase(points.begin()+pointSelected);
			update();
		}
	}
}

// Callback pour le mouvement de la souris
void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
	mousePos = event->pos();
	emit MouseMoved();
	if (event->buttons() & Qt::LeftButton)
	{
		if (pointSelected >= 0)
		{
			points[pointSelected] = Point(convertXY(event->pos().x(), event->pos().y()));
			update();
		}
	}

	/*if (event->buttons() == Qt::MidButton)
	{
		if (!mouseLook)
		{
			tmpMousePos = event->pos();
			tmpRotValue = rotValue;
			mouseLook = true;
		}


		rotValue.setX(tmpRotValue.x() + (tmpMousePos.x() - event->pos().x()) * 0.2);
		rotValue.setY(tmpRotValue.y() + (tmpMousePos.y() - event->pos().y()) * -0.2);
	}
	else
	{
		// turn off mouse look
		mouseLook = false;
	}*/
}

// Callback pour la relâche de la souris
void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
	if (event->button() == Qt::LeftButton && pointSelected >= 0)
	{
		update();
	}
}

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

void GLWidget::drawFacesWithID(vector<Face> faces)
{
	int nbFaces = faces.size();
	if (nbFaces == 0)
		return;
	glColor3f(150.0f, 150.0f, 150.0f);
	for (int k = 0; k < faces.size(); k++)
	{
		Point pts[6];
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 3; i++)
		{
			pts[i*2] = getPointfromID(points, getSidefromID(sides, faces[k].sidesID[i]).pLow);
			pts[i*2+1] = getPointfromID(points, getSidefromID(sides, faces[k].sidesID[i]).pHigh);
			glVertex3f(pts[i * 2].coord.x(), pts[i * 2].coord.y(), pts[i * 2].coord.z());
			glVertex3f(pts[i * 2 + 1].coord.x(), pts[i * 2 + 1].coord.y(), pts[i * 2 + 1].coord.z());
		}
		glEnd();
	}
}

// Dessiner des côtés à partir des points
/*void GLWidget::drawLines(vector<Point> points, vector<Side> sides)
{
	int nbPoints = sides.size();
	if (nbPoints == 0)
		return;
	glColor3f(150.0f, 150.0f, 150.0f);
	glBegin(GL_LINES);
	for (int i = 0; i < nbPoints; i++) 
	{
		int indexP1 = getPointIndex(points, sides[i].pLow);
		int indexP2 = getPointIndex(points, sides[i].pHigh);
		glVertex3f(points.at(indexP1).coord.x(), points.at(indexP1).coord.y(), points.at(indexP1).coord.z());
		glVertex3f(points.at(indexP2).coord.x(), points.at(indexP2).coord.y(), points.at(indexP2).coord.z());
	}
	glEnd();
}*/

void GLWidget::drawLinesFromSides(vector<Side> sides)
{
	int nbSides = sides.size();
	if (nbSides == 0)
		return;
	glColor3f(0.0f, 0.0f, 150.0f);
	for (int i = 0; i < nbSides; i++)
	{
		glBegin(GL_LINES);
		glVertex3f(sides[i].points[0].coord.x(), sides[i].points[0].coord.y(), sides[i].points[0].coord.z());
		glVertex3f(sides[i].points[1].coord.x(), sides[i].points[1].coord.y(), sides[i].points[1].coord.z());
		glEnd();
	}
}

void  GLWidget::drawLinesFromPoints(vector<Point> pts)
{
	int nbSides = pts.size();
	if (nbSides == 0)
		return;
	glColor3f(150.0f, 150.0f, 0.0f);
	glBegin(GL_LINES);
	for (int i = 0; i < nbSides; i+=2)
	{
		glVertex3f(pts[i].coord.x(), pts[i].coord.y(), pts[i].coord.z());
		glVertex3f(pts[i + 1].coord.x(), pts[i + 1].coord.y(), pts[i + 1].coord.z());
	}
	glEnd();
}

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
void GLWidget::drawPoints(vector<Point> points)
{
	int nbPoints = points.size();
	if (nbPoints == 0)
		return;
	glColor3f(1, 1, 1);
	glPointSize(POINT_SIZE);
	glBegin(GL_POINTS);
	for (int i = 0; i < nbPoints; i++)
		glVertex3f(points[i].coord.x(), points[i].coord.y(), points[i].coord.z());
	glEnd();
}

void drawConvexHullFromSides()
{
	int nbSides = _convexHull.size();
	if (nbSides == 0)
		return;
	glColor3f(0.0f, 150.0f, 0.0f);

	list<Side>::iterator itE;
	list<Side>::iterator fromEdge = _convexHull.end();
	list<Side>::iterator toEdge = _convexHull.begin();

	for (itE = _convexHull.begin(); itE != _convexHull.end();itE++)
	{
		glBegin(GL_LINES);
		glVertex3f(itE->points[0].coord.x(), itE->points[0].coord.y(), itE->points[0].coord.z());
		glVertex3f(itE->points[1].coord.x(), itE->points[1].coord.y(), itE->points[1].coord.z());
		glEnd();
	}
}

void drawConvexHull3D(vector<face> faces)
{
	int nbPoints = faces.size();
	if (nbPoints == 0)
		return;
	glColor3f(150.0f, 150.0f, 150.0f);
	for (int i = 0; i < faces.size(); i++)
	{
			glBegin(GL_LINES);
			glVertex3f(faces[i].I[0], faces[i].I[0], faces[i].I[0]);
			glVertex3f(faces[i].I[1], faces[i].I[1], faces[i].I[1]);
			glEnd();
			glBegin(GL_LINES);
			glVertex3f(faces[i].I[1], faces[i].I[1], faces[i].I[1]);
			glVertex3f(faces[i].I[2], faces[i].I[2], faces[i].I[2]);
			glEnd();
			glBegin(GL_LINES);
			glVertex3f(faces[i].I[2], faces[i].I[2], faces[i].I[2]);
			glVertex3f(faces[i].I[0], faces[i].I[0], faces[i].I[0]);
			glEnd();
	}
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
		update();
		break;

	case Qt::Key_Right:
		m_phi -= 2.0f;
		update();
		break;

	case Qt::Key_Up:
		m_theta += 2.0f;
		update();
		break;

	case Qt::Key_Down:
		m_theta -= 2.0f;
		update();
		break;
	}
}

void GLWidget::resetData() 
{
	points.clear();
	sides.clear();
	faces.clear();
	resetGlobalID();
}

void GLWidget::resetCamera() {
	m_theta = 180.0f;
	m_phi = 0.0f;
	QApplication::setOverrideCursor(Qt::PointingHandCursor);
}

float a = 0.005f;
float b = -0.0025f;

void GLWidget::movePoints(vector<Point> &pts) {
	for (int i = 0; i < pts.size(); i++)
	{
		float x_old = pts[i].coord.x(); float y_old = pts[i].coord.y();
		if (i % 2 == 0) {
			pts[i].coord.setX(x_old * cos(a) - y_old * sin(a));
			pts[i].coord.setY(x_old * sin(a) + y_old * cos(a));
		}
		else {
			pts[i].coord.setX(x_old * cos(b) - y_old * sin(b));
			pts[i].coord.setY(x_old * sin(b) + y_old * cos(b));
		}

	}
	//a += 0.005f;
}