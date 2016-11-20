#include "stdafx.h"
#include "GLWidget.h"

#include <QDebug>
#include <QKeyEvent>
#include <math.h>
#include "AlgoMath.h"

extern vector<Point> points;

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

	range = 100.0;
	m_aspectRatio = double(width) / double(height);
	//gluOrtho2D(0, width, height, 0);
	if (width <= height)
		glOrtho(-range, range, -range / m_aspectRatio, range / m_aspectRatio, range, -range);
	else
		glOrtho(-range * m_aspectRatio, range * m_aspectRatio, -range, range, range, -range);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	qDebug() << width << " " << height;
}

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

	//drawGridandAxis();

	// Points
	drawLines(TriangulationSimple(points));
	//drawLinesStrip(EnvelopeJarvis(points));
	//drawPoly(GrahamScan(points));
	drawPoints(points);

	glPopMatrix();
}

// Callback pour les click de la souris
void GLWidget::mousePressEvent(QMouseEvent *event)
{
	pointSelected = findNearestPoint(event->pos());
	if (event->buttons() & Qt::LeftButton)
	{
		if (pointSelected == -1) {
			points.push_back(Point(QVector3D((int)(((float)event->pos().x() - screenW / 2) / 4.55), (int)(((float)event->pos().y() - screenH / 2) / 4.55), 0)));
			//points.push_back(Point(QVector3D((int)((float)event->pos().x() * 2.0 * range / screenW - range), (int)((float)-event->pos().y() * 2.0 * range / screenH + range), 0)));
			qDebug() << ((float)event->pos().x() - screenW / 2) / 4.55 << " " << ((float)event->pos().y() - screenH / 2) / 4.55;

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
			points[pointSelected] = Point(QVector3D((int)(((float)event->pos().x() - screenW / 2) / 4.55), (int)(((float)event->pos().y() - screenH / 2) / 4.55), 0));
			update();
		}
	}

	// is the middle mouse button down?
	if (event->buttons() == Qt::MidButton)
	{
		// was it already down? If not, store the current coords
		if (!mouseLook)
		{
			tmpMousePos = event->pos();
			tmpRotValue = rotValue;
			mouseLook = true;
		}

		// update the rotation values depending on the relative mouse position
		rotValue.setX(tmpRotValue.x() + (tmpMousePos.x() - event->pos().x()) * 0.2);
		rotValue.setY(tmpRotValue.y() + (tmpMousePos.y() - event->pos().y()) * -0.2);
	}
	else
	{
		// turn off mouse look
		mouseLook = false;
	}
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
		QPoint d = QPoint((p.x() - screenW / 2) / 4.55, (p.y()-screenH / 2) / 4.55) - QPoint(points[i].coord.x(), points[i].coord.y());
		if (sqrt(pow(d.x(), 2) + pow(d.y(), 2) <= POINT_SIZE*4))
			return i;
	}
	return -1;
}

void GLWidget::drawGridandAxis()
{
	// Grid

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

void GLWidget::drawLinesStrip(vector<QVector3D> pts)
{
	int nbPoints = pts.size();
	if (nbPoints == 0)
		return;
	glColor3f(150.0f, 150.0f, 0);
	glLineWidth(6);
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < nbPoints; i++) 
	{
		glVertex3f(pts[i].x(), pts[i].y(), pts[i].z());
	}
	glEnd();
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
void GLWidget::drawLines(vector<Face> faces)
{
	int nbPoints = faces.size();
	if (nbPoints == 0)
		return;
	glColor3f(150.0f, 150.0f, 150.0f);
	for (int i = 0; i < faces.size(); i++)
	{
		if ((faces[i].points.size() == 3)) {
			glBegin(GL_LINES);
			/*glVertex3f(sides[i].points[0].coord.x() , sides[i].points[0].coord.y(), sides[i].points[0].coord.z());
			glVertex3f(sides[i].points[1].coord.x(), sides[i].points[1].coord.y(), sides[i].points[1].coord.z());*/
			glVertex3f(faces[i].points[0].coord.x(), faces[i].points[0].coord.y(), faces[i].points[0].coord.z());
			glVertex3f(faces[i].points[1].coord.x(), faces[i].points[1].coord.y(), faces[i].points[1].coord.z());
			/*for (int i = 0; i <= 2; i++)
			{
			glVertex3f(faces[i].points[i].coord.x(), faces[i].points[i].coord.y(), faces[i].points[i].coord.z());
			}*/
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

void GLWidget::drawPoly(vector<Point> points)
{
	int nbPoints = points.size();
	if (nbPoints == 0)
		return;
	glColor3f(150.0f, 0.0f, 150.0f);
	glLineWidth(2);
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < nbPoints; i++) {
		glVertex3f(points[i].coord.x(), points[i].coord.y(), points[i].coord.z());
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

void GLWidget::resetCamera() {
	m_theta = 180.0f;
	m_phi = 0.0f;
}