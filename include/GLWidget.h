#pragma once
#include <QtOpenGL>
#include "GL/glu.h"

#include "AlgoMath.h"

#include <QGLWidget>
#include <QVector2D>
#include <vector>

using namespace std;

#define POINT_SIZE 10

const float pi = 3.141592653f;
const float twoPi = 2.0f * pi;
const float piBy2 = 0.5f * pi;
const float degToRad = pi / 180.0f;
const float radToDeg = 180.0f / pi;


class GLWidget : public QOpenGLWidget
{
	Q_OBJECT
public:
	GLWidget(QWidget *parent);
	void initializeGL();
	void resizeGL(int width, int height);
	void paintGL();
	~GLWidget();

	
	void changeModeEnvelop(int mode) { modeEnvelop = mode; }
	void changeModeTriangulation(int mode);

protected:
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	int findNearestPoint(QPoint p);

	virtual void keyPressEvent(QKeyEvent* e);

	void drawGridandAxes();
	QVector3D convertXY(int X, int Y);
	// Dessiner des c�t�s � partir des points
	void drawPoints(vector<Point> points, QVector3D color);
	//void drawLines(vector<Point> points, vector<Side> sides);
	void drawLinesFromSides(vector<Side> sides);
	void drawLinesFromPoints(vector<Point> pts);
	void drawFacesWithID(vector<Face> faces, bool stipple);
	void drawFaces(vector<Face> faces);
	void drawPoly(vector<Point> pts, QVector3D color, float width);
	void movePoints(vector<Point> &pts);
	void move3DPoints(vector<QVector3D> &pts);
	void recalculateDelaunay(vector<Point> pts);

public slots:
	void timeOutSlot();
	void setVoronoi(int f) { showVoronoi = f == 0 ? false : true; }
	void setMovePoint(int m) { movePoint = m == 0 ? false : true; }
	void setGrid(int g) { showGrid = g == 0 ? false : true; }
	void setEnvelop3D(int e) { showEnvelop3D = e == 0 ? false : true; qDebug() << "here"; }
	void resetData();
	void resetCamera();

signals:
	void MouseMoved();

private:
	QPoint mousePos;
	QTimer *t_Timer;
	QColor bgColor;

	float m_theta; // Rotation x-axis
	float m_phi; // Rotation  y-axis
	double range;
	float m_aspectRatio;
	bool mouseLook;
	QPoint rotValue;
	QPoint tmpMousePos;
	QPoint tmpRotValue;

	int screenW;
	int screenH;

	vector<Point> points;
	vector<Side> sides;
	vector<Face> faces;
	unsigned int current_id_points = 0;
	unsigned int current_id_sides = 0;
	unsigned int current_id_faces = 0;

	int pointSelected;

	int modeEnvelop = 0;
	int modeTriangulation = 3;
	bool showVoronoi = false;
	bool movePoint = false;
	bool showEnvelop3D = false;
	bool showGrid = false;
};





