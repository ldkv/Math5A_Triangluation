#pragma once
#include <QtOpenGL>
#include "GL/glu.h"

#include "AlgoMath.h"

#include <vector>

using namespace std;

#define POINT_SIZE 10




#include <QGLWidget>
#include <QVector2D>

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

	void resetCamera();

protected:
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	int findNearestPoint(QPoint p);

	virtual void keyPressEvent(QKeyEvent* e);

	// Dessiner la grille et les axes
	void drawGridAxes();
	// Dessiner des côtés à partir des points
	void drawPoints(vector<Point> points);
	void drawLines(vector<Point> points, vector<Side> sides);
	void drawPoly(vector<Point> pts, QVector3D color, float width);

public slots:
	void timeOutSlot();

signals:
	void MouseMoved();

private:
	QPoint mousePos;
	QTimer *t_Timer;
	QColor bgColor;

	float m_theta; // Rotation x-axis
	float m_phi; // Rotation  y-axis
	float m_aspectRatio;
	bool mouseLook;
	QPoint rotValue;
	QPoint tmpMousePos;
	QPoint tmpRotValue;

	int screenW;
	int screenH;
	double range;

	vector<Point> points;
	vector<Side> sides;
	vector<Face> faces;
	unsigned int current_id_points = 0;
	unsigned int current_id_sides = 0;
	unsigned int current_id_faces = 0;

	int pointSelected;
};





