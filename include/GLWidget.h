#pragma once
#include <QtOpenGL>
#include <QGLWidget>
#include <QVector2D>
#include <vector>

#include "GL/glu.h"
#include "AlgoMath.h"

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
	~GLWidget();
	void initializeGL();
	void resizeGL(int width, int height);
	void paintGL();
	void changeModeEnvelop(int mode) { modeEnvelop = mode; needUpdate = true; }
	void changeModeTriangulation(int mode);
	QString labelTimer[5];

protected:
	// Les événements Qt de souris et du clavier
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	virtual void keyPressEvent(QKeyEvent* e);

public slots:
	void timeOutSlot();
	// Fonctions pour mettre à jour les paramètres de l'UI
	void setVoronoi(int f) { showVoronoi = f == 0 ? false : true; }
	void setMovePoint(int m) { movePoint = m == 0 ? false : true; }
	void setGrid(int g) { showGrid = g == 0 ? false : true; }
	void setEnvelop3D(int e) { showEnvelop3D = e == 0 ? false : true; }
	// Réinitialiser les données
	void resetData();
	// Réinitialiser le caméra au paramètres par défaut
	void resetCamera();

signals:
	// Signal Qt pour mettre à jour les labels de Timers
	void labelChanged(int);

private:
	// Conversion de coordonnées d'écran à coordonnées de la scène OPENGL
	QVector3D convertXY(int X, int Y);
	// Chercher du point (dans la nuage existante) la plus proche de la souris
	int findNearestPoint(QPoint p);
	// Dessiner la grille et les axes
	void drawGridandAxes();
	// Dessiner des points
	void drawPoints(vector<Point> points, QVector3D color);
	// Dessiner des côtés à partir des couples de points
	void drawLinesFromPoints(vector<Point> pts);
	// Dessiner des côtés à partir de la structure Side
	void drawLinesFromSides(vector<Side> sides);
	// Dessiner des triangles
	void drawFaces(vector<Face> faces);
	// Dessiner des triangles en utilisant les IDs
	void drawFacesWithID(vector<Face> faces, bool stipple);
	// Dessiner d'un polygone à partir d'un ensemble des points
	void drawPoly(vector<Point> pts, QVector3D color, float width);
	// Animer automatiquement des points
	void movePoints(vector<Point> &pts);
	// Animer automatiquement des points en 3D
	void move3DPoints(vector<QVector3D> &pts);
	// Réinitialiser la structure des données de Delaunay à partir d'une nuage de points
	void recalculateDelaunay(vector<Point> pts);

	// Les fonctions pour l'envelope 3D
	QVector3D randomVector(float offset);
	void resetPoints3D();
	void drawConvexHull(vector<int> ids, vector<Point> pts);
	void drawConvexHull(vector<int> ids, vector<QVector3D> pts);
	void drawPointsch(vector<QVector3D> points);
	void drawConvexHullFromSides();
	QPoint mousePos;
	QTimer *t_Timer;

	// Les paramètres de caméra OPENGL
	float m_theta;	// Rotation x-axis
	float m_phi;	// Rotation  y-axis
	double range;
	float m_aspectRatio;
	bool mouseLook;
	QPoint rotValue;
	QPoint tmpMousePos;
	QPoint tmpRotValue;
	int screenW;
	int screenH;

	// Les données
	vector<Point> points;
	vector<Side> sides;
	vector<Face> faces;
	int pointSelected = -1;

	// Les paramètres de l'UI
	int modeEnvelop = 0;
	int modeTriangulation = 3;
	bool needUpdate = false;
	bool showVoronoi = false;
	bool movePoint = false;
	bool showEnvelop3D = false;
	bool showGrid = false;
};