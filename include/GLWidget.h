#pragma once
#include <QtOpenGL>
#include "GL/glu.h"

#include "AlgoMath.h"

using namespace std;

#define POINT_SIZE 10

class GLWidget : public QOpenGLWidget
{
	Q_OBJECT
public:
	GLWidget(QWidget *parent);
	void initializeGL();
	void resizeGL(int width, int height);
	void paintGL();
	~GLWidget();

protected:
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);

public slots:
	void timeOutSlot();

signals:
	void MouseMoved();

private:
	QPoint mousePos;
	QTimer *t_Timer;
	QColor bgColor;
};

