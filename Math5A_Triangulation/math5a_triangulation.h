#ifndef MATH5A_TRIANGULATION_H
#define MATH5A_TRIANGULATION_H

#include <QtWidgets/QMainWindow>
#include "ui_math5a_triangulation.h"
#include "GLWidget.h"

class Math5A_Triangulation : public QMainWindow
{
	Q_OBJECT

public:
	Math5A_Triangulation(QWidget *parent = 0);
	~Math5A_Triangulation();

private slots:

private:
	Ui::triangulationForm ui;
	GLWidget *glScene;
};

#endif // MATH5A_TRIANGULATION_H
