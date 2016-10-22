#include "stdafx.h"
#include "math5a_triangulation.h"

Math5A_Triangulation::Math5A_Triangulation(QWidget *parent)
	: QMainWindow(parent), ui(Ui::triangulationForm())
{
	ui.setupUi(this);
	glScene = new GLWidget(this);
	ui.centralLayout->addWidget(glScene);
	connect(ui.resetCam, SIGNAL(clicked()), this, SLOT(resetCamera()));
}

Math5A_Triangulation::~Math5A_Triangulation()
{

}

void Math5A_Triangulation::resetCamera()
{
	glScene->resetCamera();
	QApplication::setOverrideCursor(Qt::PointingHandCursor);
}
