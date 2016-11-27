#include "stdafx.h"
#include "math5a_triangulation.h"

Math5A_Triangulation::Math5A_Triangulation(QWidget *parent)
	: QMainWindow(parent), ui(Ui::triangulationForm())
{
	ui.setupUi(this);
	glScene = new GLWidget(this);
	ui.centralLayout->addWidget(glScene);

	// Créer les raccourcis clavier
	// Ctrl + Q pour quitter
	QShortcut *shortcut = new QShortcut(QKeySequence("Ctrl+Q"), this);
	QObject::connect(shortcut, SIGNAL(activated()), this, SLOT(quit()));
	// Ctrl + R pour reset le caméra
	shortcut = new QShortcut(QKeySequence("Ctrl+R"), this);
	QObject::connect(shortcut, SIGNAL(activated()), this, SLOT(resetCamera()));

	// Connect signals
	connect(ui.rbJarvis, SIGNAL(clicked()), this, SLOT(modeEnvelop()));
	connect(ui.rbGrahamScan, SIGNAL(clicked()), this, SLOT(modeEnvelop()));
	connect(ui.rbNoneEnv, SIGNAL(clicked()), this, SLOT(modeEnvelop()));

	connect(ui.rbTriSimple, SIGNAL(clicked()), this, SLOT(modeTriangulation()));
	connect(ui.rbTriDelaunay, SIGNAL(clicked()), this, SLOT(modeTriangulation()));
	connect(ui.rbVoronoi, SIGNAL(clicked()), this, SLOT(modeTriangulation()));
	connect(ui.rbNoneTri, SIGNAL(clicked()), this, SLOT(modeTriangulation()));
	connect(ui.cbFlipping, SIGNAL(stateChanged(int)), glScene, SLOT(setFlipping(int)));

	connect(ui.bResetData, SIGNAL(clicked()), glScene, SLOT(resetData()));
	connect(ui.bResetCam, SIGNAL(clicked()), glScene, SLOT(resetCamera()));
	connect(ui.bQuit, SIGNAL(clicked()), this, SLOT(quit()));
}

Math5A_Triangulation::~Math5A_Triangulation()
{
	delete[] glScene;
}

void Math5A_Triangulation::modeEnvelop() 
{
	if (ui.rbJarvis->isChecked())
	{
		glScene->changeModeEnvelop(1);
	}
	else if (ui.rbGrahamScan->isChecked())
	{
		glScene->changeModeEnvelop(2);
	}
	else
		glScene->changeModeEnvelop(0);
}

void Math5A_Triangulation::modeTriangulation()
{
	if (ui.rbTriSimple->isChecked())
		glScene->changeModeTriangulation(1);
	else if (ui.rbTriDelaunay->isChecked())
		glScene->changeModeTriangulation(2);
	else if (ui.rbVoronoi->isChecked())
		glScene->changeModeTriangulation(3);
	else
		glScene->changeModeTriangulation(0);
}

// Quitter
void Math5A_Triangulation::quit()
{
	qApp->quit();
}