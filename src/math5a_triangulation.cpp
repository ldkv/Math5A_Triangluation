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
	QObject::connect(shortcut, SIGNAL(activated()), glScene, SLOT(resetCamera()));
	// Ctrl + D pour reset le caméra
	shortcut = new QShortcut(QKeySequence("Ctrl+D"), this);
	QObject::connect(shortcut, SIGNAL(activated()), glScene, SLOT(resetData()));

	// Signal pour mettre à jour les lablels des Timers
	connect(glScene, SIGNAL(labelChanged(int)), this, SLOT(updateLabels(int)));

	// Connect signals
	connect(ui.rbJarvis, SIGNAL(clicked()), this, SLOT(modeEnvelop()));
	connect(ui.rbGrahamScan, SIGNAL(clicked()), this, SLOT(modeEnvelop()));
	connect(ui.rbNoneEnv, SIGNAL(clicked()), this, SLOT(modeEnvelop()));

	connect(ui.rbTriSimple, SIGNAL(clicked()), this, SLOT(modeTriangulation()));
	connect(ui.rbFlipping, SIGNAL(clicked()), this, SLOT(modeTriangulation())); 
	connect(ui.rbTriDelaunay, SIGNAL(clicked()), this, SLOT(modeTriangulation()));
	connect(ui.rbNoneTri, SIGNAL(clicked()), this, SLOT(modeTriangulation()));
	connect(ui.cbVoronoi, SIGNAL(stateChanged(int)), glScene, SLOT(setVoronoi(int)));
	connect(ui.cbMovePoint, SIGNAL(stateChanged(int)), glScene, SLOT(setMovePoint(int)));
	connect(ui.cbShowEnvelop3D, SIGNAL(stateChanged(int)), glScene, SLOT(setEnvelop3D(int)));
	connect(ui.cbShowGrid, SIGNAL(stateChanged(int)), glScene, SLOT(setGrid(int)));
	
	connect(ui.bResetData, SIGNAL(clicked()), glScene, SLOT(resetData()));
	connect(ui.bResetCam, SIGNAL(clicked()), glScene, SLOT(resetCamera()));
	connect(ui.bQuit, SIGNAL(clicked()), this, SLOT(quit()));
}

Math5A_Triangulation::~Math5A_Triangulation()
{
	delete[] glScene;
}

void Math5A_Triangulation::updateLabels(int label)
{
	switch (label)
	{
	case 0:
		ui.laTimeJarvis->setText(glScene->labelTimer[label]);
		break;
	case 1:
		ui.laTimeGraham->setText(glScene->labelTimer[label]);
		break;
	case 2:
		ui.laTimeTriSimple->setText(glScene->labelTimer[label]);
		break;
	case 3:
		ui.laTimeFlipping->setText(glScene->labelTimer[label]);
		break;
	case 4:
		ui.laTimeDelaunay->setText(glScene->labelTimer[label]);
		break;
	default:
		break;
	}
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
	else if (ui.rbFlipping->isChecked())
		glScene->changeModeTriangulation(2);
	else if (ui.rbTriDelaunay->isChecked())
		glScene->changeModeTriangulation(3);
	else
		glScene->changeModeTriangulation(0);
}

// Quitter
void Math5A_Triangulation::quit()
{
	qApp->quit();
}