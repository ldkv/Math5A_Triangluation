#include "stdafx.h"
#include "math5a_triangulation.h"

Math5A_Triangulation::Math5A_Triangulation(QWidget *parent)
	: QMainWindow(parent), ui(Ui::triangulationForm())
{
	ui.setupUi(this);
	glScene = new GLWidget(this);
	ui.centralLayout->addWidget(glScene);
	//connect(ui.config1, SIGNAL(clicked()), this, SLOT(newForm()));
}

Math5A_Triangulation::~Math5A_Triangulation()
{

}
