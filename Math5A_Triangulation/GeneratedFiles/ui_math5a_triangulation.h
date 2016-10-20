/********************************************************************************
** Form generated from reading UI file 'math5a_triangulation.ui'
**
** Created by: Qt User Interface Compiler version 5.7.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MATH5A_TRIANGULATION_H
#define UI_MATH5A_TRIANGULATION_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_triangulationForm
{
public:
    QWidget *centralWidget;
    QHBoxLayout *horizontalLayout_2;
    QHBoxLayout *centralLayout;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;
    QDockWidget *dockWidget;
    QWidget *dockWidgetContents_2;
    QVBoxLayout *verticalLayout_2;
    QVBoxLayout *verticalLayout;
    QPushButton *config3;
    QPushButton *config2;
    QPushButton *config1;
    QSpacerItem *horizontalSpacer;

    void setupUi(QMainWindow *triangulationForm)
    {
        if (triangulationForm->objectName().isEmpty())
            triangulationForm->setObjectName(QStringLiteral("triangulationForm"));
        triangulationForm->resize(1510, 966);
        centralWidget = new QWidget(triangulationForm);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        horizontalLayout_2 = new QHBoxLayout(centralWidget);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        centralLayout = new QHBoxLayout();
        centralLayout->setSpacing(6);
        centralLayout->setObjectName(QStringLiteral("centralLayout"));

        horizontalLayout_2->addLayout(centralLayout);

        triangulationForm->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(triangulationForm);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1510, 31));
        triangulationForm->setMenuBar(menuBar);
        mainToolBar = new QToolBar(triangulationForm);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        triangulationForm->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(triangulationForm);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        triangulationForm->setStatusBar(statusBar);
        dockWidget = new QDockWidget(triangulationForm);
        dockWidget->setObjectName(QStringLiteral("dockWidget"));
        dockWidgetContents_2 = new QWidget();
        dockWidgetContents_2->setObjectName(QStringLiteral("dockWidgetContents_2"));
        verticalLayout_2 = new QVBoxLayout(dockWidgetContents_2);
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setContentsMargins(11, 11, 11, 11);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(6);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        config3 = new QPushButton(dockWidgetContents_2);
        config3->setObjectName(QStringLiteral("config3"));

        verticalLayout->addWidget(config3);

        config2 = new QPushButton(dockWidgetContents_2);
        config2->setObjectName(QStringLiteral("config2"));

        verticalLayout->addWidget(config2);

        config1 = new QPushButton(dockWidgetContents_2);
        config1->setObjectName(QStringLiteral("config1"));

        verticalLayout->addWidget(config1);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        verticalLayout->addItem(horizontalSpacer);


        verticalLayout_2->addLayout(verticalLayout);

        dockWidget->setWidget(dockWidgetContents_2);
        triangulationForm->addDockWidget(static_cast<Qt::DockWidgetArea>(1), dockWidget);

        retranslateUi(triangulationForm);

        QMetaObject::connectSlotsByName(triangulationForm);
    } // setupUi

    void retranslateUi(QMainWindow *triangulationForm)
    {
        triangulationForm->setWindowTitle(QApplication::translate("triangulationForm", "Math5A_Triangulation", 0));
        config3->setText(QApplication::translate("triangulationForm", "Config1222", 0));
        config2->setText(QApplication::translate("triangulationForm", "Config2", 0));
        config1->setText(QApplication::translate("triangulationForm", "Config3", 0));
    } // retranslateUi

};

namespace Ui {
    class triangulationForm: public Ui_triangulationForm {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MATH5A_TRIANGULATION_H
