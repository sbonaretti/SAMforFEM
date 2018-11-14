/*
 * ISTB - University of Bern, Serena Bonaretti
 */

// QT includes
#include <qapplication.h>
#include "MainWindow.h"

int main( int argc, char** argv )
{
  
	// QT Stuff
  QApplication app( argc, argv );

  MainWindow mainwindow;
  app.setMainWidget(&mainwindow);
  mainwindow.show();

  return app.exec();
}