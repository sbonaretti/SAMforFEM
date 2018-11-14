/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <PointReaderXyz.h>

#include <QProcess.h>
#include <QFile.h>
#include <QTextStream.h>

namespace points{

	// constructor
	PointReaderXyz::PointReaderXyz(){
	}
	
	// destructor
	PointReaderXyz::~PointReaderXyz(){
	}

	// overwritten virtual function
	void PointReaderXyz::Update(){
		
		QFile inFile(_fileName);
		inFile.open(QIODevice::ReadOnly | QIODevice::Text);
		QTextStream file(&inFile);
		
		int index = 0;
		while (!file.atEnd()){
			double coord[3];
			file >> coord[0]; file >> coord[1]; file >> coord[2];
			_points->InsertNextPoint(coord); 
		}
		inFile.close();
	
	}

}