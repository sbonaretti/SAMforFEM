/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <PointReaderXyzId.h>

#include <QProcess.h>
#include <QFile.h>
#include <QTextStream.h>

namespace points{

	// constructor
	PointReaderXyzId::PointReaderXyzId(){
	}
	
	// destructor
	PointReaderXyzId::~PointReaderXyzId(){
	}

	// overwritten virtual function
	void PointReaderXyzId::Update(){
		
		QFile inFile(_fileName);
		inFile.open(QIODevice::ReadOnly | QIODevice::Text);
		QTextStream file(&inFile);
		
		_points = vtkPoints::New();
		_idVector = vtkDoubleArray::New();
		while (!file.atEnd()){
			double coord[3];
			file >> coord[0]; file >> coord[1]; file >> coord[2];
			_points->InsertNextPoint(coord);
			int id; file >> id;
			_idVector->InsertNextValue(id);
		}

		inFile.close();
	
	}

}