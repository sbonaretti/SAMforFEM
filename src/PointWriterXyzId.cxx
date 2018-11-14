/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <PointWriterXyzId.h>

namespace points{

	// constructor
	PointWriterXyzId::PointWriterXyzId(){
	}
	
	// destructor
	PointWriterXyzId::~PointWriterXyzId(){
	}

	// overwritten virtual function
	void PointWriterXyzId::Update(){
		
		QFile outFile(_fileName.toAscii().data());
		outFile.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream writeFile(&outFile);
		for (int a=0; a<_points->GetNumberOfPoints(); a++){
			double coord[3];
			_points->GetPoint(a, coord);
			if (a == _points->GetNumberOfPoints()-1){
				writeFile << coord[0]; writeFile << ' ';
				writeFile << coord[1]; writeFile << ' ';
				writeFile << coord[2]; writeFile << ' ';
				writeFile << _idVector->GetValue(a);
			}
			else {
				writeFile << coord[0]; writeFile << ' ';
				writeFile << coord[1]; writeFile << ' ';
				writeFile << coord[2]; writeFile << ' ';
				writeFile << _idVector->GetValue(a); writeFile << endl;
			}

		}
		outFile.close();
	
	}

}

