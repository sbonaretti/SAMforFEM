/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <PointWriterXyzVector.h>

namespace points{

	// constructor
	PointWriterXyzVector::PointWriterXyzVector(){
	}
	
	// destructor
	PointWriterXyzVector::~PointWriterXyzVector(){
	}

	// overwritten virtual function
	void PointWriterXyzVector::Update(){
		
		QFile outFile(_fileName.toAscii().data());
		outFile.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream writeFile(&outFile);
		for (int a=0; a<_points->GetNumberOfPoints(); a++){
			double coord[3];
			_points->GetPoint(a, coord);
			if (a == _points->GetNumberOfPoints()-1){
				writeFile << coord[0]; writeFile << endl;
				writeFile << coord[1]; writeFile << endl;
				writeFile << coord[2]; writeFile;
			}
			else {
				writeFile << coord[0]; writeFile << endl;
				writeFile << coord[1]; writeFile << endl;
				writeFile << coord[2]; writeFile << endl;
			}

		}
		outFile.close();
	
	}

}