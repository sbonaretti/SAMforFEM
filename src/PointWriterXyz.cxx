/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <PointWriterXyz.h>

namespace points{

	// constructor
	PointWriterXyz::PointWriterXyz(){
	}
	
	// destructor
	PointWriterXyz::~PointWriterXyz(){
	}

	// overwritten virtual function
	void PointWriterXyz::Update(){
		
		QFile outFile(_fileName.toAscii().data());
		outFile.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream writeFile(&outFile);
		for (int a=0; a<_points->GetNumberOfPoints(); a++){
			double coord[3];
			_points->GetPoint(a, coord);
			if (a == _points->GetNumberOfPoints()-1){
				writeFile << coord[0]; writeFile << ' ';
				writeFile << coord[1]; writeFile << ' ';
				writeFile << coord[2]; writeFile;
			}
			else {
				writeFile << coord[0]; writeFile << ' ';
				writeFile << coord[1]; writeFile << ' ';
				writeFile << coord[2]; writeFile << endl;
			}

		}
		outFile.close();
	
	}

}