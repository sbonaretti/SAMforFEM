/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <VnlReaderMatrix.h>

#include <iostream>
#include <vector>

#include <QFile>
#include <QString>
#include <QTextStream>


namespace vnl{
	
	// constructor
	VnlReaderMatrix::VnlReaderMatrix(){
	}

	// destructor
	VnlReaderMatrix::~VnlReaderMatrix(){
	}
	
	// overwritten virtual function
	void VnlReaderMatrix::Update(){

		QFile file(_fileName);
		file.open(QIODevice::ReadOnly | QIODevice::Text);
		QTextStream readFile(&file);

		std::vector<double> matrixTemp;
		
		while( !readFile.atEnd() ){
			double temp; 
			readFile >> temp; 
			matrixTemp.push_back(temp);	
		}

		int r = matrixTemp.at(0);
		int c = matrixTemp.at(1);

		_matrix.set_size(r,c);

		int count = 2;
		for (int i=0; i<r; i++){
			for (int j=0; j<c; j++){
				_matrix(i,j) = matrixTemp.at(count);
				count ++;
			}
		}

		file.close();
	}

}