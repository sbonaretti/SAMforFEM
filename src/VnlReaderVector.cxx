/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <VnlReaderVector.h>

#include <vector>

#include <QFile>
#include <QString>
#include <QTextStream>


namespace vnl{
	
	// constructor
	VnlReaderVector::VnlReaderVector(){
	}

	// destructor
	VnlReaderVector::~VnlReaderVector(){
	}
	
	// overwritten virtual function
	void VnlReaderVector::Update(){

		QFile file(_fileName);
		file.open(QIODevice::ReadOnly | QIODevice::Text);
		QTextStream readFile(&file);

		std::vector<double> vectorTemp;
		
		while( !readFile.atEnd() ){
			double temp; 
			readFile >> temp; vectorTemp.push_back(temp);	
		}

		_vector.set_size(vectorTemp.size());
		for (int i=0; i<vectorTemp.size(); i++){
			_vector(i) = vectorTemp.at(i);
		}		

		file.close();
	}
}