/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <VnlReaderEValues.h>

#include <vector>

#include <QFile>
#include <QString>
#include <QTextStream>


namespace vnl{
	
	// constructor
	VnlReaderEValues::VnlReaderEValues(){
	}

	// destructor
	VnlReaderEValues::~VnlReaderEValues(){
	}
	
	// overwritten virtual function
	void VnlReaderEValues::Update(){

		QFile file(_fileName);
		file.open(QIODevice::ReadOnly | QIODevice::Text);
		QTextStream readFile(&file);

		std::vector<double> eValuesTemp;
		std::vector<double> normalizedEValuesTemp;

		while( !readFile.atEnd() ){
			double temp; 
			readFile >> temp;
			readFile >> temp; eValuesTemp.push_back(temp);
			readFile >> temp; normalizedEValuesTemp.push_back(temp);		
		}

		_eValues.set_size(eValuesTemp.size());
		_normalizedEValues.set_size(normalizedEValuesTemp.size());
		for (int i=0; i<eValuesTemp.size(); i++){
			_eValues(i) = eValuesTemp.at(i);
			_normalizedEValues(i) = normalizedEValuesTemp.at(i);
		}		

		file.close();
	}
}