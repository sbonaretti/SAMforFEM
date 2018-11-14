/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <VnlWriterEValues.h>

#include <QFile>
#include <QString>
#include <QTextStream>


namespace vnl{
	
	// constructor
	VnlWriterEValues::VnlWriterEValues(){
	}

	// destructor
	VnlWriterEValues::~VnlWriterEValues(){
	}
	
	// overwritten virtual function
	void VnlWriterEValues::Update(){

		QFile file(_fileName);
		file.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream writeFile(&file);

		
		for (int i=0; i<_eValues.size(); i++){

			if (i == _eValues.size()-1){
				writeFile << i+1; writeFile << " "; writeFile << _eValues(i); writeFile << " "; writeFile << _normalizedEValues(i);
			}
			else {
				writeFile << i+1; writeFile << " "; writeFile << _eValues(i); writeFile << " "; writeFile << _normalizedEValues(i); writeFile << endl;
			}
		}

		file.close();
	}
}