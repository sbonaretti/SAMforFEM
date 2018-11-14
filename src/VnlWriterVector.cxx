/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <VnlWriterVector.h>

#include <QFile>
#include <QString>
#include <QTextStream>


namespace vnl{
	
	// constructor
	VnlWriterVector::VnlWriterVector(){
	}

	// destructor
	VnlWriterVector::~VnlWriterVector(){
	}
	
	// overwritten virtual function
	void VnlWriterVector::Update(){

		QFile file(_fileName);
		file.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream writeFile(&file);

		
		for (int i=0; i<_vector.size(); i++){
			if (i == _vector.size()-1)
				writeFile << _vector(i);
			else 
				writeFile << _vector(i) << ' ';
		}
		file.close();
	}
}