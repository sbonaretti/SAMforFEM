/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <VnlWriterMatrix.h>

#include <QFile>
#include <QString>
#include <QTextStream>


namespace vnl{
	
	// constructor
	VnlWriterMatrix::VnlWriterMatrix(){
	}

	// destructor
	VnlWriterMatrix::~VnlWriterMatrix(){
	}
	
	// overwritten virtual function
	void VnlWriterMatrix::Update(){

		QFile file(_fileName);
		file.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream writeFile(&file);

		int r = _matrix.rows();
		int c = _matrix.cols();

		// writes number of rows and number of columns
		writeFile << r; writeFile << ' '; writeFile << c; writeFile << endl;


		for (int i=0; i<r; i++){
			for (int j=0; j<c; j++){
				
				if (i == r-1 && j == c-1){
					writeFile << _matrix(i,j);
				}
				else {
					writeFile << _matrix(i,j); writeFile << ' ';
				}
			}
		}
		file.close();
	}

	void VnlWriterMatrix::MatrixShapeUpdate(){

		QFile file(_fileName);
		file.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream writeFile(&file);

		int r = _matrix.rows();
		int c = _matrix.cols();

		// writes number of rows and number of columns
		writeFile << r; writeFile << ' '; writeFile << c; writeFile << endl;


		for (int i=0; i<r; i++){
			for (int j=0; j<c; j++){
				
				if (i == r-1 && j == c-1){
					writeFile << _matrix(i,j);
				}
				else if (i != r-1 && j == c-1){
					writeFile << _matrix(i,j) << endl;
				}
				else {
					writeFile << _matrix(i,j); writeFile << ' ';
				}
			}
		}
		file.close();
	
	
	}

}