/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <MeshWriterAnsys.h>

#include <vtkGenericCell.h>
#include <vtkIdList.h>

#include <QFile>
#include <QString>
#include <QTextStream>


namespace mesh{
	
	// constructor
	MeshWriterAnsys::MeshWriterAnsys(){
	}

	// destructor
	MeshWriterAnsys::~MeshWriterAnsys(){
	}
	
	// overwritten virtual function
	void MeshWriterAnsys::Update(){
		
		QFile outFile(_fileName);
		outFile.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream ansysFile(&outFile);

		// write file header (as in the morpher volume output file)
		ansysFile << "/COM,ADOMSH - 12.0" << endl;
		ansysFile << "/PREP7" << endl;
		ansysFile << "/NOPR" << endl;
		ansysFile << "/TITLE," << endl;
		ansysFile << "*IF,_CDRDOFF,EQ,1,THEN     !if solid model was read in" << endl;
		ansysFile << "_CDRDOFF=             !reset flag, numoffs already performed" << endl;
		ansysFile << "*ELSE              !offset database for the following FE model" << endl;
		ansysFile << "NUMOFF,NODE,   " << _mesh->GetNumberOfPoints() << endl;
		ansysFile << "NUMOFF,ELEM,   " << _mesh->GetNumberOfCells() << endl;
		ansysFile << "NUMOFF,TYPE,       1" << endl;
		ansysFile << "*ENDIF" << endl;
		ansysFile << "*SET,_BUTTON ,  1.000000000000" << endl;
		ansysFile << "*SET,_RETURN ,  0.000000000000" << endl;
		ansysFile << "*SET,_STATUS ,  0.000000000000" << endl;
		ansysFile << "*SET,_UIQR   ,  1.000000000000" << endl;
		ansysFile << "DOF,DELETE" << endl;
		ansysFile << "ET,       1,200" << endl;
		ansysFile << "KEYOP,       1, 1,        9" << endl;

		// mesh
		if (_meshFlag == true){
			// for nodes: 3i8,6e20.13; for elements: 19i8 - used for Morpher.exe
		
			// write node header
			ansysFile << "NBLOCK,6,SOLID,     " << _mesh->GetNumberOfPoints() << ",     " << _mesh->GetNumberOfPoints() << endl;
			ansysFile << "(3i8,6e20.13)" << endl;
			
			// write nodes
			for (int i=0; i<_mesh->GetNumberOfPoints(); i++) {
				
				ansysFile << QString("%1").arg(i+1, 8, 10);
				ansysFile << QString("%1").arg(0, 8, 10);
				ansysFile << QString("%1").arg(0, 8, 10);
				
				double coord[3];
				_mesh->GetPoint(i, coord);
				ansysFile << QString("%1").arg(coord[0], 20, 'E',13);
				ansysFile << QString("%1").arg(coord[1], 20, 'E',13);
				ansysFile << QString("%1").arg(coord[2], 20, 'E',13) << endl;
			}

			// write element header
			ansysFile << "N,R5.3,LOC,       -1," << endl;
			ansysFile << "EBLOCK,19,SOLID,     " << _mesh->GetNumberOfCells() << ",     " << _mesh->GetNumberOfCells() << endl;
			ansysFile << "(19i8)" << endl;

			// write elements
			vtkGenericCell* cell = vtkGenericCell::New();
			vtkIdList* idList;
			for (int i=0; i<_mesh->GetNumberOfCells(); i++) {
				
				// get the cell
				_mesh->GetCell(i, cell);
				
				// get the 10 tetra vertices
				idList = cell->GetPointIds();

				ansysFile << QString("%1").arg(1, 8, 10);
				ansysFile << QString("%1").arg(1, 8, 10);
				ansysFile << QString("%1").arg(1, 8, 10);
				ansysFile << QString("%1").arg(1, 8, 10);
				ansysFile << QString("%1").arg(0, 8, 10);
				ansysFile << QString("%1").arg(0, 8, 10);
				ansysFile << QString("%1").arg(0, 8, 10);
				ansysFile << QString("%1").arg(0, 8, 10);
				ansysFile << QString("%1").arg(10, 8, 10);
				ansysFile << QString("%1").arg(0, 8, 10);
				ansysFile << QString("%1").arg(i+1, 8, 10);

				ansysFile << QString("%1").arg(idList->GetId(0)+1, 8, 10);
				ansysFile << QString("%1").arg(idList->GetId(1)+1, 8, 10);
				ansysFile << QString("%1").arg(idList->GetId(2)+1, 8, 10);
				ansysFile << QString("%1").arg(idList->GetId(3)+1, 8, 10);
				ansysFile << QString("%1").arg(idList->GetId(4)+1, 8, 10);
				ansysFile << QString("%1").arg(idList->GetId(5)+1, 8, 10);
				ansysFile << QString("%1").arg(idList->GetId(6)+1, 8, 10);
				ansysFile << QString("%1").arg(idList->GetId(7)+1, 8, 10) << endl;;
				ansysFile << QString("%1").arg(idList->GetId(8)+1, 8, 10);
				ansysFile << QString("%1").arg(idList->GetId(9)+1, 8, 10) << endl;

			}
		}
		
		// write file footer
		ansysFile << "      -1" << endl;

		// close file
		outFile.close();

		
	}

}