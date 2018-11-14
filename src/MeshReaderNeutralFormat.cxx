/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <MeshReaderNeutralFormat.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

#include <QFile>
#include <QTextStream>
#include <QString>


namespace mesh{

	// constructor
	MeshReaderNeutralFormat::MeshReaderNeutralFormat(){
	}

	// destructor
	MeshReaderNeutralFormat::~MeshReaderNeutralFormat(){
	}
	
	// overwritten virtual function
	
	void MeshReaderNeutralFormat::Update(){

	
		QFile file(_fileName); // open file to read
		if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) 
			cerr << "Error: temp.inp file could not be opened" << endl;
		QTextStream in(&file);

		vtkPoints* points = vtkPoints::New();
		vtkCellArray* polys = vtkCellArray::New();
			

		// nodes
		QString line;
		in >> line; // number of nodes
		int nOfNodes = atof(line);
		
		for (int i=0; i<nOfNodes; i++){
			double coord[3];
			in >> coord[0]; in >> coord[1]; in >> coord[2];
			points->InsertNextPoint(coord);
		}	
		
		_mesh->SetPoints(points);

		
		// elements
		in >> line; // number of elements
		int nOfElements = atof(line);

		for (int i=0; i<nOfElements; i++) {

			double element[10];
			in >> line; // number 1 in the row - unuseful
			in >> element[0]; 
			in >> element[1]; 
			in >> element[2]; 
			in >> element[3]; 
			in >> element[4]; 
			in >> element[5]; 
			in >> element[6];
			in >> element[7]; 
			in >> element[8]; 
			in >> element[9]; 

			polys->InsertNextCell(10);
			polys->InsertCellPoint(element[0]-1);
			polys->InsertCellPoint(element[1]-1); 
			polys->InsertCellPoint(element[2]-1);
			polys->InsertCellPoint(element[3]-1);
			polys->InsertCellPoint(element[4]-1);
			polys->InsertCellPoint(element[5]-1);
			polys->InsertCellPoint(element[6]-1);
			polys->InsertCellPoint(element[7]-1);
			polys->InsertCellPoint(element[8]-1);
			polys->InsertCellPoint(element[9]-1);
		}

		_mesh->SetPolys(polys);
		
		file.close(); // close file

		// cleaning up
		points->Delete();
		polys->Delete();

		std::cout << "nodes: " << _mesh->GetNumberOfPoints() << std::endl;
		std::cout << "elements: " << _mesh->GetNumberOfCells() << std::endl;

	}
}