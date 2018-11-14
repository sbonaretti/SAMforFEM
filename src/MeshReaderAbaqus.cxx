/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <MeshReaderAbaqus.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkIntArray.h>

#include <vtkGenericCell.h>
#include <vtkGeometryFilter.h>

#include <QFile>
#include <QTextStream>
#include <QString>


namespace mesh{

	// constructor
	MeshReaderAbaqus::MeshReaderAbaqus(){

		_elementType = 10;  //change also _nodeID in AbaqusFileCreation.h
		_abaqusElementType = ("C3D10");

	}

	// destructor
	MeshReaderAbaqus::~MeshReaderAbaqus(){
	}
	
	// overwritten virtual function
	
	void MeshReaderAbaqus::Update(){

		// read file
		QFile file(_fileName);
		if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) 
			cerr << "Error: mesh file could not be opened" << endl;
		QTextStream in(&file);

		QString line;
		
		vtkPoints* points = vtkPoints::New();
		vtkIntArray* nodeID = vtkIntArray::New();
		vtkCellArray* polys = vtkCellArray::New();
		int flag = 0;
		double coord[3];
		int nOfCoord = 3;
		int nOfNodes = 0;
		int nOfElements= 0;

		
		// number of nodes and elements
		while( !in.atEnd() ){
			
			in >> line; // evaluates the first element of the line, not the entire line (i.e. for the second node or element: line == 2,)
			
			if( line.contains("*Node") ){
				flag = 1;
				in >> line; // takes "1,"
			}
			
			if( line.contains("*Element") ){
				flag = 2;
				line = in.readLine(); // takes the whole *Element line
				in >> line; // takes "1,"
			}

			
			// read nodes, find number of nodes and put them in array
			if (flag == 1){
				line.remove(QChar (','));
				int temp = line.toInt();
				
				if (nOfNodes == 0){
					
					nOfNodes = 1;
					
					for (int a=0; a<nOfCoord; a++){
						char comma;
						if (a == nOfCoord-1){ //end of row
							in >> coord[a];
						}
						else {
							in >> coord[a];
							in >> comma;
						}
					}
					points->InsertNextPoint(coord);
				}
							
				else if (nOfNodes != 0 && temp == nOfNodes+1){
					nOfNodes = temp;

					for (int a=0; a<nOfCoord; a++){
						char comma;
						if (a == nOfCoord-1){ //end of row
							in >> coord[a];
						}
						else {
							in >> coord[a];
							in >> comma;
						}
					}
					points->InsertNextPoint(coord);
				
				}
			}
		
		
		// read elements and put them in the array
			if (flag == 2){
				
				line.remove(QChar (','));
				int temp = line.toInt();

				if (nOfElements == 0){
					
					nOfElements = 1;					
					
					int temp[10];
					for (int a=0; a<_elementType; a++){
						char comma;
						if (a == _elementType-1) //end of row
							in >> temp[a];
						else {
							in >> temp[a];
							in >> comma;
						}
					}
					polys->InsertNextCell(10);
					polys->InsertCellPoint(temp[0]-1);
					polys->InsertCellPoint(temp[1]-1); // to have normals pointing outside
					polys->InsertCellPoint(temp[2]-1);
					polys->InsertCellPoint(temp[3]-1);
					polys->InsertCellPoint(temp[4]-1);
					polys->InsertCellPoint(temp[5]-1); //7
					polys->InsertCellPoint(temp[6]-1); //5
					polys->InsertCellPoint(temp[7]-1); //6
					polys->InsertCellPoint(temp[8]-1);
					polys->InsertCellPoint(temp[9]-1);
				}

				else if (nOfElements != 0 && temp == nOfElements+1){
					
					nOfElements = temp;

					int temp[10];
					for (int a=0; a<_elementType; a++){
						char comma;
						if (a == _elementType-1) //end of row
							in >> temp[a];
						else {
							in >> temp[a];
							in >> comma;
						}
					}
					polys->InsertNextCell(10);
					polys->InsertCellPoint(temp[0]-1);
					polys->InsertCellPoint(temp[1]-1); // to have normals pointing outside
					polys->InsertCellPoint(temp[2]-1);
					polys->InsertCellPoint(temp[3]-1);
					polys->InsertCellPoint(temp[4]-1);
					polys->InsertCellPoint(temp[5]-1); //7
					polys->InsertCellPoint(temp[6]-1); //5
					polys->InsertCellPoint(temp[7]-1); //6
					polys->InsertCellPoint(temp[8]-1);
					polys->InsertCellPoint(temp[9]-1);
				}
			}
		}
		
		file.close();

		_mesh->SetPoints(points);
		_mesh->SetPolys(polys);

		// cleaning up
		points->Delete();
		nodeID->Delete();
		polys->Delete();

		std::cout << "nodes: " << _mesh->GetNumberOfPoints() << std::endl;
		std::cout << "elements: " << _mesh->GetNumberOfCells() << std::endl;

	}

}