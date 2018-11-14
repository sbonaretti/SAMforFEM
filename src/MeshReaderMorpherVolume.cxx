/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <MeshReaderMorpherVolume.h>
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
	MeshReaderMorpherVolume::MeshReaderMorpherVolume(){
	}

	// destructor
	MeshReaderMorpherVolume::~MeshReaderMorpherVolume(){
	}
	
	// overwritten virtual function
	
	void MeshReaderMorpherVolume::Update(){

		// read file
		QFile file(_fileName);
		if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) 
			cerr << "Error: mesh file could not be opened" << endl;
		QTextStream in(&file);
		
		int flag = 0;
		vtkPoints* points = vtkPoints::New();
		vtkIntArray* nodeID = vtkIntArray::New();
		int nOfNodes = 0;
		int nodeIndex = 0;
		vtkCellArray* polys = vtkCellArray::New();
		int nOfElements = 0;
		int elementIndex = 0;
		
		while( !in.atEnd() ){
			
			QString line;
			in >> line; // evaluates the first element of each line (first element, not the entire line (i.e. for the second node or element: line == 2,))
			
			if( line.contains("NBLOCK") ){
				flag = 1;
				in >> nOfNodes; // takes the following line
				//std::cout << nOfNodes << std::endl;
				in >> line; 
				in >> line; 
				in >> line; 

			}

			if( line.contains("EBLOCK") ){
				flag = 2;
				in >> line;
				in >> nOfElements;
				//std::cout << nOfElements << std::endl;
				line = in.readLine(); // takes the whole *Element line
				in >> line; // takes the following line
			}
		
			// read nodes and put them in _mesh
				if (flag == 1){
					if (nodeIndex < nOfNodes){
						if (nodeIndex == 0){
							int number;
							in >> number;	//std::cout << number << std::endl;
							in >> number;	//std::cout << number << std::endl;
							in >> number;	//std::cout << number << std::endl;
							double coord[3];
							in >> coord[0]; in >> coord[1]; in >> coord[2];	//std::cout << coord[0] << ' '  << coord[1] << ' ' << coord[2] << std::endl;
							points->InsertNextPoint(coord);
							nodeIndex ++;
						}
						else {
							int number;
							in >> number;	//std::cout << number << std::endl;
							in >> number;	//std::cout << number << std::endl;
							double coord[3];
							in >> coord[0]; in >> coord[1]; in >> coord[2];	//std::cout << coord[0] << ' '  << coord[1] << ' ' << coord[2] << std::endl;
							points->InsertNextPoint(coord);
							nodeIndex ++;
						}
					}
				}
			
			// read elements and put them in _mesh
			if (flag == 2){
				if (elementIndex < nOfElements){
					int temp[10]; 
					if (elementIndex == 0){
						int number;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
											
						in >> temp[0]; //std::cout << temp[0] << std::endl;
						in >> temp[1]; //std::cout << temp[1] << std::endl;
						in >> temp[2]; //std::cout << temp[2] << std::endl;
						in >> temp[3]; //std::cout << temp[3] << std::endl;
						in >> temp[4]; //std::cout << temp[4] << std::endl;
						in >> temp[5]; //std::cout << temp[5] << std::endl;
						in >> temp[6]; //std::cout << temp[6] << std::endl;
						in >> temp[7]; //std::cout << temp[7] << std::endl;
						in >> temp[8]; //std::cout << temp[8] << std::endl;
						in >> temp[9]; //std::cout << temp[9] << std::endl;
						
						polys->InsertNextCell(10);
						polys->InsertCellPoint(temp[0]-1);
						polys->InsertCellPoint(temp[1]-1); // to have normals pointing outside
						polys->InsertCellPoint(temp[2]-1);
						polys->InsertCellPoint(temp[3]-1);
						polys->InsertCellPoint(temp[4]-1);
						polys->InsertCellPoint(temp[5]-1);//7
						polys->InsertCellPoint(temp[6]-1);//5
						polys->InsertCellPoint(temp[7]-1);//6
						polys->InsertCellPoint(temp[8]-1);
						polys->InsertCellPoint(temp[9]-1);

						elementIndex ++;
					}
		
					else {
						int number;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
						in >> number; //std::cout << number << std::endl;
												
						in >> temp[0]; //std::cout << temp[0] << std::endl;
						in >> temp[1]; //std::cout << temp[1] << std::endl;
						in >> temp[2]; //std::cout << temp[2] << std::endl;
						in >> temp[3]; //std::cout << temp[3] << std::endl;
						in >> temp[4]; //std::cout << temp[4] << std::endl;
						in >> temp[5]; //std::cout << temp[5] << std::endl;
						in >> temp[6]; //std::cout << temp[6] << std::endl;
						in >> temp[7]; //std::cout << temp[7] << std::endl;
						in >> temp[8]; //std::cout << temp[8] << std::endl;
						in >> temp[9]; //std::cout << temp[9] << std::endl;
						//std::cout << temp[0] << ' ' << temp[1] << ' ' << temp[2] << ' ' << temp[3] << ' ' << temp[4] << ' ' << temp[5] << ' ' << temp[6] << ' ' << temp[7] << ' ' << temp[8] << ' ' << temp[9] << std::endl;
						
						polys->InsertNextCell(10);
						polys->InsertCellPoint(temp[0]-1);
						polys->InsertCellPoint(temp[1]-1); 
						polys->InsertCellPoint(temp[2]-1);
						polys->InsertCellPoint(temp[3]-1);
						polys->InsertCellPoint(temp[4]-1);
						polys->InsertCellPoint(temp[5]-1);
						polys->InsertCellPoint(temp[6]-1);
						polys->InsertCellPoint(temp[7]-1);
						polys->InsertCellPoint(temp[8]-1);
						polys->InsertCellPoint(temp[9]-1);

						elementIndex ++;
					}

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