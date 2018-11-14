/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <MeshReaderMorpherSurface.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkIntArray.h>

#include <QFile>
#include <QTextStream>
#include <QString>


namespace mesh{

	// constructor
	MeshReaderMorpherSurface::MeshReaderMorpherSurface(){
	}

	// destructor
	MeshReaderMorpherSurface::~MeshReaderMorpherSurface(){
	}
	
	// overwritten virtual function
	
	void MeshReaderMorpherSurface::Update(){
			
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
		_mesh->Allocate();

		while( !in.atEnd() ){
			
			QString line;
			in >> line; // evaluates the first element of each line (first element, not the entire line (i.e. for the second node or element: line == 2,))
			
			if( line.contains("NBLOCK") ){
				flag = 1;
				in >> nOfNodes; // takes the following line
				// std::cout << nOfNodes << std::endl;
				in >> line; 
				in >> line; 
				in >> line; 

			}

			if( line.contains("EBLOCK") ){
				flag = 2;
				in >> line;
				in >> nOfElements;
				// std::cout << nOfElements << std::endl;
				line = in.readLine(); // takes the whole *Element line
				in >> line; // takes the following line
			}
			
			// read nodes and put them in _mesh
			if (flag == 1){
				if (nodeIndex < nOfNodes){
					if (nodeIndex == 0){
						int number;
						in >> number;	//std::cout << number << std::endl;
						nodeID->InsertNextValue(number);
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						double coord[3];
						in >> coord[0]; in >> coord[1]; in >> coord[2];	//std::cout << coord[0] << ' '  << coord[1] << ' ' << coord[2] << std::endl;
						points->InsertNextPoint(coord);
						nodeIndex ++;
					}
					else {
						nodeID->InsertNextValue(line.toInt());
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
					int temp1[3]; int temp2[3];
					if (elementIndex == 0){
						int number;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
											
						in >> temp1[0]; //std::cout << temp1[0] << std::endl;
						in >> temp1[1]; //std::cout << temp1[1] << std::endl;
						in >> temp1[2]; //std::cout << temp1[2] << std::endl;
						
						polys->InsertNextCell(3);
						for (int i=0; i<nodeID->GetSize(); i++){
							if (nodeID->GetValue(i) == temp1[0])
								temp2[0] = i;
							if (nodeID->GetValue(i) == temp1[1])
								temp2[1] = i;
							if (nodeID->GetValue(i) == temp1[2])
								temp2[2] = i;
						}
						//std::cout << temp2[0] << ' '  << temp2[1] << ' ' << temp2[2] << std::endl;
						polys->InsertCellPoint(temp2[0]);
						polys->InsertCellPoint(temp2[2]); // to have normals pointing outside
						polys->InsertCellPoint(temp2[1]);
						elementIndex ++;
					}
		
					else {
						int number;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
						in >> number;	//std::cout << number << std::endl;
												
						in >> temp1[0]; 
						in >> temp1[1]; 
						in >> temp1[2];

						polys->InsertNextCell(3);
						for (int i=0; i<nodeID->GetSize(); i++){
							if (nodeID->GetValue(i) == temp1[0])
								temp2[0] = i;
							if (nodeID->GetValue(i) == temp1[1])
								temp2[1] = i;
							if (nodeID->GetValue(i) == temp1[2])
								temp2[2] = i;
						}
						//std::cout << nodeIDs[0] << ' '  << nodeIDs[1] << ' ' << nodeIDs[2] << std::endl;
						polys->InsertCellPoint(temp2[0]);
						polys->InsertCellPoint(temp2[2]); // to have normals pointing outside
						polys->InsertCellPoint(temp2[1]);
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