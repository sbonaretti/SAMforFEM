/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <MeshWriterAbaqus.h>

#include <vtkGenericCell.h>
#include <vtkIdList.h>

#include <QFile>
#include <QString>
#include <QTextStream>


namespace mesh{
	
	// constructor
	MeshWriterAbaqus::MeshWriterAbaqus(){

		_elementType = 10;  //change also _nodeID in AbaqusFileCreation.h
		_abaqusElementType = ("C3D10");
		_poissonRatio = 0.3;
		_loadPoint1 = 0;
		_loadPoint2 = 0;
	}

	// destructor
	MeshWriterAbaqus::~MeshWriterAbaqus(){
	}
	
	// overwritten virtual function
	void MeshWriterAbaqus::Update(){

		QFile outFile(_fileName);
		outFile.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream abaqusFile(&outFile);
			
		// header
		abaqusFile << "*Heading"  << endl;	
		abaqusFile << "**Abaqus input file" << endl;
		// defining part
		if (_youngModulusElemFlag == true || _youngModulusNodeFlag == true){
			abaqusFile << "**PARTS" << endl;
			abaqusFile << "*Part, name=PART-1" << endl;
		}

		
		/************ MESH ************/
		if (_meshFlag == true){ 
			
			// nodes
			abaqusFile << "**" << endl;
			abaqusFile << "*Node" << endl;
			for (int i=0; i<_mesh->GetNumberOfPoints(); i++){
				double coord[3];
				_mesh->GetPoint(i, coord);
				abaqusFile << i+1 << ", " << coord[0] <<  ", " << coord[1] << ", " << coord[2] << endl;
			}

			if (_youngModulusElemFlag == false &&  _youngModulusNodeFlag == false){
				// elements
				vtkGenericCell* cell = vtkGenericCell::New();
				vtkIdList* idList;
				abaqusFile << "*Element, type=" << _abaqusElementType << endl; // check the kind of element in abaqus
				for (int i=0; i<_mesh->GetNumberOfCells(); i++){
					_mesh->GetCell(i, cell);
					idList = cell->GetPointIds();
					abaqusFile << i+1;
					for (int a=0; a<_elementType; a++){
						abaqusFile << ", " << idList->GetId(a)+1;
						if (a == _elementType-1)
							abaqusFile << endl;
					}
				}
			}
		}

		
		/************ YOUNG'S MODULUS FOR ELEMENTS ************/
		if (_youngModulusElemFlag == true){

			// elements
			vtkGenericCell* cell = vtkGenericCell::New();
			vtkIdList* idList;
			for (int i=0; i<_mesh->GetNumberOfCells(); i++){
				abaqusFile << "*Element, type=" << _abaqusElementType << ", elset=P_SET_"<< i+1 << endl; // check the kind of element in abaqus
				_mesh->GetCell(i, cell);
				idList = cell->GetPointIds();
				abaqusFile << i+1;
				for (int a=0; a<_elementType; a++){
					abaqusFile << ", " << idList->GetId(a)+1;
					if (a == _elementType-1)
						abaqusFile << endl;
				}
			}

			// part
			abaqusFile << "**" << endl;
			for (int i=0; i<_mesh->GetNumberOfCells(); i++){
				abaqusFile << "*Solid Section, elset=P_SET_" << i+1 << ", material=MATRL_" << i+1 << endl;
				abaqusFile << "1.," << endl;
			}
			abaqusFile << "*End Part" << endl;
			
			// assembly
			abaqusFile << "**" << endl;
			abaqusFile << "** ASSEMBLY" << endl;
			abaqusFile << "*Assembly, name=Assembly" << endl;
			abaqusFile << "*Instance, name=PART-1-1, part=PART-1" << endl;
			abaqusFile << "*End Instance" << endl;
			
			
			/************ LOAD POINT APPLICATION ************/
			if (_loadFlag == true){

				// load application point 1
				double pt[3]; // input point
				_loadPoints->GetPoint(0,pt);
				std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
				findClosestNode (pt, _loadPoint1);

				// load application point 2
				_loadPoints->GetPoint(1,pt);
				std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
				findClosestNode (pt,_loadPoint2);
				
				abaqusFile << "**" << endl;
				abaqusFile << "*Nset, nset=_PickedSet1, internal, instance=PART-1-1" << endl;
				abaqusFile << _loadPoint1 << "," << endl; // node ID
			}
	
			/************ BOUNDARY CONDITION APPLICATION POINT ************/
			if (_bcFlag == true){
	
				
				// bc application point 1
				double pt[3]; // input point
				int bcApplicationPoint1;
				_bcPoints->GetPoint(0,pt);
				std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
				findClosestNode (pt, bcApplicationPoint1);
				
				// bc application point 2
				int bcApplicationPoint2;
				_bcPoints->GetPoint(1,pt);
				std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
				findClosestNode (pt,bcApplicationPoint2);

				// bc application point 3
				int bcApplicationPoint3;
				_bcPoints->GetPoint(2,pt);
				std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
				findClosestNode (pt,bcApplicationPoint3);
								
				// bc point application
				abaqusFile << "*Nset, nset=_PickedSet5, internal, instance=PART-1-1" << endl;
				abaqusFile << bcApplicationPoint1 << "," << endl; // node ID	
				abaqusFile << "*Nset, nset=_PickedSet6, internal, instance=PART-1-1" << endl;
				abaqusFile << bcApplicationPoint2 << "," << endl; // node ID	
				abaqusFile << "*Nset, nset=_PickedSet7, internal, instance=PART-1-1" << endl;
				abaqusFile << bcApplicationPoint3 << "," << endl; // node ID	
			}

			abaqusFile << "*End Assembly" << endl;

			// materials
			abaqusFile << "**" << endl;
			abaqusFile << "** MATERIALS" << endl;
			for (int i=0; i<_mesh->GetNumberOfCells(); i++){
				abaqusFile << "**" << endl;
				abaqusFile << "*Material, name=MATRL_" << i+1 << endl;
				abaqusFile << "*Elastic, Type=ISO" << endl; 
				abaqusFile << _youngModulus->GetValue(i) * 1000 << ", " << _poissonRatio << endl; // *1000 = from GPa to MPa
			}
		}


		// step
		abaqusFile << "** ----------------------------------------------------------------" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "** STEP" << endl;
		abaqusFile << "*Step, name=Step-1" << endl;
		abaqusFile << "*Static" << endl;
		abaqusFile << "1., 1., 1e-05, 1." << endl;
		
		/************ LOAD POINT APPLICATION ************/
		if (_loadFlag == true){
				
			// getting the coordinates from the index
			double loadPoint1[3];
			_mesh->GetPoint(_loadPoint1-1, loadPoint1);
			double loadPoint2[3];
			_mesh->GetPoint(_loadPoint2-1, loadPoint2);
			
			// decomposing force magnitude along the 3 force directions
			int nOfCoord = 3;
			double diff[3];
			for (int a=0; a<nOfCoord; a++)
				diff[a] = loadPoint1[a] - loadPoint2[a];
			double norm = 0.0;
			for (int a=0; a<nOfCoord; a++)
				norm += diff[a]*diff[a];
			norm = sqrt (norm);
			for (int a=0; a<nOfCoord; a++)
				diff[a] /= norm;
			
			double loadComponents[3];
			double force = _loadMagnitude->GetValue(0);
			for (int a=0; a<nOfCoord; a++)
				loadComponents[a]=force*diff[a];

			abaqusFile << "**" << endl;
			abaqusFile << "** LOAD - Concentrated force" << endl;
			abaqusFile << "*Cload, op=NEW" << endl;
			abaqusFile << "_PickedSet1, 1, " << loadComponents[0] << endl;
			abaqusFile << "_PickedSet1, 2, " << loadComponents[1] << endl;
			abaqusFile << "_PickedSet1, 3, " << loadComponents[2] << endl;
		}

		
		/************ BOUNDARY CONDITION APPLICATION POINT ************/
		if (_bcFlag == true){
			
			abaqusFile << "**" << endl;
			abaqusFile << "** BOUNDARY CONDITIONS" << endl;
			abaqusFile << "**" << endl;
			abaqusFile << "** Name: BC-1 Type: Displacement/Rotation" << endl;
			abaqusFile << "*Boundary" << endl;
			abaqusFile << "_PickedSet5, 1, 1" << endl;
			abaqusFile << "_PickedSet5, 3, 3" << endl;
			abaqusFile << "** Name: BC-2 Type: Displacement/Rotation" << endl;
			abaqusFile << "*Boundary" << endl;
			abaqusFile << "_PickedSet6, 1, 1" << endl;
			abaqusFile << "_PickedSet6, 2, 2" << endl;
			abaqusFile << "_PickedSet6, 3, 3" << endl;
			abaqusFile << "_PickedSet6, 4, 4" << endl;
			abaqusFile << "_PickedSet6, 5, 5" << endl;
			abaqusFile << "_PickedSet6, 6, 6" << endl;
			abaqusFile << "** Name: BC-3 Type: Displacement/Rotation" << endl;
			abaqusFile << "*Boundary" << endl;
			abaqusFile << "_PickedSet7, 2, 2" << endl;
		}
	
				
		
		/************ YOUNG'S MODULUS FOR NODES ************/
		if (_youngModulusNodeFlag == true){

			// elements
			vtkGenericCell* cell = vtkGenericCell::New();
			vtkIdList* idList;
			abaqusFile << "*Element, type=" << _abaqusElementType << endl; // check the kind of element in abaqus
			for (int i=0; i<_mesh->GetNumberOfCells(); i++){
				_mesh->GetCell(i, cell);
				idList = cell->GetPointIds();
				abaqusFile << i+1;
				for (int a=0; a<_elementType; a++){
					abaqusFile << ", " << idList->GetId(a)+1;
					if (a == _elementType-1)
						abaqusFile << endl;
				}
			}

			// part
			abaqusFile << "**" << endl;
			abaqusFile << "*Elset, elset=_PickedSet4, internal" << endl;
			for (int i=0; i<_mesh->GetNumberOfCells(); i++){
				if (i == _mesh->GetNumberOfCells()-1)
					abaqusFile << _mesh->GetNumberOfCells() << endl;
				else 
					abaqusFile << i+1 << ", ";
				if (i%10 == 0)
					abaqusFile << endl;

			}
			abaqusFile << "** Section: Section-1" << endl;
			abaqusFile << "*Solid Section, elset=_PickedSet4, material=Material-1" << endl;
			abaqusFile << "1.," << endl;
			abaqusFile << "*End Part" << endl;
		
			// assembly
			abaqusFile << "**" << endl;
			abaqusFile << "** ASSEMBLY" << endl;
			abaqusFile << "*Assembly, name=Assembly" << endl;
			abaqusFile << "*Instance, name=PART-1-1, part=PART-1" << endl;
			abaqusFile << "*End Instance" << endl;
			abaqusFile << "**" << endl;

			abaqusFile << "*End Assembly" << endl;

			// materials
			abaqusFile << "**" << endl;
			abaqusFile << "** MATERIALS" << endl;
			abaqusFile << "*Material, name=Material-1" << endl;
			abaqusFile << "*Elastic" << endl;
			double range[2];
			_youngModulus->GetRange(range);
			abaqusFile << range[0] << ", " << _poissonRatio << ", " <<  range[0] << endl;
			abaqusFile << range[1] << ", " << _poissonRatio << ", " <<  range[1] << endl;
			abaqusFile << "** ----------------------------------------------------------------" << endl;
			abaqusFile << "*INITIAL CONDITION, TYPE=TEMPERATURE" << endl;
			for (int i=0; i<_youngModulus->GetNumberOfTuples(); i++)
				abaqusFile << "PART-1-1." << i+1 << ", " << _youngModulus->GetValue(i) << endl;
		}

		// file end
		abaqusFile << "**" << endl;
		abaqusFile << "** OUTPUT REQUESTS" << endl;
		abaqusFile << "**" << endl; 
		abaqusFile << "*Restart, write, frequency=0" << endl;
		abaqusFile << "**" << endl; 
		abaqusFile << "** FIELD OUTPUT: F-Output-1" << endl;
		abaqusFile << "**" << endl; 
		abaqusFile << "*Output, field, variable=PRESELECT" << endl;
		abaqusFile << "**" << endl; 
		abaqusFile << "** HISTORY OUTPUT: H-Output-1" << endl;
		abaqusFile << "**" << endl; 
		abaqusFile << "*Output, history, variable=PRESELECT" << endl;
		abaqusFile << "*End Step" << endl;

		// close file
		outFile.close();

		}

		
		
		void MeshWriterAbaqus::findClosestNode (double point[3], int &minIndex){
	
		// find the absolute distances between the input point and all the mesh nodes
		double min = 0.0; //int minIndex;
		int nOfCoord = 3;
		
		for (int i=0; i<_mesh->GetNumberOfPoints(); i++){
			double coord[3];
			_mesh->GetPoint(i, coord);
			
			double distance = 0.0;
			for (int a=0; a<nOfCoord; a++)
				distance += ((coord[a]-(point[a])) * (coord[a]-(point[a])));
				distance = sqrt(distance);
			if (i==0){
				min = distance;
				minIndex=i+1; // in abaqus node id starts from 1
			}
			else{
				if (distance < min){
					min = distance;
					minIndex=i+1; // in abaqus node id starts from 1
					}
				}
		}
	}

}