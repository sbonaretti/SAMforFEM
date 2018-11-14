/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <FemAbaqusInpWriter.h>

#include <vtkGenericCell.h>
#include <vtkIdList.h>

#include <QFile>
#include <QTextStream>

namespace fem{

	// constructor
	FemAbaqusInpWriter::FemAbaqusInpWriter(){

		_elementType = 10;
		_abaqusElementType = ("C3D10");
		_poissonRatio = 0.3;
	}

	// destructor
	FemAbaqusInpWriter::~FemAbaqusInpWriter(){
	}

	
	// member functions
	void FemAbaqusInpWriter::WriteNodeWalkingInp(){


		QFile outFile(_fileName);
		outFile.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream abaqusFile(&outFile);
			
		// header
		abaqusFile << "*Heading" << endl;
		abaqusFile << "*Preprint, echo=NO, model=NO, history=NO, contact=NO" << endl;
		abaqusFile << "**" << endl;

		//part
		abaqusFile << "** PART" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "*Part, name=FEMUR" << endl;

		// nodes
		abaqusFile << "**" << endl;
		abaqusFile << "*Node" << endl;
		for (int i=0; i<_mesh->GetNumberOfPoints(); i++){
			double coord[3];
			_mesh->GetPoint(i, coord);
			abaqusFile << i+1 << ", " << coord[0] <<  ", " << coord[1] << ", " << coord[2] << endl;
		}
		
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
				if (a == _elementType-1){
					abaqusFile << endl;
				}
			}
		}
		
		// Nset
		abaqusFile << "*Nset, nset=FULLBONE, generate" << endl;
		abaqusFile << "1, " << _mesh->GetNumberOfPoints() << ", 1" << endl;
		// Elset
		abaqusFile << "*Elset, elset=FULLBONE_ELEM,generate " << endl;
		abaqusFile << "1, " << _mesh->GetNumberOfCells() << ", 1" << endl;
		
		// section
		abaqusFile << "**Section: FEMURSECTION" << endl;
		abaqusFile << "*Solid Section, elset=FULLBONE_ELEM, material=BONE" << endl; 
		abaqusFile << "," << endl;
		abaqusFile << "*End Part" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "**" << endl;
		
		// assembly
		abaqusFile << "*Assembly, name=ASSEMBLY" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "**" << endl;
		
		// instance
		abaqusFile << "*Instance, name=FEMURPART, part=FEMUR" << endl;
		abaqusFile << "*End Instance" << endl; 
		
		// Nset
		abaqusFile << "*Nset, nset=FULLBONE, instance=FEMURPART" << endl;
		abaqusFile << "1, " << _mesh->GetNumberOfPoints() << ", 1" << endl;
		// Elset
		abaqusFile << "*Elset, elset=FULLBONE_ELEM, instance=FEMURPART" << endl;
		abaqusFile << "1, " << _mesh->GetNumberOfCells() << ", 1" << endl;
		
		// Application points and coordinate system tranformation
		abaqusFile << "*Nset, nset=P0, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(0) << "," << endl;
		abaqusFile << "*Nset, nset=P1, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(1) << "," << endl; 
		abaqusFile << "*Nset, nset=P2, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(2) << "," << endl; 
		abaqusFile << "*Nset, nset=F1, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(3) << "," << endl;
		abaqusFile << "*Nset, nset=F3, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(4) << "," << endl;
		abaqusFile << "*Nset, nset=F5, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(5) << "," << endl;
		abaqusFile << "*Nset, nset=_T-LOCAL_BC, internal" << endl;
		abaqusFile << "P0," << endl;
		abaqusFile << "F3," << endl;
		abaqusFile << "F5," << endl;
		abaqusFile << "*Transform, nset=_T-LOCAL_BC" << endl; 
		double A[3];
		_boundaryConditions->GetPoint(5,A);
		double B[3];
		_boundaryConditions->GetPoint(6,B);
		abaqusFile << A[0] << "," << A[1] << "," << A[2] << "," << B[0] << "," << B[1] << "," << B[2] << endl;
		abaqusFile << "*Nset, nset=_T-LOCAL_FCS, internal" << endl; 
		abaqusFile << "P1," << endl;
		abaqusFile << "P2," << endl;
		abaqusFile << "*Transform, nset=_T-LOCAL_FCS" << endl;
		double C[3];
		_boundaryConditions->GetPoint(7,C);
		double D[3];
		_boundaryConditions->GetPoint(8,D);
		abaqusFile << C[0] << "," << C[1] << "," << C[2] << "," << D[0] << "," << D[1] << "," << D[2] << endl;
		abaqusFile << "*End Assembly" << endl;
		
		// material properties
		abaqusFile << "*Material, name=BONE" << endl;
		abaqusFile << "*Elastic " << endl;
		abaqusFile << "50, 0.3,50," << endl;
		abaqusFile << "30000, 0.3,30000," << endl;
		abaqusFile << "**PREDEFINED CONDITIONS" << endl;
		abaqusFile << "** Name: HU_predefined_pseudo_calibrated   Type: Temperature" << endl;
		abaqusFile << "*Initial Conditions, type=TEMPERATURE" << endl;
		for (int i=0; i<_mesh->GetNumberOfPoints(); i++){
			abaqusFile << "FEMURPART." << i+1 << "," << _mechProp->GetValue(i) * 1000 << endl; 
		}
		abaqusFile << "**" << endl;
		abaqusFile << "**" << endl;
		
		// step
		abaqusFile << "*Step, name=STEP-1, nlgeom=no" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "*Static" << endl;
		abaqusFile << "1., 1., 1e-05, 1." << endl;
		abaqusFile << "**" << endl;
		
		// loads
		abaqusFile << "**LOADS" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "** Name: hipContact Type:Concentrated force" << endl; 
		double force[3];
		_boundaryConditions->GetPoint(0, force);
		abaqusFile << "*Cload, op=NEW" << endl;
		abaqusFile << "P0,1," << force[0] << endl;
		abaqusFile << "** Name: abductor Type:Concentrated force" << endl; 
		_boundaryConditions->GetPoint(1, force);
		abaqusFile << "*Cload, op=NEW" << endl;
		abaqusFile << "P1, 1," << -force[2] << endl;
		abaqusFile << "P1, 2," << -force[0] << endl;
		abaqusFile << "P1, 3," << +force[1] << endl;
		abaqusFile << "** Name: tfsProx Type:Concentrated force" << endl; 
		_boundaryConditions->GetPoint(2, force);
		abaqusFile << "*Cload, op=NEW" << endl; 
		abaqusFile << "P1, 1," << -force[2] << endl;
		abaqusFile << "P1, 2," << -force[0] << endl;
		abaqusFile << "P1, 3," << +force[1] << endl;
		abaqusFile << "** Name: tfsDist Type:Concentrated force" << endl; 
		_boundaryConditions->GetPoint(3, force);
		abaqusFile << "*Cload, op=NEW" << endl; 
		abaqusFile << "P1, 1," << -force[2] << endl;
		abaqusFile << "P1, 2," << -force[0] << endl;
		abaqusFile << "P1, 3," << +force[1] << endl;
		abaqusFile << "** Name: vastLat Type:Concentrated force" << endl; 
		_boundaryConditions->GetPoint(4, force);
		abaqusFile << "*Cload, op=NEW" << endl; 
		abaqusFile << "P2, 1," << -force[2] << endl;
		abaqusFile << "P2, 2," << -force[0] << endl;
		abaqusFile << "P2, 3," << +force[1] << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "**" << endl;

		// boundary conditions
		abaqusFile << "** BOUNDARY CONDITIONS" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "** Name: head Type: Displacement/Rotation" << endl;
		abaqusFile << "*Boundary" << endl;
		abaqusFile << "P0, 2, 2" << endl;
		abaqusFile << "P0, 3, 3" << endl;
		abaqusFile << "** Name: knee Type: Displacement/Rotation" << endl;
		abaqusFile << "*Boundary" << endl;
		abaqusFile << "F3, 1, 1" << endl;
		abaqusFile << "F3, 2, 2" << endl;
		abaqusFile << "F3, 3, 3" << endl;
		abaqusFile << "** Name: lateral Type: Displacement/Rotation" << endl;
		abaqusFile << "*Boundary" << endl;
		abaqusFile << "F5, 3, 3" << endl;
		abaqusFile << "**" << endl;

		// output
		abaqusFile << "** OUTPUT REQUESTS" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "*Restart, write, frequency=0" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "** FIELD OUTPUT: F-Output-1" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "*Output, field" << endl;
		abaqusFile << "*Node Output" << endl;
		//abaqusFile << "CF, NT, RF, U" << endl; // concentrated force, nodal temperature, reaction forces, displacements
		abaqusFile << "NT, U" << endl;
		abaqusFile << "*Element Output, directions=YES" << endl;
		//abaqusFile << "LE, PE, PEEQ, PEMAG, S" << endl; // all logaritmic strain components, all plastic strain components, equivalent plastic strain, plastic strain magnitude, stresses, all stress components
		abaqusFile << "S, MISES, E, EP" << endl; 
		abaqusFile << "*Contact Output" << endl;
		abaqusFile << "CDISP, CSTRESS" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "** HISTORY OUTPUT: H-Output-1" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "*Output, history, variable=PRESELECT" << endl;
		abaqusFile << "** OUTPUT REQUESTS" << endl;
		abaqusFile << "*End Step" << endl;

		// close file
		outFile.close();

	}

	void FemAbaqusInpWriter::WriteNodeFallingInp(){
	
		QFile outFile(_fileName);
		outFile.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream abaqusFile(&outFile);
			
		// header
		abaqusFile << "*Heading" << endl;
		abaqusFile << "*Preprint, echo=NO, model=NO, history=NO, contact=NO" << endl;
		abaqusFile << "**" << endl;

		//part
		abaqusFile << "** PART" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "*Part, name=FEMUR" << endl;

		// nodes
		abaqusFile << "**" << endl;
		abaqusFile << "*Node" << endl;
		for (int i=0; i<_mesh->GetNumberOfPoints(); i++){
			double coord[3];
			_mesh->GetPoint(i, coord);
			abaqusFile << i+1 << ", " << coord[0] <<  ", " << coord[1] << ", " << coord[2] << endl;
		}
		
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
				if (a == _elementType-1){
					abaqusFile << endl;
				}
			}
		}
		
		// Nset
		abaqusFile << "*Nset, nset=FULLBONE, generate" << endl;
		abaqusFile << "1, " << _mesh->GetNumberOfPoints() << ", 1" << endl;
		// Elset
		abaqusFile << "*Elset, elset=FULLBONE_ELEM,generate " << endl;
		abaqusFile << "1, " << _mesh->GetNumberOfCells() << ", 1" << endl;
		
		// section
		abaqusFile << "**Section: FEMURSECTION" << endl;
		abaqusFile << "*Solid Section, elset=FULLBONE_ELEM, material=BONE" << endl; 
		abaqusFile << "," << endl;
		abaqusFile << "*End Part" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "**" << endl;
		
		// assembly
		abaqusFile << "*Assembly, name=ASSEMBLY" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "**" << endl;
		
		// instance
		abaqusFile << "*Instance, name=FEMURPART, part=FEMUR" << endl;
		abaqusFile << "*End Instance" << endl; 
		
		// Nset
		abaqusFile << "*Nset, nset=FULLBONE, instance=FEMURPART" << endl;
		abaqusFile << "1, " << _mesh->GetNumberOfPoints() << ", 1" << endl;
		// Elset
		abaqusFile << "*Elset, elset=FULLBONE_ELEM, instance=FEMURPART" << endl;
		abaqusFile << "1, " << _mesh->GetNumberOfCells() << ", 1" << endl;
		
		// Application points and coordinate system tranformation
		abaqusFile << "*Nset, nset=P0, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(0) << "," << endl;
		abaqusFile << "*Nset, nset=P1, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(1) << "," << endl; 
		abaqusFile << "*Nset, nset=P2, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(2) << "," << endl; 
		abaqusFile << "*Nset, nset=F1, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(3) << "," << endl;
		abaqusFile << "*Nset, nset=F3, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(4) << "," << endl;
		abaqusFile << "*Nset, nset=F5, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(5) << "," << endl;
		
		abaqusFile << "*Nset, nset=_T-LOCAL_BC, internal" << endl;
		abaqusFile << "P0," << endl;
		abaqusFile << "P1," << endl;
		abaqusFile << "F1," << endl;
		abaqusFile << "F3," << endl;
		abaqusFile << "*Transform, nset=_T-LOCAL_BC" << endl; 
		double C[3];
		_boundaryConditions->GetPoint(1,C);
		double E[3];
		_boundaryConditions->GetPoint(2,E);
		abaqusFile << C[0] << "," << C[1] << "," << C[2] << "," << E[0] << "," << E[1] << "," << E[2] << endl;
		
		/*
		abaqusFile << "*Nset, nset=_T-LOCAL_FCS, internal" << endl; 
		abaqusFile << "P1," << endl;
		abaqusFile << "P2," << endl;
		abaqusFile << "*Transform, nset=_T-LOCAL_FCS" << endl;
		double C[3];
		_boundaryConditions->GetPoint(7,C);
		double D[3];
		_boundaryConditions->GetPoint(8,D);
		abaqusFile << C[0] << "," << C[1] << "," << C[2] << "," << D[0] << "," << D[1] << "," << D[2] << endl;
		*/
		abaqusFile << "*End Assembly" << endl;
		
		// material properties
		abaqusFile << "*Material, name=BONE" << endl;
		abaqusFile << "*Elastic " << endl;
		abaqusFile << "50, 0.3,50," << endl;
		abaqusFile << "30000, 0.3,30000," << endl;
		abaqusFile << "**PREDEFINED CONDITIONS" << endl;
		abaqusFile << "** Name: HU_predefined_pseudo_calibrated   Type: Temperature" << endl;
		abaqusFile << "*Initial Conditions, type=TEMPERATURE" << endl;
		for (int i=0; i<_mesh->GetNumberOfPoints(); i++){
			abaqusFile << "FEMURPART." << i+1 << "," << _mechProp->GetValue(i) * 1000 << endl; 
		}
		abaqusFile << "**" << endl;
		abaqusFile << "**" << endl;
		
		// step
		abaqusFile << "*Step, name=STEP-1, nlgeom=no" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "*Static" << endl;
		abaqusFile << "1., 1., 1e-05, 1." << endl;
		abaqusFile << "**" << endl;
		
		// loads
		abaqusFile << "**LOADS" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "** Name: lateralForce Type:Concentrated force" << endl; 
		double force[3];
		_boundaryConditions->GetPoint(0, force);
		abaqusFile << "*Cload, op=NEW" << endl;
		abaqusFile << "P0, 1," << -force[0] << endl;
		abaqusFile << "P0, 2," << -force[1] << endl;
		abaqusFile << "P0, 3," << force[2] << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "**" << endl;

		// boundary conditions
		abaqusFile << "** BOUNDARY CONDITIONS" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "** Name: P1 Type: Displacement/Rotation" << endl;
		abaqusFile << "*Boundary" << endl;
		abaqusFile << "P1, 1, 1" << endl;
		abaqusFile << "P1, 2, 2" << endl;
		abaqusFile << "P1, 3, 3" << endl;
		abaqusFile << "** Name: F3 Type: Displacement/Rotation" << endl;
		abaqusFile << "*Boundary" << endl;
		abaqusFile << "F3, 1, 1" << endl;
		abaqusFile << "F3, 2, 2" << endl;
		abaqusFile << "F3, 3, 3" << endl;
		abaqusFile << "** Name: F5 Type: Displacement/Rotation" << endl;
		abaqusFile << "*Boundary" << endl;
		abaqusFile << "F5, 1, 1" << endl;
		abaqusFile << "F5, 2, 2" << endl;
		abaqusFile << "F5, 3, 3" << endl;
		abaqusFile << "**" << endl;

		// output
		abaqusFile << "** OUTPUT REQUESTS" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "*Restart, write, frequency=0" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "** FIELD OUTPUT: F-Output-1" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "*Output, field" << endl;
		abaqusFile << "*Node Output" << endl;
		//abaqusFile << "CF, NT, RF, U" << endl; // concentrated force, nodal temperature, reaction forces, displacements
		abaqusFile << "NT, U" << endl;
		abaqusFile << "*Element Output, directions=YES" << endl;
		//abaqusFile << "LE, PE, PEEQ, PEMAG, S" << endl; // all logaritmic strain components, all plastic strain components, equivalent plastic strain, plastic strain magnitude, stresses, all stress components
		abaqusFile << "S, MISES, E, EP" << endl; 
		abaqusFile << "*Contact Output" << endl;
		abaqusFile << "CDISP, CSTRESS" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "** HISTORY OUTPUT: H-Output-1" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "*Output, history, variable=PRESELECT" << endl;
		abaqusFile << "** OUTPUT REQUESTS" << endl;
		abaqusFile << "*End Step" << endl;

		// close file
		outFile.close();
	}

	void FemAbaqusInpWriter::WriteNodeStandingInp(){

		QFile outFile(_fileName);
		outFile.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream abaqusFile(&outFile);
			
		// header
		abaqusFile << "*Heading" << endl;
		abaqusFile << "*Preprint, echo=NO, model=NO, history=NO, contact=NO" << endl;
		abaqusFile << "**" << endl;

		//part
		abaqusFile << "***** PART SECTION *****" << endl;
		abaqusFile << "*Part, name=FEMUR" << endl;
		abaqusFile << "**" << endl;

		// nodes
		abaqusFile << "*Node" << endl;
		for (int i=0; i<_mesh->GetNumberOfPoints(); i++){
			double coord[3];
			_mesh->GetPoint(i, coord);
			abaqusFile << i+1 << ", " << coord[0] <<  ", " << coord[1] << ", " << coord[2] << endl;
		}
		
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
				if (a == _elementType-1){
					abaqusFile << endl;
				}
			}
		}
		abaqusFile << "**" << endl;
		
		// Nset generation
		abaqusFile << "*Nset, nset=FULLBONE, generate" << endl;
		abaqusFile << "1, " << _mesh->GetNumberOfPoints() << ", 1" << endl;
		abaqusFile << "**" << endl;
		
		// Elset generation
		abaqusFile << "*Elset, elset=FULLBONE_ELEM,generate " << endl;
		abaqusFile << "1, " << _mesh->GetNumberOfCells() << ", 1" << endl;
		abaqusFile << "**" << endl;
		
		// section
		abaqusFile << "*Solid Section, elset=FULLBONE_ELEM, material=BONE" << endl; 
		abaqusFile << "**" << endl;
		abaqusFile << "*End Part" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "**" << endl;
		
		// assembly
		abaqusFile << "***** ASSEMBLY SECTION *****" << endl;
		abaqusFile << "*Assembly, name=ASSEMBLY" << endl;
		abaqusFile << "**" << endl;
				
		// instance
		abaqusFile << "*Instance, name=FEMURPART, part=FEMUR" << endl;
		abaqusFile << "*End Instance" << endl;
		abaqusFile << "**" << endl;
		
		// Nset assignment
		abaqusFile << "*Nset, nset=FULLBONE_NODE, instance=FEMURPART" << endl;
		abaqusFile << "1, " << _mesh->GetNumberOfPoints() << ", 1" << endl;
		abaqusFile << "**" << endl;
		
		// Elset assignment
		abaqusFile << "*Elset, elset=FULLBONE_ELEM, instance=FEMURPART" << endl;
		abaqusFile << "1, " << _mesh->GetNumberOfCells() << ", 1" << endl;
		abaqusFile << "**" << endl;
		
		// Application points  
		abaqusFile << "*Nset, nset=P0, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(0) << "," << endl;
		abaqusFile << "*Nset, nset=P1, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(1) << "," << endl; 
		abaqusFile << "*Nset, nset=P2, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(2) << "," << endl; 
		abaqusFile << "*Nset, nset=F1, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(3) << "," << endl;
		abaqusFile << "*Nset, nset=F3, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(4) << "," << endl;
		abaqusFile << "*Nset, nset=F5, instance=FEMURPART" << endl;
		abaqusFile << _appPoints->GetValue(5) << "," << endl;
		abaqusFile << "**" << endl;

		// cordinate system transformation
		abaqusFile << "*Nset, nset=_T-LOCAL_BC, internal" << endl;
		abaqusFile << "P0, F3, F5" << endl;
		abaqusFile << "*Transform, nset=_T-LOCAL_BC" << endl; 
		double A[3];
		_boundaryConditions->GetPoint(1,A);
		double B[3];
		_boundaryConditions->GetPoint(2,B);
		abaqusFile << A[0] << "," << A[1] << "," << A[2] << "," << B[0] << "," << B[1] << "," << B[2] << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "*Nset, nset=_T-LOCAL_FCS, internal" << endl; 
		abaqusFile << "P1, P2" << endl;
		abaqusFile << "*Transform, nset=_T-LOCAL_FCS" << endl;
		double C[3];
		_boundaryConditions->GetPoint(3,C);
		double D[3];
		_boundaryConditions->GetPoint(4,D);
		abaqusFile << C[0] << "," << C[1] << "," << C[2] << "," << D[0] << "," << D[1] << "," << D[2] << endl;
		
		abaqusFile << "*End Assembly" << endl;
		abaqusFile << "**" << endl;
		
		// material properties
		abaqusFile << "***** MATERIAL SECTION ***** " << endl;
		abaqusFile << "*Material, name=BONE" << endl;
		abaqusFile << "*Elastic " << endl;
		abaqusFile << "50, 0.3,50," << endl;
		abaqusFile << "30000, 0.3,30000," << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "*Initial Conditions, type=TEMPERATURE" << endl;
		for (int i=0; i<_mesh->GetNumberOfPoints(); i++){
			abaqusFile << "FEMURPART." << i+1 << "," << _mechProp->GetValue(i) * 1000 << endl; 
		}
		abaqusFile << "**" << endl;
		abaqusFile << "**" << endl;
		
		// step
		abaqusFile << "***** STEP SECTION *****" << endl;
		abaqusFile << "*Step, name=STEP-1, nlgeom=no" << endl;
		abaqusFile << "*Static" << endl;
		abaqusFile << "1., 1., 1e-05, 1." << endl;
		abaqusFile << "**" << endl;
		
		// loads
		abaqusFile << "**LOADS" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "** Name: hipContact Type:Concentrated force" << endl; 
		double force[3];
		_boundaryConditions->GetPoint(0, force);
		abaqusFile << "*Cload, op=NEW" << endl;
		abaqusFile << "P0,1," << force[0] << endl;
		abaqusFile << "**" << endl;

		// boundary conditions
		abaqusFile << "** BOUNDARY CONDITIONS" << endl;
		abaqusFile << "** Name: head Type: Displacement/Rotation" << endl;
		abaqusFile << "*Boundary" << endl;
		abaqusFile << "P0, 2, 2" << endl;
		abaqusFile << "P0, 3, 3" << endl;
		abaqusFile << "** Name: knee Type: Displacement/Rotation" << endl;
		abaqusFile << "*Boundary" << endl;
		abaqusFile << "F3, 1, 1" << endl;
		abaqusFile << "F3, 2, 2" << endl;
		abaqusFile << "F3, 3, 3" << endl;
		abaqusFile << "** Name: lateral Type: Displacement/Rotation" << endl;
		abaqusFile << "*Boundary" << endl;
		abaqusFile << "F5, 3, 3" << endl;
		abaqusFile << "**" << endl;

		// output
		abaqusFile << "** OUTPUT REQUESTS" << endl;
		abaqusFile << "*Restart, write, frequency=0" << endl;
		abaqusFile << "*Output, field" << endl;
		abaqusFile << "*Node Output" << endl;
		//abaqusFile << "CF, NT, RF, U" << endl; // concentrated force, nodal temperature, reaction forces, displacements
		abaqusFile << "NT, U" << endl;
		abaqusFile << "*Element Output, directions=YES" << endl;
		//abaqusFile << "LE, PE, PEEQ, PEMAG, S" << endl; // all logaritmic strain components, all plastic strain components, equivalent plastic strain, plastic strain magnitude, stresses, all stress components
		abaqusFile << "S, MISES, E, EP" << endl; 
		abaqusFile << "*Contact Output" << endl;
		abaqusFile << "CDISP, CSTRESS" << endl;
		abaqusFile << "**" << endl;
		abaqusFile << "*Output, history, variable=PRESELECT" << endl;
		abaqusFile << "*End Step" << endl;

		// close file
		outFile.close();

	
	}
}
