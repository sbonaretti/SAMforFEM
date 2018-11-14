/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <FemAssignerElements.h>

#include <MeshReaderAbaqus.h>
#include <ImageHandler.h>

#include <vtkGenericCell.h>
#include <vtkIdList.h>

using namespace image;
using namespace mesh;


namespace fem{

	// constructor
	FemAssignerElements::FemAssignerElements(){
	}

	// destructor
	FemAssignerElements::~FemAssignerElements(){
	}

	
	// overwritten virtual functions
	void FemAssignerElements::Update(){
		
		
		// nodes assignment
		vtkDoubleArray* rAsh = vtkDoubleArray::New();
		vtkDoubleArray* rApp = vtkDoubleArray::New(); // for each node (for the calculation of E)
		vtkDoubleArray* nodesYoungModulus = vtkDoubleArray::New();
		
		for (int i=0; i<_greyLevel->GetNumberOfTuples(); i++){
			
			// rhoAsh
			double rhoAsh = _greyLevel->GetValue(i);
			rAsh->InsertNextValue(rhoAsh);
			
			// rhoApp
			double rhoApp = rhoAsh/0.60;  
			rApp->InsertNextValue(rhoApp);
			
			// young's modulus [GPa]
			double E = _assignmentLawOne * pow(rhoApp,_assignmentLawTwo); //[mg/mm3]
			if (E < 0.050) 
				E = 0.050;
			nodesYoungModulus->InsertNextValue(E);
		}

		double range[2];
		rAsh->GetRange(range);
		std::cout << "node rho ash: " << range[0] << ' ' << range[1] << std::endl;

		rApp->GetRange(range);
		std::cout << "node rho app: " << range[0] << ' ' << range[1] << std::endl;

		nodesYoungModulus->GetRange(range);
		std::cout << "node young modulus: " << range[0] << ' ' << range[1] << std::endl;
	
			
				
		// element assignment as average on the 10 nodes		
		std::cout << "number of cells: " << _mesh->GetNumberOfCells() << std::endl; 	
		vtkGenericCell* cell = vtkGenericCell::New();
		vtkIdList* idList;
		
		// apparent density
		_rApp = vtkDoubleArray::New(); // for each element (for the GetRhoApp() )
		for (int i=0; i<_mesh->GetNumberOfCells(); i++){
		
			// get the element
			_mesh->GetCell(i, cell);
			idList = cell->GetPointIds();
			
			double Eelement = 0.0;
			for (int a=0; a<idList->GetNumberOfIds(); a++){
				Eelement += rApp->GetValue(idList->GetId(a));
			}
			Eelement /= idList->GetNumberOfIds();
			if (Eelement < 0.050 || Eelement != Eelement) // Eelement != Eelement detects if the variable is nan (not a number)
				Eelement = 0.050; // the min bone Young's modulus is considered to be 0.050GPa 
			_rApp->InsertNextValue( Eelement);

		}
		_rApp->GetRange(range);
		std::cout << "element apparent density: " << range[0] << ' ' << range[1] << std::endl;
		
		// Young's modulus		
		for (int i=0; i<_mesh->GetNumberOfCells(); i++){
		
			// get the element
			_mesh->GetCell(i, cell);
			idList = cell->GetPointIds();
			
			double Eelement = 0.0;
			for (int a=0; a<idList->GetNumberOfIds(); a++){
				Eelement += nodesYoungModulus->GetValue(idList->GetId(a));
			}
			Eelement /= idList->GetNumberOfIds();
			if (Eelement < 0.050 || Eelement != Eelement) // Eelement != Eelement detects if the variable is nan (not a number)
				Eelement = 0.050; // the min bone Young's modulus is considered to be 0.050GPa 
			_youngModulus->InsertNextValue( Eelement);

		}
		_youngModulus->GetRange(range);
		std::cout << "element young modulus: " << range[0] << ' ' << range[1] << std::endl;

		
	
	}
}