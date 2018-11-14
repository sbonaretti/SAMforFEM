/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <FemAssignerNodes.h>

#include <MeshReaderAbaqus.h>
#include <ImageHandler.h>

#include <vtkPoints.h>
#include <vtkIdList.h>

using namespace image;
using namespace mesh;


namespace fem{

	// constructor
	FemAssignerNodes::FemAssignerNodes(){
	}

	// destructor
	FemAssignerNodes::~FemAssignerNodes(){
	}

	
	// overwritten virtual functions
	void FemAssignerNodes::Update(){

		std::cout << std::endl;
		std::cout << "-- Update from FemAssignerNodes --" << std::endl;

		_rAsh = vtkDoubleArray::New();
		_rApp = vtkDoubleArray::New();
		
		for (int i=0; i<_greyLevel->GetNumberOfTuples(); i++){
			
			// rhoAsh
			double rhoAsh = _greyLevel->GetValue(i);
			_rAsh->InsertNextValue(rhoAsh);
			
			// rhoApp
			double rhoApp = rhoAsh/0.60;
			_rApp->InsertNextValue(rhoApp);
			
			// young's modulus [GPa]
			double E = _assignmentLawOne * pow(rhoApp,_assignmentLawTwo); //[mg/mm3]
			
			// ben's old
			//if (E < 0.500 || E != E)
			//	E = 0.500;

			// ben's new 
			//if (E < 0.010 || E != E)
			//	E = 0.010;
			
			// my old
			if (E < 0.050 || E != E) // E != E detects if the variable is nan (not a number)
				E = 0.050;
			_youngModulus->InsertNextValue(E);
		}

		double range[2];
		_rAsh->GetRange(range);
		std::cout << "rho ash: " << range[0] << ' ' << range[1] << std::endl;

		_rApp->GetRange(range);
		std::cout << "rho app: " << range[0] << ' ' << range[1] << std::endl;
		
		_youngModulus->GetRange(range);
		std::cout << "young modulus: " << range[0] << ' ' << range[1] << std::endl;

	}
}

	