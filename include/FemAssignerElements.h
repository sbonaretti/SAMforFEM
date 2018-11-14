/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef FEMASSIGNERELEMENTS_H
#define FEMASSIGNERELEMENTS_H

#include <FemAssigner.h>


namespace fem{

	/**
	* Assignement of the Young's modulus to the mesh elements.
	* ! Call GreyLevelAssignmentUpdate() before Update() 
	*/


	class FemAssignerElements : public FemAssigner {
	
	public:

		// constructor/destructor
		FemAssignerElements();
		~FemAssignerElements();

			
		// overwritten virtual function
		/**
		* Execute the assignment
		*/		
		void Update();

	};

}
#endif //FEMASSIGNERELEMENTS_H 
