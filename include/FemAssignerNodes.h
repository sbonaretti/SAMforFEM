/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef FEMASSIGNERNODES_H
#define FEMASSIGNERNODES_H

#include <FemAssigner.h>


namespace fem{

	/**
	* Assignement of the Young's modulus to the mesh nodes.
	* ! Call GreyLevelAssignmentUpdate() before Update() 
	*/


	class FemAssignerNodes : public FemAssigner {
	
	public:

		// constructor/destructor
		FemAssignerNodes();
		~FemAssignerNodes();

		// overwritten virtual function
		/**
		* Execute the assignment
		*/		
		void Update();

	};

}
#endif //FEMASSIGNERNODES_H 
