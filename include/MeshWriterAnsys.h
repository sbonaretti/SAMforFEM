/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef MESHWRITERANSYS_H
#define MESHWRITERANSYS_H

#include <MeshWriter.h>

namespace mesh{

	/**
	* Writes the .cdb file containing a quadratic tetra volume mesh and used by the Ansys Morpher 
	* 
	* In each element, nodes are written in the following order: 0 1 2 3 4 5 6 7 8 9 (and not: 0 1 2 3 4 7 5 6 8 9)
	*/

	class MeshWriterAnsys: public MeshWriter{

	public:
		
		// constructor/destructor
		MeshWriterAnsys();
		~MeshWriterAnsys();

		// overwritten virtual function
		void Update();

	};
}
#endif //MESHWRITERANSYS_H