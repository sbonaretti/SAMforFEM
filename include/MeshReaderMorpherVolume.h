/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef MESHREADERMORPHERVOLUME_H
#define MESHREADERMORPHERVOLUME_H

#include <MeshReader.h>

namespace mesh{
	
	/**
	* Reads the .cdb file containing a quadratic tetra volume mesh and used by the Ansys Morpher 
	* 
	* In each element, nodes are written in the following order: 0 1 2 3 4 5 6 7 8 9 (and not: 0 1 2 3 4 7 5 6 8 9)
	*/

	class MeshReaderMorpherVolume: public MeshReader{
	
	public: 

		// constructor/destructor
		MeshReaderMorpherVolume();
		~MeshReaderMorpherVolume();
		
		// overwritten virtual function
		void Update();
		
	};
}
#endif //MESHREADERMORPHERVOLUME_H 