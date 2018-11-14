/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef MESHREADERMORPHERSURFACE_H
#define MESHREADERMORPHERSURFACE_H

#include <MeshReader.h>

namespace mesh{

	/**
	* Reads the .cdb file created by the Ansys Morpher containing a surface mesh
	*/
	class MeshReaderMorpherSurface: public MeshReader{
	
	public: 

		// constructor/destructor
		MeshReaderMorpherSurface();
		~MeshReaderMorpherSurface();
		
		// overwritten virtual function
		void Update();
		
	};
}
#endif //MESHREADERMORPHERSURFACE_H 