/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef MESHREADERABAQUS_H
#define MESHREADERABAQUS_H

#include <MeshReader.h>

namespace mesh{
	
	/**
	* Reads the .cdb file containing a quadratic tetra volume mesh and used by the Ansys Morpher 
	* 
	* In each element, nodes are read in the following order: 0 1 2 3 4 7 5 6 8 9 
	*/

	class MeshReaderAbaqus: public MeshReader{
	
	public: 

		// constructor/destructor
		MeshReaderAbaqus();
		~MeshReaderAbaqus();
		
		// overwritten virtual function
		void Update();

	protected:
		int _elementType;
		QString _abaqusElementType;
		
	};
}
#endif //MESHREADERABAQUS_H 