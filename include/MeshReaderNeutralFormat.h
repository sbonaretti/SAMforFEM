/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef MESHREADERNEUTRALFORMAT_H
#define MESHREADERNEUTRALFORMAT_H

#include <MeshReader.h>

namespace mesh{

	/**
	* Reads the .inp file created by NetGen containing a quadratic tetra volume mesh
	*/
	
	class MeshReaderNeutralFormat: public MeshReader{
	
	public: 

		// constructor/destructor
		MeshReaderNeutralFormat();
		~MeshReaderNeutralFormat();
		
		// overwritten virtual function
		void Update();
		
	};
}
#endif //MESHREADERNEUTRALFORMAT_H 