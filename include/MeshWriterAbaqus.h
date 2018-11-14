/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef MESHWRITERABAQUS_H
#define MESHWRITERABAQUS_H

#include <MeshWriter.h>

namespace mesh{

	/**
	* Writes the .inp file containing a quadratic tetra volume mesh
	* 
	* In each element, nodes are written in the following order: 0 1 2 3 4 5 6 7 8 9 
	*/

	class MeshWriterAbaqus: public MeshWriter{

	public:
		
		// constructor/destructor
		MeshWriterAbaqus();
		~MeshWriterAbaqus();

		// overwritten virtual function
		void Update();

	private:
		void findClosestNode (double point[3], int &minIndex);

	protected:
		int _elementType;
		QString _abaqusElementType;
		double _poissonRatio; 
		int _loadPoint1; // index of the mesh node closest to the input point
		int _loadPoint2;
		
				
	};
}
#endif //MESHWRITERABAQUS_H