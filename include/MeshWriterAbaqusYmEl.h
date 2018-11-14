/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef MESHWRITERABAQUSYMEL_H
#define MESHWRITERABAQUSYMEL_H

#include <MeshWriter.h>

#include <vtkDoubleArray.h>

namespace mesh{

	/**
	* Writes the .inp file containing a quadratic tetra volume mesh
	* 
	* In each element, nodes are written in the following order: 0 1 2 3 4 5 6 7 8 9 
	*/


	class MeshWriterAbaqusYmEl: public MeshWriter{

	public:
		
		// constructor/destructor
		MeshWriterAbaqusYmEl();
		~MeshWriterAbaqusYmEl();

		void SetYoungModulus(vtkDoubleArray* youngModulus) {_youngModulus = youngModulus;}

		// overwritten virtual function
		void Update();

	protected:
		int _elementType;
		QString _abaqusElementType;
		vtkDoubleArray* _youngModulus;
		double _poissonRatio; 
		

			
	};
}
#endif //MESHWRITERABAQUSYMEL_H