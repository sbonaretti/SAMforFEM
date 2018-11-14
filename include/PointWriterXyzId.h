
/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef POINTWRITERXYZID_H
#define POINTWRITERXYZID_H

#include <PointWriter.h>

#include <vtkDoubleArray.h>

namespace points{

	/**
	* Writes the .txt file containing N points and the vector one in N x 3 form 
	*/

	class PointWriterXyzId: public PointWriter{
		
	public:

		// constructor/destructor
		PointWriterXyzId();
		~PointWriterXyzId();

		// member functions
		void SetIdVector(vtkDoubleArray* vector) {_idVector = vector;}

		// overwritten virtual function
		void Update();
		
	private:

		vtkDoubleArray* _idVector;

	};
}
#endif //POINTWRITERXYZID_H