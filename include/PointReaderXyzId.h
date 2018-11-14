/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef POINTREADERXYZID_H
#define POINTREADERXYZID_H

#include <PointReader.h>

#include <vtkDoubleArray.h>

namespace points{
	class PointReaderXyzId: public PointReader{
		
	public:

		// constructor/destructor
		PointReaderXyzId();
		~PointReaderXyzId();
		
		// member function
		vtkDoubleArray* GetIdVector() {return _idVector;}
	
		// overwritten virtual function
		void Update();

	private:

		vtkDoubleArray* _idVector;

	};
}
#endif //POINTREADERXYZID_H

