/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef POINTREADERXYZ_H
#define POINTREADERXYZ_H

#include <PointReader.h>

namespace points{
	class PointReaderXyz: public PointReader{
		
	public:

		// constructor/destructor
		PointReaderXyz();
		~PointReaderXyz();

		// overwritten virtual function
		void Update();
		
	};
}
#endif //POINTREADERXYZ_H

