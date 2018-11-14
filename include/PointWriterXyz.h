/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef POINTWRITERXYZ_H
#define POINTWRITERXYZ_H

#include <PointWriter.h>

namespace points{

	/**
	* Writes the .txt file containing N points in N x 3 form 
	*/

	class PointWriterXyz: public PointWriter{
		
	public:

		// constructor/destructor
		PointWriterXyz();
		~PointWriterXyz();

		// overwritten virtual function
		void Update();
		
	};
}
#endif //POINTWRITERXYZ_H

