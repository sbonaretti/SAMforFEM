/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef POINTWRITERXYZVECTOR_H
#define POINTWRITERXYZVECTOR_H

#include <PointWriter.h>

namespace points{

	/**
	* Writes the .txt file containing N points in 3N x 1 form 
	*/

	class PointWriterXyzVector: public PointWriter{
		
	public:

		// constructor/destructor
		PointWriterXyzVector();
		~PointWriterXyzVector();

		// overwritten virtual function
		void Update();
		
	};
}
#endif //POINTWRITERXYZVECTOR_H