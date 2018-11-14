/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef VNLREADERVECTOR_H
#define VNLREADERVECTOR_H

#include <VnlReader.h>

#include <vnl/vnl_vector.h>

namespace vnl{

	/**
	* Reads the vector from a text file
	*/

	class VnlReaderVector: public VnlReader{

	public:
		
		// constructor/destructor
		VnlReaderVector();
		~VnlReaderVector();

		// accessor
		/**
		* Gets the vector
		*/
		vnl_vector<double> GetVnlVector() {return _vector;}
		
		// overwritten virtual function
		void Update();

	protected:
		vnl_vector<double> _vector;
				
	};
}
#endif //VNLREADERVECTOR_H