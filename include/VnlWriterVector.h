/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef VNLWRITERVECTOR_H
#define VNLWRITERVECTOR_H

#include <VnlWriter.h>

#include <vnl/vnl_vector.h>

namespace vnl{

	/**
	* Writes the vnl vector into a text file
	*/

	class VnlWriterVector: public VnlWriter{

	public:
		
		// constructor/destructor
		VnlWriterVector();
		~VnlWriterVector();

		// accessor
		/**
		* Sets the vnl vector
		*/
		void SetVnlVector(vnl_vector<double> vector) {_vector = vector;}


		// overwritten virtual function
		void Update();

	protected:
		vnl_vector<double> _vector;

					
	};
}
#endif //VNLWRITERVECTOR_H