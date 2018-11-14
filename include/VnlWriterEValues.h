/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef VNLWRITEREVALUES_H
#define VNLWRITEREVALUES_H

#include <VnlWriter.h>

#include <vnl/vnl_vector.h>

namespace vnl{

	/**
	* Writes the vnl id, eigenvalues and normalized eigenvalues into a text file
	*/

	class VnlWriterEValues: public VnlWriter{

	public:
		
		// constructor/destructor
		VnlWriterEValues();
		~VnlWriterEValues();

		// accessor
		/**
		* Sets the eigenvalues vector
		*/
		void SetVnlEValues(vnl_vector<double> eValues) {_eValues = eValues;}
		/**
		* Sets the normalized eigenvalues vector
		*/
		void SetVnlNormalizedEValues(vnl_vector<double> normalizedEValues) {_normalizedEValues = normalizedEValues;}


		// overwritten virtual function
		void Update();

	protected:
		vnl_vector<double> _eValues;
		vnl_vector<double> _normalizedEValues;

					
	};
}
#endif //VNLWRITEREVALUES_H