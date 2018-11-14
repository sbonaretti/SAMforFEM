/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef VNLREADEREVALUES_H
#define VNLREADEREVALUES_H

#include <VnlReader.h>

#include <vnl/vnl_vector.h>

namespace vnl{

	/**
	* Reads the vnl id, eigenvalues and normalized eigenvalues from a text file
	*/

	class VnlReaderEValues: public VnlReader{

	public:
		
		// constructor/destructor
		VnlReaderEValues();
		~VnlReaderEValues();

		// accessor
		/**
		* Gets the eigenvalues vector
		*/
		vnl_vector<double> GetVnlEValues() {return _eValues;}
		/**
		* Gets the normalized eigenvalues vector
		*/
		vnl_vector<double> GetVnlNormalizedEValues() {return _normalizedEValues;}


		// overwritten virtual function
		void Update();

	protected:
		vnl_vector<double> _eValues;
		vnl_vector<double> _normalizedEValues;

					
	};
}
#endif //VNLREADEREVALUES_H