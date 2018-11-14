/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef STATISTICSENTROPY_H
#define STATISTICSENTROPY_H

#include <vnl/vnl_vector.h>


namespace statistics{

	/**
	* Class for the calculation of the entropy
	*/

	class StatisticsEntropy{

	public:

		// constructor/destructor
		StatisticsEntropy();
		~StatisticsEntropy();

		// accessors
		/**
		* Sets the vector
		*/
		void SetVector (vnl_vector<double> vector) {_vector = vector;}
		/**
		* Gets the entropy
		*/
		double GetEntropy(){return _entropy;}

		// methods
		/**
		* Calculates the Shannon Entropy (- sum(pi*log(pi)))
		*/
		void CalculateShannonEntropy();

	
	protected:
		vnl_vector<double> _vector;
		double _entropy;

		};
}
#endif //STATISTICSENTROPY_H



		

