/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <StatisticsEntropy.h>

#include <iostream>

namespace statistics{

	// constructor
	StatisticsEntropy::StatisticsEntropy(){
		_entropy = 0.0;
	}
	// destructor
	StatisticsEntropy::~StatisticsEntropy(){
	}

		// methods
	void StatisticsEntropy::CalculateShannonEntropy(){

		// normalize data
		double max =  _vector.max_value();
		for (int i=0; i<_vector.size(); i++){
			_vector(i) = _vector(i) / max;
		}

		// calculate the entropy
		for (int i=0; i<_vector.size(); i++){
			if (_vector(i) != 0){
				_entropy = _entropy + _vector(i) * std::log(_vector(i));
			}
		}
		_entropy = - _entropy;
			
	}
}